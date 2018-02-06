##################################################################
########## Effects analysis models using metafor #################
##################################################################
library(metafor)

# # Tests and comments from the paper of metafor (2010)
# gdata = acute_gene_tables_raw[["10891"]] # PGC1 in acute response
# gdata = longterm_gene_tables[["1282"]] # COL1 gene
# gdata = longterm_gene_tables[["4168"]] # A negative example
# gdata = acute_gene_tables[["5166"]] # survived rep but not meta
# # gdata = acute_gene_tables[["11326"]] # gene with significant modifiers
# # another test: make a dataset with large effects
# # gdata$yi = rnorm(nrow(gdata),mean=5)
# plot(gdata[,"yi"],gdata[,"vi"],pch=as.numeric(as.factor(gdata$tissue)))
# # knha - a correction that accounts for the uncertainty in the random effect
# res0 = rma.mv(yi,vi,random = ~ 1|gse,data=gdata,subset = (tissue=="muscle"))
# # explanation of the result above:
# #   mu - the average effect is 0.067, the CI contains zero
# res = rma.mv(yi,vi,mods = ~  training + time, random= ~1|gse ,data=gdata[gdata$tissue=="muscle",])
# # CIs of the anova stats:
# confint(res0)
# # Forest plot - very informative
# forest(res)
# # difference in tau before and after using moderators - teaches us about the
# # percentage of explained variance due to the moderators
# # when the test for residuals (QE) is signficant - we may be missing additional
# # moderators
# # We can use predict to get expected effects for new moderators:
# # Currently does not work because factors should be transformed into dummy vars.
# # Also pages 18-19 show nice figures and analysis of the predictions.
# predict(res,newmods = as.matrix(data.frame("",time=5,"",training="endurance",tissue="muscle")), addx = TRUE)
# predict(res,newmods = as.matrix(gdata))
# # look at the fitted values
# predict(res)
# # Separate by tissue
# res2 = rma(yi=yi,vi=vi,mods = ~ training + time ,data=gdata,knha=T, subset = (tissue=="muscle"))
# res1 = rma(yi=yi,vi=vi,mods = ~ training + time ,data=gdata,knha=T, subset = (tissue=="blood"))
# # Residual analysis for detecting outlier datasets
# barplot(as.numeric(rstudent(res)$z));abline(-2,0);abline(2,0)
# # Residual analysis is informative but not enough
# # case deletion diagnostics are informatice as well
# plot(influence(res))
# # funnel plots
# funnel(res0)
# # radial plots: useful for consistency analysis
# # can be used only for models without moderators
# radial(res0)
# # qq plots for the standardized residuals
# qqnorm(res0,main="random")
# qqnorm(res,main="mixed, both tissues")
# qqnorm(res1,main="blood")
# qqnorm(res2,main="muscle")
# # tests for publication bias
# regtest(res0,predictor="vi",model="lm")
# regtest(res2,predictor="vi",model="lm")
# # anova tests
# anova(res0,res)
# 
# # # Compare to the rmeta package
# # install.packages('rmeta')
# # library(rmeta)
# # blood_gdata = gdata[gdata$tissue=="blood",]
# # rmeta_res0 = meta.summaries(d=blood_gdata$yi,se = blood_gdata$vi)
# # summary(rmeta_res0)[[3]]
# # plot(rmeta_res0)
# # funnelplot(rmeta_res0)
# # gdata$yi
# # 
# # # Tests on simulated data
# # n = 5;effect = 3; hetero=0.5; effect2=6
# # vi = rep(0.5,n)
# # yi = rnorm(n,sd=vi) + effect + rnorm(n,sd=hetero)
# # rma(yi,vi)$pval
# # dummy = 1:n
# # rma.mv(yi,vi,random=~ 1|dummy)$pval
# # # vs.
# # rma(c(yi,rnorm(n,sd=vi) + rnorm(n,sd=hetero)),c(vi,vi))$pval
# # # vs.
# # ref = factor(c(rep(1,n),rep(2,n)))
# # rma.mv(c(yi,rnorm(n,sd=vi) + rnorm(n,sd=hetero)),c(vi,vi),random = ~1|ref)$pval
# # 
# # # random effects
# # beta1 = 1 ; beta2 = 0.5
# # vi = c(vi,vi)
# # yi = c(yi+beta1,yi+beta2)
# # ref = factor(c(rep(1,n),rep(2,n)))
# # rma(yi,vi)
# # rma.mv(yi,vi,rand = ~1|ref)


##################################################################
########## Downstream analyses and display items #################
##################################################################

# Get average effect per time point x training
get_subset_forest_plot<-function(gdata,tissue="all",training="all",sortby = "time",labelsby=c("gse","time"),...){
  if(!tissue=="all"){
    gdata = gdata[gdata$tissue==tissue,]
  }
  if(!training=="all"){
    gdata = gdata[gdata$training==training,]
  }
  ord = order(gdata[,sortby])
  gdata = gdata[ord,]
  res = rma.mv(yi,vi,random = ~ 1|gse,data=gdata)
  res0 = rma(yi,vi,data=gdata)
  rtest = regtest(res0)
  print(rtest)
  forest(res,slab=apply(gdata[,labelsby],1,paste,collapse=','),...)
}

get_gene_weighted_avg_pattern <-function(gdata){
  gdata$wt = 1/(sqrt(gdata$vi))
  res = by(gdata, paste(gdata$tissue,gdata$training,gdata$time,sep=','), function(x) weighted.mean(x$yi, x$wt),simplify = T)
  v = as.numeric(res);names(v) = names(res)
  return(v)
}
reorder_weighted_avg_matrix<-function(x){
  cn = colnames(x)
  cn = sapply(cn,function(x)strsplit(x,split=',')[[1]])
  ord = order(cn[1,],cn[2,],as.numeric(cn[3,]))
  x = x[,ord]
  return(x)
}
# geneset - a list of gene symbols
print_drem_matrices<-function(x,dir_name,gene_conv = entrez2symbol,min_time_points=3,geneset=NULL){
  cn = colnames(x)
  cn = sapply(cn,function(x)strsplit(x,split=',')[[1]])
  tb = table(cn[1,],cn[2,])
  for(tissue in rownames(tb)){
    xx = x[,cn[1,]==tissue]
    times = as.numeric(cn[3,cn[1,]==tissue])
    trs = cn[2,cn[1,]==tissue]
    for(tr in unique(trs)){
      if(sum(trs==tr)<min_time_points){next}
      curr_xx = xx[,trs==tr]
      curr_times = times[trs==tr]
      colnames(curr_xx) = paste(curr_times,"h",sep='')
      rownames(curr_xx) = gene_conv[rownames(curr_xx)]
      if(!is.null(geneset)){curr_xx = curr_xx[intersect(rownames(curr_xx),geneset),]}
      fname = paste(dir_name,tissue,"_",tr,'.txt',sep='')
      curr_xx = cbind(rownames(curr_xx),curr_xx)
      colnames(curr_xx)[1] = "UID"
      write.table(curr_xx,file=fname,quote=F,sep="\t",col.names = T,row.names = F)
    }
  }
}
plot_gene_pattern<-function(x,errs = NULL,tosmooth=T,min_time_points=3,
                            tr2col=c(untrained="red",endurance="blue",resistance="green",both="purple"),
                            legend.x="topright",legend.cex=0.8,legend.ncol=2,
                            mfrow = c(1,2),main_prefix = "acute",y_lim_add = 0.5,y_lim_min = 1.5,
                            ...){
  inds = !grepl("yoga|treatment",names(x))
  x = x[inds]
  if(!is.null(errs)){errs = errs[inds]}
  cn = names(x)
  cn = sapply(cn,function(x)strsplit(x,split=',')[[1]])
  tb = table(cn[1,],cn[2,])
  if(!is.null(mfrow)){par(mfrow=mfrow)}
  for(tissue in rownames(tb)){
    xx = x[cn[1,]==tissue]
    if(!is.null(errs)){yy = errs[cn[1,]==tissue]}
    if(length(xx)<min_time_points){next}
    times = as.numeric(cn[3,cn[1,]==tissue])
    trs = cn[2,cn[1,]==tissue]
    ylim = c(min(xx)-y_lim_add,max(xx)+y_lim_add)
    ylim[1] = min(-y_lim_min,ylim[1])
    ylim[2] = max(y_lim_min,ylim[2])
    plot(xx,x=times,col="white",las=2,main=paste(main_prefix,tissue),xlim=c(-1,max(times)),ylim=ylim,ylab="log2 Fold Change",xlab="time")
    for(tr in unique(trs)){
      if(tosmooth){xx1 = c(0,smooth(xx[trs==tr]))}
      else{xx1 = c(0,xx[trs==tr])}
      if(length(xx1)<min_time_points){next}
      tt1 = c(-1,times[trs==tr])
      lines(xx1,x=tt1,type='l',lwd=2,col=tr2col[tr],pch=20)
      if(!is.null(errs)){
        yy1 = c(0.01,yy[trs==tr])
        arrows(tt1, xx1-yy1, tt1, xx1+yy1, length=0.05, angle=90, code=3,col=tr2col[tr])
      }
    }
    legend(x=legend.x,legend = unique(trs),fill=tr2col[unique(trs)],cex=legend.cex,ncol=legend.ncol)  
  }
}
library(cluster)
pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
# x: rows are genes
perform_gene_clustering<-function(x,fun=pam1,standardize_genes=T,num_pcs=5){
  xx = x
  if(standardize_genes){xx = t(apply(x,1,function(x)(x-mean(x))/(sd(x))))}
  if(num_pcs>1){xx = prcomp(xx)$x[,1:num_pcs]}
  gs.pam.xx <- clusGap(xx, FUN = fun, K.max = 10, B = 20)
  num_c = maxSE(gs.pam.xx[[1]][,3],gs.pam.xx[[1]][,4])
  cc = fun(xx,num_c)$cluster
  return(cc)
}
cluster_homogeneities<-function(x,cl,...){
  hs = c()
  for(cc in unique(cl)){
    cc = as.character(cc)
    if(sum(cl==cc)==1){
      hs[cc] = 1
    }
    else{
      rhos = cor(t(x[cl==cc,]),...)
      hs[cc] = mean(rhos[lower.tri(rhos)])
    }
  }
  names(hs) = as.character(unique(cl))
  return(hs)
}
standardize_rows = function(x){
  xx = t(apply(x,1,function(x)(x-mean(x))/(sd(x))))
  return(x)
}
options(expressions=10000)
cluster_genes_by_homogeneity <- function(x,thr=0.5,func=kmeans,...){
  if(is.null(dim(x))){
    return (1)
  }
  if(nrow(x)<2){
    return(1)
  }
  v = rep(1,nrow(x));names(v) = rownames(x)
  hs = cluster_homogeneities(x,v)
  if(hs>=thr){
    return(v)
  }
  k=2
  while(all(hs<thr) && k < nrow(x)){
    v = as.numeric(func(x,k)[[1]])
    names(v) = rownames(x)
    hs = cluster_homogeneities(x,v)
    k = k+1
  }
  if(k==nrow(x)){
    v = 1:nrow(x);names(x) = rownames(x)
    return(v)
  }
  for(cl in names(hs)){
    if(hs[cl]>=thr){next}
    newx = x[v==cl,]
    newx_v = as.numeric(cluster_genes_by_homogeneity(newx,thr,func))
    newx_v = as.numeric(newx_v) + max(as.numeric(v))
    names(newx_v) = rownames(newx)
    v[names(newx_v)] = newx_v
    # print(table(v))
    # print("###")
    # print(table(newx_v))
    # print("###")
    # print(cluster_homogeneities(newx,newx_v))
  }
  v = as.numeric(as.factor(v))
  names(v) = rownames(x)
  return(v)
}

# mm = standardize_rows(m)
# vv = cluster_genes_by_homogeneity(mm,0.5,pam1)
# sort(table(vv))
# cluster_homogeneities(mm,vv)

#########################################################

##################################################################
############## Mixed effect models using lme4 ####################
##################################################################
library(ggplot2);library(grid);library(lme4)
# These functions were used to test methods for analyzing
# fold change data at the subject level.

# Mixed-effects analysis functions
get_dummy_subject_randomized_vector<-function(times,subjects){
  dummy_times = times
  for(su in unique(subjects)){
    inds = which(subjects==su)
    if(length(inds)<=1){next}
    curr_t = sample(times[inds])
    dummy_times[inds] = curr_t
  }
  return(dummy_times)
}
run_mixed_effect_model1<-function(x,d1,return_pvals=T){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  lmer_obj_0 = lmer(x ~  factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj = lmer(x ~  ordered(timev) + factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj_int1 = list()
  lmer_obj_int1[["tissue"]] =  lmer(x ~  ordered(timev) * factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj_int1[["training"]] =  lmer(x ~  ordered(timev) * factor(training) + factor(tissue) + factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj_int1[["sex"]] =  lmer(x ~  ordered(timev) * factor(sex) + factor(tissue) + factor(training) + (1|dataset/subjs) ,REML=F,data=d)
  a1 = get_pairwise_anova_aic_bic_p(lmer_obj,lmer_obj_0)
  int_a = sapply(lmer_obj_int1,get_pairwise_anova_aic_bic_p,obj2=lmer_obj)
  int_a_0 = sapply(lmer_obj_int1,get_pairwise_anova_aic_bic_p,obj2=lmer_obj_0)
  if(return_pvals){
    pvals = unlist(c(a1[3],int_a[3,],int_a_0[3,]))
    return(pvals)
  }
  return(list(lmer_obj=lmer_obj,lmer_obj_0=lmer_obj_0,interaction_lmers = lmer_obj_int1,pvals=pvals))
}
run_mixed_effect_model2<-function(x,d1,return_pvals=T,empirical_pval=T,...){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  lmer_obj_0 = lmer(x ~  factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj = lmer(x ~  ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  a1 = get_pairwise_anova_aic_bic_p(lmer_obj,lmer_obj_0)
  if(empirical_pval){
    rand_scores = get_lmer_model_empirical_null(
      x ~  ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
      d,RMEL=F)
    emp_p = (1+sum(rand_scores[,2]<=a1[2]))/(1+nrow(rand_scores))
    return(emp_p)
  }
  if(return_pvals){return(as.numeric(a1))}
  return(list(lmer_obj=lmer_obj,lmer_obj_0=lmer_obj_0))
}
get_mixed_effect_model_time_empirical_p<-function(x,frm1,frm0,d1,reps=1000,min_reps=100,statistic="Chisq",thr=0.05,...){
  min_reps = min(min_reps,reps)
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  rand_scores = c()
  m1 = lmer(frm1,d,REML=F)
  m0 = lmer(frm0,d,REML=F)
  anova_real = as.numeric(anova(m1,m0)[2,statistic])
  d_copy=d
  for(j in 1:reps){
    d_copy$timev = get_dummy_subject_randomized_vector(d$timev,d$subjs)
    rand_m1 = lmer(frm1,d_copy,REML=F,control=lmerControl(check.rankX =  "silent.drop.cols"))
    rand_m0 = lmer(frm0,d_copy,REML=F,control=lmerControl(check.rankX =  "silent.drop.cols"))
    anova_rand = as.numeric(anova(rand_m1,rand_m0)[2,statistic])
    rand_scores[j] = anova_rand
    if(j>=min_reps){
      emp_p = (1+sum(rand_scores>=anova_real))/(j+1)
      if(emp_p>thr){break}
      print(emp_p)
    }
  }
  return(emp_p)
}
get_pairwise_anova_aic_bic_p<-function(obj1,obj2){
  return(anova(obj1,obj2)[2,c(2,3,8)])
}
get_mixed_effect_model_time_apprx_p<-function(x,frm1,frm0,d1,statistic="Pr(>Chisq)",...){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  m1 = lmer(frm1,d,REML=F)
  m0 = lmer(frm0,d,REML=F)
  return(as.numeric(anova(m1,m0)[2,statistic]))
}
get_mixed_effect_model_time_apprx_stat_diff<-function(x,frm1,frm0,d1,statistic="Pr(>Chisq)",...){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  m1 = lmer(frm1,d,REML=F)
  m0 = lmer(frm0,d,REML=F)
  an =anova(m1,m0)[,statistic]
  return(an[1]-an[2])
}

get_stats_df_chisq<-function(stats,ps,range,use_log=T,func=abs){
  mse_minus_log = c()
  if(use_log){log_ps = -log(ps,base=10)}
  for(j in range[-1]){
    curr_ps = pchisq(stats,df=j,lower.tail = F)
    curr_ps[curr_ps<min(ps)] = min(ps)
    mse_minus_log[as.character(j)] = mean(func(ps-curr_ps))
    if(use_log){
      curr_ps = -log(curr_ps,base=10)
      mse_minus_log[as.character(j)] = mean(func(log_ps-curr_ps))
    }
    plot(curr_ps,log_ps,main=j);abline(0,1)
  }
  return(range[which(mse_minus_log==min(mse_minus_log))[1]])
}
get_stats_df_chisq2<-function(stats,ps,range){
  comparison_ps = c()
  for(j in range[-1]){
    curr_ps = pchisq(stats,df=j,lower.tail = F)
    curr_ps[curr_ps<min(ps)] = min(ps)
    p = wilcox.test(ps,curr_ps,alternative = "less",paired=T)$p.value
    print(p)
    comparison_ps[j]=p
    if (p < 1e-50){break}
  }
  plot(comparison_ps)
  return(j)
}
get_stats_df_chisq3<-function(stats,ps,range){
  comparison_ps = c()
  for(j in range[-1]){
    curr_ps = pchisq(stats,df=j,lower.tail = F)
    curr_ps[curr_ps<min(ps)] = min(ps)
    if(all(curr_ps>=ps)){break}
  }
  return(j)
}
################ Visualization of subject specific data ################
get_profile_longi_line_plot<-function(x,time,subj,subj_groups=NULL){
  d = data.frame(profile=x,time=time,subj=subj)
  if(is.null(subj_groups)){
    p <- ggplot(data = d, aes(x = time, y = profile, group = subj))
  }
  else{
    p <- ggplot(data = d, aes(x = time, y = profile, group = subj,colour=subj_groups))
  }
  #p + geom_line()
  return(p)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




