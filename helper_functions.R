library(e1071)
# Helper functions

simplify_tissue_info<-function(subject_col2tissue){
  subject_col2tissue_slim = subject_col2tissue
  subject_col2tissue_slim[
    grepl("muscl",subject_col2tissue,ignore.case = T) |
      grepl("vastus",subject_col2tissue,ignore.case = T) |
      grepl("bicep",subject_col2tissue,ignore.case = T)
    ] = "muscle"
  subject_col2tissue_slim[
    grepl("adipo",subject_col2tissue,ignore.case = T) |
      grepl("fat",subject_col2tissue,ignore.case = T)
    ] = "fat"
  subject_col2tissue_slim[
    grepl("blood",subject_col2tissue,ignore.case = T) |
      grepl("cytes",subject_col2tissue,ignore.case = T) |
      grepl("pbmc",subject_col2tissue,ignore.case = T) |
      grepl("phil",subject_col2tissue,ignore.case = T)
    ] = "blood"
  return(subject_col2tissue_slim)
}

simplify_sex_info<-function(sexinfo){
  newv = as.character(sexinfo)
  names(newv) = names(sexinfo)
  newv[grepl(": f",newv)] = "female"
  newv[newv!="" & newv!="female"] = "male"
  newv[newv==""] = NA
  return(newv)
}

simplify_training_type<-function(metadata){
  raw_training_data = as.character(metadata$Training.program.during.experiment.type)
  dataset_subgroup = as.character(metadata$Study.subgroup)
  sample_is_endurance = grepl("endur",raw_training_data,ignore.case = T)
  sample_is_resistance = grepl("resis",raw_training_data,ignore.case = T)| 
    grepl("strength",raw_training_data,ignore.case = T)
  table(sample_is_endurance,sample_is_resistance)
  # The short description of the training program
  sample2training_type = rep("other",nrow(metadata))
  sample2training_type[sample_is_endurance] = "endurance"
  sample2training_type[sample_is_resistance] = "resistance"
  sample2training_type[sample_is_endurance & sample_is_resistance] = "both"
  sample2training_type[grepl("untrained",dataset_subgroup,ignore.case = T)] = "untrained"
  sample2training_type[grepl("control",raw_training_data,ignore.case = T)] = "untrained"
  sample2training_type[grepl("no training",raw_training_data,ignore.case = T)] = "untrained"
  sample2training_type[grepl("yoga",raw_training_data,ignore.case = T)] = "yoga"
  table(sample2training_type)
  names(sample2training_type) = metadata[,1]
  return(sample2training_type)
}
get_subject_info_from_gsms<-function(subj,gsms,metadata,colname,sample2subject){
  subj_gsms = gsms[sample2subject[gsms]==subj]
  currdata = as.character(metadata[subj_gsms,colname])
  return(unique(currdata)[1])
}

# Sept 2017 time windows: selected to have at least 3 cohorts in a tp
# The method below gives:
# blood muscle
# 1     17      5
# 4      4      8
# 24     6      7
simplify_time_acute<-function(tt){
  tt[tt<=2.5] = 1
  tt[tt>2.5 & tt<=5] = 4
  tt[tt>=20] = 24
  return(tt)
}
simplify_time_longterm_blood<-function(tt){
  tt[tt<=200] = 100
  return(tt)
}
simplify_time_longterm_muscle<-function(tt){
  tt[tt<=150] = 100
  tt[tt>150] = 200
  return(tt)
}
simplify_time_longterm<-function(tt,type="muscle"){
  if(type=="muscle"){
    return (simplify_time_longterm_muscle(tt))
  }
  if(type=="blood"){
    return(simplify_time_longterm_blood(tt))
  }
  tt[tt==70.1 | tt==70.2 | tt == 70] = 80
  tt[tt==84] = 80
  tt[tt>100 & tt<200] = 150
  return(tt)
}

# # Before Sept 2017
# simplify_time_acute<-function(tt){
#   tt[tt<=1] = 0.5
#   tt[tt>=2.5 & tt<=5] = 4
#   tt[tt==20] = 24
#   tt[tt==72 | tt==96] = 48
#   return(tt)
# }


####################### Working with GPL tables ####################
merge_two_id_lists<-function(l1,l2){
  lnew = lapply(1:length(l1),function(x,y,z)union(y[[x]],z[[x]]),y=l1,z=l2)
  return(lnew)
}
split_gpl_table_entry <-function(x){strsplit(x,split="( /// )|,|;",perl=T)}
map_gpl_rows_to_entrez<-function(x,converters){
  row2entrez=lapply(1:nrow(x),function(x)c())
  # look for entrez column
  is_entrez = which(grepl(colnames(x),pattern="entrez",ignore.case=T))
  if(length(is_entrez)>0){
    curr_entrez = as.character(x[,is_entrez[1]])
    row2entrez = split_gpl_table_entry(curr_entrez)
  }
  for(conv_fields in names(converters)){
    curr_cols = which(grepl(colnames(x),pattern=conv_fields,ignore.case=T))
    print(conv_fields)
    if(length(curr_cols)>0){
      for (j in curr_cols){
        currv = as.character(x[,j])
        currl = split_gpl_table_entry(currv)
        currl_vals = unique(unlist(currl))
        curr_conv = converters[[conv_fields]][currl_vals]
        curr_conv = curr_conv[!is.na(names(curr_conv))]
        currl = lapply(currl,function(x,y)unique(unlist(y[x])),y=curr_conv)
        row2entrez = merge_two_id_lists (row2entrez,currl)
      }
    }
  }
  return(row2entrez)
}
get_num_probe_genes <- function(x){
  if(is.null(x)){return(0)}
  return(length(na.omit((x))))
}
reverse_mapping_list<-function(l){
  newl = list()
  for(currname in names(l)){
    currvals = l[[currname]]
    for(v in currvals){
      v = as.character(v)
      newl[[v]] = union(newl[[v]],currname)
    }
  }
  return(newl)
}

########### Functions for transforming into genes ##############
# Convert a matrix x with probe ids as names into a new
# matrix of genes. 
# h is the list that maps each gene to its probe ids
# f is the function to apply, default is mean
probe2genes_conv<-function(x,h,f=mean, genes=NULL,...){
  if (is.null(genes)){genes = names(h)}
  y = sapply(genes,applyFunctionOnGroup,h=h,f=f,x=x,...)
  y = t(as.matrix(data.frame(y,check.names=F)))
  return (y)
}
# Apply a function f on values in the vector x (or a matrix)
# that are under the names of h[[g]]
# For example if h maps genes to probes then we can apply
# the mean value of the probes by setting f=mean
# We assume that all mapped probes are necessarily represented as 
# the names in x.
applyFunctionOnGroup<-function(g,h,x,f,intsct=F,...){
  curr_names = h[[g]]
  if (is.null(dim(x))){return (f(x[curr_names]))}
  if (length(curr_names)==1){return (sapply(x[curr_names,],f))}
  return (apply(x[curr_names,],2,f,...))
}

# clean a mapping list by the row names in m
clean_mapping<-function(m,h){
  newh = lapply(h,function(x,y)intersect(x,y),y=rownames(m))
  names(newh) = names(h)
  newh = newh[sapply(newh,length)>0]
  return(newh)
}

# A function that collects data on the NA cells in a matrix
collect_na_stats<-function(x){
  na_x = is.na(x)
  stats = list()
  stats[['percent NAs']] = sum(c(na_x))/length(na_x)
  stats[['num rows with NAs']] = sum(apply(na_x,1,any))
  stats[['num cols with NAs']] = sum(apply(na_x,2,any))
  stats[['row NA counts']] = rowSums(na_x)
  stats[['col NA counts']] = colSums(na_x)
  return(stats)
}
# A function that transforms a probe level data into a gene level matrix
# OPTIONAL: set out_path to save the new gene matrix to an RData file
# ... additional parameters for probe2genes_conv
transform_matrix_into_genes<-function(gsms_mat,gpl,gpl_mappings_entrez2probes,out_path=NULL,...){
  mode(gsms_mat) = 'numeric'
  gsms_mat_na_stats = collect_na_stats(gsms_mat)
  curr_mapping = gpl_mappings_entrez2probes[[gpl]]
  curr_mapping = clean_mapping(gsms_mat,curr_mapping)
  print('difference before and after cleaning the gpl mapping list:')
  print(paste(length(gpl_mappings_entrez2probes[[gpl]]),length(curr_mapping)))
  print('merging genes by averaging probes')
  entrez_mat = probe2genes_conv(gsms_mat,curr_mapping,na.rm=T)
  print('done')
  print(dim(entrez_mat))
  entrez_mat_na_stats = collect_na_stats(entrez_mat)
  return(list(gsms_mat_na_stats=gsms_mat_na_stats,entrez_mat=entrez_mat,entrez_mat_na_stats=entrez_mat_na_stats))
}

########################## Working with GEOquery data ##########################
get_gsm_data_size<-function(gsm){
  if(class(gsm)!="GSM"){return(0)}
  return(nrow(Table(gsm)))
}

########################## Extract data from our data contrainers ##########################
get_submatrix<-function(m,cols){
  cols = intersect(colnames(m),cols)
  if(length(cols)==0){return (NULL)}
  return (m[,cols])
}
# assumption: all gsms are of the same platform
get_data_matrix_from_matrix_lists<-function(m_list,gsms){
  l = lapply(m_list,get_submatrix,cols=gsms)
  l = l[!sapply(l,is.null)]
  l = l[sapply(l,length)>0]
  if(length(l)==0){return(NULL)}
  newm = c()
  for(mm in l){newm = cbind(newm,mm)}
  if(length(newm)<2 || is.null(newm)){return (NULL)}
  if(is.null(dim(newm))){newm = matrix(newm,ncol=1)}
  gsms = intersect(gsms,colnames(newm))
  if(length(gsms)==0){return(NULL)}
  return(newm[,gsms])
}

######################### Analysis and statistical tests ##########################
get_time_from_subj_names<-function(s){
  arr = strsplit(s,split="_")[[1]]
  return(as.numeric(arr[length(arr)]))
}

# x,y,d all have the same sample set in the same order
run_lm<-function(y,x,...){
  obj = lm(y~.,data=data.frame(y,x,check.names=F))
  return(obj)
}
leave_study_out <-function(y,x,d,func=run_lm,...){
  preds_list = list()
  for(batch in unique(d)){
    print(batch)
    inds = which(d==batch)
    trx = x[-inds,];try = y[-inds]
    m = func(try,trx)
    p = predict(m,data.frame(x[inds,],check.names=F))
    preds_list[[batch]] = cbind(p,y[inds])
    colnames(preds_list[[batch]]) = c("predictions","real.labels")
    if(!all(names(p)==rownames(x)[inds])){
      print ("ERROR in prediction")
      break
    }
  }
  return(preds_list)
}
library(randomForest)
run_rf<-function(y,x,...){randomForest(x,y,...)}
# Assumes that all input other than x are factors
run_anova_on_time_and_subject<-function(x,subjs,timev){
  aov_obj = aov(x~ timev + Error(subjs/timev))
  pval = summary(aov_obj)[[2]][[1]][1,5]
  return(list(aov_obj=aov_obj,pval=pval))
}

# Assumes that all input other than x are factors or numbers
run_anova_on_time_age_sex_and_subject<-function(x,subjs,timev,sexv=NULL,agev=NULL,get_pvalue=T){
  has_sex = !is.null(sexv) && !any(sexv=="") && length(levels(sexv))>1
  has_age = !is.null(agev) && !any(agev=="") && (length(levels(agev))>1 || class(agev)=="numeric")
  if(has_sex && has_age){
    aov_obj = aov(x~timev*sexv*agev + Error(subjs))
    pval = get_timev_pval(summary(aov_obj)[[2]][[1]])
    return(list(aov_obj=aov_obj,pval=pval))
  }
  if(has_sex){
    aov_obj = aov(x~timev*sexv + Error(subjs))
    pval = get_timev_pval(summary(aov_obj)[[2]][[1]])
    return(list(aov_obj=aov_obj,pval=pval))
  }
  if(has_age){
    aov_obj = aov(x~timev*agev + Error(subjs))
    pval = get_timev_pval(summary(aov_obj)[[2]][[1]])
    return(list(aov_obj=aov_obj,pval=pval))
  }
  return(run_anova_on_time_and_subject(x,subjs,timev))
}
get_timev_pval<-function(results_table){
  rownames(results_table) = gsub(rownames(results_table),pattern=" ",replace="")
  return(results_table["timev",5])
}
# x= data_matrix[1,curr_samples];subjs = factor(curr_subjects);timev=factor(curr_times,ordered=T);sexv=factor(curr_sex)
# ttt = run_anova_on_time_age_sex_and_subject(data_matrix[1,curr_samples],curr_subjects,curr_times,curr_sex,NULL,F)
# ttt_sum = summary(ttt)
# get_anova_p(ttt)

run_paired_test <- function(x,subjs,timev,times=NULL,func=t.test,...){
  if(is.null(times)){
    times = sort(unique(as.numeric(timev)))[1:2]
  }
  names(x) = subjs
  x1 = x[timev==times[1]];x2=x[timev==times[2]]
  x2 = x2[names(x1)]
  #print(all(names(x1)==names(x2)))
  return(func(x1,x2,paired=T,...))
}
# install.packages('lme4')
library(lme4)
run_lmer4_anova<-function(x,subjs,timev,sexv=NULL,agev=NULL,get_pvalue=T){
  has_sex = !is.null(sexv) && !any(sexv=="") && length(levels(sexv))>1
  has_age = !is.null(agev) && !any(agev=="") && (length(levels(agev))>1 || class(agev)=="numeric")
  if(has_sex && has_age){
    lm2 = lmer(x ~timev + sexv + agev + (1|subjs),REML=F)
    lm1 = lmer(x ~ 1 + sexv + agev + (1|subjs),REML=F)
  }
  else if(has_sex){
    lm2 = lmer(x ~ timev + sexv + (1|subjs),REML=F)
    lm1 = lmer(x ~ 1 + sexv + (1|subjs),REML=F)
  }
  else if(has_age){
    lm2 = lmer(x ~ timev + agev + (1|subjs),REML=F)
    lm1 = lmer(x ~ 1 + agev + (1|subjs),REML=F)
  }
  else{
    lm2 = lmer(x ~ timev + (1|subjs) ,REML=F)
    lm1 = lmer(x ~ 1 + (1|subjs),REML=F)
  }
  anova_obj = anova(lm1,lm2,refit=F)
  pval = anova_obj[2,8]
  if(get_pvalue){return(pval)}
  return(list(lm1=lm1,lm2=lm2,anova_obj=anova_obj,pval=pval))
}

# #####################################################################
# # QA of the functions above
# nsubj = 10
# nt = 2
# mean_diff = 4
# subj = factor(rep(1:nsubj,nt))
# timev = c(sapply(1:nt,function(x)rep(x,nsubj)))
# subj2sex = c(rep("F",nsubj/2),rep("M",nsubj/2)) # Balanced sexes
# sexv = factor(subj2sex[subj])
# # Effect due to time only
# x = rnorm(nsubj*nt)
# for(j in 2:nt){x[timev==j] = x[timev==j]+mean_diff}
# p1 = run_paired_ttest(x,subjs,timev)$p.value
# p2 = run_anova_on_time_and_subject(x,subjs,timev)$pval
# p3 = run_anova_on_time_age_sex_and_subject(x,subjs,timev)$pval
# print(abs(p1-p2))
# print(abs(p2-p3))
# # Effect due to sex only
# x = rnorm(nsubj*nt)
# for(j in unique(sexv)[-1]){x[sexv==j] = x[sexv==j]+mean_diff}
# p1 = run_paired_ttest(x,subjs,timev)$p.value
# p2 = run_anova_on_time_and_subject(x,subjs,timev)$pval
# p3 = run_anova_on_time_age_sex_and_subject(x,subjs,timev,sexv=sexv)$pval
# p4 = run_lmer4_anova(x,subjs,timev,sexv = sexv)
# summary(p3[[1]])
# print(abs(p1-p2))
# print(abs(p2-p3))

#####################################################################
# We currently take the first repeat when subjects have more than
# a single time series data
get_fold_changes_vs_baseline<-function(x,subjs,timev,baseline=NULL,func = function(a,b){a-b},metadata=NULL){
  if(is.null(baseline)){baseline = sort(unique(timev))[1]}
  print(baseline)
  newx = c()
  for(subj in unique(subjs)){
    inds = subjs==subj
    if(sum(inds)<=1){next}
    inds_t = timev[inds]
    currx = x[,inds]
    #print(subj);print(timev[inds])
    base_t_ind = which(inds_t==baseline)
    if(ncol(currx)>2){
      currx_diff = apply(currx[,-base_t_ind],2,func,b=currx[,base_t_ind])
    }
    else{
      currx_diff = as.matrix(func(currx[,-base_t_ind],currx[,base_t_ind]),ncol=1)
    }
    colnames(currx_diff) = paste(subj,inds_t[-base_t_ind],sep="_")
    newx = cbind(newx,currx_diff)
  }
  rownames(newx) = rownames(x)
  return(newx)
}

################ GO enrichment analysis ################
# GO enrichment
library(topGO)
run_topgo_enrichment_fisher<-function(genesOfInterest,geneUniverse,go_term_size=10,go_max_size=1000,go_dags=c("BP","MF"),...){
  l = list()
  if(class(genesOfInterest)=="character"){
    l[["set1"]] = genesOfInterest
  }
  else{
    l = genesOfInterest
  }
  res = c()
  for(nn in names(l)){
    genesOfInterest = l[[nn]]
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    for(type in go_dags){
      myGOdata <- new("topGOdata", description="My project", ontology=c("BP"), nodeSize=go_term_size,
                  allGenes=geneList, annot = annFUN.org, mapping="org.Hs.eg.db", ID = "entrez")
      allGOs = usedGO(myGOdata)
      resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")
      allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = 
                       "resultFisher", ranksOf = "classicFisher", topNodes = length(score(resultFisher)))
      setname = rep(nn,nrow(allRes))
      typev = rep(type,nrow(allRes))
      allRes = cbind(setname,typev,allRes)
      res = rbind(res,allRes)
    }
  }
  go_pvals = as.numeric(res[,ncol(allRes)])
  go_pvals[is.na(go_pvals)] = 1e-30
  go_qvals = p.adjust(go_pvals,method='fdr')
  res = cbind(res,go_qvals)
  return(res)
}
extract_top_go_results<-function(res,qval=0.1,maxsize=2000){
  res = res[,-which(colnames(res)=="typev")]
  res = res[res$Annotated<maxsize,]
  res = unique(res)
  res[res$GO.ID=="GO:0002532",]
  new_qvals = p.adjust(res$classicFisher,method='fdr')
  #plot(new_qvals,res$go_qvals);abline(0,1)
  res$go_qvals = new_qvals
  return(res[res$go_qvals<=qval,])
}
get_most_sig_enrichments_by_groups <- function(res){
  gs = unique(as.character(res[,1]))
  m = c()
  for(g in gs){
    res0 = res[res[,1]==g,]
    ps = as.numeric(res0$classicFisher)
    m = rbind(m,res0[ps==min(ps),])
  }
  return(m)
}

# GSEA
library(fgsea)
fgsea_wrapper <- function(pathways,scores,nperm=2000,run_nperm=1000,...){
  num_runs = nperm/run_nperm
  l = list()
  for(j in 1:num_runs){
    l[[j]] = fgsea(pathways,scores,nperm = run_nperm,...)
  }
  emp_pvals = sapply(l,function(x)x$pval)
  emp_pvals = emp_pvals*run_nperm
  min_to_add = min(emp_pvals)
  emp_pvals = emp_pvals-min_to_add
  new_pvals = rowSums(emp_pvals)+min_to_add
  new_pvals = new_pvals/nperm
  new_qvals = p.adjust(new_pvals,method='fdr')
  res = l[[1]]
  res[,"pval"] = new_pvals
  res[,"padj"] = new_qvals
  return(res)
}

################ Visualization ################
require(ggplot2)
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

# library(pbkrtest)
# KRSumFun <- function(object, objectDrop, ...) {
#   krnames <- c("ndf","ddf","Fstat","p.value","F.scaling")
#   r <- if (missing(objectDrop)) {
#     setNames(rep(NA,length(krnames)),krnames)
#   } else {
#     krtest <- KRmodcomp(object,objectDrop)
#     unlist(krtest$stats[krnames])
#   }
#   attr(r,"method") <- c("Kenward-Roger via pbkrtest package")
#   r
# }

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
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

#### Classification analysis methods
################## Implementation of wrappers for base classifiers ####################################
# These are used within some multilabel algorithms, most notably the naive binary relevance method.
# We implement below wrappers for binary and non binary data. These methods
# First run feature selection, record the selected clustering and use them to learn the
# classifier. 
# The non-binary data analysis is based on the CMA package in R.
# For binary data we used Fisher's exact test.
# The main reason for these wrappers: the need to save running time.
#source("https://bioconductor.org/biocLite.R")
#biocLite("CMA")
library(CMA);library(randomForest)
featureSelectionClassifier<-function(x,y,classification_function=randomForest,
        features_inds_to_keep=c(),fsmethod="limma",numFeatures=100,...){
  selection = GeneSelection(x,y,method=fsmethod,trace=F)
  selected_inds = toplist(selection,k=numFeatures,show=F)$index
  selected_inds = union(selected_inds,features_inds_to_keep)
  #print ("features selected")
  classifier = classification_function(x[,selected_inds],y,...)
  #print (paste("classifier learned",class(classifier)))
  obj = list(inds=selected_inds,classifier=classifier)
  class(obj) = "featureSelectionClassifier"
  return (obj)
}

predict.featureSelectionClassifier<-function(obj,newx,...){
  print ("in fs predict")
  args = c(list(obj$classifier,newx[,obj$inds]),...)
  return (do.call(predict,args=args))
}

get_hg_pval<-function(x,y,...){
  tab = table(x,y)
  if (length(tab) !=4){return (1)}
  return (fisher.test(tab,...)$p.value)
}

featureSelectionClassifier_binary_data<-function(x,y,
        classification_function=svm,numFeatures=250,features_to_keep=c(),...){
  f_pvals = apply(x,2,get_hg_pval,y=y,alternative = 'greater')
  names(f_pvals) = colnames(x)
  selected_features1 = names(sort(f_pvals)[1:numFeatures])
  f_pvals = apply(x,2,get_hg_pval,y=y,alternative = 'less')
  names(f_pvals) = colnames(x)
  selected_features2 = names(sort(f_pvals)[1:numFeatures])
  selected_features = union(selected_features1,selected_features2)
  #print (paste(length(selected_features)," features selected"))
  classifier = classification_function(as.matrix(x[,selected_features]),y,...)
  #print (paste("classifier learned",class(classifier)))
  obj = list(inds=selected_features,classifier=classifier)
  class(obj) = "featureSelectionClassifier"
  return (obj)
}

calcAupr <- function(pred, gs,roc=F,useAbs=T) {
  ord.idx = NULL
  if (useAbs){
    ord.idx <- order(abs(pred), decreasing = T)
  }
  else{
    ord.idx <- order(pred, decreasing = T)
  }
  prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx))) #also known as positive predictive value
  rec  <- cumsum(gs[ord.idx]) / sum(gs)                     #also know as true positive rate
  fpr  <- cumsum(gs[ord.idx] == 0) / (length(gs) - sum(gs)) #false positive rate
  prec <- c(prec[1], prec)
  rec <- c(0, rec)
  fpr <- c(0, fpr)
  if (!roc){
    aupr <- areaUnderCurve(rec, prec)
    return (aupr)
  }
  #plot(fpr,rec,type='l')
  auroc <- areaUnderCurve(fpr, rec)
  return(auroc)
}

areaUnderCurve <- function(x, y) {
  dx <- diff(x)
  my <- y[1:(length(y) - 1)] + diff(y) / 2
  return(sum(dx * my))
}

leave_study_out2 <-function(y,x,d,func=run_lm,pred_args=list(),...){
  preds_list = list()
  for(batch in unique(d)){
    #print(batch)
    inds = which(d==batch)
    trx = x[-inds,];try = y[-inds]
    m = func(y=try,x=trx,...)
    p = NULL
    try({
      #print(max(inds));print(dim(x))
      args = c(list(m,x[inds,]),pred_args)
      p = do.call(predict,args)
    })
    if(is.null(p)){
      args = c(list(m,data.frame(x[inds,],check.names=F)),pred_args)
      p = do.call(predict,args)
    }
    preds_list[[batch]] = list(preds=p,real=y[inds])
  }
  return(preds_list)
}

simple_sampling_based_learner<-function(x,y,d,func = featureSelectionClassifier,reps=10,max_num_samples = 1000,...){
  y_table = table(y)
  min_class = names(y_table[y_table==min(y_table)])
  min_class_size = y_table[y_table==min(y_table)]
  positives = rownames(x)[y==min_class]
  datasets = unique(d[positives])
  negatives = rownames(x)[y!=min_class & is.element(d[rownames(x)],set=datasets)]
  bgcs = setdiff(rownames(x),positives);bgcs = setdiff(bgcs,negatives)
  classifiers = list()
  samp_size = min(c(length(positives),max_num_samples))
  for (i in 1:reps){
    newbgcs = NULL
    if(length(bgcs)>0){
      newbgcs = sample(bgcs)[1:min(samp_size,length(bgcs))]
    }
    pos_sample = sample(positives)[1:samp_size]
    neg_sample = NULL
    if(length(negatives)>0){
      neg_sample = sample(negatives)[1:min(samp_size,length(negatives))]
    }
    #print(paste("number of positives, negatives, and bgcs:"))
    #print(c(length(pos_sample),length(neg_sample),length(newbgcs)))
    curr_samples = c(pos_sample,neg_sample,newbgcs)
    curr_inds = is.element(rownames(x),curr_samples)
    newy = y[curr_inds]
    newx = x[curr_inds,]
    classifiers[[i]] = func(newx,as.factor(newy),...)
  }
  obj = classifiers
  class(obj) = "simple_sampling_based_learner"
  return(obj)
}
predict.simple_sampling_based_learner<-function(obj,x,...){
  preds = c()
  counts = 0
  for (i in 1:length(obj)){
    if (length(preds)==0){
      try({
        preds = predict(obj[[i]],x,...);
        counts = counts + 1
      })
    }
    else{
      try({
        preds = preds + predict(obj[[i]],x,...);
        counts = counts + 1
      })
    }
  }
  preds = preds / counts
  return (preds)
}

# Supervised analysis
# A configuration specifies the classification test.
# Each config has the following information in a list:
# The label name
# The features matrix name
# The sample set
# Logical: add covariates or not
# Classifier
# Classifier Arguments
svm_config = list(pred_args = list(probability=T),numFeatures=500,
                  classification_function=svm,probability=T,
                  kernel="linear")
rf_config = list(pred_args = list(type="vote"),numFeatures=500,
                 classification_function=randomForest)
##################################################################
##################################################################
##################################################################
############# Functions for analysis of configurations ###########
##################################################################

run_lso_tests_on_a_configuration_set<-function(configurations,classification_xs,classification_ys,d,x_cov){
  lso_results = list()
  for(cname in names(configurations)){
    config = configurations[[cname]]
    y_names = config$sample_set
    y = classification_ys[[config$yname]][y_names]
    x = classification_xs[[config$xname]][y_names,]
    # In case one of the classes has less than 2 datasets exclude it:
    tab = table(y,d[y_names])
    print(tab)
    classes_to_exclude = rownames(tab)[rowSums(tab>0)<2]
    if(length(classes_to_exclude)>0){
      print(paste("excluded classes in LSO: ",classes_to_exclude,collapse= " "))
      y_names = y_names[!is.element(y,set=classes_to_exclude)]
      y = classification_ys[[config$yname]][y_names]
      y = as.factor(as.character(y))
      names(y) = y_names
      x = classification_xs[[config$xname]][y_names,]
    }
    features_inds_to_keep=c()
    if(config$include_covs=="TRUE"){
      curr_covs = x_cov[,setdiff(colnames(x_cov),config$yname)]
      curr_covs = model.matrix(~.,curr_covs)[,-1]
      new_y_names = intersect(rownames(x),rownames(curr_covs))
      x = x[new_y_names,];y=y[new_y_names]
      x = cbind(curr_covs[new_y_names,],x)
      y_names = new_y_names
      features_inds_to_keep = 1:ncol(curr_covs)
    }
    curr_lso_test_args = c(list(y,x,d[y_names],func=featureSelectionClassifier,
                                features_inds_to_keep=features_inds_to_keep),config$classification_args)
    lso = do.call(leave_study_out2,curr_lso_test_args)
    lso_results[[cname]] = lso
  }
  return(lso_results)
}

get_standard_classification_performance_scores_for_results_list<-function(lso_results){
  performance_scores = list()
  for(cname in names(lso_results)){
    lso = lso_results[[cname]]
    lso_preds = lapply(lso,function(x)x$preds)
    if(! "matrix" %in% class(lso_preds[[1]])){
      lso_preds = lapply(lso,function(x)attr(x$preds,"probabilities"))
    }
    class_names = colnames(lso_preds[[1]])
    m = c()
    for(mm in lso_preds){m = rbind(m,mm[,class_names])}
    colnames(m) = class_names
    lso_preds = m
    aucs = c()
    for(class_name in class_names){
      lso_y = unlist(sapply(lso,function(x)x$real))
      class_inds = lso_y==class_name
      lso_y = rep(0,length(lso_y))
      lso_y[class_inds]=1
      lso_p = lso_preds[,class_name]
      roc = calcAupr(lso_p,as.numeric(lso_y),roc=T)
      proc = as.numeric(pROC::auc(lso_y, lso_p))
      aupr = calcAupr(lso_p,as.numeric(lso_y),roc=F)
      aucs = rbind(aucs,c(roc=roc,aupr=aupr,pROC=proc))
    }
    rownames(aucs) = class_names
    acc = lso_y == (lso_p>=0.5)
    acc = sum(acc)/length(acc)
    performance_scores[[cname]] = list(aucs=aucs,acc=acc)
  }
  return(performance_scores)
}

get_subject_performance_scores_from_results_list<-function(lso_results,name2subj=sample2subj,paired_tests_analysis=T,x_cov){
  subject_performance_scores = list()
  for(cname in names(lso_results)){
    subject_performance_scores[[cname]]=list()
    lso = lso_results[[cname]]
    lso_preds = lapply(lso,function(x)x$preds)
    if(! "matrix" %in% class(lso_preds[[1]])){
      lso_preds = lapply(lso,function(x)attr(x$preds,"probabilities"))
    }
    class_names = colnames(lso_preds[[1]])
    m = c()
    for(mm in lso_preds){m = rbind(m,mm[,class_names])}
    colnames(m) = class_names
    lso_preds = m
    curr_subjects = name2subj[rownames(lso_preds)]
    lso_y = unlist(sapply(lso,function(x)x$real))
    class1 = colnames(lso_preds)[1]
    if(paired_tests_analysis){
      get_profile_longi_line_plot(lso_preds[,class1],lso_y,curr_subjects) + 
        geom_line()
      xx1 = c();xx2=c()
      for(subj in unique(curr_subjects)){
        subj_inds = curr_subjects==subj
        if(sum(subj_inds)!=2){next} 
        curr_y = lso_y[subj_inds]
        if(length(unique(curr_y))!=2){next}
        curr_p = lso_preds[subj_inds,class1]
        xx1 = c(xx1,curr_p[curr_y==class1])
        xx2 = c(xx2,curr_p[curr_y!=class1])
      }
      curr_p = NA
      if(length(xx1)>0){curr_p = wilcox.test(xx1,xx2,paired=T)$p.value}
      d1 = data.frame(subjs=factor(curr_subjects),timev=lso_y,
                      tissue = x_cov[rownames(lso_preds),"tissue"],
                      training=x_cov[rownames(lso_preds),"training"],
                      sex=x_cov[rownames(lso_preds),"sex"])
      mixed_effects_model_emp_pval = NA
      # Should not work when splitting the database by one of the covariates
      # TODO: fixed later
      try({ 
        mixed_effects_model_emp_pval = get_mixed_effect_model_time_empirical_p(x=lso_preds[,class1],
                                                                               frm0 = x ~ factor(tissue) + factor(training) + factor(sex) + (1|subjs),
                                                                               frm1 = x ~ ordered(timev) + factor(tissue) + factor(training) + factor(sex) + (1|subjs),
                                                                               d1=d1,reps=100,min_reps = 100)
      })
      subject_performance_scores[[cname]][["paired_tests_analysis"]] = 
        list(x1=xx1,x2=xx2,binary_wilcox_p=curr_p,mixed_effects_emp_p = mixed_effects_model_emp_pval)
    }
    subj2acc = c()
    for(subj in curr_subjects){
      subj_inds = curr_subjects==subj
      if(sum(subj_inds)<2){next}
      curr_y = lso_y[subj_inds]
      if(!grepl("time",cname) || length(unique(curr_y))==1){
        curr_y = lso_y[subj_inds] == class1
        curr_p = mean(lso_preds[subj_inds,class1])>0.5
        subj2acc[subj] = sum(curr_p==curr_y[1])
      }else{
        curr_p = as.numeric(lso_preds[subj_inds,class1]>0.5)
        subj2acc[subj] = mean(curr_p==curr_y)
      }
    }
    subject_performance_scores[[cname]][["acc"]] = mean(subj2acc,na.rm=T)
    print(paste(cname,mean(subj2acc,na.rm=T)))
  }
  return(subject_performance_scores)
}
##################################################################
##################################################################
##################################################################

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

# # Functions for replicability analysis
# # Get the data matrices for the analysis
# get_paired_ttest_matrices = function(dataset_ids,tt,sample2subj_with_dataset_number,
#                                      dataset2preprocessed_data,geneset,min_subjs = 3,toprint=T,toplot=T){
#   dataset_pvals = c();dataset_stats = c();dataset_pvals_two_tail=c()
#   for (dataset in unique(dataset_ids)){
#     all_dataset_samples = names(dataset_ids[dataset_ids==dataset])
#     time_points = setdiff(unique(tt[all_dataset_samples],na.rm=T),min(tt))
#     if(length(time_points)==0){next}
#     #if(toprint){print(time_points)}
#     for(tt_val in time_points){
#       dataset_samples = all_dataset_samples
#       curr_tt = tt[dataset_samples]
#       dataset_samples = dataset_samples[curr_tt==min(tt) | curr_tt==tt_val]
#       dataset_subjects = sample2subj_with_dataset_number[dataset_samples]
#       # all NAs == we do not have the data
#       if(all(is.na(dataset_subjects))){next}
#       tab = table(dataset_subjects)
#       paired_subjects = names(tab[tab==2])
#       if(length(paired_subjects)<min_subjs){next}
#       #print(table(curr_tt))
#       # select subjects with paired samples
#       dataset_samples = dataset_samples[is.element(dataset_subjects,set=paired_subjects)]
#       dataset_subjects = sample2subj_with_dataset_number[dataset_samples]
#       dataset_times = tt[dataset_samples]
#       # reorder  
#       ord = order(dataset_times,dataset_subjects,decreasing=F)
#       dataset_samples = dataset_samples[ord]
#       dataset_times = dataset_times[ord]
#       dataset_subjects = dataset_subjects[ord]
#       assertthat::assert_that(all(dataset_subjects[dataset_times==min(tt)]==dataset_subjects[dataset_times>min(tt)]))
#       # run t-test
#       mat = dataset2preprocessed_data[[dataset]]$gene_data[geneset,dataset_samples]
#       if(is.null(mat)){next}
#       ttests = apply(mat,1,run_paired_test,subjs = dataset_subjects,timev=dataset_times,alternative="less")
#       ttest_stats = sapply(ttests,function(x)x$statistic)
#       ttest_pvals = sapply(ttests,function(x)x$p.value)
#       # save the results
#       formatted_name = get_dataset_name_for_rep_analysis(dataset,dataset_ids,sample2tissue,sample2training_type)
#       tt_val = as.character(tt_val)
#       formatted_name = paste(tt_val,formatted_name,sep=';')
#       print(formatted_name)
#       if(toplot){hist(ttest_pvals,main=dataset)}
#       dataset_pvals = cbind(dataset_pvals,ttest_pvals)
#       colnames(dataset_pvals)[ncol(dataset_pvals)] = formatted_name
#       dataset_stats = cbind(dataset_stats,ttest_stats)
#       colnames(dataset_stats)[ncol(dataset_stats)] = formatted_name
#       dataset_pvals_two_tail = cbind(dataset_pvals_two_tail,2*pmin(1-ttest_pvals,ttest_pvals))
#       colnames(dataset_pvals_two_tail)[ncol(dataset_pvals_two_tail)] = formatted_name
#     }
#   }
#   rownames(dataset_stats) = gsub(rownames(dataset_stats),pattern = "\\.t$",replace="",perl=T)
#   return(list(dataset_stats=dataset_stats,dataset_pvals=dataset_pvals,dataset_pvals_two_tail=dataset_pvals_two_tail))
# }
# 
# get_dataset_name_for_rep_analysis<-function(nn,dataset_ids,sample2tissue,sample2training_type){
#   samp = names(which(dataset_ids==nn))[1]
#   gse = strsplit(split=';',nn)[[1]][1]
#   newname = paste(sample2tissue[samp],sample2training_type[samp],gse,sep=';')
#   return(newname)
# }


## Simple functions for cohort-level analyses
# some simple functions to analyze a time point in a dataset's matrix
get_paired_ttest_yi_vi <-function(x,sample2time,t1,t2){
  x1 = x[sample2time==t1]
  x2 = x[sample2time==t2]
  d = x2-x1
  n = length(d)
  sdd = sd(d)/sqrt(length(d))
  return(c(yi=mean(d),vi=sdd^2))
}
# # test
# x1 = rnorm(10);x2=rnorm(10);x=c(x1,x2)
# s2t = c(rep(0,10),rep(1,10))
# pp =get_paired_ttest_yi_vi(x,s2t,0,1)
# pp
# tt = t.test(x1,x2,paired = T)
# tt$estimate/tt$statistic
# tt$estimate

get_ttest_yi_vi_per_dataset<-function(mat,metadata){
  dataset_times = metadata$time[colnames(mat)]
  if(any(is.na(dataset_times))){
    mat = mat[,!is.na(dataset_times)]
    dataset_times = metadata$time[colnames(mat)]
  }
  min_time = min(dataset_times)
  other_times = setdiff(unique(dataset_times),min_time)
  times2effects = list()
  for(other_time in other_times){
    curr_mat = mat[,dataset_times==min_time | dataset_times==other_time]
    curr_subjects = metadata$subject[colnames(curr_mat)]
    subjects_to_keep = names(which(table(curr_subjects)==2))
    curr_mat = curr_mat[,is.element(curr_subjects,set = subjects_to_keep)]
    ord = order(metadata$subject[colnames(curr_mat)],metadata$time[colnames(curr_mat)])
    curr_mat = curr_mat[,ord]
    curr_times = metadata$time[colnames(curr_mat)]
    paired_test_data = t(apply(curr_mat,1,get_paired_ttest_yi_vi,sample2time=curr_times,t1=min_time,t2=other_time))
    times2effects[[as.character(other_time)]] = paired_test_data
  }
  return(times2effects)
}
get_ttest_pval_per_dataset<-function(mat,metadata){
  dataset_times = metadata$time[colnames(mat)]
  if(any(is.na(dataset_times))){
    mat = mat[,!is.na(dataset_times)]
    dataset_times = metadata$time[colnames(mat)]
  }
  min_time = min(dataset_times)
  other_times = setdiff(unique(dataset_times),min_time)
  times2pvals = c()
  for(other_time in other_times){
    curr_mat = mat[,dataset_times==min_time | dataset_times==other_time]
    curr_subjects = metadata$subject[colnames(curr_mat)]
    subjects_to_keep = names(which(table(curr_subjects)==2))
    curr_mat = curr_mat[,is.element(curr_subjects,set = subjects_to_keep)]
    ord = order(metadata$subject[colnames(curr_mat)],metadata$time[colnames(curr_mat)])
    curr_mat = curr_mat[,ord]
    curr_times = metadata$time[colnames(curr_mat)]
    paired_test_data = apply(curr_mat,1,get_paired_ttest_pval,sample2time=curr_times,t1=min_time,t2=other_time)
    times2pvals = cbind(times2pvals,paired_test_data)
    colnames(times2pvals)[ncol(times2pvals)] = as.character(other_time)
  }
  return(times2pvals)
}
# some simple functions to analyze a time point in a dataset's matrix
get_paired_ttest_pval <-function(x,sample2time,t1,t2){
  x1 = x[sample2time==t1]
  x2 = x[sample2time==t2]
  return(t.test(x1,x2,paired = T)$p.value)
}
##
