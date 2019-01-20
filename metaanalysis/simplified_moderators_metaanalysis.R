# Our algorithm for analysis of a single gene
# Input: datasets for time point t and gene g
# 1. Run meta-analysis for endurence and resistance (and both)
# 2. Run meta-analysis for controls
# 3. Get controls intercept estimation beta_c
# 4. Get current gene's overall diff expression: beta_e = max(beta_0,beta_0 + beta_endurence)
# 5. SCORE 1: abs(beta_e-beta_c)
# 6. SCORE 2: Egger's test of the exercise meta-analysis
# 7. Return SCORE 1, SCORE 2, and the betas (and their significance) from the exercise meta-analysis

setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(org.Hs.eg.db);library(metafor)
source('/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)

############################################################################
############################################################################
############################################################################
# Functions for data filtering

#' Change the time field from numeric to ordered factor.
#' @param gdata A data frame with a column "time"
#' @return A new data frame.
simplify_time_in_gdata<-function(gdata,func=binarize_feature,...){
  gdata$time = func(gdata$time,...)
  gdata$time = ordered(gdata$time)
  return(gdata)
}
binarize_feature<-function(tt,thr=5){return(as.numeric(tt <= thr))}
simplify_age_gdata<-function(gdata,func=binarize_feature,...){
  gdata$avg_age = func(as.numeric(gdata$avg_age),...)
  gdata$avg_age = ordered(gdata$avg_age)
  return(gdata)
}
simplify_training_for_exercise_analysis<-function(gdata){
  rows = !is.element(gdata$training,set=c("yoga","control","untrained"))
  gdata = gdata[rows,]
  gdata$training[gdata$training=="both"] = "resistance"
  return(gdata)
}

# Clean the datasets: from muscle, remove 1 cohort with time point zero
# From long term exclude all studies with time > 150 days
clean_acute_table <-function(gdata,tissue="muscle"){
  gdata = gdata[gdata$tissue==tissue,]
  gdata$vi = pmax(gdata$vi,1e-5)
  if(tissue=="muscle"){
    gdata = gdata[as.numeric(gdata$time) > 1,]
  }
  gdata = simplify_time_in_gdata(gdata,simplify_time_acute)
  gdata = simplify_age_gdata(gdata,thr=40)
  gdata = simplify_training_for_exercise_analysis(gdata)
  return(gdata)
}
clean_longterm_table <-function(gdata,tissue="muscle"){
  gdata = gdata[gdata$tissue==tissue,]
  gdata$vi = pmax(gdata$vi,1e-5)
  gdata = gdata[as.numeric(gdata$time) < 150 ,]
  gdata = simplify_age_gdata(gdata,thr=40)
  gdata = simplify_training_for_exercise_analysis(gdata)
  return(gdata)
}
get_untrained_table <-function(gdata,tissue="muscle"){
  gdata = gdata[gdata$tissue==tissue,]
  gdata$vi = pmax(gdata$vi,1e-5)
  rows = is.element(gdata$training,set=c("yoga","control","untrained"))
  return(gdata[rows,])
}

meta_analysis_wrapper<-function(gdata,func = rma.mv,...){
  for(rel.tol in c(1e-8,1e-7,1e-6)){
    cc=list(iter.max=10000,rel.tol=rel.tol)
    res = NULL
    try({
      res = func(yi,vi,data=gdata,...)
    },silent=TRUE)
    if(!is.null(res)){return(res)}
  }
  return(NA)
}
run_simple_re<-function(gdata){return(meta_analysis_wrapper(gdata,func=rma.uni))}
try_get_field<-function(res,fname){
  if(is.element(fname,set=names(res))){
    return(res[[fname]])
  }
  return(NA)
}
pvalue_qqplot<-function(ps,n_random=50000,...){
  x = runif(n_random)
  qqplot(y=-log(ps,10),x=-log(x,10),xlab="Theoretical",ylab="Observed",...)
  abline(0,1,lty=2,lwd=2)
}
simple_stouffer_meta_analysis<-function(gdata){
  gdata = gdata[!is.na(as.numeric(gdata$p)),]
  p = metap::sumz(as.numeric(gdata$p),weights = as.numeric(gdata$df))$p[1,1]
  return(p)
}

get_coeffs_str<-function(coeffs){
  coeffs_v = NULL
  for(i in 1:nrow(coeffs)){
    if(i==1){
      coeffs_v = paste(coeffs_v,paste("b0:",
        format(coeffs[i,"beta"],digits=2),"(",
        format(coeffs[i,"pval"],digits=2),")"
        ,sep=""),sep="")
    }
    else{
      coeffs_v = paste(coeffs_v,paste(rownames(coeffs)[i],":",
        format(coeffs[i,"beta"],digits=2),"(",
        format(coeffs[i,"pval"],digits=2),")"
        ,sep=""),sep=",")
    }
  }
  return(coeffs_v)
}

############################################################################
############################################################################
############################################################################

# Prepare the datasets for the different analyses below
# Load the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata
load("PADB_dataset_level_meta_analysis_data.RData")

# Get the cleaned, filtered datasets
datasets = list()
datasets[["acute,muscle"]] = lapply(acute_gene_tables,clean_acute_table,tissue="muscle")
datasets[["acute,blood"]] = lapply(acute_gene_tables,clean_acute_table,tissue="blood")
datasets[["longterm,muscle"]] = lapply(longterm_gene_tables,clean_longterm_table,tissue="muscle")
datasets[["longterm,blood"]] = lapply(longterm_gene_tables,clean_longterm_table,tissue="blood")

# Reshape data for replication analysis
rep_datasets = lapply(datasets,function(x)t(sapply(x,function(y)as.numeric(y$p))))
for(nn in names(rep_datasets)){
  colnames(rep_datasets[[nn]]) = 
    paste(rownames(datasets[[nn]][[1]]),datasets[[nn]][[1]]$V1,
          datasets[[nn]][[1]]$time,sep=";")
}
sapply(rep_datasets,dim)
sapply(rep_datasets,rownames)
sapply(rep_datasets,colnames)
# Get the untrained controls datasets
untrained_datasets = list()
untrained_datasets[["acute,muscle"]] = lapply(acute_gene_tables_raw,get_untrained_table,tissue="muscle")
untrained_datasets[["acute,blood"]] = lapply(acute_gene_tables_raw,get_untrained_table,tissue="blood")
untrained_datasets[["longterm,muscle"]] = lapply(longterm_gene_tables_raw,get_untrained_table,tissue="muscle")
untrained_datasets[["longterm,blood"]] = lapply(longterm_gene_tables_raw,get_untrained_table,tissue="blood")
sapply(untrained_datasets,function(x)dim(x[[1]]))

############################################################################
############################################################################
############################################################################
# Simple random effects analysis - to understand the data
simple_REs = lapply(datasets, function(x)lapply(x,run_simple_re))
simple_REs_untrained = lapply(untrained_datasets, function(x)lapply(x,run_simple_re))
simple_RE_pvals = lapply(simple_REs,function(x)sapply(x,try_get_field,fname="QMp"))
simple_RE_I2s = lapply(simple_REs,function(x)sapply(x,try_get_field,fname="I2"))
simple_RE_tau2s = lapply(simple_REs,function(x)sapply(x,try_get_field,fname="tau2"))
simple_RE_beta = lapply(simple_REs,function(x)sapply(x,try_get_field,fname="beta"))
st_pvals = lapply(datasets,function(x)lapply(x,simple_stouffer_meta_analysis))
simple_REs_untrained_beta = lapply(simple_REs_untrained,function(x)sapply(x,try_get_field,fname="beta"))
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  pvalue_qqplot(simple_RE_pvals[[nn]],main=nn,pch=20,cex=0.5)
}
for(nn in names(simple_RE_pvals)){
  r = hist(simple_RE_I2s[[nn]],plot=F)
  per = sum(simple_RE_I2s[[nn]]>70,na.rm = T)/length(simple_RE_I2s[[nn]])
  per = format(per*100,digits=3)
  cols=rep("white",length(r$mids))
  cols[r$mids > 70] = "blue"
  hist(simple_RE_I2s[[nn]],main=paste(nn,"(",per,"%)",sep=""),xlab = "I^2(%)",col=cols)
}
for(nn in names(simple_RE_pvals)){
  pvalue_qqplot(unlist(st_pvals[[nn]]),main=nn,pch=20,cex=0.5)
}
for(nn in names(simple_RE_pvals)){
  x1 = simple_RE_I2s[[nn]]
  x2 = log(simple_RE_pvals[[nn]])
  inds = !is.na(x1) & !is.na(x2)
  print(nn)
  print(cor(x1[inds],x2[inds],method="spearman"))
}
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  x1 = simple_RE_I2s[[nn]]
  x2 = simple_RE_tau2s[[nn]]
  inds = !is.na(x1) & !is.na(x2)
  print(nn)
  print(cor(x1[inds],x2[inds],method="spearman"))
  plot(x=x1,y=x2,xlab="I2(%)",ylab="Tau^2",pch=20,cex=0.2,main=nn)
}
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  x1 = simple_RE_I2s[[nn]]
  x2 = simple_RE_beta[[nn]]
  inds = !is.na(x1) & !is.na(x2)
  print(nn)
  print(cor(x1[inds],x2[inds],method="spearman"))
  plot(x=x1,y=x2,xlab="I2(%)",ylab="beta",pch=20,cex=0.2,main=nn)
}

# Example for high I2 and low tau
# TODO: revise and add some stats to the report
nn  = 2
table(simple_RE_I2s[[nn]] < 50 & abs(simple_RE_beta[[nn]]) < 0.25)
selected_is = names(which(simple_RE_I2s[[nn]] > 60 &
    simple_RE_tau2s[[nn]] < 0.1 & abs(simple_RE_beta[[nn]]) > 0.25))
selected_is = names(which(simple_RE_I2s[[nn]] > 90 &
    simple_RE_tau2s[[nn]] < 0.01))
selected_is = names(which(simple_RE_I2s[[nn]] < 25 &
    simple_RE_pvals[[nn]] < 0.00001 & abs(simple_RE_beta[[nn]]) > 0.25))
selected_is = names(which(
    simple_RE_pvals[[nn]] < 1e-8 & abs(simple_RE_beta[[nn]]) > 0.25))

selected_i = sample(selected_is)[1]
forest(simple_REs[[nn]][[selected_i]])
median(as.numeric(datasets[[nn]][[selected_i]]$p))
simple_REs[[nn]][[selected_i]]

# Filters 1 and 2 of the analysis:
to_rem1 = sapply(rep_datasets,function(x)apply(x<0.05,1,sum,na.rm=T) <= 1)
to_rem2 = list()
for(nn in names(simple_RE_beta)){
  if(grepl("acute",nn)){
    to_rem2[[nn]] = simple_RE_I2s[[nn]] < 50 & (abs(simple_RE_beta[[nn]]) < 0.25 | simple_RE_pvals[[nn]] > 0.05)
  }
  else{
    to_rem2[[nn]] = simple_RE_I2s[[nn]] < 50 & (abs(simple_RE_beta[[nn]]) < 0.1 | simple_RE_pvals[[nn]] > 0.05)
  }
  print(all(names(to_rem1[[nn]])==names(to_rem2[[nn]]),na.rm=T))
  print(table(to_rem1[[nn]] | to_rem2[[nn]])/length(to_rem2[[nn]]))
  print(sum(to_rem1[[nn]] | to_rem2[[nn]],na.rm = T))
}

rm(acute_gene_tables_raw);rm(acute_gene_tables)
rm(longterm_gene_tables_raw);rm(longterm_gene_tables)
save.image(file="workspace_before_rep_analysis.RData")

# Create the input files for meta-regression and replication analysis
meta_reg_datasets = list()
for(nn in names(simple_RE_beta)){
  print(all(names(to_rem1[[nn]])==names(to_rem2[[nn]]),na.rm=T))
  to_rem = to_rem1[[nn]] | to_rem2[[nn]]
  curr_genes = names(which(!to_rem))
  meta_reg_datasets[[nn]] = datasets[[nn]][curr_genes]
}
sapply(meta_reg_datasets,length)
meta_reg_to_mods = list(
  "acute,muscle" = c("time","training"),
  "acute,blood" = c("time","prop_males"),
  "longterm,muscle" = c("training","time","avg_age","prop_males"),
  "longterm,blood" = "training"
)
save(meta_reg_datasets,meta_reg_to_mods,
     untrained_datasets,file="meta_analysis_input.RData")

############################################################################
############################################################################
############################################################################
# Meta-regression and model selection for selected genes
library(parallel)

I2_thr = 50
AIC_diff_thr = 5
ACUTE_beta_thr = 0.25
LONGTERM_beta_thr = 0.1
P_thr = 0.0001

# Load the results (run on sherlock) instead of running the code above
load("meta_analysis_results.RData")
load("workspace_before_rep_analysis.RData")

# Algorithm for selecting genes from each meta-reg analysis
analysis2selected_genes = list()
analysis2selected_genes_stats = list()
for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]]
  pvals = sapply(analysis1,function(x)x[[1]]$mod_p)
  i2s = sapply(analysis1,function(x)x[[1]]$I2)
  aic_diffs = sapply(analysis1,function(x)unname(x[[1]]$aic_c - x$`simple:base_model`$aic_c))
  if(grepl("acute",nn)){
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>ACUTE_beta_thr))
  }
  else{
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>LONGTERM_beta_thr))  
  }
  i2s_base = simple_RE_I2s[[nn]][names(i2s)]
  pvals_base = simple_RE_pvals[[nn]][names(i2s)]
  effects_base = simple_RE_beta[[nn]][names(i2s)]
  # gene set 1: low I2 before inference, significant base p and high effects
  if(grepl("acute",nn)){
    gene_set1 = names(which(i2s_base < I2_thr & pvals_base < P_thr & abs(effects_base)>ACUTE_beta_thr))
  }
  else{
    gene_set1 = names(which(i2s_base < I2_thr & pvals_base < P_thr & abs(effects_base)>LONGTERM_beta_thr))
  }
  gene_set2 = names(which(aic_diffs <= -AIC_diff_thr & model2beta & pvals < P_thr))
  gene_set2_names = sapply(analysis1[gene_set2],function(x)names(x)[1])
  gene_set2_names = sapply(gene_set2_names,function(x)strsplit(x,split=":")[[1]][2])
  gene_set2_names[gene_set1] = "base_model"
  # is.element("10891",set=gene_set2)
  # gene_set2_names["10891"]
  # table(gene_set2_names)
  analysis2selected_genes[[nn]] = gene_set2_names
  coeffs = lapply(analysis1[names(gene_set2_names)],function(x)x[[1]]$coeffs)
  coeffs_v = sapply(coeffs,get_coeffs_str)
  m = cbind(
    unlist(names(gene_set2_names)), # entrez gene id
    unlist(entrez2symbol[names(gene_set2_names)]), # gene symbol
    unlist(gene_set2_names), # gene group
    pvals[names(gene_set2_names)], # model's p-value
    aic_diffs[names(gene_set2_names)], # AICc difference
    coeffs_v # details about the coefficients
  )
  # Correct the table above: base model should have base p-values, 
  # and a comment about the AIC difference
  base_model_genes = is.element(m[,1],set=gene_set1)
  if(sum(base_model_genes)>0){
    m[base_model_genes,4] = pvals_base[m[base_model_genes,1]]
    m[base_model_genes,5] = paste(m[base_model_genes,5],"*",sep="")
  }
  colnames(m)= c("Entrez","Symbol","Group","Model pvalue","AICc diff","Coefficients")
  analysis2selected_genes_stats[[nn]] = m
}
sapply(analysis2selected_genes,length)
analysis2selected_genes[[3]]

# Validation analyses of the selected gene sets
# 1. Replication analysis
# Load SCREEN's results
scr_path = "/Users/David/Desktop/MoTrPAC/PA_database/screen_res/"
pvals_files = list.files(scr_path)
pvals_files = pvals_files[grepl("pvals.txt",pvals_files)]
screen_results = list()
for (ff in pvals_files){
  currname = gsub("_pvals.txt","",ff)
  screen_results[[currname]] = list()
  for(m in c("bum","znormix")){
    outfile = paste(scr_path,currname,"_",m,".txt",sep="")
    screen_results[[currname]][[m]] = read.delim(outfile,row.names = 1,header=T)
  }
}
names(screen_results) = gsub("_",",",names(screen_results))
save(screen_results,file=paste(scr_path,"screen_results.RData",sep=""))
# Look at the FDR values of the selected gene sets vs. the others
par(mfrow=c(3,2),mar=c(4,4,1,1),
    cex.lab=1.1,cex.axis=1,cex=1,cex.main=1.05,lwd=0.5,mgp=c(2.1,1,0),
    las=3)
for(nn in names(analysis2selected_genes)[1:3]){
  x = screen_results[[nn]]$bum
  s = names(analysis2selected_genes[[nn]])
  s_c = setdiff(rownames(x),s)
  l1 = list();l2=list()
  for(j in 1:8){
    l1[[as.character(j+1)]] = x[s,j]
    l2[[as.character(j+1)]] = x[s_c,j]
  }
  boxplot(l1,las=2,col="blue",ylab="local fdr",main=paste(nn,": selected",sep=""),
          pch=20,ylim=c(0,1),xlab="Number of cohorts")
  abline(h = 0.2,lty=2,col="red",lwd=1.5)
  boxplot(l2,las=2,ylab="local fdr",main=paste(nn,": other",sep=""),
          pch=20,ylim=c(0,1),xlab="Number of cohorts")
  abline(h = 0.2,lty=2,col="red",lwd=1.5)
}

# 2. Look at effects in untrained
par(mfrow=c(1,3),mar=c(10,4,2,2),cex.lab=1.2,cex.axis=1.2,cex=1,cex.main=1.05)
for(nn in names(analysis2selected_genes)[1:3]){
  x = simple_REs_untrained_beta[[nn]]
  s = names(analysis2selected_genes[[nn]])
  s_c = setdiff(names(x),s)
  analysis1 = all_meta_analysis_res[[nn]][s]
  y = sapply(analysis1,function(x)max(abs(x[[1]]$coeffs[,1])))
  p = wilcox.test(x[s],y,paired=T)$p.value
  p = format(p,digits=2)
  l = list(
    untrained = abs(x[s]),exercise=y,"other genes" = abs(x[s_c])
  )
  cols = c("blue","white","red")
  boxplot(l,las=2,col=cols,horizontal=F,
          ylab="Fold change",main = paste(nn," p=",p,sep=""),
          pch=20)
}

# 3. GO enrichments
bg = unique(c(unlist(sapply(simple_REs,names))))
gs = lapply(analysis2selected_genes[1:3],names)
go_res = run_topgo_enrichment_fisher(
  gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
go_res_fdr = go_res[go_res$go_qvals < 0.1,]
table(go_res_fdr$setname)

par(mfrow=c(2,2))
for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]]
  model_names = sapply(analysis1,function(x)names(x)[1])
  table(model_names)
  ps1 = sapply(analysis1,function(x)x[[1]]$mod_p)
  i2s = sapply(analysis1,function(x)x[[1]]$I2)
  i2s_2 = simple_RE_I2s[[nn]][names(i2s)]
  aic_diffs = sapply(analysis1,function(x)x[[1]]$aic_c - x[[2]]$aic_c)
  i2_diffs = i2s-i2s_2
  # pvalue_qqplot(ps1);abline(0,1)
  plot(ps1,i2s)
  plot(i2s_2,i2s)
  plot(aic_diffs,i2_diffs)
  print(length(ps1))
  print(table(i2s<50 & ps1 < 0.0001))
}

# Examples of genes for the paper
curr_genes = analysis2selected_genes$`acute,muscle`
curr_genes[curr_genes=="base_model"]
#PGC1a
gene = "10891"
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,muscle`[[gene]]
gdata = meta_reg_datasets$`acute,muscle`[[gene]]
curr_m$slab.null = F
curr_times = rep("early",nrow(gdata))
curr_times[gdata$time==24] = "late(>12h)"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)
boxplot(yi~time,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",names=c("early","late"),ylab="Fold change")
#GJA1
gene = "2697"
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,muscle`[[gene]]
curr_m$I2
gdata = meta_reg_datasets$`acute,muscle`[[gene]]
curr_m$slab.null = F
curr_times = rep("early",nrow(gdata))
curr_times[gdata$time==24] = "late(>12h)"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)
#GJA1
gene = "85366"
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,muscle`[[gene]]
curr_m$I2
gdata = meta_reg_datasets$`acute,muscle`[[gene]]
curr_m$slab.null = F
curr_times = rep("early",nrow(gdata))
curr_times[gdata$time==24] = "late(>12h)"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)

# known_genes = c("PPARGC1A","COX1","NDUFA","PDK4","VEGFA","KDR","THY1","MYL4",
#                 "MYH1","COL1A1","ACTC1","TNNT2","GADD45G","MMP9","NR4A1")

# Prepare supplementary tables for the results
system(paste("mkdir","supp_tables"))
supp_path = paste(getwd(),"supp_tables/",sep="/")

supp_table_genes = c()
for(nn in names(analysis2selected_genes_stats)){
  m = analysis2selected_genes_stats[[nn]]
  m = cbind(rep(nn,nrow(m)),m)
  colnames(m)[1] = "Discovered in"
  supp_table_genes = rbind(supp_table_genes,m)
}
write.table(supp_table_genes,file=paste(supp_path,"supp_table_genes.txt",sep=""),
            sep="\t",quote=F,col.names = T,row.names = F)


supp_table_enrichments = go_res_fdr[,c(1:4,9)]
colnames(supp_table_enrichments)[5] = "q-value"
colnames(supp_table_enrichments)[1] = "Discovered in"
write.table(supp_table_enrichments,file=paste(supp_path,"supp_table_enrichments.txt",sep=""),
            sep="\t",quote=F,col.names = T,row.names = F)
