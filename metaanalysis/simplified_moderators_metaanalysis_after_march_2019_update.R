# Our algorithm for analysis of a single gene
# Input: datasets for time point t and gene g
# 1. Run meta-analysis for endurence and resistance (and both)
# 2. Run meta-analysis for controls
# 3. Get controls intercept estimation beta_c
# 4. Get current gene's overall diff expression: beta_e = max(beta_0,beta_0 + beta_endurence)
# 5. SCORE 1: abs(beta_e-beta_c)
# 6. SCORE 2: Egger's test of the exercise meta-analysis
# 7. Return SCORE 1, SCORE 2, and the betas (and their significance) from the exercise meta-analysis
library(org.Hs.eg.db);library(metafor)
source('/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)

setwd('/Users/David/Desktop/MoTrPAC/project_release_feb_2018/data/')
# Prepare the datasets for the different analyses below
# Load the datasets and their metadata
load("human_ge_cohort_preprocessed_db_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
acute_sample2time = sample2time
load("human_ge_cohort_preprocessed_db_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata
longterm_sample2time = sample2time
load("human_ge_cohort_preprocessed_db_gene_tables.RData")

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

# New preprocessing methods: after the March 2019 update
# New analyses:
#   Acute: time is binned to immediate (<=1h), early (1-5h), and late(>=20h)
#   Age is kept as numeric
clean_acute_table <-function(gdata,tissue="muscle",remove_untrained=T){
  if(remove_untrained){
    gdata = gdata[!is.element(gdata$training,set=c("yoga","control","untrained")),]
  }
  gdata = gdata[gdata$tissue==tissue,]
  gdata$vi = pmax(gdata$vi,1e-5)
  newtime = rep(2,nrow(gdata))
  newtime[gdata$time<=1] = 1
  newtime[gdata$time>=20 ] = 3
  gdata$time = ordered(newtime)
  gdata = simplify_training_for_exercise_analysis(gdata)
  gdata$sdd = as.numeric(gdata$sdd)
  return(gdata)
}
# Long term data preprocessing:
#   Binarize time based on 150 days
clean_longterm_table <-function(gdata,tissue="muscle",remove_untrained=T){
  if(remove_untrained){
    gdata = gdata[!is.element(gdata$training,set=c("yoga","control","untrained")),]
  }
  gdata = gdata[gdata$tissue==tissue,]
  gdata$vi = pmax(gdata$vi,1e-5)
  # time analysis
  newtime = rep(2,nrow(gdata))
  newtime[gdata$time<150] = 1
  gdata$time = ordered(newtime)
  gdata = simplify_training_for_exercise_analysis(gdata)
  gdata$sdd = as.numeric(gdata$sdd)
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
  if(is.null(dim(coeffs))){coeffs=as.matrix(coeffs,nrow=1)}
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

# Get the cleaned, filtered datasets
datasets = list()
datasets[["acute,muscle"]] = lapply(acute_gene_tables,clean_acute_table,tissue="muscle")
datasets[["acute,blood"]] = lapply(acute_gene_tables,clean_acute_table,tissue="blood")
datasets[["longterm,muscle"]] = lapply(longterm_gene_tables,clean_longterm_table,tissue="muscle")
datasets[["longterm,blood"]] = lapply(longterm_gene_tables,clean_longterm_table,tissue="blood")

# A new filter added after the update on March 2019: exclude
# genes with extremely low number of studies from each analysis
par(mfrow=c(2,2))
for(nn in names(datasets)){
  num_datasets = sapply(datasets[[nn]],function(x)length(unique(x$gse)))
  hist(num_datasets,main=nn)
  to_rem = (num_datasets/max(num_datasets)) < 0.75
  print(table(to_rem))
  datasets[[nn]] = datasets[[nn]][!to_rem] 
}
sapply(datasets,length)

# Reshape data for replication analysis
get_gene_data_for_rep_analysis<-function(gdata,exclude){
  gdata = gdata[!is.element(gdata$V1,set=exclude),]
  x = gdata$p
  names(x)=paste(gdata$V1,gdata$time,sep=";")
  return(x)
}
load("human_ge_gene_coverage_analysis.RData")
rep_datasets = list()
for(nn in names(datasets)){
  m = sapply(datasets[[nn]][all_genes],get_gene_data_for_rep_analysis,
             exclude=low_coverage_platforms)
  print(dim(m))
  rep_datasets[[nn]] = t(m)
  print(colnames(rep_datasets[[nn]]))
}
sapply(rep_datasets,dim)
sapply(rep_datasets,function(x)all(rownames(x)==all_genes))
sapply(rep_datasets,colnames)

# Get the untrained controls datasets
untrained_datasets = list()
untrained_datasets[["acute,muscle"]] = lapply(acute_gene_tables,get_untrained_table,tissue="muscle")
untrained_datasets[["acute,blood"]] = lapply(acute_gene_tables,get_untrained_table,tissue="blood")
untrained_datasets[["longterm,muscle"]] = lapply(longterm_gene_tables,get_untrained_table,tissue="muscle")
untrained_datasets[["longterm,blood"]] = lapply(longterm_gene_tables,get_untrained_table,tissue="blood")
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

# Filters 1 exclude genes with no sig p at 0.01
# before the update on March 2019 the threshold was 0.05 and 1
# Returns true if the gene has p<thr in at least num studies
rep_filter <-function(gdata,num=1,thr=0.01){
  return(sum(gdata$p<=thr,na.rm = T)>=num)
}
to_rem1 = sapply(datasets,function(x)!sapply(x,rep_filter,num=2,thr=0.05))

rm(acute_gene_tables_raw);rm(acute_gene_tables);
rm(longterm_gene_tables_raw);rm(longterm_gene_tables)
save.image(file="workspace_before_rep_analysis.RData")

# Create the input files for meta-regression and replication analysis
meta_reg_datasets = list()
for(nn in names(simple_RE_beta)){
  curr_genes = names(which(!to_rem1[[nn]]))
  meta_reg_datasets[[nn]] = datasets[[nn]][curr_genes]
}
sapply(meta_reg_datasets,length)

# After the data update of March 2019
# Not enough non-male acute muscle studies
# Not enough resistance training datasets for acute, blood
# For longterm blood we take training only as we assume that it will
# be the primary cause of differential abundance
meta_reg_to_mods = list(
  "acute,muscle" = c("time","training","avg_age"),
  "acute,blood" = c("time","prop_males","avg_age"),
  "longterm,muscle" = c("training","time","avg_age","prop_males"),
  "longterm,blood" = "training"
)
save(meta_reg_datasets,meta_reg_to_mods,
     untrained_datasets,file="meta_analysis_input.RData")

############################################################################
############################################################################
############################################################################
# Meta-regression and model selection for selected genes
library(parallel);library(metafor)

# Load the results (run on sherlock) instead of running the code above
load("meta_analysis_results.RData")
load("workspace_before_rep_analysis.RData")
load("meta_analysis_input.RData")
source('/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')

I2_thr = 50
AIC_diff_thr = 5
ACUTE_beta_thr = 0.25
LONGTERM_beta_thr = 0.1
P_thr = 1e-04

# A few helper methods to work with the meta-analysis output
get_aicc_diff<-function(x){
  if(is.element("simple:base_model",set=names(x)) &&
     is.element("aic_c",set=names(x[["simple:base_model"]]))){
    val = x[[1]]$aic_c - x$`simple:base_model`$aic_c
    if(is.null(val) && is.na(val) || is.nan(val)){return(0)}
    if(is.infinite(val)){return(-1000)}
    return(unname(val))
  }
  return(0)
}

# Algorithm for selecting genes from each meta-reg analysis
analysis2selected_genes = list()
analysis2selected_genes_stats = list()
all_pvals = c()
for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]]
  pvals = sapply(analysis1,function(x)x[[1]]$mod_p)
  all_pvals = c(all_pvals,pvals)
  i2s = simple_RE_I2s[[nn]][names(pvals)]
  # separate into two gene sets: those that passed the aic diff test vs. those that did not
  # define the set of filters
  # 1. AICc filter
  aic_diffs = sapply(analysis1,get_aicc_diff)
  genes_with_high_aic_diff = aic_diffs <= -AIC_diff_thr
  # 2. Beta filter
  if(grepl("acute",nn)){
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>ACUTE_beta_thr))
  }
  else{
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>LONGTERM_beta_thr))  
  }
  # 3. Is the top model simple base
  is_base_model = sapply(analysis1,function(x)names(x)[1] =="simple:base_model")
  # 4. Pval filter
  pval_filter = pvals <= P_thr
  # 5. I2 filter
  i2_filter = i2s <= I2_thr
  which(genes_with_high_aic_diff & i2_filter & pval_filter & model2beta)
  
  selected_aic_diff_genes = names(aic_diffs)[genes_with_high_aic_diff & model2beta & pval_filter]
  selected_base_model_genes = names(aic_diffs)[i2_filter & model2beta & pval_filter]
  selected_base_model_genes = setdiff(selected_base_model_genes,selected_aic_diff_genes)
  curr_selected_genes = union(selected_aic_diff_genes,selected_base_model_genes)
  curr_selected_genes_names = sapply(analysis1[selected_aic_diff_genes],function(x)names(x)[1])
  curr_selected_genes_names = sapply(curr_selected_genes_names,function(x)strsplit(x,split=":")[[1]][2])
  curr_selected_genes_names[selected_base_model_genes] = "base_model"
  analysis2selected_genes[[nn]] = curr_selected_genes_names
  coeffs = lapply(analysis1[names(curr_selected_genes_names)],function(x)x[[1]]$coeffs)[selected_aic_diff_genes]
  coeffs[selected_base_model_genes] = lapply(simple_REs[[nn]][selected_base_model_genes],
                                             function(x){y=cbind(x$beta,x$pval);colnames(y)=c("beta","pval");y})
  coeffs = coeffs[curr_selected_genes]
  coeffs_v = sapply(coeffs,get_coeffs_str)
  m = cbind(
    unlist(names(curr_selected_genes_names)), # entrez gene id
    unlist(entrez2symbol[names(curr_selected_genes_names)]), # gene symbol
    unlist(curr_selected_genes_names), # gene group
    pvals[names(curr_selected_genes_names)], # model's p-value
    aic_diffs[names(curr_selected_genes_names)], # AICc difference
    coeffs_v # details about the coefficients
  )
  colnames(m)= c("Entrez","Symbol","Group","Model pvalue","AICc diff","Coefficients")
  analysis2selected_genes_stats[[nn]] = m
}
sapply(analysis2selected_genes,length)

# Sanity checks and tests: our p-value threshold is lower than BY correction
max(all_pvals[p.adjust(all_pvals,method = "BY")<0.05])>1e-4

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
# par(mfrow=c(1,3),mar=c(10,4,2,2),cex.lab=1.2,cex.axis=1.2,cex=1,cex.main=0.9)
# for(nn in names(analysis2selected_genes)[1:3]){
#   x = simple_REs_untrained_beta[[nn]]
#   analyzed_genes = names(meta_reg_datasets[[nn]])
#   x = x[intersect(names(x),analyzed_genes)]
#   s = names(analysis2selected_genes[[nn]])
#   s_c = setdiff(names(x),s)
#   analysis1 = all_meta_analysis_res[[nn]][s]
#   y = sapply(analysis1,function(x)max(abs(x[[1]]$coeffs[,1])))
#   p = wilcox.test(x[s],y,paired=T)$p.value
#   p = format(p,digits=2)
#   l = list(
#     "untrained effects" = abs(x[s]),"exercise effects"=y,"other genes,\n untrained effects" = abs(x[s_c])
#   )
#   cols = c("blue","white","red")
#   boxplot(l,las=2,col=cols,horizontal=F,cex.axis=0.9,
#           ylab="Fold change",main = paste(nn," p=",p,sep=""),pch=20)
# }

# 2. alt, Look at effects in untrained - excluded untrained others
par(mar=c(8,4,2,2),cex.lab=1.2,cex.axis=1.2)
l = list();cols=c()
for(nn in names(analysis2selected_genes)[1:3]){
  x = simple_REs_untrained_beta[[nn]]
  analyzed_genes = names(meta_reg_datasets[[nn]])
  x = x[intersect(names(x),analyzed_genes)]
  s = names(analysis2selected_genes[[nn]])
  s_c = setdiff(names(x),s)
  analysis1 = all_meta_analysis_res[[nn]][s]
  y = sapply(analysis1,function(x)max(abs(x[[1]]$coeffs[,1])))
  p = wilcox.test(x[s],y,paired=T)$p.value
  p = format(p,digits=2)
  l[[paste(nn,"untrained effects",sep="\n")]] = abs(x[s])
  l[[paste(nn,"exercise effects",sep="\n")]] = y
  cols = c(cols,c("blue","red"))
}
boxplot(l,las=2,col=cols,horizontal=F,cex.axis=1.05,pch=20,ylim=c(0,1.6),ylab="Fold change")

# 3. GO enrichments
bg = unique(c(unlist(sapply(simple_REs,names))))
gs = lapply(analysis2selected_genes,names)
go_res = run_topgo_enrichment_fisher(
  gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
go_res1 = go_res[go_res$Annotated < 1500,]
go_res1$classicFisher[is.na(as.numeric(go_res1$classicFisher))] = 1e-30
go_res1 = go_res[go_res1$Significant > 3,]
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
table(go_res_fdr$setname)
gene_group_enrichments = go_res
gene_group_enrichments_fdr = go_res_fdr
get_most_sig_enrichments_by_groups(gene_group_enrichments_fdr,num=5)[,1:4]

gs = lapply(analysis2selected_genes,names)
gs = lapply(gs,function(x,y)y[x],y=unlist(entrez2symbol))

par(mfrow=c(2,2))
for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]]
  model_names = sapply(analysis1,function(x)names(x)[1])
  table(model_names)
  ps1 = sapply(analysis1,function(x)x[[1]]$mod_p)
  i2s = sapply(analysis1,function(x)x[[1]]$I2)
  i2s_2 = simple_RE_I2s[[nn]][names(i2s)]
  aic_diffs = abs(sapply(analysis1,function(x)x[[1]]$aic_c - x[[2]]$aic_c))
  hist(aic_diffs,main=nn)
}

# Examples of genes for the paper
curr_genes = analysis2selected_genes$`acute,muscle`
curr_genes[curr_genes=="base_model"]
#PGC1a
gene = "10891"
curr_genes[gene]
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,muscle`[[gene]]
gdata = meta_reg_datasets$`acute,muscle`[[gene]]
curr_m$slab.null = F
curr_times = rep("0-1h",nrow(gdata))
curr_times[gdata$time==2] = "2-5h"
curr_times[gdata$time==3] = ">20h"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
analysis1[[gene]][[1]]$mod_p
aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
dev.off()
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)
boxplot(yi~time,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",names=c("0-1h","2-5h",">20h"),ylab="Fold change")

# Selected examples, before the March 2019 update
gene = "1282"
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`longterm,muscle`[[gene]]
all_meta_analysis_res$`longterm,muscle`[[gene]]
curr_m$I2
gdata = meta_reg_datasets$`longterm,muscle`[[gene]]
curr_m$slab.null = F
curr_times = rep("< 150 days",nrow(gdata))
curr_times[gdata$time==24] = "> 150 days"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)

############################################################################
############################################################################
############################################################################
# Interpretation of the results
library(topGO)

# some helper functions for reformatting the data
get_ts<-function(gdata){
  v = as.numeric(gdata$tstat)
  ns = paste(gdata$V1,gdata$training,gdata$time,
             format(gdata$avg_age,digits=1),format(gdata$prop_males,digits=1),sep=";")
  v = tapply(v,ns,mean)
  return(v)
}
get_t_matrix_from_list<-function(tstats){
  all_ns = unique(unlist(sapply(tstats,names)))
  all_genes = names(tstats)
  m = matrix(0,nrow=length(all_genes),ncol=length(all_ns),
             dimnames = list(all_genes,all_ns))
  for(i in 1:length(all_genes)){
    m[i,names(tstats[[i]])] = tstats[[i]]
  }
  return(m)
}

# Get effect matrices - mean responses - t statistics
# Reminder: the resulting matrices may have many zeroes because
# we basically merge all available t-statistics and not all genes
# are represented in all studies.
dataset_tstats = lapply(datasets,function(x)sapply(x,get_ts))
mean_effect_matrices = lapply(dataset_tstats,get_t_matrix_from_list)
sapply(mean_effect_matrices,dim)
sapply(mean_effect_matrices,function(x)table(x==0))
for(nn in names(mean_effect_matrices)){
  colnames(mean_effect_matrices[[nn]]) = 
    gsub("endurance","E",colnames(mean_effect_matrices[[nn]]),ignore.case = T)
  colnames(mean_effect_matrices[[nn]]) = 
    gsub("resistance","R",colnames(mean_effect_matrices[[nn]]),ignore.case = T)
  colnames(mean_effect_matrices[[nn]]) = 
    gsub("GE_","",colnames(mean_effect_matrices[[nn]]),ignore.case = T)
  colnames(mean_effect_matrices[[nn]]) = 
    gsub(" ","",colnames(mean_effect_matrices[[nn]]),ignore.case = T)
}
sapply(mean_effect_matrices,colnames)

# helper function for getting the clusters
get_num_clusters_wss_kmeans<-function(data,k.max=10,wss_imp_thr=0.6){
  if(ncol(m)>2){
  data = t(scale(t(data)))
  }
  k.max = min(k.max,nrow(data)/2)
  wss <- sapply(1:k.max, 
                function(k){kmeans(data, k, nstart=100,iter.max = 500 )$tot.withinss})
  plot(1:k.max, wss,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  wss_imp_factors = wss[-1]/wss[-length(wss)]
  if(all(wss_imp_factors>wss_imp_thr)){return(1)}
  k = 2
  for(j in 2:length(wss_imp_factors)){
    if(wss_imp_factors[j] > wss_imp_thr){break}
    k = k+1
  }
  return(c(k=k))
}

gene_subgroups = list()
gene_t_patterns = list()
par(mfrow=c(4,4))
for(nn in names(analysis2selected_genes_stats)){
  curr_genes = analysis2selected_genes_stats[[nn]]
  curr_groups = curr_genes[,"Group"]
  for(gg in unique(curr_groups)){
    curr_m = as.matrix(curr_genes[gg==curr_genes[,"Group"],])
    if(ncol(curr_m)==1){curr_m=t(curr_m)}
    m = as.matrix(mean_effect_matrices[[nn]][curr_m[,"Entrez"],])
    if(ncol(m)==1){m=t(m)}
    rownames(m) = curr_m[,"Entrez"]
    
    # Merge based on the moderator's average not the t-test scores themselves
    if(gg != "base_model"){
      curr_mods = strsplit(gg,split=";")[[1]]
      curr_m_meta = sapply(colnames(m),function(x)strsplit(x,split=";")[[1]])
      rownames(curr_m_meta) = c("cohort","training","time","avg_age","prop_males")
      curr_m_meta = curr_m_meta[curr_mods,]
      if(!is.null(dim(curr_m_meta))){
        curr_m_meta = apply(curr_m_meta,2,paste,collapse=";")
      }
      m = t(apply(m,1,function(x,y)tapply(x,y,mean),y=curr_m_meta))
      if(nrow(m)>5){currk = get_num_clusters_wss_kmeans(m,10)}
      else{currk=2}
    }
    else{
      m = as.matrix(rowMeans(m),ncol=1)
      colnames(m)[1] = "base,meant"
      currk=2
    }

    print(paste(nn,gg,currk))
    if(nrow(m)>2){
      m_kmeans = kmeans(m,centers = currk)
      m_kmeans = m_kmeans$cluster
    }
    else{
      m_kmeans = rep(1,nrow(m))
      names(m_kmeans) = curr_m[,"Entrez"]
    }
    m = cbind(m,m_kmeans)
    colnames(m)[ncol(m)] = "kmeans_clusters"
    gene_t_patterns[[paste(nn,gg,sep=",")]] = m
    
    for(kk in unique(m_kmeans)){
      currname = paste(nn,gg,kk,sep=",")
      gene_subgroups[[currname]] = names(m_kmeans)[m_kmeans==kk]
    }
  }
}
sort(sapply(gene_subgroups,length))
base_model_ms = c()
for(gg in names(gene_t_patterns)[grepl("base_model",names(gene_t_patterns))]){
  currm = gene_t_patterns[[gg]]
  currm = cbind(rep(gg,nrow(currm)),currm)
  base_model_ms = rbind(base_model_ms,currm)
}
gene_t_patterns[["base_models"]] = base_model_ms

# Enrichment analysis of the new groups
bg = unique(c(unlist(sapply(simple_REs,names))))
gs = gene_subgroups
gs = gs[sapply(gs,length)>5]
sapply(gs,length)
go_res = run_topgo_enrichment_fisher(
  gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
go_res1 = go_res[go_res$Annotated < 1500,]
go_res1$classicFisher[is.na(as.numeric(go_res1$classicFisher))] = 1e-30
go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
table(go_res_fdr$setname)
gene_subgroup_enrichments = go_res
gene_subgroup_enrichments_fdr = go_res_fdr

table(as.character(gene_subgroup_enrichments_fdr[,1]))
get_most_sig_enrichments_by_groups(gene_subgroup_enrichments_fdr,num=2)[,c(1,4,8)]
gene_subgroup_enrichments_fdr[grepl("ossi",gene_subgroup_enrichments_fdr$Term),]
get_most_sig_enrichments_by_groups(gene_group_enrichments_fdr,num=3)[,1:4]

# Other enrichment analyses
library(ReactomePA)
reactome_pathways_groups = run_reactome_enrichment_analysis(
  lapply(analysis2selected_genes,names),universe=bg)
reactome_pathways_groups1 = reactome_pathways_groups[reactome_pathways_groups$Count>2,]
ps = reactome_pathways_groups1$pvalue
qs = p.adjust(ps,method="fdr")
reactome_pathways_groups1$qvalue = qs
reactome_pathways_groups_fdr = reactome_pathways_groups1[qs <= 0.1,]
table(reactome_pathways_groups_fdr[,1])

reactome_pathways_subgroups = run_reactome_enrichment_analysis(gene_subgroups,universe=bg)
reactome_pathways_subgroups1 = reactome_pathways_subgroups[reactome_pathways_subgroups$Count>2,]
ps = reactome_pathways_subgroups1$pvalue
qs = p.adjust(ps,method="fdr")
reactome_pathways_subgroups1$qvalue = qs
reactome_pathways_subgroups_fdr = reactome_pathways_subgroups1[qs <= 0.1,]
table(reactome_pathways_subgroups_fdr[,1])
get_most_sig_enrichments_by_groups(reactome_pathways_subgroups_fdr,pcol="pvalue",num = 4)[,c(1,3)]

reactome_pathways_subgroups_fdr[,1] = as.character(reactome_pathways_subgroups_fdr[,1])

save(gene_subgroups,analysis2selected_genes,
     gene_t_patterns,
  gene_subgroup_enrichments,gene_subgroup_enrichments_fdr,
     gene_group_enrichments,gene_group_enrichments_fdr,
     reactome_pathways_subgroups,reactome_pathways_subgroups_fdr,
     reactome_pathways_groups,reactome_pathways_groups_fdr,
     file="mm_results_march_2019_p1e5.RData")

# We want to examine which clusters to analyze, sort by number of enrichments
enriched_clusters_go = sort(table(as.character(gene_subgroup_enrichments_fdr$setname)))
enriched_clusters_reactome = sort(table(as.character(reactome_pathways_subgroups_fdr[,1])))
all_enriched_clusters = union(names(enriched_clusters_go),names(enriched_clusters_reactome))

# Age associated clusters
age_associated_clusters = all_enriched_clusters[grep("age",all_enriched_clusters)]
sex_associated_clusters = all_enriched_clusters[grep("male",all_enriched_clusters)]
training_associated_clusters = all_enriched_clusters[grep("train",all_enriched_clusters)]
sapply(gene_subgroups[age_associated_clusters],length)

# Heatmaps and line plots
# Analysis of the acute blood datasets
plot_with_err_bars<-function(xnames,avg,sdev,add=F,arrow_col="black",...){
  if(add){
    lines(avg,pch=19,...)
  }
  else{
    ylim = c(min(avg)-max(sdev),max(avg)+max(sdev))
    plot(avg,xaxt = "n",pch=19, type='l',ylim=ylim,...)
    axis(1, at=1:length(xnames), labels=xnames)
  }
  # hack: we draw arrows but with very special "arrowheads"
  arrows(1:length(xnames), avg-sdev, 1:length(xnames), avg+sdev,
         length=0.05, angle=90, code=3,col=arrow_col)
}
shorten_by_words<-function(x,num=3){
  if(is.na(x) || length(x)==0){return("")}
  arr =  strsplit(x,split=" ")[[1]]
  num = min(num,length(arr))
  return(paste(arr[1:num],collapse=" "))
}
library(gplots)

# Training-specific responses (indep of time)
set_name = "acute,muscle,training,1"
table_name = "acute,muscle"
table_name2 = "acute,muscle,training"
# Training-specific responses with time
set_name = "acute,muscle,time;training,2"
table_name = "acute,muscle"
table_name2 = "acute,muscle,time;training"

# plot the t-stat over all cohorts
set_genes = gene_subgroups[[set_name]]
mat = mean_effect_matrices[[table_name]][set_genes,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,apply(mat,2,sd)>0]
mat[mat>5]=5;mat[mat< -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.8)
# plot the average t-stats
mat = gene_t_patterns[[table_name2]][set_genes,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,-ncol(mat)]
mat[mat>5]=5;mat[mat< -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.8)
# get the enrichments of the set
gene_subgroup_enrichments_fdr[gene_subgroup_enrichments_fdr[,1]==set_name,]
reactome_pathways_subgroups_fdr[reactome_pathways_subgroups_fdr[,1]==set_name,]

# Sex-specific responses in longterm training
set_name = "longterm,muscle,prop_males,1"
gene_set = gene_subgroups[[set_name]]
mat = gene_t_patterns$`longterm,muscle,prop_males`[gene_set,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,-ncol(mat)]
# mat[mat>2]=2;mat[mat< -2]=-2
mat = mat[,colnames(mat)!="NaN"]
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.5)

set_name = "acute,blood,time;avg_age,1"
gene_set = gene_subgroups[[set_name]]
mat = gene_t_patterns$`acute,blood,time;avg_age`[gene_set,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,-ncol(mat)]
mat[mat>4]=4;mat[mat< -4]=-4
mat = mat[,colnames(mat)!="NaN"]
ages = sapply(colnames(mat),function(x)as.numeric(strsplit(x,split=";")[[1]][2]))
heatmap.2(mat[,order(ages)],trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.5)

set_name = "acute,blood,avg_age,1"
gene_set = gene_subgroups[[set_name]]
mat = gene_t_patterns$`acute,blood,avg_age`[gene_set,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,-ncol(mat)]
mat[mat>4]=4;mat[mat< -4]=-4
mat = mat[,colnames(mat)!="NaN"]
ages = sapply(colnames(mat),function(x)as.numeric(strsplit(x,split=";")[[1]][2]))
heatmap.2(mat[,order(ages)],trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.8)

# par(mfrow=c(1,2))
# gene_set = gene_subgroups$`acute,blood,time,1`
# mat = mean_effect_matrices$`acute,blood`[gene_set,]
# times = ordered(as.numeric(sapply(colnames(mat),function(x)strsplit(x,split=";")[[1]][3])))
# D = dummy(times,levelsToKeep = levels(times))
# D = t(t(D)/colSums(D))
# mat1 = mat %*% D
# boxplot(mat1,xlab="Time after bout (hours)",ylab="Avg t-statistic",col="blue")
# gene_set = gene_subgroups$`acute,blood,time,2`
# mat = mean_effect_matrices$`acute,blood`[gene_set,]
# times = ordered(as.numeric(sapply(colnames(mat),function(x)strsplit(x,split=";")[[1]][3])))
# D = dummy(times,levelsToKeep = levels(times))
# D = t(t(D)/colSums(D))
# mat1 = mat %*% D
# boxplot(mat1,xlab="Time after bout (hours)",ylab="Avg t-statistic",col ="red")

# Try another plot: select cluster names and plot all of their patterns
# along with selected enrichments

# The different patterns of acute muscle, time-dependent response
pref = "acute,muscle,time,\\d"
clusters = names(gene_subgroups)[grepl(pref,names(gene_subgroups))]
par(mfrow=c(2,3))
cols = c("blue","cyan","red","green","black")
names(cols) = clusters
for(set_name in sort(clusters)){
  gene_set = gene_subgroups[[set_name]]
  top_enrichments1 = reactome_pathways_subgroups_fdr[
    reactome_pathways_subgroups_fdr[,1]==set_name,"Description"]
  top_enrichments2 = gene_subgroup_enrichments_fdr[
    gene_subgroup_enrichments_fdr[,1]==set_name,"Term"]
  
  # remove embryo enrichments
  top_enrichments1 = top_enrichments1[!grepl("embryo",top_enrichments1)]
  top_enrichments2 = top_enrichments2[!grepl("embryo",top_enrichments2)]
  
  enrichment1 = top_enrichments1[1]
  if(any(grepl("muscle",top_enrichments1))){
    enrichment1 = top_enrichments1[grepl("muscle",top_enrichments1)][1]
  }
  enrichment2 = top_enrichments2[1]
  if(any(grepl("muscle",top_enrichments2))){
    enrichment2 = top_enrichments2[grepl("muscle",top_enrichments2)][1]
  }
  if(is.na(enrichment1)){enrichment1 = NULL}
  if(is.na(enrichment2)){enrichment2 = NULL}
  enrichment1 = shorten_by_words(tolower(enrichment1))
  enrichment2 = shorten_by_words(tolower(enrichment2))
  new_set_name = paste(set_name," (",length(gene_set)," genes)",sep="")
  curr_main = paste(new_set_name,enrichment1,enrichment2,sep="\n")
  curr_main = gsub("\n\n","\n",curr_main)
  print(curr_main)
  mat = gene_t_patterns$`acute,muscle,time`[gene_set,]
  rownames(mat) = unlist(entrez2symbol[rownames(mat)])
  mat = mat[,-ncol(mat)]
  plot_with_err_bars(c("0-1h","2-5h",">20h"),
                     colMeans(mat),apply(mat,2,sd),col=cols[set_name],lwd=3,
                     main=curr_main,ylab = "Mean t-statistic",xlab="Time",
                     cex.lab=1.2,cex.axis=1.2,arrow_col = cols[set_name])
}

# Write input to DREM
m = gene_t_patterns$`acute,muscle,time`[,1:3]
rownames(m) = entrez2symbol[rownames(m)]
write.table(m,file="drem/acute_muscle_time.txt",
            sep="\t",col.names = T,row.names = T,quote=F)

pref = "acute,blood,time,\\d"
clusters = names(gene_subgroups)[grepl(pref,names(gene_subgroups))]
par(mfrow=c(1,2))
for(set_name in sort(clusters)){
  gene_set = gene_subgroups[[set_name]]
  top_enrichments1 = reactome_pathways_subgroups_fdr[
    reactome_pathways_subgroups_fdr[,1]==set_name,"Description"]
  top_enrichments2 = gene_subgroup_enrichments_fdr[
    gene_subgroup_enrichments_fdr[,1]==set_name,"Term"]

  # remove embryo enrichments
  top_enrichments1 = top_enrichments1[!grepl("embryo",top_enrichments1)]
  top_enrichments2 = top_enrichments2[!grepl("embryo",top_enrichments2)]

  enrichment1 = top_enrichments1[1]
  if(any(grepl("muscle",top_enrichments1))){
    enrichment1 = top_enrichments1[grepl("muscle",top_enrichments1)][1]
  }
  enrichment2 = top_enrichments2[1]
  if(any(grepl("muscle",top_enrichments2))){
    enrichment2 = top_enrichments2[grepl("muscle",top_enrichments2)][1]
  }
  if(is.na(enrichment1)){enrichment1 = NULL}
  if(is.na(enrichment2)){enrichment2 = NULL}
  enrichment1 = shorten_by_words(tolower(enrichment1))
  enrichment2 = shorten_by_words(tolower(enrichment2))

  set_name = paste(set_name,"(",nrow(mat)," genes)",sep="")
  curr_main = paste(set_name,enrichment1,enrichment2,sep="\n")
  curr_main = gsub("\n\n","\n",curr_main)
  print(curr_main)
  mat = gene_t_patterns$`acute,blood,time`[gene_set,]
  rownames(mat) = unlist(entrez2symbol[rownames(mat)])
  mat = mat[,-ncol(mat)]
  mat[mat>6]=6;mat[mat< -6]=-6
  plot_with_err_bars(c("0-1h","2-5h",">20h"),
                     colMeans(mat),apply(mat,2,sd),col="blue",lwd=2,
                     main=curr_main,ylab = "avg t")
}

# pref = "longterm,muscle,time;avg_age,\\d"
# clusters = names(gene_subgroups)[grepl(pref,names(gene_subgroups))]
# par(mfrow=c(2,3))
# for(set_name in sort(clusters)){
#   gene_set = gene_subgroups[[set_name]]
#   top_enrichments1 = reactome_pathways_subgroups_fdr[
#     reactome_pathways_subgroups_fdr[,1]==set_name,"Description"]
#   top_enrichments2 = gene_subgroup_enrichments_fdr[
#     gene_subgroup_enrichments_fdr[,1]==set_name,"Term"]
#   
#   # remove embryo enrichments
#   top_enrichments1 = top_enrichments1[!grepl("embryo",top_enrichments1)]
#   top_enrichments2 = top_enrichments2[!grepl("embryo",top_enrichments2)]
#   
#   enrichment1 = top_enrichments1[1]
#   if(any(grepl("muscle",top_enrichments1))){
#     enrichment1 = top_enrichments1[grepl("muscle",top_enrichments1)][1]
#   }
#   enrichment2 = top_enrichments2[1]
#   if(any(grepl("muscle",top_enrichments2))){
#     enrichment2 = top_enrichments2[grepl("muscle",top_enrichments2)][1]
#   }
#   if(is.na(enrichment1)){enrichment1 = NULL}
#   if(is.na(enrichment2)){enrichment2 = NULL}
#   enrichment1 = shorten_by_words(tolower(enrichment1))
#   enrichment2 = shorten_by_words(tolower(enrichment2))
#   
#   curr_main = paste(set_name,enrichment1,enrichment2,sep="\n")
#   curr_main = gsub("\n\n","\n",curr_main)
#   print(curr_main)
#   mat = gene_t_patterns$`longterm,muscle,time;avg_age`[gene_set,]
#   rownames(mat) = unlist(entrez2symbol[rownames(mat)])
#   mat = mat[,-ncol(mat)]
#   mat[mat>6]=6;mat[mat< -6]=-6
#   plot_with_err_bars(colnames(mat),
#                      colMeans(mat),apply(mat,2,sd),col="blue",lwd=2,
#                      main=curr_main,ylab = "avg t")
# }

###############################################
###############################################
###############################################
# Prepare supplementary tables for the results
###############################################
###############################################
###############################################

system(paste("mkdir","supp_tables"))
supp_path = paste(getwd(),"supp_tables/",sep="/")
supp_file = paste(supp_path,"SupplementaryTables.xlsx",sep="/")
library(xlsx)

# 1. Tables presenting the cohorts
metadata_row_for_supp_table<-function(x,y,samp2time){
  curr_samps = samp2time[x$gsms]
  curr_samps = curr_samps[!is.na(curr_samps)]
  Nsample = length(curr_samps)
  curr_samps[curr_samps==min(samp2time)]="Pre"
  if(is.element("gene_fold_changes",set=names(y))){
    sample_info = strsplit(colnames(y$gene_fold_changes),split="_")
    sample_info_tps =  sapply(sample_info,function(x)x[length(x)])
    subjects =  sapply(sample_info,function(x)x[1])
    Nsubject = sum(sample_info_tps==sample_info_tps[1])
  }
  else{
    Nsubject=Nsample
  }
  tps = paste(unique(curr_samps),collapse=",")
  res = c("GSE"=x$gse,"Tissue"=x$tissue,"Training"=x$training,"Nsample"=Nsample,"Nsubject"=Nsubject,
          "Time_points"=tps,"Avg_age"=x$avg_age,"Prop_males"=x$prop_males)
  return(res)
}
m1 = c()
for(nn in names(acute_metadata)){
  m1 = rbind(m1,
        c(
          nn,metadata_row_for_supp_table(acute_metadata[[nn]],
                                         acute_datasets[[nn]],acute_sample2time)
        ))
}

m2 = c()
for(nn in names(longterm_datasets)){
    m2 = rbind(m2,
               c(
                 nn,metadata_row_for_supp_table(longterm_metadata[[nn]],
                                                longterm_datasets[[nn]],longterm_sample2time)
               ))
}
supp_table_1_all_cohorts = rbind(m1,m2)
meta_analysis_group=c()
for(nn in supp_table_1_all_cohorts[,1]){
  currg = ""
  for(nn2 in names(datasets)){
    if(is.element(nn,set=datasets[[nn2]][[1]][,1])){
      currg = nn2
    }
  }
  for(nn2 in names(untrained_datasets)){
    if(is.element(nn,set=untrained_datasets[[nn2]][[1]][,1])){
      currg = nn2
    }
  }
  meta_analysis_group[nn]=currg
}
supp_table_1_all_cohorts = cbind(supp_table_1_all_cohorts,meta_analysis_group)
colnames(supp_table_1_all_cohorts)[1]="Cohort_ID"
write.xlsx(supp_table_1_all_cohorts,file=supp_file,sheetName = "STable1",row.names = F)
length(unique(supp_table_1_all_cohorts[
  supp_table_1_all_cohorts[,ncol(supp_table_1_all_cohorts)]!="","GSE"]))

supp_table_genes = c()
for(nn in names(analysis2selected_genes_stats)){
  m = analysis2selected_genes_stats[[nn]]
  m = cbind(rep(nn,nrow(m)),m)
  colnames(m)[1] = "Discovered in"
  currgenes = rownames(m)
  subgroups = c()
  for(g in currgenes){
    curr_subgroups = names(gene_subgroups)[sapply(gene_subgroups,function(x,y)is.element(y,set=x),y=g)]
    curr_subgroups = paste(curr_subgroups,collapse=" and ")
    subgroups[g]=curr_subgroups
  }
  m = cbind(m,subgroups)
  supp_table_genes = rbind(supp_table_genes,m)
}
rownames(supp_table_genes)=NULL
write.xlsx(supp_table_genes,file=supp_file,sheetName = "STable2",row.names = F,append = T)

# Mean t-statistic matrices
base_model_ms = gene_t_patterns$base_models
base_model_ms = cbind(rownames(base_model_ms),
                      unlist(entrez2symbol[rownames(base_model_ms)]),base_model_ms)
colnames(base_model_ms)[1:3] = c("Entrez","Symbol","Analysis_type")
rownames(base_model_ms)=NULL
write.xlsx(base_model_ms,file=supp_file,sheetName = "STable3_base_models_ts",
           row.names = F,append = T)
sheet_counter=3
for(gg in names(gene_t_patterns)[!grepl("base_",names(gene_t_patterns))]){
  sheet_counter = sheet_counter+1
  ms = gene_t_patterns[[gg]]
  ms = cbind(rownames(ms),unlist(entrez2symbol[rownames(ms)]),ms)
  colnames(ms)[1:2] = c("Entrez","Symbol")
  gg = gsub(",",replacement = "_",gg)
  currname=paste("tstats_",gg,sep="")
  write.xlsx(ms,file=supp_file,
             sheetName = paste("STable",sheet_counter,"_",currname,sep=""),
             row.names = F,append = T)
}

# GO enrichments of subgroups
sheet_counter = sheet_counter+1
supp_table_enrichments = gene_subgroup_enrichments_fdr[,c(1:4,9:10)]
colnames(supp_table_enrichments)[6] = "q-value"
colnames(supp_table_enrichments)[5] = "Genes"
colnames(supp_table_enrichments)[1] = "Discovered in"
write.xlsx(supp_table_enrichments,file=supp_file,
           sheetName = paste("STable",sheet_counter,"_GO_enrichments",sep=""),
           row.names = F,append=T)

# Same for Reactome
sheet_counter = sheet_counter+1
supp_table_enrichments = reactome_pathways_subgroups_fdr[,c(1,3,8)]
colnames(supp_table_enrichments)[3] = "q-value"
colnames(supp_table_enrichments)[1] = "Discovered in"
write.xlsx(supp_table_enrichments,file=supp_file,
           sheetName = paste("STable",sheet_counter,"_Reactome_enrichments",sep=""),
           row.names = F,append=T)

# # Original groups (without kmeans and partitions by meta-regression type)
# # Enrichment tables: go or reactome and group or subgroup
# supp_table_enrichments = gene_group_enrichments_fdr[,c(1:4,9:10)]
# colnames(supp_table_enrichments)[6] = "q-value"
# colnames(supp_table_enrichments)[5] = "Genes"
# colnames(supp_table_enrichments)[1] = "Discovered in"
# write.table(supp_table_enrichments,
#             file=paste(supp_path,"supp_table_enrichments_go_group.txt",sep=""),
#             sep="\t",quote=F,col.names = T,row.names = F)
# 
# supp_table_enrichments = reactome_pathways_groups_fdr[,c(1,3,8)]
# colnames(supp_table_enrichments)[3] = "q-value"
# colnames(supp_table_enrichments)[1] = "Discovered in"
# write.table(supp_table_enrichments,
#             file=paste(supp_path,"supp_table_enrichments_reactome_group.txt",sep=""),
#             sep="\t",quote=F,col.names = T,row.names = F)


