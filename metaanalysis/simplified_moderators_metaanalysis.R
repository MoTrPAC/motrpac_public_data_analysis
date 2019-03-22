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

# # Before the update in March 2019
# setwd('/Users/David/Desktop/MoTrPAC/PA_database')
# # Prepare the datasets for the different analyses below
# # Load the datasets and their metadata
# load("PADB_univariate_results_and_preprocessed_data_acute.RData")
# acute_datasets = cohort_data
# acute_metadata = cohort_metadata
# acute_sample2time = sample2time
# load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
# longterm_datasets = cohort_data
# longterm_metadata = cohort_metadata
# longterm_sample2time = sample2time
# load("PADB_dataset_level_meta_analysis_data.RData")

# After the update
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

# Preprocessing before the data update on March 2019
# # Clean the datasets: from muscle, remove 1 cohort with time point zero
# # From long term exclude all studies with time > 150 days
# clean_acute_table <-function(gdata,tissue="muscle",remove_untrained=T){
#   if(remove_untrained){
#     gdata = gdata[!is.element(gdata$training,set=c("yoga","control","untrained")),]
#   }
#   gdata = gdata[gdata$tissue==tissue,]
#   gdata$vi = pmax(gdata$vi,1e-5)
#   if(tissue=="muscle"){
#     gdata = gdata[as.numeric(gdata$time) > 1,]
#   }
#   gdata = simplify_time_in_gdata(gdata,simplify_time_acute)
#   gdata = simplify_age_gdata(gdata,thr=40)
#   gdata = simplify_training_for_exercise_analysis(gdata)
#   gdata$prop_males = as.numeric(gdata$prop_males)
#   gdata$p = as.numeric(gdata$p)
#   gdata$sdd = as.numeric(gdata$sdd)
#   gdata$df = as.numeric(gdata$df)
#   gdata$tstat = as.numeric(gdata$tstat)
#   return(gdata)
# }
# clean_longterm_table <-function(gdata,tissue="muscle",remove_untrained=T){
#   if(remove_untrained){
#     gdata = gdata[!is.element(gdata$training,set=c("yoga","control","untrained")),]
#   }
#   gdata = gdata[gdata$tissue==tissue,]
#   gdata$vi = pmax(gdata$vi,1e-5)
#   gdata = gdata[as.numeric(gdata$time) < 150 ,]
#   gdata = simplify_age_gdata(gdata,thr=40)
#   gdata = simplify_training_for_exercise_analysis(gdata)
#   gdata$prop_males = as.numeric(gdata$prop_males)
#   gdata$p = as.numeric(gdata$p)
#   gdata$sdd = as.numeric(gdata$sdd)
#   gdata$df = as.numeric(gdata$df)
#   gdata$tstat = as.numeric(gdata$tstat)
#   return(gdata)
# }

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

# # Reshape data for replication analysis: before the update
# rep_datasets = lapply(datasets,function(x)t(sapply(x,function(y)as.numeric(y$p))))
# for(nn in names(rep_datasets)){
#   colnames(rep_datasets[[nn]]) = 
#     paste(rownames(datasets[[nn]][[1]]),datasets[[nn]][[1]]$V1,
#           datasets[[nn]][[1]]$time,sep=";")
# }
# sapply(rep_datasets,dim)
# sapply(rep_datasets,rownames)
# sapply(rep_datasets,colnames)

# Reshape data for replication analysis: after the update
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

# # Get the untrained controls datasets : before the update 
# untrained_datasets = list()
# untrained_datasets[["acute,muscle"]] = lapply(acute_gene_tables_raw,get_untrained_table,tissue="muscle")
# untrained_datasets[["acute,blood"]] = lapply(acute_gene_tables_raw,get_untrained_table,tissue="blood")
# untrained_datasets[["longterm,muscle"]] = lapply(longterm_gene_tables_raw,get_untrained_table,tissue="muscle")
# untrained_datasets[["longterm,blood"]] = lapply(longterm_gene_tables_raw,get_untrained_table,tissue="blood")
# sapply(untrained_datasets,function(x)dim(x[[1]]))

# Get the untrained controls datasets : after the update 
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
# for(nn in names(simple_RE_pvals)){
#   pvalue_qqplot(simple_RE_pvals[[nn]],main=nn,pch=20,cex=0.5)
# }
for(nn in names(simple_RE_pvals)){
  r = hist(simple_RE_I2s[[nn]],plot=F)
  per = sum(simple_RE_I2s[[nn]]>70,na.rm = T)/length(simple_RE_I2s[[nn]])
  per = format(per*100,digits=3)
  cols=rep("white",length(r$mids))
  cols[r$mids > 70] = "blue"
  hist(simple_RE_I2s[[nn]],main=paste(nn,"(",per,"%)",sep=""),xlab = "I^2(%)",col=cols)
}
# for(nn in names(simple_RE_pvals)){
#   pvalue_qqplot(unlist(st_pvals[[nn]]),main=nn,pch=20,cex=0.5)
# }
# for(nn in names(simple_RE_pvals)){
#   x1 = simple_RE_I2s[[nn]]
#   x2 = log(simple_RE_pvals[[nn]])
#   inds = !is.na(x1) & !is.na(x2)
#   print(nn)
#   print(cor(x1[inds],x2[inds],method="spearman"))
# }
# par(mfrow=c(2,2))
# for(nn in names(simple_RE_pvals)){
#   x1 = simple_RE_I2s[[nn]]
#   x2 = simple_RE_tau2s[[nn]]
#   inds = !is.na(x1) & !is.na(x2)
#   print(nn)
#   print(cor(x1[inds],x2[inds],method="spearman"))
#   plot(x=x1,y=x2,xlab="I2(%)",ylab="Tau^2",pch=20,cex=0.2,main=nn)
# }
# par(mfrow=c(2,2))
# for(nn in names(simple_RE_pvals)){
#   x1 = simple_RE_I2s[[nn]]
#   x2 = simple_RE_beta[[nn]]
#   inds = !is.na(x1) & !is.na(x2)
#   print(nn)
#   print(cor(x1[inds],x2[inds],method="spearman"))
#   plot(x=x1,y=x2,xlab="I2(%)",ylab="beta",pch=20,cex=0.2,main=nn)
# }

# # Example for high I2 and low tau
# # TODO: revise and add some stats to the report
# nn  = 2
# table(simple_RE_I2s[[nn]] < 50 & abs(simple_RE_beta[[nn]]) < 0.25)
# selected_is = names(which(simple_RE_I2s[[nn]] > 60 &
#     simple_RE_tau2s[[nn]] < 0.1 & abs(simple_RE_beta[[nn]]) > 0.25))
# selected_is = names(which(simple_RE_I2s[[nn]] > 90 &
#     simple_RE_tau2s[[nn]] < 0.01))
# selected_is = names(which(simple_RE_I2s[[nn]] < 25 &
#     simple_RE_pvals[[nn]] < 0.00001 & abs(simple_RE_beta[[nn]]) > 0.25))
# selected_is = names(which(
#     simple_RE_pvals[[nn]] < 1e-8 & abs(simple_RE_beta[[nn]]) > 0.25))
# selected_i = sample(selected_is)[1]
# forest(simple_REs[[nn]][[selected_i]])
# median(as.numeric(datasets[[nn]][[selected_i]]$p))
# simple_REs[[nn]][[selected_i]]

# Filters 1 exclude genes with no sig p at 0.01
# before the update on March 2019 the threshold was 0.05
rep_filter <-function(gdata,num=1,thr=0.01){
  return(sum(gdata$p<=thr,na.rm = T)>=num)
}
to_rem1 = sapply(datasets,function(x)!sapply(x,rep_filter))
sapply(to_rem1,table)
rm(acute_gene_tables_raw)
rm(acute_gene_tables)
rm(longterm_gene_tables_raw)
rm(longterm_gene_tables)
save.image(file="workspace_before_rep_analysis.RData")

# Create the input files for meta-regression and replication analysis
meta_reg_datasets = list()
for(nn in names(simple_RE_beta)){
  curr_genes = names(which(!to_rem1[[nn]]))
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
library(parallel);library(metafor)

I2_thr = 50
AIC_diff_thr = 5
ACUTE_beta_thr = 0.25
LONGTERM_beta_thr = 0.1
P_thr = 0.0001

# Load the results (run on sherlock) instead of running the code above
load("meta_analysis_results.RData")
load("workspace_before_rep_analysis.RData")
load("meta_analysis_input.RData")

# Algorithm for selecting genes from each meta-reg analysis
analysis2selected_genes = list()
analysis2selected_genes_stats = list()
for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]]
  pvals = sapply(analysis1,function(x)x[[1]]$mod_p)
  i2s = simple_RE_I2s[[nn]][names(pvals)]
  aic_diffs = unlist(sapply(analysis1,function(x)unname(x[[1]]$aic_c - x$`simple:base_model`$aic_c)))
  aic_diffs[setdiff(names(i2s),names(aic_diffs))] = 0
  aic_diffs = aic_diffs[names(i2s)]
  if(grepl("acute",nn)){
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>ACUTE_beta_thr))
  }
  else{
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>LONGTERM_beta_thr))  
  }
  
  # separate into two gene sets: those that passed the aic diff test vs. those that did not
  genes_with_high_aic_diff = names(aic_diffs)[aic_diffs <= -AIC_diff_thr]
  selected_aic_diff_genes = names(aic_diffs)[aic_diffs <= -AIC_diff_thr & model2beta & pvals <= P_thr]
  genes_without_high_aic_diff = names(aic_diffs)[aic_diffs > -AIC_diff_thr]
  genes_without_high_aic_diff_i2s = i2s[genes_without_high_aic_diff]
  genes_without_high_aic_diff_pvals = simple_RE_pvals[[nn]][genes_without_high_aic_diff]
  genes_without_high_aic_diff_betas = abs(simple_RE_beta[[nn]][genes_without_high_aic_diff])
  if(grepl("acute",nn)){
    genes_without_high_aic_diff_betas = genes_without_high_aic_diff_betas>ACUTE_beta_thr
  }
  else{
    genes_without_high_aic_diff_betas = genes_without_high_aic_diff_betas>LONGTERM_beta_thr
  }
  selected_base_model_genes = genes_without_high_aic_diff[
    genes_without_high_aic_diff_i2s < I2_thr &
      genes_without_high_aic_diff_pvals < P_thr &
      genes_without_high_aic_diff_betas
    ]
  selected_base_model_genes = selected_base_model_genes[!is.na(selected_base_model_genes)]
  
  curr_selected_genes = union(selected_aic_diff_genes,selected_base_model_genes)
  curr_selected_genes_names = sapply(analysis1[selected_aic_diff_genes],function(x)names(x)[1])
  curr_selected_genes_names = sapply(curr_selected_genes_names,function(x)strsplit(x,split=":")[[1]][2])
  curr_selected_genes_names[selected_base_model_genes] = "base_model"
  # is.element("10891",set=gene_set2)
  # curr_selected_genes_names["10891"]
  # table(curr_selected_genes_names)
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
  # # Correct the table above: base model should have base p-values, 
  # # and a comment about the AIC difference
  # base_model_genes = is.element(m[,1],set=selected_base_model_genes)
  # if(sum(base_model_genes)>0){
  #   m[base_model_genes,4] = pvals_base[m[base_model_genes,1]]
  #   m[base_model_genes,5] = paste(m[base_model_genes,5],"*",sep="")
  # }
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
par(mfrow=c(1,3),mar=c(10,4,2,2),cex.lab=1.2,cex.axis=1.2,cex=1,cex.main=0.9)
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
gs = lapply(analysis2selected_genes,names)
go_res = run_topgo_enrichment_fisher(
  gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
go_res1 = go_res[go_res$Annotated < 1500,]
go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
table(go_res_fdr$setname)
gene_group_enrichments = go_res
gene_group_enrichments_fdr = go_res_fdr

gs = lapply(analysis2selected_genes,names)
gs = lapply(gs,function(x,y)y[x],y=unlist(entrez2symbol))

# par(mfrow=c(2,2))
# for(nn in names(all_meta_analysis_res)){
#   analysis1 = all_meta_analysis_res[[nn]]
#   model_names = sapply(analysis1,function(x)names(x)[1])
#   table(model_names)
#   ps1 = sapply(analysis1,function(x)x[[1]]$mod_p)
#   i2s = sapply(analysis1,function(x)x[[1]]$I2)
#   i2s_2 = simple_RE_I2s[[nn]][names(i2s)]
#   aic_diffs = sapply(analysis1,function(x)x[[1]]$aic_c - x[[2]]$aic_c)
#   i2_diffs = i2s-i2s_2
#   # pvalue_qqplot(ps1);abline(0,1)
#   plot(ps1,i2s)
#   plot(i2s_2,i2s)
#   plot(aic_diffs,i2_diffs)
#   print(length(ps1))
#   print(table(i2s<50 & ps1 < 0.0001))
# }

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
dev.off()
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)
boxplot(yi~time,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",names=c("early","late"),ylab="Fold change")
# GJA1
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
# MYLK2
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
# Time and training type - CNRIP1
gene ="25927"
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
boxplot(yi~time+training,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",ylab="Fold change",las=1)
# 7423 VEGFB
gene ="7423"
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,blood`[[gene]]
gdata = meta_reg_datasets$`acute,blood`[[gene]]
curr_m$slab.null = F
curr_times = rep("early",nrow(gdata))
curr_times[gdata$time==24] = "late(>12h)"
curr_times[gdata$time==4] = "inter(<5h)"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)
boxplot(yi~time+training,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",ylab="Fold change",las=1)

# Males-dependent example
gene = "7188"
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,blood`[[gene]]
gdata = meta_reg_datasets$`acute,blood`[[gene]]
curr_m$slab.null = F
curr_times = rep("early",nrow(gdata))
curr_times[gdata$time==24] = "late(>12h)"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,blood`
aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)
boxplot(yi~time,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",ylab="Fold change")
table(gdata$time,gdata$prop_males)
boxplot(yi~prop_males,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",ylab="Fold change")

# known_genes = c("PPARGC1A","COX1","NDUFA","PDK4","VEGFA","KDR","THY1","MYL4",
#                 "MYH1","COL1A1","ACTC1","TNNT2","GADD45G","MMP9","NR4A1")


############################################################################
############################################################################
############################################################################
# Interpretation of the results
library(topGO)

# Get effect matrices - mean responses - t statistics
mean_effect_matrices = lapply(datasets,function(x)t(sapply(x,function(y)as.numeric(y$tstat))))
for(nn in names(mean_effect_matrices)){
  colnames(mean_effect_matrices[[nn]]) = 
    paste(rownames(datasets[[nn]][[1]]),datasets[[nn]][[1]]$training,
          datasets[[nn]][[1]]$time,sep=";")
  colnames(mean_effect_matrices[[nn]]) = gsub("endurance","E",
                                              colnames(mean_effect_matrices[[nn]]),ignore.case = T)
  colnames(mean_effect_matrices[[nn]]) = gsub("resistance","R",
                                              colnames(mean_effect_matrices[[nn]]),ignore.case = T)
}

gene_subgroups = list()
gene_t_patterns = list()
for(nn in names(analysis2selected_genes_stats)){
  curr_genes = analysis2selected_genes_stats[[nn]]
  curr_groups = curr_genes[,"Group"]
  for(gg in unique(curr_groups)){
    curr_m = as.matrix(curr_genes[gg==curr_genes[,"Group"],])
    if(ncol(curr_m)==1){curr_m=t(curr_m)}
    m = as.matrix(mean_effect_matrices[[nn]][curr_m[,"Entrez"],])
    if(ncol(m)==1){m=t(m)}
    rownames(m) = curr_m[,"Entrez"]
    
    # cluster based on the moderator's average not the t-test scores themselves
    if(gg != "base_model"){
      curr_mods = strsplit(gg,split=";")[[1]]
      curr_gdata = datasets[[nn]][[1]]
      for(j in names(curr_gdata)){
        if(is.numeric(curr_gdata[[j]])){curr_gdata[[j]] = format(curr_gdata[[j]],digits=2)}
      }
      curr_mods = curr_gdata[,curr_mods]
      if(!is.null(dim(curr_mods))){
        curr_mods = apply(curr_mods,1,paste,collapse=";")
      }
      m = t(apply(m,1,function(x,y)tapply(x,y,mean),y=curr_mods))
    }
    else{
      m = matrix(rowMeans(m),ncol=1)
      rownames(m) = curr_m[,"Entrez"]
      colnames(m) = "base_model_mean_t"
    }
    
    currk = 2
    # if(nrow(curr_m)>50){currk=4}
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
    for(j in 1:ncol(m)){m[,j]=as.numeric(format(m[,j],digits=3))}
    gene_t_patterns[[paste(nn,gg,sep=",")]] = m
    
    for(kk in unique(m_kmeans)){
      currname = paste(nn,gg,kk,sep=",")
      gene_subgroups[[currname]] = names(m_kmeans)[m_kmeans==kk]
    }
  }
}
sapply(gene_subgroups,length)
gene_subgroups[["longterm,muscle,prop_males,1"]]
sapply(gene_t_patterns,colnames)
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
go_res = run_topgo_enrichment_fisher(
  gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
go_res1 = go_res[go_res$Annotated < 1500,]
go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
table(go_res_fdr$setname)
gene_subgroup_enrichments = go_res
gene_subgroup_enrichments_fdr = go_res_fdr

get_most_sig_enrichments_by_groups(gene_subgroup_enrichments_fdr,num=2)
get_most_sig_enrichments_by_groups(gene_group_enrichments_fdr,num=3)

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
get_most_sig_enrichments_by_groups(reactome_pathways_subgroups_fdr,pcol="pvalue",num = 2)

save(gene_subgroup_enrichments,gene_subgroup_enrichments_fdr,
     gene_group_enrichments,gene_group_enrichments_fdr,
     reactome_pathways_subgroups,reactome_pathways_subgroups_fdr,
     reactome_pathways_groups,reactome_pathways_groups_fdr,
     file="topGO_res_jan_2019.RData")


# Heatmaps
library(gplots)
gene_set = gene_subgroups$`acute,muscle,training,2`
ord = order(datasets$`acute,muscle`[[1]]$training,
            datasets$`acute,muscle`[[1]]$time)
mat = mean_effect_matrices$`acute,muscle`[gene_set,ord]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat[mat>5]=5;mat[mat< -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.9)

gene_set = gene_subgroups$`acute,muscle,training,1`
ord = order(datasets$`acute,muscle`[[1]]$training,
            datasets$`acute,muscle`[[1]]$time)
mat = mean_effect_matrices$`acute,muscle`[gene_set,ord]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat[mat>5]=5;mat[mat< -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.9)

gene_set = gene_subgroups$`acute,muscle,time,1`
ord = order(datasets$`acute,muscle`[[1]]$time,
            datasets$`acute,muscle`[[1]]$training)
mat = mean_effect_matrices$`acute,muscle`[gene_set,ord]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat[mat>5]=5;mat[mat< -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered)

gene_set = gene_subgroups$`longterm,muscle,base_model,2`
ord = order(datasets$`longterm,muscle`[[1]]$time,
            datasets$`longterm,muscle`[[1]]$training)
mat = mean_effect_matrices$`longterm,muscle`[gene_set,ord]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat[mat>5]=5;mat[mat< -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.7)

gene_set = gene_subgroups$`longterm,muscle,base_model,1`
ord = order(datasets$`longterm,muscle`[[1]]$time,
            datasets$`longterm,muscle`[[1]]$training)
mat = mean_effect_matrices$`longterm,muscle`[gene_set,ord]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat[mat>5]=5;mat[mat< -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.7)

gene_set = gene_subgroups$`acute,blood,prop_males,1`
ord = order(datasets$`acute,blood`[[1]]$prop_males,
            datasets$`acute,blood`[[1]]$time)
mat = mean_effect_matrices$`acute,blood`[gene_set,ord]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat[mat>5]=5;mat[mat< -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.7)

gene_set = gene_subgroups$`acute,blood,time,1`
ord = order(datasets$`acute,blood`[[1]]$time)
mat = mean_effect_matrices$`acute,blood`[gene_set,ord]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat[mat>5]=5;mat[mat < -5]=-5
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered)

# Analysis of the acute blood datasets
plot_with_err_bars<-function(xnames,avg,sdev,add=F,...){
  if(add){
    lines(avg,pch=19,...)
  }
  else{
    plot(avg,xaxt = "n",pch=19, type='l',...)
    axis(1, at=1:length(xnames), labels=xnames)
  }
  # hack: we draw arrows but with very special "arrowheads"
  arrows(1:length(xnames), avg-sdev, 1:length(xnames), avg+sdev,
         length=0.05, angle=90, code=3)
}

par(mfrow=c(1,2))
gene_set = gene_subgroups$`acute,blood,time,1`
ord = order(datasets$`acute,blood`[[1]]$time)
mat = mean_effect_matrices$`acute,blood`[gene_set,ord]
times = ordered(as.numeric(sapply(colnames(mat),function(x)strsplit(x,split=";")[[1]][3])))
D = dummy(times,levelsToKeep = levels(times))
D = t(t(D)/colSums(D))
mat1 = mat %*% D
boxplot(mat1,xlab="Time after bout (hours)",ylab="Avg t-statistic",col="blue",
        main="Down-regulated (196 genes)")
gene_set = gene_subgroups$`acute,blood,time,2`
ord = order(datasets$`acute,blood`[[1]]$time)
mat = mean_effect_matrices$`acute,blood`[gene_set,ord]
times = ordered(as.numeric(sapply(colnames(mat),function(x)strsplit(x,split=";")[[1]][3])))
D = dummy(times,levelsToKeep = levels(times))
D = t(t(D)/colSums(D))
mat1 = mat %*% D
boxplot(mat1,xlab="Time after bout (hours)",ylab="Avg t-statistic",col ="red",
       main= "Up-regulated (99 genes)")

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
          "Time_points"=tps)
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


