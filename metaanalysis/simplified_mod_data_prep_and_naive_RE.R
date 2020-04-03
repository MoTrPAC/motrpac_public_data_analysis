library(org.Hs.eg.db);library(metafor)
source('~/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)

# setwd('~/Desktop/MoTrPAC/project_release_feb_2018/data/')
setwd('~/Desktop/MoTrPAC/project_release_feb_2018/revision_feb_2020/')

# Prepare the datasets for the different analyses below
# Load the datasets and their metadata
load("human_ge_cohort_preprocessed_db_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
acute_sample2time = sample2time
acute_sample_meta = sample_metadata
sample2sex_a = sample2sex
load("human_ge_cohort_preprocessed_db_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata
longterm_sample2time = sample2time
acute_sample2time = sample2time
longterm_sample_meta = sample_metadata
sample2sex_l = sample2sex
load("human_ge_cohort_preprocessed_db_gene_tables.RData")

# set the output dir for RData files
# out_dir = "" # for the WD
# out_dir = "~/Desktop/MoTrPAC/project_release_feb_2018/data/"
out_dir = "~/Desktop/MoTrPAC/project_release_feb_2018/revision_feb_2020/"
out_dir_figs = paste(out_dir,"figures/",sep="")
try(dir.create(out_dir_figs))

############################################################################
############################################################################
############################################################################
# Get some stats for the main paper

# GSM-based table: subjects, tissues, sex
gses_to_remove = c("GSE8668","GSE18966","GSE83578","GSE28498","GSE41914","GSE43856","GSE51835")
complete_sample_table = c()
complete_sample_table_excluded_from_meta_analysis = c()
for(ac in names(acute_metadata)){
  if(acute_metadata[[ac]]$gse %in% gses_to_remove){next}
  if(acute_metadata[[ac]]$tissue == "fat"){next}
  curr_gsms = acute_metadata[[ac]]$gsms
  curr_subjects = acute_sample_meta$subject[curr_gsms]
  curr_subjects = paste(acute_metadata[[ac]]$gse,curr_subjects,sep=";")
  curr_sex = acute_sample_meta$sex[curr_gsms]
  curr_tissue = rep(acute_metadata[[ac]]$tissue,length(curr_gsms))
  curr_training = rep(acute_metadata[[ac]]$training,length(curr_gsms))
  curr_type = rep("acute",length(curr_gsms))
  m = cbind(curr_gsms,curr_subjects,curr_sex,curr_tissue,curr_training,curr_type)
  if(length(unique(acute_metadata[[ac]]$times))<2){
    complete_sample_table_excluded_from_meta_analysis = rbind(
      complete_sample_table_excluded_from_meta_analysis, m
    )
  }
  else{
    complete_sample_table = rbind(complete_sample_table,m)
  }
}
for(lo in names(longterm_metadata)){
  if(longterm_metadata[[lo]]$gse %in% gses_to_remove){next}
  if(longterm_metadata[[lo]]$tissue == "fat"){next}
  curr_gsms = longterm_metadata[[lo]]$gsms
  curr_subjects = longterm_sample_meta$subject[curr_gsms]
  curr_subjects = paste(longterm_metadata[[lo]]$gse,curr_subjects,sep=";")
  curr_sex = longterm_sample_meta$sex[curr_gsms]
  curr_tissue = rep(longterm_metadata[[lo]]$tissue,length(curr_gsms))
  curr_training = rep(longterm_metadata[[lo]]$training,length(curr_gsms))
  curr_type = rep("longterm",length(curr_gsms))
  m = cbind(curr_gsms,curr_subjects,curr_sex,curr_tissue,curr_training,curr_type)
  if(length(unique(longterm_metadata[[lo]]$times))<2){
    complete_sample_table_excluded_from_meta_analysis = rbind(
      complete_sample_table_excluded_from_meta_analysis, m
    )
  }
  else{
    complete_sample_table = rbind(complete_sample_table,m)
  }
}
dim(complete_sample_table)
print("sample based counts:")
print("table of NAs in sex:")
na_inds = is.na(complete_sample_table[,3])
unique(complete_sample_table[na_inds,4:6])
print(table(na_inds,complete_sample_table[,"curr_type"]))
print("tissue vs. study type:")
print(table(complete_sample_table[,"curr_tissue"],complete_sample_table[,"curr_type"]))
print("tissue counts:")
tissues_count_info = unique(complete_sample_table[,c(1,4)])
print(table(tissues_count_info[,2]))
# should fit the numbers in Supp Table 1
dim(complete_sample_table[complete_sample_table[,4]!="fat",])
# This should fit the count of non-fat excluded datasets (362 as of May 2019):
dim(complete_sample_table_excluded_from_meta_analysis[
  complete_sample_table_excluded_from_meta_analysis[,4]!="fat",
  ])

# Read in the sample-level metadata and get the counts
metadata_file = 'GEO_sample_metadata.xlsx'
complete_sample_table_metaanalysis = complete_sample_table[
  !(rownames(complete_sample_table) %in% 
      rownames(complete_sample_table_excluded_from_meta_analysis)),
]
dim(complete_sample_table_metaanalysis)

# Read the sample metadata - these are the sample sets to be analyzed
# The loaded RData objects contain more datasets than what is used for the
# different meta-analyses
library(readxl)
acute_metadata_raw = data.frame(read_xlsx(metadata_file,sheet=1))
acute_metadata_raw[is.na(acute_metadata_raw)] = ""
longterm_metadata_raw = data.frame(read_xlsx(metadata_file,sheet=2))
longterm_metadata_raw[is.na(longterm_metadata_raw)] = ""

# Compare the sex imputation results
added_names = setdiff(names(sample2sex_a),names(sample2sex_l))
sample2sex = sample2sex_l
sample2sex[added_names] = sample2sex_a[added_names]
length(sample2sex)
table(sample2sex)
length(unique(names(sample2sex)))
table(sample2sex)
table(sample2sex[rownames(complete_sample_table)],complete_sample_table[,3])
complete_sample_table[,3] = sample2sex[complete_sample_table[,1]]

rawsex = acute_metadata_raw$Gender
names(rawsex) = acute_metadata_raw$GSM
rawsex[longterm_metadata_raw$GSM] = longterm_metadata_raw$Gender
length(rawsex)

# compare the names
length(setdiff(names(rawsex),names(sample2sex)))
length(setdiff(names(rawsex),names(rawsex)))

rawsex = rawsex[names(sample2sex)]
table(sample2sex,rawsex)
table(sample2sex[complete_sample_table[,1]])
table(is.na(sample2sex[complete_sample_table[,1]]))
table(rawsex[complete_sample_table[,1]])
length(sample2sex)
dim(complete_sample_table)
table(complete_sample_table[,4])
length(unique(complete_sample_table[,2]))
imputed_sex = names(rawsex)[rawsex==""]
imputed_sex = imputed_sex[!is.na(sample2sex[imputed_sex])]
imputed_sex = intersect(imputed_sex,rownames(complete_sample_table))
table(unique(complete_sample_table[,2:4])[,2])
table(unique(complete_sample_table[imputed_sex,2:4])[,2])
table(is.na(unique(complete_sample_table[,2:4])[,2]))
dim(
  complete_sample_table[setdiff(rownames(complete_sample_table),imputed_sex),])
table(complete_sample_table[setdiff(rownames(complete_sample_table),imputed_sex),3])

# sample_level_meta = rbind(acute_metadata[,c("GSM","GSE","Subject.id","Tissue","Gender","Numeric_Age")],
#                           longterm_metadata[,c("GSM","GSE","Subject.id","Tissue","Gender","Numeric_Age")])
# dim(sample_level_meta)
# sample_level_meta_in_metaanalysis = sample_level_meta[
#   is.element(sample_level_meta[,1],
#              set=complete_sample_table[,1]),]
# dim(sample_level_meta_in_metaanalysis)
# subject_ids = paste(sample_level_meta_in_metaanalysis[,"GSE"],
#                     sample_level_meta_in_metaanalysis[,"Subject.id"])
# sample_level_meta_in_metaanalysis = cbind(sample_level_meta_in_metaanalysis,subject_ids)
# length(unique(subject_ids))
# sex_info = unique(sample_level_meta_in_metaanalysis[,c("subject_ids","Gender")])
# table(sex_info$Gender!="")
# table(grepl("f",sex_info$Gender,ignore.case = T))
# acute_metadata[1:5,1:5]

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
clean_acute_table <-function(gdata,tissue="muscle",remove_untrained=T,gses_to_remove = NULL){
  if(remove_untrained){
    gdata = gdata[!is.element(gdata$training,set=c("yoga","control","untrained")),]
  }
  gdata = gdata[gdata$tissue==tissue,]
  if(!is.null(gses_to_remove)){
    gdata = gdata[! gdata$gse %in% gses_to_remove,]
  }
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
# Binarize time based on 150 days
clean_longterm_table <-function(gdata,tissue="muscle",remove_untrained=T,gses_to_remove=NULL){
  if(remove_untrained){
    gdata = gdata[!is.element(gdata$training,set=c("yoga","control","untrained")),]
  }
  gdata = gdata[gdata$tissue==tissue,]
  if(!is.null(gses_to_remove)){
    gdata = gdata[!gdata$gse %in% gses_to_remove,]
  }
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
datasets[["acute,blood"]] = lapply(acute_gene_tables,clean_acute_table,tissue="blood",
                                   gses_to_remove = c("GSE8668","GSE18966","GSE83578","GSE28498","GSE41914","GSE43856","GSE51835"))
datasets[["longterm,muscle"]] = lapply(longterm_gene_tables,clean_longterm_table,tissue="muscle")
datasets[["longterm,blood"]] = lapply(longterm_gene_tables,clean_longterm_table,tissue="blood")

# A new filter added after the update on March 2019: exclude
# genes with extremely low number of studies from each analysis
# par(mfrow=c(2,2))
for(nn in names(datasets)){
  num_datasets = sapply(datasets[[nn]],function(x)length(unique(x$gse)))
  # hist(num_datasets,main=nn)
  to_rem = (num_datasets/max(num_datasets)) < 0.75
  print(table(to_rem))
  datasets[[nn]] = datasets[[nn]][!to_rem]
}

# Reshape data for replication analysis
get_gene_data_for_rep_analysis<-function(gdata,exclude){
  gdata = gdata[!is.element(gdata$V1,set=exclude),]
  x = gdata$p
  names(x)=paste(gdata$V1,gdata$time,sep=";")
  return(x)
}
load("~/Desktop/MoTrPAC/project_release_feb_2018/data/human_ge_gene_coverage_analysis.RData")
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
# write down the p-value matrices
for(nn in names(rep_datasets)){
  pvals_m = rep_datasets[[nn]]
  pvals_m[is.na(pvals_m)] = 0.5
  nn_file = paste(
    gsub(",","_",nn),"_pvals_matrix.txt",sep=""
  )
  write.table(pvals_m,file=nn_file,sep="\t",row.names = T,col.names = T,quote = F)
}


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
simple_REs_untrained_beta = lapply(simple_REs_untrained,function(x)sapply(x,try_get_field,fname="beta"))

pdf(paste0(out_dir_figs,"Supp_Figure1B.pdf"))
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  pvalue_qqplot(simple_RE_pvals[[nn]],n_random = 10000,
                main=nn,pch=20,cex=0.5,
                cex.main=1.7,cex.lab=1.4)
}
dev.off()

pdf(paste0(out_dir_figs,"Supp_Figure1A.pdf"))
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  r = hist(simple_RE_I2s[[nn]],plot=F)
  per = sum(simple_RE_I2s[[nn]]>70,na.rm = T)/length(simple_RE_I2s[[nn]])
  per = format(per*100,digits=3)
  cols=rep("white",length(r$mids))
  cols[r$mids > 70] = "blue"
  hist(simple_RE_I2s[[nn]],main=paste(nn,"(",per,"%)",sep=""),
       xlab = "I^2(%)",col=cols,cex.main=1.7,cex.lab=1.4)
}
dev.off()

# Filters 1 exclude genes with no sig p at 0.01
# before the update on March 2019 the threshold was 0.05 and 1
# Returns true if the gene has p<thr in at least num studies
rep_filter <-function(gdata,num=1,thr=0.01){
  return(sum(gdata$p<=thr,na.rm = T)>=num)
}
to_rem1 = sapply(datasets,function(x)!sapply(x,rep_filter,num=2,thr=0.05))

rm(acute_gene_tables);rm(longterm_gene_tables)
save.image(file=paste(out_dir,"workspace_before_rep_analysis.RData",sep=""))

# Create the input files for meta-regression and replication analysis
meta_reg_datasets = list()
for(nn in names(simple_RE_beta)){
  curr_genes = names(which(!to_rem1[[nn]]))
  meta_reg_datasets[[nn]] = datasets[[nn]][curr_genes]
}
sapply(meta_reg_datasets,length)

# Original moderators
# meta_reg_to_mods = list(
#   "acute,muscle" = c("time","training","avg_age"),
#   "acute,blood" = c("time","prop_males","avg_age"),
#   "longterm,muscle" = c("training","time","avg_age","prop_males"),
#   "longterm,blood" = "training"
# )

# Revision on Feb 2020
# Not enough non-male acute muscle studies
# No resistance training and enough time points for acute, blood
# For longterm blood we take training only as we assume that it will
# be the primary cause of differential abundance
meta_reg_to_mods = list(
  "acute,muscle" = c("time","training","avg_age"),
  "acute,blood" = c("prop_males","avg_age"),
  "longterm,muscle" = c("training","time","avg_age","prop_males"),
  "longterm,blood" = "training"
)
save(meta_reg_datasets,meta_reg_to_mods,
     untrained_datasets,file=paste(out_dir,"meta_analysis_input.RData",sep=""))

############################################################################
############################################################################
############################################################################
############################################################################