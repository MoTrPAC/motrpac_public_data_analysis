# This script loads the raw profiles and the metadata and
# creates standardized datasets and metadata. The preprocessed 
# database is then saved in an easy to use objects. These are kept in
# RData files, listed below.
# The process below takes care of:
# 1. Getting the gene data matrix per study
# 2. Partition of studies to cohorts
# 3. Sex imputation
# 4. Analysis of platform gene overlap
# 5. Computing summary statistics per time point
# 6. Transforming the data to gene tables for the meta-analysis

###############################################
###############################################
############### Config workspace ##############
###############################################
###############################################
# Paths
WD = "/Users/David/Desktop/MoTrPAC/project_release_feb_2018/data"
SCRIPTS = "/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/"
# Assumption: these files is in the working dir (WD)
# our annotated excel file with acute samples in sheet1
# and longterm samples in sheet2. The format of these tables is the same and the basic id is the GSM id
# of a sample.
metadata_file = 'GEO_sample_metadata.xlsx' 
# The output expression database from human_ge_data_download.R
raw_data_output_obj = 'human_ge_profiles_db.RData'
# The output expression matrices from rnaseq_data_retreival.R
rnaseq_matrices_obj = "rnaseq_matrices.RData"

# Analysis constants
# specifies that we will work on transcriptomics
# in the future the analysis can be adapted to specifying
# different omics but as of March 2019 there are not enough
# datasets for a meta-analysis of other omics data.
OMIC_TYPE = "mRNA" 
MIN_NROW = 5000 # minimal number of rows in a datamatrix of a specific dataset
# Where to keep the cohort-based data and metadata
OUT_FILE_ACUTE = "human_ge_cohort_preprocessed_db_acute.RData"
OUT_FILE_LONGTERM = "human_ge_cohort_preprocessed_db_longterm.RData"
# Keep the results of the comparison of the different platforms by
# their gene coverage
GENE_FILTER_ANALYSIS = "human_ge_gene_coverage_analysis.RData"
# For the meta-analysis the basic data unit is a "gene table":
# the table of summary statistics and moderators for each gene.
# In these tables we put a row for each summary statistics. That is,
# if we have several time points then each one will have a row (as a result
# of comparing it to the base pre time point)
OUT_FILE_GENE_TABLES = "human_ge_cohort_preprocessed_db_gene_tables.RData"

# Set the working directory and load the gene expression database
setwd(WD)
library('xlsx');library(corrplot)
source(paste(SCRIPTS,'ge_download_preprocessing_helper_functions.R',sep=''))
load(raw_data_output_obj)
load(rnaseq_matrices_obj)
load('gpl_mappings_to_entrez.RData')

# Comments about the acute metadata
# Time series:
#   -1 means pre-treatment
#   Time is measured in hours
# Training - columns "Training program during experiment type"
#   The type is generally endurance or resistance
#   Untrained == controls
#   There can also be: training with treatment (e.g., LPS)
#   Column "Study subgroup" contains some additional raw information from the study description
# We exclude samples without time info - happens due to some acute/longterm mixed
# datasets such as GSE28392. We also exclude samples without subject ids, and samples 
# whose intervention includes a treatment (very few).
# Study subgroup represetns different training types or treatments within a study - 
# i.e., different interventions within a study.

# Comments about the long-term metadata
# Time series:
#   0 means pre-treatment
#   Time is measured in days
# Training - columns "Training program during experiment type"
#   The type is generally endurance or resistance
#   Untrained == controls
#   There can also be: training with treatment (e.g., LPS), or both endurance and resistance
#   Column "Study subgroup" contains some additional raw information from the study description
# We exclude samples without time info - happens due to some acute/longterm mixed
# datasets such as GSE28392. We also exclude samples without subject ids, and samples 
# whose intervention includes a treatment (very few).
# Study subgroup represetns different training types or treatments within a study - 
# i.e., different interventions within a study.

###############################################
###############################################
########## Preprocessing functions ############
###############################################
###############################################
#' 
#' Clean data from a metadata sheet: solve issues with duplications and samples without subject ids or time point information.
#' @param metadata A data frame. A table of samples (rows) vs. their annotated data (columns)
#' @description If there are duplications keep the last sample. We observed that in some datasets replicates that were done due to low quality appear LAST.
#' @return A metadata data frame in the same format as the input.
clean_raw_metadata<-function(metadata){
  # remove rows without a GSM id
  metadata = metadata[as.character(metadata[,1])!="",]
  # Solve duplications in the GSMs (happened in one dataset with several GEO ids), only relevant for acute
  gsm_duplications = names(which(table(metadata[,1])>1))
  to_keep = rep(T,nrow(metadata))
  # If there are duplications keep the last sample:
  for(curr_gsm in gsm_duplications){
    curr_ids = which(metadata[,1]==curr_gsm)
    if(length(curr_ids)<2){next}
    to_keep[curr_ids[-1]]=F
  }
  metadata = metadata[to_keep,];rownames(metadata) = metadata[,1]
  # Exclude samples without time info - happens due to some acute/longterm mixed datasets such as GSE28392
  stand_time_col = colnames(metadata)[grepl("standard",colnames(metadata),ignore.case = T)&grepl("time",colnames(metadata),ignore.case = T)]
  sample2time = as.numeric(as.character(metadata[,stand_time_col]))
  metadata = metadata[!is.na(sample2time) & sample2time!="",]
  # exclude samples without subject ids - we should contact the authors
  sample2subject = as.character(metadata[,"Subject.id"])
  names(sample2subject) = metadata[,1]
  samples_without_subject = names(which(is.na(sample2subject)|sample2subject==""))
  # get these problematic datasets and remove these samples from the  metadata table
  GSEs_without_subjects = as.character(metadata[samples_without_subject,"GSE"])
  metadata = metadata[!is.element(rownames(metadata),set=samples_without_subject),]
  return(metadata)
}

#' 
#' Simplify the training type information.
#' @param metadata A data frame. A table of samples (rows) vs. their annotated data (columns)
#' @description Parse the training type info into endurance, resistance, both, yoga, or untrained. 
#' @return A metadata data frame in the same format as the input.
simplify_training_type<-function(metadata){
  raw_training_data = as.character(metadata$Training.program.during.experiment.type)
  dataset_subgroup = as.character(metadata$Study.subgroup)
  sample_is_endurance = grepl("endur",raw_training_data,ignore.case = T)
  sample_is_resistance = grepl("resis",raw_training_data,ignore.case = T)| 
    grepl("strength",raw_training_data,ignore.case = T)
  sample_has_treatment = grepl("treatment",raw_training_data,ignore.case = T)
  # The short description of the training program
  sample2training_type = rep("other",nrow(metadata))
  sample2training_type[sample_is_endurance] = "endurance"
  sample2training_type[sample_is_resistance] = "resistance"
  sample2training_type[sample_is_endurance & sample_is_resistance] = "both"
  sample2training_type[grepl("untrained",raw_training_data,ignore.case = T)] = "untrained"
  sample2training_type[grepl("no training",raw_training_data,ignore.case = T)] = "untrained"
  sample2training_type[grepl("no exercise",raw_training_data,ignore.case = T)] = "untrained"
  sample2training_type[grepl("yoga",raw_training_data,ignore.case = T)] = "yoga"
  sample2training_type[sample_has_treatment] = paste(sample2training_type[sample_has_treatment],"treatment",sep="_")
  sample2training_type[raw_training_data==""] = ""
  names(sample2training_type) = metadata[,1]
  return(sample2training_type)
}
#' 
#' Simplify the tissue information.
#' @param metadata A data frame. A table of samples (rows) vs. their annotated data (columns)
#' @description Parse the tissue info into muscle, blood, or fat. 
#' @return A metadata data frame in the same format as the input.
simplify_tissue_info<-function(subject_col2tissue){
  subject_col2tissue_slim = as.character(subject_col2tissue)
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
#' 
#' Simplify the sex information.
#' @param metadata A data frame. A table of samples (rows) vs. their annotated data (columns)
#' @description Parse the training type info into male or female.
#' @return A metadata data frame in the same format as the input.
simplify_sex_info<-function(sexinfo){
  sexinfo = tolower(as.character(sexinfo))
  newv = as.character(sexinfo)
  names(newv) = names(sexinfo)
  newv[grepl("f",newv)] = "female"
  newv[newv!="" & newv!="female"] = "male"
  newv[newv==""] = NA
  return(newv)
}
#' 
#' A wrapper of the functions above that simplifies the information of a subject,
#' @param metadata A data frame. A table of samples (rows) vs. their annotated data (columns)
#' @return A list. Each element is a feature of the samples. Output includes sex, age, subject ids, replicates, training, and time.
get_simplified_sample_information<-function(metadata){
  # The short description of the training program
  sample2training_type = simplify_training_type(metadata)
  # # Sanity checks: should be all false
  # table(is.na(sample2training_type)|sample2training_type=="")
  # table(is.na(metadata$GSE)|metadata$GSE=="")
  # Study ids are primarily based on pubmed data
  study_ids = as.character(metadata$pmid)
  study_ids[study_ids==""] = as.character(metadata[study_ids=="","GSE"])
  names(study_ids) = metadata[,1]
  # Subject names as given in the original GEO datasets
  sample2subject = as.character(metadata[,"Subject.id"])
  names(sample2subject) = metadata[,1]
  # Dataset id -  a level below study id. It separates samples by:
  # GSE id, tissue, platform, training, study subgroup
  # The dataset ids contain the above information, separated by ';'
  dataset_ids = paste(metadata$GSE,metadata$Tissue,
                      metadata$platform_id,sample2training_type,
                      metadata$Study.subgroup,sep=';')
  names(dataset_ids) = metadata[,1]
  # Time series by sample
  stand_time_col = colnames(metadata)[grepl("standard",colnames(metadata),
                ignore.case = T)&grepl("time",colnames(metadata),ignore.case = T)]
  sample2time = as.numeric(as.character(metadata[,stand_time_col]))
  names(sample2time) = metadata[,1]
  table(sample2time)
  # Other important data
  sample2tissue = simplify_tissue_info(tolower(as.character(metadata$Tissue)))
  names(sample2tissue) = metadata[,1]
  sample2age = tolower(as.character(metadata$Age))
  sample2age = gsub(sample2age,pattern = "age: ",replace="")
  names(sample2age) = metadata[,1]
  names(sample2tissue) = metadata[,1]
  sample2sex = simplify_sex_info(tolower(as.character(metadata$Gender)))
  names(sample2sex) = metadata[,1]
  sample2replicate_info = as.character(metadata[,"Replicate.info"])
  names(sample2replicate_info) = metadata[,1]
  sample_metadata = list(age=sample2age,sex=sample2sex,replicates=sample2replicate_info,
                         time=sample2time,tissue=sample2tissue,subject=sample2subject,
                         training=sample2training_type,dataset = dataset_ids)
  return(sample_metadata)
}
#' This function analyzes the samples in analysis_samples and prepares their
#' gene expression data. We use their metadata to get the study ids and check
#' the samples time points, mapping to subject ids, and information about replicates if
#' available.
#' Unfortunately, to to the complexity of trascriptomics data we have different sources
#' of information. For microarrays data can be fRMA normalized, RMA normalized, or based
#' on data given in the GSE series matrix. For RNA-seq data we use either data preprocessed
#' by the recount database or we preprocess the read count data given in the raw data
#' of the GSEs. 
#' The preferance of data sources is as followes:
#'   microarrays: fRMA, RMA, series matrix
#'   rnaseq: recount, other information
#' For microarray data we use the GPL information objects to map probes to genes.
#' A gene's expression profile is computed by the mean of its probes.
#' For RNAseq data we assume that the data are already given in entrez genes, so
#' no averaging is required.
preprocess_expression_data<-function(metadata,analysis_samples,sample_metadata,
    CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes,rnaseq_matrices){
  dataset_ids = sample_metadata$dataset
  sample2time = sample_metadata$time
  sample2replicate_info = sample_metadata$replicates
  sample2subject = sample_metadata$subject
  dataset2preprocessed_data = list()
  for (dataset in unique(dataset_ids[analysis_samples])){
    if(is.element(dataset,set=names(dataset2preprocessed_data))){next}
    assign("last.warning", NULL, envir = baseenv())
    dataset_samples = names(which(dataset_ids==dataset))
    # cut by samples that have time points and are not repeats
    dataset_samples = dataset_samples[!is.na(sample2time[dataset_samples])]
    dataset_samples = dataset_samples[sample2replicate_info[dataset_samples] == ""]
    dataset_samples = dataset_samples[sample2subject[dataset_samples]!=""]
    if(length(dataset_samples)==0){next}
    curr_platforms = as.character(metadata[dataset_samples,"platform_id"])
    platform = curr_platforms[1]
    if(metadata[dataset_samples[1],"Type"]=="SRA"){platform = "recount"}
    
    # for convinience we order the samples by the time
    dataset_samples = dataset_samples[order(sample2time[dataset_samples])]
    
    # Order: see explanation above
    frma_mat = get_data_matrix_from_matrix_lists(CEL_frma_profiles,dataset_samples)
    rma_mat = get_data_matrix_from_matrix_lists(CEL_rma_profiles,dataset_samples)
    gse_mat = get_data_matrix_from_matrix_lists(gse_matrices,dataset_samples)
    data_matrix = frma_mat
    data_source_description = "fRMA"
    if(is.null(data_matrix)|| is.null(dim(data_matrix)) ||
       length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
      data_source_description = "RMA"
      data_matrix = rma_mat
    }
    if(is.null(data_matrix)|| is.null(dim(data_matrix)) || 
       length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
      data_source_description = "GSE series matrix or Recount"
      data_matrix = gse_mat
    }
    if(is.null(data_matrix)|| is.null(dim(data_matrix)) || 
       length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
      gse = as.character(metadata[dataset_samples,"GSE"][1])
      if(is.element(gse,set=names(rnaseq_matrices))){
        data_source_description = "RNAseq from raw data"
        dataset_samples = intersect(dataset_samples,colnames(rnaseq_matrices[[gse]]))
        data_matrix = rnaseq_matrices[[gse]][,dataset_samples]
      }
    }
    if(is.null(data_matrix) || is.null(dim(data_matrix)) ||
       length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
      data_source_description = "GE data unavailable"
    }
    print (paste("**** Analyzing the gene expression data of",dataset,"****"))
    print(paste("data_source_description:",data_source_description))
    if(data_source_description == "GE data unavailable"){next}
    
    na_rows = apply(is.na(data_matrix),1,any) | apply(is.nan(data_matrix),1,any)
    print(table(na_rows))
    data_matrix = data_matrix[!na_rows,]
    if(is.null(data_matrix)|| is.null(dim(data_matrix)) ||
       length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
      print("Data matrix has too many rows with NA or NaN values:")
      print(dataset)
      print("Skipping for now, address later")
      next
    }
    if(nrow(data_matrix) < MIN_NROW){
      print("Data matrix has less rows than the minimal number specified:")
      print(dataset)
      print("Skipping")
      next
    }
    # Analyze the data matrix
    # 1. Check if the data was logged and correct if needed
    if(max(data_matrix,na.rm=T)>50){
      data_matrix = log(data_matrix+1,base=2)
      data_matrix[is.na(data_matrix)] = min(data_matrix,na.rm=T)
    }
    # 2. get current features
    if(length(warnings())>0){print("Got some warnings, Please debug:")}
    print("Adding a new dataset to the preprocessed data containers:")
    print(paste("****",dataset,"****"))
    
    # For fold change calculations exclude subjects without a baseline or subjects with baseline only
    curr_times = sample2time[dataset_samples]
    curr_table = table(sample2subject[dataset_samples],curr_times)
    dataset_samples_for_fchange = dataset_samples
    if(is.null(dim(curr_table))||sd(curr_times,na.rm=T)==0){
      dataset_samples_for_fchange = c()
    }
    else{
      subjects_to_remove = rownames(curr_table)[curr_table[,1]==0 | rowSums(curr_table)<=1]
      if(length(subjects_to_remove)>0){
        curr_subjects = sample2subject[dataset_samples]
        samples_to_exclude = dataset_samples[is.element(curr_subjects,set=subjects_to_remove)]
        dataset_samples_for_fchange = setdiff(dataset_samples_for_fchange,samples_to_exclude)
      }
    }
    
    # Get the data matrices for further analysis
    # Transform to genes (entrez or symbols)
    if(data_source_description == "RNAseq from raw data"){
      genes_data_matrix = data_matrix
    }
    else{
      genes_data_matrix_obj = transform_matrix_into_genes(
        data_matrix,gpl = platform,gpl_mappings_entrez2probes)
      genes_data_matrix = genes_data_matrix_obj$entrez_mat
    }
    # exclude genes with NAs
    if( data_source_description != "RNAseq from raw data" &&
      sum(genes_data_matrix_obj$entrez_mat_na_stats[["row NA counts"]]==0)>10000){
      genes_data_matrix = genes_data_matrix[
        genes_data_matrix_obj$entrez_mat_na_stats[["row NA counts"]]==0,]
    }
    
    # add the results to the data containers
    dataset2preprocessed_data[[dataset]] = list()
    dataset2preprocessed_data[[dataset]][["probe_data"]] = data_matrix
    dataset2preprocessed_data[[dataset]][["gene_data"]] = genes_data_matrix
    
    # release memory
    rm(data_matrix);rm(genes_data_matrix);rm(genes_data_matrix_obj);gc()
  }
  return(dataset2preprocessed_data)
}
# We currently take the first repeat when subjects have more than
# a single time series data
get_fold_changes_vs_baseline<-function(x,subjs,timev,baseline=NULL,func = function(a,b){a-b},metadata=NULL){
  if(is.null(baseline)){baseline = sort(unique(timev))[1]}
  print(paste("baseline time point is:", baseline))
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
## Simple functions for cohort-level analyses
# some simple functions to analyze a time point in a dataset's matrix
get_paired_ttest_yi_vi <-function(x,sample2time,t1,t2){
  x1 = x[sample2time==t1]
  x2 = x[sample2time==t2]
  d = x2-x1
  n = length(d)
  sdd = sd(d)/sqrt(length(d))
  meand = mean(d)
  return(c(yi=meand,vi=sdd^2,tstat = meand/sdd,sdd=sdd,df=length(x1)-1))
}
# # test
# x1 = rnorm(10);x2=rnorm(10)+1;x=c(x1,x2)
# s2t = c(rep(0,10),rep(1,10))
# pp = get_paired_ttest_yi_vi(x,s2t,0,1)
# pp
# tt = t.test(x1,x2,paired = T)
# abs(tt$statistic) - abs(pp[1]/pp["sdd"])
# tt$estimate/tt$statistic
# tt$estimate
# tt$p.value
# 2*pt(tt$statistic,df=9)

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
  tt = t.test(x1,x2,paired = T)
  return(tt$p.value)
}
##
# Gene selection and display items
get_matrix_p_adjust<-function(x,q=0.1,...){
  v = c(x)
  v = v[!is.na(x)]
  vq = p.adjust(v,...)
  thr = max(v[vq<=q])
  return(thr)
}
###############################################
###############################################
######### Acute data preprocessing ############
###############################################
###############################################

# Sheet 1 has the acute samples metadata
metadata = clean_raw_metadata(read.xlsx2(file=metadata_file,sheetIndex=1))
sample_metadata = get_simplified_sample_information(metadata)
dim(metadata)
# Select the relevant samples
# Get the current metadata
curr_rows = rep(F,nrow(metadata))
names(curr_rows) = rownames(metadata)
# Of the selected rows use the selected omics and exercise types
if(OMIC_TYPE == "mRNA"){curr_rows = (metadata[,3] == "RNA" | metadata[,3] == "SRA")}
print(table(curr_rows))
analysis_samples = rownames(metadata)[curr_rows]

dataset2preprocessed_data = preprocess_expression_data(metadata,analysis_samples,sample_metadata,
    CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes,rnaseq_matrices)

# test_samples = analysis_samples[is.element(metadata[analysis_samples,"GSE"],set=c("GSE71972","GSE87749"))]
# dataset2preprocessed_data = preprocess_expression_data(metadata,test_samples,sample_metadata,
#     CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes,rnaseq_matrices)

# get datasets with missing matrices
all_gses = as.character(unique(metadata$GSE[is.element(metadata$Type,set=c("RNA","SRA"))]))
for(gse in all_gses){
  if(sum(grepl(gse,names(dataset2preprocessed_data)))>0){next}
  ind = which(grepl(gse,names(gse_matrices)))
  if(length(ind)>0){
    curr_obj = gse_matrices[[ind[1]]]
    if(class(curr_obj)=="matrix" && nrow(curr_obj) < MIN_NROW){next} 
  }
  print(gse)
}

# Encode the dataset ids and create their metadata info
cohort_ids = paste("GE_A_",1:length(dataset2preprocessed_data),sep='')
cohort_info = sapply(names(dataset2preprocessed_data),function(x)strsplit(x,split=';')[[1]])
cohort_metadata = list()
for(j in 1:length(dataset2preprocessed_data)){
  c_id = cohort_ids[j]
  d_id = names(dataset2preprocessed_data)[j]
  gsms = colnames(dataset2preprocessed_data[[d_id]]$gene_data)
  c_info = cohort_info[[j]]
  gse = c_info[1]
  full_tissue = c_info[2]
  gpl = c_info[3]
  training = c_info[4]
  additional_info = NA
  if(length(c_info)>4){additional_info = c_info[5]}
  tissue = simplify_tissue_info(full_tissue)
  curr_times = sample_metadata$time[gsms]
  cohort_metadata[[c_id]] = list(gsms=gsms,tissue=tissue,training=training,gse=gse,gpl=gpl,
                                 full_tissue=full_tissue,additional_info=additional_info,times=curr_times)
  if(length(unique(curr_times))<2){next}
  gene_fchanges = get_fold_changes_vs_baseline(dataset2preprocessed_data[[j]]$gene_data,sample_metadata$subject[gsms],curr_times)
  dataset2preprocessed_data[[j]][["gene_fchanges"]] = gene_fchanges
}
sapply(cohort_metadata,function(x)c(x$gse,length(unique(x$times))))
cohort_data = dataset2preprocessed_data
names(cohort_data) = cohort_ids
sample2time=sample_metadata$time;sample2sex=sample_metadata$sex;sample2age=sample_metadata$age
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_ACUTE)

sapply(cohort_data,function(x)dim(x$gene_data))
sum(sapply(cohort_data,function(x)ncol(x$gene_data)))

### Add fold changes and t-tests #######
# compute ttest p-values, yi's and vi's
for(j in 1:length(cohort_data)){
  if(length(cohort_data[[j]])<3){next}
  res1 = get_ttest_yi_vi_per_dataset(cohort_data[[j]]$gene_data,sample_metadata)
  res2 = get_ttest_pval_per_dataset(cohort_data[[j]]$gene_data,sample_metadata)
  for(nn in names(res1)){
    res1[[nn]] = cbind(res1[[nn]],res2[,nn])
    colnames(res1[[nn]])[ncol(res1[[nn]])] = "p"
  }
  cohort_data[[j]][["time2ttest_stats"]] = res1
}
sapply(cohort_data,function(x)colnames(x$time2ttest_stats[[1]]))
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_ACUTE)

###############################################
###############################################
####### longterm data preprocessing ###########
###############################################
###############################################

# Sheet 2 has the longterm samples metadata
metadata = clean_raw_metadata(read.xlsx2(file=metadata_file,sheetIndex=2))
sample_metadata = get_simplified_sample_information(metadata)

# Select the relevant samples
# Get the current metadata
curr_rows = rep(F,nrow(metadata))
names(curr_rows) = rownames(metadata)
# Of the selected rows use the selected omics and exercise types
if(OMIC_TYPE == "mRNA"){curr_rows = (metadata[,3] == "RNA" | metadata[,3] == "SRA")}
print(table(curr_rows))
analysis_samples = rownames(metadata)[curr_rows]

dataset2preprocessed_data = preprocess_expression_data(metadata,analysis_samples,sample_metadata,
      CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes,rnaseq_matrices)

# test_samples = analysis_samples[is.element(metadata[analysis_samples,"GSE"],set=c("GSE117525"))]
# test_dataset2preprocessed_data = preprocess_expression_data(metadata,test_samples,sample_metadata,
#     CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes,rnaseq_matrices)

# get datasets with missing matrices
all_gses = as.character(unique(metadata$GSE[is.element(metadata$Type,set=c("RNA","SRA"))]))
for(gse in all_gses){
  if(sum(grepl(gse,names(dataset2preprocessed_data)))>0){next}
  ind = which(grepl(gse,names(gse_matrices)))
  if(length(ind)>0){
    curr_obj = gse_matrices[[ind[1]]]
    if(class(curr_obj)=="matrix" && nrow(curr_obj) < MIN_NROW){next}
  }
  print(gse)
}

# Encode the dataset ids and create their metadata info
cohort_ids = paste("GE_L_",1:length(dataset2preprocessed_data),sep='')
cohort_info = sapply(names(dataset2preprocessed_data),function(x)strsplit(x,split=';')[[1]])
cohort_metadata = list()
for(j in 1:length(dataset2preprocessed_data)){
  c_id = cohort_ids[j]
  d_id = names(dataset2preprocessed_data)[j]
  gsms = colnames(dataset2preprocessed_data[[d_id]]$gene_data)
  c_info = cohort_info[[j]]
  gse = c_info[1]
  full_tissue = c_info[2]
  gpl = c_info[3]
  training = c_info[4]
  additional_info = NA
  if(length(c_info)>4){additional_info = c_info[5]}
  tissue = simplify_tissue_info(full_tissue)
  curr_times = sample_metadata$time[gsms]
  cohort_metadata[[c_id]] = list(gsms=gsms,tissue=tissue,training=training,gse=gse,gpl=gpl,
                                 full_tissue=full_tissue,additional_info=additional_info,times=curr_times)
  if(length(unique(curr_times))<2){next}
  gene_fchanges = get_fold_changes_vs_baseline(dataset2preprocessed_data[[j]]$gene_data,sample_metadata$subject[gsms],curr_times)
  dataset2preprocessed_data[[j]][["gene_fchanges"]] = gene_fchanges
}
cohort_data = dataset2preprocessed_data
names(cohort_data) = cohort_ids
sample2time=sample_metadata$time;sample2sex=sample_metadata$sex;sample2age=sample_metadata$age
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_LONGTERM)

# Fix: added on March 2019 - some datasets have samples with all NA/NAN values
for(j in 1:length(cohort_data)){
  curr_matrix = cohort_data[[j]]$gene_data
  curr_matrix = curr_matrix[,colSums(is.na(curr_matrix)|is.nan(curr_matrix)) < (nrow(curr_matrix)/2)]
  curr_matrix = curr_matrix[!apply(is.na(curr_matrix) | is.nan(curr_matrix),1,all),]
  cohort_metadata[[j]]$gsms = colnames(curr_matrix)
  cohort_metadata[[j]]$times = sample2time[colnames(curr_matrix)]
  cohort_data[[j]]$gene_data = curr_matrix
}

sapply(cohort_data,function(x)dim(x$gene_data))
sum(sapply(cohort_data,function(x)ncol(x$gene_data)))

### Add fold changes and t-tests #######
# compute ttest p-values, yi's and vi's
for(j in 1:length(cohort_data)){
  if(length(cohort_data[[j]])<3){next}
  res1 = get_ttest_yi_vi_per_dataset(cohort_data[[j]]$gene_data,sample_metadata)
  res2 = get_ttest_pval_per_dataset(cohort_data[[j]]$gene_data,sample_metadata)
  for(nn in names(res1)){
    res1[[nn]] = cbind(res1[[nn]],res2[,nn])
    colnames(res1[[nn]])[ncol(res1[[nn]])] = "p"
  }
  cohort_data[[j]][["time2ttest_stats"]] = res1
}
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_LONGTERM)

# At this point we have two RData files with all cohorts and all data
# in our annotated resource.
# For the meta-analysis we need to take care of three tasks:
# 1. Exclude platforms with extremely low gene coverage
# 2. Impute sex
# 3. Create gene tables

###############################################
###############################################
# Analysis of gene coverage: get the worst datasets #
###############################################
###############################################

load(OUT_FILE_LONGTERM)
ge_matrices1 = lapply(cohort_data,function(x)x$gene_data)
load(OUT_FILE_ACUTE)
ge_matrices2 = lapply(cohort_data,function(x)x$gene_data)
for(nn in names(ge_matrices2)){
  ge_matrices1[[nn]] = ge_matrices2[[nn]]
}

get_list_intersect<-function(l){
  x = l[[1]]
  for(i in 1:length(l)){
    x = intersect(x,l[[i]])
  }
  return(x)
}
gene_lists = lapply(ge_matrices1,rownames)
all_genes = get_list_intersect(gene_lists)
# remove the worst dataset at a time
included_datasets = names(ge_matrices1)
low_coverage_platforms = c()
while(length(included_datasets)>1){
  curr_sizes = c()
  for(j in 1:length(included_datasets)){
    curr_sizes[j] = length(get_list_intersect(gene_lists[included_datasets[-j]]))
  }
  ind = which(curr_sizes==max(curr_sizes))[1]
  excluded = included_datasets[ind]
  included_datasets = included_datasets[-ind]
  all_genes = get_list_intersect(gene_lists[included_datasets])
  print(paste(length(all_genes),excluded))
  low_coverage_platforms = c(low_coverage_platforms,excluded)
  if(length(all_genes)>10000){break}
}

# Save two objects: 
# The set of genes covered in most cases
# The set of cohorts from low coverage plarforms
save(all_genes,low_coverage_platforms,file=GENE_FILTER_ANALYSIS)

###############################################
###############################################
######## Sex imputation using ML flow #########
###############################################
###############################################
stats_matrix = c()
load(OUT_FILE_LONGTERM)
ge_matrices1 = lapply(cohort_data,function(x)x$gene_data)
gses1 = sapply(cohort_metadata, function(x)x$gse)
tissue1 = sapply(cohort_metadata, function(x)x$tissue)
sex1 = sample2sex
for(nn in names(cohort_metadata)){
  stats_matrix = rbind(stats_matrix,c(nn,ncol(cohort_data[[nn]]$gene_data),cohort_metadata[[nn]]$tissue))
}
longterm_gsms = unique(unlist(sapply(ge_matrices1,colnames)))
load(OUT_FILE_ACUTE)
ge_matrices2 = lapply(cohort_data,function(x)x$gene_data)
gses2 = sapply(cohort_metadata, function(x)x$gse)
tissue2 = sapply(cohort_metadata, function(x)x$tissue)
sex2 = sample2sex
acute_gsms = unique(unlist(sapply(ge_matrices2,colnames)))
for(nn in names(cohort_metadata)){
  stats_matrix = rbind(stats_matrix,c(nn,ncol(cohort_data[[nn]]$gene_data),cohort_metadata[[nn]]$tissue))
}
stats_matrix = as.data.frame(stats_matrix)
stats_matrix = cbind(stats_matrix,grepl("GE_L",stats_matrix[,1]))
stats_matrix$V2 = as.numeric(as.character(stats_matrix $V2))
aggregate(stats_matrix$V2, by=list(Category=stats_matrix$V3,L=stats_matrix$`grepl("GE_L", stats_matrix[, 1])`), FUN=sum)
length(unique(c(gses1,gses2)))

# merge the gene expression data before the analysis below
for(nn in names(ge_matrices2)){
  ge_matrices1[[nn]] = ge_matrices2[[nn]]
}
gses1[names(gses2)]=gses2
tissue1[names(tissue2)]=tissue2
sex1[names(sex2)] = sex2

# Define our y target: sex
y=sex1
table(y)

# Merge the gene expression profiles
x = c(); z = c(); y2=c()
for(nn in setdiff(names(ge_matrices1),low_coverage_platforms)){
  currx = ge_matrices1[[nn]][all_genes,]
  x = cbind(x,currx)
  currz = rep(gses1[nn],ncol(currx))
  z[colnames(currx)] = currz
  y2[colnames(currx)] = rep(tissue1[nn],ncol(currx))
}

samp_names = intersect(names(y),colnames(x))
samp_names = samp_names[y2[samp_names]!="fat"]
y = y[samp_names];y2 = y2[samp_names]
x = x[,samp_names];z = z[samp_names]
# some useful stats of the data
length(unique(union(acute_gsms,longterm_gsms)))
table(y[acute_gsms])
table(y2[acute_gsms])
table(y2[longterm_gsms])

# Define the set of samples with missing sex information
missing_set = is.na(y)
# table(missing_set)

library(GenomicRanges)
library(Homo.sapiens)
library(DESeq2)
geneRanges <- function(db=Homo.sapiens, column="ENTREZID"){
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}
entrez_gr = geneRanges(Homo.sapiens, column="ENTREZID")
entrez_gr = as.data.frame(entrez_gr)
selected_genes = entrez_gr[entrez_gr$seqnames=="chrY" | entrez_gr$seqnames=="chrX","ENTREZID"]

# # normalize x using rank transform
# quntx = normalizeQuantiles(x)
# quntx_sex = quntx[intersect(selected_genes,rownames(x)),]
# boxplot(quntx[,sample(1:ncol(x))[1:20]])
# dim(quntx_sex)

# Try another normalization
getRankedBasedProfile<-function(x,fs=NULL){
  if (!is.null(fs)){
    x = x[fs]
  }
  # get the ranks in decreasing order
  N = length(x)
  rs = N-rank(x,ties.method = "average")+1
  w_rs = rs*exp(-rs/N)
  names(w_rs) = names(rs)
  return (w_rs)
}
newx = apply(x,2,getRankedBasedProfile)
newx_sex = newx[intersect(selected_genes,rownames(x)),]
boxplot(newx[,sample(1:ncol(x))[1:20]])
dim(newx_sex)

inds = !missing_set 
table(inds)
source("/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions_classification.R")
# lso_res_qx_sex = leave_study_out2(as.factor(y[inds]),
#                 t(quntx_sex[,inds]),z[inds],
#                 func = svm,class.weights=c("female"=10,"male"=1),
#                 pred_args=list(probability=T),probability=T,kernel="linear")
lso_res_newx_sex = leave_study_out2(as.factor(y[inds]),
                t(newx_sex[,inds]),z[inds],
                func = svm,class.weights=c("female"=10,"male"=1),
                pred_args=list(probability=T),probability=T,kernel="linear")

preds = c()
lso_res = lso_res_newx_sex
for(nn in names(lso_res)){
  preds = rbind(preds,cbind(
    as.numeric(lso_res[[nn]]$real == "female"),
    attr(lso_res[[nn]]$preds,"probabilities")[,"female"]))
}
library(pROC)
boxplot(preds[,2]~preds[,1])
calcAupr(preds[,2],preds[,1])
calcAupr(preds[,2],preds[,1],roc = T,useAbs = T)
as.numeric(pROC::auc(preds[,1], preds[,2]))

# Impute sex and add it to the files
tr = newx_sex[,!missing_set]
svm_model = svm(x=t(tr),y=as.factor(y[!missing_set]),kernel="linear",class.weights=c("female"=10,"male"=1))
te = newx_sex[,missing_set]
preds = predict(svm_model,t(te))
table(preds)

load(OUT_FILE_ACUTE)
inds = intersect(names(preds),names(sample2sex))
sample2sex[inds] = as.character(preds[inds])
table(preds[inds])
table(sample2sex)
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_ACUTE)

load(OUT_FILE_LONGTERM)
inds = intersect(names(preds),names(sample2sex))
sample2sex[inds] = as.character(preds[inds])
table(preds[inds])
table(sample2sex)
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_LONGTERM)

###############################################
###############################################
### Reshape the data for the meta-analysis ####
###############################################
###############################################

#' Put the moderator information (covariates) of a dataset in a single object
#' with: sex,age,tissue, training type, gse
#' Assumption: "other" without a description for a training type == "untrained"
#' @param metadata A named list. Each item represents the metadata of a cohort.
#' @return A named vector with the study id, tissue, training type, avg age, and proportion of males.
get_dataset_moderators<-function(metadata){
  arrs = t(sapply(metadata,function(x)c(
    x$gse,x$tissue,x$training,x$avg_age,x$age_sd,x$male_prop,length(x$gsms))))
  colnames(arrs) = c("gse","tissue","training","avg_age","age_sd","prop_males","N")
  return(arrs)
}

#' Put all data of a single gene in a single object
#' @param gene A character or an index. The gene to analyze.
#' @param dataset_effects. A list with the summary statistics.
#' @param moderators. The study covariates to add to each data frame.
#' @return A data frame with the information for the gene.
#' @description Some datasets, RNAseq likely have zero variance for some genes, this creates NAs and NANs in the summary statistics. Whenever the sdd is zero we put zero in yi and tstat but keep the p-value as NaN
get_gene_table<-function(gene,dataset_effects,moderators){
  m = c()
  for(nn in names(dataset_effects)){
    if(! gene %in% rownames(dataset_effects[[nn]][[1]])){next}
    mm = sapply(dataset_effects[[nn]],function(u,v)u[v,],v=gene)
    for(j in 1:ncol(mm)){
      m = rbind(m,c(nn,colnames(mm)[j],moderators[nn,],mm[,j]))
    }
  }
  if(is.null(m)){return(NULL)}
  if(is.null(dim(m))){m = as.matrix(m,nrow=1)}
  colnames(m)[2]="time"
  m = data.frame(m,stringsAsFactors=F)
  # Change variables from character to numeric 
  m = transform(m,yi=as.numeric(yi),vi=as.numeric(vi),time=as.numeric(time),
                avg_age = as.numeric(avg_age),N=as.numeric(N),p=as.numeric(p),
                df=as.numeric(df),prop_males=as.numeric(prop_males),age_sd=as.numeric(age_sd),
                tstat = as.numeric(tstat))
  # Some datasets, RNAseq likely have zero variance for some genes,
  # this creates NAs and NANs in the summary statistics. 
  # Whenever the sdd is zero we put zero in yi and tstat but keep the p-value as NaN
  nan_inds = m$sdd==0
  m[nan_inds,c("yi","tstat")] = 0
  m[nan_inds,"sdd"] = 1e-6 # put some low number instead of zero
  return(m)
}

load(OUT_FILE_LONGTERM)
metadata = clean_raw_metadata(read.xlsx2(file=metadata_file,sheetIndex=2))
sample_metadata = get_simplified_sample_information(metadata)
for(nn in names(cohort_metadata)){
  samps = cohort_metadata[[nn]]$gsms
  curr_ages = as.numeric(metadata[samps,]$Numeric_Age)
  curr_raw_ages = as.character(metadata[samps,]$Age)
  is_male = sample2sex[samps] == "male"
  # print(table(is_male))
  is_male = is_male[!is.na(is_male)]
  curr_p = sum(is_male)/length(is_male)
  cohort_metadata[[nn]]$avg_age = mean(curr_ages,na.rm=T)
  if(is.nan(cohort_metadata[[nn]]$avg_age)){break}
  cohort_metadata[[nn]]$age_sd = sd(curr_ages,na.rm=T)
  if(any(grepl("(\\±|\\+)+",curr_raw_ages))){
    currsd = strsplit(as.character(curr_raw_ages[1]),split="(\\±|\\+)+")[[1]]
    currsd = gsub("\\D+$","",currsd[[length(currsd)]])
    currsd = gsub("^\\D+","",currsd[[length(currsd)]])
    cohort_metadata[[nn]]$age_sd = currsd
  }
  cohort_metadata[[nn]]$male_prop = curr_p
  print(paste(cohort_metadata[[nn]]$avg_age,cohort_metadata[[nn]]$age_sd,cohort_metadata[[nn]]$male_prop))
}

all_covered_genes = unique(unlist(sapply(cohort_data,function(x)rownames(x$gene_data))))
moderators = get_dataset_moderators(cohort_metadata)
data_datasets_effects = lapply(cohort_data,function(x)x$time2ttest_stats)
data_datasets_effects = data_datasets_effects[sapply(data_datasets_effects,length)>0]
gene_tables = lapply(all_covered_genes,get_gene_table,
                              dataset_effects=data_datasets_effects,moderators=moderators)
names(gene_tables) = all_covered_genes
longterm_gene_tables = gene_tables
rm(gene_tables);gc()

load(OUT_FILE_ACUTE)
metadata = clean_raw_metadata(read.xlsx2(file=metadata_file,sheetIndex=1))
sample_metadata = get_simplified_sample_information(metadata)
for(nn in names(cohort_metadata)){
  samps = cohort_metadata[[nn]]$gsms
  curr_ages = as.numeric(metadata[samps,]$Numeric_Age)
  curr_raw_ages = as.character(metadata[samps,]$Age)
  is_male = sample2sex[samps] == "male"
  # print(table(is_male))
  is_male = is_male[!is.na(is_male)]
  curr_p = sum(is_male)/length(is_male)
  cohort_metadata[[nn]]$avg_age = mean(curr_ages,na.rm=T)
  if(is.nan(cohort_metadata[[nn]]$avg_age)){break}
  cohort_metadata[[nn]]$age_sd = sd(curr_ages,na.rm=T)
  if(any(grepl("(\\±|\\+)+",curr_raw_ages))){
    currsd = strsplit(as.character(curr_raw_ages[1]),split="(\\±|\\+)+")[[1]]
    currsd = gsub("\\D+$","",currsd[[length(currsd)]])
    currsd = gsub("^\\D+","",currsd[[length(currsd)]])
    cohort_metadata[[nn]]$age_sd = currsd
    print(paste(curr_raw_ages[1],currsd))
  }
  cohort_metadata[[nn]]$male_prop = curr_p
  # print(paste(cohort_metadata[[nn]]$avg_age,cohort_metadata[[nn]]$age_sd,cohort_metadata[[nn]]$male_prop))
}

curr_datasets = !sapply(cohort_data,function(x)is.null(x[["time2ttest_stats"]]))
all_covered_genes = unique(unlist(sapply(cohort_data[curr_datasets],
                                         function(x)rownames(x$gene_data))))
gene_sets = lapply(cohort_data,function(x)rownames(x$gene_data))
moderators = get_dataset_moderators(cohort_metadata)
data_datasets_effects = lapply(cohort_data,function(x)x$time2ttest_stats)
data_datasets_effects = data_datasets_effects[sapply(data_datasets_effects,length)>0]
gene_tables = lapply(all_covered_genes,get_gene_table,
                     dataset_effects=data_datasets_effects,moderators=moderators)
names(gene_tables) = all_covered_genes
acute_gene_tables = gene_tables

# # for QA do a for loop
# gene_tables = list()
# for(gene in all_covered_genes){
#   gene_tables[[gene]] = get_gene_table(gene,dataset_effects=data_datasets_effects,moderators=moderators)
#   print(gene)
#   for(nn in names(gene_sets)){
#     if(is.element(gene,set=gene_sets[[nn]])){
#       print(nn)
#     }
#   }
# }

save(longterm_gene_tables,acute_gene_tables,file=OUT_FILE_GENE_TABLES)



