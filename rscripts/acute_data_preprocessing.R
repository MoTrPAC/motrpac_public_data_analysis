# This script loads the raw profiles and the metadata
# creates standardized datasets and metadata and saves
# it in easy to use objects
setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library('xlsx');library('GEOquery');library(corrplot)
source('repos/motrpac/helper_functions.R')

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

###############################################
###############################################
############# Load the data ###################
###############################################
###############################################

metadata_file = 'GEO_sample_metadata.xlsx'
# Sheet 1 has the acute samples metadata
metadata = read.xlsx2(file=metadata_file,sheetIndex=1)
metadata = metadata[as.character(metadata[,1])!="",]
gsm_duplications = names(which(table(metadata[,1])>1))
# Solve duplications in the GSMs (happened in one dataset with several GEO ids)
to_keep = rep(T,nrow(metadata))
for(curr_gsm in gsm_duplications){
  curr_ids = which(metadata[,1]==curr_gsm)
  to_keep[curr_ids[-1]]=F
}
metadata = metadata[to_keep,]
rownames(metadata) = metadata[,1]

# Exclude samples without time info - happens due to some acute/longterm mixed datasets
# datasets such as GSE28392
sample2time = sample2time = as.numeric(as.character(metadata$Acute..Standardized.Time..hours...1.is.baseline.))
metadata = metadata[!is.na(sample2time) & sample2time!="",]
print(dim(metadata))
# exclude samples without subject ids - we should contact the authors
sample2subject = as.character(metadata[,"Subject.id"])
names(sample2subject) = metadata[,1]
samples_without_subject = names(which(is.na(sample2subject)|sample2subject==""))
# get these problematic datasets and remove these samples from the  metadata table
GSEs_without_subjects = as.character(metadata[samples_without_subject,"GSE"])
table(GSEs_without_subjects)
metadata = metadata[!is.element(rownames(metadata),set=samples_without_subject),]
print(dim(metadata))

# Get sample information
# The short description of the training program
sample2training_type = simplify_training_type(metadata)
table(sample2training_type)
# Sanity checks
table(is.na(sample2training_type)|sample2training_type=="")
table(is.na(metadata$GSE)|metadata$GSE=="")
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
dataset_ids = paste(metadata$GSE,metadata$Tissue,metadata$platform_id,sample2training_type,metadata$Study.subgroup,sep=';')
names(dataset_ids) = metadata[,1]
sort(table(dataset_ids))[1:5]
# Time series by sample
sample2time = as.numeric(as.character(metadata$Acute..Standardized.Time..hours...1.is.baseline.))
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

# Check for irregularities in subject ids
# Printed tables - do not have the exact same values
# Datasets with some zeroes are okay - these are datasets with
# some missing time points per subject. Not a problem.
for(dataset in unique(dataset_ids)){
  curr_samples = names(which(dataset_ids==dataset))
  curr_samples_rep_info = sample2replicate_info[curr_samples]
  # Keep samples that are not repeats and have times
  curr_samples = curr_samples[curr_samples_rep_info == ""]
  curr_times = sample2time[curr_samples]
  curr_samples = curr_samples[curr_times !=""]
  curr_samples = curr_samples[sample2subject[curr_samples]!=""]
  curr_subjects = sample2subject[curr_samples]
  curr_times = sample2time[curr_samples]
  curr_table = table(curr_subjects,curr_times)
  if(length(curr_table)==0){next}
  # Exclude for this test, subjects without the baseline time point
  curr_table = curr_table[curr_table[,1]>0,]
  if(all(c(curr_table)==c(curr_table)[1])){next}
  print(dataset)
  print(curr_table)
  #break
}

sample_metadata = list(age=sample2age,sex=sample2sex,
                       time=sample2time,tissue=sample2tissue,subject=sample2subject,training=sample2training_type,dataset = dataset_ids)
save(sample_metadata,file="PADB_sample_metadata_acute.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############### Set constants #################
###############################################
###############################################

# Analysis constants
# set constants
# Set constants for the analysis
OMIC_TYPE = "mRNA" # mRNA, methylation
EXERCISE_TYPE = "both" # endurance, resistance, both, or other (yoga)
MIN_NROW = 5000 # minimal number of rows in a datamatrix of a specific dataset
OUT_FILE = "PADB_univariate_results_and_preprocessed_data_acute.RData"

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############### Preprocessing #################
###############################################
###############################################

# load the gene expression database
load('gpl_mappings_to_entrez.RData')
load("PA_database_profiles.RData")
load(OUT_FILE)

# Select the relevant samples
# Get the current metadata
curr_rows = rep(F,nrow(metadata))
names(curr_rows) = rownames(metadata)
# Of the selected rows use the selected omics and exercise types
if(OMIC_TYPE == "mRNA"){curr_rows = (metadata[,3] == "RNA" | metadata[,3] == "SRA")}
print(table(curr_rows))
analysis_samples = rownames(metadata)[curr_rows]

# We start with a simple fold-change exploratory analysis
# The loop below gets the gene and probe matrices for our datasets
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
  
  # Order: fRMA, RMA, GSE object, GSM profiles
  frma_mat = get_data_matrix_from_matrix_lists(CEL_frma_profiles,dataset_samples)
  rma_mat = get_data_matrix_from_matrix_lists(CEL_rma_profiles,dataset_samples)
  gse_mat = get_data_matrix_from_matrix_lists(gse_matrices,dataset_samples)
  print(dim(gse_mat))
  data_matrix = frma_mat
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
    print("Not enough samples in fRMA object, use RMA:")
    print(dataset)
    data_matrix = rma_mat
  }
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
    print("Not enough samples in RMA object, use GSEs:")
    print(dataset)
    data_matrix = gse_mat
  }
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
    print("Not enough samples in GSE object:")
    print(dataset)
    print("Skipping for now, address later")
    next
  }
  na_rows = apply(is.na(data_matrix),1,any) | apply(is.nan(data_matrix),1,any)
  print(table(na_rows))
  data_matrix = data_matrix[!na_rows,]
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
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
    data_matrix = log(data_matrix,base=2)
    data_matrix[is.na(data_matrix)] = min(data_matrix,na.rm=T)
  }
  
  # 2. get current features
  if(length(warnings())>0){print("Got some warnings, Please debug:")}
  
  print("Adding a new dataset to the preprocessed data containers:")
  print(paste("****************",dataset))
  
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
  genes_data_matrix_obj = transform_matrix_into_genes(data_matrix,gpl = platform,gpl_mappings_entrez2probes)
  genes_data_matrix = genes_data_matrix_obj$entrez_mat
  # exclude genes with NAs
  if(sum(genes_data_matrix_obj$entrez_mat_na_stats[["row NA counts"]]==0)>10000){
    genes_data_matrix = genes_data_matrix[genes_data_matrix_obj$entrez_mat_na_stats[["row NA counts"]]==0,]
  }
  # get the fold change matrices
  gene_fold_changes = NULL; probe_fold_changes = NULL
  if(length(dataset_samples_for_fchange)>0){
    curr_subjects = sample2subject[dataset_samples_for_fchange]
    # probe_fold_changes = get_fold_changes_vs_baseline(data_matrix[,dataset_samples_for_fchange],curr_subjects,curr_times[dataset_samples_for_fchange])
    gene_fold_changes = get_fold_changes_vs_baseline(genes_data_matrix[,dataset_samples_for_fchange],curr_subjects,curr_times[dataset_samples_for_fchange])
    corrplot(cor(gene_fold_changes),order="hclust")
  }
  
  # add the results to the data containers
  dataset2preprocessed_data[[dataset]] = list()
  dataset2preprocessed_data[[dataset]][["probe_data"]] = data_matrix
  dataset2preprocessed_data[[dataset]][["gene_data"]] = genes_data_matrix
  # dataset2preprocessed_data[[dataset]][["probe_fold_changes"]] = probe_fold_changes
  dataset2preprocessed_data[[dataset]][["gene_fold_changes"]] = gene_fold_changes
  
  # release memory and save
  rm(data_matrix);rm(genes_data_matrix);rm(genes_data_matrix_obj);gc()
  save(dataset2preprocessed_data,file=OUT_FILE)
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
  cohort_metadata[[c_id]] = list(gsms=gsms,tissue=tissue,training=training,gse=gse,gpl=gpl,
                                 full_tissue=full_tissue,additional_info=additional_info)
  curr_times = sample2time[gsms]
  if(length(unique(curr_times))<2){next}
  gene_fchanges = get_fold_changes_vs_baseline(dataset2preprocessed_data[[j]]$gene_data,sample2subject[gsms],curr_times)
  dataset2preprocessed_data[[j]][["gene_fchanges"]] = gene_fchanges
}
names(dataset2preprocessed_data) = cohort_ids
cohort_data = dataset2preprocessed_data
save(cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE)

rm(CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes,gpl_mappings_to_entrez,gpl_tables)
gc()

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
########### Fold changes and t-tests ##########
###############################################
###############################################
load(OUT_FILE)
load("PADB_sample_metadata_acute.RData")
acute_metadata = get(load("PADB_sample_metadata_acute.RData"))

# # test
# res1 = get_ttest_yi_vi_per_dataset(cohort_data[[1]]$gene_data,acute_metadata)
# mat1 = cohort_data[[1]]$gene_fchanges
# res2 = rowMeans(mat1[,grepl(colnames(mat1),pattern='2.5')])
# max(abs(res1$`2.5`[,"yi"]-res2))

# compute ttest p-values, yi's and vi's
for(j in 1:length(cohort_data)){
  if(length(cohort_data[[j]])<3){next}
  res1 = get_ttest_yi_vi_per_dataset(cohort_data[[j]]$gene_data,acute_metadata)
  res2 = get_ttest_pval_per_dataset(cohort_data[[j]]$gene_data,acute_metadata)
  for(nn in names(res1)){
    res1[[nn]] = cbind(res1[[nn]],res2[,nn])
    colnames(res1[[nn]]) = c("yi","vi","p")
  }
  cohort_data[[j]][["time2ttest_stats"]] = res1
}
table(sapply(cohort_data,function(x)is.null(x$gene_fold_changes)))
save(cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE)

###############################################
###############################################
#################### End ######################
###############################################
###############################################



