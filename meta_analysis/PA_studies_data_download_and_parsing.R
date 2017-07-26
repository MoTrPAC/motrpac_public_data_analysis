setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library('xlsx');library('GEOquery')
source('helper_functions.R')

metadata_file = 'GEO_sample_metadata.xlsx'
# Before June 2017 all samples where in one sheet.
# To make the database more convinient the acute and long-term annotations were separated.
# Nevertheless the old sheet has all ids and can be used to download the data.
metadata = read.xlsx2(file=metadata_file,sheetIndex=3)
metadata = metadata[as.character(metadata[,1])!="",]
gsm_duplications = names(which(table(metadata[,1])>1))
# Solve duplications in the GSMs (happened in one dataset with several GEO ids)
to_keep = rep(T,nrow(metadata))
for(curr_gsm in gsm_duplications){
  curr_ids = which(metadata[,1]==curr_gsm)
  to_keep[curr_ids[-1]]=F
}
table(to_keep)
metadata = metadata[to_keep,]
rownames(metadata) = metadata[,1]

GEO_destdir = "GEO"

# Some stats
has_acute = as.character(metadata$Acute..Standardized.Time..hours...1.is.baseline.)!=""
has_longterm = as.character(metadata$Long.term..Standardized.Time..days.)!=""
table(has_acute,has_longterm)
table(has_acute)
table(has_longterm)
sample_is_endurance = grepl("endur",metadata$Training.program.during.experiment.type,ignore.case = T)
sample_is_resistance = grepl("resis",metadata$Training.program.during.experiment.type,ignore.case = T)| 
  grepl("strength",metadata$Training.program.during.experiment.type,ignore.case = T)
table(sample_is_endurance,sample_is_resistance)
# Moreover, this preprocessing determines the analysis
# General comments:
#   1. Any cross validation should be done at the study level
#   2. Meta-analysis can be done at the dataset level
#   3. Dataset ids for "acute" data analysis: use the training program in the dataset name

# Study ids are primarily based on pubmed data
study_ids = as.character(metadata[,5])
study_ids[study_ids==""] = as.character(metadata[study_ids=="",2])
names(study_ids) = metadata[,1]
# Subject names as given in the original GEO datasets
sample2subject = as.character(metadata[,"Subject.id"])
names(sample2subject) = metadata[,1]
# The short description of the training program
sample2training_program = as.character(metadata$Training.program.during.experiment.type)
names(sample2training_program) = metadata[,1]
# Dataset id -  a level below study id. It separates samples by:
# GSE id, tissue, platform
# The dataset ids contain the above information, separated by ';'
dataset_ids = paste(metadata$GSE,metadata$Tissue,metadata$platform_id,sep=';')
names(dataset_ids) = metadata[,1]
sort(table(dataset_ids))[1:5]
# Dataset ids for "acute" data analysis: use the training program in the dataset name
# This is done for acute only because of the acute-long mixed datasets
dataset_ids_acute = paste(dataset_ids,sample2training_program,sep=";")
names(dataset_ids_acute) = metadata[,1]
dataset_ids_longterm = dataset_ids
names(dataset_ids_longterm) = metadata[,1]
time_col_acute = which(grepl("acute", colnames(metadata),ignore.case = T) & 
                         grepl("standardized", colnames(metadata),ignore.case = T))
time_col_longterm = which(grepl("long", colnames(metadata),ignore.case = T) & 
                            grepl("standardized", colnames(metadata),ignore.case = T))
# Other important data
sample2tissue = tolower(as.character(metadata$Tissue))
names(sample2tissue) = metadata[,1]
sample2within_study_group = as.character(metadata$Exercise.group.)
names(sample2within_study_group) = metadata[,1]
sample2age = tolower(as.character(metadata$Age))
sample2age = gsub(sample2age,pattern = "age: ",replace="")
names(sample2age) = metadata[,1]
sample2sex = tolower(as.character(metadata$Gender))
names(sample2sex) = metadata[,1]
sample2replicate_info = as.character(metadata[,"Replicate.info"])
names(sample2replicate_info) = metadata[,1]

# (this marks the start/end of a subsection of the analysis)
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# Check for irregularities in subject ids
# THIS CODE MUST RUN WITHOUT ANY OUTPUT BEFORE MOVING ON
for(dataset in unique(dataset_ids_acute)){
  curr_samples = names(which(dataset_ids_acute==dataset))
  curr_samples_rep_info = as.character(metadata[curr_samples,"Replicate.info"])
  # Keep samples that are not repeats and have times
  curr_samples = curr_samples[curr_samples_rep_info == ""]
  curr_times = as.character(metadata[curr_samples,time_col_acute])
  curr_samples = curr_samples[curr_times !=""]
  curr_samples = curr_samples[sample2subject[curr_samples]!=""]
  curr_subjects = sample2subject[curr_samples]
  curr_times = as.character(metadata[curr_samples,time_col_acute])
  curr_table = table(curr_subjects,curr_times)
  if(length(curr_table)==0){next}
  # Exclude for this test, subjects without the baseline time point
  curr_table = curr_table[curr_table[,1]>0,]
  if(all(c(curr_table)==c(curr_table)[1])){next}
  print(dataset)
  print(curr_table)
  break
}
for(dataset in unique(dataset_ids_longterm)){
  curr_samples = names(which(dataset_ids_longterm==dataset))
  curr_samples_rep_info = as.character(metadata[curr_samples,"Replicate.info"])
  # Keep samples that are not repeats and have times
  curr_samples = curr_samples[curr_samples_rep_info == ""]
  curr_times = as.character(metadata[curr_samples,time_col_longterm])
  curr_samples = curr_samples[curr_times !=""]
  curr_samples = curr_samples[sample2subject[curr_samples]!=""]
  curr_subjects = sample2subject[curr_samples]
  curr_times = as.character(metadata[curr_samples,time_col_longterm])
  curr_table = table(curr_subjects,curr_times)
  if(length(curr_table)==0){next}
  # Exclude for this test, subjects without the baseline time point
  curr_table = curr_table[curr_table[,1]>0,]
  if(all(curr_table==c(curr_table)[1])){next}
  print(dataset)
  print(curr_table)
  break
}

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# Download and save the expression data

####### Get the data by taking the series matrices
geo_ids = unique(as.character(metadata[,2]))
geo_ids = geo_ids[grepl("^GSE",geo_ids,perl=T) | grepl("^GDS",geo_ids,perl=T)]
# download the pheno and meta data
geo_objs = list()
gse_fails = c()
for (gse in geo_ids){
  try({
    geo_objs[[gse]] = getGEO(gse,destdir = GEO_destdir,GSEMatrix = T,getGPL = F)
  })
}
# failed downloads:
failed_gses = setdiff(geo_ids,names(geo_objs))
failed_gses_data = metadata[is.element(metadata[,2],set=failed_gses),]
table(failed_gses_data[,3])
# see the object classes
sapply(geo_objs,function(x)sapply(x,class)) # all should be ExpressionSet

####### Go directly to the GSMs
gsms = unique(metadata[,1])
sra_gsms = metadata[grepl("SRA",metadata[,3]),1]
# download the database
gsm_objs = list()
for(gsm in gsms){
	gsm_objs[[gsm]] = ""
	try({gsm_objs[[gsm]] = getGEO(gsm,destdir = GEO_destdir)})
	ind = length(gsm_objs)
	if (ind %% 10 ==0 ){print(ind)}
}
failed_gsms = names(which(sapply(gsm_objs,class)=="character"))
failed_gsms_data = metadata[is.element(metadata[,1],set=failed_gsms),]
failed_gsms_data[,3]
table(sapply(gsm_objs,class))

####### Get and analyze the supplementary CEL files
supplementary_CEL_files = c()
for(gsm in names(gsm_objs)){
  if(gsm=="" || class(gsm_objs[[gsm]])=="character"){next}
  supp_files = Meta(gsm_objs[[gsm]])$supplementary_file
  supp_files = supp_files[grepl("CEL",supp_files)]
  if(length(supp_files)<1){next}
  supplementary_CEL_files[gsm] = supp_files[1]
}
# download the CEL files
for(gsm in names(supplementary_CEL_files)){
  dstfile = paste(GEO_destdir,"/",gsm,".CEL",sep="")
  if(is.element(dstfile,set=paste(GEO_destdir,list.files(GEO_destdir),sep='/'))){next}
  print(gsm)
  # try downloading each file several times
  download_success = F
  for (i in 1:10){
    if(download_success){break}
    try({
      download.file(supplementary_CEL_files[gsm],destfile = dstfile)
      download_success = T
    })
  }
}
# Run preprocessing of the CELs using frma
#source("https://bioconductor.org/biocLite.R")
#biocLite("frma")
#biocLite('hugene10stv1frmavecs')
#biocLite("oligo")
#biocLite("hugene.1.0.st.v1frmavecs")
#biocLite("hgu133plus2frmavecs")
#biocLite("annotate",dependencies=T)
#biocLite("hgu219frmavecs");biocLite("pd.hg.u219");biocLite("hgu219.db")
# COMMENT: frma failed on hgu219frmavecs - no such package 
library('frma'); library('affy');library('oligo');library('hugene.1.0.st.v1frmavecs')
library("hgu133plus2frmavecs"); library("hgu219frmavecs")
curr_gsms = names(supplementary_CEL_files)
curr_files = paste(GEO_destdir,"/",curr_gsms,".CEL",sep="")
names(curr_files) = curr_gsms
curr_gse_and_platform = paste(metadata[curr_gsms,2],metadata[curr_gsms,"platform_id"],sep="_")
names(curr_gse_and_platform) = curr_gsms
CEL_frma_profiles = list(); CEL_rma_profiles = list()
batches = unique(curr_gse_and_platform)
for(CEL_batch in batches){
  print(CEL_batch)
  if(is.element(CEL_batch,set=names(CEL_frma_profiles))){next}
  # try using the oligo or the affy packages
  gsm_lst = curr_files[curr_gse_and_platform==CEL_batch]
  #metadata[names(gsm_lst),2][1]
  frma_success = F
  try({
    affybatch = read.celfiles(filenames = gsm_lst)
    rma_obj = rma(affybatch)
    frma_obj = frma(affybatch,target="core")
    frma_success = T
  })
  if  (!frma_success){
    try({
      affybatch = ReadAffy(filenames = gsm_lst)
      frma_obj = frma(affybatch,target="core")
      frma_success = T
    })
  }
  if(!frma_success){
    print ("fRMA ERROR")
    rma_data = exprs(rma_obj)
    CEL_rma_profiles[[CEL_batch]] = rma_data
    next
  }
  frma_data = exprs(frma_obj)
  rma_data = exprs(rma_obj)
  print(cor(frma_data[,1:3],rma_data[,1:3]))
  # Compare these new calculations to the gsm object
  gsm = colnames(frma_data)[1]
  gsm = gsub(gsm,pattern = ".CEL$",replace='')
  gsm_original_data = Table(gsm_objs[[gsm]])[,2]
  names(gsm_original_data) = Table(gsm_objs[[gsm]])[,1]
  frma_v = frma_data[,1]
  shared_probes = intersect(names(gsm_original_data),names(frma_v))
  #plot(frma_v[shared_probes],gsm_original_data[shared_probes])
  print(cor(frma_v[shared_probes],gsm_original_data[shared_probes],method="spearman"))
  print(cor(rma_data[shared_probes,1],gsm_original_data[shared_probes],method="spearman"))
  CEL_frma_profiles[[CEL_batch]] = frma_data
  CEL_rma_profiles[[CEL_batch]] = rma_data
}

save(geo_objs,gsm_objs,CEL_frma_profiles,CEL_rma_profiles,file="PA_database_profiles.RData")

# Run this code and resave
remove_cel_from_colnames<-function(cols){return(gsub(cols,pattern = ".CEL$",replace=''))}
for(i in 1:length(CEL_rma_profiles)){
  colnames(CEL_rma_profiles[[i]]) = remove_cel_from_colnames(colnames(CEL_rma_profiles[[i]]))
}
for(i in 1:length(CEL_frma_profiles)){
  colnames(CEL_frma_profiles[[i]]) = remove_cel_from_colnames(colnames(CEL_frma_profiles[[i]]))
}
gse_matrices = list()
for(nn in names(geo_objs)){
  ll = geo_objs[[nn]]
  for(nn2 in names(ll)){
    gse_matrices[[paste(nn,nn2,sep=";")]] = exprs(ll[[nn2]])
  }
}
length(gse_matrices)
sapply(gse_matrices,dim)

# Get platform data
platforms = unique(metadata$platform_id)
gpl_objs = list()
for(platform in platforms){
  try({
    gpl_objs[[platform]] = getGEO(platform,destdir = GEO_destdir)
  })
}
gpl_tables = lapply(gpl_objs,Table)
has_entrez = sapply(gpl_tables,function(x)any(grepl("entrez gene",colnames(x),ignore.case = T)))
has_gene = sapply(gpl_tables,function(x)any(grepl("gene",colnames(x),ignore.case = T)))
num_gene_cols = sapply(gpl_tables,function(x)sum(grepl("gene",colnames(x),ignore.case = T)))
num_probes = sapply(gpl_tables,nrow)

save(gse_matrices,gsm_objs,gpl_tables,CEL_frma_profiles,CEL_rma_profiles,file="PA_database_profiles.RData")

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# Use the GPL data to map probes to entrez ids

# Step 1: for each row in the gpl data table try to map the probe
# to the relevant entrez gene ids.
# Here we use the bioconductor org.Hs.eg.db package to map the following
# ids into entrez genes: genebank, gene symbol, gene name, ensembl gene, and unigene.
# The relevant columns in the gpl matrix are those that contain the information
# of the ids above or a direct mapping to entrez genes.
# For each row we go over all relevant columns, map them to entrez ids
# and keep the union.
converters = list()
# note that this genebank mapping contains all refseq ids
converters[["GB_ACC"]] = as.list(org.Hs.egACCNUM2EG)
converters[["GB_LIST"]] = converters[["GB_ACC"]]
converters[["symbol"]] = as.list(org.Hs.egSYMBOL2EG)
converters[["gene"]] = converters[["symbol"]]
converters[["ensembl"]] = as.list(org.Hs.egENSEMBL2EG)
converters[["unigene"]] = as.list(org.Hs.egUNIGENE2EG)

# correct the version issue in the GB_ACC columns
all_gpls = names(gpl_objs)
for (gpl in all_gpls){
  curr_col = which(colnames(gpl_objs[[gpl]])=="GB_ACC")
  if(length(curr_col)>0){
    gpl_objs[[gpl]][,"GB_OLD_ACC"] = gpl_objs[[gpl]][,"GB_ACC"]
    gpl_objs[[gpl]][,"GB_ACC"] = gsub(gpl_objs[[gpl]][,"GB_ACC"],pattern="\\.\\d+$",perl=T,replace="")
  }
}

sapply(gpl_tables,dim)
gpl_mappings_to_entrez = lapply(gpl_tables,map_gpl_rows_to_entrez,converters=converters)
all(names(gpl_mappings_to_entrez)==names(gpl_tables))
sapply(gpl_mappings_to_entrez,function(x)sum(sapply(x,is.null)))
sapply(gpl_mappings_to_entrez,function(x)sum(sapply(is.na(x),any)))
sapply(gpl_mappings_to_entrez,length) == sapply(gpl_tables,nrow)
sapply(gpl_mappings_to_entrez,function(x)sum(sapply(x,get_num_probe_genes)==1))
gpl_mappings_to_entrez[["GPL6947"]]

# # Possible improvement
# # Borrow data across the illumina platforms using the genebank data
# # We observed that in the illumina platforms some rows are mapped only using the 
# # genbank ids. On the other hand, these are not consistent across the platforms.
# # We therefore obtain a unified mapping of genebank ids to entrez
# # and then use it to fill gaps.
# # Practically it added a few hundrend annotations in each GPL
# manufacts = c()
# for (gpl in all_gpls){
#   gpl_obj = gpl_objs[[gpl]]
#   manufacts[gpl] = ""
#   try({manufacts[gpl] = Meta(gpl_obj)$manufacturer})
# }
# illu_gpls = names(manufacts[grepl("Illumin",manufacts)])
# all_gbs_to_entrez = list()
# for(gpl in illu_gpls){
#   if(!is.element("GB_ACC",set=colnames(gpl_tables[[gpl]]))){next}
#   curr_gbs = as.character(gpl_tables[[gpl]][,"GB_ACC"])
#   curr_ez = gpl_mappings_to_entrez[[gpl]]
#   for (j in 1:length(curr_gbs)){
#     all_gbs_to_entrez[[curr_gbs[j]]] = union(all_gbs_to_entrez[[curr_gbs[j]]],curr_ez[[j]])
#   }
# }
# table(sapply(all_gbs_to_entrez,length))
# count_corrections1 = 0
# count_corrections2 = 0
# for(gpl in illu_gpls){
#   if(!is.element("GB_ACC",set=colnames(gpl_tables[[gpl]]))){next}
#   curr_gbs = gpl_tables[[gpl]][,"GB_ACC"]
#   curr_ez = gpl_mappings_to_entrez[[gpl]]
#   for (j in 1:length(curr_gbs)){
#     new_ezs = union(curr_ez[[j]],all_gbs_to_entrez[[curr_gbs[j]]])
#     new_ezs = new_ezs[!is.na(new_ezs)]
#     if(length(new_ezs)>length(curr_ez[[j]])){
#       count_corrections1 = count_corrections1+1
#       if(length(curr_ez[[j]])==0){count_corrections2 = count_corrections2+1}
#       print(paste(c("#####",curr_ez[[j]]),collapse=';'))
#       print(all_gbs_to_entrez[[curr_gbs[j]]])
#       curr_ez[[j]] = new_ezs
#     }
#   }
#   gpl_mappings_to_entrez[[gpl]] = curr_ez
# }

# Step 3: use the correct rows from the GPL table as the probe ids.
# This is done by downloading a sample GSM profile for each gpl.
# We then go over the gpl table and seek the correct column that fits the ids
# in the names of the rows in the downloaded profile
gsms_and_platforms = as.matrix(metadata[,c("GSM","platform_id")])
for(gpl in all_gpls){
  curr_gsms = gsms_and_platforms[gsms_and_platforms[,2]==gpl,1]
  sample_gsm = curr_gsms[1]
  gsm_obj = Table(getGEO(sample_gsm))
  if(length(gsm_obj)==0){next}
  v = as.numeric(gsm_obj[,"VALUE"])
  names(v) = gsm_obj[,1]
  curr_table = gpl_tables[[gpl]]
  num_appearing = c()
  num_missing = c()
  for(j in 1:ncol(curr_table)){
    num_appearing[j] = length(intersect(names(v),curr_table[,j]))
    num_missing[j] = length(setdiff(names(v),curr_table[,j]))
  }
  curr_col = which(num_appearing==max(num_appearing))[1]
  print(num_appearing)
  #if(length(curr_cols)>1){
  #       problematic_rows = which(!apply(curr_table[,curr_cols],1,function(x)all(x==x[1])))
  #       curr_table[problematic_rows,curr_cols]
  #}
  names(gpl_mappings_to_entrez[[gpl]]) = curr_table[,curr_col]
}

# Step 4: we now have a mapping of probes to entrez genes,
# get the reverse mapping
gpl_mappings_entrez2probes = list()
for(gpl in all_gpls){
  gpl_mappings_entrez2probes[[gpl]] = reverse_mapping_list(gpl_mappings_to_entrez[[gpl]])
}
sapply(gpl_mappings_entrez2probes,length)

# Step 5: take unique mappings
gpl2probe2is_single_gene = lapply(gpl_mappings_to_entrez,function(x)sapply(x,get_num_probe_genes)==1)
gpl2probe2num_genes = lapply(gpl_mappings_to_entrez,function(x)sapply(x,get_num_probe_genes))
# percentage of probes that were mapped into a single gene
sapply(gpl2probe2is_single_gene,sum)/sapply(gpl2probe2num_genes,function(x)sum(x>0))
gpl_unique_mappings = list()
for(gpl in names(gpl_mappings_to_entrez)){
  l = gpl_mappings_to_entrez[[gpl]]
  newl = rep(NA,length(l))
  names(newl) = names(l)
  newl[gpl2probe2is_single_gene[[gpl]]] = unlist(gpl_mappings_to_entrez[gpl2probe2is_single_gene[[gpl]]])
  gpl_unique_mappings[[gpl]] = newl
}
sapply(gpl_unique_mappings,function(x)length(unique(x)))

save(gpl_mappings_to_entrez,gpl_mappings_entrez2probes,gpl_unique_mappings,file='gpl_mappings_to_entrez.RData')

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

# Tests: the data we need from gsm objects
gsm = gsm_objs[[1]]
gsm_meta = Meta(gsm)
platform = gsm_meta$platform_id
gsm_data = Table(gsm)

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# Get additional data not given in GEO directly
# load the database
load('gpl_mappings_to_entrez.RData')
load("PA_database_profiles.RData")

table(metadata$Type)
genomic_gsms = rownames(metadata)[metadata$Type=="genomic"]
SRA_mRNA_gsms = rownames(metadata)[metadata$Type=="SRA"]
mRNA_microarray_gsms = rownames(metadata)[metadata$Type=="RNA"]
no_gsm_obj = names(gsm_objs)[sapply(gsm_objs,class)=="character"]
metadata[no_gsm_obj,1:5]
for(gsm in no_gsm_obj){
  gsm_objs[[gsm]] = ""
  try({gsm_objs[[gsm]] = getGEO(gsm,destdir = GEO_destdir)})
  ind = length(gsm_objs)
  if (ind %% 10 ==0 ){print(ind)}
}
gsm_data_size = sapply(gsm_objs,get_gsm_data_size)
table(gsm_data_size)
table(gsm_data_size[SRA_mRNA_gsms])
gsm_no_data = names(which(gsm_data_size==0))
table(metadata[gsm_no_data,]$Type)

## Look at the genomic profiles
#gsm_genomic = genomic_gsms[1]
#gsm_objs[[gsm_genomic]]

# Get SRA data
# install recount
# source("https://bioconductor.org/biocLite.R")
# biocLite("recount")
# # code that we need:
# ## Find a project of interest
# project_info <- abstract_search('GSE32465')
# ## Download the gene-level RangedSummarizedExperiment data
# download_study(project_info$project)
# ## Load the data
# load(file.path(project_info$project, 'rse_gene.Rdata'))
# ## View GEO ids
# colData(rse_gene)$geo_accession
# dim(rse_gene)
# rse <- scale_counts(rse_gene)
# rse_fpkm = fpkm(rse)
# # For mapping SRXs to SRPs
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("SRAdb")
# biocLite("DESeq2")

library(recount); library(SRAdb);library(DESeq2)
sqlfile <- getSRAdbFile()
sra_connection <- dbConnect(SQLite(), "SRAmetadb.sqlite")
sra_has_profile_supp = c()
for(sra_gsm in SRA_mRNA_gsms){
  curr_obj = gsm_objs[[sra_gsm]]
  curr_obj_meta = Meta(curr_obj)
  supp_files = names(curr_obj_meta)[grepl(names(curr_obj_meta),pattern="supplementary_file")]
  curr_supp_files = curr_obj_meta[supp_files]
  sra_has_profile_supp[sra_gsm] = any(grepl("pkm",unlist(curr_supp_files)))
  curr_gse = as.character(metadata[sra_gsm,"GSE"])
}

sra_has_profile_recount = c()
covered_srps = c()
for(sra_gsm in SRA_mRNA_gsms){
  curr_obj = gsm_objs[[sra_gsm]]
  curr_obj_meta = Meta(curr_obj)
  curr_relation = curr_obj_meta$relation
  sra_relation = curr_relation[grepl("SRA:",curr_relation)]
  srx_id = strsplit(sra_relation,split="term=")[[1]][2]
  conversion <- sraConvert(srx_id, sra_con = sra_connection)
  sra_studies = as.character(conversion[,"study"])
  
  for(sra_study in sra_studies){
    if(is.element(sra_study,set=names(covered_srps))){
      sra_has_profile_recount[sra_gsm]=covered_srps[sra_study]
      break
    }
    covered_srps[sra_study] = F
    try({
      download_study(sra_study)
      sra_has_profile_recount[sra_gsm] = T
      covered_srps[sra_study]=T
    })
  }
}

table(sra_has_profile_recount)
table(sra_has_profile_supp)
table(sra_has_profile_recount | sra_has_profile_supp)

# Get the data from recount, keep it in the GSE matrices
# Map recount's rows into entrez and add it to the gpls
recount_covered_srps = names
library(org.Hs.eg.db)
ensembl2entrez <- as.list(org.Hs.egENSEMBL2EG)
all_entrez = unique(unlist(ensembl2entrez))
table(sapply(ensembl2entrez,length))
recount_entrez2rows = NULL
for(srp in recount_covered_srps){
  download_study(srp)
  load(file.path(srp, 'rse_gene.Rdata'))
  curr_gsms = colData(rse_gene)$geo_accession
  dim(rse_gene)
  rse <- scale_counts(rse_gene)
  dds <- DESeqDataSet(rse,~1)
  rse_fpkm = fpkm(dds)
  colnames(rse_fpkm) = curr_gsms
  if(is.null(recount_entrez2rows)){
    curr_ensembl_ids = rownames(rse_fpkm)
    curr_ensemble_gene = sapply(curr_ensembl_ids,function(x)strsplit(x,split='\\.')[[1]][1])
    gene2enseble_rows = reverse_mapping_list(curr_ensemble_gene)
    entrez2enseble_rows = list()
    for(e in all_entrez){entrez2enseble_rows[[e]]=""}
    for(j in 1:length(gene2enseble_rows)){
      rows = gene2enseble_rows[[j]]
      ensmbl = names(gene2enseble_rows)[j]
      es = ensembl2entrez[[ensmbl]]
      for(e in es){entrez2enseble_rows[[e]] = c(entrez2enseble_rows[[e]],rows)}
    }
    entrez2enseble_rows = lapply(entrez2enseble_rows,function(x)unique(setdiff(x,"")))
    entrez2enseble_rows = entrez2enseble_rows[sapply(entrez2enseble_rows,length)>0]
    recount_entrez2rows = entrez2enseble_rows
  }
  gse_matrices[[srp]] = rse_fpkm
}
sapply(gse_matrices[recount_covered_srps],dim)
gpl_mappings_entrez2probes[["recount"]] = recount_entrez2rows
save(gpl_mappings_to_entrez,gpl_mappings_entrez2probes,gpl_unique_mappings,file='gpl_mappings_to_entrez.RData')
save(gse_matrices,gsm_objs,gpl_tables,CEL_frma_profiles,CEL_rma_profiles,file="PA_database_profiles.RData")

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
#!!!!!!!!!!!!!!!! THIS CODE IS MOVING GRADUALLY TO SEPARATE SCRIPTS !!!!!!!!!!!!!!!!!!!!!!!


# The analysis starts here
# load the gene expression database
load('gpl_mappings_to_entrez.RData')
load("PA_database_profiles.RData")

# Set constants for the analysis
EXPERIMENT_TYPE = "longterm" # acute, long term
OMIC_TYPE = "mRNA" # mRNA, methylation
EXERCISE_TYPE = "both" # endurance, resistance, both, or other (yoga)
MIN_NROW = 5000 # minimal number of rows in a datamatrix of a specific dataset
OUT_FILE_ACUTE = "PADB_univariate_results_and_preprocessed_data.RData"
OUT_FILE_LONGTERM = "PADB_univariate_results_and_preprocessed_data_longterm.RData"
OUT_FILE = OUT_FILE_LONGTERM
ANALYSIS_OUT_FILE_ACUTE = "PADB_statistical_analysis_results_acute.RData"
ANALYSIS_OUT_FILE_LONGTERM = "PADB_statistical_analysis_results_longterm.RData"

# Get the current metadata
curr_rows = rep(F,nrow(metadata))
names(curr_rows) = rownames(metadata)
# Take all rows of the selected experiment type
if(EXPERIMENT_TYPE=="acute"){
  analysis_dataset_ids = dataset_ids_acute
  time_col = which(grepl("acute", colnames(metadata),ignore.case = T) & 
                     grepl("standardized", colnames(metadata),ignore.case = T))
  curr_rows[!is.na(metadata[,time_col]) & metadata[,time_col] != ""] = T
}
if(EXPERIMENT_TYPE!="acute"){
  analysis_dataset_ids = dataset_ids_longterm
  time_col = which(grepl("long", colnames(metadata),ignore.case = T) & 
                     grepl("standardized", colnames(metadata),ignore.case = T))
  curr_rows[!is.na(metadata[,time_col]) & metadata[,time_col] != ""] = T
}
# Of the selected rows use the selected omics and exercise types
if(OMIC_TYPE == "mRNA"){
  curr_rows = curr_rows & (metadata[,3] == "RNA" | metadata[,3] == "SRA") 
}

# Try this:
analysis_dataset_ids = dataset_ids_acute

# Select the sample set for this analysis
analysis_samples = rownames(metadata)[curr_rows]

# Get t-test p-values and statistics
# Load existing data
load(OUT_FILE) # Acute data
# OR recreate
dataset2preprocessed_data = list()
dataset2univar_pvals = list()
for (dataset in unique(analysis_dataset_ids[analysis_samples])){
  assign("last.warning", NULL, envir = baseenv())
  dataset_samples = names(which(analysis_dataset_ids==dataset))
  if(length(dataset_samples)==0){next}
  # cut by samples that have time points and are not repeats
  dataset_samples = dataset_samples[metadata[dataset_samples,time_col] != ""]
  dataset_samples = dataset_samples[metadata[dataset_samples,"Replicate.info"] == ""]
  if(length(dataset_samples)==0){next}
  curr_platforms = as.character(metadata[dataset_samples,"platform_id"])
  platform = curr_platforms[1]
  if(metadata[dataset_samples[1],"Type"]=="SRA"){
    platform = "recount"
  }
  #if(platform!="recount"){next}
  entry_name = dataset
  print(entry_name)
  if(is.element(entry_name,set=names(dataset2univar_pvals))){next}
  
  # Filter the current sample set:
  curr_samples = dataset_samples
  curr_samples = curr_samples[sample2replicate_info[curr_samples] == ""]
  curr_samples = curr_samples[sample2subject[curr_samples]!=""]
  curr_times = as.character(metadata[curr_samples,time_col])
  curr_samples = curr_samples[curr_times !=""]
  
  # Order: fRMA, RMA, GSE object, GSM profiles
  frma_mat = get_data_matrix_from_matrix_lists(CEL_frma_profiles,curr_samples)
  rma_mat = get_data_matrix_from_matrix_lists(CEL_rma_profiles,curr_samples)
  gse_mat = get_data_matrix_from_matrix_lists(gse_matrices,curr_samples)
  print(dim(gse_mat))
  
  data_matrix = frma_mat
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(curr_samples)){
    print("Not enough samples in fRMA object, use RMA:")
    print(dataset)
    data_matrix = rma_mat
  }
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(curr_samples)){
    print("Not enough samples in RMA object, use GSEs:")
    print(dataset)
    data_matrix = gse_mat
  }
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(curr_samples)){
    print("Not enough samples in GSE object:")
    print(dataset)
    print("Skipping for now, address later")
    next
  }
  
  na_rows = apply(is.na(data_matrix),1,any) | apply(is.nan(data_matrix),1,any)
  print(table(na_rows))
  data_matrix = data_matrix[!na_rows,]
  
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(curr_samples)){
    print("Data matrix has too many rows with NA or NaN values:")
    print(dataset)
    print("Skipping for now, address later")
    next
  }

  #######################################
  # Analyze the data matrix
  # TODO before getting the stats below:
  # 1. Check if the data was logged
  if(max(data_matrix,na.rm=T)>50){
    data_matrix = log(data_matrix,base=2)
    data_matrix[is.na(data_matrix)] = min(data_matrix,na.rm=T)
  }
  par(mfrow=c(1,1))
  hist(data_matrix,breaks=100)
  # 2. Check if there are negative values
  # 4. check exercise groups
  ######################################
  
  # Reorder the columns
  curr_samples = curr_samples[order(metadata[curr_samples,time_col])]
  curr_times = as.numeric(as.character(metadata[curr_samples,time_col]))
  curr_subjects = factor(sample2subject[curr_samples])
  curr_ages = factor(sample2age[curr_samples],ordered=T)
  curr_sex = factor(sample2sex[curr_samples])
  numeric_ages = as.numeric(as.character(curr_ages))
  
  if(!any(is.na(numeric_ages))){curr_ages = numeric_ages}
  if(length(unique(curr_ages))<2){curr_ages=NULL}
  if(length(unique(curr_sex))<2){curr_sex=NULL}
  
  if(sd(curr_times)==0){
    print("Dataset has only one time point:")
    print(dataset)
    print("Skipping")
    next
  }
  
  if(nrow(data_matrix) < MIN_NROW){
    print("Data matrix has less rows than the minimal number specified:")
    print(dataset)
    print("Skipping")
    next
  }
  
  if(length(warnings())>0){
    print("Got some warnings, breaking. Please debug:")
   # break
  }
  
  print("Adding a new dataset to the preprocessed data containers:")
  print(paste("****************",entry_name))
  #next
  
  #if(length(unique(curr_sex))<2){curr_sex=NULL}else{break}
  #next
  #if(length(unique(curr_times))>2){break}else{next}
  
  # # Skipped the probe p-values for now - to save some running time
  # curr_paired_tests = apply(data_matrix,1,run_anova_on_time_age_sex_and_subject,
  #                         subjs=curr_subjects,timev=curr_times,sexv=curr_sex,agev=curr_ages)
  # probe_paired_test_pvals = sapply(curr_paired_tests,function(x)x$pval)
  # hist(probe_paired_test_pvals,main=dataset)
  #sapply(apply(data_matrix[1:10,],1,run_paired_ttest,subjs=curr_subjects,timev=curr_times),function(x)x$p.value)
  probe_fold_changes = get_fold_changes_vs_baseline(data_matrix,curr_subjects,as.numeric(curr_times))
  # Transform to genes (entrez or symbols) and rerun
  genes_data_matrix_obj = transform_matrix_into_genes(data_matrix,gpl = platform,gpl_mappings_entrez2probes)
  genes_data_matrix = genes_data_matrix_obj$entrez_mat
  # exclude genes with NAs
  if(sum(genes_data_matrix_obj$entrez_mat_na_stats[["row NA counts"]]==0)>10000){
    genes_data_matrix = genes_data_matrix[genes_data_matrix_obj$entrez_mat_na_stats[["row NA counts"]]==0,]
  }
  print(table(is.na(genes_data_matrix)))
  # curr_paired_tests = apply(genes_data_matrix,1,run_anova_on_time_age_sex_and_subject,
  #                           subjs=curr_subjects,timev=curr_times,sexv=curr_sex,agev=curr_ages)
  # gene_paired_test_pvals = sapply(curr_paired_tests,function(x)x$pval)
  # par(mfrow=c(1,3))
  # hist(gene_paired_test_pvals)
  gene_paired_tests_lme4_pvals = apply(genes_data_matrix,1,run_lmer4_anova,
                             subjs=curr_subjects,timev=curr_times,sexv=curr_sex,agev=curr_ages)
  hist(gene_paired_tests_lme4_pvals)
  # plot(gene_paired_tests_lme4_pvals,gene_paired_test_pvals)
  # Sys.sleep(30)
  #sapply(apply(data_matrix[1:10,],1,run_paired_ttest,subjs=curr_subjects,timev=curr_times),function(x)x$p.value)
  gene_fold_changes = get_fold_changes_vs_baseline(genes_data_matrix,curr_subjects,as.numeric(curr_times))

  # add the results to the data containers
  dataset2univar_pvals[[entry_name]] = list()
  #dataset2univar_pvals[[entry_name]][["genes anova"]] = gene_paired_test_pvals
  dataset2univar_pvals[[entry_name]][["genes lme4 anova"]] = gene_paired_tests_lme4_pvals
  #dataset2univar_pvals[[entry_name]][["probes anova"]] = probe_paired_test_pvals
  dataset2preprocessed_data[[entry_name]] = list()
  dataset2preprocessed_data[[entry_name]][["probe_data"]] = data_matrix
  dataset2preprocessed_data[[entry_name]][["gene_data"]] = genes_data_matrix
  dataset2preprocessed_data[[entry_name]][["time"]] = curr_times
  dataset2preprocessed_data[[entry_name]][["subject"]] = curr_subjects
  dataset2preprocessed_data[[entry_name]][["gsms"]] = colnames(data_matrix)
  dataset2preprocessed_data[[entry_name]][["age"]] = curr_ages
  dataset2preprocessed_data[[entry_name]][["sex"]] = curr_sex
  dataset2preprocessed_data[[entry_name]][["probe_fold_changes"]] = probe_fold_changes
  dataset2preprocessed_data[[entry_name]][["gene_fold_changes"]] = gene_fold_changes
  
  # release memory and save
  rm(curr_paired_tests);rm(data_matrix);
  rm(genes_data_matrix);rm(genes_data_matrix_obj)
  gc()
  save(dataset2preprocessed_data,dataset2univar_pvals,file=OUT_FILE)
}

############################# ANALYSIS STARTS HERE ######################
load(OUT_FILE_LONGTERM)
load(OUT_FILE_ACUTE)

# 0. Preprocessing
# Dataset information
cols = names(dataset2univar_pvals)
arrs = strsplit(cols,split=';')
combined_gses = sapply(arrs,function(x)x[1])
is_endurance = sapply(arrs,function(x)any(grepl("endur",x,ignore.case = T)))
is_resistance = sapply(arrs,function(x)any(grepl("resis",x,ignore.case = T)|grepl("strength",x,ignore.case = T)))
dataset_is_endurance = c()
dataset_is_resistance = c()
for(dataset in names(dataset2preprocessed_data)){
  currmat = dataset2preprocessed_data[[dataset]]$gene_data
  currmeta = metadata[colnames(currmat),]
  training_program = as.character(currmeta$Training.program.during.experiment.type)
  print(dataset);print(unique(training_program))
  dataset_is_endurance[dataset]=F;dataset_is_resistance[dataset]=F
  if(any(grepl("endur",training_program,ignore.case = T))){dataset_is_endurance[dataset]=T}
  if(any(grepl("aerob",training_program,ignore.case = T))){dataset_is_endurance[dataset]=T}
  if(any(grepl("resistance",training_program,ignore.case = T))){dataset_is_resistance[dataset]=T}
  if(any(grepl("streng",training_program,ignore.case = T))){dataset_is_resistance[dataset]=T}
}
############ TODO::::: #############
# for longterm - need to be corrected
is_endurance = dataset_is_endurance
is_resistance = dataset_is_resistance
#################################
tissue = rep("muscle",length(is_endurance))
tissue[sapply(arrs,function(x)any(grepl("blood",x,ignore.case = T)))] = "blood"
tissue[sapply(arrs,function(x)any(grepl("PBMC",x,ignore.case = T)))] = "blood"
tissue[sapply(arrs,function(x)any(grepl("cytes",x,ignore.case = T)))] = "blood"
tissue[sapply(arrs,function(x)any(grepl("phil",x,ignore.case = T)))] = "blood"
cols[tissue!="muscle"]
names(tissue) = cols

exercise_type = is_resistance
names(exercise_type) = cols
exercise_type[is_resistance]="resistance"
exercise_type[is_endurance]="endurance"
exercise_type[is_endurance&is_resistance]="both"
exercise_type[!is_endurance&!is_resistance]="other"
table(exercise_type)

# Analyze the gene output: get the gene sets
gene_sets = lapply(dataset2preprocessed_data,function(x)rownames(x$gene_fold_changes))
sapply(gene_sets,length)
inters = gene_sets[[1]]
for(l in gene_sets){inters = intersect(inters,l)}
print(length(inters))

#########################################################################
# 1. REPLICABILITY: get the gene p-values
gene_combined_pvalues = c()
for(dataset in names(dataset2univar_pvals)){
  currp = dataset2univar_pvals[[dataset]][[1]]
  currp = currp[inters]
  gene_combined_pvalues = cbind(gene_combined_pvalues,currp)
}
par(mfrow=c(3,3))
apply(gene_combined_pvalues[,1:9],2,hist)
colnames(gene_combined_pvalues) = names(dataset2univar_pvals)
sort(apply(gene_combined_pvalues<0.05,1,sum),decreasing=T)[1:10]
gene_qvals = apply(gene_combined_pvalues,2,p.adjust,method='fdr')
# By analysis by tissue and exercise type
rep_scores = c()
for(tt in unique(tissue)){
  for(ee in unique(exercise_type)){
    inds = tissue==tt & exercise_type==ee
    if(sum(inds)<2){next}
    qmat = gene_qvals[,inds]
    # use SCREEN later
    rep_scores = cbind(rep_scores,rowSums(qmat<0.2)/ncol(qmat))
    colnames(rep_scores)[ncol(rep_scores)] = paste(tt,ee,sum(inds),sep=';')
  }
}
library(corrplot)
q_corrs = cor(gene_qvals)
colnames(q_corrs) = tissue
rownames(q_corrs) = exercise_type
par(mfrow=c(1,1))
corrplot(q_corrs,tl.cex=1,order='hclust')
cor(rep_scores)
rep_gene_sets = apply(rep_scores,2,function(x)names(x)[x>0.499]) # acute has more datasets
rep_gene_sets = apply(rep_scores,2,function(x)names(x)[x>0.6]) # long term has fewer
sapply(rep_gene_sets,length)

# GO enrichment
library(topGO)
set2go_results = list()
for (i in 1:length(rep_gene_sets)){
  geneUniverse <- rownames(rep_scores)
  genesOfInterest = rep_gene_sets[[i]]
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  myGOdata <- new("topGOdata", description="My project", ontology="BP", nodeSize=5,
                allGenes=geneList, annot = annFUN.org, mapping="org.Hs.eg.db", ID = "entrez")
  resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")
  allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = 
                     "resultFisher", ranksOf = "classicFisher", topNodes = length(score(resultFisher)))
  go_pvals = as.numeric(allRes[,ncol(allRes)])
  go_qvals = p.adjust(go_pvals,method='fdr')
  allRes = cbind(allRes,go_qvals)
  set2go_results[[i]] = allRes[go_qvals<0.1 & allRes$Annotated<1000,]
}

save(gene_combined_pvalues,rep_scores,gene_qvals,rep_gene_sets,set2go_results,file=ANALYSIS_OUT_FILE_ACUTE)

# # Pathway enrichment (reactome)
# # https://bioconductor.org/packages/devel/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html#pathway-enrichment-analysis
# library(ReactomePA)
# sapply(rep_gene_sets, length)
# de <- rep_gene_sets[[3]]
# # simple enrichment
# x <- enrichPathway(gene=rep_gene_sets[[3]],pvalueCutoff=1, readable=T,universe = rownames(rep_scores),minGSSize=3,qvalueCutoff = 1,pAdjustMethod="BH")
# barplot(x, showCategory=15)
# ## gsea
# #x <- gsePathway(sort(rep_scores[,i],decreasing=T), organism = "human", exponent = 1, nPerm = 1000,
# #                minGSSize = 3, maxGSSize = 500, pvalueCutoff = 1,
# #                pAdjustMethod = "none", verbose = TRUE, seed = FALSE, by = "DOSE")


#########################################################################

#########################################################################
# 2. Combined gene data
gene_subject_combined_ge_matrix = c()
sample_metadata = list()
for(dataset in names(dataset2preprocessed_data)){
  currmat = dataset2preprocessed_data[[dataset]]$gene_data
  currmat = currmat[inters,]
  curr_datasets = rep(dataset,ncol(currmat))
  names(curr_datasets) = colnames(currmat)
  gene_subject_combined_ge_matrix = cbind(gene_subject_combined_ge_matrix,currmat)
  curr_subjects = dataset2preprocessed_data[[dataset]]$subject
  sample_metadata$dataset = c(sample_metadata$dataset,curr_datasets)
  sample_metadata$subject = c(sample_metadata$subject,curr_subjects)
  sample_metadata$time = c(sample_metadata$time,dataset2preprocessed_data[[dataset]]$time)
  sample_metadata$tissue = c(sample_metadata$tissue,rep(tissue[dataset],ncol(currmat)))
  sample_metadata$exercise_type = c(sample_metadata$exercise_type,rep(exercise_type[dataset],ncol(currmat)))
}
dim(gene_subject_combined_ge_matrix)
library(preprocessCore)
gene_subject_combined_ge_matrix_q = normalize.quantiles.robust(gene_subject_combined_ge_matrix)
gene_data_pca = prcomp(t(gene_subject_combined_ge_matrix_q),scale=T,center=T,retx=T)
gene_data_PCs = data.frame(gene_data_pca$x[,1:200])
gene_data_PCs = cbind(gene_data_PCs,sample_metadata$tissue,sample_metadata$exercise_type)
run_example = leave_study_out(sample_metadata$time,gene_data_PCs,study_ids[colnames(gene_subject_combined_ge_matrix)],run_rf)
cor(t(sapply(run_example,colMeans)))
tmp = t(sapply(run_example,colMeans))
plot(tmp[,1],tmp[,2])
sapply(run_example,function(x)cor.test(x[,1],x[,2],method="spearman")$p.value)

#########################################################################
# 3. Get the gene fold changes
gene_fold_change_combine_mat = c()
subject_metadata = list()
sort(rowSums(gene_fold_change_combine_mat>1),decreasing=T)[1:10]
sort(rowSums(gene_fold_change_combine_mat< (-1)),decreasing=T)[1:10]

for(dataset in names(dataset2preprocessed_data)){
  currmat = dataset2preprocessed_data[[dataset]]$gene_fold_changes
  currmat = currmat[inters,]
  curr_times = as.numeric(sapply(colnames(currmat),get_time_from_subj_names))
  if(any(is.na(curr_times))){break}
  curr_cols = paste(dataset,colnames(currmat),sep=";")
  curr_gsms = colnames(dataset2preprocessed_data[[dataset]]$probe_data)
  curr_study = study_ids[curr_gsms[1]]
  curr_study = rep(curr_study,ncol(currmat))
  dataset_arr = strsplit(dataset,split=';')[[1]]
  curr_training = rep(dataset_arr[length(dataset_arr)],ncol(currmat))
  colnames(currmat) = curr_cols
  gene_fold_change_combine_mat = cbind(gene_fold_change_combine_mat,currmat)
  subject_metadata$study = c(subject_metadata$study,curr_study)
  subject_metadata$time = c(subject_metadata$time,curr_times)
  subject_metadata$tissue = c(subject_metadata$tissue,rep(tissue[dataset],ncol(currmat)))
  subject_metadata$exercise_type = c(subject_metadata$exercise_type,rep(exercise_type[dataset],ncol(currmat)))
}

# PCA + lm tests
gene_fold_change_pca = prcomp(t(gene_fold_change_combine_mat),scale=T,center=T,retx=T)
plot(gene_fold_change_pca)
gene_fold_change_PCs = data.frame(gene_fold_change_pca$x[,1:20])
gene_fold_change_PCs = cbind(gene_fold_change_PCs,subject_metadata$tissue,subject_metadata$exercise_type)
curr_t = subject_metadata$time
simple_lm = lm(curr_t~.,data=gene_fold_change_PCs)
summary(simple_lm)
library(randomForest)
simpleRF = randomForest(x = gene_fold_change_PCs,y=curr_t)
dummy_times = sample(curr_t)
run_example = leave_study_out(curr_t,gene_fold_change_PCs,subject_metadata$study,run_rf)
sapply(run_example,function(x)cor.test(x[,1],x[,2],method="spearman")$p.value)
cor(t(sapply(run_example,colMeans)))

