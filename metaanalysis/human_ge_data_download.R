# This script downloads, collates, and normalizes raw expression profiles

# Paths
WD = "/Users/David/Desktop/MoTrPAC/project_release_feb_2018/data/"
SCRIPTS = "/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/"
GEO_destdir = paste(WD,"GEO/",sep="")
metadata_file = 'GEO_sample_metadata.xlsx'
raw_data_output_obj = 'human_ge_profiles_db.RData'

# Configurate the session
setwd(WD)
library('xlsx');library('GEOquery')
source(paste(SCRIPTS,'ge_download_preprocessing_helper_functions.R',sep=''))
# Run preprocessing of the CELs using frma
# source("https://bioconductor.org/biocLite.R")
# biocLite("frma")
# biocLite('hugene10stv1frmavecs')
# biocLite("oligo")
# biocLite("hugene.1.0.st.v1frmavecs")
# biocLite("hgu133plus2frmavecs")
# biocLite("annotate",dependencies=T)
# biocLite("hgu219frmavecs");biocLite("pd.hg.u219");biocLite("hgu219.db")
# COMMENT: frma failed on hgu219frmavecs - no such package 
library('frma'); library('affy');library('oligo');library('hugene.1.0.st.v1frmavecs')
library("hgu133plus2frmavecs")
library("hgu219frmavecs")
# # Get SRA data
# # install recount
# source("https://bioconductor.org/biocLite.R")
# biocLite("recount")
# # code that we need:
# ## Find a project of interest
# library(recount)
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
library(org.Hs.eg.db)

acute_metadata = read.xlsx2(file=metadata_file,sheetIndex=1)
longterm_metadata = read.xlsx2(file=metadata_file,sheetIndex=2)

# Solve duplications in the GSMs of acute (happened in one dataset with several GEO ids)
gsm_duplications = names(which(table(acute_metadata[,1])>1))
to_keep = rep(T,nrow(acute_metadata))
for(curr_gsm in gsm_duplications){
  curr_ids = which(acute_metadata[,1]==curr_gsm)
  to_keep[curr_ids[-1]]=F
}
acute_metadata = acute_metadata[to_keep,]

# Merge first three columns and keep in a table
metadata = as.matrix(unique(rbind(acute_metadata[,1:4],longterm_metadata[,1:4])))
rownames(metadata) = metadata[,1]

# Download and save the expression data
# Get the data by taking the series matrices
all_gse_ids = unique(metadata[,2])
all_gse_ids = all_gse_ids[grepl("^GSE",all_gse_ids,perl=T) | grepl("^GDS",all_gse_ids,perl=T)]
# download the pheno and meta data
geo_objs = list()
gse_fails = c()
for (gse in all_gse_ids){
  try({
    geo_objs[[gse]] = getGEO(gse,destdir = GEO_destdir,GSEMatrix = T,getGPL = F)
  })
}
# failed downloads:
failed_gses = setdiff(all_gse_ids,names(geo_objs))
# Test
length(geo_objs)
table(unlist(sapply(geo_objs,function(x)unlist(sapply(x,class))))) # all should be ExpressionSet

####### Go directly to the GSMs
gsms = metadata[,1]
sra_gsms = union(acute_metadata[grepl("SRA",acute_metadata[,3]),1],
                 longterm_metadata[grepl("SRA",longterm_metadata[,3]),1])
# download the database
gsm_objs = list()
all_gsms_in_destdir = list.files(GEO_destdir)
for(gsm in gsms){
  if(is.element(gsm,set=names(gsm_objs))){next}
	gsm_objs[[gsm]] = ""
	try({gsm_objs[[gsm]] = getGEO(gsm,destdir = GEO_destdir)})
	ind = length(gsm_objs)
	if (ind %% 100 ==0 ){print(ind)}
}
length(gsm_objs)
failed_gsms = names(which(sapply(gsm_objs,class)=="character"))

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
  for (i in 1:20){
    if(download_success){break}
    try({
      download.file(supplementary_CEL_files[gsm],destfile = dstfile)
      download_success = T
    })
  }
}

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
  if(!frma_success){
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

# Remove the 'CEL' suffix from the names
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
sapply(gse_matrices[grepl("GSE47969",names(gse_matrices))],dim)
sapply(gse_matrices[grepl("GSE58250",names(gse_matrices))],dim)
sapply(gse_matrices[grepl("GSE19062",names(gse_matrices))],dim)
sapply(gse_matrices[grepl("GSE117525",names(gse_matrices))],dim)
sapply(gse_matrices,dim)

# Get platform data
platforms = unique(metadata[,'platform_id'])
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
save(gse_matrices,gsm_objs,gpl_tables,CEL_frma_profiles,CEL_rma_profiles,file=raw_data_output_obj)

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
# all(names(gpl_mappings_to_entrez)==names(gpl_tables))
# sapply(gpl_mappings_to_entrez,function(x)sum(sapply(x,is.null)))
# sapply(gpl_mappings_to_entrez,function(x)sum(sapply(is.na(x),any)))
# sapply(gpl_mappings_to_entrez,length) == sapply(gpl_tables,nrow)
# sapply(gpl_mappings_to_entrez,function(x)sum(sapply(x,get_num_probe_genes)==1))
# gpl_mappings_to_entrez[["GPL6947"]]

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

# Get additional data not given in GEO directly
genomic_gsms = rownames(metadata)[metadata[,'Type']=="genomic"]
SRA_mRNA_gsms = rownames(metadata)[metadata[,'Type']=="SRA"]
mRNA_microarray_gsms = rownames(metadata)[metadata[,'Type']=="RNA"]
no_gsm_obj = names(gsm_objs)[sapply(gsm_objs,class)=="character"]
no_gsm_obj = no_gsm_obj[no_gsm_obj!=""]
for(gsm in no_gsm_obj){
  gsm_objs[[gsm]] = ""
  try({gsm_objs[[gsm]] = getGEO(gsm,destdir = GEO_destdir)})
  ind = length(gsm_objs)
  if (ind %% 10 ==0 ){print(ind)}
}
gsm_data_size = sapply(gsm_objs,get_gsm_data_size)
gsm_no_data = setdiff(names(which(gsm_data_size==0)),"")
table(metadata[gsm_no_data,"Type"])

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
recount_covered_srps = names(which(covered_srps))
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
save(gse_matrices,gsm_objs,gpl_tables,CEL_frma_profiles,CEL_rma_profiles,file=raw_data_output_obj)

# Manual correction: for dataset GSE117525 exclude the samples from the RMA data
CEL_rma_profiles = CEL_rma_profiles[!grepl("GSE117525",names(CEL_rma_profiles))]
save(gse_matrices,gsm_objs,gpl_tables,CEL_frma_profiles,CEL_rma_profiles,file=raw_data_output_obj)

# Manually added datasets
load("gpl_mappings_to_entrez.RData")
which(sapply(gpl_mappings_entrez2probes,length)==0)
# "GPL20880"
pl = getGEO("GPL20880",destdir = GEO_destdir)
pl = pl@dataTable@table
entrez2probe = pl[,1]
names(entrez2probe) = pl[,3]
entrez2probe = as.list(entrez2probe)
gpl_mappings_entrez2probes[["GPL20880"]] = entrez2probe
save(gpl_mappings_to_entrez,gpl_mappings_entrez2probes,gpl_unique_mappings,file='gpl_mappings_to_entrez.RData')
# GPL14613 - miRNA - skipped
# GPL16791 - sequencing - skipped
# GPL16770 - miRNA - skipped
# GPL10850 - miRNA - skipped
# GPL7731 - miRNA - skipped
# GPL11154 - Sequencing
# GPL8227 - miRNA








