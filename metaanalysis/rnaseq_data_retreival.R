
rna_seq_datasets = c("GSE71972","GSE87749","GSE107934","GSE97084")

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
acute_metadata = read.xlsx2(file=metadata_file,sheetIndex=1)
longterm_metadata = read.xlsx2(file=metadata_file,sheetIndex=2)

# load human data for annotation
library(org.Hs.eg.db)
xx <- as.list(org.Hs.egENSEMBL2EG)
ensembl2entrez = xx
entrez2ensembl = as.list(org.Hs.egENSEMBL)
# remove nas
entrez2ensembl = lapply(entrez2ensembl,function(x)x[!is.na(x)])
entrez2ensembl = entrez2ensembl[!sapply(entrez2ensembl,length)==0]

# source("https://bioconductor.org/biocLite.R")
# biocLite("Homo.sapiens")
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

entrez_count_mat_to_fpkm <- function(m){
  # first get the ranges of the entrez genes in the matrix m
  gr = subset(entrez_gr,is.element(ENTREZID,set=rownames(m)))
  # subset the data to those with GRs - some may be lost here
  m = m[gr$ENTREZID,]
  # format the data in what DESeq expects
  se <- SummarizedExperiment(list(counts=m), colData=colnames(m))
  dds <- DESeqDataSet(se, ~ 1)
  rowRanges(dds) <- gr
  newm = fpkm(dds)
  newm = log(newm+1,base=2)
  return(newm)
}

rnaseq_matrices = list()

###########################
# GSE97084
###########################
gse = "GSE97084"
curr_meta = acute_metadata[acute_metadata$GSE==gse,]
if(nrow(curr_meta)==0){
  curr_meta = longterm_metadata[longterm_metadata$GSE==gse,]
}
count_matrix1 = read.delim("./raw_exp_data/GSE97084_GeneCount_raw.tsv",stringsAsFactors = F)
rownames(count_matrix1) = count_matrix1$GeneID
count_matrix1 = count_matrix1[,-c(1:6)]
count_matrix2 = read.delim("./raw_exp_data/GSE97084_GeneCount_raw_2.tsv",stringsAsFactors = F)
rownames(count_matrix2) = count_matrix2$GeneID
count_matrix2 = count_matrix2[,-c(1:6)]
count_matrix2 = count_matrix2[,-which(colnames(count_matrix2)=="X")]
all(rownames(count_matrix1)==rownames(count_matrix2))
count_matrix = cbind(count_matrix1,count_matrix2)
dim(count_matrix)
# fix the column name issue
ids1 = as.character(curr_meta$Subject.id)
times1 = as.character(curr_meta$Long.term..Standardized.Time..days.)

ids2 = colnames(count_matrix)
ids2 = gsub(ids2,pattern = "s_",replacement = "")
ids2 = gsub(ids2,pattern = ".180_GeneCount",replacement = "")
ids2 = gsub(ids2,pattern = "\\w$",replacement = "")

ids2_to_gsm = c()
for(id in unique(ids1)){
  inds1 = ids1==id
  inds2 = ids2==id
  ns1 = as.character(curr_meta[inds1,]$GSM)
  ns2 = colnames(count_matrix)[inds2]
  ns2 = ns2[order(ns2)]
  ns1 = ns1[order(times1[inds1])]
  ids2_to_gsm[ns2]=ns1
}
colnames(count_matrix) = ids2_to_gsm[colnames(count_matrix)]
count_matrix = probe2genes_conv(count_matrix,entrez2ensembl,f=sum)
# check NAs
fpkm_matrix = entrez_count_mat_to_fpkm(count_matrix)
rnaseq_matrices[["GSE97084"]] = fpkm_matrix
save(rnaseq_matrices,file="rnaseq_matrices.RData")

###########################
# GSE71972
###########################
gse = "GSE71972"
curr_meta = acute_metadata[acute_metadata$GSE==gse,]
if(nrow(curr_meta)==0){
  curr_meta = longterm_metadata[longterm_metadata$GSE==gse,]
}
gse_dir = paste("raw_exp_data/",gse,"_RAW/",sep="")
gse_files = list.files(gse_dir)
gse_files = gse_files[grepl("post_star",gse_files)]
gse_data = c()
for(f in gse_files){
  currd = read.table(paste(gse_dir,f,sep=""),stringsAsFactors = F)
  currx = currd[,2]
  names(currx) = currd[,1]
  curr_gsm = strsplit(f,split="_")[[1]][1]
  if(length(gse_data)>0){
    print(all(rownames(gse_data)==names(currx)))
  }
  gse_data = cbind(gse_data,currx)
  colnames(gse_data)[ncol(gse_data)]=curr_gsm
}
all_ens = unique(unlist(entrez2ensembl))
dim(gse_data)
gse_data = gse_data[intersect(rownames(gse_data),all_ens),]
curr_entrez2ensembl = lapply(entrez2ensembl, intersect,y=rownames(gse_data))
curr_entrez2ensembl = curr_entrez2ensembl[sapply(curr_entrez2ensembl,length)>0]
dim(gse_data)
count_matrix = probe2genes_conv(gse_data,curr_entrez2ensembl,f=sum)
fpkm_matrix = entrez_count_mat_to_fpkm(count_matrix)
rnaseq_matrices[["GSE71972"]] = fpkm_matrix
save(rnaseq_matrices,file="rnaseq_matrices.RData")

###########################
# GSE107934
###########################
gse="GSE107934"
curr_meta = acute_metadata[acute_metadata$GSE==gse,]
if(nrow(curr_meta)==0){
  curr_meta = longterm_metadata[longterm_metadata$GSE==gse,]
}
gse_dir = paste("raw_exp_data/",gse,"_RAW/",sep="")
gse_files = list.files(gse_dir)
gse_data = c()
for(f in gse_files){
  currd = read.table(paste(gse_dir,f,sep=""),stringsAsFactors = F)
  currx = currd[,2]
  names(currx) = currd[,1]
  curr_gsm = strsplit(f,split="_")[[1]][1]
  if(length(gse_data)>0){
    print(all(rownames(gse_data)==names(currx)))
  }
  gse_data = cbind(gse_data,currx)
  colnames(gse_data)[ncol(gse_data)]=curr_gsm
}
dim(gse_data)
gse_data = gse_data[intersect(rownames(gse_data),all_ens),]
curr_entrez2ensembl = lapply(entrez2ensembl, intersect,y=rownames(gse_data))
curr_entrez2ensembl = curr_entrez2ensembl[sapply(curr_entrez2ensembl,length)>0]
dim(gse_data)
count_matrix = probe2genes_conv(gse_data,curr_entrez2ensembl,f=sum)
fpkm_matrix = entrez_count_mat_to_fpkm(count_matrix)
dim(fpkm_matrix)
rnaseq_matrices[["GSE107934"]] = fpkm_matrix
save(rnaseq_matrices,file="rnaseq_matrices.RData")

###########################
# GSE87749
###########################
gse = "GSE87749"
curr_meta = acute_metadata[acute_metadata$GSE==gse,]
if(nrow(curr_meta)==0){
  curr_meta = longterm_metadata[longterm_metadata$GSE==gse,]
}
rownames(curr_meta) = as.character(curr_meta$GSM)
count_matrix = read.delim("./raw_exp_data/GSE87748_RNAseq_Readcoverage_v3.txt",
                          stringsAsFactors = F,row.names = 1)
count_matrix[1:5,1:5]

# fix the column name issue
ids1 = as.character(curr_meta$Subject.id)
ids1 = sapply(ids1,function(x)strsplit(x,split=" ")[[1]][2])
times1 = as.character(curr_meta$Acute..Standardized.Time..hours...1.is.baseline.)

cols2 = colnames(count_matrix)
cols2 = gsub(cols2,pattern = "^X",replacement = "")
ids2 = sapply(cols2,function(x)strsplit(x,split="_")[[1]][1])
times2 = sapply(cols2,function(x)strsplit(x,split="_")[[1]][3])
times2[times2=="Basal"]=-1
times2[times2=="240"]=4
names(times2) = colnames(count_matrix)

times2 = as.numeric(times2)
times1 = as.numeric(times1)

# Create a mapping from the data sample names to GSMs
ids2_to_gsm = c()
for(id in unique(ids1)){
  inds1 = ids1==id
  inds2 = ids2==id
  curr_t1 = times1[inds1]
  curr_t2 = times2[inds2]
  ns1 = as.character(curr_meta[inds1,]$GSM)
  ns2 = colnames(count_matrix)[inds2]
  ns2 = ns2[order(times2[inds2])]
  ns1 = ns1[order(times1[inds1])]
  # A sanity check: make sure that the mapped times from
  # both annotations map to the exact same results
  print(all(curr_meta[ns1,]$Acute..Standardized.Time..hours...1.is.baseline. ==
              times2[ns2]))
  ids2_to_gsm[ns2]=ns1
}
colnames(count_matrix) = ids2_to_gsm[colnames(count_matrix)]
count_matrix = probe2genes_conv(count_matrix,entrez2ensembl,f=sum)
dim(count_matrix)
table(rowSums(is.na(count_matrix))) # we see that NAs are in complete rows
count_matrix = count_matrix[!apply(is.na(count_matrix),1,any),]
dim(count_matrix)
fpkm_matrix = entrez_count_mat_to_fpkm(count_matrix)
dim(fpkm_matrix)
rnaseq_matrices[[gse]] = fpkm_matrix
save(rnaseq_matrices,file="rnaseq_matrices.RData")


