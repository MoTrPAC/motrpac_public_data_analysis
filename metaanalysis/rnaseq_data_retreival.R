
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
source("https://bioconductor.org/biocLite.R")
biocLite("Homo.sapiens")
library(GenomicRanges)
library(Homo.sapiens)
geneRanges <- function(db=Homo.sapiens, column="ENTREZID"){
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
}
entrez_gr = geneRanges(Homo.sapiens, column="ENTREZID")

count_mat_to_fpkm <- function(x,gr){
  se <- SummarizedExperiment(list(counts=m), colData=DataFrame(sample=1:4))
  dds <- DESeqDataSet(se, ~ 1)
  
}


# GSE97084
count_matrix = read.delim("./raw_exp_data/GSE97084_GeneCount_raw_2.tsv",stringsAsFactors = F)
rownames(count_matrix) = count_matrix$GeneID
count_matrix = count_matrix[,-c(1:6)]
count_matrix[1:5,1:5]
count_matrix = probe2genes_conv(count_matrix,entrez2ensembl,f=sum)

