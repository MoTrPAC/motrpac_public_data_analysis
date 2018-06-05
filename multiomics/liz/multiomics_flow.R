setwd("/Users/David/Desktop/multiomics/liz/arabidopsis")
# Examples for differential analysis
library(ggplot2);library(grid);library(lme4);library(DESeq2)

# Get sample data on gene level
gene_counts = read.csv("gene_counts.txt", sep="",
                       stringsAsFactors = TRUE,row.names = NULL)
rownames(gene_counts)=gene_counts[,2];gene_counts=gene_counts[,-c(1:2)]
row_metadata = as.matrix(gene_counts[,1:5])
chrs = sapply(row_metadata[,1],function(x)unique(strsplit(x,split=';')[[1]]))
chrs = paste("chr",chrs,sep='')
starts = sapply(row_metadata[,2],function(x)as.numeric(strsplit(x,split=';')[[1]]))
min_starts = sapply(starts,min)
ends = sapply(row_metadata[,3],function(x)as.numeric(strsplit(x,split=';')[[1]]))
sizes = sapply(row_metadata[,5],function(x)as.numeric(strsplit(x,split=';')[[1]]))
gr0 <- GRanges(Rle(chrs),IRanges(min_starts, width=sizes))
count_data = as.matrix(gene_counts[,-c(1:5)])
# Transform to FPKM
se <- SummarizedExperiment(list(counts=count_data), colData=DataFrame(sample=1:ncol(count_data)))
dds <- DESeqDataSet(se, ~ 1)
rowRanges(dds) <- gr0
gene_fpkm = fpkm(dds);rownames(gene_fpkm)=rownames(count_data)
colnames(gene_fpkm) = gsub(colnames(gene_fpkm),pattern=".all_trimmed.hisat2.bam",replace='')
colnames(gene_fpkm) = gsub(colnames(gene_fpkm),pattern="^X",replace='')
# get the metadata
sample_data = t(sapply(colnames(gene_fpkm),function(x)strsplit(x,split="_|h")[[1]]))
sample_data[sample_data[,2]==""|sample_data[,2]=="arvest",2] = "pre"
sample_data[sample_data[,1]==""|sample_data[,2]=="pre",1] = "-2"
colnames(sample_data) = c("time","type","repeat")
sample_data = data.frame(sample_data)
sample_data$time = as.numeric(as.character(sample_data$time))
sort(sample_data[,1])
# For our tests get a small dataset: log and filter by expression level and variance
gene_fpkm[gene_fpkm==0] = min(gene_fpkm[gene_fpkm>0])
gene_fpkm = log(gene_fpkm,base=2)
hist(gene_fpkm)
sds = apply(gene_fpkm,1,sd)
gene_fpkm = gene_fpkm[sds >= sort(sds,decreasing = T)[1000],]
############## Mixed effect models using lme4 ####################
# d1 is the metadata
# x is the expression value
# frm0 is the null hypothesis formula
# frm1 is the alt hypothesis formula
get_mixed_effect_model_time_apprx_p<-function(x,frm1,frm0,d1,statistic="Pr(>Chisq)",...){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  m1 = lmer(frm1,d,REML=F)
  m0 = lmer(frm0,d,REML=F)
  return(as.numeric(anova(m1,m0)[2,statistic]))
}
get_mixed_effect_model_time_apprx_stat_diff<-function(x,frm1,frm0,d1,statistic="Pr(>Chisq)",...){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  m1 = lmer(frm1,d,REML=F)
  m0 = lmer(frm0,d,REML=F)
  an =anova(m1,m0)[,statistic]
  return(an[1]-an[2])
}
pvals = apply(gene_fpkm,1,get_mixed_effect_model_time_apprx_p,
  frm1=x~factor(time)+factor(type)+1|repeat.,frm0=x~factor(type)+1|repeat.,
  d1 = sample_data)
hist(pvals)
table(p.adjust(pvals)<0.1)


