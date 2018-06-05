#setwd("/Users/David/Desktop/multiomics/")

# Libraries
# install.packages("devtools") 
# library(devtools)
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomeInfoDbData",dependencies=T)
# devtools::install_github("PMBio/MOFA", subdir="MOFAtools")
library(MOFAtools);library(MultiAssayExperiment)
# For the GCloud:
# library("googleCloudStorageR")

#system('gsutil cp gs://motrpac-portal-projects/ipop/iPOP_MultiOmics.RData .')

# Look at the ipop dataset
load("iPOP_MultiOmics.RData")
# ls()
# dim(rnaseq.df)
# rownames(rnaseq.df)
# colnames(rnaseq.df)
# dim(ck.df)
# dim(clinic.df)
# clinic.df[1:2,]
# table(clinic.df$SubjectID)
# colnames(clinic.df)

# Format the data and run MOFA

# # Look at the clinical data
# subjects = clinic.df$SubjectID
# times = clinic.df$CollectionDate
# times = sapply(times,function(x)strsplit(x,split=" ")[[1]][1])
# times = strptime(times,format = "%m/%d/%Y")
# table(is.na(times))
# status = clinic.df$CL4
# rownames(clinic.df) = clinic.df$SampleID

get_subj2first_ind<-function(x,col1="SampleID",col2="CL4"){
  subjects = x[,col1]
  status = x[,col2]
  subj2first_id=c()
  for(i in 1:length(subjects)){
    if(is.element(subjects[i],set=names(subj2first_id)) || status[i]!="Healthy"){next}
    subj2first_id[subjects[i]]=i
  }
  return(subj2first_id)
}
remove_cols<-function(x,regs = c("SampleID","SubjectID","CollectionDate","CL4")){
  inds = which(is.element(colnames(x),set=regs))
  return(x[,-inds])
}
get_cols<-function(x,regs = c("SampleID","SubjectID","CollectionDate","CL4")){
  inds = which(is.element(colnames(x),set=regs))
  return(x[,inds])
}

omics_data = list(
  "rnaseq" = rnaseq.df,
  "prot" = prot.df,
  "metab" = metb.df,
  "cytok" = ck.df,
  "clinic" = clinic.df
)
omics_data2meta = lapply(omics_data,get_cols)
omics_data = lapply(omics_data, remove_cols)
sapply(omics_data,dim)
omics_data = lapply(omics_data,function(x){x=as.matrix(x);mode(x)="numeric";x})
sapply(omics_data,class)
sapply(omics_data,mode)
sapply(omics_data,mean)
boxplot(t(omics_data$rnaseq[sample(1:500)[1:100],]))
omics_data$clinic[1:5,1:5]
clinic.df[1:5,1:5]

####### Preprocess the data #########
# Comments from the paper:
# RNA-seq
# - removed the genes with average read counts over all samples smaller than 1
# - samples with average read counts over all filtered genes smaller than 0.5 are filtered out
# - RData contains transcriptome in three formats:
#     raw read count data, variance stabilizing transformed (VST) data, and log transformed
# Metabolites and protein MS intensities were log10 transformed
# Transcripts were log2(n+1) transformed
# Only microbial taxa that were present (>0) or microbial genes having >1% abundance in > half
# of entire collection (> 400) were used
# Transcriptome, metabolome and cytokines were normalized based on size factor by log(X + 0.5)
# and converted into the log space for downstream analyses.
sapply(omics_data,mean)
# correct the metabolic range using log10 
omics_data$metab = log(omics_data$metab+0.5,base=10)
omics_data$cytok = log(omics_data$cytok+0.5,base=2)
sapply(omics_data,function(x)table(is.na(x)))

# med_h_data=list()
# for(nn in names(omics_data)){
#   x1 = omics_data[[nn]]
#   x2 = omics_data2meta[[nn]]
#   rownames(x1) = x2$SampleID
#   h_inds = x2$CL4=="Healthy"
#   h_x1 = x1[h_inds,]
#   h_x2 = x2[h_inds,]
#   med_h_x1 = c()
#   for(ss in unique(h_x2$SubjectID)){
#     currinds = h_x2$SubjectID==ss
#     if(sum(currinds)>1){
#       v = apply(h_x1[currinds,],2,median,na.rm=T)
#     }
#     else{
#       v = h_x1[currinds,]
#     }
#     med_h_x1 = rbind(med_h_x1,v)
#     rownames(med_h_x1)[nrow(med_h_x1)] = ss
#   }
#   med_h_data[[nn]] = med_h_x1
# }
# sapply(med_h_data, dim)
# sapply(med_h_data,function(x)table(is.na(x)))
# 
# shared_subjects = rownames(med_h_data[[1]])
# for(i in 1:length(med_h_data)){shared_subjects=intersect(shared_subjects,rownames(med_h_data[[i]]))}
# med_h_data_shared = lapply(med_h_data,function(x)x[shared_subjects,])
# sapply(med_h_data_shared, dim)
# sapply(med_h_data_shared,median)

# # Time series
# infections_data=list()
# for(nn in names(omics_data)){
#   x1 = omics_data[[nn]]
#   x2 = omics_data2meta[[nn]]
#   times = strptime(x2$CollectionDate,format = "%m/%d/%Y")
#   ord = order(x2$SubjectID,times)
#   x2 = x2[ord,];x1=x1[ord,];times=times[ord]
#   h_inds = x2$CL4=="Healthy"
#   inf_inds = x2$CL4 == "Infection"
#   inf_x1 = c()
#   for(i in 1:nrow(x1)){
#     if(!inf_inds[i]){next}
#     curr_subj = x2$SubjectID[i]
#     last_h_ind = h_inds & times < times[i] & x2$SubjectID==curr_subj
#     if(sum(h_inds & times < times[i] & x2$SubjectID==curr_subj,na.rm = T)==0){next}
#     last_h = max(which(last_h_ind))
#     v = x1[i,]-x1[last_h,]
#     curr_name = paste(curr_subj,x2$CollectionDate[i],sep=';')
#     inf_x1 = rbind(inf_x1,v)
#     rownames(inf_x1)[nrow(inf_x1)] = curr_name
#   }
#   infections_data[[nn]] = inf_x1
# }
# sapply(infections_data, dim)
# sapply(infections_data,function(x)table(is.na(x)))
# 
# shared_subjects = rownames(infections_data[[1]])
# for(i in 1:length(infections_data)){shared_subjects=intersect(shared_subjects,rownames(infections_data[[i]]))}
# infections_data_shared = lapply(infections_data,function(x)x[shared_subjects,])
# sapply(infections_data_shared, dim)
# sapply(infections_data_shared,median)
# corrs1 = cor(infections_data_shared$metab,infections_data_shared$prot,method = "spearman")
# library(gplots)
# heatmap.2(corrs1,trace="none")

# All data
for (nn in names(omics_data)){
  rownames(omics_data[[nn]]) = omics_data2meta[[nn]]$SampleID
}
sapply(omics_data, rownames)
shared_subjects = rownames(omics_data[[1]])
for(i in 1:length(omics_data)){shared_subjects=intersect(shared_subjects,rownames(omics_data[[i]]))}
omics_data_shared = lapply(omics_data,function(x)x[shared_subjects,])
sapply(omics_data_shared, dim)

########## MOFA ############
# http://htmlpreview.github.io/?https://github.com/bioFAM/MOFA/blob/master/MOFAtools/vignettes/MOFA_example_CLL.html
# clinic_ind = which(names(med_h_data_shared)=="clinic")
# mofa_multidata = lapply(med_h_data_shared,t)
# sapply(mofa_multidata,class)
# sapply(mofa_multidata,mode)
# hdf_name = "/ipop_healthy_mofa_res.hdf5"
# workspace_name = "ipop_analysis_workspace.RData"
# 
# clinic_ind = which(names(infections_data_shared)=="clinic")
# mofa_multidata = lapply(infections_data_shared,t)
# sapply(mofa_multidata,class)
# sapply(mofa_multidata,mode)
# hdf_name = "/ipop_infections_mofa_res.hdf5"
# workspace_name = "ipop_analysis_infections_workspace.RData"

clinic_ind = which(names(omics_data_shared)=="clinic")
mofa_multidata = lapply(omics_data_shared,t)
sapply(mofa_multidata,class)
sapply(mofa_multidata,mode)
hdf_name = "/ipop_allsamples_mofa_res.hdf5"
workspace_name = "ipop_analysis_allsamples_workspace.RData"

#sapply(mofa_multidata,colnames)
MAE = MultiAssayExperiment(experiments = mofa_multidata[-clinic_ind],
                           colData = t(mofa_multidata[[clinic_ind]]))
MOFAobject <- createMOFAobject(MAE)

# define I/O
DirOptions <- list(
  "dataDir" = tempdir(), # Folder to store the input matrices as .txt files, it can be a simple temporary folder
  "outFile" = paste(getwd(),hdf_name,sep='') # Output file of the model (use hdf5 extension)
)
DataOptions <- getDefaultDataOptions()
# Define Model Options
ModelOptions <- getDefaultModelOptions(MOFAobject)
# ModelOptions$likelihood <- c("gaussian","poisson")
# Define Training Options
TrainOptions <- getDefaultTrainOptions()
TrainOptions$maxiter <- 200
# Prepare MOFA object for training
MOFAobject <- prepareMOFA(MOFAobject, DirOptions, DataOptions, ModelOptions,TrainOptions)
# train a MOFA model
MOFAobject <- runMOFA(MOFAobject, DirOptions)
# save workspace
save.image(file=workspace_name)
system(paste('gsutil cp ',hdf_name, ' gs://motrpac-portal-projects/ipop/',sep=''))
system(paste('gsutil cp ',workspace_name,' gs://motrpac-portal-projects/ipop/',sep=''))

# Load and look at solutions
load('ipop_analysis_allsamples_workspace.RData')
load('ipop_analysis_infections_workspace.RData')

# Factors and features
# which factors are active in which view
r2 <- calculateVarianceExplained(MOFAobject)
# look at the top weight features in an active view
plotWeightsHeatmap(MOFAobject, "prot", factors=1:5, show_colnames=F)
plotWeights(MOFAobject, view = "rnaseq", factor = 1)
plotTopWeights(MOFAobject, "rnaseq", 1)
# look at the top features in the original dataset
plotDataHeatmap(MOFAobject, view="rnaseq", 
                factor=1, features=20, show_rownames=FALSE)

# Visualize samples along factors
rownames(mofa_multidata[[clinic_ind]])
plotFactorScatter(MOFAobject, 
      factors = 1:2, color_by = "LDL")
# overview of factor pairs
plotFactorScatters(MOFAobject, factors = 1:4, color_by = "LDL")
# use factors to cluster samples
h <- clusterSamples(MOFAobject, k=2, factors=1)
plotFactorScatters(MOFAobject, factors=1:3, color_by=h$cluster)

# CF analysis using CausalImpact
install.packages("CausalImpact")
library("CausalImpact")

