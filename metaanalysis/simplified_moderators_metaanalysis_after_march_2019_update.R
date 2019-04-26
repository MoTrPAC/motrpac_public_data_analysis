library(org.Hs.eg.db);library(metafor)
source('/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)

setwd('/Users/David/Desktop/MoTrPAC/project_release_feb_2018/data/')
# Prepare the datasets for the different analyses below
# Load the datasets and their metadata
load("human_ge_cohort_preprocessed_db_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
acute_sample2time = sample2time
acute_sample_meta = sample_metadata
load("human_ge_cohort_preprocessed_db_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata
longterm_sample2time = sample2time
acute_sample2time = sample2time
longterm_sample_meta = sample_metadata
load("human_ge_cohort_preprocessed_db_gene_tables.RData")

# Get some stats
# Total number of samples
tot_samples = length(union(
  unique(unlist(sapply(acute_metadata, function(x)x$gsms))),
  unique(unlist(sapply(longterm_metadata, function(x)x$gsms)))
))
# subjects, tissues, sex
complete_sample_table = c()
for(ac in names(acute_metadata)){
  curr_gsms = acute_metadata[[ac]]$gsms
  curr_subjects = acute_sample_meta$subject[curr_gsms]
  curr_subjects = paste(acute_metadata[[ac]]$gse,curr_subjects,sep=";")
  curr_sex = acute_sample_meta$sex[curr_gsms]
  curr_tissue = rep(acute_metadata[[ac]]$tissue,length(curr_gsms))
  curr_training = rep(acute_metadata[[ac]]$training,length(curr_gsms))
  curr_type = rep("acute",length(curr_gsms))
  m = cbind(curr_gsms,curr_subjects,curr_sex,curr_tissue,curr_training,curr_type)
  complete_sample_table = rbind(complete_sample_table,m)
}
for(lo in names(longterm_metadata)){
  curr_gsms = longterm_metadata[[lo]]$gsms
  curr_subjects = longterm_sample_meta$subject[curr_gsms]
  curr_subjects = paste(longterm_metadata[[lo]]$gse,curr_subjects,sep=";")
  curr_sex = longterm_sample_meta$sex[curr_gsms]
  curr_tissue = rep(longterm_metadata[[lo]]$tissue,length(curr_gsms))
  curr_training = rep(longterm_metadata[[lo]]$training,length(curr_gsms))
  curr_type = rep("longterm",length(curr_gsms))
  m = cbind(curr_gsms,curr_subjects,curr_sex,curr_tissue,curr_training,curr_type)
  complete_sample_table = rbind(complete_sample_table,m)
}
apply(complete_sample_table,2,table)

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

# Reshape data for replication analysis
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
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  pvalue_qqplot(simple_RE_pvals[[nn]],main=nn,pch=20,cex=0.5)
}

for(nn in names(simple_RE_pvals)){
  r = hist(simple_RE_I2s[[nn]],plot=F)
  per = sum(simple_RE_I2s[[nn]]>70,na.rm = T)/length(simple_RE_I2s[[nn]])
  per = format(per*100,digits=3)
  cols=rep("white",length(r$mids))
  cols[r$mids > 70] = "blue"
  hist(simple_RE_I2s[[nn]],main=paste(nn,"(",per,"%)",sep=""),xlab = "I^2(%)",col=cols)
}

# Filters 1 exclude genes with no sig p at 0.01
# before the update on March 2019 the threshold was 0.05 and 1
# Returns true if the gene has p<thr in at least num studies
rep_filter <-function(gdata,num=1,thr=0.01){
  return(sum(gdata$p<=thr,na.rm = T)>=num)
}
to_rem1 = sapply(datasets,function(x)!sapply(x,rep_filter,num=2,thr=0.05))

rm(acute_gene_tables);rm(longterm_gene_tables)
save.image(file="workspace_before_rep_analysis.RData")

# Create the input files for meta-regression and replication analysis
meta_reg_datasets = list()
for(nn in names(simple_RE_beta)){
  curr_genes = names(which(!to_rem1[[nn]]))
  meta_reg_datasets[[nn]] = datasets[[nn]][curr_genes]
}
sapply(meta_reg_datasets,length)

# After the data update of March 2019
# Not enough non-male acute muscle studies
# Not enough resistance training datasets for acute, blood
# For longterm blood we take training only as we assume that it will
# be the primary cause of differential abundance
meta_reg_to_mods = list(
  "acute,muscle" = c("time","training","avg_age"),
  "acute,blood" = c("time","prop_males","avg_age"),
  "longterm,muscle" = c("training","time","avg_age","prop_males"),
  "longterm,blood" = "training"
)
save(meta_reg_datasets,meta_reg_to_mods,
     untrained_datasets,file="meta_analysis_input.RData")

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# Meta-regression and model selection for selected genes
library(parallel);library(metafor)

# Load the results (run on sherlock) instead of running the code above
load("meta_analysis_results.RData")
load("workspace_before_rep_analysis.RData")
load("meta_analysis_input.RData")
source('/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')

I2_thr = 50
AIC_diff_thr = 5
ACUTE_beta_thr = 0.25
LONGTERM_beta_thr = 0.1
P_thr = 1e-04

# A few helper methods to work with the meta-analysis output
get_aicc_diff<-function(x){
  if(is.element("simple:base_model",set=names(x)) &&
     is.element("aic_c",set=names(x[["simple:base_model"]]))){
    val = x[[1]]$aic_c - x$`simple:base_model`$aic_c
    if(is.null(val) && is.na(val) || is.nan(val)){return(0)}
    if(is.infinite(val)){return(-1000)}
    return(unname(val))
  }
  if(is.element("base2:base_model",set=names(x)) &&
     is.element("aic_c",set=names(x[["base2:base_model"]]))){
    val = x[[1]]$aic_c - x$`base2:base_model`$aic_c
    if(is.null(val) && is.na(val) || is.nan(val)){return(0)}
    if(is.infinite(val)){return(-1000)}
    return(unname(val))
  }
  return(0)
}
get_simple_model_beta<-function(x){
  if(is.element("simple:base_model",set=names(x))){
    return(x[["simple:base_model"]]$coeffs[1,1])
  }
  if(is.element("base2:base_model",set=names(x))){
    return(x[["base2:base_model"]]$coeffs[1,1])
  }
  return(0)
}

# Algorithm for selecting genes from each meta-reg analysis
analysis2selected_genes = list()
analysis2selected_genes_stats = list()
all_pvals = c()
for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]]
  pvals = sapply(analysis1,function(x)x[[1]]$mod_p)
  all_pvals = c(all_pvals,pvals)
  i2s = simple_RE_I2s[[nn]][names(pvals)]
  # separate into two gene sets: those that passed the aic diff test vs. those that did not
  # define the set of filters
  # 1. AICc filter
  aic_diffs = sapply(analysis1,get_aicc_diff)
  genes_with_high_aic_diff = aic_diffs <= -AIC_diff_thr
  # 2. Beta filter
  if(grepl("acute",nn)){
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>ACUTE_beta_thr))
  }
  else{
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>LONGTERM_beta_thr))  
  }
  # 3. Is the top model simple base or is the AICc diff not large enough
  is_base_model = sapply(analysis1,function(x)names(x)[1] =="simple:base_model") | 
    !genes_with_high_aic_diff
  # 3.1 For base models make sure we get the correct beta value
  if(grepl("acute",nn)){
    model2beta[is_base_model] = sapply(analysis1[is_base_model],
         function(x)get_simple_model_beta(x)>ACUTE_beta_thr)
  }
  else{
    model2beta[is_base_model] = sapply(analysis1[is_base_model],
         function(x)get_simple_model_beta(x)>LONGTERM_beta_thr)  
  }
  # 4. Pval filter
  pval_filter = pvals <= P_thr
  # 5. I2 filter
  i2_filter = i2s <= I2_thr
  which(genes_with_high_aic_diff & i2_filter & pval_filter & model2beta)
  
  selected_aic_diff_genes = names(aic_diffs)[genes_with_high_aic_diff & model2beta & pval_filter]
  selected_base_model_genes = names(aic_diffs)[i2_filter & model2beta & pval_filter]
  selected_base_model_genes = setdiff(selected_base_model_genes,selected_aic_diff_genes)
  curr_selected_genes = union(selected_aic_diff_genes,selected_base_model_genes)
  curr_selected_genes_names = sapply(analysis1[selected_aic_diff_genes],function(x)names(x)[1])
  curr_selected_genes_names = sapply(curr_selected_genes_names,function(x)strsplit(x,split=":")[[1]][2])
  curr_selected_genes_names[selected_base_model_genes] = "base_model"
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
  colnames(m)= c("Entrez","Symbol","Group","Model pvalue","AICc diff","Coefficients")
  analysis2selected_genes_stats[[nn]] = m
}
sapply(analysis2selected_genes,length)

# Sanity checks and tests: our p-value threshold is lower than BY correction
max(all_pvals[p.adjust(all_pvals,method = "BY")<0.05])>1e-4

############################################################################
############################################################################
############################################################################
# Representation of the selected gene sets as bipartite graphs

bipartite_graphs = list()
for(nn in names(analysis2selected_genes_stats)){
  curr_edges = c()
  curr_m = analysis2selected_genes_stats[[nn]]
  for(i in 1:nrow(curr_m)){
    curr_gene = curr_m[i,"Symbol"]
    curr_entrez = curr_m[i,"Entrez"]
    curr_group = curr_m[i,"Group"]
    curr_coeffs = all_meta_analysis_res[[nn]][[curr_entrez]][[1]]$coeffs[,c("beta","pval")]
    if(is.null(dim(curr_coeffs))){curr_coeffs=matrix(curr_coeffs,nrow=1)}
    curr_coeffs[,2] = -log(curr_coeffs[,2],base=10)
    if(grepl("base_m",curr_group)){
      curr_edges = rbind(curr_edges,c(curr_gene,curr_entrez,"base_model",curr_coeffs))
      next
    }
    curr_gene_edges = c()
    for(j in 1:nrow(curr_coeffs)){
      curr_gene_edges = rbind(curr_gene_edges,
                              c(curr_gene,curr_entrez,
                                rownames(curr_coeffs)[j],curr_coeffs[j,c("beta","pval")]))
    }
    curr_edges = rbind(curr_edges,curr_gene_edges)
  }
  bipartite_graphs[[nn]] = curr_edges
}
sapply(bipartite_graphs,dim)

# Fix the feature name issue
for(nn in names(bipartite_graphs)){
  m = bipartite_graphs[[nn]]
  m[,3] = gsub(m[,3],pattern="time.L",replacement = "Time-linear")
  m[,3] = gsub(m[,3],pattern="time.Q",replacement = "Time-Q")
  m[,3] = gsub(m[,3],pattern="avg_age",replacement = "Age")
  m[,3] = gsub(m[,3],pattern="prop_males",replacement = "Sex")
  m[,3] = gsub(m[,3],pattern="trainingresistance",replacement = "Training-RE")
  m[,3] = gsub(m[,3],pattern="intrcpt",replacement = "b0")
  colnames(m) = c("Entrez","Symbol","Group","Effect","-log_P")
  bipartite_graphs[[nn]] = m
}

# Create interaction networks for each analysis, summarizing the number 
# of detected genes.
get_gene_set_by_feature_name<-function(m,fname,up=T){
  cnames = colnames(m)
  m = m[sapply(m[,3],grepl,name1),]
  if(is.null(dim(m))){
    m = matrix(m,nrow=1)
    colnames(m) = cnames
  }
  effects = as.numeric(m[,"Effect"])
  if(up){
    return(m[effects>0,2])
  }
  return(m[effects<0,2])
}
gene_overlaps = list()
gene_sets_per_cov = list()
for(nn in names(bipartite_graphs)){
  m = bipartite_graphs[[nn]]
  m = m[m[,3]!="b0",]
  # m = m[m[,3]!="base_model",]
  if(nrow(m)==0){next}
  covered_genes1 = unique(m[,1])
  # m = m[abs(as.numeric(m[,4]))>0.05,]
  covered_genes2 = unique(m[,1])
  setdiff(covered_genes1,covered_genes2)
  curr_features = unique(m[,3])
  if(length(curr_features)==1){next}
  curr_names = sort(c(paste(curr_features,",Up",sep=""),
                    paste(curr_features,",Down",sep="")))
  overlap_m = matrix(0,nrow=length(curr_names),ncol=length(curr_names),
                     dimnames = list(curr_names,curr_names))
  for(name1 in curr_names){
    gene_sets_per_cov[[paste(nn,name1,sep=",")]] = 
      get_gene_set_by_feature_name(m,name1,grepl(",Up",name1))
  }
  for(name1 in curr_names){
    set1 = gene_sets_per_cov[[paste(nn,name1,sep=",")]]
    for(name2 in curr_names){
      set2 = gene_sets_per_cov[[paste(nn,name2,sep=",")]]
      overlap_m[name1,name2] = length(intersect(set1,set2))
      overlap_m[name2,name1] = overlap_m[name1,name2]
    }
  }
  overlap_m = overlap_m[rowSums(overlap_m)>0,rowSums(overlap_m)>0]
  gene_overlaps[[nn]] = overlap_m
}
intersect(gene_sets_per_cov$`acute,muscle,Time-linear,Down`,
          gene_sets_per_cov$`acute,muscle,Time-Q,Up`)

library(corrplot)
corrplot(gene_overlaps[[1]],is.corr = F,method="number",type = "upper",cl.length = 5,
         cl.cex = 1.2,cl.ratio = 0.3,bg = "gray",mar=c(1, 0, 1, 0))
corrplot(gene_overlaps[[2]],is.corr = F,method="number",type = "upper",cl.length = 5,
         cl.cex = 1.2,cl.ratio = 0.3,bg = "gray",mar=c(1, 0, 1, 0))
corrplot(gene_overlaps[[3]],is.corr = F,method="number",type = "upper",cl.length = 5,
         cl.cex = 1.2,cl.ratio = 0.3,bg = "gray",mar=c(1, 0, 1, 0))

sort(sapply(gene_sets_per_cov,length))

bg = unique(c(unlist(sapply(all_meta_analysis_res,names))))
gs = gene_sets_per_cov
gs = gs[sapply(gs,length)>10]

go_res = run_topgo_enrichment_fisher(
  gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
go_res1 = go_res[go_res$Annotated < 1500,]
go_res1$classicFisher[is.na(as.numeric(go_res1$classicFisher))] = 1e-30
go_res1 = go_res[go_res1$Significant > 2,]
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
table(as.character(go_res_fdr$setname))
get_most_sig_enrichments_by_groups(go_res_fdr,num=2)[,1:4]
go_enrichments_by_cov_fdr = go_res_fdr

reactome_pathways_by_cov = run_reactome_enrichment_analysis(gs,universe=bg)
reactome_pathways_by_cov1 = reactome_pathways_by_cov[reactome_pathways_by_cov$Count>2,]
ps = reactome_pathways_by_cov1$pvalue
qs = p.adjust(ps,method="fdr")
reactome_pathways_by_cov1$qvalue = qs
reactome_pathways_by_cov_fdr = reactome_pathways_by_cov1[qs <= 0.1,]
table(reactome_pathways_by_cov_fdr[,1])
get_most_sig_enrichments_by_groups(reactome_pathways_by_cov_fdr,pcol="pvalue",num = 2)[,c(1,3)]
reactome_pathways_by_cov_fdr[,1] = as.character(reactome_pathways_by_cov_fdr[,1])

############################################################################
############################################################################
############################################################################
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
# par(mfrow=c(1,3),mar=c(10,4,2,2),cex.lab=1.2,cex.axis=1.2,cex=1,cex.main=0.9)
# for(nn in names(analysis2selected_genes)[1:3]){
#   x = simple_REs_untrained_beta[[nn]]
#   analyzed_genes = names(meta_reg_datasets[[nn]])
#   x = x[intersect(names(x),analyzed_genes)]
#   s = names(analysis2selected_genes[[nn]])
#   s_c = setdiff(names(x),s)
#   analysis1 = all_meta_analysis_res[[nn]][s]
#   y = sapply(analysis1,function(x)max(abs(x[[1]]$coeffs[,1])))
#   p = wilcox.test(x[s],y,paired=T)$p.value
#   p = format(p,digits=2)
#   l = list(
#     "untrained effects" = abs(x[s]),"exercise effects"=y,"other genes,\n untrained effects" = abs(x[s_c])
#   )
#   cols = c("blue","white","red")
#   boxplot(l,las=2,col=cols,horizontal=F,cex.axis=0.9,
#           ylab="Fold change",main = paste(nn," p=",p,sep=""),pch=20)
# }

# 2. alt, Look at effects in untrained - excluded untrained others
par(mar=c(8,4,2,2),cex.lab=1.2,cex.axis=1.2)
l = list();cols=c()
for(nn in names(analysis2selected_genes)[1:3]){
  x = simple_REs_untrained_beta[[nn]]
  analyzed_genes = names(meta_reg_datasets[[nn]])
  x = x[intersect(names(x),analyzed_genes)]
  s = names(analysis2selected_genes[[nn]])
  s_c = setdiff(names(x),s)
  analysis1 = all_meta_analysis_res[[nn]][s]
  y = sapply(analysis1,function(x)max(abs(x[[1]]$coeffs[,1])))
  p = wilcox.test(x[s],y,paired=T)$p.value
  p = format(p,digits=2)
  l[[paste(nn,"untrained effects",sep="\n")]] = abs(x[s])
  l[[paste(nn,"exercise effects",sep="\n")]] = y
  cols = c(cols,c("blue","red"))
}
dev.off()
par(mar=c(8,8,8,8))
names(l) = c("ac,mu,untr","ac,mu,exercise","ac,bl,untr","ac,bl,exercise",
             "lo,mu,untr","lo,mu,exercise")
boxplot(l,las=2,col=cols,horizontal=T,cex.axis=1,
        pch=20,ylim=c(0,1.6),xlab="Fold change")

# # 3. GO enrichments: not a must
# bg = unique(c(unlist(sapply(simple_REs,names))))
# gs = lapply(analysis2selected_genes,names)
# go_res = run_topgo_enrichment_fisher(
#   gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
# go_res1 = go_res[go_res$Annotated < 1500,]
# go_res1$classicFisher[is.na(as.numeric(go_res1$classicFisher))] = 1e-30
# go_res1 = go_res[go_res1$Significant > 3,]
# go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
# go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
# go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
# table(go_res_fdr$setname)
# gene_group_enrichments = go_res
# gene_group_enrichments_fdr = go_res_fdr
# get_most_sig_enrichments_by_groups(gene_group_enrichments_fdr,num=5)[,1:4]

############################################################################
############################################################################
############################################################################
# Some figures

dev.off()
par(mfrow=c(2,2))
for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]]
  model_names = sapply(analysis1,function(x)names(x)[1])
  table(model_names)
  ps1 = sapply(analysis1,function(x)x[[1]]$mod_p)
  i2s = sapply(analysis1,function(x)x[[1]]$I2)
  i2s_2 = simple_RE_I2s[[nn]][names(i2s)]
  aic_diffs = abs(sapply(analysis1,function(x)x[[1]]$aic_c - x[[2]]$aic_c))
  hist(aic_diffs,main=nn,xlab = "AICc difference")
}

# Examples of genes for the paper
curr_genes = analysis2selected_genes$`acute,muscle`
curr_genes[curr_genes=="base_model"]
#PGC1a
gene = "10891"
curr_genes[gene]
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,muscle`[[gene]]
gdata = meta_reg_datasets$`acute,muscle`[[gene]]
curr_m$slab.null = F
curr_times = rep("0-1h",nrow(gdata))
curr_times[gdata$time==2] = "2-5h"
curr_times[gdata$time==3] = ">20h"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
analysis1[[gene]][[1]]$mod_p
aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
dev.off()
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)
boxplot(yi~time,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",names=c("0-1h","2-5h",">20h"),ylab="Fold change")

# Other genes: LPL, CPT1B, SMAD3, ACTN3, VEGFA, FOXO1 and IL6R
genes = c("1375","4023","4088","89","7422","2308","3570")
par(mfrow=c(2,2))
for(gene in genes[5:7]){
  gene_name = entrez2symbol[[gene]]
  curr_m = simple_REs$`acute,muscle`[[gene]]
  gdata = meta_reg_datasets$`acute,muscle`[[gene]]
  curr_m$slab.null = F
  curr_times = rep("0-1h",nrow(gdata))
  curr_times[gdata$time==2] = "2-5h"
  curr_times[gdata$time==3] = ">20h"
  curr_m$slab = paste(gdata$training,curr_times,sep=",")
  analysis1 = all_meta_analysis_res$`acute,muscle`
  aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
  forest(curr_m,main=paste(gene_name),annotate = T)
}

# RPL24 in acute blood
gene = "3094"
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,blood`[[gene]]
gdata = meta_reg_datasets$`acute,blood`[[gene]]
curr_m$slab.null = F
curr_times = rep("0-1h",nrow(gdata))
curr_times[gdata$time==2] = "2-5h"
curr_times[gdata$time==3] = ">20h"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,blood`
analysis1[[gene]][[1]]
aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
forest(curr_m,main=paste(gene_name),annotate = T)


#RXRA
gene = "81848"
curr_genes[gene]
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`acute,muscle`[[gene]]
gdata = meta_reg_datasets$`acute,muscle`[[gene]]
curr_m$slab.null = F
curr_times = rep("0-1h",nrow(gdata))
curr_times[gdata$time==2] = "2-5h"
curr_times[gdata$time==3] = ">20h"
curr_m$slab = paste(gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
analysis1[[gene]][[1]]$mod_p
aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
dev.off()
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)
par(mar=c(6,6,6,6))
boxplot(yi~time,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time",names=c("0-1h","2-5h",">20h"),ylab="Fold change")

# Selected examples, before the March 2019 update
gene = "1282"
gene_name = entrez2symbol[[gene]]
curr_m = simple_REs$`longterm,muscle`[[gene]]
all_meta_analysis_res$`longterm,muscle`[[gene]]
curr_m$I2
gdata = meta_reg_datasets$`longterm,muscle`[[gene]]
curr_m$slab.null = F
curr_times = rep("",nrow(gdata))
curr_times[gdata$time==2] = ",> 150 days"
curr_m$slab = paste(gdata$training,curr_times,sep="")
analysis1 = all_meta_analysis_res$`acute,muscle`
forest(curr_m,main=paste(gene_name," all cohorts"),annotate = T)

############################################################################
############################################################################
############################################################################
# Interpretation of the results: defining subgroups by clustering mean patterns

# Some helper functions for reformatting the data
get_ts<-function(gdata){
  v = as.numeric(gdata$tstat)
  ns = paste(gdata$V1,gdata$training,gdata$time,
             format(gdata$avg_age,digits=1),format(gdata$prop_males,digits=1),sep=";")
  v = tapply(v,ns,mean)
  return(v)
}
get_t_matrix_from_list<-function(tstats){
  all_ns = unique(unlist(sapply(tstats,names)))
  all_genes = names(tstats)
  m = matrix(0,nrow=length(all_genes),ncol=length(all_ns),
             dimnames = list(all_genes,all_ns))
  for(i in 1:length(all_genes)){
    m[i,names(tstats[[i]])] = tstats[[i]]
  }
  return(m)
}

# Get effect matrices - mean responses - t statistics
# Reminder: the resulting matrices may have many zeroes because
# we basically merge all available t-statistics and not all genes
# are represented in all studies.
dataset_tstats = lapply(datasets,function(x)sapply(x,get_ts))
mean_effect_matrices = lapply(dataset_tstats,get_t_matrix_from_list)
sapply(mean_effect_matrices,dim)
sapply(mean_effect_matrices,function(x)table(x==0))
for(nn in names(mean_effect_matrices)){
  colnames(mean_effect_matrices[[nn]]) = 
    gsub("endurance","E",colnames(mean_effect_matrices[[nn]]),ignore.case = T)
  colnames(mean_effect_matrices[[nn]]) = 
    gsub("resistance","R",colnames(mean_effect_matrices[[nn]]),ignore.case = T)
  colnames(mean_effect_matrices[[nn]]) = 
    gsub("GE_","",colnames(mean_effect_matrices[[nn]]),ignore.case = T)
  colnames(mean_effect_matrices[[nn]]) = 
    gsub(" ","",colnames(mean_effect_matrices[[nn]]),ignore.case = T)
}
sapply(mean_effect_matrices,colnames)

# helper function for getting the clusters
library(cluster)
process_t_matrix<-function(data,thr1=1,thr2=5){
  data[data > -thr1 & data < thr1] = 0
  data[data > thr2] = thr2
  data[data < -thr2] = -thr2
  data = data[!apply(data==0,1,all),]
  return(data)
}
get_num_clusters_wss_kmeans<-function(data,k.max=10,wss_imp_thr=0.7){
  for(j in 1:ncol(data)){
    if(sd(data[,j])>0){
      tmp1 = mean(data[,j])
      tmp2 = sd(data[,j])
      data[,j] = (data[,j]-tmp1)/tmp2
    }
  }
  k.max = min(k.max,nrow(data)/2)
  k.max = min(k.max,length(unique(apply(data,1,paste,collapse="")))-1)
  wss <- sapply(1:k.max,
                function(k){kmeans(data, k, nstart=100,iter.max = 200)$tot.withinss})
  plot(1:k.max, wss,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  wss_imp_factors = wss[-1]/wss[-length(wss)]
  if(all(wss_imp_factors>wss_imp_thr)){return(1)}
  k = 2
  for(j in 2:length(wss_imp_factors)){
    if(wss_imp_factors[j] > wss_imp_thr){break}
    k = k+1
  }
  return(c(k=k))
}
# get_num_clusters_wss_kmeans(process_t_matrix(mmm))

gene_subgroups = list()
gene_t_patterns = list()
par(mfrow=c(2,2))
set.seed(123)
for(nn in names(analysis2selected_genes_stats)){
  curr_genes = analysis2selected_genes_stats[[nn]]
  curr_groups = curr_genes[,"Group"]
  for(gg in unique(curr_groups)){
    curr_m = as.matrix(curr_genes[gg==curr_genes[,"Group"],])
    if(ncol(curr_m)==1){curr_m=t(curr_m)}
    m = as.matrix(mean_effect_matrices[[nn]][curr_m[,"Entrez"],])
    if(ncol(m)==1){m=t(m)}
    rownames(m) = curr_m[,"Entrez"]
    
    if(nrow(m)<=5){
      m_kmeans = rep(1,nrow(m))
      names(m_kmeans) = curr_m[,"Entrez"]
      m = as.matrix(rowMeans(m),ncol=1)
      colnames(m)[1] = "mean_t"
    }
    else{
      # Merge based on the moderator's average not the t-test scores themselves
      if(gg != "base_model"){
        curr_mods = strsplit(gg,split=";")[[1]]
        curr_m_meta = sapply(colnames(m),function(x)strsplit(x,split=";")[[1]])
        rownames(curr_m_meta) = c("cohort","training","time","avg_age","prop_males")
        curr_m_meta = curr_m_meta[curr_mods,]
        if(!is.null(dim(curr_m_meta))){
          curr_m_meta = apply(curr_m_meta,2,paste,collapse=";")
        }
        m = t(apply(m,1,function(x,y)tapply(x,y,mean),y=curr_m_meta))
        if(nrow(m)>5){
          m_processed = process_t_matrix(m)
          # April 2019: change the clustering threshold
          # for large matrices default is fine
          # for smaller ones we need lower thresholds
          currk = get_num_clusters_wss_kmeans(m_processed,15)
          if(nrow(m)<50){
            currk = get_num_clusters_wss_kmeans(m_processed,15,wss_imp_thr = 0.6)
          }
        }
        else{currk=1}
        m_kmeans = kmeans(m_processed,centers = currk)
        m_kmeans = m_kmeans$cluster
      }
      else{
        m = as.matrix(rowMeans(m),ncol=1)
        colnames(m)[1] = "base,mean_t"
        currk=2
        m_kmeans = kmeans(m,centers = currk)
        m_kmeans = m_kmeans$cluster
      }
    }
    
    m = cbind(m,m_kmeans)
    colnames(m)[ncol(m)] = "kmeans_clusters"
    gene_t_patterns[[paste(nn,gg,sep=",")]] = m
    currk = length(unique(m_kmeans))
    print(paste(nn,gg,currk))
    print(table(m_kmeans))
    
    for(kk in unique(m_kmeans)){
      currname = paste(nn,gg,kk,sep=",")
      gene_subgroups[[currname]] = names(m_kmeans)[m_kmeans==kk]
    }
  }
}
sort(sapply(gene_subgroups,length))
base_model_ms = c()
for(gg in names(gene_t_patterns)[grepl("base_model",names(gene_t_patterns))]){
  print(gg)
  currm = gene_t_patterns[[gg]]
  currm = cbind(rep(gg,nrow(currm)),currm)
  base_model_ms = rbind(base_model_ms,currm)
}
gene_t_patterns[["base_models"]] = base_model_ms

# Enrichment analysis of the new groups
bg = unique(c(unlist(sapply(simple_REs,names))))
gs = gene_subgroups
gs = gs[sapply(gs,length)>5]
sapply(gs,length)
go_res = run_topgo_enrichment_fisher(
  gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
go_res1 = go_res[go_res$Annotated < 1500,]
go_res1$classicFisher[is.na(as.numeric(go_res1$classicFisher))] = 1e-30
go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
table(go_res_fdr$setname)
gene_subgroup_enrichments = go_res
gene_subgroup_enrichments_fdr = go_res_fdr

table(as.character(gene_subgroup_enrichments_fdr[,1]))
get_most_sig_enrichments_by_groups(gene_subgroup_enrichments_fdr,num=2)[,c(1,4,8)]
gene_subgroup_enrichments_fdr[grepl("ossi",gene_subgroup_enrichments_fdr$Term),]
get_most_sig_enrichments_by_groups(gene_group_enrichments_fdr,num=3)[,1:4]

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
# subgroups
reactome_pathways_subgroups = run_reactome_enrichment_analysis(gene_subgroups,universe=bg)
reactome_pathways_subgroups1 = reactome_pathways_subgroups[reactome_pathways_subgroups$Count>2,]
ps = reactome_pathways_subgroups1$pvalue
qs = p.adjust(ps,method="fdr")
reactome_pathways_subgroups1$qvalue = qs
reactome_pathways_subgroups_fdr = reactome_pathways_subgroups1[qs <= 0.1,]
table(reactome_pathways_subgroups_fdr[,1])
reactome_pathways_subgroups_fdr[,1] = as.character(reactome_pathways_subgroups_fdr[,1])

# We want to examine which clusters to analyze, sort by number of enrichments
enriched_clusters_go = sort(table(as.character(gene_subgroup_enrichments_fdr$setname)))
enriched_clusters_reactome = sort(table(as.character(reactome_pathways_subgroups_fdr[,1])))
all_enriched_clusters = union(names(enriched_clusters_go),names(enriched_clusters_reactome))

# Some stats about the clusters for the paper
length(gene_subgroups)
hist(sapply(gene_subgroups,length))
table(sapply(gene_subgroups,length) > 5)
large_clusters = names(which(sapply(gene_subgroups,length) > 5))
intersect(large_clusters,all_enriched_clusters)
length(intersect(large_clusters,all_enriched_clusters))
setdiff(all_enriched_clusters,large_clusters)

# Look at the top enrichments
get_most_sig_enrichments_by_groups(reactome_pathways_subgroups_fdr,pcol="pvalue",num = 2)[,c(1,3)]
get_most_sig_enrichments_by_groups(gene_subgroup_enrichments_fdr,num = 2)

# Age associated clusters
age_associated_clusters = all_enriched_clusters[grep("age",all_enriched_clusters)]
sex_associated_clusters = all_enriched_clusters[grep("male",all_enriched_clusters)]
training_associated_clusters = all_enriched_clusters[grep("train",all_enriched_clusters)]
sapply(gene_subgroups[age_associated_clusters],length)

# Heatmaps and line plots
# Analysis of the acute blood datasets
plot_with_err_bars<-function(xnames,avg,sdev,add=F,arrow_col="black",...){
  if(add){
    lines(avg,pch=19,...)
  }
  else{
    ylim = c(min(avg)-max(sdev),max(avg)+max(sdev))
    plot(avg,xaxt = "n",pch=19, type='l',ylim=ylim,...)
    axis(1, at=1:length(xnames), labels=xnames)
  }
  # hack: we draw arrows but with very special "arrowheads"
  arrows(1:length(xnames), avg-sdev, 1:length(xnames), avg+sdev,
         length=0.05, angle=90, code=3,col=arrow_col)
}
shorten_by_words<-function(x,num=5){
  if(is.na(x) || length(x)==0){return("")}
  arr =  strsplit(x,split="\\s|−")[[1]]
  num = min(num,length(arr))
  return(paste(arr[1:num],collapse=" "))
}
library(gplots)
hclust_func<-function(x){return(hclust(x,method = "ward.D2"))}
# Plot all clusters with enrichments
# pdf("supp_tables/all_heatmaps.pdf")
for(set_name in all_enriched_clusters){
  if(grepl("acute,muscle,time,",set_name)){next}
  arr = strsplit(set_name,split=",")[[1]]
  table_name = paste(arr[1:2],collapse=",")
  table_name2 = paste(arr[1:3],collapse=",")
  set_genes = gene_subgroups[[set_name]]
  mat = gene_t_patterns[[table_name2]][set_genes,]
  rownames(mat) = unlist(entrez2symbol[rownames(mat)])
  mat = mat[,-ncol(mat)]
  mat[mat>4]=4;mat[mat< -4]=-4
  if(is.null(dim(mat)) || nrow(mat)<2){next}
  
  if((grepl(",time;",set_name)||grepl(",time,",set_name))&& grepl("acute",set_name)){
    colnames(mat) = gsub("^1","0-1h;",colnames(mat))
    colnames(mat) = gsub("^2","2-5h;",colnames(mat))
    colnames(mat) = gsub("^3",">20h;",colnames(mat))
  }
  if((grepl(",time;",set_name)||grepl(",time,",set_name))&& !grepl("acute",set_name)){
    colnames(mat) = gsub("^1","",colnames(mat))
    colnames(mat) = gsub("^2","> 150 days;",colnames(mat))
  }
  curr_enrichments = "";curr_reactome="";curr_go=""
  colnames(mat) = gsub(";",", ",colnames(mat))
  mat = mat[,colnames(mat)!="NaN"]
  # get the enrichments of the set
  curr_go = gene_subgroup_enrichments_fdr[gene_subgroup_enrichments_fdr[,1]==set_name,]
  if(nrow(curr_go)>2){curr_go = curr_go[1:2,c(4,8)]}
  else{curr_go = curr_go[,c(4,8)]}
  curr_reactome = reactome_pathways_subgroups_fdr[reactome_pathways_subgroups_fdr[,1]==set_name,]
  if(nrow(curr_reactome)>2){curr_reactome = curr_reactome[1:2,c(3,6)]}
  else{curr_reactome = curr_reactome[,c(3,6)]}
  if(!is.null(dim(curr_go))){colnames(curr_go) = colnames(curr_reactome)}
  curr_enrichments = rbind(curr_go,curr_reactome)
  curr_enrichments = curr_enrichments[order(curr_enrichments[,2]),]
  num_enrichments = min(3,nrow(curr_enrichments))
  curr_enrichments[,2] = format(as.numeric(curr_enrichments[,2]),digits=2)
  curr_enrichments[,1] = sapply(as.character(curr_enrichments[,1]),shorten_by_words)
  curr_enrichments = curr_enrichments[1:num_enrichments,]
  curr_main = paste(set_name,"\n",sep="")
  print(dim(curr_enrichments))
  for(j in 1:nrow(curr_enrichments)){
    curr_e = paste(curr_enrichments[j,1]," (p=",curr_enrichments[j,2],")",sep="")
    curr_main = paste(curr_main,curr_e,"\n",sep="")
  }
  curr_main = gsub("\\.\\.\\.","",curr_main)
  cex_genes = 0.4
  if(nrow(mat)<40){cex_genes = 0.7}
  if(nrow(mat)<20){cex_genes = 1}
  if(nrow(mat)<10){cex_genes = 1}
  pdf(paste("supp_tables/",gsub(",|;","_",set_name),".pdf",sep=""))
  par(cex.main=0.7)
  heatmap.2(mat,trace = "none",scale = "none",Colv = T,col=bluered,
            cexRow = cex_genes,main=curr_main,
            Rowv = T,srtCol=45,hclustfun = hclust_func,density.info="none",
            key.title = NA,keysize = 1.1,key.xlab = "t-statistic",
            key.par = list("cex.axis"=1.1),margins = c(10,10))
  dev.off()
}
# dev.off()

# Training-specific responses (indep of time)

# # plot the t-stat over all cohorts
# set_genes = gene_subgroups[[set_name]]
# mat = mean_effect_matrices[[table_name]][set_genes,]
# rownames(mat) = unlist(entrez2symbol[rownames(mat)])
# mat = mat[,apply(mat,2,sd)>0]
# mat[mat>5]=5;mat[mat< -5]=-5
# heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 0.8)

# Specifically for longterm muscle base models:
makeRects <- function(m,lwd=2){
  coords = expand.grid(nrow(m):1, 1:ncol(m))[m,]
  xl=coords[,2]-0.49
  yb=coords[,1]-0.49
  xr=coords[,2]+0.49
  yt=coords[,1]+0.49
  rect(xl,yb,xr,yt,border="black",lwd=lwd)
}

set_name = "longterm,muscle,base_model,2"
table_name = "longterm,muscle"
set_genes = gene_subgroups[[set_name]]
mat = mean_effect_matrices[[table_name]][set_genes,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,apply(mat,2,sd)>0]
mat = mat[,apply(mat==0,2,sum)/nrow(mat) < 0.5]
mat = mat[apply(mat==0,1,sum)/ncol(mat) < 0.5,]
mat[mat>4]=4;mat[mat< -4]=-4
colnames(mat) = NULL;
# rownames(mat)=rep("",nrow(mat)) # comment this out for down-regulates
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 1.2,
          add.expr={makeRects(mat==0)},Rowv = F,margins = c(5,8))
# get the top enrichments of the set
curr_gos = gene_subgroup_enrichments_fdr[gene_subgroup_enrichments_fdr[,1]==set_name,]
curr_pathways = reactome_pathways_subgroups_fdr[reactome_pathways_subgroups_fdr[,1]==set_name,]








# Try another plot: select cluster names and plot all of their patterns
# along with selected enrichments










# The different patterns of acute muscle, time-dependent response
pref = "acute,muscle,time,\\d"
clusters = names(gene_subgroups)[grepl(pref,names(gene_subgroups))]
par(mfrow=c(2,2),cex.main=1.2)
cols = c("blue","black","red","green")
names(cols) = clusters
for(set_name in sort(clusters)){
  gene_set = gene_subgroups[[set_name]]
  top_enrichments1 = reactome_pathways_subgroups_fdr[
    reactome_pathways_subgroups_fdr[,1]==set_name,"Description"]
  top_enrichments2 = gene_subgroup_enrichments_fdr[
    gene_subgroup_enrichments_fdr[,1]==set_name,"Term"]
  
  # remove embryo enrichments
  top_enrichments1 = top_enrichments1[!grepl("embryo",top_enrichments1)]
  top_enrichments2 = top_enrichments2[!grepl("embryo",top_enrichments2)]
  
  enrichment1 = top_enrichments1[1]
  if(any(grepl("muscle",top_enrichments1))){
    enrichment1 = top_enrichments1[grepl("muscle",top_enrichments1)][1]
  }
  enrichment2 = top_enrichments2[1]
  if(any(grepl("muscle",top_enrichments2))){
    enrichment2 = top_enrichments2[grepl("muscle",top_enrichments2)][1]
  }
  if(length(top_enrichments1)> 1 && (is.na(enrichment2) || length(enrichment2)==0)){
    enrichment2 = top_enrichments1[2]
  }
  if(is.na(enrichment1)){enrichment1 = NULL}
  if(is.na(enrichment2)){enrichment2 = NULL}
  enrichment1 = shorten_by_words(tolower(enrichment1))
  enrichment2 = shorten_by_words(tolower(enrichment2))
  new_set_name = paste(set_name," (",length(gene_set)," genes)",sep="")
  curr_main = paste(new_set_name,enrichment1,enrichment2,sep="\n")
  curr_main = gsub("\n\n","\n",curr_main)
  print(curr_main)
  mat = gene_t_patterns$`acute,muscle,time`[gene_set,]
  print(paste(set_name,is.element("10891",set=rownames(mat))))
  rownames(mat) = unlist(entrez2symbol[rownames(mat)])
  mat = mat[,-ncol(mat)]
  plot_with_err_bars(c("0-1h","2-5h",">20h"),
                     colMeans(mat),apply(mat,2,sd),col=cols[set_name],lwd=3,
                     main=curr_main,ylab = "Mean t-statistic",xlab="Time",cex.main=1.1,
                     cex.lab=1.2,cex.axis=1.3,arrow_col = cols[set_name])
  abline(h = 0,lty=2,col="black")
}

# Write input to DREM
m = gene_t_patterns$`acute,muscle,time`[,1:3]
rownames(m) = entrez2symbol[rownames(m)]
write.table(m,file="drem/acute_muscle_time.txt",
            sep="\t",col.names = T,row.names = T,quote=F)
# Write input for cytoscape
m = gene_subgroups
m = m[grepl("acute,muscle,time,",names(m))]
mm = c()
for(cluster_name in names(m)){
  mm = rbind(mm,
             cbind(m[[cluster_name]],rep(cluster_name,length(m[[cluster_name]]))))
}
write.table(mm,file="supp_tables/acute_muscle_time_subgroups.txt",
            sep="\t",col.names = F,row.names = F,quote=F)

pref = "acute,blood,time,\\d"
clusters = names(gene_subgroups)[grepl(pref,names(gene_subgroups))]
par(mfrow=c(1,3))
cols = c("blue","red","green")
names(cols) = clusters
for(set_name in sort(clusters)[1:3]){
  gene_set = gene_subgroups[[set_name]]
  top_enrichments1 = reactome_pathways_subgroups_fdr[
    reactome_pathways_subgroups_fdr[,1]==set_name,"Description"]
  top_enrichments2 = gene_subgroup_enrichments_fdr[
    gene_subgroup_enrichments_fdr[,1]==set_name,"Term"]

  # remove embryo enrichments
  top_enrichments1 = top_enrichments1[!grepl("embryo",top_enrichments1)]
  top_enrichments2 = top_enrichments2[!grepl("embryo",top_enrichments2)]

  enrichment1 = top_enrichments1[1]
  if(any(grepl("muscle",top_enrichments1))){
    enrichment1 = top_enrichments1[grepl("muscle",top_enrichments1)][1]
  }
  enrichment2 = top_enrichments2[1]
  if(any(grepl("muscle",top_enrichments2))){
    enrichment2 = top_enrichments2[grepl("muscle",top_enrichments2)][1]
  }
  if(length(top_enrichments1)> 1 && (is.na(enrichment2) || length(enrichment2)==0)){
    enrichment2 = top_enrichments1[2]
  }
  if(is.na(enrichment1)){enrichment1 = NULL}
  if(is.na(enrichment2)){enrichment2 = NULL}
  enrichment1 = shorten_by_words(tolower(enrichment1))
  enrichment2 = shorten_by_words(tolower(enrichment2))
  new_set_name = paste(set_name," (",length(gene_set)," genes)",sep="")
  curr_main = paste(new_set_name,enrichment1,enrichment2,sep="\n")
  curr_main = gsub("\n\n","\n",curr_main)
  print(curr_main)
  mat = gene_t_patterns$`acute,blood,time`[gene_set,]
  rownames(mat) = unlist(entrez2symbol[rownames(mat)])
  mat = mat[,-ncol(mat)]
  plot_with_err_bars(c("0-1h","2-5h",">20h"),arrow_col = cols[set_name],
                     colMeans(mat),apply(mat,2,sd),lwd=2,xlab = "Time",
                     main=curr_main,ylab = "Mean t-statistic",col=cols[set_name])
}

# Write input to DREM
m = gene_t_patterns$`acute,blood,time`[,1:3]
rownames(m) = entrez2symbol[rownames(m)]
write.table(m,file="drem/acute_blood_time.txt",
            sep="\t",col.names = T,row.names = T,quote=F)
# Write input for cytoscape
m = gene_subgroups
m = m[grepl("acute,blood,time,",names(m))]
mm = c()
for(cluster_name in names(m)){
  mm = rbind(mm,
             cbind(m[[cluster_name]],rep(cluster_name,length(m[[cluster_name]]))))
}
write.table(mm,file="supp_tables/acute_blood_time_subgroups.txt",
            sep="\t",col.names = F,row.names = F,quote=F)

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
options(java.parameters = "-Xmx2g" )

# 1. Tables presenting the cohorts
longterm_meta_tmp = c()
for(gd in longterm_gene_tables){
  longterm_meta_tmp = rbind(longterm_meta_tmp,gd[,1:9])
  longterm_meta_tmp = unique(longterm_meta_tmp)
  if(length(setdiff(names(longterm_metadata),longterm_meta_tmp[,1]))==0){break}
  print(setdiff(names(longterm_metadata),longterm_meta_tmp[,1]))
}
rownames(longterm_meta_tmp) = longterm_meta_tmp[,1]
acute_meta_tmp = c()
for(gd in acute_gene_tables){
  acute_meta_tmp = rbind(acute_meta_tmp,gd[,1:9])
  acute_meta_tmp = unique(acute_meta_tmp)
  if(length(setdiff(names(acute_metadata),acute_meta_tmp[,1]))==0){break}
  print(setdiff(names(acute_metadata),acute_meta_tmp[,1]))
}
acute_meta_tmp = unique(acute_meta_tmp[,-2])
rownames(acute_meta_tmp) = acute_meta_tmp[,1]
cohorts_without_time_points = union(setdiff(names(longterm_metadata),longterm_meta_tmp[,1]),
                                    setdiff(names(acute_metadata),acute_meta_tmp[,1]))

metadata_row_for_supp_table<-function(cohort_name,cohort_metadata,
                                      cohort_info,sample2subject){
  samp2time = cohort_metadata[[nn]]$times
  curr_samps = samp2time[cohort_metadata[[cohort_name]]$gsms]
  curr_samps = curr_samps[!is.na(curr_samps)]
  Nsample = length(curr_samps)
  curr_samps[curr_samps==min(samp2time)]="Pre"
  Nsubject=length(unique(sample2subject[cohort_metadata[[cohort_name]]$gsms]))
  tps = paste(unique(curr_samps),collapse=",")
  res = c(cohort_name,"GSE"=cohort_metadata[[cohort_name]]$gse,
          "Tissue"=cohort_metadata[[cohort_name]]$tissue,
          "Training"=cohort_metadata[[cohort_name]]$training,
          "Nsample"=Nsample,"Nsubject"=Nsubject,
          "Time_points"=tps,"Avg_age"=cohort_info[cohort_name,"avg_age"],
          "Prop_males"=cohort_info[cohort_name,"prop_males"],
          "Additional Info"=cohort_metadata[[cohort_name]]$additional_info)
  return(res)
}
m1 = c()
for(nn in names(acute_metadata)){
  curr_info = metadata_row_for_supp_table(
    nn,acute_metadata,acute_meta_tmp,acute_sample_meta$subject)
  m1 = rbind(m1,curr_info)
}
for(nn in names(longterm_datasets)){
  curr_info = metadata_row_for_supp_table(
    nn,longterm_metadata,longterm_meta_tmp,longterm_sample_meta$subject)
  m1 = rbind(m1,curr_info)
}
supp_table_1_all_cohorts = m1
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
length(unique(supp_table_1_all_cohorts[
  supp_table_1_all_cohorts[,ncol(supp_table_1_all_cohorts)]!="","GSE"]))

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

sheet_counter=3
sapply(bipartite_graphs,function(x)unique(x[,3]))
for(nn in names(bipartite_graphs)){
  m = bipartite_graphs[[nn]]
  m[,3] = gsub(m[,3],pattern="time.L",replacement = "Time-linear")
  m[,3] = gsub(m[,3],pattern="time.Q",replacement = "Time-Q")
  m[,3] = gsub(m[,3],pattern="avg_age",replacement = "Age")
  m[,3] = gsub(m[,3],pattern="prop_males",replacement = "Sex")
  m[,3] = gsub(m[,3],pattern="trainingresistance",replacement = "Training-RE")
  m[,3] = gsub(m[,3],pattern="intrcpt",replacement = "b0")
  colnames(m) = c("Entrez","Symbol","Group","Effect","-log_P")
  write.xlsx(m,file=supp_file,
               sheetName = paste("STable",sheet_counter,"_",nn,sep=""),
               row.names = F,append = T,col.names = T)
  write.table(m,sep="\t",
             file = paste("supp_tables/STable",sheet_counter,"_",nn,".txt",sep=""),
             row.names = F,col.names = T,quote=F)
  sheet_counter = sheet_counter+1
}

# GO enrichments of subgroups
supp_table_enrichments = gene_subgroup_enrichments_fdr[,c(1:4,9:10)]
colnames(supp_table_enrichments)[6] = "q-value"
colnames(supp_table_enrichments)[5] = "Genes"
colnames(supp_table_enrichments)[1] = "Discovered in"
write.xlsx(supp_table_enrichments,file=supp_file,
           sheetName = paste("STable",sheet_counter,"_GO_enrichments",sep=""),
           row.names = F,append=T)

# Same for Reactome
sheet_counter = sheet_counter+1
supp_table_enrichments = reactome_pathways_subgroups_fdr[,c(1,3,8,9)]
colnames(supp_table_enrichments)[3] = "q-value"
colnames(supp_table_enrichments)[1] = "Discovered in"
write.xlsx(supp_table_enrichments,file=supp_file,
           sheetName = paste("STable",sheet_counter,"_Reactome_enrichments",sep=""),
           row.names = F,append=T)

# save the workspace
save.image(file="meta_analysis_interpretation_results.RData")

# # Other QCs
# # Look at the time-associated genes and their L/Q p-values
# set_name="acute,muscle,time"
# gene_set = rownames(gene_t_patterns[[set_name]])
# coeffes = lapply(all_meta_analysis_res$`acute,muscle`[gene_set],
#                  function(x)x[[1]]$coeffs)
# pvals = sapply(coeffes,function(x)x[,"pval"])
# hist(pvals[1,],breaks=100)
# hist(pvals[2,],breaks=100)
# hist(pvals[3,],breaks=100)
# 
# set_name="acute,blood,time"
# gene_set = rownames(gene_t_patterns[[set_name]])
# length(gene_set)
# coeffes = lapply(all_meta_analysis_res$`acute,blood`[gene_set],
#                  function(x)x[[1]]$coeffs)
# pvals = sapply(coeffes,function(x)x[,"pval"])
# hist(pvals[1,],breaks=100)
# hist(pvals[2,],breaks=100)
# hist(pvals[3,],breaks=100)

