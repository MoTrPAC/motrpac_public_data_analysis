#### Compare databases of expression profiles ###
# # Load the newly created database
# load("/Users/David/Desktop/MoTrPAC/project_release_feb_2018/data/human_ge_profiles_db.RData")
# gse_matrices_new = gse_matrices
# gsm_objs_new = gsm_objs
# gpl_tables_new = gpl_tables
# CEL_frma_profiles_new = CEL_frma_profiles
# CEL_rma_profiles_new = CEL_rma_profiles
# # Load the old db
# load("/Users/David/Desktop/MoTrPAC/PA_database/PA_database_profiles.RData")
# # Differences in RMA matrices
# x1 = CEL_rma_profiles
# x2 = CEL_rma_profiles_new
# assertthat::are_equal(length(x1),length(x2))
# setdiff(names(x1),names(x2))
# all(is.element(names(x1),names(x2)))
# for(nn in intersect(names(x1),names(x2))){
#   xx1 = x1[[nn]]
#   xx2 = x2[[nn]]
#   if (!all(dim(xx1)==dim(xx2))){
#     print(dim(xx1));print(dim(xx2))
#     break
#   }
#   try({print(all(xx1==xx2,na.rm = T))})
# }
# # Differences in the GPL mappings: none observed
# x1 = gpl_tables
# x2 = gpl_tables_new
# assertthat::are_equal(length(x1),length(x2))
# all(is.element(names(x1),names(x2)))
# for(nn in intersect(names(x1),names(x2))){
#   xx1 = x1[[nn]]
#   xx2 = x2[[nn]]
#   if(length(xx1)==0 && length(xx2)==0){next}
#   try({print(all(xx1==xx2,na.rm = T))})
# }
# # Differences in the GSE sets
# x1 = gse_matrices
# x2 = gse_matrices_new
# assertthat::are_equal(length(x1),length(x2))
# setdiff(names(x1),names(x2))
# all(is.element(names(x1),names(x2)))
# for(nn in intersect(names(x1),names(x2))){
#   xx1 = x1[[nn]]
#   xx2 = x2[[nn]]
#   if (dim(xx1)!=dim(xx2)){
#     print(dim(xx1));print(dim(xx1))
#     next
#   }
#   try({print(all(xx1==xx2,na.rm = T))})
# }
# # Differences in the GSM set:
# # 72 GSMs from the GSE44818 dataset were removed during the manual inspection process
# x1 = gsm_objs
# x2 = gsm_objs_new
# assertthat::are_equal(length(x1),length(x2))
# setdiff(names(x1),names(x2))
# all(is.element(names(x1),names(x2)))
# for(nn in intersect(names(x1),names(x2))){
#   xx1 = x1[[nn]]
#   xx2 = x2[[nn]]
#   try({print(all(xx1==xx2,na.rm = T))})
# }

############ preprocessed matrices ###########
load("/Users/David/Desktop/MoTrPAC/project_release_feb_2018/data/human_ge_cohort_preprocessed_db_acute.RData")
cohort_data_new = cohort_data
cohort_metadata_new = cohort_metadata
sample2time_new = sample2time
sample2sex_new = sample2sex
load("/Users/David/Desktop/MoTrPAC/PA_database/PADB_univariate_results_and_preprocessed_data_acute.RData")
x1 = cohort_data_new
x2 = cohort_data
assertthat::are_equal(length(x1),length(x2))
setdiff(names(x1),names(x2))
for(nn in names(x1)){
  print (paste("**** analyzing: ",nn))
  xx1 = x1[[nn]];xx2=x2[[nn]]
  #print(all(xx1$gene_data==xx2$gene_data))
  print(max(abs(xx1$probe_data-xx2$probe_data)))
  #if(max(abs(xx1$gene_data-xx2$gene_data))>1e-10){break}
  print(table(is.infinite(xx1$probe_data)))
  print("####")
  print(table(is.infinite(xx2$probe_data)))
}
for(nn in names(which(sapply(cohort_data,length)>2))){
  print (paste("**** analyzing: ",nn))
  xx1 = x1[[nn]]$time2ttest_stats;xx2=x2[[nn]]$time2ttest_stats
  print(quantile(abs(xx1[[1]][,3]-xx2[[1]][,3])))
}
x1 = sample2time_new
x2 = sample2time
assertthat::are_equal(length(x1),length(x2))
all(x1==x2)
# Longterm data comparison
load("/Users/David/Desktop/MoTrPAC/project_release_feb_2018/data/human_ge_cohort_preprocessed_db_longterm.RData")
cohort_data_new = cohort_data
cohort_metadata_new = cohort_metadata
sample2time_new = sample2time
sample2sex_new = sample2sex
load("/Users/David/Desktop/MoTrPAC/PA_database/PADB_univariate_results_and_preprocessed_data_longterm.RData")
x1 = cohort_data_new
x2 = cohort_data
assertthat::are_equal(length(x1),length(x2))
setdiff(names(x1),names(x2))
for(nn in names(which(sapply(cohort_data,length)>2))){
  print (paste("**** analyzing: ",nn))
  xx1 = x1[[nn]]$time2ttest_stats;xx2=x2[[nn]]$time2ttest_stats
  print(quantile(abs(xx1[[1]][,3]-xx2[[1]][,3]),na.rm=T))
}
for(nn in names(x1)){
  xx1 = x1[[nn]];xx2=x2[[nn]]
  print(all(xx1$gene_data==xx2$gene_data))
  print(all(xx1$probe_data==xx2$probe_data))
  if(!all(xx1$gene_data==xx2$gene_data)){
    meta1 = cohort_metadata[[nn]]
    meta2 = cohort_metadata_new[[nn]]
    print(all(meta1$gsms==meta2$gsms))
    cohort_metadata[[nn]]
    print(max(abs(xx1$probe_data-xx2$probe_data)))
    break
  }
}
x1 = sample2time_new
x2 = sample2time
assertthat::are_equal(length(x1),length(x2))
all(x1==x2)


