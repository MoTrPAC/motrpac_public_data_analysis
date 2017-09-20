simplify_tissue_info<-function(subject_col2tissue){
  subject_col2tissue_slim = subject_col2tissue
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

simplify_sex_info<-function(sexinfo){
  newv = as.character(sexinfo)
  names(newv) = names(sexinfo)
  newv[grepl(": f",newv)] = "female"
  newv[newv!="" & newv!="female"] = "male"
  newv[newv==""] = NA
  return(newv)
}

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
get_subject_info_from_gsms<-function(subj,gsms,metadata,colname,sample2subject){
  subj_gsms = gsms[sample2subject[gsms]==subj]
  currdata = as.character(metadata[subj_gsms,colname])
  return(unique(currdata)[1])
}

# Sept 2017 time windows: selected to have at least 3 cohorts in a tp
# The method below gives:
# blood muscle
# 1     17      5
# 4      4      8
# 24     6      7
simplify_time_acute<-function(tt){
  tt[tt<2] = 1
  tt[tt>=2 & tt<=5] = 4
  tt[tt>=20] = 24
  return(tt)
}
simplify_time_longterm_blood<-function(tt){
  tt[tt<=200] = 100
  return(tt)
}
simplify_time_longterm_muscle<-function(tt){
  tt[tt<=150] = 100
  tt[tt>150] = 200
  return(tt)
}
simplify_time_longterm<-function(tt,type="muscle"){
  if(type=="muscle"){
    return (simplify_time_longterm_muscle(tt))
  }
  if(type=="blood"){
    return(simplify_time_longterm_blood(tt))
  }
  tt[tt==70.1 | tt==70.2 | tt == 70] = 80
  tt[tt==84] = 80
  tt[tt>100 & tt<200] = 150
  return(tt)
}

# # Before Sept 2017
# simplify_time_acute<-function(tt){
#   tt[tt<=1] = 0.5
#   tt[tt>=2.5 & tt<=5] = 4
#   tt[tt==20] = 24
#   tt[tt==72 | tt==96] = 48
#   return(tt)
# }


####################### Working with GPL tables ####################
merge_two_id_lists<-function(l1,l2){
  lnew = lapply(1:length(l1),function(x,y,z)union(y[[x]],z[[x]]),y=l1,z=l2)
  return(lnew)
}
split_gpl_table_entry <-function(x){strsplit(x,split="( /// )|,|;",perl=T)}
map_gpl_rows_to_entrez<-function(x,converters){
  row2entrez=lapply(1:nrow(x),function(x)c())
  # look for entrez column
  is_entrez = which(grepl(colnames(x),pattern="entrez",ignore.case=T))
  if(length(is_entrez)>0){
    curr_entrez = as.character(x[,is_entrez[1]])
    row2entrez = split_gpl_table_entry(curr_entrez)
  }
  for(conv_fields in names(converters)){
    curr_cols = which(grepl(colnames(x),pattern=conv_fields,ignore.case=T))
    print(conv_fields)
    if(length(curr_cols)>0){
      for (j in curr_cols){
        currv = as.character(x[,j])
        currl = split_gpl_table_entry(currv)
        currl_vals = unique(unlist(currl))
        curr_conv = converters[[conv_fields]][currl_vals]
        curr_conv = curr_conv[!is.na(names(curr_conv))]
        currl = lapply(currl,function(x,y)unique(unlist(y[x])),y=curr_conv)
        row2entrez = merge_two_id_lists (row2entrez,currl)
      }
    }
  }
  return(row2entrez)
}
get_num_probe_genes <- function(x){
  if(is.null(x)){return(0)}
  return(length(na.omit((x))))
}
reverse_mapping_list<-function(l){
  newl = list()
  for(currname in names(l)){
    currvals = l[[currname]]
    for(v in currvals){
      v = as.character(v)
      newl[[v]] = union(newl[[v]],currname)
    }
  }
  return(newl)
}

########### Functions for transforming into genes ##############
# Convert a matrix x with probe ids as names into a new
# matrix of genes. 
# h is the list that maps each gene to its probe ids
# f is the function to apply, default is mean
probe2genes_conv<-function(x,h,f=mean, genes=NULL,...){
  if (is.null(genes)){genes = names(h)}
  y = sapply(genes,applyFunctionOnGroup,h=h,f=f,x=x,...)
  y = t(as.matrix(data.frame(y,check.names=F)))
  return (y)
}
# Apply a function f on values in the vector x (or a matrix)
# that are under the names of h[[g]]
# For example if h maps genes to probes then we can apply
# the mean value of the probes by setting f=mean
# We assume that all mapped probes are necessarily represented as 
# the names in x.
applyFunctionOnGroup<-function(g,h,x,f,intsct=F,...){
  curr_names = h[[g]]
  if (is.null(dim(x))){return (f(x[curr_names]))}
  if (length(curr_names)==1){return (sapply(x[curr_names,],f))}
  return (apply(x[curr_names,],2,f,...))
}

# clean a mapping list by the row names in m
clean_mapping<-function(m,h){
  newh = lapply(h,function(x,y)intersect(x,y),y=rownames(m))
  names(newh) = names(h)
  newh = newh[sapply(newh,length)>0]
  return(newh)
}

# A function that collects data on the NA cells in a matrix
collect_na_stats<-function(x){
  na_x = is.na(x)
  stats = list()
  stats[['percent NAs']] = sum(c(na_x))/length(na_x)
  stats[['num rows with NAs']] = sum(apply(na_x,1,any))
  stats[['num cols with NAs']] = sum(apply(na_x,2,any))
  stats[['row NA counts']] = rowSums(na_x)
  stats[['col NA counts']] = colSums(na_x)
  return(stats)
}
# A function that transforms a probe level data into a gene level matrix
# OPTIONAL: set out_path to save the new gene matrix to an RData file
# ... additional parameters for probe2genes_conv
transform_matrix_into_genes<-function(gsms_mat,gpl,gpl_mappings_entrez2probes,out_path=NULL,...){
  mode(gsms_mat) = 'numeric'
  gsms_mat_na_stats = collect_na_stats(gsms_mat)
  curr_mapping = gpl_mappings_entrez2probes[[gpl]]
  curr_mapping = clean_mapping(gsms_mat,curr_mapping)
  print('difference before and after cleaning the gpl mapping list:')
  print(paste(length(gpl_mappings_entrez2probes[[gpl]]),length(curr_mapping)))
  print('merging genes by averaging probes')
  entrez_mat = probe2genes_conv(gsms_mat,curr_mapping,na.rm=T)
  print('done')
  print(dim(entrez_mat))
  entrez_mat_na_stats = collect_na_stats(entrez_mat)
  return(list(gsms_mat_na_stats=gsms_mat_na_stats,entrez_mat=entrez_mat,entrez_mat_na_stats=entrez_mat_na_stats))
}

########################## Working with GEOquery data ##########################
get_gsm_data_size<-function(gsm){
  if(class(gsm)!="GSM"){return(0)}
  return(nrow(Table(gsm)))
}

########################## Extract data from our data contrainers ##########################
get_submatrix<-function(m,cols){
  cols = intersect(colnames(m),cols)
  if(length(cols)==0){return (NULL)}
  return (m[,cols])
}
# assumption: all gsms are of the same platform
get_data_matrix_from_matrix_lists<-function(m_list,gsms){
  l = lapply(m_list,get_submatrix,cols=gsms)
  l = l[!sapply(l,is.null)]
  l = l[sapply(l,length)>0]
  if(length(l)==0){return(NULL)}
  newm = c()
  for(mm in l){newm = cbind(newm,mm)}
  if(length(newm)<2 || is.null(newm)){return (NULL)}
  if(is.null(dim(newm))){newm = matrix(newm,ncol=1)}
  gsms = intersect(gsms,colnames(newm))
  if(length(gsms)==0){return(NULL)}
  return(newm[,gsms])
}

######################### Analysis and statistical tests ##########################

run_paired_test <- function(x,subjs,timev,times=NULL,func=t.test,...){
  if(is.null(times)){
    times = sort(unique(as.numeric(timev)))[1:2]
  }
  names(x) = subjs
  x1 = x[timev==times[1]];x2=x[timev==times[2]]
  x2 = x2[names(x1)]
  #print(all(names(x1)==names(x2)))
  return(func(x1,x2,paired=T,...))
}
# install.packages('lme4')
library(lme4)
run_lmer4_anova<-function(x,subjs,timev,sexv=NULL,agev=NULL,get_pvalue=T){
  has_sex = !is.null(sexv) && !any(sexv=="") && length(levels(sexv))>1
  has_age = !is.null(agev) && !any(agev=="") && (length(levels(agev))>1 || class(agev)=="numeric")
  if(has_sex && has_age){
    lm2 = lmer(x ~timev + sexv + agev + (1|subjs),REML=F)
    lm1 = lmer(x ~ 1 + sexv + agev + (1|subjs),REML=F)
  }
  else if(has_sex){
    lm2 = lmer(x ~ timev + sexv + (1|subjs),REML=F)
    lm1 = lmer(x ~ 1 + sexv + (1|subjs),REML=F)
  }
  else if(has_age){
    lm2 = lmer(x ~ timev + agev + (1|subjs),REML=F)
    lm1 = lmer(x ~ 1 + agev + (1|subjs),REML=F)
  }
  else{
    lm2 = lmer(x ~ timev + (1|subjs) ,REML=F)
    lm1 = lmer(x ~ 1 + (1|subjs),REML=F)
  }
  anova_obj = anova(lm1,lm2,refit=F)
  pval = anova_obj[2,8]
  if(get_pvalue){return(pval)}
  return(list(lm1=lm1,lm2=lm2,anova_obj=anova_obj,pval=pval))
}

# #####################################################################
# # QA of the functions above
# nsubj = 10
# nt = 2
# mean_diff = 4
# subj = factor(rep(1:nsubj,nt))
# timev = c(sapply(1:nt,function(x)rep(x,nsubj)))
# subj2sex = c(rep("F",nsubj/2),rep("M",nsubj/2)) # Balanced sexes
# sexv = factor(subj2sex[subj])
# # Effect due to time only
# x = rnorm(nsubj*nt)
# for(j in 2:nt){x[timev==j] = x[timev==j]+mean_diff}
# p1 = run_paired_ttest(x,subjs,timev)$p.value
# p2 = run_anova_on_time_and_subject(x,subjs,timev)$pval
# p3 = run_anova_on_time_age_sex_and_subject(x,subjs,timev)$pval
# print(abs(p1-p2))
# print(abs(p2-p3))
# # Effect due to sex only
# x = rnorm(nsubj*nt)
# for(j in unique(sexv)[-1]){x[sexv==j] = x[sexv==j]+mean_diff}
# p1 = run_paired_ttest(x,subjs,timev)$p.value
# p2 = run_anova_on_time_and_subject(x,subjs,timev)$pval
# p3 = run_anova_on_time_age_sex_and_subject(x,subjs,timev,sexv=sexv)$pval
# p4 = run_lmer4_anova(x,subjs,timev,sexv = sexv)
# summary(p3[[1]])
# print(abs(p1-p2))
# print(abs(p2-p3))

#####################################################################
# We currently take the first repeat when subjects have more than
# a single time series data
get_fold_changes_vs_baseline<-function(x,subjs,timev,baseline=NULL,func = function(a,b){a-b},metadata=NULL){
  if(is.null(baseline)){baseline = sort(unique(timev))[1]}
  print(baseline)
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

################ GO enrichment analysis ################
# GO enrichment
library(topGO)
run_topgo_enrichment_fisher<-function(genesOfInterest,geneUniverse,go_term_size=10,go_max_size=1000,go_dags=c("BP","MF"),...){
  l = list()
  if(class(genesOfInterest)=="character"){
    l[["set1"]] = genesOfInterest
  }
  else{
    l = genesOfInterest
  }
  res = c()
  for(nn in names(l)){
    genesOfInterest = l[[nn]]
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    for(type in go_dags){
      myGOdata <- new("topGOdata", description="My project", ontology=c("BP"), nodeSize=go_term_size,
                  allGenes=geneList, annot = annFUN.org, mapping="org.Hs.eg.db", ID = "entrez")
      allGOs = usedGO(myGOdata)
      resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")
      allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = 
                       "resultFisher", ranksOf = "classicFisher", topNodes = length(score(resultFisher)))
      setname = rep(nn,nrow(allRes))
      typev = rep(type,nrow(allRes))
      allRes = cbind(setname,typev,allRes)
      res = rbind(res,allRes)
    }
  }
  go_pvals = as.numeric(res[,ncol(allRes)])
  go_pvals[is.na(go_pvals)] = 1e-30
  go_qvals = p.adjust(go_pvals,method='fdr')
  res = cbind(res,go_qvals)
  return(res)
}
extract_top_go_results<-function(res,qval=0.1,maxsize=2000){
  res = res[,-which(colnames(res)=="typev")]
  res = res[res$Annotated<maxsize,]
  res = unique(res)
  res[res$GO.ID=="GO:0002532",]
  new_qvals = p.adjust(res$classicFisher,method='fdr')
  #plot(new_qvals,res$go_qvals);abline(0,1)
  res$go_qvals = new_qvals
  return(res[res$go_qvals<=qval,])
}
get_most_sig_enrichments_by_groups <- function(res,num=1){
  gs = unique(as.character(res[,1]))
  m = c()
  for(g in gs){
    res0 = res[res[,1]==g,]
    ps = as.numeric(res0$classicFisher)
    thr = sort(ps,decreasing = F)[min(num,length(ps))]
    m = rbind(m,res0[ps<=thr,])
  }
  return(m)
}

# GSEA
library(fgsea)
fgsea_wrapper <- function(pathways,scores,nperm=2000,run_nperm=1000,...){
  num_runs = nperm/run_nperm
  l = list()
  for(j in 1:num_runs){
    l[[j]] = fgsea(pathways,scores,nperm = run_nperm,...)
  }
  emp_pvals = sapply(l,function(x)x$pval)
  emp_pvals = emp_pvals*run_nperm
  min_to_add = min(emp_pvals)
  emp_pvals = emp_pvals-min_to_add
  new_pvals = rowSums(emp_pvals)+min_to_add
  new_pvals = new_pvals/nperm
  new_qvals = p.adjust(new_pvals,method='fdr')
  res = l[[1]]
  res[,"pval"] = new_pvals
  res[,"padj"] = new_qvals
  return(res)
}

## Simple functions for cohort-level analyses
# some simple functions to analyze a time point in a dataset's matrix
get_paired_ttest_yi_vi <-function(x,sample2time,t1,t2){
  x1 = x[sample2time==t1]
  x2 = x[sample2time==t2]
  d = x2-x1
  n = length(d)
  sdd = sd(d)/sqrt(length(d))
  return(c(yi=mean(d),vi=sdd^2))
}
# # test
# x1 = rnorm(10);x2=rnorm(10);x=c(x1,x2)
# s2t = c(rep(0,10),rep(1,10))
# pp =get_paired_ttest_yi_vi(x,s2t,0,1)
# pp
# tt = t.test(x1,x2,paired = T)
# tt$estimate/tt$statistic
# tt$estimate

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
  return(t.test(x1,x2,paired = T)$p.value)
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
