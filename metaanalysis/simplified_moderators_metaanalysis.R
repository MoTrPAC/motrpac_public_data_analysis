library(org.Hs.eg.db);library(metafor)
source('~/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)

# setwd('~/Desktop/MoTrPAC/project_release_feb_2018/data/')
out_dir = "~/Desktop/MoTrPAC/project_release_feb_2018/revision_feb_2020/"
out_dir_figs = paste(out_dir,"figures/",sep="")
try(dir.create(out_dir_figs))
out_dir_rdata = paste(out_dir,"rdata/",sep="")
try(dir.create(out_dir_rdata))
try(dir.create(paste(out_dir,"supp_tables/",sep="")))
setwd(out_dir)

# Meta-regression and model selection for selected genes
library(parallel);library(metafor)

# Load the results (run on sherlock) instead of running the code above
load("meta_analysis_results.RData")
load("workspace_before_rep_analysis.RData")
load("meta_analysis_input.RData")
source('~/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')

I2_thr = 50
AIC_diff_thr = 5
ACUTE_beta_thr = 0.1
LONGTERM_beta_thr = 0.1
P_thr = 1e-03

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
  i2s[is.na(i2s)] = 100
  # separate into two gene sets: those that passed the aic diff test vs. those that did not
  # define the set of filters
  # 1. AICc filter
  aic_diffs = sapply(analysis1,get_aicc_diff)
  genes_with_high_aic_diff = aic_diffs <= -AIC_diff_thr
  # 2. Beta filter
  curr_effect_thr = ACUTE_beta_thr
  if(grepl("acute",nn)){
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>ACUTE_beta_thr))
  }
  else{
    model2beta = sapply(analysis1,function(x)any(abs(x[[1]]$coeffs[,1])>LONGTERM_beta_thr))  
    curr_effect_thr = LONGTERM_beta_thr
  }
  # 3. Is the top model simple base or is the AICc diff not large enough
  is_base_model = sapply(analysis1,function(x)names(x)[1] =="simple:base_model") | 
    !genes_with_high_aic_diff
  # 3.1 For base models make sure we get the correct beta value
  model2beta[is_base_model] = sapply(analysis1[is_base_model],
           function(x)get_simple_model_beta(x)>curr_effect_thr) 
  # 4. Pval filter
  pval_filter = pvals <= P_thr
  # 5. I2 filter
  i2_filter = i2s <= I2_thr
  
  # Selection process: use the model selection only for the muscle datasets
  if(grepl("blood",nn)){
    curr_selected_genes = i2_filter & pval_filter
    curr_selected_genes = curr_selected_genes & 
      abs(simple_RE_beta[[nn]][names(curr_selected_genes)]) > curr_effect_thr
    curr_selected_genes = names(curr_selected_genes)[curr_selected_genes]
    curr_selected_genes_names = rep("base_model",length(curr_selected_genes))
    names(curr_selected_genes_names) = curr_selected_genes
    coeffs = lapply(simple_REs[[nn]][curr_selected_genes],
                       function(x){y=cbind(x$beta,x$pval);colnames(y)=c("beta","pval");y})
    coeffs_v = sapply(coeffs,get_coeffs_str)
  }
  else{
    
    selected_aic_diff_genes = names(aic_diffs)[genes_with_high_aic_diff & model2beta & pval_filter]
    selected_base_model_genes = names(aic_diffs)[i2_filter & model2beta & pval_filter]
    selected_base_model_genes = setdiff(selected_base_model_genes,selected_aic_diff_genes)
    selected_base_model_genes = selected_base_model_genes[!is.na(selected_base_model_genes)]
    
    curr_selected_genes = union(selected_aic_diff_genes,selected_base_model_genes)
    curr_selected_genes_names = sapply(analysis1[selected_aic_diff_genes],function(x)names(x)[1])
    curr_selected_genes_names = sapply(curr_selected_genes_names,function(x)strsplit(x,split=":")[[1]][2])
    curr_selected_genes_names[selected_base_model_genes] = "base_model"
    coeffs = lapply(analysis1[names(curr_selected_genes_names)],function(x)x[[1]]$coeffs)[selected_aic_diff_genes]
    coeffs[selected_base_model_genes] = lapply(simple_REs[[nn]][selected_base_model_genes],
                                               function(x){y=cbind(x$beta,x$pval);colnames(y)=c("beta","pval");y})
    coeffs = coeffs[curr_selected_genes]
    coeffs_v = sapply(coeffs,get_coeffs_str)
  }
  
  analysis2selected_genes[[nn]] = curr_selected_genes_names
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

# Sanity checks and tests: our p-value threshold is similar to the FDR adjusted one
fdr0.01_thr = max(all_pvals[p.adjust(all_pvals,method = "fdr")<0.01])

# Overlaps
names(analysis2selected_genes_stats)
intersect(analysis2selected_genes_stats[[1]][,2],analysis2selected_genes_stats[[2]][,2])
intersect(analysis2selected_genes_stats[[1]][,2],analysis2selected_genes_stats[[3]][,2])
intersect(analysis2selected_genes_stats[[3]][,2],analysis2selected_genes_stats[[2]][,2])
intersect(analysis2selected_genes_stats[[3]][,2],analysis2selected_genes_stats[[4]][,2])

# Add GSEA enrichments: use the inferred beta values

gsea_input_scores = list()
for(nn in names(simple_RE_beta)){
  currbetas = simple_RE_beta[[nn]]
  curri2 = simple_RE_I2s[[nn]]
  curri2 = curri2/100
  currps = simple_RE_pvals[[nn]]
  r_currbeta = sign(currbetas) * rank(abs(currbetas))
  r_curri2 = sign(currbetas) * rank(1-curri2)
  gsea_input_scores[[nn]] = cbind(currbetas,curri2,r_currbeta,r_curri2)
}

# # redo the lonterm muscle analysis without the overlap with the acute analysis
# remove_by_gse<-function(gtable,gses){
#   gtable = gtable[! (gtable$gse %in% gses),]
#   return(gtable)
# }
# longterm_reduced_datasets = lapply(datasets[["longterm,muscle"]],remove_by_gse,gses = c(
#   "GSE59088","GSE27285","GSE28998","GSE28392","GSE41769","GSE45426","GSE106865","GSE107934",
#   "GSE19062","GSE43219","GSE43856","GSE87749"
# ))
# simple_re_longterm_reduced = lapply(longterm_reduced_datasets,run_simple_re)
simple_re_longterm_reduced_beta = 
  sapply(simple_re_longterm_reduced,try_get_field,fname="beta")
simple_re_longterm_reduced_i2 = 
  sapply(simple_re_longterm_reduced,try_get_field,fname="I2")/100
gsea_input_scores[["longterm_overlap_removed"]] = cbind(
  simple_re_longterm_reduced_beta,
  simple_re_longterm_reduced_i2,
  sign(simple_re_longterm_reduced_beta) * rank(simple_re_longterm_reduced_beta),
  sign(simple_re_longterm_reduced_beta) * rank(simple_re_longterm_reduced_i2)
)

library(fgsea)
gsea_reactome_results = list()
for(nn in names(gsea_input_scores)){
  # beta GSEA
  currscores = gsea_input_scores[[nn]][,3]
  currscores = na.omit(currscores)
  currscores = sample(currscores)
  pathways = reactomePathways(names(currscores))
  pathways = pathways[sapply(pathways,length)>10]
  pathways = pathways[sapply(pathways,length)<200]
  all_p_genes = unique(unlist(pathways))
  fgsea_res = fgsea(pathways,currscores[all_p_genes],nperm=50000)
  fgsea_res = fgsea_res[order(fgsea_res$padj),]
  gsea_reactome_results[[paste("beta",nn,sep="_")]] = fgsea_res
  # I2 GSEA
  currscores = gsea_input_scores[[nn]][,4]
  currscores = na.omit(currscores)
  currscores = sample(currscores)
  pathways = reactomePathways(names(currscores))
  pathways = pathways[sapply(pathways,length)>10]
  pathways = pathways[sapply(pathways,length)<200]
  all_p_genes = unique(unlist(pathways))
  fgsea_res = fgsea(pathways,currscores[all_p_genes],nperm=50000)
  fgsea_res = fgsea_res[order(fgsea_res$padj),]
  gsea_reactome_results[[paste("I2",nn,sep="_")]] = fgsea_res
}


all_gsea_ps = unlist(sapply(gsea_reactome_results,function(x)x$pval))
hist(all_gsea_ps)
fdr_BY_thr = max(all_gsea_ps[p.adjust(all_gsea_ps,method="BY") < 0.1],na.rm=T)
gsea_reactome_results_fdr = lapply(gsea_reactome_results,
                               function(x)x[x$pval<fdr_BY_thr,])
sapply(gsea_reactome_results_fdr,nrow)
lapply(gsea_reactome_results_fdr,function(x)x[1:3,1:5])

save(gsea_reactome_results,gsea_reactome_results_fdr,gsea_input_scores,
     file = paste(out_dir_rdata,"gsea_reactome_results.RData",sep=""))

# merge all gsea results into one table
all_fdr_corrected_gsea_results = c()
for(nn in names(gsea_reactome_results_fdr)){
  m = gsea_reactome_results_fdr[[nn]]
  m = as.data.frame(m)
  m$analysis = nn
  all_fdr_corrected_gsea_results = rbind(all_fdr_corrected_gsea_results,m)
}
all_fdr_corrected_gsea_results$leadingGenes = sapply(all_fdr_corrected_gsea_results$leadingEdge,
  function(x,y)y[x],y=entrez2symbol)
all_fdr_corrected_gsea_results$leadingEdge = sapply(all_fdr_corrected_gsea_results$leadingEdge,
                                                    paste,collapse=";")
all_fdr_corrected_gsea_results$leadingGenes = sapply(all_fdr_corrected_gsea_results$leadingGenes,
                                                    paste,collapse=";")
all_fdr_corrected_gsea_results = as.matrix(all_fdr_corrected_gsea_results)
all_fdr_corrected_gsea_results = all_fdr_corrected_gsea_results[,
      c("analysis","pathway","padj","NES","pval","ES","size","leadingGenes","leadingEdge")]
write.table(all_fdr_corrected_gsea_results,file=paste(supp_path,"STable9.txt",sep="")
            ,row.names = F,quote=F,sep="\t")

# GSEA plots
num_pathways = 6
library(gridExtra)
library(grid)
for(nn in names(gsea_input_scores)){
  currscores = gsea_input_scores[[nn]][,3]
  currscores = na.omit(currscores)
  currscores = sample(currscores)
  curr_res = gsea_reactome_results[[paste("beta",nn,sep="_")]]
  inds1 = which(curr_res[["NES"]]>0)
  inds2 = which(curr_res[["NES"]]<0)
  curr_selected_pathways = as.data.frame(curr_res[,1])[
    c(inds1[1:(num_pathways/2)],inds2[1:(num_pathways/2)]),1]
  pathways = reactomePathways(names(currscores))
  all_p_genes = unique(unlist(pathways))
  curr_f_name = paste0(out_dir_figs,"gsea_top6_",gsub(",","_",nn),".pdf")
  pdf(curr_f_name)
  currpathways = pathways[curr_selected_pathways]
  curr_res = as.data.frame(curr_res)
  pl = plotGseaTable(currpathways,
                currscores[all_p_genes],curr_res,gseaParam=0.5,
                colwidths = c(5, 3, 0.8, 0, 1.2),render=T)
  dev.off()
  write.table(t(t(names(currpathways))),row.names = F,
              col.names = F,quote=F)
}

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
table(bipartite_graphs$`acute,blood`[,3])

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
  if(nrow(m)==0){next}
  curr_features = unique(m[,3])
  # if(length(curr_features)==1){next}
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

# dev.off()
pdf(paste(out_dir_figs,"Supp_Figure4AB.pdf"))
par(mfrow=c(1,2))
library(corrplot)
for(j in c("longterm,muscle","acute,muscle")){
  corrplot(gene_overlaps[[j]],is.corr = F,method="number",type = "upper",cl.length = 5,
           cl.cex = 0.8,cl.ratio = 0.3,bg = "gray",mar=c(1, 0, 1, 0),number.cex=0.5)
}
dev.off()

bg = unique(c(unlist(sapply(all_meta_analysis_res,names))))
gs = gene_sets_per_cov
gs = gs[sapply(gs,length)>10]

# go_res = run_topgo_enrichment_fisher(
#   gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
# go_res1 = go_res[go_res$Annotated < 1500,]
# go_res1$classicFisher[is.na(as.numeric(go_res1$classicFisher))] = 1e-30
# go_res1 = go_res[go_res1$Significant > 2,]
# go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
# go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
# go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
# table(as.character(go_res_fdr$setname))
# get_most_sig_enrichments_by_groups(go_res_fdr,num=2)[,1:4]
# go_enrichments_by_cov_fdr = go_res_fdr
# 
# reactome_pathways_by_cov = run_reactome_enrichment_analysis(gs,universe=bg)
# reactome_pathways_by_cov1 = reactome_pathways_by_cov[reactome_pathways_by_cov$Count>2,]
# ps = reactome_pathways_by_cov1$pvalue
# qs = p.adjust(ps,method="fdr")
# reactome_pathways_by_cov1$qvalue = qs
# reactome_pathways_by_cov_fdr = reactome_pathways_by_cov1[qs <= 0.1,]
# table(reactome_pathways_by_cov_fdr[,1])
# get_most_sig_enrichments_by_groups(reactome_pathways_by_cov_fdr,pcol="pvalue",num = 2)[,c(1,3)]
# reactome_pathways_by_cov_fdr[,1] = as.character(reactome_pathways_by_cov_fdr[,1])
# 
# save(gs,bg,go_enrichments_by_cov_fdr,reactome_pathways_by_cov_fdr,
#      file = paste(out_dir_rdata,"covariate_sets_enrichments.RData",sep=""))
# table(as.character(go_enrichments_by_cov_fdr[,1]))
# table(as.character(reactome_pathways_by_cov_fdr[,1]))
# reactome_pathways_by_cov_fdr[grepl("blood",reactome_pathways_by_cov_fdr[,1]),]
# go_enrichments_by_cov_fdr[grepl("blood",go_enrichments_by_cov_fdr[,1]),]

############################################################################
############################################################################
############################################################################
# Validation analyses of the selected gene sets
# TODO: rerun screen

# # 1. Replication analysis
# # Load SCREEN's results
# scr_path = "~/Desktop/MoTrPAC/PA_database/screen_res/"
# pvals_files = list.files(scr_path)
# pvals_files = pvals_files[grepl("pvals.txt",pvals_files)]
# screen_results = list()
# for (ff in pvals_files){
#   currname = gsub("_pvals.txt","",ff)
#   screen_results[[currname]] = list()
#   for(m in c("bum","znormix")){
#     outfile = paste(scr_path,currname,"_",m,".txt",sep="")
#     screen_results[[currname]][[m]] = read.delim(outfile,row.names = 1,header=T)
#   }
# }
# names(screen_results) = gsub("_",",",names(screen_results))
# save(screen_results,file=paste(scr_path,"screen_results.RData",sep=""))
# # Look at the FDR values of the selected gene sets vs. the others
# par(mfrow=c(3,2),mar=c(4,4,1,1),
#     cex.lab=1.1,cex.axis=1,cex=1,cex.main=1.05,lwd=0.5,mgp=c(2.1,1,0),
#     las=3)
# for(nn in names(analysis2selected_genes)[1:3]){
#   x = screen_results[[nn]]$bum
#   s = names(analysis2selected_genes[[nn]])
#   s_c = setdiff(rownames(x),s)
#   l1 = list();l2=list()
#   for(j in 1:8){
#     l1[[as.character(j+1)]] = x[s,j]
#     l2[[as.character(j+1)]] = x[s_c,j]
#   }
#   boxplot(l1,las=2,col="blue",ylab="local fdr",main=paste(nn,": selected",sep=""),
#           pch=20,ylim=c(0,1),xlab="Number of cohorts")
#   abline(h = 0.2,lty=2,col="red",lwd=1.5)
#   boxplot(l2,las=2,ylab="local fdr",main=paste(nn,": other",sep=""),
#           pch=20,ylim=c(0,1),xlab="Number of cohorts")
#   abline(h = 0.2,lty=2,col="red",lwd=1.5)
# }

# 2. Look at effects in untrained - excluded untrained others

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

pdf(paste(out_dir_figs,"Figure2E.pdf"))
par(mar=c(4,8,4,4))
names(l) = c("untrained","exercise","untrained","exercise",
             "untrained","exercise")
cols = c("red","red","cyan","cyan","green","green")
par(cex.lab=2.3,cex.axis=1.6)
boxplot(l,las=2,col=cols,horizontal=T,pch=20,ylim=c(0,2))
legend(x=0.5,y=4.5,c("muscle, acute","blood, acute","muscle, long-term"),
       fil=c("red","cyan","green"),cex=1.6)
dev.off()

# # 3. GO enrichments of the whole sets
bg = unique(c(unlist(sapply(simple_REs,names))))
gs = lapply(analysis2selected_genes,names)
go_res = run_topgo_enrichment_fisher(
  gs,bg,go_dags = "BP",go_term_size = 20,go_max_size = 200)
go_res1 = go_res[go_res$Annotated < 1500,]
go_res1$classicFisher[is.na(as.numeric(go_res1$classicFisher))] = 1e-30
go_res1 = go_res[go_res1$Significant > 3,]
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
go_res1$go_qvals = p.adjust(as.numeric(go_res1$classicFisher),method='fdr')
go_res_fdr = go_res1[go_res1$go_qvals < 0.1,]
table(go_res_fdr$setname)
gene_group_enrichments = go_res
gene_group_enrichments_fdr = go_res_fdr
get_most_sig_enrichments_by_groups(gene_group_enrichments_fdr,num=5)[,1:4]

save(gs,bg,gene_group_enrichments_fdr,
     file = paste(out_dir_rdata,"large_gene_sets_enrichments.RData",sep=""))
table(as.character(gene_group_enrichments_fdr[,1]))
gene_group_enrichments_fdr[grepl("blood",gene_group_enrichments_fdr[,1]),]

############################################################################
############################################################################
############################################################################
# Some figures

pdf(paste(out_dir_figs,"Supp_Figure3.pdf"))
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
dev.off()

# Examples of genes for the paper

# Acute, muscle:
curr_genes = analysis2selected_genes$`acute,muscle`
#PGC1a
gene = "10891"

# #BHLHE40
# gene = "8553"
# # SMAD3
# gene = "4088"
# # FOXO1
# gene = "2308"
# # ID1
# gene = "3397"
# # VEGFA
# gene = "7422"
# # HES1
# gene = "3280"
# NNT
# gene = "23530"

gene_name = entrez2symbol[[gene]]
gdata = meta_reg_datasets$`acute,muscle`[[gene]]
gdata = gdata[order(gdata$time),]
curr_times = rep("0-1h",nrow(gdata))
curr_times[gdata$time==2] = "2-5h"
curr_times[gdata$time==3] = ">20h"
gdata$training = gsub("endurance","EN",gdata$training)
gdata$training = gsub("resistance","RE",gdata$training)
gdata$V1 = gsub("GE_","",gdata$V1)
slabels = paste(gdata$V1,gdata$training,curr_times,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
analysis1[[gene]][[1]]$mod_p
aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
pdf(paste(out_dir_figs,"Figure2B.pdf"))
forest(x = gdata$yi,sei = gdata$sdd,slab = slabels,showweights = F,
       main=gene_name,annotate = F,cex = 1.3,xlab = "",
       pch = 19,top=0.1,col="darkgray",cex.axis=1.4)
dev.off()

pdf(paste(out_dir_figs,"Figure2C.pdf"))
par(cex.lab=1.5,cex.axis=1.4,mar = c(6,6,6,6))
boxplot(yi~time,data=gdata,
        main=paste(gene_name,": by time (aic diff:",format(aic_diff,digits=3),")",sep=""),
        xlab="Time window",names=c("0-1h","2-5h",">20h"),ylab="")
dev.off()

# Other genes: LPL, CPT1B, SMAD3, ACTN3, VEGFA, FOXO1 and IL6R
genes = c("1375","4023","4088","89","7422","2308","3570")
# Validated genes
library("metaviz")
validated_genes = c(
  "SMAD3" = "4088",
  "ID1" = "3397",
  "NR4A1" = "3164",
  "HES1" = "3280",
  "SCN2B" = "6327",
  "SLC25A25" = "114789",
  "PPARGC1A" = "10891",
  "NNT" = "23530",
  "SH3KBP1" = "30011"
)
for(gene in validated_genes){
  gene_name = entrez2symbol[[gene]]
  gdata = meta_reg_datasets$`acute,muscle`[[gene]]
  gdata = gdata[order(gdata$time),]
  curr_times = rep("0-1h",nrow(gdata))
  curr_times[gdata$time==2] = "2-5h"
  curr_times[gdata$time==3] = ">20h"
  gdata$training = gsub("endurance","EN",gdata$training)
  gdata$training = gsub("resistance","RE",gdata$training)
  gdata$V1 = gsub("GE_","",gdata$V1)
  slabels = paste(gdata$V1,gdata$training,curr_times,sep=",")
  analysis1 = all_meta_analysis_res$`acute,muscle`
  analysis1[[gene]][[1]]$mod_p
  pdf(paste(out_dir_figs,gene_name,".pdf",sep=""))
  forest(x = gdata$yi,sei = gdata$sdd,slab = slabels,showweights = F,
         main=gene_name,annotate = F,cex = 1.3,xlab = "",
         pch = 19,top=0.1,col="darkgray",cex.axis=1.4)
  dev.off()
}

# # Other blood genes
# genes = c("3482","9172","57105","1178","5004","5005")
# par(mfrow=c(2,2))
# for(gene in genes[5:6]){
#   gene_name = entrez2symbol[[gene]]
#   curr_m = simple_REs$`acute,blood`[[gene]]
#   gdata = meta_reg_datasets$`acute,blood`[[gene]]
#   curr_m$slab.null = F
#   curr_times = rep("0-1h",nrow(gdata))
#   curr_times[gdata$time==2] = "2-5h"
#   curr_times[gdata$time==3] = ">20h"
#   curr_m$slab = paste(gdata$training,curr_times,sep=",")
#   analysis1 = all_meta_analysis_res$`acute,blood`
#   analysis1[[gene]][[1]]
#   aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
#   forest(curr_m,main=paste(gene_name),annotate = T)
# }

# COL4A1 in longterm muscle
gene = "1282"
gene_name = entrez2symbol[[gene]]
all_meta_analysis_res$`longterm,muscle`[[gene]]
gdata = meta_reg_datasets$`longterm,muscle`[[gene]]
gdata = gdata[order(gdata$time),]
curr_times = rep("",nrow(gdata))
curr_times[gdata$time==2] = ",> 150 days"
gdata$V1 = gsub("GE_","",gdata$V1)
gdata$training = gsub("endurance","EN",gdata$training)
gdata$training = gsub("resistance","RE",gdata$training)
slabels = paste(gdata$training,curr_times,sep="")
slabels = paste(gdata$V1,slabels,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
pdf(paste(out_dir_figs,"Figure2D.pdf"))
forest(x = gdata$yi,sei = gdata$sdd,slab = slabels,showweights = F,
       main=gene_name,annotate = F,cex = 1.3,xlab = "",
       pch = 19,top=0.1,col="darkgray",cex.axis=1.4)
dev.off()

# # Plots with sex info - longterm muscle
# genes = c("8897","567","50","11217")
# par(mfrow=c(2,2))
# for(gene in genes){
#   gene_name = entrez2symbol[[gene]]
#   curr_m = simple_REs$`longterm,muscle`[[gene]]
#   all_meta_analysis_res$`longterm,muscle`[[gene]]
#   curr_m$I2
#   gdata = meta_reg_datasets$`longterm,muscle`[[gene]]
#   curr_m$slab.null = F
#   curr_times = rep("",nrow(gdata))
#   curr_times[gdata$time==2] = ",> 150 days"
#   curr_m$slab = paste("prop males: ",format(gdata$prop_males,digits=2),sep="")
#   analysis1 = all_meta_analysis_res$`longterm,muscle`
#   aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
#   currmain=paste(gene_name,": by sex (aic diff:",format(aic_diff,digits=3),")",sep="")
#   forest(curr_m,main=currmain,annotate = T)
# }

# # Genes in both muscle analyses: DDN, ID3 and PLEKHO1
# genes = c("23109","3399","51177")
# par(mfrow=c(3,1),mar=c(3,0,3,0))
# for(gene in genes){
#   gene_name = entrez2symbol[[gene]]
#   curr_m = simple_REs$`longterm,muscle`[[gene]]
#   all_meta_analysis_res$`longterm,muscle`[[gene]]
#   curr_m$I2
#   gdata = meta_reg_datasets$`longterm,muscle`[[gene]]
#   curr_m$slab.null = F
#   curr_m$slab = paste(gdata$training)
#   analysis1 = all_meta_analysis_res$`longterm,muscle`
#   aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
#   currmain=paste("LT,",gene_name," (aic diff:",format(aic_diff,digits=3),")",sep="")
#   forest(curr_m,main=currmain,annotate = T,addfit = F,showweights = T,top=0,cex=0.7,
#          transf=transf.ztor, digits=3, alim=c(-1,1))
# }

# par(mfrow=c(3,1),mar=c(3,0,3,0))
# for(gene in genes){
#   gene_name = entrez2symbol[[gene]]
#   curr_m = simple_REs$`acute,muscle`[[gene]]
#   all_meta_analysis_res$`acute,muscle`[[gene]]
#   curr_m$I2
#   gdata = meta_reg_datasets$`acute,muscle`[[gene]]
#   curr_m$slab.null = F
#   curr_m$slab = paste(gdata$training)
#   analysis1 = all_meta_analysis_res$`acute,muscle`
#   aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
#   currmain=paste("Acute,",gene_name," (aic diff:",format(aic_diff,digits=3),")",sep="")
#   forest(curr_m,main=currmain,annotate = T,addfit = F,showweights = T,top=0,cex=0.7,
#          transf=transf.ztor, digits=3, alim=c(-1,1))
# }


############################################################################
############################################################################
############################################################################
# Interpretation of the results: defining subgroups by clustering mean patterns
library(parallel)
# Some helper functions for reformatting the data
get_ts<-function(gdata){
  v = as.numeric(gdata$tstat)
  ns = paste(gdata$V1,gdata$training,gdata$time,
             format(gdata$avg_age,digits=1),format(gdata$prop_males,digits=1),sep=";")
  v = tapply(v,ns,mean)
  return(v)
}
get_t_matrix_from_list<-function(tstats){
  if("matrix" %in% class(tstats)){
    if(ncol(tstats)>nrow(tstats)){
      tstats = t(tstats)
    }
    return(tstats)
  }
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
dataset_tstats = mclapply(datasets,function(x)sapply(x,get_ts),mc.cores = 4)
mean_effect_matrices = mclapply(dataset_tstats,get_t_matrix_from_list,mc.cores = 4)
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
get_num_clusters_wss_kmeans<-function(data,k.max=10,wss_imp_thr=0.7,seed = 1){
  set.seed(seed)
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
          # assumption: we do not have more than 5 patterns per cluster
          currk = min(currk,5)
        }
        else{currk=1}
        m_kmeans = kmeans(m_processed,centers = currk)
        m_kmeans = m_kmeans$cluster
      }
      else{
        m = as.matrix(rowMeans(m),ncol=1)
        colnames(m)[1] = "base,mean_t"
        m_kmeans = as.numeric(simple_RE_beta[[nn]][rownames(m)]>0)+1
        names(m_kmeans) = rownames(m)
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

save(gs,bg,gene_subgroup_enrichments_fdr,
     reactome_pathways_subgroups_fdr,
     file = paste(out_dir_rdata,"gene_subgroups_enrichments.RData",sep=""))
table(as.character(gene_subgroup_enrichments_fdr[,1]))
gene_subgroup_enrichments_fdr[grepl("blood",gene_subgroup_enrichments_fdr[,1]),]
gene_subgroup_enrichments_fdr[grepl("avg_age",gene_subgroup_enrichments_fdr[,1]),]
table(as.character(reactome_pathways_subgroups_fdr[,1]))
reactome_pathways_subgroups_fdr[grepl("blood",reactome_pathways_subgroups_fdr[,1]),]
reactome_pathways_subgroups_fdr[grepl("avg_age",reactome_pathways_subgroups_fdr[,1]),]

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
for(set_name in names(gene_subgroups)){
  # if(grepl("acute,muscle,time,",set_name)){next}
  arr = strsplit(set_name,split=",")[[1]]
  table_name = paste(arr[1:2],collapse=",")
  table_name2 = paste(arr[1:3],collapse=",")
  set_genes = gene_subgroups[[set_name]]
  if(length(set_genes)<3){next}
  mat = gene_t_patterns[[table_name2]][set_genes,]
  mat = mat[,-ncol(mat)]
  if(grepl("base_",set_name)){
    mat = mean_effect_matrices[[table_name]][set_genes,]
  }
  curr_gene_names = entrez2symbol[rownames(mat)]
  missing_names = sapply(curr_gene_names, length) == 0
  curr_gene_names[missing_names] = rownames(mat)[missing_names]
  rownames(mat) = unlist(curr_gene_names)
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
  curr_main = paste(set_name,"\n",sep="")
  if(num_enrichments > 0){
    curr_enrichments = curr_enrichments[1:num_enrichments,]
    print(dim(curr_enrichments))
    for(j in 1:nrow(curr_enrichments)){
      curr_e = paste(curr_enrichments[j,1]," (p=",curr_enrichments[j,2],")",sep="")
      curr_main = paste(curr_main,curr_e,"\n",sep="")
    }
  }
  curr_main = gsub("\\.\\.\\.","",curr_main)
  cex_genes = 0.4
  if(nrow(mat)<40){cex_genes = 0.7}
  if(nrow(mat)<20){cex_genes = 1}
  if(nrow(mat)<10){cex_genes = 1.5}
  pdf(paste(out_dir_figs,gsub(",|;","_",set_name),".pdf",sep=""))
  par(cex.main=0.7)
  heatmap.2(mat,trace = "none",scale = "none",Colv = T,col=bluered,
            cexRow = cex_genes,main=curr_main,
            Rowv = T,srtCol=45,hclustfun = hclust_func,density.info="none",
            key.title = NA,keysize = 1.1,key.xlab = "t-statistic",
            key.par = list("cex.axis"=1.1),margins = c(10,10))
  # heatmap.2(mat,trace = "none",scale = "none",Colv = T,col=bluered,
  #           cexRow = cex_genes,main="",key = F,
  #           Rowv = T,srtCol=45,hclustfun = hclust_func,density.info="none",
  #           key.title = NA,keysize = 1.1,key.xlab = "t-statistic",
  #           key.par = list("cex.axis"=1.1),margins = c(10,10))
  dev.off()
}

# Specifically for longterm muscle base models:
makeRects <- function(m,lwd=2){
  coords = expand.grid(nrow(m):1, 1:ncol(m))[m,]
  xl=coords[,2]-0.49
  yb=coords[,1]-0.49
  xr=coords[,2]+0.49
  yt=coords[,1]+0.49
  rect(xl,yb,xr,yt,border="black",lwd=lwd)
}

set_name = "longterm,muscle,base_model,1"
table_name = "longterm,muscle"
set_genes = gene_subgroups[[set_name]]
length(set_genes)
mat = mean_effect_matrices[[table_name]][set_genes,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,apply(mat,2,sd)>0]
mat = mat[,apply(mat==0,2,sum)/nrow(mat) < 0.5]
mat = mat[apply(mat==0,1,sum)/ncol(mat) < 0.5,]
mat[mat>4]=4;mat[mat< -4]=-4
colnames(mat) = NULL;
pdf(paste(out_dir_figs,"Figure4B.pdf"))
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 1.2,
          add.expr={makeRects(mat==0)},Rowv = F,margins = c(5,8))
dev.off()
# get the top enrichments of the set
curr_gos = gene_subgroup_enrichments_fdr[gene_subgroup_enrichments_fdr[,1]==set_name,]
curr_pathways = reactome_pathways_subgroups_fdr[reactome_pathways_subgroups_fdr[,1]==set_name,]
# get node file for cytoscape
lonterm_muscle_base_model_effects = cbind(
  set_genes,
  unlist(entrez2symbol[set_genes]),
  sapply(all_meta_analysis_res$`longterm,muscle`[set_genes],
         function(x)x[[1]]$coeffs[1,1])
)

set_names = c("longterm,muscle,base_model,1",
              "longterm,muscle,base_model,2")
table_name = "longterm,muscle"
set_genes = unlist(gene_subgroups[set_names])
length(set_genes)
mat = mean_effect_matrices[[table_name]][set_genes,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,apply(mat,2,sd)>0]
mat = mat[,apply(mat==0,2,sum)/nrow(mat) < 0.5]
mat = mat[apply(mat==0,1,sum)/ncol(mat) < 0.5,]
mat[mat>4]=4;mat[mat< -4]=-4
colnames(mat) = NULL;
pdf(paste(out_dir_figs,"Figure4A.pdf"))
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 1.2,
          add.expr={makeRects(mat==0)},Rowv = F,margins = c(5,8))
dev.off()
# get the top enrichments of the set
curr_gos = gene_subgroup_enrichments_fdr[gene_subgroup_enrichments_fdr[,1]==set_name,]
curr_pathways = reactome_pathways_subgroups_fdr[reactome_pathways_subgroups_fdr[,1]==set_name,]

lonterm_muscle_base_model_effects = rbind(
  lonterm_muscle_base_model_effects, 
  cbind(
    set_genes,
    unlist(entrez2symbol[set_genes]),
    sapply(all_meta_analysis_res$`longterm,muscle`[set_genes],
           function(x)x[[1]]$coeffs[1,1])
  )
)
colnames(lonterm_muscle_base_model_effects) = c("entrez","Symbol","fchange")
write.table(lonterm_muscle_base_model_effects,
            file=paste(out_dir_figs,"Figure4_longterm_genes_nodes.txt"),
            quote = F,row.names = F,col.names = T,sep="\t")

###################################################################
# The different patterns of acute muscle, time-dependent response
###################################################################
pref = "acute,muscle,time,\\d"
clusters = names(gene_subgroups)[grepl(pref,names(gene_subgroups))]
cols = c("blue","black","red","green")
names(cols) = clusters

pdf(paste(out_dir_figs,"Figure3A.pdf"))
par(mfrow=c(2,2),cex.main=1.2)
for(set_name in sort(clusters)[c(3,4,1,2)]){
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
  curr_main = "" # consider keeping the main instead
  plot_with_err_bars(c("0-1h","2-5h",">20h"),
                     colMeans(mat),apply(mat,2,sd),col=cols[set_name],lwd=3,
                     main=curr_main,ylab = "Mean t-statistic",xlab="Time",cex.main=1.1,
                     cex.lab=1.2,cex.axis=1.3,arrow_col = cols[set_name])
  abline(h = 0,lty=2,col="black")
}
dev.off()

# redo the same figure, keep the labels
pdf(paste(out_dir_figs,"Figure3A_with_text.pdf"))
par(mfrow=c(2,2),cex.main=1.2)
for(set_name in sort(clusters)[c(3,4,1,2)]){
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
dev.off()

# Write input to DREM
m = gene_t_patterns$`acute,muscle,time`[,1:3]
rownames(m) = entrez2symbol[rownames(m)]
colnames(m) = c("1h","4h","24h")
m = rbind(colnames(m),m)
write.table(m,file="./drem/drem_acute_muscle_time.txt",
            sep="\t",col.names = F,row.names = T,quote=F)
# Write input for cytoscape
m = gene_subgroups
m = m[grepl("acute,muscle,time,",names(m))]
mm = c()
for(cluster_name in names(m)){
  mm = rbind(mm,
             cbind(m[[cluster_name]],rep(cluster_name,length(m[[cluster_name]]))))
}
write.table(mm,file=paste(out_dir_figs,"acute_muscle_time_subgroups.txt",sep=""),
            sep="\t",col.names = F,row.names = F,quote=F)

# pref = "acute,blood,time,\\d"
# clusters = names(gene_subgroups)[grepl(pref,names(gene_subgroups))]
# par(mfrow=c(1,3))
# cols = c("blue","red","green")
# names(cols) = clusters
# for(set_name in sort(clusters)[1:3]){
#   gene_set = gene_subgroups[[set_name]]
#   top_enrichments1 = reactome_pathways_subgroups_fdr[
#     reactome_pathways_subgroups_fdr[,1]==set_name,"Description"]
#   top_enrichments2 = gene_subgroup_enrichments_fdr[
#     gene_subgroup_enrichments_fdr[,1]==set_name,"Term"]
# 
#   # remove embryo enrichments
#   top_enrichments1 = top_enrichments1[!grepl("embryo",top_enrichments1)]
#   top_enrichments2 = top_enrichments2[!grepl("embryo",top_enrichments2)]
# 
#   enrichment1 = top_enrichments1[1]
#   if(any(grepl("muscle",top_enrichments1))){
#     enrichment1 = top_enrichments1[grepl("muscle",top_enrichments1)][1]
#   }
#   enrichment2 = top_enrichments2[1]
#   if(any(grepl("muscle",top_enrichments2))){
#     enrichment2 = top_enrichments2[grepl("muscle",top_enrichments2)][1]
#   }
#   if(length(top_enrichments1)> 1 && (is.na(enrichment2) || length(enrichment2)==0)){
#     enrichment2 = top_enrichments1[2]
#   }
#   if(is.na(enrichment1)){enrichment1 = NULL}
#   if(is.na(enrichment2)){enrichment2 = NULL}
#   enrichment1 = shorten_by_words(tolower(enrichment1))
#   enrichment2 = shorten_by_words(tolower(enrichment2))
#   new_set_name = paste(set_name," (",length(gene_set)," genes)",sep="")
#   curr_main = paste(new_set_name,enrichment1,enrichment2,sep="\n")
#   curr_main = gsub("\n\n","\n",curr_main)
#   print(curr_main)
#   mat = gene_t_patterns$`acute,blood,time`[gene_set,]
#   rownames(mat) = unlist(entrez2symbol[rownames(mat)])
#   mat = mat[,-ncol(mat)]
#   plot_with_err_bars(c("0-1h","2-5h",">20h"),arrow_col = cols[set_name],
#                      colMeans(mat),apply(mat,2,sd),lwd=2,xlab = "Time",
#                      main=curr_main,ylab = "Mean t-statistic",col=cols[set_name])
# }

# # Write input to DREM
# m = gene_t_patterns$`acute,blood,time`[,1:3]
# rownames(m) = entrez2symbol[rownames(m)]
# write.table(m,file="drem/acute_blood_time.txt",
#             sep="\t",col.names = T,row.names = T,quote=F)
# # Write input for cytoscape
# m = gene_subgroups
# m = m[grepl("acute,blood,time,",names(m))]
# mm = c()
# for(cluster_name in names(m)){
#   mm = rbind(mm,
#              cbind(m[[cluster_name]],rep(cluster_name,length(m[[cluster_name]]))))
# }
# write.table(mm,file="supp_tables/acute_blood_time_subgroups.txt",
#             sep="\t",col.names = F,row.names = F,quote=F)

###############################################
###############################################
###############################################
# Look at the covariate corrleations
###############################################
###############################################
###############################################
library(corrplot);library(lme4);library(lmerTest)

get_assoc_matrix<-function(gdata,ys = c("avg_age","prop_males","N")){
  n = length(ys)
  newgdata = c()
  alltimes = sort(unique(gdata$time))
  for(cohort in unique(gdata$V1)){
    currinds = which(gdata$V1==cohort)
    v = gdata[currinds[1],-c(1:2)]
    for(tt in alltimes){
      if(is.element(tt,gdata[currinds,"time"])){
        v = c(v,1)
      }
      else{
        v = c(v,0)
      }
      names(v)[length(v)] = paste("time_window",tt,sep="")
    }
    newgdata = rbind(newgdata,v)
  }
  rownames(newgdata) = NULL
  mode(newgdata) = "numeric"
  gdata = data.frame(newgdata)
  p = matrix(0,nrow=n,ncol = ncol(gdata),dimnames = list(ys,colnames(gdata)))
  betas = matrix(1,nrow=n,ncol = ncol(gdata),dimnames = list(ys,colnames(gdata)))
  rhos = matrix(1,nrow=n,ncol = ncol(gdata),dimnames = list(ys,colnames(gdata)))
  rhosp = matrix(0,nrow=n,ncol = ncol(gdata),dimnames = list(ys,colnames(gdata)))
  r2s = matrix(1,nrow=n,ncol = ncol(gdata),dimnames = list(ys,colnames(gdata)))
  for(y in ys){
    for(x in colnames(p)){
      if(x==y){next}
      form = as.formula(paste(y,"~",x))
      model = summary(lm(form,data=gdata))
      currp = model$coefficients[2,4]
      currbeta = model$coefficients[2,1]
      rhoTest = cor.test(gdata[,y],gdata[,x],method="spearman")
      p[y,x] = currp
      betas[y,x] = currbeta
      rhos[y,x] = rhoTest$estimate
      rhosp[y,x] = rhoTest$p.value
      r2s[y,x] = model$r.squared
    }
  }
  return(list(rhos=rhos,rhosp=rhosp,lmmp=p,lmmb=betas,r2s=r2s))
}

pdf(paste(out_dir_figs,"Supp_Figure1.pdf"))
par(mfrow=c(2,2))
for(nn in names(meta_reg_datasets)[c(1,3)]){
  dataset = meta_reg_datasets[[nn]]
  dataset_sizes = sapply(dataset,nrow)
  selected_dataset = which(dataset_sizes == max(dataset_sizes))[1]
  gdata = dataset[[selected_dataset]]
  gdata = gdata[,c("V1","time","training","avg_age","prop_males","N")]
  gdata$training = as.numeric(grepl("resis",gdata$training))
  gdata = gdata[,!apply(gdata,2,function(x)length(unique(x))==1)]
  gdata = gdata[!apply(is.na(gdata), 1,any),]
  cov_eval = get_assoc_matrix(gdata)
  corrplot(cov_eval$rhos[1:3,1:3],p.mat = cov_eval$rhosp[1:3,1:3],
           main=paste(nn," - rho (Spearman)"),mar=c(0,2,4,0),cex.main=1.5,tl.cex = 1.4)
  corrplot(cov_eval$r2s,p.mat = cov_eval$lmmp,main=paste(nn," - lm test"),
           mar=c(0,2,4,0),cl.lim = c(0,1),cex.main=1.5,tl.cex = 1.4)
  print(cov_eval$r2s)
  print(median(cov_eval$r2s))
}
dev.off()

###############################################
###############################################
###############################################
# Prepare supplementary tables for the results
###############################################
###############################################
###############################################

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

system(paste("mkdir","supp_tables"))
supp_path = paste(getwd(),"supp_tables/",sep="/")
supp_file = paste(supp_path,"SupplementaryTables.xlsx",sep="/")
library(readxl)

# Parse the information from the datasets used in the meta-analysis
meta_reg_datasets_cohort_names = 
  lapply(meta_reg_datasets,function(x)unique(unlist(lapply(x,function(y)y$V1))))
untrained_datasets_cohort_names = 
  lapply(untrained_datasets,function(x)unique(unlist(lapply(x,function(y)y$V1))))


get_gene_table_with_cohort_info<-function(cohort_name,gtables){
  for(g in gtables){
    if(cohort_name %in% g[,1]){
      v = g[which(g[,1]==cohort_name)[1],]
      names(v) = colnames(g)
      return((v))
    }
  }
  return(NULL)
}

metadata_row_for_supp_table<-function(cohort_name,cohort_metadata,
                                      gtable_row,sample2subject){
  samp2time = cohort_metadata[[cohort_name]]$times
  curr_samps = samp2time[cohort_metadata[[cohort_name]]$gsms]
  curr_samps = curr_samps[!is.na(curr_samps)]
  Nsample = length(curr_samps)
  curr_samps[curr_samps==min(samp2time)]="Pre"
  Nsubject=length(unique(sample2subject[cohort_metadata[[cohort_name]]$gsms]))
  tps = paste(unique(curr_samps),collapse=",")
  res = c(cohort_name,"GSE"= cohort_metadata[[cohort_name]]$gse,
          "Tissue" = cohort_metadata[[cohort_name]]$tissue,
          "Training"= cohort_metadata[[cohort_name]]$training,
          "Nsample"= Nsample,"Nsubject" = Nsubject,
          "Time_points"= tps,
          "Avg_age"= as.numeric(gtable_row["avg_age"]),
          "Prop_males" = as.numeric(gtable_row["prop_males"]),
          "Additional Info" = cohort_metadata[[cohort_name]]$additional_info)
  return(res)
}
m1 = c()
# Add the meta-anaysis used cohorts
for(analysis_name in names(meta_reg_datasets_cohort_names)){
  curr_metadata = acute_metadata
  curr_sample_meta = acute_sample_meta
  if(grepl("longterm",analysis_name)){
    curr_metadata = longterm_metadata
    curr_sample_meta = longterm_sample_meta
  }
  for(cohort_name in meta_reg_datasets_cohort_names[[analysis_name]]){
    curr_gtable_row = get_gene_table_with_cohort_info(
      cohort_name,meta_reg_datasets[[analysis_name]]) 
    curr_info = metadata_row_for_supp_table(
      cohort_name,curr_metadata,curr_gtable_row,curr_sample_meta$subject)
    curr_info["Analysis_group"] = analysis_name
    names(curr_info)[1] = "Cohort_ID"
    m1 = rbind(m1,curr_info)
  }
}
# add the untrained cohorts
for(analysis_name in names(untrained_datasets_cohort_names)){
  curr_metadata = acute_metadata
  curr_sample_meta = acute_sample_meta
  if(grepl("longterm",analysis_name)){
    curr_metadata = longterm_metadata
    curr_sample_meta = longterm_sample_meta
  }
  for(cohort_name in untrained_datasets_cohort_names[[analysis_name]]){
    curr_gtable_row = get_gene_table_with_cohort_info(
      cohort_name,untrained_datasets[[analysis_name]]) 
    curr_info = metadata_row_for_supp_table(
      cohort_name,curr_metadata,curr_gtable_row,curr_sample_meta$subject)
    curr_info["Analysis_group"] = paste(analysis_name,"untrained",sep=",")
    names(curr_info)[1] = "Cohort_ID"
    m1 = rbind(m1,curr_info)
  }
}
# add the unused datasets
for(nn in names(acute_metadata)){
  if(nn %in% m1[,1]){next}
  curr_gtable_row = get_gene_table_with_cohort_info(nn,acute_gene_tables) 
  if(is.null(curr_gtable_row)){
    curr_gtable_row = c("avg_age" = NA, "prop_males" = NA)
  }
  curr_info = metadata_row_for_supp_table(
    nn,acute_metadata,curr_gtable_row,acute_sample_meta$subject)
  curr_info["Analysis_group"] = "Not used in meta-analysis"
  names(curr_info)[1] = "Cohort_ID"
  m1 = rbind(m1,curr_info)
}

for(nn in names(longterm_datasets)){
  if(nn %in% m1[,1]){next}
  curr_gtable_row = get_gene_table_with_cohort_info(nn,longterm_gene_tables) 
  if(is.null(curr_gtable_row)){
    curr_gtable_row = c("avg_age" = NA, "prop_males" = NA)
  }
  curr_info = metadata_row_for_supp_table(
    nn,longterm_metadata,curr_gtable_row,longterm_sample_meta$subject)
  curr_info["Analysis_group"] = "Not used in meta-analysis"
  names(curr_info)[1] = "Cohort_ID"
  m1 = rbind(m1,curr_info)
}

supp_table_1_all_cohorts = m1
# write.xlsx(supp_table_1_all_cohorts,file=supp_file,sheetName = "STable1",row.names = F)
write.table(supp_table_1_all_cohorts,file=paste(supp_path,"STable1.txt",sep="")
            ,row.names = F,quote=F,sep="\t")
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
# write.xlsx(supp_table_genes,file=supp_file,sheetName = "STable2",row.names = F,append = T)
write.table(supp_table_genes,file=paste(supp_path,"STable2.txt",sep="")
            ,row.names = F,quote=F,sep="\t")

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
  # write.xlsx(m,file=supp_file,
  #              sheetName = paste("STable",sheet_counter,"_",nn,sep=""),
  #              row.names = F,append = T,col.names = T)
  write.table(m,sep="\t",
             file = paste("supp_tables/STable",sheet_counter,"_",nn,".txt",sep=""),
             row.names = F,col.names = T,quote=F)
  sheet_counter = sheet_counter+1
}

# GO enrichments of subgroups
supp_table_enrichments = gene_subgroup_enrichments_fdr[,c(1:4,9:10,8)]
colnames(supp_table_enrichments)[6] = "q-value"
colnames(supp_table_enrichments)[5] = "Genes"
colnames(supp_table_enrichments)[1] = "Discovered in"
# write.xlsx(supp_table_enrichments,file=supp_file,
#            sheetName = paste("STable",sheet_counter,"_GO_enrichments",sep=""),
#            row.names = F,append=T)
write.table(supp_table_enrichments,file=paste(supp_path,"STable7.txt",sep="")
            ,row.names = F,quote=F,sep="\t")


# Same for Reactome
sheet_counter = sheet_counter+1
supp_table_enrichments = reactome_pathways_subgroups_fdr[,c(1,3,8,9,6)]
colnames(supp_table_enrichments)[3] = "q-value"
colnames(supp_table_enrichments)[1] = "Discovered in"
# write.xlsx(supp_table_enrichments,file=supp_file,
#            sheetName = paste("STable",sheet_counter,"_Reactome_enrichments",sep=""),
#            row.names = F,append=T)
write.table(supp_table_enrichments,file=paste(supp_path,"STable8.txt",sep="")
            ,row.names = F,quote=F,sep="\t")

# merge all gsea results into one table
all_fdr_corrected_gsea_results = c()
for(nn in names(gsea_reactome_results)){
  m = gsea_reactome_results[[nn]]
  m = as.data.frame(m)
  m$analysis = nn
  all_fdr_corrected_gsea_results = rbind(all_fdr_corrected_gsea_results,m)
}
all_fdr_corrected_gsea_results$leadingEdge = sapply(all_fdr_corrected_gsea_results$leadingEdge,
                                                    paste,collapse=";")
all_fdr_corrected_gsea_results = as.matrix(all_fdr_corrected_gsea_results)
all_fdr_corrected_gsea_results = all_fdr_corrected_gsea_results[,
          c(ncol(all_fdr_corrected_gsea_results),1:(ncol(all_fdr_corrected_gsea_results)-1))]
write.table(all_fdr_corrected_gsea_results,file=paste(supp_path,"STable9.txt",sep="")
            ,row.names = F,quote=F,sep="\t")

# save the workspace
save.image(file=paste(out_dir,"meta_analysis_interpretation_results.RData",sep=""))

# ###############################################
# ###############################################
# ###############################################
# # Prepare gene-csv files for the portal
# ###############################################
# ###############################################
# ###############################################
# load("meta_analysis_interpretation_results.RData")
# 
# write_to_bucket<-function(m,fname,bucket){
#   write.csv(m,row.names = F,quote = T,file=fname)
#   system(paste("~/google-cloud-sdk/bin/gsutil cp",fname,bucket))
#   system(paste("rm",fname))
# }
# 
# for(nn in names(meta_reg_datasets)){
#   bucket = paste("gs://bic_data_analysis/meta_analysis_human/",nn,"/",sep="")
#   bucket = gsub(",","_",bucket)
#   genes = names(meta_reg_datasets[[nn]])
#   curr_gene_table = c()
#   for(gene in genes){
#     gname = entrez2symbol[[gene]]
#     gdata = meta_reg_datasets[[nn]][[gene]]
#     base_model = simple_REs[[nn]][[gene]]
#     if(length(base_model)==1 && is.na(base_model)){next}
#     i2 = base_model$I2
#     tau2 = base_model$tau2
#     p_model = base_model$pval
#     selected = is.element(gene,names(analysis2selected_genes[[nn]]))
# 
#     if(grepl("acute",nn)){
#       curr_times = rep("0-1h",nrow(gdata))
#       curr_times[gdata$time==2] = "2-5h"
#       curr_times[gdata$time==3] = ">20h"
#     }
#     else{
#       curr_times = rep("0-150 days",nrow(gdata))
#       curr_times[gdata$time==2] = ">150 days"
#     }
#     
#     gdata$time = curr_times
#     meta_reg_analysis = all_meta_analysis_res[[nn]][[gene]]
#     aic_diff = meta_reg_analysis[[1]]$aic_c - meta_reg_analysis$`simple:base_model`$aic_c
#     
#     selected_model_name = names(meta_reg_analysis)[1]
#     if(grepl("base_model",selected_model_name)){
#       selected_model_name = ":base_model (simple RE)" # add ':' just for the split later
#     }
#     else{
#       p_model = meta_reg_analysis[[1]]$mod_p
#     }
#     
#     selected_model_name = strsplit(selected_model_name,":")[[1]][2]
#     selected_model_name = gsub("avg_","",selected_model_name)
#     selected_model_name = gsub(";",",",selected_model_name)
#     
#     curr_gene_table = rbind(curr_gene_table,
#       c(gname,gene,selected_model_name,i2,tau2,p_model,aic_diff,selected))
#     
#     curr_gfile = paste(gname,".csv",sep="")
#     write_to_bucket(gdata,curr_gfile,bucket)
#   }
#   colnames(curr_gene_table) = c("Symbol","Entrez","SelectedModel",
#                                 "I2","Tau2","P-value","AICcDiff","Selected?")
#   
#   fname = "gene_stats.csv"
#   write_to_bucket(curr_gene_table,fname,bucket)
# }


