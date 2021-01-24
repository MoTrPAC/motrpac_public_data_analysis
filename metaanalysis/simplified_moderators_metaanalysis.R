#' This is the main analysis script that takes the model
#' selection analysis results, applied filters to define the set
#' of genes and then run multiple analyses for interpretation;
#' This code also generates most of the display items for the paper.
#' One exception is the histograms and qq-plots for supplementary figuren 2,
#' these are created by the script simplified_mod_data_prep_and_naive_RE.R

###########################################################
# Load libraries and the data
library(org.Hs.eg.db);library(parallel);
library(metafor);library(data.table);library(fgsea)
source('~/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)
symbol2entrez = as.list(org.Hs.egSYMBOL2EG)

# this is the main directory - as you can see from the names, our preliminary
# data collection was done in 2018, but the full analysis, including adding new
# datasets and mitigating some issues related to the low number of studies, were
# addressed in early 2020.
out_dir = "~/Desktop/MoTrPAC/project_release_feb_2018/revision_feb_2020/"
out_dir_figs = paste(out_dir,"figures/",sep="")
try(dir.create(out_dir_figs))
out_dir_rdata = paste(out_dir,"rdata/",sep="")
try(dir.create(out_dir_rdata))
try(dir.create(paste(out_dir,"supp_tables/",sep="")))
setwd(out_dir)

# Load the results (run on stanford's hpc) instead of running the code above
load("meta_analysis_results.RData")
# same obj as above but with R2 scores
load("./all_meta_analysis_res_withR2.RData")
load("workspace_before_rep_analysis.RData")
load("meta_analysis_input.RData")

###########################################################
# Set thresholds for gene selection
I2_thr = 50
AIC_diff_thr = 5
ACUTE_beta_thr = 0.1
LONGTERM_beta_thr = 0.1
P_thr = 1e-03

###########################################################
# Helper functions for the analyses below

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

# Analysis for correlation or association among the covariates in an
# analysis. The results are presented in Supp Fig 3
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

###########################################################
# Algorithm for selecting genes from each meta-regression analysis
analysis2selected_genes = list()
analysis2selected_genes_stats = list()
all_pvals = c()
for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]]
  pvals = sapply(analysis1,function(x)x[[1]]$mod_p)
  # get the R2s of the selected models
  R2s = c()
  ngenes_with_R2 = 0
  for(gene in names(pvals)){
    obj = analysis1[[gene]]
    R2s[gene] = NA
    if(!is.na(obj[[1]]$R2) && !is.null(obj[[1]]$R2)){
      R2s[gene] = analysis1[[gene]][[1]]$R2
    }
    for(jj in 1:length(obj)){
      if(!is.na(obj[[jj]]$R2) && !is.null(obj[[jj]]$R2)){
        ngenes_with_R2 = ngenes_with_R2 + 1
        break
      }
    }
  }
  
  all_pvals = c(all_pvals,pvals)
  i2s = simple_RE_I2s[[nn]][names(pvals)]
  i2s[is.na(i2s)] = 100
  # separate into two gene sets: 
  # those that passed the aic diff test vs. those that did not
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
  model2beta[is_base_model] = sapply(analysis1[is_base_model],function(x)get_simple_model_beta(x)>curr_effect_thr) 
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
    coeffs_v, # details about the coefficients
    R2s[curr_selected_genes]
  )
  colnames(m)= c("Entrez","Symbol","Group","Model pvalue","AICc diff","Coefficients","R2")
  analysis2selected_genes_stats[[nn]] = m
}
sapply(analysis2selected_genes,length)

# Sanity checks and tests: our p-value threshold is similar to the FDR adjusted one
fdr0.01_thr = max(all_pvals[p.adjust(all_pvals,method = "fdr")<0.01])

# Examine some Overlaps
intersect(analysis2selected_genes_stats[[1]][,2],analysis2selected_genes_stats[[2]][,2])
intersect(analysis2selected_genes_stats[[1]][,2],analysis2selected_genes_stats[[3]][,2])

###########################################################
# Add GSEA enrichments: use the inferred beta values
# use GTEx data to filter genes that are not expressed in tissues
gtex_mean_tpm = fread("./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                      skip=2,sep="\t",data.table = F,stringsAsFactors = F)
gtex_entrez = symbol2entrez[gtex_mean_tpm[,2]]
inds = sapply(gtex_entrez,length)>1
gtex_entrez[inds] = sapply(gtex_entrez[inds],function(x)x[1])
gtex_entrez = unlist(gtex_entrez)
gtex_mean_tpm = gtex_mean_tpm[gtex_mean_tpm[,2] %in% names(gtex_entrez),]
gtex_mean_tpm[,1] = gtex_entrez[gtex_mean_tpm[,2]]
colSums(gtex_mean_tpm > 5 )
colnames(gtex_mean_tpm)[grepl("blood",colnames(gtex_mean_tpm),ignore.case = T)]
colnames(gtex_mean_tpm)[grepl("musc",colnames(gtex_mean_tpm),ignore.case = T)]
blood_gtex_genes = gtex_mean_tpm[gtex_mean_tpm[, "Whole Blood"] > 1,1]
muscle_gtex_genes = gtex_mean_tpm[gtex_mean_tpm[, "Muscle - Skeletal"] > 1,1]
length(blood_gtex_genes);length(muscle_gtex_genes)

gsea_input_scores = list()
for(nn in names(simple_RE_beta)){
  currbetas = simple_RE_beta[[nn]]
  curri2 = simple_RE_I2s[[nn]]
  curri2 = curri2/100
  currps = simple_RE_pvals[[nn]]
  r_currbeta = rank(currbetas)
  r_curri2 =  rank(sign(currbetas)*curri2)
  gsea_input_scores[[nn]] = cbind(currbetas,curri2,r_currbeta,r_curri2)
  if(grepl("muscle",nn)){
    inds = names(currbetas) %in% muscle_gtex_genes
  }
  else{
    inds = names(currbetas) %in% blood_gtex_genes
  }
  gsea_input_scores[[nn]] = gsea_input_scores[[nn]][inds,]
}
sapply(gsea_input_scores,dim)

# redo the lonterm muscle analysis without the overlap with the acute analysis
remove_by_gse<-function(gtable,gses){
  gtable = gtable[! (gtable$gse %in% gses),]
  return(gtable)
}
longterm_reduced_datasets = lapply(datasets[["longterm,muscle"]],
   remove_by_gse,gses = c("GSE59088","GSE27285","GSE28998","GSE28392",
        "GSE41769","GSE45426","GSE106865","GSE107934","GSE19062","GSE43219","GSE43856","GSE87749"))
simple_re_longterm_reduced = lapply(longterm_reduced_datasets,run_simple_re)
simple_re_longterm_reduced_beta = 
  sapply(simple_re_longterm_reduced,try_get_field,fname="beta")
simple_re_longterm_reduced_i2 = 
  sapply(simple_re_longterm_reduced,try_get_field,fname="I2")/100
gsea_input_scores[["longterm_overlap_removed"]] = cbind(
  simple_re_longterm_reduced_beta,
  simple_re_longterm_reduced_i2,
  rank(simple_re_longterm_reduced_beta),
  rank(sign(simple_re_longterm_reduced_beta)*simple_re_longterm_reduced_i2)
)
inds = names(simple_re_longterm_reduced_beta) %in% muscle_gtex_genes
gsea_input_scores[["longterm_overlap_removed"]] = gsea_input_scores[["longterm_overlap_removed"]][inds,]
sapply(gsea_input_scores,dim)

gsea_reactome_results = list()
for(nn in names(gsea_input_scores)){
  # beta GSEA
  currscores = gsea_input_scores[[nn]][,1]
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
  pathways = pathways[sapply(pathways,length)>20]
  pathways = pathways[sapply(pathways,length)<200]
  all_p_genes = unique(unlist(pathways))
  fgsea_res = fgsea(pathways,currscores[all_p_genes],nperm=50000)
  fgsea_res = fgsea_res[order(fgsea_res$padj),]
  gsea_reactome_results[[paste("I2",nn,sep="_")]] = fgsea_res
}

all_gsea_ps = unlist(sapply(gsea_reactome_results,function(x)x$pval))
fdr_BY_thr = max(all_gsea_ps[p.adjust(all_gsea_ps,method="BY") < 0.1],na.rm=T)
gsea_reactome_results_fdr = lapply(gsea_reactome_results,
                               function(x)x[x$pval<fdr_BY_thr,])
sapply(gsea_reactome_results_fdr,nrow)
lapply(gsea_reactome_results_fdr,function(x)x[1:3,c(1,3,5,7)])

save(gsea_reactome_results,gsea_reactome_results_fdr,gsea_input_scores,
     muscle_gtex_genes,blood_gtex_genes,
     file = paste(out_dir_rdata,"gsea_reactome_results.RData",sep=""))
save(muscle_gtex_genes,blood_gtex_genes,gsea_input_scores,
     file = paste(out_dir_rdata,"gsea_reactome_input.RData",sep=""))

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
write.table(all_fdr_corrected_gsea_results,file=paste(supp_path,"STable2.txt",sep="")
            ,row.names = F,quote=F,sep="\t")

# GSEA plots for Supp Figure 1
num_pathways = 6
library(gridExtra);library(grid)
for(nn in names(gsea_input_scores)){
  currscores = gsea_input_scores[[nn]][,1]
  currscores = na.omit(currscores)
  currscores = sample(currscores)
  curr_res = gsea_reactome_results[[paste("beta",nn,sep="_")]]
  inds1 = which(curr_res[["NES"]]>0)
  inds2 = which(curr_res[["NES"]]<0)
  inds = c(inds1[1:(num_pathways/2)],inds2[1:(num_pathways/2)])
  curr_selected_pathways = as.data.frame(curr_res[,1])[inds,1]
  pathways = reactomePathways(names(currscores))
  all_p_genes = unique(unlist(pathways))
  curr_f_name = paste0(out_dir_figs,"gsea_top6_",gsub(",","_",nn),".pdf")
  # pdf(curr_f_name)
  currpathways = pathways[curr_selected_pathways]
  curr_res = as.data.frame(curr_res)
  # dev.off()
  pl = plotGseaTable(currpathways,
                currscores[all_p_genes],curr_res,gseaParam=0.5,
                colwidths = c(5, 3, 0.8, 0, 1.2),render=T)
  dev.off()
  write.table(as.matrix(curr_res[inds,c(1,5,3)]),row.names = F,
              col.names = F,quote=F,sep="\t")
}

# venn diagrams for the enriched pathway overlap - for figure 2
library(gplots)llibrary(VennDiagram)
l = lapply(gsea_reactome_results_fdr[c(1,3,9,7)],function(x)x[[1]])
names(l) = c("acute,muscle","acute,blood","long-term,muscle","long-term,blood")
venn(l)
vp = venn.diagram(l,filename = NULL,fill = 2:5, alpha = 0.3,
                  cex = 1.4,cat.cex=1.6)
pdf(paste0(out_dir_figs,"gsea_venn.pdf"))
grid.draw(vp);dev.off()
dev.off()
grid.draw(vp);dev.off()
intersect_all = l[[1]]
for(ll in l){intersect_all = intersect(intersect_all,ll)}
write.table(t(t(intersect_all)),quote=F,row.names = F)
muscle_intersect = intersect(l$`acute,muscle`,l$`long-term,muscle`)
write.table(t(t(muscle_intersect)),quote=F,row.names = F)

# venn diagrams for the gene set overlap, for supp fig 6
l = lapply(analysis2selected_genes,names)
vp = venn.diagram(l,filename = NULL,fill = 2:5, alpha = 0.3,
                  cex = 1.4,cat.cex=1.6)
pdf(paste0(out_dir_figs,"genes_venn.pdf"))
grid.draw(vp);dev.off()
dev.off()
grid.draw(vp);dev.off()
intersect_all = l[[1]]
for(ll in l){intersect_all = intersect(intersect_all,ll)}
entrez2symbol[intersect(l$`acute,muscle`,l$`longterm,muscle`)]

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

# Create interaction networks for each analysis, 
# summarizing the number of detected genes.
get_gene_set_by_feature_name<-function(m,fname,up=T){
  cnames = colnames(m)
  m = m[sapply(m[,3],grepl,name1),]
  if(is.null(dim(m))){
    m = matrix(m,nrow=1)
    colnames(m) = cnames
  }
  effects = as.numeric(m[,"Effect"])
  if(up){return(m[effects>0,2])}
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

# Plot overlap matrices for Supp Fig 4
pdf(paste0(out_dir_figs,"Supp_Figure4AB.pdf"))
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

# Some enrichment analyses: commented out becuase this is slow
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
# Look at effects in untrained - excluded untrained others
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

# Create boxplots for Figure 3C
pdf(paste0(out_dir_figs,"Figure3C.pdf"))
par(mar=c(4,8,4,4))
names(l) = c("untrained","exercise","untrained","exercise",
             "untrained","exercise")
cols = c("red","red","cyan","cyan","green","green")
par(cex.lab=2.3,cex.axis=1.6)
boxplot(l,las=2,col=cols,horizontal=T,pch=20,ylim=c(0,2))
legend(x=0.7,y=4.8,c("muscle, acute","blood, acute","muscle, long-term"),
       fil=c("red","cyan","green"),cex=1.6)
dev.off()

# p-values, comparing the distributions in the plot
# acute, muscle (1.525566e-88)
wilcox.test(l[[1]],l[[2]],paired=T)$p.value
# acute_blood (6.222399e-08)
wilcox.test(l[[3]],l[[4]],paired=T)$p.value
# longterm muscle (3.110914e-68)
wilcox.test(l[[5]],l[[6]],paired=T)$p.value

# create figure src for the paper
l_for_src1 = cbind(l[[1]],l[[2]])
colnames(l_for_src1) = c("untrained","exercise")
write.table(l_for_src1,paste0(out_dir,"supp_tables/fig3c_src_acute_muscle.txt"),
                sep="\t",quote=F,row.names = F,col.names = T)
l_for_src2 = cbind(l[[3]],l[[4]])
colnames(l_for_src2) = c("untrained","exercise")
write.table(l_for_src2,paste0(out_dir,"supp_tables/fig3c_src_acute_blood.txt"),
            sep="\t",quote=F,row.names = F,col.names = T)
l_for_src3 = cbind(l[[5]],l[[6]])
colnames(l_for_src3) = c("untrained","exercise")
write.table(l_for_src3,paste0(out_dir,"supp_tables/fig3c_src_longterm_muscle.txt"),
            sep="\t",quote=F,row.names = F,col.names = T)

# 3. GO enrichments of the whole sets
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

############################################################################
# Examples of genes for the paper (appear in Figures 2 and 3 and Supp Fig 5)

# Our version for metafor's forest plot
gdata_forest_plot<-function(gdata,col.cex=0.8,plot.cex=0.8,title=""){
  ind1 = max(gdata$yi+3*gdata$sdd)
  ind2 = min(gdata$yi-3*gdata$sdd)
  annot_pos = ind2 - c(5,4:1)
  forest(x = gdata$yi,sei = gdata$sdd, xlim=c(min(annot_pos)-2, ind1), 
         slab = gdata$V1,
         ilab=cbind(gdata$N,gdata$training,round(gdata$avg_age,digits = 1), 
                    round(gdata$prop_males,digits = 2)*100,gdata$time),
         ilab.xpos=annot_pos, cex=plot.cex,
         ylim=c(-1, nrow(gdata)+3),
         xlab="Fold change", mlab="", 
         psize=1, header="Cohort",showweights = F,annotate=F,
         main = title)
  text(annot_pos, nrow(gdata)+2, c("N","Type","Age","%M","Time"),cex=col.cex)
  return(NULL)
}

#### Acute, muscle:
curr_genes = analysis2selected_genes$`acute,muscle`
#PGC1a
gene = "10891"
gene_name = entrez2symbol[[gene]]
gdata = meta_reg_datasets$`acute,muscle`[[gene]]
gdata = gdata[order(gdata$time),]
curr_times = rep("0-1h",nrow(gdata))
curr_times[gdata$time==2] = "2-5h"
curr_times[gdata$time==3] = ">20h"
gdata$time = curr_times
gdata$training = gsub("endurance","EE",gdata$training)
gdata$training = gsub("resistance","RE",gdata$training)
gdata$V1 = gsub("GE_","",gdata$V1)
# add src info for the paper
write.table(gdata,paste0(out_dir,"supp_tables/pcg1a_gdata_source.txt"),
            sep="\t",quote=F,row.names = T,col.names = T)

analysis1 = all_meta_analysis_res$`acute,muscle`
analysis1[[gene]][[1]]$mod_p
aic_diff = analysis1[[gene]][[1]]$aic_c - analysis1[[gene]]$`simple:base_model`$aic_c
png(paste0(out_dir_figs,"Figure2A.png"), units="px", width=1600, height=1600, res=300)
gdata_forest_plot(gdata,col.cex = 0.85)
dev.off()

# Other genes: LPL, CPT1B, SMAD3, ACTN3, VEGFA, FOXO1 and IL6R
genes = c("1375","4023","4088","89","7422","2308","3570")
# Validated genes
library("metaviz")
validated_genes = c(
  "SMAD3" = "4088","ID1" = "3397","NR4A1" = "3164",
  "HES1" = "3280","SCN2B" = "6327","SLC25A25" = "114789",
  "PPARGC1A" = "10891","MRPL34" = "64981","SH3KBP1" = "30011"
)
for(gene in validated_genes){
  gene_name = entrez2symbol[[gene]]
  gdata = meta_reg_datasets$`acute,muscle`[[gene]]
  gdata = gdata[order(gdata$time,decreasing = F),]
  curr_times = rep("0-1h",nrow(gdata))
  curr_times[gdata$time==2] = "2-5h"
  curr_times[gdata$time==3] = ">20h"
  gdata$time = curr_times
  gdata$training = gsub("endurance","EE",gdata$training)
  gdata$training = gsub("resistance","RE",gdata$training)
  gdata$V1 = gsub("GE_","",gdata$V1)
  slabels = paste(gdata$training,curr_times,sep=",")
  analysis1 = all_meta_analysis_res$`acute,muscle`
  analysis1[[gene]][[1]]$mod_p
  png(paste0(out_dir_figs,gene_name,".png",sep=""), 
      units="px", width=1600, height=1600, res=300)
  gdata_forest_plot(gdata,col.cex = 0.85,plot.cex = 0.85,title=gene_name)
  dev.off()
  write.table(gdata,paste0(out_dir,"supp_tables/",gene_name,"_gdata_source.txt"),
              sep="\t",quote=F,row.names = T,col.names = T)
}

# COL4A1 in longterm muscle
gene = "1282"
gene_name = entrez2symbol[[gene]]
all_meta_analysis_res$`longterm,muscle`[[gene]]
gdata = meta_reg_datasets$`longterm,muscle`[[gene]]
gdata = gdata[order(gdata$time),]
curr_times = rep("< 150 d",nrow(gdata))
curr_times[gdata$time==2] = ",> 150 d"
gdata$time = curr_times
gdata$V1 = gsub("GE_","",gdata$V1)
gdata$training = gsub("endurance","EE",gdata$training)
gdata$training = gsub("resistance","RE",gdata$training)
slabels = paste(gdata$training,curr_times,sep="")
# slabels = paste(gdata$V1,slabels,sep=",")
analysis1 = all_meta_analysis_res$`acute,muscle`
png(paste0(out_dir_figs,"Figure3C.png"), units="px", width=1600, height=1600, res=300)
gdata_forest_plot(gdata,col.cex = 0.85)
dev.off()
write.table(gdata,paste0(out_dir,"supp_tables/",gene_name,"_gdata_source.txt"),
            sep="\t",quote=F,row.names = T,col.names = T)

############################################################################
# Interpretation of the results: defining subgroups by clustering mean patterns

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

# Run clustering analysis on the selected genes
gene_subgroups = list();gene_t_patterns = list()
set.seed(123) # for reproducibility
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

# some stats: acute and age genes
unique(unlist(gene_subgroups[
  grepl("age",names(gene_subgroups)) & grepl("acute",names(gene_subgroups))]))
unique(unlist(gene_subgroups[
  grepl("age",names(gene_subgroups)) & grepl("long",names(gene_subgroups))]))
unique(unlist(gene_subgroups[
  grepl("males",names(gene_subgroups)) & grepl("long",names(gene_subgroups))]))

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

# Some stats about the clustering solution
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
    colnames(mat) = gsub("^1","<150d",colnames(mat))
    colnames(mat) = gsub("^2",">150d",colnames(mat))
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
  
  # some specific adjustments
  if(set_name == "longterm,muscle,time;avg_age,1"){
    colnames(mat) = sapply(colnames(mat),function(x)strsplit(x,split=", ")[[1]][1])
    heatmap.2(mat,trace = "none",scale = "none",Colv = T,col=bluered,
              cexRow = 1.2,main=curr_main,cexCol = 1.15,
              Rowv = T,srtCol=60,hclustfun = hclust_func,density.info="none",
              key = F,margins = c(10,10))
    dev.off()
  }
  if(set_name == "acute,muscle,time;avg_age,1"){
    colnames(mat) = gsub(" ","",colnames(mat))
    colnames(mat) = gsub(",,",", ",colnames(mat))
    heatmap.2(mat,trace = "none",scale = "none",Colv = T,col=bluered,
              cexRow = 1.2,main="",cexCol = 1.15,
              Rowv = T,srtCol=60,hclustfun = hclust_func,density.info="none",
              key = F,margins = c(10,10))
    dev.off()
  }
  if(set_name == "longterm,muscle,prop_males,1"){
    colnames(mat) = as.character(round(100*as.numeric(colnames(mat))))
    heatmap.2(mat,trace = "none",scale = "none",Colv = T,col=bluered,
              cexRow = 0.2,main="",cexCol = 1.25,
              Rowv = T,srtCol=60,hclustfun = hclust_func,density.info="none",
              key = F,margins = c(10,10))
    dev.off()
  }
  else{
    heatmap.2(mat,trace = "none",scale = "none",Colv = T,col=bluered,
              cexRow = cex_genes,main=curr_main,
              Rowv = T,srtCol=45,hclustfun = hclust_func,density.info="none",
              key.title = NA,keysize = 1.1,key.xlab = "t-statistic",
              key.par = list("cex.axis"=1.1),margins = c(10,10))
  }
  
  write.table(mat,paste0(out_dir,"supp_tables/",gsub(",","_",set_name),"_t_table.txt"),
    sep="\t",row.names=T,col.names=T,quote=F)
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

# get the top enrichments of the set
set_name = "longterm,muscle,base_model,1"
set_names = c("longterm,muscle,base_model,1",
              "longterm,muscle,base_model,2")
table_name = "longterm,muscle"
set_genes = unlist(gene_subgroups[set_names])
mat = mean_effect_matrices[[table_name]][set_genes,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,apply(mat,2,sd)>0]
mat = mat[,apply(mat==0,2,sum)/nrow(mat) < 0.5]
mat = mat[apply(mat==0,1,sum)/ncol(mat) < 0.5,]
mat[mat>4]=4;mat[mat< -4]=-4
colnames(mat) = NULL;
pdf(paste(out_dir_figs,"Figure5A.pdf"))
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,cexRow = 1.2,
          add.expr={makeRects(mat==0)},Rowv = F,margins = c(5,8))
dev.off()
write.table(mat,paste0(out_dir,"supp_tables/",gsub(",","_",set_name),"_t_table.txt"),
            sep="\t",row.names=T,col.names=T,quote=F)
# get the top enrichments of the set
curr_gos = gene_subgroup_enrichments_fdr[gene_subgroup_enrichments_fdr[,1]==set_name,]
curr_pathways = reactome_pathways_subgroups_fdr[reactome_pathways_subgroups_fdr[,1]==set_name,]

set_name = "acute,blood,base_model,2"
table_name = "acute,blood"
set_genes = unlist(gene_subgroups[set_name])
mat = mean_effect_matrices[[table_name]][set_genes,]
rownames(mat) = unlist(entrez2symbol[rownames(mat)])
mat = mat[,apply(mat,2,sd)>0]
mat = mat[,apply(mat==0,2,sum)/nrow(mat) < 0.5]
mat = mat[apply(mat==0,1,sum)/ncol(mat) < 0.5,]
mat[mat>4]=4;mat[mat< -4]=-4
matnames = sapply(strsplit(colnames(mat),split=";"),function(x)x)[-2,]
matnames[1,] = gsub("A_","ID:",matnames[1,])
matnames[2,] = paste0(matnames[2,],"hr")
matnames[3,] = paste0(matnames[3,],"y")
matnames[4,] = paste0(as.numeric(matnames[4,])*100,"% M")
matnames = apply(matnames[1:2,],2,paste,collapse=", ")
colnames(mat) = matnames
pdf(paste(out_dir_figs,"Figure3D.pdf"))
heatmap.2(mat,trace = "none",scale = "none",Colv = F,col=bluered,
          cexRow = cex_genes,main="",
          Rowv = F,srtCol=25,hclustfun = hclust_func,density.info="none",
          key.title = NA,keysize = 1.2,key.xlab = "t-statistic",
          add.expr={makeRects(mat==0)},
          key.par = list("cex.axis"=1),margins = c(10,10))
dev.off()
write.table(mat,paste0(out_dir,"supp_tables/",gsub(",","_",set_name),"_t_table.txt"),
            sep="\t",row.names=T,col.names=T,quote=F)

lonterm_muscle_base_model_effects = rbind(
  lonterm_muscle_base_model_effects, 
  cbind(set_genes,
    unlist(entrez2symbol[set_genes]),
    sapply(all_meta_analysis_res$`longterm,muscle`[set_genes],
           function(x)x[[1]]$coeffs[1,1]))
)
colnames(lonterm_muscle_base_model_effects) = c("entrez","Symbol","fchange")
write.table(lonterm_muscle_base_model_effects,
            file=paste(out_dir_figs,"Figure4_longterm_genes_nodes.txt"),
            quote = F,row.names = F,col.names = T,sep="\t")

###################################################################
# The different patterns of acute muscle, time-dependent response
# This analysis is important for all that is presented in Figure 4

pref = "acute,muscle,time,\\d"
clusters = names(gene_subgroups)[grepl(pref,names(gene_subgroups))]
cols = c("blue","black","red","green")
names(cols) = clusters

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
  print(paste("stats for cluster:",set_name))
  write.table(t(c(colMeans(mat),apply(mat,2,sd))),sep="\t",row.names = F,col.names = F)
}
dev.off()

# Write input to DREM - 
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
  mm = rbind(mm,cbind(m[[cluster_name]],rep(cluster_name,length(m[[cluster_name]]))))
}
write.table(mm,file=paste(out_dir_figs,"acute_muscle_time_subgroups.txt",sep=""),
            sep="\t",col.names = F,row.names = F,quote=F)

###################################################################################
# Look at the covariate corrleations

pdf(paste(out_dir_figs,"Supp_Figure3.pdf"))
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

###################################################################################
# Prepare supplementary tables 
# reload the raw input data for the meta-analyses
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
  colnames(m)[colnames(m)=="R2"] = "R2(if applicable)"
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
write.table(supp_table_genes,file=paste(supp_path,"STable3_v2.txt",sep="")
            ,row.names = F,quote=F,sep="\t")
hist(as.numeric(supp_table_genes[,8]))
table(is.na(as.numeric(supp_table_genes[,8])))
sum(as.numeric(supp_table_genes[,8])>60,na.rm=T)/sum(!is.na(as.numeric(supp_table_genes[,8])))

sheet_counter=4
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
sheet_counter = sheet_counter+1
supp_table_enrichments = gene_subgroup_enrichments_fdr[,c(1:4,9:10,8)]
colnames(supp_table_enrichments)[6] = "q-value"
colnames(supp_table_enrichments)[5] = "Genes"
colnames(supp_table_enrichments)[1] = "Discovered in"
# write.xlsx(supp_table_enrichments,file=supp_file,
#            sheetName = paste("STable",sheet_counter,"_GO_enrichments",sep=""),
#            row.names = F,append=T)
write.table(supp_table_enrichments,file=paste(supp_path,"STable8.txt",sep="")
            ,row.names = F,quote=F,sep="\t")

# Same for Reactome
sheet_counter = sheet_counter+1
supp_table_enrichments = reactome_pathways_subgroups_fdr[,c(1,3,8,9,6)]
colnames(supp_table_enrichments)[3] = "q-value"
colnames(supp_table_enrichments)[1] = "Discovered in"
# write.xlsx(supp_table_enrichments,file=supp_file,
#            sheetName = paste("STable",sheet_counter,"_Reactome_enrichments",sep=""),
#            row.names = F,append=T)
write.table(supp_table_enrichments,file=paste(supp_path,"STable9.txt",sep="")
            ,row.names = F,quote=F,sep="\t")

# save the workspace
save.image(file=paste(out_dir,"meta_analysis_interpretation_results.RData",sep=""))

###################################################################################
# Prepare gene-csv files for the data portal
# load the entire session with all results from the analyses above
load("meta_analysis_interpretation_results.RData")
library(data.table)
out_dir_gcp = paste0(out_dir,"gcp/")
dir.create(out_dir_gcp,showWarnings = F)
for(nn in names(meta_reg_datasets)){
  nn_out_dir_gcp = paste0(out_dir_gcp,gsub(",","_",nn),"/")
  dir.create(nn_out_dir_gcp,showWarnings = F)
  genes = names(meta_reg_datasets[[nn]])
  curr_gene_table = c()
  for(gene in genes){
    gname = entrez2symbol[[gene]]
    gdata = meta_reg_datasets[[nn]][[gene]]
    base_model = simple_REs[[nn]][[gene]]
    if(length(base_model)==1 && is.na(base_model)){next}
    i2 = base_model$I2
    tau2 = base_model$tau2
    p_model = base_model$pval
    selected = is.element(gene,names(analysis2selected_genes[[nn]]))

    if(grepl("acute",nn)){
      curr_times = rep("0-1h",nrow(gdata))
      curr_times[gdata$time==2] = "2-5h"
      curr_times[gdata$time==3] = ">20h"
    }
    else{
      curr_times = rep("0-150 days",nrow(gdata))
      curr_times[gdata$time==2] = ">150 days"
    }

    gdata$time = curr_times
    meta_reg_analysis = all_meta_analysis_res[[nn]][[gene]]
    aic_c = meta_reg_analysis[[1]]$aic_c
    aic_diff = meta_reg_analysis[[1]]$aic_c - meta_reg_analysis$`simple:base_model`$aic_c
    i2_selected_model = NA
    if("I2" %in% names(meta_reg_analysis[[1]])){
      i2_selected_model = meta_reg_analysis[[1]]$I2
    }

    selected_model_name = names(meta_reg_analysis)[1]
    if(grepl("base_model",selected_model_name)){
      selected_model_name = ":base_model (simple RE)" # add ':' just for the split later
    }
    else{
      p_model = meta_reg_analysis[[1]]$mod_p
    }

    selected_model_name = strsplit(selected_model_name,":")[[1]][2]
    selected_model_name = gsub("avg_","",selected_model_name)
    selected_model_name = gsub(";",",",selected_model_name)
    curr_gene_table = rbind(curr_gene_table,
      c(gname,gene,i2,tau2,selected_model_name,
        p_model,i2_selected_model,aic_c,aic_diff,
        selected))

    curr_gfile = paste0(nn_out_dir_gcp,gname,".csv")
    fwrite(gdata,file=curr_gfile,sep=",",quote = F)
  }
  colnames(curr_gene_table) = c("Symbol","Entrez ID",
                                "I2(base model)","Tau2(base model)",
                                "Selected model name","P-value",
                                "I2(selected model - if available)",
                                "AICc(selected model)",
                                "AICcDiff","Selected?")
  fname = paste0(out_dir_gcp,gsub(",","_",nn),"_gene_stats.csv")
  fwrite(curr_gene_table,file=fname,sep=",",quote=F)
}

###################################################################################
# For the paper revision: examine blood vs. muscle using the base models
x1 = simple_RE_beta$`acute,muscle`
x2 = simple_RE_beta$`acute,blood`
x12_shared_genes = intersect(names(x1),names(x2))
x12_shared_genes = x12_shared_genes[!is.na(x1[x12_shared_genes])]
x12_shared_genes = x12_shared_genes[!is.na(x2[x12_shared_genes])]
plot(x = x1[x12_shared_genes],y = x2[x12_shared_genes],pch=20,cex=0.7,
     xlab = "Muscle",ylab="Blood",main="log2 fold change",col="gray")
cor.test(x1[x12_shared_genes],x2[x12_shared_genes])
cor.test(x1[x12_shared_genes],x2[x12_shared_genes],method="spearman")
median(abs(x1),na.rm=T);median(abs(x2),na.rm=T)
table(abs(x1)>0.1);table(abs(x2)>0.1)

x1 = simple_RE_I2s$`acute,muscle`
x2 = simple_RE_I2s$`acute,blood`
x12_shared_genes = intersect(names(x1),names(x2))
x12_shared_genes = x12_shared_genes[!is.na(x1[x12_shared_genes])]
x12_shared_genes = x12_shared_genes[!is.na(x2[x12_shared_genes])]
plot(x = x1[x12_shared_genes],y = x2[x12_shared_genes],pch=20,cex=0.7,
     xlab = "Muscle",ylab="Blood",main="I2",col="gray")
cor.test(x1[x12_shared_genes],x2[x12_shared_genes])
cor.test(x1[x12_shared_genes],x2[x12_shared_genes],method="spearman")
cor.test(x1[x12_shared_genes],x2[x12_shared_genes],method="spearman")$p.value
median(x1,na.rm=T);median(x2,na.rm=T)

x1 = simple_RE_beta$`acute,muscle`
x2 = simple_RE_beta$`longterm,muscle`
x12_shared_genes = intersect(names(x1),names(x2))
x12_shared_genes = x12_shared_genes[!is.na(x1[x12_shared_genes])]
x12_shared_genes = x12_shared_genes[!is.na(x2[x12_shared_genes])]
plot(x1[x12_shared_genes],x2[x12_shared_genes],pch=20,cex=0.7)
cor.test(x1[x12_shared_genes],x2[x12_shared_genes])
cor.test(x1[x12_shared_genes],x2[x12_shared_genes],method="spearman")

x1 = simple_RE_beta$`longterm,blood`
x2 = simple_RE_beta$`longterm,muscle`
x12_shared_genes = intersect(names(x1),names(x2))
x12_shared_genes = x12_shared_genes[!is.na(x1[x12_shared_genes])]
x12_shared_genes = x12_shared_genes[!is.na(x2[x12_shared_genes])]
plot(x1[x12_shared_genes],x2[x12_shared_genes],pch=20,cex=0.7)
cor.test(x1[x12_shared_genes],x2[x12_shared_genes])
cor.test(x1[x12_shared_genes],x2[x12_shared_genes],method="spearman")

x1 = simple_RE_beta$`longterm,blood`
x2 = simple_RE_beta$`acute,blood`
x12_shared_genes = intersect(names(x1),names(x2))
x12_shared_genes = x12_shared_genes[!is.na(x1[x12_shared_genes])]
x12_shared_genes = x12_shared_genes[!is.na(x2[x12_shared_genes])]
plot(x1[x12_shared_genes],x2[x12_shared_genes],pch=20,cex=0.7)
cor.test(x1[x12_shared_genes],x2[x12_shared_genes])
cor.test(x1[x12_shared_genes],x2[x12_shared_genes],method="spearman")

# examine aicc diffs
all_aicc_diffs = lapply(all_meta_analysis_res,function(x)sapply(x,get_aicc_diff))
sapply(all_aicc_diffs,length)



