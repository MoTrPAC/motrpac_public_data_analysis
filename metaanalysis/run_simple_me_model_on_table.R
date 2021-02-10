######################################################
# Load libraries
.libPaths("~/R/packages") # add our local installed R packages
libraries = c("metafor","optparse")
for(lib in libraries){
  load_worked = F
  try({
    suppressPackageStartupMessages(library(lib,character.only = T))
    load_worked = T
  })
  if(!load_worked){
    try({
      install.packages(lib)
      suppressPackageStartupMessages(library(lib,character.only = T))
      load_worked = T
    })
  }
  if(!load_worked){
    print(paste("Cannot run the script without the following packages:",paste(libraries,collapse=",")))
    q("no",status=1)
  }
}

######################################################
# Some helper functions
meta_analysis_wrapper<-function(gdata,func = rma.mv,...){
  res=NULL
  try({
    res = func(yi,vi,data=gdata,...)
  },silent=TRUE)
  if(!is.null(res)){
    return(res)
  }
  
  for(rel.tol in c(1e-8,1e-7,1e-6)){
    cc=list(iter.max=10000,rel.tol=rel.tol)
    res = NULL
    try({
      res = func(yi,vi,data=gdata,control=cc,...)
    },silent=TRUE)
    if(!is.null(res)){return(res)}
  }
  return(NA)
}


gdata_forest_plot<-function(gdata,model,fulltable=T,col.cex=0.8,
                            plot.cex=0.8,min_interval_size = 1,
                            max_interval_size = 10){
  gdata$V1 = gsub("GE_","",gdata$V1)
  gdata$training = gsub("endurance","endur",gdata$training)
  gdata$training = gsub("resistance","resist",gdata$training)
  gdata$training = gsub("_treatment","",gdata$training)
  ind1 = max(gdata$yi+3*gdata$sdd)
  ind2 = min(gdata$yi-3*gdata$sdd)
  #print(paste(ind1,ind2))
  if(ind1 - ind2 < min_interval_size){
    diff = ind1 - ind2 - min_interval_size
    ind1 = ind1 + diff/2
    ind2 = ind2 - diff/2
  }
  if(ind1 - ind2 > max_interval_size){
    diff = ind1 - ind2 - max_interval_size
    ind1 = ind1 - (diff/2)
    ind2 = ind2 + (diff/2)
    #print(paste(ind1,ind2))
  }
  
  I2 = format(model$I2,digits=4)
  p = format(model$pval,digits=2)
  beta = format(model$b[1,1],digits=2)
  re_bound = paste0("[",
       format(model$ci.lb,digits=3),
       ",",
       format(model$ci.ub,digits=3),
       "]")
  main = paste0("RE: P-value=",p,", beta=",beta,", I2=",I2)
  main_col = "black"
  if(I2 > 50){
    main_col = "red"
    main = paste0(main,"*")
  }
  if(I2 > 75){
    main_col = "darkred"
    main = paste0(main,"*")
  }
  
  annot_pos = ind2
  ilab = NULL
  slab = 1:nrow(gdata)
  slab_name = "Table row"
  if(fulltable){
    annot_pos = ind2 - c(5,4:1)
    ilab = cbind(gdata$N,gdata$training,round(gdata$avg_age,digits = 1), 
                 round(gdata$prop_males,digits = 2)*100,gdata$time)
    slab = gdata$V1
    slab_name = "Cohort"
  }
  plot_xlim = c(min(annot_pos)-2, ind1)
  #print(plot_xlim)
  
  forest(x = gdata$yi,sei = gdata$sdd, 
         xlim = plot_xlim,
         alim = c(ind2,ind1),
         slab = slab,
         ilab = ilab,
         ilab.xpos = annot_pos, 
         header=slab_name,
         cex=plot.cex,
         ylim=c(-2, nrow(gdata)+3),
         xlab="log2 fold change", 
         psize = 1,
         top = 3,
         digits = 1,
         showweights = F,annotate=F,
         fonts = "Helvetica",
         col = "darkblue",
         pch = 20, steps = 10,
         main = main,col.main=main_col,cex.main=1
  )
  if(fulltable){
    text(annot_pos, nrow(gdata)+2, c("N","Type","Age","%M","Time"),cex=col.cex)
  }
  ### add summary estimate from the random-effects model to forest plot
  poly_col = "gray"
  if(I2 > 50){
    poly_col = "red"
  }
  if(I2 > 75){
    poly_col = "darkred"
  }
  addpoly(model,annotate = F,col=poly_col)
  
  return(NULL)
}

# # Genes to test:
# gene = "5166" #PDK4
# gene= "10891" # PCG1a
# gname = entrez2symbol[[gene]]
# nn = "acute,blood"
# nn = "longterm,muscle"
# nn = "acute,muscle"
# gdata = meta_reg_datasets[[nn]][[gene]]
# nrow(gdata)
# 
# model0 = meta_analysis_wrapper(gdata,func = rma)
# model = meta_analysis_wrapper(gdata,random= ~V1|gse, struct="CS")
# if(is.null(model)){model=model0}
# model$I2 = model0$I2
# 
# png("~/Desktop/test_forest.png",
#     width = 4.5, height = 4, units = 'in', res = 400, pointsize=8)
# par(mar=c(5,3,3,2))
# gdata_forest_plot(gdata,model)
# dev.off()


######################################################
# Parse the input parameters
option_list <- list(
  make_option("--csv", default="COL4A.csv",
              help="Input csv file for meta-analysis"),
  make_option("--forest",default=1,
              help="Numeric, which forest plot should be created:\n0 = none\n1 = the metafor default\n2 = our version with the extended table"),
  make_option("--verbose",default=T,
              help="TRUE: should the script print additional output other than the ME results?"),
  make_option("--out",default="out",
              help="Name of the output png file for the forest plot, the script will add a suffix:forest.png")
)
opt <- parse_args(OptionParser(option_list=option_list))
if(opt$verbose){
  print("Running simple mixed effects meta-analysis, input parameters:")
  for(param in names(opt)){
    print(paste("Parameter:",param," - ", opt[[param]]))
  }
}

######################################################
# Check the input
if(is.na(opt$csv) || !file.exists(opt$csv)){
  print("Input csv file does not exist")
  q("no",status=1)
}

d = read.csv(opt$csv,stringsAsFactors = F)
if(opt$verbose){
  print("completed reading the input csv file, head:")
  print(head(d))
}

if(any(!(c("V1","gse","yi","vi") %in% colnames(d)))){
  print("One of the following columns is not in d: V1, gse, yi, or vi")
  q("no",status=1)
}

opt$forest = as.numeric(opt$forest)
if(is.na(opt$forest) || !(opt$forest %in% 0:2)){
  print("Wrong specification of the forest parameters, rerun the script with -h or --help")
  q("no",status=1)
}

######################################################
# run the analysis
me_result0 = meta_analysis_wrapper(d,func = rma)
me_result = meta_analysis_wrapper(d,random= ~V1|gse, struct="CS")
if(is.null(me_result)){me_result=me_result0}
me_result$I2 = me_result0$I2
fitscores = fitstats(me_result)
if(opt$verbose){
  print("Completed fitting the mixed effects model:")
  print(summary(me_result))
}

results_summary_stats = c(
  "beta" = me_result$beta,
  "beta_se" = me_result$se,
  "CI lower bound:" = me_result$ci.lb,
  "CI upper bound:" = me_result$ci.ub,
  "model p-value" = me_result$pval,
  "sigma2" = me_result$sigma2,
  "tau2" = me_result$tau2,
  "I2" = me_result$I2,
  "QE" = me_result$QE,
  "Qp (heterogeneity p-value)" = me_result$QEp,
  "AICc" = fitscores[5,1]
)

######################################################
# Create the forest plot
if(opt$forest == 0){q("no")}
if(opt$verbose){
  print("creating a forest plot, input parameter is:")
  print(opt$forest)
  print("output file is:")
  paste0(getwd(),"/",opt$out,".forest.png")
}

plot_pointsize = 10
if(nrow(d)>20){
  plot_pointsize = 8
}

png(paste0(getwd(),"/",opt$out,".forest.png"),
    width = 4.5, height = 4, units = 'in', 
    res = 400, pointsize = plot_pointsize)
par(mar=c(5,3,3,2))
if(opt$forest == 1){
  tmp = gdata_forest_plot(d,me_result,F)
}
if(opt$forest == 2){
  tmp = gdata_forest_plot(d,me_result)
}
tmp2 = dev.off()

# print("################################")
# print("warnings:")
# print(warnings())
print("################################")
print("Mixed effects model summary stats:")
write.table(t(t(results_summary_stats)),col.names=F,quote=F,sep="\t")



