######################################################
# Load libraries
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
  if(!is.null(res)){return(res)}
  
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

gdata_forest_plot<-function(gdata,model,fulltable=T,col.cex=0.8,plot.cex=0.8,title=""){
  gdata$V1 = gsub("GE_","",gdata$V1)
  gdata$training = gsub("endurance","endur",gdata$training)
  gdata$training = gsub("resistance","resist",gdata$training)
  ind1 = max(gdata$yi+3*gdata$sdd)
  ind2 = min(gdata$yi-3*gdata$sdd)
  
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
  
  forest(x = gdata$yi,sei = gdata$sdd, xlim=c(min(annot_pos)-2, ind1), 
         slab = slab,
         ilab = ilab,
         ilab.xpos = annot_pos, 
         cex=plot.cex,
         ylim=c(-2, nrow(gdata)+3),
         xlab="Fold change", 
         psize=1, header=slab_name,
         showweights = F,annotate=F,
         fonts = "Helvetica",
         col = "darkblue", pch = 20, steps = 10,
         main = title
  )
  if(fulltable){
    text(annot_pos, nrow(gdata)+2, c("N","Type","Age","%M","Time"),cex=col.cex)
  }
  ### add summary estimate from the random-effects model to forest plot
  addpoly(model,col="gray")
  return(NULL)
}

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
me_result = meta_analysis_wrapper(d,random= ~V1|gse, struct="CS")
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

png(paste0(getwd(),"/",opt$out,".forest.png"),
    width = 4, height = 4, units = 'in', res = 300,pointsize=10)
if(opt$forest == 1){
  tmp = gdata_forest_plot(d,me_result,F)
}
if(opt$forest == 2){
  tmp = gdata_forest_plot(d,me_result)
}
tmp2 = dev.off()

print("################################")
print("Mixed effects model summary stats:")
write.table(t(t(results_summary_stats)),col.names=F,quote=F,sep="\t")



