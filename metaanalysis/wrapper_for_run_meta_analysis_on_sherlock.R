# This script is a wrapper to run run_meta_analysis_on_sherlock.R
# We run each metanalysis separately and break each gene list to a group.
# The goal is to create MANY small jobs that run in parallel.
# The assumptions:
#   we have meta_analysis_input.RData in the working dir
#   we have run_meta_analysis_on_sherlock.R in the working dir
# TODO: make the input above command line arguments

# Helper functions for running in sherlock
get_sh_default_prefix<-function(err="",log="",time="6:00:00"){
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=",time,sep=""),
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH -c 1",
      "#SBATCH --mem=2000",
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "#SBATCH -x sh-113-15",
      "",
      "module load R"
    )
  )
}
run_command<-function(cmd,out_path,name,batch_script_func=get_sh_default_prefix,...){
  err_path = paste(out_path,name,".err",sep="")
  log_path = paste(out_path,name,".log",sep="")
  curr_sh_file = paste(out_path,name,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,cmd)
  system(paste("sbatch",curr_sh_file))
}
print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}
get_my_jobs<-function(){
  tmp = paste("tmp",abs(rnorm(1)),sep='')
  system(paste("sacct > ",tmp),wait = T)
  jobs = readLines(tmp)[-c(1:2)]
  jobs = jobs[!grepl("^\\d+\\.",jobs)]
  jobs = t(sapply(jobs,function(x)strsplit(x,split="\\s+")[[1]]))
  rownames(jobs) = as.character(jobs[,1])
  jobs = jobs[,-1]
  jobs = jobs[jobs[,1]!="bash",]
  system(paste("rm",tmp))
  new_jobs = rownames(jobs)[jobs[,5]=="RUNNING" | jobs[,5]=="PENDING"]
  if(length(new_jobs)==1){
    return(matrix(jobs[new_jobs,],nrow=1))
  }
  return(jobs[new_jobs,])
}
get_job_id<-function(x){return(x[1])}
wait_for_job<-function(jobs_before=NULL,waittime=30,max_wait=6000){
  Sys.sleep(waittime)
  curr_jobs = get_my_jobs()
  if(length(curr_jobs)==0){return(NULL)}
  new_jobs = rownames(curr_jobs)[curr_jobs[,5]=="RUNNING" | curr_jobs[,5]=="PENDING"]
  if(length(new_jobs)==0){return(matrix(NA,0,0))}
  print(paste("new added jobs are: ",new_jobs))
  i = 1
  while(T){
    Sys.sleep(waittime)
    curr_jobs = get_my_jobs()
    if(length(curr_jobs)==0){return(NULL)}
    new_jobs = rownames(curr_jobs)[curr_jobs[,5]=="RUNNING" | curr_jobs[,5]=="PENDING"]
    if(length(new_jobs)==0){return(matrix(NA,0,0))}
    i=i+1
    if(i>max_wait){break}
  }
}


# Load the input
load("meta_analysis_input.RData")
num_genes_per_job = 200
for (i in 1:length(meta_reg_datasets)){
  curr_length = length(meta_reg_datasets[[i]])
  print(paste("############ analyzing dataset number:",i))
  start = 1
  while(start < curr_length){
    end = start+num_genes_per_job-1
    
    cmd = paste(
      "Rscript run_meta_analysis_on_sherlock.R", 
      "~/motrpac_metaanalysis/v2_march_2019/meta_analysis/",
      1,i,start,end
    )
    run_command(cmd,paste(getwd(),"/",sep=""),
                paste("run_",i,"_",start,"_",end,sep=""))
    
    start = end+1
    # if(end>=100){break} # for QA
  }
  # break # for QA
}

wait_for_job()
# read the output and put all results in one object

all_meta_analysis_res = list()
for (i in 1:length(meta_reg_datasets)){
  curr_length = length(meta_reg_datasets[[i]])
  nn = names(meta_reg_datasets)[i]
  if(!is.element(nn,set=names(all_meta_analysis_res))){
    all_meta_analysis_res[[nn]] = list()
  }
  newnn = gsub(",","_",nn)
  start = 1
  while(start < curr_length){
    end = start+num_genes_per_job-1
    end = min(end,curr_length)
    curr_file = paste("meta_analysis_results",newnn,start,end,".RData",sep="_")
    curr_rdata = curr_file
    analysis_res = NULL
    try({load(curr_rdata)})
    if(is.null(analysis_res)){
      print("#### failed loading data file:")
      print(curr_file)
      start = end+1
      next
    }
    for(gg in names(analysis_res)){
      all_meta_analysis_res[[nn]][[gg]] = analysis_res[[gg]]
    }
    system(paste("rm",curr_file))
    start = end+1
  }
}
print("merged all rdata files, results size is:")
sapply(all_meta_analysis_res,length)
save(all_meta_analysis_res,file="meta_analysis_results.RData")

# clear the working directory
system("rm *.log")
system("rm *.err")
for(i in 1:length(all_meta_analysis_res)){
  system(paste("rm run_",i,"_*",sep=""))
}



