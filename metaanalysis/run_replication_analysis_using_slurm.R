
# Input dir: should have some p-value matrices.
# These are files of the format *pvals.txt
curr_dir = "/home/users/davidama/motrpac_metaanalysis/rep_analysis/"
setwd(curr_dir)
pvals_files = list.files()
pvals_files = pvals_files[grepl("pvals.txt",pvals_files)]

get_sh_prefix_one_node_specify_cpu_and_mem<-function(err="",log="",Ncpu=4,mem_size=32000,time="24:00:00"){
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=",time,sep=""),
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH --nodes=1",
      paste("#SBATCH -c",Ncpu),
      paste("#SBATCH --mem=",mem_size,sep=""),
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "#SBATCH -x sh-113-15",
      "",
      "module load R"
    )
  )
}
print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}
run_sbatch_rscript_command<-function(cmd,out_path,name,batch_script_func=get_sh_default_prefix,...){
  err_path = paste(out_path,name,".err",sep="")
  log_path = paste(out_path,name,".log",sep="")
  curr_sh_file = paste(out_path,name,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,cmd)
  system(paste("sbatch",curr_sh_file))
}

for (ff in pvals_files){
  for(m in c("bum","znormix")){
    currname = gsub("_pvals.txt","",ff)
    currname = paste(currname,"_",m,sep="")
    outfile = paste(curr_dir,currname,".txt",sep="")
    curr_cmd = paste("Rscript /home/users/davidama/repos/screen/R/cmd_line_runnable.R",
                     ff,"1",m,outfile,"path=/home/users/davidama/R/packages/ emEps=1e-5 nH=10000 minP=0.001")
    if(m=="bum"){
      curr_cmd = paste("Rscript /home/users/davidama/repos/screen/R/cmd_line_runnable.R ",
                       ff,"1",m,outfile,"path=/home/users/davidama/R/packages/ emEps=1e-3 nH=10000 minP=0.001")
    }
    run_sbatch_rscript_command(curr_cmd,curr_dir,currname,get_sh_prefix_one_node_specify_cpu_and_mem)
  }
}

# Go over the results, save them into an RData file and do some comparisons
screen_results = list()
for (ff in pvals_files){
  currname = gsub("_pvals.txt","",ff)
  screen_results[[currname]] = list()
  screen_results[[currname]][["pvals"]] = read.delim(ff,row.names = 1,header=T)
  for(m in c("bum","znormix")){
    outfile = paste(curr_dir,currname,"_",m,".txt",sep="")
    screen_results[[currname]][[m]] = read.delim(outfile,row.names = 1,header=T)
  }
}

save(screen_results,file="screen_results.RData")

# comp_screen_tables<-function(x1,x2,thr=0.2){
#   res=list()
#   for(j in 1:ncol(x1)){
#     res[[j]] = table(x1[,j]<thr,x2[,j]<thr)
#   }
#   return(res)
# }
# 
# for(nn in names(screen_results)){
#   print(comp_screen_tables(screen_results[[nn]][[2]],screen_results[[nn]][[3]]))
# }
# 
# nn = names(screen_results)[4]
# x1 = screen_results[[nn]][[2]][,6]
# x2 = screen_results[[nn]][[3]][,6]
# names(x1) = rownames(screen_results[[nn]][[2]])
# names(x2) = rownames(screen_results[[nn]][[3]])
# all(names(x1)==names(x2))
# sort(abs(x1-x2),decreasing = T)[1:10]



