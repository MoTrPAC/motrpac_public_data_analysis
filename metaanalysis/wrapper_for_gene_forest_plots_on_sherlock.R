.libPaths("~/R/packages") # add our local installed R packages
library(optparse)
option_list <- list(
  make_option("--input_dir", default="acute_blood/",
              help="Input dir"),
  make_option("--ind1",default=1,
              help="index 1 on sorted csv files"),
  make_option("--ind2",default=10,
              help="index 2 on sorted csv files")
)
opt <- parse_args(OptionParser(option_list=option_list))

args = commandArgs(trailingOnly=TRUE)
input_dir = opt$input_dir
print(paste("Working on input dir:",input_dir))
pngs_dir = paste0(input_dir,"/png/")
out_dir = paste0(input_dir,"/out/")
pngs_dir = gsub("//","/",pngs_dir)
out_dir = gsub("//","/",out_dir)
dir.create(pngs_dir,showWarnings = F)
dir.create(out_dir,showWarnings = F)

csv_files = list.files(input_dir,full.names = T)
csv_files =  csv_files [ grepl(".csv$",csv_files)]
csv_files = sort(csv_files)
opt$ind1 = max(1,opt$ind1)
opt$ind2 = min(opt$ind2,length(csv_files))

for (i in opt$ind1:opt$ind2){
  curr_gene = gsub(input_dir,"",csv_files[i])
  csv_files[i] = gsub("//","/",csv_files[i])
  print(paste("Analyzing:",csv_files[i]))
  cmd = paste(
    "Rscript",
    "copy_of_run_simple_me_model_on_table.R",
    "--csv", csv_files[i],
    "--out", paste0(pngs_dir,curr_gene),
    "--forest 2",
    "--verbose FALSE",
    ">",paste0(out_dir,curr_gene,".out")
  )
  system(cmd)
  #print(warnings())
}









