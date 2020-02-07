
# This is a simple code for comparing different versions of the data created by
# the different scripts in the data_collection folder

p1 = "~/Desktop/MoTrPAC/project_release_feb_2018/revision_feb_2020/"
p2 = "~/Desktop/MoTrPAC/project_release_feb_2018/data/"
load(paste(p1,"human_ge_cohort_preprocessed_db_gene_tables.RData",sep=""))
a1 = acute_gene_tables
l1 = longterm_gene_tables
load(paste(p2,"human_ge_cohort_preprocessed_db_gene_tables.RData",sep=""))
a2 = acute_gene_tables
l2 = longterm_gene_tables

shared_a = intersect(names(a1),names(a2))
for(gene in sample(shared_a)[1:3000]){
  print (gene)
  m1 = a1[[gene]]
  m2 = a2[[gene]]
  for(j in c(2,6:ncol(m1))){
    m1[[j]] = as.numeric(m1[[j]])
    m1[[j]] = round(m1[[j]],digits = 6)
    m2[[j]] = as.numeric(m2[[j]])
    m2[[j]] = round(m2[[j]],digits = 6)
  }
  m1[,2] = as.integer(m1[,2])
  m2[,2] = as.integer(m2[,2])
  rownames(m1) = apply(m1[,c(2:5,ncol(m1))],1,paste,collapse=",")
  rownames(m1) = gsub(" ","",rownames(m1))
  rownames(m2) = apply(m2[,c(2:5,ncol(m2))],1,paste,collapse=",")
  rownames(m2) = gsub(" ","",rownames(m2))
  stopifnot(all(rownames(m1) %in% rownames(m2)))
  # print(setdiff(rownames(m1),rownames(m2)))
  m2 = m2[rownames(m1),]
  stopifnot(sum(m2[,-c(1,8)]!=m1[,-c(1,8)],na.rm = T)==0)
  stopifnot(all(is.na(m1[is.na(m2)])))
  stopifnot(all(is.na(m2[is.na(m1)])))
}

shared_l = intersect(names(l1),names(l2))
for(gene in sample(shared_l)[1:3000]){
  print (gene)
  m1 = l1[[gene]]
  m2 = l2[[gene]]
  for(j in c(2,6:ncol(m1))){
    m1[[j]] = as.numeric(m1[[j]])
    m1[[j]] = round(m1[[j]],digits = 6)
    m2[[j]] = as.numeric(m2[[j]])
    m2[[j]] = round(m2[[j]],digits = 6)
  }
  m1[,2] = as.integer(m1[,2])
  m2[,2] = as.integer(m2[,2])
  rownames(m1) = apply(m1[,c(2:5,ncol(m1))],1,paste,collapse=",")
  rownames(m1) = gsub(" ","",rownames(m1))
  rownames(m2) = apply(m2[,c(2:5,ncol(m2))],1,paste,collapse=",")
  rownames(m2) = gsub(" ","",rownames(m2))
  stopifnot(all(rownames(m1) %in% rownames(m2)))
  # print(setdiff(rownames(m1),rownames(m2)))
  m2 = m2[rownames(m1),]
  stopifnot(sum(m2[,-c(1,8)]!=m1[,-c(1,8)],na.rm = T)==0)
  stopifnot(all(is.na(m1[is.na(m2)])))
  stopifnot(all(is.na(m2[is.na(m1)])))
}


