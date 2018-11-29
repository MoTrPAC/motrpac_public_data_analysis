library(SmartSVA)
library(corrplot)
library(lme4)

run_wilcox<-function(x,y,...){
  x1 = x[y==y[1]]
  x2 = x[y!=y[1]]
  return(wilcox.test(x1,x2,...)$p.value)
}

run_test<-function(x,y,func=t.test,...){
  x1 = x[y==y[1]]
  x2 = x[y!=y[1]]
  return(func(x1,x2,...)$p.value)
}
#######################################
#######################################
#######################################
# Metabolomics 

# Some preprocessing of the rat pilot data (plasma)
input_dataset = "/Users/David/Desktop/MoTrPAC/data/pilot/pilot_rat_plasma/CAS Pilot Rat Plasma Compiled.txt"
data = read.delim(input_dataset,stringsAsFactors = F)

hmdb = data$HMDB..
hmdb_table = table(hmdb)
selected_hmdbs = names(hmdb_table)[hmdb_table>1]
selected_hmdbs = selected_hmdbs[grepl("^HMDB",selected_hmdbs)]
selected_rows = is.element(hmdb,set=selected_hmdbs)

X = data[selected_rows,c(paste("R",1:10,sep=""),paste("S",1:10,sep=""))]
meta1 = data[selected_rows,2:6]
table(meta1$HMDB..)
sort(meta1$HMDB)
meta2 = sapply(meta1$Full.Name,function(x)strsplit(x,split="_|\\s")[[1]])
site = sapply(meta2,function(x)x[1])
method = sapply(meta2,function(x)x[2])
tt = table(site,meta1$HMDB..)
selected_sites = c("Broad","UM")
selected_hmdbs = colnames(tt)[colSums(tt[selected_sites,]>0)==2]
rows2 = is.element(site,set=selected_sites) & is.element(meta1$HMDB..,set=selected_hmdbs)
X = X[rows2,]
meta1 = meta1[rows2,]
site = site[rows2]
method = method[rows2]
tt1 = table(site,meta1$HMDB..)
tt2 = table(method,meta1$HMDB..)
rowSums(tt2==2)
selected_method = "HILIC"
selected_hmdbs = colnames(tt2)[tt2[selected_method,]==2]
rows3 = is.element(method,set=selected_method) & is.element(meta1$HMDB..,set=selected_hmdbs)
X = X[rows3,]
meta1 = meta1[rows3,]
site=site[rows3];hmdb=meta1$HMDB..
newX = c();newS = c()
for(ss in unique(site)){
  inds = site==ss
  currX = X[inds,]
  rownames(currX) = hmdb[inds]
  colnames(currX) = paste(colnames(currX),ss,sep=";")
  currX = as.matrix(currX)
  ord = order(rownames(currX))
  currX = currX[ord,]
  if(length(newX)>0 && !(all(rownames(newX)==rownames(currX)))){
    print("rownames do not fit")
  }
  print(mean(currX,na.rm = T))
  newX = cbind(newX,currX)
  newS = c(newS,rep(ss,ncol(currX)))
}

# Analyze the data
# Correlations
newX = log(newX)
corrplot(cor(newX),tl.cex = 0.8)
# PCA
par(mfrow=c(1,2))
pca = prcomp(t(newX))
plot(pca)
cumsum(pca$sdev)/sum(pca$sdev)
cols = as.factor(newS)
cols2 = as.factor(grepl("^R",rownames(pcax)))
pcax = pca$x
plot(PC1~PC2,data=pcax,main="Blood plasma metabolites (25), explained_variance=0.63",col=cols,lwd=4,pch=20)
legend("center",legend = unique(cols),fill=unique(cols),cex=2)
boxplot(PC7~cols2,data=pcax,col=c("red","blue"),names=c("S","R (exercise)"),ylab="PC7")

d1 = data.frame(pcax,site=cols,treatment=cols2,subject = rep(1:10,4))
summary(lme(PC2~treatment,random = ~1|site + 1|subject,data=d1))

pca_ps = c()
for(j in 1:ncol(pcax)){
  pca_ps[j] = run_wilcox(pcax[,j],newS,paired=T)
}

library(PairedData)
currd = data.frame(pc=pcax[,1],group=newS)
broad <- subset(currd,  group == "Broad", pc,
                 drop = TRUE)
# subset weight data after treatment
um <- subset(currd,  group == "UM", pc,
                drop = TRUE)
# Plot paired data
library(PairedData)
pd <- paired(broad, um)
plot(pd, type = "profile",main="PC1: paired sample difference") + theme_bw()

# SVA: reconfirms PC1
pheno=rep("R",ncol(newX))
pheno[grepl("^S",colnames(newX))] = "S"
mdat= data.frame(pheno,newS)
mod = model.matrix(~as.factor(pheno), data=mdat)
mod0 = model.matrix(~1,data=mdat)
n.sv = num.sv(newX,mod,method="leek")
svobj = sva(newX,mod,mod0,n.sv=n.sv)
cor(pcax[,1],svobj$sv[,1])

Y = t(newX)
corrplot(cor(Y[1:20,],Y[21:40,]),order="hclust")
diag(cor(Y[1:20,],Y[21:40,]))

pvals = c()
for(ss in unique(newS)){
  currX = newX[,newS==ss]
  ord = order(colnames(currX))
  currX = currX[,ord]
  currpheno = grepl("^R",colnames(currX))
  ps = apply(currX,1,run_wilcox,y=currpheno)
  pvals = cbind(pvals,ps)
}
cor(pvals,method="spearman")
table(pvals[,1]<0.1,pvals[,2]<0.1)

# Adjust for SVA's factors
newXresids = c()
for(j in 1:nrow(newX)){
  l = lm(newX[j,]~svobj$sv[,1])
  newXresids = rbind(newXresids,l$residuals)
}
dim(newXresids)
hist(abs(newX-newXresids))

newX = newXresids
corrplot(cor(newXresids),order="hclust")
Y = t(newXresids)
corrplot(cor(Y[1:20,],Y[21:40,]),order="hclust")
diag(cor(Y[1:20,],Y[21:40,]))

pca = prcomp(t(newX))
pcax = pca$x
plot(PC1~PC2,data=pcax)
text(PC1~PC2,data=pcax,labels = rownames(pcax),pos=4)
pvals = c()
for(ss in unique(newS)){
  currX = newX[,newS==ss]
  ord = order(colnames(currX))
  currX = currX[,ord]
  currpheno = grepl("^R",colnames(currX))
  ps = apply(currX,1,run_wilcox,y=currpheno)
  pvals = cbind(pvals,ps)
}
cor(pvals,method="spearman")
table(pvals[,1]<0.1,pvals[,2]<0.1)

hmdb2name = unique(cbind(data$HMDB..,data$Featrue.name))
hmdb2name = hmdb2name[!is.na(hmdb2name[,1]) & grepl("hmdb",hmdb2name[,1],ignore.case = T),]
rownames(hmdb2name) = hmdb2name[,1]
our_pc = pca$rotation[,1]
our_pc = cbind(names(our_pc),hmdb2name[names(our_pc),2],our_pc)
write.table(our_pc,row.names = F,quote=F,sep=",")

#######################################
#######################################
#######################################
# RNA-seq 

# Some preprocessing of the rat pilot data (plasma)
load("/Users/David/Desktop/MoTrPAC/data/pilot/rna_seq/gene_matrices.RData")
mlab_gcounts_matrix = log(mlab_gcounts_matrix+1,base=2)
fpkm_matrix = log(fpkm_matrix+1,base=2)
boxplot(mlab_gcounts_matrix)
boxplot(fpkm_matrix)

corrplot(cor(fpkm_matrix),order="hclust")
corrplot(cor(mlab_gcounts_matrix),order="hclust")

m1 = fpkm_matrix[shared,]
m2 = mlab_gcounts_matrix[shared,]
colnames(m2) = gsub("_ReadsPerGene","",colnames(m2))
corrplot(cor(m1,m2[,colnames(m1)]))

library(preprocessCore)
m1_qnorm1 = cbind(normalize.quantiles.robust(m1[,1:10]),normalize.quantiles.robust(m1[,11:20]))
boxplot(m1_qnorm1)
colSums(m1_qnorm1==0)
colnames(m1_qnorm1) = colnames(m1)
rownames(m1_qnorm1) = rownames(m1)
corrplot(cor(m1_qnorm1))

rnaseq_pca = prcomp(t(m1_qnorm1))
plot(rnaseq_pca)
cumsum(rnaseq_pca$sdev)/sum(rnaseq_pca$sdev)
cols = rep("Muscle",ncol(m1))
cols[grepl("^A_",colnames(m1))] = "Adipose"
cols = as.factor(cols)
pcax = rnaseq_pca$x
plot(PC2~PC1,data=pcax,main="RNA-seq (20k), explained_variance=0.41",col=cols,lwd=4,pch=20)
legend("center",legend = unique(cols),fill=unique(cols),cex=2)
outlier_sample = which(pcax[,"PC2"]< -50)
colnames(m1)

x_muscle = m1_qnorm1[,11:20]
x_adipose = m1_qnorm1[,c(1:8,10)]
muscle_genes = apply(x_muscle>1,1,any)
adipose_genes = apply(x_adipose>1,1,any)
x_muscle = x_muscle[muscle_genes,]
x_adipose = x_adipose[adipose_genes,]
dim(x_muscle)
dim(x_adipose)
is_running = grepl(colnames(x_muscle),pattern = "^\\D_R")
muscle_ps = apply(x_muscle,1,run_test,y=is_running,paired=F,func=t.test)
is_running = grepl(colnames(x_adipose),pattern = "^\\D_R")
adipose_ps = apply(x_adipose,1,run_test,y=is_running,paired=F,func=t.test)
muscle_qs = p.adjust(muscle_ps,method='fdr')
adipose_qs = p.adjust(adipose_ps,method='fdr')
hist(muscle_ps)
hist(adipose_ps)
muscle_diff_genes = rownames(x_muscle)[muscle_ps < 0.01]
adipose_diff_genes = rownames(x_adipose)[adipose_ps < 0.01]
intersect(muscle_diff_genes,adipose_diff_genes)
all_diff_genes = union(muscle_diff_genes,adipose_diff_genes)
corrs1 = cor(t(m1_qnorm1[all_diff_genes,11:20]))
heatmap(corrs1)
tmp2 = m1_qnorm1[all_diff_genes,c(1:8,10)]
tmp2 = tmp2[apply(tmp2,1,sd)>0,]
corrs2 = cor(t(tmp2))
library(gplots)
heatmap.2(corrs2[adipose_diff_genes,adipose_diff_genes],trace="none",labels=NULL)
heatmap.2(corrs1[adipose_diff_genes,adipose_diff_genes],trace="none",labels=NULL)




