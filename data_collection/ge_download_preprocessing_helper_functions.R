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

