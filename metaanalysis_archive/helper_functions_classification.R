library(e1071);library(randomForest)
################## Implementation of wrappers for base classifiers ####################################
# These are used within some multilabel algorithms, most notably the naive binary relevance method.
# We implement below wrappers for binary and non binary data. These methods
# First run feature selection, record the selected clustering and use them to learn the
# classifier. 
# The non-binary data analysis is based on the CMA package in R.
# For binary data we used Fisher's exact test.
# The main reason for these wrappers: the need to save running time.
#source("https://bioconductor.org/biocLite.R")
#biocLite("CMA")
library(CMA);library(randomForest)
featureSelectionClassifier<-function(x,y,classification_function=randomForest,
                                     features_inds_to_keep=c(),fsmethod="limma",numFeatures=100,...){
  selection = GeneSelection(x,y,method=fsmethod,trace=F)
  selected_inds = toplist(selection,k=numFeatures,show=F)$index
  selected_inds = union(selected_inds,features_inds_to_keep)
  #print ("features selected")
  classifier = classification_function(x[,selected_inds],y,...)
  #print (paste("classifier learned",class(classifier)))
  obj = list(inds=selected_inds,classifier=classifier)
  class(obj) = "featureSelectionClassifier"
  return (obj)
}

predict.featureSelectionClassifier<-function(obj,newx,...){
  print ("in fs predict")
  args = c(list(obj$classifier,newx[,obj$inds]),...)
  return (do.call(predict,args=args))
}

get_hg_pval<-function(x,y,...){
  tab = table(x,y)
  if (length(tab) !=4){return (1)}
  return (fisher.test(tab,...)$p.value)
}

featureSelectionClassifier_binary_data<-function(x,y,
                                                 classification_function=svm,numFeatures=250,features_to_keep=c(),...){
  f_pvals = apply(x,2,get_hg_pval,y=y,alternative = 'greater')
  names(f_pvals) = colnames(x)
  selected_features1 = names(sort(f_pvals)[1:numFeatures])
  f_pvals = apply(x,2,get_hg_pval,y=y,alternative = 'less')
  names(f_pvals) = colnames(x)
  selected_features2 = names(sort(f_pvals)[1:numFeatures])
  selected_features = union(selected_features1,selected_features2)
  #print (paste(length(selected_features)," features selected"))
  classifier = classification_function(as.matrix(x[,selected_features]),y,...)
  #print (paste("classifier learned",class(classifier)))
  obj = list(inds=selected_features,classifier=classifier)
  class(obj) = "featureSelectionClassifier"
  return (obj)
}

calcAupr <- function(pred, gs,roc=F,useAbs=T) {
  ord.idx = NULL
  if (useAbs){
    ord.idx <- order(abs(pred), decreasing = T)
  }
  else{
    ord.idx <- order(pred, decreasing = T)
  }
  prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx))) #also known as positive predictive value
  rec  <- cumsum(gs[ord.idx]) / sum(gs)                     #also know as true positive rate
  fpr  <- cumsum(gs[ord.idx] == 0) / (length(gs) - sum(gs)) #false positive rate
  prec <- c(prec[1], prec)
  rec <- c(0, rec)
  fpr <- c(0, fpr)
  if (!roc){
    aupr <- areaUnderCurve(rec, prec)
    return (aupr)
  }
  #plot(fpr,rec,type='l')
  auroc <- areaUnderCurve(fpr, rec)
  return(auroc)
}

areaUnderCurve <- function(x, y) {
  dx <- diff(x)
  my <- y[1:(length(y) - 1)] + diff(y) / 2
  return(sum(dx * my))
}

leave_study_out2 <-function(y,x,d,func=run_lm,pred_args=list(),...){
  preds_list = list()
  for(batch in unique(d)){
    #print(batch)
    inds = which(d==batch)
    trx = x[-inds,];try = y[-inds]
    m = func(y=try,x=trx,...)
    p = NULL
    try({
      #print(max(inds));print(dim(x))
      args = c(list(m,x[inds,]),pred_args)
      p = do.call(predict,args)
    })
    if(is.null(p)){
      args = c(list(m,data.frame(x[inds,],check.names=F)),pred_args)
      p = do.call(predict,args)
    }
    preds_list[[batch]] = list(preds=p,real=y[inds])
  }
  return(preds_list)
}

simple_sampling_based_learner<-function(x,y,d,func = featureSelectionClassifier,reps=10,max_num_samples = 1000,...){
  y_table = table(y)
  min_class = names(y_table[y_table==min(y_table)])
  min_class_size = y_table[y_table==min(y_table)]
  positives = rownames(x)[y==min_class]
  datasets = unique(d[positives])
  negatives = rownames(x)[y!=min_class & is.element(d[rownames(x)],set=datasets)]
  bgcs = setdiff(rownames(x),positives);bgcs = setdiff(bgcs,negatives)
  classifiers = list()
  samp_size = min(c(length(positives),max_num_samples))
  for (i in 1:reps){
    newbgcs = NULL
    if(length(bgcs)>0){
      newbgcs = sample(bgcs)[1:min(samp_size,length(bgcs))]
    }
    pos_sample = sample(positives)[1:samp_size]
    neg_sample = NULL
    if(length(negatives)>0){
      neg_sample = sample(negatives)[1:min(samp_size,length(negatives))]
    }
    #print(paste("number of positives, negatives, and bgcs:"))
    #print(c(length(pos_sample),length(neg_sample),length(newbgcs)))
    curr_samples = c(pos_sample,neg_sample,newbgcs)
    curr_inds = is.element(rownames(x),curr_samples)
    newy = y[curr_inds]
    newx = x[curr_inds,]
    classifiers[[i]] = func(newx,as.factor(newy),...)
  }
  obj = classifiers
  class(obj) = "simple_sampling_based_learner"
  return(obj)
}
predict.simple_sampling_based_learner<-function(obj,x,...){
  preds = c()
  counts = 0
  for (i in 1:length(obj)){
    if (length(preds)==0){
      try({
        preds = predict(obj[[i]],x,...);
        counts = counts + 1
      })
    }
    else{
      try({
        preds = preds + predict(obj[[i]],x,...);
        counts = counts + 1
      })
    }
  }
  preds = preds / counts
  return (preds)
}

# Supervised analysis
# A configuration specifies the classification test.
# Each config has the following information in a list:
# The label name
# The features matrix name
# The sample set
# Logical: add covariates or not
# Classifier
# Classifier Arguments
svm_config = list(pred_args = list(probability=T),numFeatures=500,
                  classification_function=svm,probability=T,
                  kernel="linear")
rf_config = list(pred_args = list(type="vote"),numFeatures=500,
                 classification_function=randomForest)


##################################################################
############# Functions for analysis of configurations ###########
##################################################################

run_lso_tests_on_a_configuration_set<-function(configurations,classification_xs,classification_ys,d,x_cov){
  lso_results = list()
  for(cname in names(configurations)){
    config = configurations[[cname]]
    y_names = config$sample_set
    y = classification_ys[[config$yname]][y_names]
    x = classification_xs[[config$xname]][y_names,]
    # In case one of the classes has less than 2 datasets exclude it:
    tab = table(y,d[y_names])
    print(tab)
    classes_to_exclude = rownames(tab)[rowSums(tab>0)<2]
    if(length(classes_to_exclude)>0){
      print(paste("excluded classes in LSO: ",classes_to_exclude,collapse= " "))
      y_names = y_names[!is.element(y,set=classes_to_exclude)]
      y = classification_ys[[config$yname]][y_names]
      y = as.factor(as.character(y))
      names(y) = y_names
      x = classification_xs[[config$xname]][y_names,]
    }
    features_inds_to_keep=c()
    if(config$include_covs=="TRUE"){
      curr_covs = x_cov[,setdiff(colnames(x_cov),config$yname)]
      curr_covs = model.matrix(~.,curr_covs)[,-1]
      new_y_names = intersect(rownames(x),rownames(curr_covs))
      x = x[new_y_names,];y=y[new_y_names]
      x = cbind(curr_covs[new_y_names,],x)
      y_names = new_y_names
      features_inds_to_keep = 1:ncol(curr_covs)
    }
    curr_lso_test_args = c(list(y,x,d[y_names],func=featureSelectionClassifier,
                                features_inds_to_keep=features_inds_to_keep),config$classification_args)
    lso = do.call(leave_study_out2,curr_lso_test_args)
    lso_results[[cname]] = lso
  }
  return(lso_results)
}

get_standard_classification_performance_scores_for_results_list<-function(lso_results){
  performance_scores = list()
  for(cname in names(lso_results)){
    lso = lso_results[[cname]]
    lso_preds = lapply(lso,function(x)x$preds)
    if(! "matrix" %in% class(lso_preds[[1]])){
      lso_preds = lapply(lso,function(x)attr(x$preds,"probabilities"))
    }
    class_names = colnames(lso_preds[[1]])
    m = c()
    for(mm in lso_preds){m = rbind(m,mm[,class_names])}
    colnames(m) = class_names
    lso_preds = m
    aucs = c()
    for(class_name in class_names){
      lso_y = unlist(sapply(lso,function(x)x$real))
      class_inds = lso_y==class_name
      lso_y = rep(0,length(lso_y))
      lso_y[class_inds]=1
      lso_p = lso_preds[,class_name]
      roc = calcAupr(lso_p,as.numeric(lso_y),roc=T)
      proc = as.numeric(pROC::auc(lso_y, lso_p))
      aupr = calcAupr(lso_p,as.numeric(lso_y),roc=F)
      aucs = rbind(aucs,c(roc=roc,aupr=aupr,pROC=proc))
    }
    rownames(aucs) = class_names
    acc = lso_y == (lso_p>=0.5)
    acc = sum(acc)/length(acc)
    performance_scores[[cname]] = list(aucs=aucs,acc=acc)
  }
  return(performance_scores)
}

get_subject_performance_scores_from_results_list<-function(lso_results,name2subj=sample2subj,paired_tests_analysis=T,x_cov){
  subject_performance_scores = list()
  for(cname in names(lso_results)){
    subject_performance_scores[[cname]]=list()
    lso = lso_results[[cname]]
    lso_preds = lapply(lso,function(x)x$preds)
    if(! "matrix" %in% class(lso_preds[[1]])){
      lso_preds = lapply(lso,function(x)attr(x$preds,"probabilities"))
    }
    class_names = colnames(lso_preds[[1]])
    m = c()
    for(mm in lso_preds){m = rbind(m,mm[,class_names])}
    colnames(m) = class_names
    lso_preds = m
    curr_subjects = name2subj[rownames(lso_preds)]
    lso_y = unlist(sapply(lso,function(x)x$real))
    class1 = colnames(lso_preds)[1]
    if(paired_tests_analysis){
      get_profile_longi_line_plot(lso_preds[,class1],lso_y,curr_subjects) + 
        geom_line()
      xx1 = c();xx2=c()
      for(subj in unique(curr_subjects)){
        subj_inds = curr_subjects==subj
        if(sum(subj_inds)!=2){next} 
        curr_y = lso_y[subj_inds]
        if(length(unique(curr_y))!=2){next}
        curr_p = lso_preds[subj_inds,class1]
        xx1 = c(xx1,curr_p[curr_y==class1])
        xx2 = c(xx2,curr_p[curr_y!=class1])
      }
      curr_p = NA
      if(length(xx1)>0){curr_p = wilcox.test(xx1,xx2,paired=T)$p.value}
      d1 = data.frame(subjs=factor(curr_subjects),timev=lso_y,
                      tissue = x_cov[rownames(lso_preds),"tissue"],
                      training=x_cov[rownames(lso_preds),"training"],
                      sex=x_cov[rownames(lso_preds),"sex"])
      mixed_effects_model_emp_pval = NA
      # Should not work when splitting the database by one of the covariates
      # TODO: fixed later
      try({ 
        mixed_effects_model_emp_pval = get_mixed_effect_model_time_empirical_p(x=lso_preds[,class1],
                                                                               frm0 = x ~ factor(tissue) + factor(training) + factor(sex) + (1|subjs),
                                                                               frm1 = x ~ ordered(timev) + factor(tissue) + factor(training) + factor(sex) + (1|subjs),
                                                                               d1=d1,reps=100,min_reps = 100)
      })
      subject_performance_scores[[cname]][["paired_tests_analysis"]] = 
        list(x1=xx1,x2=xx2,binary_wilcox_p=curr_p,mixed_effects_emp_p = mixed_effects_model_emp_pval)
    }
    subj2acc = c()
    for(subj in curr_subjects){
      subj_inds = curr_subjects==subj
      if(sum(subj_inds)<2){next}
      curr_y = lso_y[subj_inds]
      if(!grepl("time",cname) || length(unique(curr_y))==1){
        curr_y = lso_y[subj_inds] == class1
        curr_p = mean(lso_preds[subj_inds,class1])>0.5
        subj2acc[subj] = sum(curr_p==curr_y[1])
      }else{
        curr_p = as.numeric(lso_preds[subj_inds,class1]>0.5)
        subj2acc[subj] = mean(curr_p==curr_y)
      }
    }
    subject_performance_scores[[cname]][["acc"]] = mean(subj2acc,na.rm=T)
    print(paste(cname,mean(subj2acc,na.rm=T)))
  }
  return(subject_performance_scores)
}


# x,y,d all have the same sample set in the same order
run_lm<-function(y,x,...){
  obj = lm(y~.,data=data.frame(y,x,check.names=F))
  return(obj)
}
leave_study_out <-function(y,x,d,func=run_lm,...){
  preds_list = list()
  for(batch in unique(d)){
    print(batch)
    inds = which(d==batch)
    trx = x[-inds,];try = y[-inds]
    m = func(try,trx)
    p = predict(m,data.frame(x[inds,],check.names=F))
    preds_list[[batch]] = cbind(p,y[inds])
    colnames(preds_list[[batch]]) = c("predictions","real.labels")
    if(!all(names(p)==rownames(x)[inds])){
      print ("ERROR in prediction")
      break
    }
  }
  return(preds_list)
}
library(randomForest)
run_rf<-function(y,x,...){randomForest(x,y,...)}