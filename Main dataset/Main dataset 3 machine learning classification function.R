library(dplyr)
library(limma)
library(randomForest)
library(e1071)
library(MASS)
#library(class)
#library(ElemStatLearn)
library(tidyr)
library(ggplot2)
library(edgeR)

otu_normalized<-normalization_f(otu_common_combine) ## main dataset
otu_normalized<-normalization_f(otu_common_species_order_mystery) ## mystery dataset

city_sampling<-function(data,fold){
  city<-unique(data$city)
  sample_table<-as.data.frame(matrix(NA,nrow=length(unique(data$city)),ncol=3))
  colnames(sample_table)<-c('city','number_of_samples','sampling')
  sample_table$city<-city
  sample_table$number_of_samples<-table(data$city)
  sample_table$sampling<-round(sample_table$number_of_samples/fold)
  samples<-NULL
  for(i in 1:length(sample_table$city)){
    samples<-c(samples,sample(which(data$city==sample_table[i,'city']),sample_table[i,'sampling']))
  }
  return(samples)
}
test_index<-city_sampling(otu_normalized,5)
otu_normalized_training<-otu_normalized[-test_index,]
otu_normalized_test<-otu_normalized[test_index,]

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


RF_CV_test<-function(data_training,data_test,iter){
  time_start<-Sys.time()
  n<-nrow(data_training)
  RF_error_CV<-numeric(iter)
  RF_error_test<-as.data.frame(matrix(NA,nrow=nrow(data_test),ncol=iter+3))
  colnames(RF_error_test)<-c('True',1:iter,'mode','check')
  RF_error_test$True<-data_test$city
  rf_importance<-NULL
  for(i in 1:iter){
    set.seed(i)
    atrain<-data_training[-i,]
    atest<-data_training[i,]
    rf_fit<-randomForest(factor(city)~.,data=atrain,ntree=1000,type='classifier')
    imp<-importance(rf_fit)
    rf_importance<-cbind(rf_importance,imp)
    rf_pred_cv<-predict(rf_fit,atest)
    RF_error_test[,i+1]<-predict(rf_fit,data_test)
    rf_check<-table(rf_pred_cv,truth=atest$city)
    ## error rate
    if(rownames(rf_check)[which(rf_check==1)]==colnames(rf_check)){RF_error_CV[i]<-0}else{RF_error_CV[i]<-1}
  }
  RF_error_test$mode<-apply(RF_error_test[,2:iter+1],1,getmode)
  RF_error_test$check<-ifelse(RF_error_test$True==RF_error_test$mode,1,0)
  ## importance score
  rf_importance_mean<-as.data.frame(apply(rf_importance,1,mean))
  colnames(rf_importance_mean)[1]<-'importance_mean'
  rf_importance_mean$taxon<-row.names(rf_importance_mean)
  rf_importance_mean<-rf_importance_mean[order(rf_importance_mean$importance_mean),]
  rf_importance_mean$taxon<-factor(rf_importance_mean$taxon,levels=rf_importance_mean$taxon)
  plot_importance_score<-ggplot(data=rf_importance_mean,aes(importance_mean,taxon,col='red'))+
    geom_point()+
    xlab('importance score')+
    ylab('features')+
    theme(legend.position="none")
  RF_output<-list(RF_error_CV=RF_error_CV,
                  RF_error_test=RF_error_test,
                  plot_importance_score=plot_importance_score)
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  return(RF_output)
}

SVM_CV_test<-function(data_training,data_test,iter){
  time_start<-Sys.time()
  n<-nrow(data_training)
  svm_error_cv<-numeric(iter)
  svm_error_test<-as.data.frame(matrix(NA,nrow=nrow(data_test),ncol=iter+3))
  colnames(svm_error_test)<-c('True',1:iter,'mode','check')
  svm_error_test$True<-data_test$city
  for(i in 1:iter){
    set.seed(i)
    atrain<-data_training[-i,]
    atest<-data_training[i,]
    svm_fit<-best.svm(factor(city)~.,data=atrain)
    svm_pred_cv<-predict(svm_fit,atest)
    svm_error_test[,i+1]<-predict(svm_fit,data_test)
    svm_check<-table(svm_pred_cv,truth=atest$city)
    ## error rate
    if(rownames(svm_check)[which(svm_check==1)]==colnames(svm_check)){svm_error_cv[i]<-0}else{svm_error_cv[i]<-1}
  }
  svm_error_test$mode<-apply(svm_error_test[,2:iter+1],1,getmode)
  svm_error_test$check<-ifelse(svm_error_test$True==svm_error_test$mode,1,0)
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  svm_output<-list(svm_error_cv=svm_error_cv,svm_error_test=svm_error_test)
  return(svm_output)
}

LDA_CV_test<-function(data_training,data_test,iter){
  time_start<-Sys.time()
  n<-nrow(data_training)
  lda_error_cv<-numeric(iter)
  lda_error_test<-as.data.frame(matrix(NA,nrow=nrow(data_test),ncol=iter+3))
  colnames(lda_error_test)<-c('True',1:iter,'mode','check')
  lda_error_test$True<-data_test$city
  for(i in 1:iter){
    set.seed(i)
    atrain<-data_training[-i,]
    atest<-data_training[i,]
    lda_fit<-lda(factor(city)~.,data=atrain)
    lda_pred_cv<-predict(lda_fit,atest)
    lda_pred_test<-predict(lda_fit,data_test)
    lda_error_test[,i+1]<-lda_pred_test$class
    lda_check<-table(lda_pred_cv$class,truth=atest$city)
    ## error rate
    if(rownames(lda_check)[which(lda_check==1)]==colnames(lda_check)){lda_error_cv[i]<-0}else{lda_error_cv[i]<-1}
  }
  lda_error_test$mode<-apply(lda_error_test[,2:iter+1],1,getmode)
  lda_error_test$check<-ifelse(lda_error_test$True==lda_error_test$mode,1,0)
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  lda_output<-list(lda_error_cv=lda_error_cv,lda_error_test=lda_error_test)
  return(lda_output)
}

RF_test_result<-RF_CV_test(otu_normalized_training,otu_normalized_test,nrow(otu_normalized_training))
SVM_test_result<-SVM_CV_test(otu_normalized_training,otu_normalized_test,nrow(otu_normalized_training))
LDA_test_result<-LDA_CV_test(otu_normalized_training,otu_normalized_test,nrow(otu_normalized_training))

