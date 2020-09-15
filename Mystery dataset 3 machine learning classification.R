library(limma)
library(randomForest)
library(e1071)
library(MASS)
#library(class)
#library(ElemStatLearn)
library(tidyr)
library(ggplot2)
library(edgeR)

####################################### normalization function ############################################
####################################### normalization function ############################################
####################################### normalization function ############################################
normalization_f<-function(data){
  data_0<-data[,which(colSums(data)!=0)]
  normalization_t<-data_0
  d0<-DGEList(normalization_t)
  group<-rep(unique(substr(colnames(normalization_t),1,3)),as.vector(table(substr(colnames(normalization_t),1,3))))
  mm<-model.matrix(~0 + group)
  y<-voom(d0, mm, plot = T)$E
  otu_normalized<-as.data.frame(t(y))
  otu_normalized$city<-substr(row.names(otu_normalized),1,3)
  names(otu_normalized)<-make.names(names(otu_normalized))
  return(otu_normalized)
}

##################################### Random Forest function ##############################################
##################################### Random Forest function ##############################################
##################################### Random Forest function ##############################################
RF_CV<-function(data,iter){
  time_start<-Sys.time()
  n<-nrow(data)
  RF_error<-numeric(iter)
  rf_importance<-NULL
  rf_predict<-data.frame()
  for(i in 1:iter){
    set.seed(i)
    k<-sample(n,1)
    atrain<-data[-k,]
    atest<-data[k,]
    rf_fit<-randomForest(factor(city)~.,data=atrain,ntree=1000,type='classifier')
    imp<-importance(rf_fit)
    rf_importance<-cbind(rf_importance,imp)
    rf_pred<-predict(rf_fit,atest)
    rf_check<-table(rf_pred,truth=atest$city)
    ## error rate
    if(rownames(rf_check)[which(rf_check==1)]==colnames(rf_check)){RF_error[i]<-0}else{RF_error[i]<-1}
    rf_predict[i,1]<-row.names(data)[k]
    rf_predict[i,2]<-rf_pred
  }
  rf_error_rate<-mean(RF_error)
  colnames(rf_predict)<-c('sample','predict')
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
  ## overall error rate
  rf_predict$city<-NA
  rf_predict$check<-NA
  for(i in 1:nrow(rf_predict)){
    rf_predict$city[i]<-strsplit(rf_predict[i,'sample'],split='_')[[1]][1]
    if(rf_predict$predict[i]==substr(rf_predict$city[i],1,3)){rf_predict$check[i]<-1}
    else{rf_predict$check[i]<-0}
  }
  ## error rate for each city
  rf_error_city<-data.frame()
  rf_city_order<-unique(rf_predict$city)[order(unique(rf_predict$city))]
  for(i in rf_city_order){
    rf_error_city[i,1]<-i
    ## error rate
    rf_error_city[i,2]<-1-mean(rf_predict[which(rf_predict$city==i),'check'])
  }
  colnames(rf_error_city)<-c('city','error_rate')
  rf_error_city$method<-'RF'
  rf_error_city<-rf_error_city[order(rf_error_city$error_rate),]
  plot_error_city<-ggplot(data=rf_error_city,aes(x=reorder(city,error_rate),y=error_rate,fill=city))+
    geom_bar(stat='identity')+
    xlab('city')+
    ylab('error rate')
  RF_output<-list(rf_error_rate=rf_error_rate,
                  plot_importance_score=plot_importance_score,
                  rf_predict=rf_predict,
                  rf_error_city=rf_error_city,
                  plot_error_city=plot_error_city)
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  return(RF_output)
}

######################################### SVM function ##################################################
######################################### SVM function ##################################################
######################################### SVM function ##################################################
SVM_CV<-function(data,iter){
  time_start<-Sys.time()
  n<-nrow(data)
  SVM_error<-numeric(iter)
  svm_predict<-data.frame()
  for(i in 1:iter){
    set.seed(i)
    k<-sample(n,1)
    atrain<-data[-k,]
    atest<-data[k,]
    svm_fit<-best.svm(factor(city)~.,data=atrain)
    svm_pred<-predict(svm_fit,atest)
    svm_check<-table(svm_pred,truth=atest$city)
    ## error rate
    if(rownames(svm_check)[which(svm_check==1)]==colnames(svm_check)){SVM_error[i]<-0}else{SVM_error[i]<-1}
    svm_predict[i,1]<-row.names(data)[k]
    svm_predict[i,2]<-svm_pred
  }
  colnames(svm_predict)<-c('sample','predict')
  svm_error_rate<-mean(SVM_error)
  ## overall error rate
  svm_predict$city<-NA
  svm_predict$check<-NA
  for(i in 1:nrow(svm_predict)){
    svm_predict$city[i]<-strsplit(svm_predict[i,'sample'],split='_')[[1]][1]
    if(svm_predict$predict[i]==substr(svm_predict$city[i],1,3)){svm_predict$check[i]<-1}
    else{svm_predict$check[i]<-0}
  }
  ## error rate for each city
  svm_error_city<-data.frame()
  svm_city_order<-unique(svm_predict$city)[order(unique(svm_predict$city))]
  for(i in svm_city_order){
    svm_error_city[i,1]<-i
    ## error rate
    svm_error_city[i,2]<-1-mean(svm_predict[which(svm_predict$city==i),'check'])
  }
  colnames(svm_error_city)<-c('city','error_rate')
  svm_error_city$method<-'SVM'
  svm_error_city<-svm_error_city[order(svm_error_city$error_rate),]
  plot_error_city<-ggplot(data=svm_error_city,aes(x=reorder(city,error_rate),y=error_rate,fill=city))+
    geom_bar(stat='identity')+
    xlab('city')+
    ylab('error rate')
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  SVM_output<-list(svm_error_rate=svm_error_rate,svm_predict=svm_predict,svm_error_city=svm_error_city,plot_error_city=plot_error_city)
  return(SVM_output)
}

######################################### LDA function ##################################################
######################################### LDA function ##################################################
######################################### LDA function ##################################################
LDA_CV<-function(data,iter){
  time_start<-Sys.time()
  n<-nrow(data)
  LDA_error<-numeric(iter)
  lda_predict<-data.frame()
  for(i in 1:iter){
    set.seed(i)
    k<-sample(n,1)
    atrain<-data[-k,]
    atest<-data[k,]
    lda_fit<-lda(factor(city)~.,data=atrain)
    lda_pred<-predict(lda_fit,atest)
    lda_check<-table(lda_pred$class,truth=atest$city)
    if(rownames(lda_check)[which(lda_check==1)]==colnames(lda_check)){LDA_error[i]<-0}else{LDA_error[i]<-1}
    lda_predict[i,1]<-row.names(data)[k]
    lda_predict[i,2]<-lda_pred$class
  }
  colnames(lda_predict)<-c('sample','predict')
  lda_error_rate<-mean(LDA_error)
  ## overall error rate
  lda_predict$city<-NA
  lda_predict$check<-NA
  for(i in 1:nrow(lda_predict)){
    lda_predict$city[i]<-strsplit(lda_predict[i,'sample'],split='_')[[1]][1]
    if(lda_predict$predict[i]==substr(lda_predict$city[i],1,3)){lda_predict$check[i]<-1}
    else{lda_predict$check[i]<-0}
  }
  ## error rate for each city
  lda_error_city<-data.frame()
  lda_city_order<-unique(lda_predict$city)[order(unique(lda_predict$city))]
  for(i in lda_city_order){
    lda_error_city[i,1]<-i
    ## error rate
    lda_error_city[i,2]<-1-mean(lda_predict[which(lda_predict$city==i),'check'])
  }
  colnames(lda_error_city)<-c('city','error_rate')
  lda_error_city$method<-'LDA'
  lda_error_city<-lda_error_city[order(lda_error_city$error_rate),]
  plot_error_city<-ggplot(data=lda_error_city,aes(x=reorder(city,error_rate),y=error_rate,fill=city))+
    geom_bar(stat='identity')+
    xlab('city')+
    ylab('error rate')
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  LDA_output<-list(lda_error_rate=lda_error_rate,lda_predict=lda_predict,lda_error_city=lda_error_city,plot_error_city=plot_error_city)
  return(LDA_output)
}

############################################# mystery dataset ############################################
############################################# mystery dataset ############################################
############################################# mystery dataset ############################################
## species-family
otu_common_species_family_mystery<-rbind(otu_common_species_mystery,otu_common_family_mystery)
otu_normalized<-normalization_f(otu_common_species_family_mystery)
result_species_family_mystery_RF<-RF_CV(otu_normalized,1000)
result_species_family_mystery_SVM<-SVM_CV(otu_normalized,1000)
result_species_family_mystery_LDA<-LDA_CV(otu_normalized,1000)

## species-order
otu_common_species_order_mystery<-rbind(otu_common_species_mystery,otu_common_order_mystery)
otu_normalized<-normalization_f(otu_common_species_order_mystery)
result_species_order_mystery_RF<-RF_CV(otu_normalized,1000)
result_species_order_mystery_SVM<-SVM_CV(otu_normalized,1000)
result_species_order_mystery_LDA<-LDA_CV(otu_normalized,1000)

importance_plot<-result_species_order_mystery_RF$plot_importance_score
print(importance_plot)
tiff(paste0( "C:/Users/41344/Desktop/Datta/CAMDA 2019/Data/Data of mystery dataset/Plot for mystery/","plot_importance_score",".tif"), height=1400, width=2100, res=250, units="px"); print(importance_plot); dev.off()


## family-order
otu_common_family_order_mystery<-rbind(otu_common_family_mystery,otu_common_order_mystery)
otu_normalized<-normalization_f(otu_common_family_order_mystery)
result_family_order_mystery_RF<-RF_CV(otu_normalized,1000)
result_family_order_mystery_SVM<-SVM_CV(otu_normalized,1000)
result_family_order_mystery_LDA<-LDA_CV(otu_normalized,1000)

## combine
otu_common_combine_mystery<-rbind(otu_common_species_mystery,otu_common_family_mystery,otu_common_order_mystery)
otu_normalized<-normalization_f(otu_common_combine_mystery)
result_common_combine_mystery_RF<-RF_CV(otu_normalized,1000)
result_common_combine_mystery_SVM<-SVM_CV(otu_normalized,1000)
result_common_combine_mystery_LDA<-LDA_CV(otu_normalized,1000)

## species/family/order
## species 
otu_normalized<-normalization_f(otu_common_species_mystery)
result_species_mystery_RF<-RF_CV(otu_normalized,1000)
result_species_mystery_SVM<-SVM_CV(otu_normalized,1000)
result_species_mystery_LDA<-LDA_CV(otu_normalized,1000)
## family
otu_normalized<-normalization_f(otu_common_family_mystery)
result_family_mystery_RF<-RF_CV(otu_normalized,1000)
result_family_mystery_SVM<-SVM_CV(otu_normalized,1000)
result_family_mystery_LDA<-LDA_CV(otu_normalized,1000)
## order
otu_normalized<-normalization_f(otu_common_order_mystery)
result_order_mystery_RF<-RF_CV(otu_normalized,1000)
result_order_mystery_SVM<-SVM_CV(otu_normalized,1000)
result_order_mystery_LDA<-LDA_CV(otu_normalized,1000)

## plot three - based on species-order
error_city<-rbind(result_species_order_mystery_RF$rf_error_city,result_species_order_mystery_SVM$svm_error_city,result_species_order_mystery_LDA$lda_error_city)
error_city$method<-as.factor(error_city$method)
error_city$method<-factor(error_city$method, levels=c('RF','SVM','LDA'), labels=c('RF','SVM','LDA')) 
error_city$mean<-0
for(i in 1:nrow(error_city)){
  error_city[i,'mean']<-mean(error_city[which(error_city$city==error_city[i,'city']),2])
}
error_city$city<-factor(error_city$city,levels=error_city[order(error_city[1:8,4]),'city'],labels=error_city[order(error_city[1:8,4]),'city'])
plot_three<-ggplot(data=error_city,aes(city,error_rate,fill=city))+geom_bar(stat='identity')+ylab('error rate')+facet_grid(.~method)
print(plot_three)

tiff(paste0( "C:/Users/41344/Desktop/Datta/CAMDA 2019/Data/Data of mystery dataset/Plot for mystery/","error_rate_three",".tif"), height=1400, width=2100, res=140, units="px"); print(plot_three); dev.off()

############################# data saved #############################
############################# data saved #############################
############################# data saved #############################
############################# data saved #############################
save(list=c('otu_common_species_order_mystery'),
     file='C:/Users/41344/Desktop/Datta/CAMDA 2019/Data/Data of mystery dataset/Mystery dataset 4 data for PCA.RData')
save(list=c('otu_common_species_order_mystery'),
     file='C:/Users/41344/Desktop/Datta/CAMDA 2019/Data/Data of mystery dataset/Mystery dataset 5 data for PCoA.RData')
save(list=c('otu_common_species_order_mystery'),
     file='C:/Users/41344/Desktop/Datta/CAMDA 2019/Data/Data of mystery dataset/Mystery dataset 6 data for ANCOM.RData')



