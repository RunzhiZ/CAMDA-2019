load('Main dataset 2 OTU tables based on common features.RData')
library(limma)
library(randomForest)
library(e1071)
library(MASS)
#library(class)
#library(ElemStatLearn)
library(tidyr)
library(ggplot2)
library(edgeR)

#### abundance ####
otu_all_species$abundance<-apply(otu_all_species[,-303],1,sum)
otu_all_family$abundance<-apply(otu_all_family[,-303],1,sum)
otu_all_order$abundance<-apply(otu_all_order[,-303],1,sum)

######## species ########
## common species ##
## species in at least 15 cities
otu_common_species_15<-otu_all_species[names(which(apply(all_species_count[,-17],1,sum)>=15)),1:302]
## species in at least 14 cities 
otu_common_species_14<-otu_all_species[names(which(apply(all_species_count[,-17],1,sum)>=14)),1:302]
## species in at least 13 cities
otu_common_species_13<-otu_all_species[names(which(apply(all_species_count[,-17],1,sum)>=13)),1:302]
## species in at least 12 cities
otu_common_species_12<-otu_all_species[names(which(apply(all_species_count[,-17],1,sum)>=12)),1:302]
## species in at least 11 cities
otu_common_species_11<-otu_all_species[names(which(apply(all_species_count[,-17],1,sum)>=11)),1:302]
## species in at least 10 cities
otu_common_species_10<-otu_all_species[names(which(apply(all_species_count[,-17],1,sum)>=10)),1:302]
## species in at least 9 cities
otu_common_species_9<-otu_all_species[names(which(apply(all_species_count[,-17],1,sum)>=9)),1:302]
## species in at least 8 cities
otu_common_species_8<-otu_all_species[names(which(apply(all_species_count[,-17],1,sum)>=8)),1:302]

######## family ########
## common family ##
## family in at least 15 cities
otu_common_family_15<-otu_all_family[names(which(apply(all_family_count[,-17],1,sum)>=15)),1:302]
## family in at least 14 cities
otu_common_family_14<-otu_all_family[names(which(apply(all_family_count[,-17],1,sum)>=14)),1:302]
## family in at least 13 cities
otu_common_family_13<-otu_all_family[names(which(apply(all_family_count[,-17],1,sum)>=13)),1:302]
## family in at least 12 cities
otu_common_family_12<-otu_all_family[names(which(apply(all_family_count[,-17],1,sum)>=12)),1:302]
## family in at least 11 cities
otu_common_family_11<-otu_all_family[names(which(apply(all_family_count[,-17],1,sum)>=11)),1:302]
## family in at least 10 cities
otu_common_family_10<-otu_all_family[names(which(apply(all_family_count[,-17],1,sum)>=10)),1:302]
## family in at least 9 cities
otu_common_family_9<-otu_all_family[names(which(apply(all_family_count[,-17],1,sum)>=9)),1:302]
## family in at least 8 cities
otu_common_family_8<-otu_all_family[names(which(apply(all_family_count[,-17],1,sum)>=8)),1:302]

######## order ########
## common order ##
## order in at least 15 cities
otu_common_order_15<-otu_all_order[names(which(apply(all_order_count[,-17],1,sum)>=15)),1:302]
## order in at least 14 cities
otu_common_order_14<-otu_all_order[names(which(apply(all_order_count[,-17],1,sum)>=14)),1:302]
## order in at least 13 cities
otu_common_order_13<-otu_all_order[names(which(apply(all_order_count[,-17],1,sum)>=13)),1:302]
## order in at least 12 cities
otu_common_order_12<-otu_all_order[names(which(apply(all_order_count[,-17],1,sum)>=12)),1:302]
## order in at least 11 cities
otu_common_order_11<-otu_all_order[names(which(apply(all_order_count[,-17],1,sum)>=11)),1:302]
## order in at least 10 cities
otu_common_order_10<-otu_all_order[names(which(apply(all_order_count[,-17],1,sum)>=10)),1:302]
## order in at least 9 cities
otu_common_order_9<-otu_all_order[names(which(apply(all_order_count[,-17],1,sum)>=9)),1:302]
## order in at least 8 cities
otu_common_order_8<-otu_all_order[names(which(apply(all_order_count[,-17],1,sum)>=8)),1:302]

##################################### normalization function ############################################
##################################### normalization function ############################################
##################################### normalization function ############################################
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
    atrain<-data[-i,]
    atest<-data[i,]
    rf_fit<-randomForest(factor(city)~.,data=atrain,ntree=1000,type='classifier')
    imp<-importance(rf_fit)
    rf_importance<-cbind(rf_importance,imp)
    rf_pred<-predict(rf_fit,atest)
    rf_check<-table(rf_pred,truth=atest$city)
    ## error rate
    if(rownames(rf_check)[which(rf_check==1)]==colnames(rf_check)){RF_error[i]<-0}else{RF_error[i]<-1}
    rf_predict[i,1]<-row.names(data)[i]
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
    xlab('Importance score')+
    ylab('Features')+
    theme(legend.position="none")
  ## overall error rate
  rf_predict$city<-NA
  rf_predict$check<-NA
  for(i in 1:nrow(rf_predict)){
    rf_predict$city[i]<-strsplit(rf_predict[i,'sample'],split='_')[[1]][1]
    if(rf_predict$predict[i]==rf_predict$city[i]){rf_predict$check[i]<-1}
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
    xlab('City')+
    ylab('Error rate')
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

##################################### SVM function ##################################################
##################################### SVM function ##################################################
##################################### SVM function ##################################################
SVM_CV<-function(data,iter){
  time_start<-Sys.time()
  n<-nrow(data)
  SVM_error<-numeric(iter)
  svm_predict<-data.frame()
  for(i in 1:iter){
    set.seed(i)
    atrain<-data[-i,]
    atest<-data[i,]
    svm_fit<-best.svm(factor(city)~.,data=atrain)
    svm_pred<-predict(svm_fit,atest)
    svm_check<-table(svm_pred,truth=atest$city)
    ## error rate
    if(rownames(svm_check)[which(svm_check==1)]==colnames(svm_check)){SVM_error[i]<-0}else{SVM_error[i]<-1}
    svm_predict[i,1]<-row.names(data)[i]
    svm_predict[i,2]<-svm_pred
  }
  colnames(svm_predict)<-c('sample','predict')
  svm_error_rate<-mean(SVM_error)
  ## overall error rate
  svm_predict$city<-NA
  svm_predict$check<-NA
  for(i in 1:nrow(svm_predict)){
    svm_predict$city[i]<-strsplit(svm_predict[i,'sample'],split='_')[[1]][1]
    if(svm_predict$predict[i]==svm_predict$city[i]){svm_predict$check[i]<-1}
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
    xlab('City')+
    ylab('Error rate')
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  SVM_output<-list(svm_error_rate=svm_error_rate,svm_predict=svm_predict,svm_error_city=svm_error_city,plot_error_city=plot_error_city)
  return(SVM_output)
}

##################################### LDA function ##################################################
##################################### LDA function ##################################################
##################################### LDA function ##################################################
LDA_CV<-function(data,iter){
  time_start<-Sys.time()
  n<-nrow(data)
  LDA_error<-numeric(iter)
  lda_predict<-data.frame()
  for(i in 1:iter){
    set.seed(i)
    atrain<-data[-i,]
    atest<-data[i,]
    lda_fit<-lda(factor(city)~.,data=atrain)
    lda_pred<-predict(lda_fit,atest)
    lda_check<-table(lda_pred$class,truth=atest$city)
    if(rownames(lda_check)[which(lda_check==1)]==colnames(lda_check)){LDA_error[i]<-0}else{LDA_error[i]<-1}
    lda_predict[i,1]<-row.names(data)[i]
    lda_predict[i,2]<-lda_pred$class
  }
  colnames(lda_predict)<-c('sample','predict')
  lda_error_rate<-mean(LDA_error)
  ## overall error rate
  lda_predict$city<-NA
  lda_predict$check<-NA
  for(i in 1:nrow(lda_predict)){
    lda_predict$city[i]<-strsplit(lda_predict[i,'sample'],split='_')[[1]][1]
    if(lda_predict$predict[i]==lda_predict$city[i]){lda_predict$check[i]<-1}
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
    xlab('City')+
    ylab('Error rate')
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  LDA_output<-list(lda_error_rate=lda_error_rate,lda_predict=lda_predict,lda_error_city=lda_error_city,plot_error_city=plot_error_city)
  return(LDA_output)
}

##################################### ubiquity #################################################
##################################### ubiquity #################################################
##################################### ubiquity #################################################
uniquity_species_rank<-otu_all_species[order(-otu_all_species$ubiquity),]
uniquity_family_rank<-otu_all_family[order(-otu_all_family$ubiquity),]
uniquity_order_rank<-otu_all_order[order(-otu_all_order$ubiquity),]
ll<-c(10,20,30,50,100,150)
#### RF ####
#### RF ####
#### RF ####
#### RF ####
otu_rf_ubi<-as.data.frame(matrix(nrow=6,ncol=3))
colnames(otu_rf_ubi)<-c('species','family','order')
row.names(otu_rf_ubi)<-paste0('top_',ll)
## species
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_species_rank[1:ll[j],1:302])
  otu_rf_ubi[j,'species']<-RF_CV(otu_normalized,nrow(otu_normalized))$rf_error_rate
}
## family
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_family_rank[1:ll[j],1:302])
  otu_rf_ubi[j,'family']<-RF_CV(otu_normalized,nrow(otu_normalized))$rf_error_rate
}

## order
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_order_rank[1:ll[j],1:302])
  otu_rf_ubi[j,'order']<-RF_CV(otu_normalized,nrow(otu_normalized))$rf_error_rate
}

#### svm ####
#### svm ####
#### svm ####
#### svm ####
otu_svm_ubi<-as.data.frame(matrix(nrow=6,ncol=3))
colnames(otu_svm_ubi)<-c('species','family','order')
row.names(otu_svm_ubi)<-paste0('top_',ll)
## species
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_species_rank[1:ll[j],1:302])
  otu_svm_ubi[j,'species']<-SVM_CV(otu_normalized,nrow(otu_normalized))$svm_error_rate
}
## family
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_family_rank[1:ll[j],1:302])
  otu_svm_ubi[j,'family']<-SVM_CV(otu_normalized,nrow(otu_normalized))$svm_error_rate
}
## order
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_order_rank[1:ll[j],1:302])
  otu_svm_ubi[j,'order']<-SVM_CV(otu_normalized,nrow(otu_normalized))$svm_error_rate
}

#### lda ####
#### lda ####
#### lda ####
#### lda ####
otu_lda_ubi<-as.data.frame(matrix(nrow=6,ncol=3))
colnames(otu_lda_ubi)<-c('species','family','order')
row.names(otu_lda_ubi)<-paste0('top_',ll)
## species
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_species_rank[1:ll[j],1:302])
  otu_lda_ubi[j,'species']<-LDA_CV(otu_normalized,nrow(otu_normalized))$lda_error_rate
}
## family
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_family_rank[1:ll[j],1:302])
  otu_lda_ubi[j,'family']<-LDA_CV(otu_normalized,nrow(otu_normalized))$lda_error_rate
}
## order
for(j in 1:6){
  otu_normalized<-normalization_f(uniquity_order_rank[1:ll[j],1:302])
  otu_lda_ubi[j,'order']<-LDA_CV(otu_normalized,nrow(otu_normalized))$lda_error_rate
}

############################################## at least x cities ##########################################
############################################## at least x cities ##########################################
############################################## at least x cities ##########################################
############################################## at least x cities ##########################################
## species
ll_species<-c('otu_common_species',
              'otu_common_species_15',
              'otu_common_species_14',
              'otu_common_species_13',
              'otu_common_species_12',
              'otu_common_species_11',
              'otu_common_species_10',
              'otu_common_species_9',
              'otu_common_species_8')
## family
ll_family<-c('otu_common_family',
             'otu_common_family_15',
             'otu_common_family_14',
             'otu_common_family_13',
             'otu_common_family_12',
             'otu_common_family_11',
             'otu_common_family_10',
             'otu_common_family_9',
             'otu_common_family_8')
## order
ll_order<-c('otu_common_order',
            'otu_common_order_15',
            'otu_common_order_14',
            'otu_common_order_13',
            'otu_common_order_12',
            'otu_common_order_11',
            'otu_common_order_10',
            'otu_common_order_9',
            'otu_common_order_8')

## RF
otu_rf_city<-as.data.frame(matrix(nrow=9,ncol=3))
colnames(otu_rf_city)<-c('species','family','order')
row.names(otu_rf_city)<-paste0('city_',16:8)
## SVM
otu_svm_city<-as.data.frame(matrix(nrow=9,ncol=3))
colnames(otu_svm_city)<-c('species','family','order')
row.names(otu_svm_city)<-paste0('city_',16:8)
## LDA
otu_lda_city<-as.data.frame(matrix(nrow=9,ncol=3))
colnames(otu_lda_city)<-c('species','family','order')
row.names(otu_lda_city)<-paste0('city_',16:8)

#### random forest ####
#### random forest ####
#### random forest ####
#### random forest ####
## species
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_species[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_rf_city[j,'species']<-RF_CV(otu_normalized,nrow(otu_normalized))$rf_error_rate
}
## family
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_family[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_rf_city[j,'family']<-RF_CV(otu_normalized,nrow(otu_normalized))$rf_error_rate
}
## order
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_order[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_rf_city[j,'order']<-RF_CV(otu_normalized,nrow(otu_normalized))$rf_error_rate
}

#### SVM ####
#### SVM ####
#### SVM ####
#### SVM ####
## species
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_species[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_svm_city[j,'species']<-SVM_CV(otu_normalized,nrow(otu_normalized))$svm_error_rate
}
## family
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_family[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_svm_city[j,'family']<-SVM_CV(otu_normalized,nrow(otu_normalized))$svm_error_rate
}
## order
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_order[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_svm_city[j,'order']<-SVM_CV(otu_normalized,nrow(otu_normalized))$svm_error_rate
}

#### lda ####
#### lda ####
#### lda ####
#### lda ####
## species
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_species[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_lda_city[j,'species']<-LDA_CV(otu_normalized,nrow(otu_normalized))$lda_error_rate
}
## family
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_family[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_lda_city[j,'family']<-LDA_CV(otu_normalized,nrow(otu_normalized))$lda_error_rate 
}  
## order
for(j in 1:9){
  data_for_normalization<-eval(parse(text=ll_order[j]))
  otu_normalized<-normalization_f(data_for_normalization)
  otu_lda_city[j,'order']<-LDA_CV(otu_normalized,nrow(otu_normalized))$lda_error_rate 
}  

######## combine ########
######## combine ########
######## combine ########
######## combine ########
otu_common_combine<-rbind(otu_common_species,otu_common_family,otu_common_order)

otu_common_family_order<-rbind(otu_common_family,otu_common_order)
otu_common_species_order<-rbind(otu_common_species,otu_common_order)
otu_common_species_family<-rbind(otu_common_species,otu_common_family)

##################################################### combine #################################################
##################################################### combine #################################################
##################################################### combine #################################################
##################################################### combine #################################################
#尝试不做normalization
#aa_test<-otu_common_combine
#aa_test<-aa_test[,which(colSums(aa_test)!=0)]
#aa_test<-as.data.frame(t(aa_test))
#aa_test$city<-substr(row.names(aa_test),1,3)
#otu_combine_RF_result_unnormalized<-RF_CV(aa_test,1000)
#otu_combine_SVM_result_unnormalzied<-SVM_CV(aa_test,1000)
#otu_combine_LDA_result_unnormalzied<-LDA_CV(aa_test,1000)
otu_normalized<-normalization_f(otu_common_combine)
otu_combine_RF_result<-RF_CV(otu_normalized,nrow(otu_normalized))
otu_combine_SVM_result<-SVM_CV(otu_normalized,nrow(otu_normalized))
otu_combine_LDA_result<-LDA_CV(otu_normalized,nrow(otu_normalized))

#################################################### species-family ###########################################
#################################################### species-family ###########################################
#################################################### species-family ###########################################
otu_normalized<-normalization_f(otu_common_species_family)
otu_species_family_RF_result<-RF_CV(otu_normalized,nrow(otu_normalized))
otu_species_family_SVM_result<-SVM_CV(otu_normalized,nrow(otu_normalized))
otu_species_family_LDA_result<-LDA_CV(otu_normalized,nrow(otu_normalized))

#################################################### species-order ###########################################
#################################################### species-order ###########################################
#################################################### species-order ###########################################
otu_normalized<-normalization_f(otu_common_species_order)
otu_species_order_RF_result<-RF_CV(otu_normalized,nrow(otu_normalized))
otu_species_order_SVM_result<-SVM_CV(otu_normalized,nrow(otu_normalized))
otu_species_order_LDA_result<-LDA_CV(otu_normalized,nrow(otu_normalized))

#################################################### family-order ###########################################
#################################################### family-order ###########################################
#################################################### family-order ###########################################
otu_normalized<-normalization_f(otu_common_family_order)
otu_family_order_RF_result<-RF_CV(otu_normalized,nrow(otu_normalized))
otu_family_order_SVM_result<-SVM_CV(otu_normalized,nrow(otu_normalized))
otu_family_order_LDA_result<-LDA_CV(otu_normalized,nrow(otu_normalized))


############################## plot for three methods ##############################
############################## plot for three methods ##############################
############################## plot for three methods ##############################
############################## plot for three methods ##############################
error_city<-rbind(otu_combine_RF_result$rf_error_city,otu_combine_SVM_result$svm_error_city,otu_combine_LDA_result$lda_error_city)
error_city$method<-as.factor(error_city$method)
error_city$method<-factor(error_city$method, levels=c('RF','SVM','LDA'), labels=c('RF','SVM','LDA'))
error_city$mean<-0
for(i in 1:nrow(error_city)){
  error_city[i,'mean']<-mean(error_city[which(error_city$city==error_city[i,'city']),2])
}
error_city$city<-factor(error_city$city,levels=error_city[order(error_city[1:16,4]),'city'],labels=error_city[order(error_city[1:16,4]),'city'])
plot_three<-ggplot(data=error_city,aes(city,error_rate,fill=city))+geom_bar(stat='identity')+ylab('error rate')+facet_grid(.~method)
print(plot_three)
tiff(paste0( "plot_three",".tif"), height=1400, width=2100, res=140, units="px"); print(plot_three); dev.off()

#### error rate plot city ####
error_rate_table_city<-as.data.frame(matrix(NA,nrow=81,ncol=4))
colnames(error_rate_table_city)<-c('error_rate','N','level','method')
error_rate_table_city[,1]<-c(otu_rf_city[,1],otu_rf_city[,2],otu_rf_city[,3],
                             otu_svm_city[,1],otu_svm_city[,2],otu_svm_city[,3],
                             otu_lda_city[,1],otu_lda_city[,2],otu_lda_city[,3])
error_rate_table_city[,2]<-rep(16:8,9)
error_rate_table_city[,3]<-rep(c('species','family','order'),times=3,each=9)
error_rate_table_city[,4]<-rep(c('RF','SVM','LDA'),each=27)
error_rate_table_city[,3]<-factor(error_rate_table_city[,3],levels=c('species','family','order'))
error_rate_table_city[,4]<-factor(error_rate_table_city[,4],levels=c('RF','SVM','LDA'))

plot_error_rate_table_city<-ggplot(data=error_rate_table_city,aes(x=N,y=error_rate))+
  geom_line()+facet_grid(level~method)+theme_bw()+xlab('N')+ylab('Error rate')
print(plot_error_rate_table_city)
tiff(paste0( "error_rate_city",".tif"), height=1400, width=1400, res=250, units="px");print(plot_error_rate_table_city);dev.off()

#### error rate plot ubi ####
error_rate_table_ubi<-as.data.frame(matrix(NA,nrow=54,ncol=4))
colnames(error_rate_table_ubi)<-c('error_rate','M','level','method')
error_rate_table_ubi[,1]<-c(otu_rf_ubi[,1],otu_rf_ubi[,2],otu_rf_ubi[,3],
                            otu_svm_ubi[,1],otu_svm_ubi[,2],otu_svm_ubi[,3],
                            otu_lda_ubi[,1],otu_lda_ubi[,2],otu_lda_ubi[,3])
error_rate_table_ubi[,2]<-rep(c(10,20,30,50,100,150),9)
error_rate_table_ubi[,3]<-rep(c('species','family','order'),times=3,each=6)
error_rate_table_ubi[,4]<-rep(c('RF','SVM','LDA'),each=18)
error_rate_table_ubi[,2]<-factor(error_rate_table_ubi[,2],levels=c('10','20','30','50','100','150'))
error_rate_table_ubi[,3]<-factor(error_rate_table_ubi[,3],levels=c('species','family','order'))
error_rate_table_ubi[,4]<-factor(error_rate_table_ubi[,4],levels=c('RF','SVM','LDA'))

plot_error_rate_table_ubi<-ggplot(data=error_rate_table_ubi,aes(x=M,y=error_rate))+
  geom_path()+facet_grid(level~method)+theme_bw()+xlab('M')+ylab('Error rate')
print(plot_error_rate_table_ubi)
tiff(paste0( "error_rate_ubi",".tif"), height=1400, width=1400, res=250, units="px");print(plot_error_rate_table_ubi);dev.off()

#### sequencing depth ####
otu_common_combine_0<-otu_common_combine[,which(colSums(otu_common_combine)!=0)]
otu_count<-as.data.frame(cbind(apply(otu_common_combine_0,2,sum),substr(colnames(otu_common_combine_0),1,3)))
colnames(otu_count)<-c('count','city')
otu_count$count<-as.numeric(otu_count$count)
otu_avg_count<-otu_count %>% group_by(city) %>% summarise(avg_count=mean(count))
error_city_dup<-error_city[!duplicated(error_city$city), ]
sequencing_depth<-inner_join(otu_avg_count,error_city_dup,by=c('city'='city'))

plot_sequencing_depth<-ggplot(data=sequencing_depth,aes(x=avg_count,y=mean,col=city))+geom_point(size=3)+xlab('Count')+ylab('Error rate')
print(plot_sequencing_depth)
tiff(paste0( "plot_sequencing_depth",".tif"), height=1400, width=1400, res=230, units="px");print(plot_sequencing_depth);dev.off()


#### details of the results of classification ####
PRED_RF<-otu_combine_RF_result$rf_predict
PRED_SVM<-otu_combine_SVM_result$svm_predict
PRED_LDA<-otu_combine_LDA_result$lda_predict
PRED_TABLE<-cbind(PRED_RF[,c(1,3,2)],PRED_SVM[,2],PRED_LDA[,2])
colnames(PRED_TABLE)[3:5]<-c('predict_RF','predict_SVM','predict_LDA')


#### Importance score ####
RF_importance<-otu_combine_RF_result$plot_importance_score
tiff(paste0( "RF_importance_main",".tif"), height=1400, width=1400, res=230, units="px");print(RF_importance);dev.off()

#### data for the following analyses ####
save(list=c('otu_common_combine',
            'otu_common_species',
            'otu_common_family',
            'otu_common_order'),
     file='otu_common_combine.RData')
