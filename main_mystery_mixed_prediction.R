load('otu_common_combine.RData')
load('Mystery dataset for prediction.RData')
library(limma)
library(randomForest)
library(e1071)
library(MASS)
#library(class)
#library(ElemStatLearn)
library(tidyr)
library(ggplot2)
library(edgeR)

mystery_common_species<-Reduce(intersect,list(row.names(otu_species_Brisbane),
                                             row.names(otu_species_Doha),
                                             row.names(otu_species_Kiev),
                                             row.names(otu_species_Oslo),
                                             row.names(otu_species_Paris),
                                             row.names(otu_species_Rio_de_Janerio),
                                             row.names(otu_species_Santiago_de_Chile),
                                             row.names(otu_species_Vienna)))
mystery_common_family<-Reduce(intersect,list(row.names(otu_family_Brisbane),
                                             row.names(otu_family_Doha),
                                             row.names(otu_family_Kiev),
                                             row.names(otu_family_Oslo),
                                             row.names(otu_family_Paris),
                                             row.names(otu_family_Rio_de_Janerio),
                                             row.names(otu_family_Santiago_de_Chile),
                                             row.names(otu_family_Vienna)))
mystery_common_order<-Reduce(intersect,list(row.names(otu_order_Brisbane),
                                             row.names(otu_order_Doha),
                                             row.names(otu_order_Kiev),
                                             row.names(otu_order_Oslo),
                                             row.names(otu_order_Paris),
                                             row.names(otu_order_Rio_de_Janerio),
                                             row.names(otu_order_Santiago_de_Chile),
                                             row.names(otu_order_Vienna)))
common_species<-intersect(mystery_common_species,row.names(otu_common_species))
common_family<-intersect(mystery_common_family,row.names(otu_common_family))
common_order<-intersect(mystery_common_order,row.names(otu_common_order))

otu_common_Brisbane<-rbind(otu_family_Brisbane[common_family,-ncol(otu_family_Brisbane)],otu_order_Brisbane[common_order,-ncol(otu_family_Brisbane)])
otu_common_Doha<-rbind(otu_family_Doha[common_family,-ncol(otu_family_Doha)],otu_order_Doha[common_order,-ncol(otu_family_Doha)])
otu_common_Kiev<-rbind(otu_family_Kiev[common_family,-ncol(otu_family_Kiev)],otu_order_Kiev[common_order,-ncol(otu_family_Kiev)])
otu_common_Oslo<-rbind(otu_family_Oslo[common_family,-ncol(otu_family_Oslo)],otu_order_Oslo[common_order,-ncol(otu_family_Oslo)])
otu_common_Paris<-rbind(otu_family_Paris[common_family,-ncol(otu_family_Paris)],otu_order_Paris[common_order,-ncol(otu_family_Paris)])
otu_common_Rio_de_Janerio<-rbind(otu_family_Rio_de_Janerio[common_family,-ncol(otu_family_Rio_de_Janerio)],otu_order_Rio_de_Janerio[common_order,-ncol(otu_family_Rio_de_Janerio)])
otu_common_Santiago_de_Chile<-rbind(otu_family_Santiago_de_Chile[common_family,-ncol(otu_family_Santiago_de_Chile)],otu_order_Santiago_de_Chile[common_order,-ncol(otu_family_Santiago_de_Chile)])
otu_common_Vienna<-rbind(otu_family_Vienna[common_family,-ncol(otu_family_Vienna)],otu_order_Vienna[common_order,-ncol(otu_family_Vienna)])

otu_common_main<-rbind(otu_common_family[common_family,],otu_common_order[common_order,])

otu_common_mystery<-cbind(otu_common_Brisbane,
                          otu_common_Doha,
                          otu_common_Kiev,
                          otu_common_Oslo,
                          otu_common_Paris,
                          otu_common_Rio_de_Janerio,
                          otu_common_Santiago_de_Chile,
                          otu_common_Vienna)

otu_common_all<-cbind(otu_common_main,otu_common_mystery)
otu_common_all_0<-otu_common_all[,-which(colSums(otu_common_all)==0)]

#### normalization ####
normalization_t<-otu_common_all_0
d0<-DGEList(normalization_t)
group<-c(substr(colnames(normalization_t[,1:298]),1,3),rep('mystery',55))
mm<-model.matrix(~0 + group)
y<-voom(d0, mm, plot = T)$E
otu_normalized<-as.data.frame(t(y))
otu_normalized$city<-c(substr(colnames(normalization_t[,1:298]),1,3),rep('mystery',55))
names(otu_normalized)<-make.names(names(otu_normalized))
rm(d0,mm,normalization_t,y,group)


error_rate<-as.data.frame(matrix(0,nrow=1000,ncol=3))
colnames(error_rate)<-c('RF','SVM','LDA')

for(i in 1:1000){
set.seed(i)
sample_Brisbane<-colnames(otu_common_Brisbane)[sample(1:ncol(otu_common_Brisbane),floor(ncol(otu_common_Brisbane)/2))]
sample_Doha<-colnames(otu_common_Doha)[sample(1:ncol(otu_common_Doha),floor(ncol(otu_common_Doha)/2))]
sample_Kiev<-colnames(otu_common_Kiev)[sample(1:ncol(otu_common_Kiev),floor(ncol(otu_common_Kiev)/2))]
sample_Oslo<-colnames(otu_common_Oslo)[sample(1:ncol(otu_common_Oslo),floor(ncol(otu_common_Oslo)/2))]
sample_Paris<-colnames(otu_common_Paris)[sample(1:ncol(otu_common_Paris),floor(ncol(otu_common_Paris)/2))]
sample_Rio_de_Janerio<-colnames(otu_common_Rio_de_Janerio)[sample(1:ncol(otu_common_Rio_de_Janerio),floor(ncol(otu_common_Rio_de_Janerio)/2))]
sample_Santiago_de_Chile<-colnames(otu_common_Santiago_de_Chile)[sample(1:ncol(otu_common_Santiago_de_Chile),floor(ncol(otu_common_Santiago_de_Chile)/2))]
sample_Vienna<-colnames(otu_common_Vienna)[sample(1:ncol(otu_common_Vienna),floor(ncol(otu_common_Vienna)/2))]

sample_test<-Reduce(union,list(sample_Brisbane,
                                sample_Doha,
                                sample_Kiev,
                                sample_Oslo,
                                sample_Paris,
                                sample_Rio_de_Janerio,
                                sample_Santiago_de_Chile,
                                sample_Vienna))
rm(sample_Brisbane,
   sample_Doha,
   sample_Kiev,
   sample_Oslo,
   sample_Paris,
   sample_Rio_de_Janerio,
   sample_Santiago_de_Chile,
   sample_Vienna)

otu_common_train<-otu_normalized[-which(colnames(otu_common_all_0) %in% sample_test),]
otu_common_test<-otu_normalized[which(colnames(otu_common_all_0) %in% sample_test),]

rm(sample_test)

combine_rf_fit<-randomForest(factor(city)~.,data=otu_common_train,ntree=1000,type='classifier')
predict_RF<-predict(combine_rf_fit,otu_common_test)
error_rate_RF<-1-length(which(predict_RF=='mystery'))/length(predict_RF)
error_rate[i,'RF']<-error_rate_RF

combine_svm_fit<-best.svm(factor(city)~.,data=otu_common_train)
predict_SVM<-predict(combine_svm_fit,otu_common_test)
error_rate_SVM<-1-length(which(predict_SVM=='mystery'))/length(predict_SVM)
error_rate[i,'SVM']<-error_rate_SVM

combine_lda_fit<-lda(factor(city)~.,data=otu_common_train)
predict_LDA<-predict(combine_lda_fit,otu_common_test)$class
error_rate_LDA<-1-length(which(predict_LDA=='mystery'))/length(predict_LDA)
error_rate[i,'LDA']<-error_rate_LDA

rm(otu_common_train,otu_common_test)
rm(combine_rf_fit,combine_svm_fit,combine_lda_fit,predict_RF,predict_SVM,predict_LDA,error_rate_RF,error_rate_SVM,error_rate_LDA)
}



