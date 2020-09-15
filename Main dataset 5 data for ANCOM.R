load('otu_common_combine.RData')
library(limma)
library(edgeR)
library(exactRankTests)
library(nlme)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

############## normalization function #################
############## normalization function #################
############## normalization function #################
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

#### normalization ####
#### normalization ####
#### normalization ####
otu_common_combine_0<-otu_common_combine[,which(colSums(otu_common_combine)!=0)]
otu_common_combine_normalized<-normalization_f(otu_common_combine_0)
data_for_ANCOM<-otu_common_combine_normalized

#### data for ANCOM ####
#### data for ANCOM ####
#### data for ANCOM ####
otu_ANCOM<-cbind(row.names(data_for_ANCOM),data_for_ANCOM)
colnames(otu_ANCOM)[1]<-'Sample.ID'
otu_ANCOM[,'Sample.ID']<-as.character(otu_ANCOM[,'Sample.ID'])
meta_ANCOM<-as.data.frame(cbind(otu_ANCOM[,'Sample.ID'],otu_ANCOM[,'city']))
colnames(meta_ANCOM)<-c('Sample.ID','city')
otu_ANCOM<-otu_ANCOM[,-ncol(otu_ANCOM)]

#### ANCOM function from ####
#### ANCOM function ####
#### ANCOM function ####
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
}

ANCOM.main <- function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}



#### multiple comparison ####
#### multiple comparison ####
#### multiple comparison ####
city_all<-unique(meta_ANCOM$city)
all_combn<-combn(city_all,2)
analysis_ANCOM<-as.data.frame(matrix(nrow=ncol(all_combn),ncol=ncol(data_for_ANCOM)-1))
colnames(analysis_ANCOM)<-colnames(data_for_ANCOM)[-ncol(data_for_ANCOM)]
analysis_ANCOM$comparison<-NA
for(i in 1:ncol(all_combn)){
  analysis_ANCOM[i,'comparison']<-paste0(all_combn[1,i],'-',all_combn[2,i])
  otu_test<-otu_ANCOM[substr(row.names(otu_ANCOM),1,3) %in% c(as.vector(all_combn[1,i]),as.vector(all_combn[2,i])),]
  meta_test<-meta_ANCOM[meta_ANCOM[,'city'] %in% c(as.vector(all_combn[1,i]),as.vector(all_combn[2,i])),]
  comparison_test=ANCOM.main(OTUdat=otu_test,
                             Vardat=meta_test,
                             adjusted=F,
                             repeated=F,
                             main.var="city",
                             adj.formula=NULL,
                             longitudinal=F,
                             repeat.var=NULL,
                             multcorr=2,
                             sig=0.05,
                             prev.cut=0.90)
  detected<-comparison_test$W.taxa
  analysis_ANCOM[i,-ncol(analysis_ANCOM)]<-0
  analysis_ANCOM[i,as.vector(detected[which(detected$detected_0.8==TRUE),'otu.names'])]<-1
}

row.names(analysis_ANCOM)<-analysis_ANCOM$comparison

#### For each city ####
analysis_ANCOM_city<-analysis_ANCOM
analysis_ANCOM_city$city_1<-NA
analysis_ANCOM_city$city_2<-NA
for(i in 1:nrow(analysis_ANCOM_city)){
  analysis_ANCOM_city$city_1[i]<-strsplit(analysis_ANCOM_city$comparison,split='-')[[i]][1]
  analysis_ANCOM_city$city_2[i]<-strsplit(analysis_ANCOM_city$comparison,split='-')[[i]][2]
}

ANCOM_city<-function(analysis_ANCOM_city,city,error_rate){
  ANCOM_tcity<-analysis_ANCOM_city[which(analysis_ANCOM_city$city_1==city|analysis_ANCOM_city$city_2==city),]
  ANCOM_tcity$city_2[which(ANCOM_tcity$city_2==city)]<-ANCOM_tcity$city_1[which(ANCOM_tcity$city_2==city)]
  ANCOM_tcity$city_1<-city
  ANCOM_tcity$comparison<-paste(ANCOM_tcity$city_1,ANCOM_tcity$city_2,sep='-')
  row.names(ANCOM_tcity)<-ANCOM_tcity$comparison
  plot_order_ANCOM<-ANCOM_tcity[,order(apply(ANCOM_tcity[,1:25],2,sum),decreasing = T)]
  ANCOM_plot<-Heatmap(t(plot_order_ANCOM),column_title=paste0('Error rate of ',city,': ',error_rate,'%'),
                      col = colorRamp2(c(0,0.5,1), brewer.pal(n=3, name="Blues")),
                      cluster_rows=F,
                      cluster_columns=F,
                      show_heatmap_legend=F)
  return(ANCOM_plot)
}

ANCOM_AKL<-ANCOM_city(analysis_ANCOM_city,'AKL',26.19)
ANCOM_BER<-ANCOM_city(analysis_ANCOM_city,'BER',12.69)
ANCOM_BOG<-ANCOM_city(analysis_ANCOM_city,'BOG',8.89)
ANCOM_HAM<-ANCOM_city(analysis_ANCOM_city,'HAM',16.67)
ANCOM_HGK<-ANCOM_city(analysis_ANCOM_city,'HGK',2.08)
ANCOM_ILR<-ANCOM_city(analysis_ANCOM_city,'ILR',5.56)
ANCOM_LON<-ANCOM_city(analysis_ANCOM_city,'LON',6.06)
ANCOM_MAR<-ANCOM_city(analysis_ANCOM_city,'MAR',13.33)
ANCOM_NYC<-ANCOM_city(analysis_ANCOM_city,'NYC',5.13)
ANCOM_OFA<-ANCOM_city(analysis_ANCOM_city,'OFA',8.33)
ANCOM_PXO<-ANCOM_city(analysis_ANCOM_city,'PXO',18.33)
ANCOM_SAC<-ANCOM_city(analysis_ANCOM_city,'SAC',3.70)
ANCOM_SAO<-ANCOM_city(analysis_ANCOM_city,'SAO',25.00)
ANCOM_SOF<-ANCOM_city(analysis_ANCOM_city,'SOF',33.33)
ANCOM_STO<-ANCOM_city(analysis_ANCOM_city,'STO',15.79)
ANCOM_TOK<-ANCOM_city(analysis_ANCOM_city,'TOK',6.67)


tiff(paste0( "ANCOM_AKL",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_AKL); dev.off()
tiff(paste0( "ANCOM_BER",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_BER); dev.off()
tiff(paste0( "ANCOM_BOG",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_BOG); dev.off()
tiff(paste0( "ANCOM_HAM",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_HAM); dev.off()
tiff(paste0( "ANCOM_HGK",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_HGK); dev.off()
tiff(paste0( "ANCOM_ILR",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_ILR); dev.off()
tiff(paste0( "ANCOM_LON",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_LON); dev.off()
tiff(paste0( "ANCOM_MAR",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_MAR); dev.off()
tiff(paste0( "ANCOM_NYC",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_NYC); dev.off()
tiff(paste0( "ANCOM_OFA",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_OFA); dev.off()
tiff(paste0( "ANCOM_PXO",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_PXO); dev.off()
tiff(paste0( "ANCOM_SAC",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_SAC); dev.off()
tiff(paste0( "ANCOM_SAO",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_SAO); dev.off()
tiff(paste0( "ANCOM_SOF",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_SOF); dev.off()
tiff(paste0( "ANCOM_STO",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_STO); dev.off()
tiff(paste0( "ANCOM_TOK",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_TOK); dev.off()

save.image(file='Main dataset 5 data for ANCOM.RData')

