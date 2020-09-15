load('otu_common_species_order.RData')
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
otu_common_species_order_mystery_normalized<-normalization_f(otu_common_species_order_mystery)
data_for_ANCOM<-otu_common_species_order_mystery_normalized

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
  plot_order_ANCOM<-ANCOM_tcity[,order(apply(ANCOM_tcity[,1:23],2,sum),decreasing = T)]
  ANCOM_plot<-Heatmap(t(plot_order_ANCOM),column_title=paste0('Error rate of ',city,': ',error_rate,'%'),
                      column_title_gp = gpar(fill = "white", col = "black", border = "white"),
                      col = colorRamp2(c(0,0.5,1), brewer.pal(n=3, name="Blues")),
                      cluster_rows=F,
                      cluster_columns=F,
                      show_heatmap_legend=F)
  return(ANCOM_plot)
}

ANCOM_Bri<-ANCOM_city(analysis_ANCOM_city,'Bri',61.11)
ANCOM_Doh<-ANCOM_city(analysis_ANCOM_city,'Doh',88.89)
ANCOM_Kie<-ANCOM_city(analysis_ANCOM_city,'Kie',72.22)
ANCOM_Osl<-ANCOM_city(analysis_ANCOM_city,'Osl',8.33)
ANCOM_Par<-ANCOM_city(analysis_ANCOM_city,'Par',16.67)
ANCOM_Rio<-ANCOM_city(analysis_ANCOM_city,'Rio',13.89)
ANCOM_San<-ANCOM_city(analysis_ANCOM_city,'San',16.67)
ANCOM_Vie<-ANCOM_city(analysis_ANCOM_city,'Vie',33.33)

tiff(paste0( "Mystery_","ANCOM_Bri",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_Bri); dev.off()
tiff(paste0( "Mystery_","ANCOM_Doh",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_Doh); dev.off()
tiff(paste0( "Mystery_","ANCOM_Kie",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_Kie); dev.off()
tiff(paste0( "Mystery_","ANCOM_Osl",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_Osl); dev.off()
tiff(paste0( "Mystery_","ANCOM_Par",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_Par); dev.off()
tiff(paste0( "Mystery_","ANCOM_Rio",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_Rio); dev.off()
tiff(paste0( "Mystery_","ANCOM_San",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_San); dev.off()
tiff(paste0( "Mystery_","ANCOM_Vie",".tif"), height=1400, width=1400, res=300, units="px"); print(ANCOM_Vie); dev.off()

save.image('Mystery dataset 5 data for ANCOM.RData')





