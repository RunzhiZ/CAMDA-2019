load('Main dataset 1 OTU tables for all samples.RData')
############################################## Obtain the common species, families and orders ###############################
l_species<-list(as.vector(AKL_species),
                as.vector(BER_species),
                as.vector(BOG_species),
                as.vector(HAM_species),
                as.vector(HGK_species),
                as.vector(ILR_species),
                as.vector(LON_species),
                as.vector(MAR_species),
                as.vector(NYC_species),
                as.vector(OFA_species),
                as.vector(PXO_species),
                as.vector(SAC_species),
                as.vector(SAO_species),
                as.vector(SOF_species),
                as.vector(STO_species),
                as.vector(TOK_species))

l_family<-list(as.vector(AKL_family),
               as.vector(BER_family),
               as.vector(BOG_family),
               as.vector(HAM_family),
               as.vector(HGK_family),
               as.vector(ILR_family),
               as.vector(LON_family),
               as.vector(MAR_family),
               as.vector(NYC_family),
               as.vector(OFA_family),
               as.vector(PXO_family),
               as.vector(SAC_family),
               as.vector(SAO_family),
               as.vector(SOF_family),
               as.vector(STO_family),
               as.vector(TOK_family))

l_order<-list(as.vector(AKL_order),
              as.vector(BER_order),
              as.vector(BOG_order),
              as.vector(HAM_order),
              as.vector(HGK_order),
              as.vector(ILR_order),
              as.vector(LON_order),
              as.vector(MAR_order),
              as.vector(NYC_order),
              as.vector(OFA_order),
              as.vector(PXO_order),
              as.vector(SAC_order),
              as.vector(SAO_order),
              as.vector(SOF_order),
              as.vector(STO_order),
              as.vector(TOK_order))
  
common_species<-Reduce(intersect,l_species) ## 7
common_family<-Reduce(intersect,l_family) ## 9
common_order<-Reduce(intersect,l_order) ## 9
all_species<-Reduce(union,l_species)
all_family<-Reduce(union,l_family)
all_order<-Reduce(union,l_order)

##################### calculate the number of feature #################
city_all<-c('AKL','BER','BOG','HAM','HGK','ILR','LON','MAR','NYC','OFA','PXO','SAC','SAO','SOF','STO','TOK')
city_species<-data.frame()
for(i in city_all){
  city_species[which(city_all==i),1]<-i
  city_species[which(city_all==i),2]<-eval(parse(text=paste0('length(',i,'_species)')))
}
colnames(city_species)<-c('city','number_species')

city_family<-data.frame()
for(i in city_all){
  city_family[which(city_all==i),1]<-i
  city_family[which(city_all==i),2]<-eval(parse(text=paste0('length(',i,'_family)')))
}
colnames(city_family)<-c('city','number_family')

city_order<-data.frame()
for(i in city_all){
  city_order[which(city_all==i),1]<-i
  city_order[which(city_all==i),2]<-eval(parse(text=paste0('length(',i,'_order)')))
}
colnames(city_order)<-c('city','number_order')

city_feature<-cbind(city_species,city_family$number_family,city_order$number_order)
colnames(city_feature)[c(3,4)]<-c('number_family','number_order')

################################################################### Species ##########################################################################
##
##
##
##
#### AKL
otu_species_AKL<-as.data.frame(matrix(nrow=length(AKL_species),ncol=14))
row.names(otu_species_AKL)<-AKL_species
for(i in 1:14){
  colnames(otu_species_AKL)[i]<-paste0('AKL_',i)
  aa<-eval(parse(text=paste0('AKL_',i)))
  for(j in AKL_species){
    if(length(which(aa$Species %in% j))==0){otu_species_AKL[j,i]<-0}
    else{otu_species_AKL[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_AKL$ubiquity<-0
for(i in 1:nrow(otu_species_AKL)){
  otu_species_AKL[i,'ubiquity']<-length(which(otu_species_AKL[i,-ncol(otu_species_AKL)]!=0))/(ncol(otu_species_AKL)-1)
}

#### BER
otu_species_BER<-as.data.frame(matrix(nrow=length(BER_species),ncol=21))
row.names(otu_species_BER)<-BER_species
for(i in 1:21){
  colnames(otu_species_BER)[i]<-paste0('BER_',i)
  aa<-eval(parse(text=paste0('BER_',i)))
  for(j in BER_species){
    if(length(which(aa$Species %in% j))==0){otu_species_BER[j,i]<-0}
    else{otu_species_BER[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_BER$ubiquity<-0
for(i in 1:nrow(otu_species_BER)){
  otu_species_BER[i,'ubiquity']<-length(which(otu_species_BER[i,-ncol(otu_species_BER)]!=0))/(ncol(otu_species_BER)-1)
}

#### BOG
otu_species_BOG<-as.data.frame(matrix(nrow=length(BOG_species),ncol=15))
row.names(otu_species_BOG)<-BOG_species
for(i in 1:ncol(otu_species_BOG)){
  colnames(otu_species_BOG)[i]<-paste0('BOG_',i)
  aa<-eval(parse(text=paste0('BOG_',i)))
  for(j in BOG_species){
    if(length(which(aa$Species %in% j))==0){otu_species_BOG[j,i]<-0}
    else{otu_species_BOG[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_BOG$ubiquity<-0
for(i in 1:nrow(otu_species_BOG)){
  otu_species_BOG[i,'ubiquity']<-length(which(otu_species_BOG[i,-ncol(otu_species_BOG)]!=0))/(ncol(otu_species_BOG)-1)
}

#### HAM
otu_species_HAM<-as.data.frame(matrix(nrow=length(HAM_species),ncol=16))
row.names(otu_species_HAM)<-HAM_species
for(i in 1:ncol(otu_species_HAM)){
  colnames(otu_species_HAM)[i]<-paste0('HAM_',i)
  aa<-eval(parse(text=paste0('HAM_',i)))
  for(j in HAM_species){
    if(length(which(aa$Species %in% j))==0){otu_species_HAM[j,i]<-0}
    else{otu_species_HAM[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_HAM$ubiquity<-0
for(i in 1:nrow(otu_species_HAM)){
  otu_species_HAM[i,'ubiquity']<-length(which(otu_species_HAM[i,-ncol(otu_species_HAM)]!=0))/(ncol(otu_species_HAM)-1)
}

#### HGK
otu_species_HGK<-as.data.frame(matrix(nrow=length(HGK_species),ncol=18))
row.names(otu_species_HGK)<-HGK_species
for(i in c(1:8,10:18)){
  colnames(otu_species_HGK)[i]<-paste0('HGK_',i)
  aa<-eval(parse(text=paste0('HGK_',i)))
  for(j in HGK_species){
    if(length(which(aa$Species %in% j))==0){otu_species_HGK[j,i]<-0}
    else{otu_species_HGK[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_HGK$V9<-NULL
otu_species_HGK$ubiquity<-0
for(i in 1:nrow(otu_species_HGK)){
  otu_species_HGK[i,'ubiquity']<-length(which(otu_species_HGK[i,-ncol(otu_species_HGK)]!=0))/(ncol(otu_species_HGK)-1)
}

#### ILR
otu_species_ILR<-as.data.frame(matrix(nrow=length(ILR_species),ncol=24))
row.names(otu_species_ILR)<-ILR_species
for(i in 1:ncol(otu_species_ILR)){
  colnames(otu_species_ILR)[i]<-paste0('ILR_',i)
  aa<-eval(parse(text=paste0('ILR_',i)))
  for(j in ILR_species){
    if(length(which(aa$Species %in% j))==0){otu_species_ILR[j,i]<-0}
    else{otu_species_ILR[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_ILR$ubiquity<-0
for(i in 1:nrow(otu_species_ILR)){
  otu_species_ILR[i,'ubiquity']<-length(which(otu_species_ILR[i,-ncol(otu_species_ILR)]!=0))/(ncol(otu_species_ILR)-1)
}

#### LON
otu_species_LON<-as.data.frame(matrix(nrow=length(LON_species),ncol=24))
row.names(otu_species_LON)<-LON_species
for(i in c(1,3:6,8:24)){
  colnames(otu_species_LON)[i]<-paste0('LON_',i)
  aa<-eval(parse(text=paste0('LON_',i)))
  for(j in LON_species){
    if(length(which(aa$Species %in% j))==0){otu_species_LON[j,i]<-0}
    else{otu_species_LON[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_LON$V2<-NULL
otu_species_LON$V7<-NULL
otu_species_LON$ubiquity<-0
for(i in 1:nrow(otu_species_LON)){
  otu_species_LON[i,'ubiquity']<-length(which(otu_species_LON[i,-ncol(otu_species_LON)]!=0))/(ncol(otu_species_LON)-1)
}

#### MAR
otu_species_MAR<-as.data.frame(matrix(nrow=length(MAR_species),ncol=10))
row.names(otu_species_MAR)<-MAR_species
for(i in 1:ncol(otu_species_MAR)){
  colnames(otu_species_MAR)[i]<-paste0('MAR_',i)
  aa<-eval(parse(text=paste0('MAR_',i)))
  for(j in MAR_species){
    if(length(which(aa$Species %in% j))==0){otu_species_MAR[j,i]<-0}
    else{otu_species_MAR[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_MAR$ubiquity<-0
for(i in 1:nrow(otu_species_MAR)){
  otu_species_MAR[i,'ubiquity']<-length(which(otu_species_MAR[i,-ncol(otu_species_MAR)]!=0))/(ncol(otu_species_MAR)-1)
}

#### NYC
otu_species_NYC<-as.data.frame(matrix(nrow=length(NYC_species),ncol=26))
row.names(otu_species_NYC)<-NYC_species
for(i in 1:ncol(otu_species_NYC)){
  colnames(otu_species_NYC)[i]<-paste0('NYC_',i)
  aa<-eval(parse(text=paste0('NYC_',i)))
  for(j in NYC_species){
    if(length(which(aa$Species %in% j))==0){otu_species_NYC[j,i]<-0}
    else{otu_species_NYC[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_NYC$ubiquity<-0
for(i in 1:nrow(otu_species_NYC)){
  otu_species_NYC[i,'ubiquity']<-length(which(otu_species_NYC[i,-ncol(otu_species_NYC)]!=0))/(ncol(otu_species_NYC)-1)
}

#### OFA
otu_species_OFA<-as.data.frame(matrix(nrow=length(OFA_species),ncol=20))
row.names(otu_species_OFA)<-OFA_species
for(i in 1:ncol(otu_species_OFA)){
  colnames(otu_species_OFA)[i]<-paste0('OFA_',i)
  aa<-eval(parse(text=paste0('OFA_',i)))
  for(j in OFA_species){
    if(length(which(aa$Species %in% j))==0){otu_species_OFA[j,i]<-0}
    else{otu_species_OFA[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_OFA$ubiquity<-0
for(i in 1:nrow(otu_species_OFA)){
  otu_species_OFA[i,'ubiquity']<-length(which(otu_species_OFA[i,-ncol(otu_species_OFA)]!=0))/(ncol(otu_species_OFA)-1)
}

#### PXO
otu_species_PXO<-as.data.frame(matrix(nrow=length(PXO_species),ncol=20))
row.names(otu_species_PXO)<-PXO_species
for(i in 1:ncol(otu_species_PXO)){
  colnames(otu_species_PXO)[i]<-paste0('PXO_',i)
  aa<-eval(parse(text=paste0('PXO_',i)))
  for(j in PXO_species){
    if(length(which(aa$Species %in% j))==0){otu_species_PXO[j,i]<-0}
    else{otu_species_PXO[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_PXO$ubiquity<-0
for(i in 1:nrow(otu_species_PXO)){
  otu_species_PXO[i,'ubiquity']<-length(which(otu_species_PXO[i,-ncol(otu_species_PXO)]!=0))/(ncol(otu_species_PXO)-1)
}

#### SAC
otu_species_SAC<-as.data.frame(matrix(nrow=length(SAC_species),ncol=18))
row.names(otu_species_SAC)<-SAC_species
for(i in 1:ncol(otu_species_SAC)){
  colnames(otu_species_SAC)[i]<-paste0('SAC_',i)
  aa<-eval(parse(text=paste0('SAC_',i)))
  for(j in SAC_species){
    if(length(which(aa$Species %in% j))==0){otu_species_SAC[j,i]<-0}
    else{otu_species_SAC[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_SAC$ubiquity<-0
for(i in 1:nrow(otu_species_SAC)){
  otu_species_SAC[i,'ubiquity']<-length(which(otu_species_SAC[i,-ncol(otu_species_SAC)]!=0))/(ncol(otu_species_SAC)-1)
}

#### SAO
otu_species_SAO<-as.data.frame(matrix(nrow=length(SAO_species),ncol=24))
row.names(otu_species_SAO)<-SAO_species
for(i in 1:ncol(otu_species_SAO)){
  colnames(otu_species_SAO)[i]<-paste0('SAO_',i)
  aa<-eval(parse(text=paste0('SAO_',i)))
  for(j in SAO_species){
    if(length(which(aa$Species %in% j))==0){otu_species_SAO[j,i]<-0}
    else{otu_species_SAO[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_SAO$ubiquity<-0
for(i in 1:nrow(otu_species_SAO)){
  otu_species_SAO[i,'ubiquity']<-length(which(otu_species_SAO[i,-ncol(otu_species_SAO)]!=0))/(ncol(otu_species_SAO)-1)
}

#### SOF
otu_species_SOF<-as.data.frame(matrix(nrow=length(SOF_species),ncol=10))
row.names(otu_species_SOF)<-SOF_species
for(i in 1:ncol(otu_species_SOF)){
  colnames(otu_species_SOF)[i]<-paste0('SOF_',i)
  aa<-eval(parse(text=paste0('SOF_',i)))
  for(j in SOF_species){
    if(length(which(aa$Species %in% j))==0){otu_species_SOF[j,i]<-0}
    else{otu_species_SOF[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_SOF$ubiquity<-0
for(i in 1:nrow(otu_species_SOF)){
  otu_species_SOF[i,'ubiquity']<-length(which(otu_species_SOF[i,-ncol(otu_species_SOF)]!=0))/(ncol(otu_species_SOF)-1)
}

#### STO
otu_species_STO<-as.data.frame(matrix(nrow=length(STO_species),ncol=20))
row.names(otu_species_STO)<-STO_species
for(i in 1:ncol(otu_species_STO)){
  colnames(otu_species_STO)[i]<-paste0('STO_',i)
  aa<-eval(parse(text=paste0('STO_',i)))
  for(j in STO_species){
    if(length(which(aa$Species %in% j))==0){otu_species_STO[j,i]<-0}
    else{otu_species_STO[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_STO$ubiquity<-0
for(i in 1:nrow(otu_species_STO)){
  otu_species_STO[i,'ubiquity']<-length(which(otu_species_STO[i,-ncol(otu_species_STO)]!=0))/(ncol(otu_species_STO)-1)
}

#### TOK
otu_species_TOK<-as.data.frame(matrix(nrow=length(TOK_species),ncol=25))
row.names(otu_species_TOK)<-TOK_species
for(i in 1:ncol(otu_species_TOK)){
  colnames(otu_species_TOK)[i]<-paste0('TOK_',i)
  aa<-eval(parse(text=paste0('TOK_',i)))
  for(j in TOK_species){
    if(length(which(aa$Species %in% j))==0){otu_species_TOK[j,i]<-0}
    else{otu_species_TOK[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_TOK$ubiquity<-0
for(i in 1:nrow(otu_species_TOK)){
  otu_species_TOK[i,'ubiquity']<-length(which(otu_species_TOK[i,-ncol(otu_species_TOK)]!=0))/(ncol(otu_species_TOK)-1)
}

#### check ubiquity
check_ubiquity<-function(otu,common){
  aa<-otu[which(row.names(otu) %in% common),]
  return(aa)
}
check_ubiquity(otu_species_AKL,common_species)
check_ubiquity(otu_species_BER,common_species)
check_ubiquity(otu_species_BOG,common_species)
check_ubiquity(otu_species_HAM,common_species)
check_ubiquity(otu_species_HGK,common_species)
check_ubiquity(otu_species_ILR,common_species)
check_ubiquity(otu_species_LON,common_species)
check_ubiquity(otu_species_MAR,common_species)
check_ubiquity(otu_species_NYC,common_species)
check_ubiquity(otu_species_OFA,common_species)
check_ubiquity(otu_species_PXO,common_species)
check_ubiquity(otu_species_SAC,common_species)
check_ubiquity(otu_species_SAO,common_species)
check_ubiquity(otu_species_SOF,common_species)
check_ubiquity(otu_species_STO,common_species)
check_ubiquity(otu_species_TOK,common_species)

rm(aa,i,j)

## ubiquity
otu_species_ubiquity<-as.data.frame(matrix(nrow=length(common_species),ncol=16))
row.names(otu_species_ubiquity)<-common_species
colnames(otu_species_ubiquity)<-c('AKL','BER','BOG','HAM',
                                  'HGK','ILR','LON','MAR',
                                  'NYC','OFA','PXO','SAC',
                                  'SAO','SOF','STO','TOK')
for(i in c('AKL','BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
  eval(parse(text=paste0('otu_species_ubiquity$',i,'<-','check_ubiquity(otu_species_',i,',common_species)$ubiquity')))
}

## common species cbind - all samples
otu_common_species<-otu_species_AKL[which(row.names(otu_species_AKL) %in% common_species),-ncol(otu_species_AKL)]
for(i in c('BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
aa<-eval(parse(text=paste0('otu_species_',i)))
bb<-aa[which(row.names(aa) %in% common_species),-ncol(aa)]
otu_common_species<-cbind(otu_common_species,bb)
}

#### all species
otu_all_species<-as.data.frame(matrix(nrow=length(all_species),ncol=302))
row.names(otu_all_species)<-all_species
colnames(otu_all_species)<-c(paste0('AKL','_',c(1:14)),
                             paste0('BER','_',c(1:21)),
                             paste0('BOG','_',c(1:15)),
                             paste0('HAM','_',c(1:16)),
                             paste0('HGK','_',c(1:8,10:18)),
                             paste0('ILR','_',c(1:24)),
                             paste0('LON','_',c(1,3:6,8:24)),
                             paste0('MAR','_',c(1:10)),
                             paste0('NYC','_',c(1:26)),
                             paste0('OFA','_',c(1:20)),
                             paste0('PXO','_',c(1:20)),
                             paste0('SAC','_',c(1:18)),
                             paste0('SAO','_',c(1:24)),
                             paste0('SOF','_',c(1:10)),
                             paste0('STO','_',c(1:20)),
                             paste0('TOK','_',c(1:25)))

for(i in colnames(otu_all_species)){
  aa<-eval(parse(text=i))
  for(j in all_species){
    if(length(which(aa$Species %in% j))==0){otu_all_species[j,i]<-0}
    else{otu_all_species[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_all_species$ubiquity<-NA
for(i in 1:nrow(otu_all_species)){
  otu_all_species[i,'ubiquity']<-(302-length(which(otu_all_species[i,1:302]==0)))/302
}
otu_all_species[common_species,'ubiquity']

################################################################# family ###########################################################################
##
##
##
##
#### AKL
otu_family_AKL<-as.data.frame(matrix(nrow=length(AKL_family),ncol=14))
row.names(otu_family_AKL)<-AKL_family
for(i in 1:14){
  colnames(otu_family_AKL)[i]<-paste0('AKL_',i)
  aa<-eval(parse(text=paste0('AKL_',i)))
  for(j in AKL_family){
    if(length(which(aa$Family %in% j))==0){otu_family_AKL[j,i]<-0}
    else{otu_family_AKL[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_AKL$ubiquity<-0
for(i in 1:nrow(otu_family_AKL)){
  otu_family_AKL[i,'ubiquity']<-length(which(otu_family_AKL[i,-ncol(otu_family_AKL)]!=0))/(ncol(otu_family_AKL)-1)
}

#### BER
otu_family_BER<-as.data.frame(matrix(nrow=length(BER_family),ncol=21))
row.names(otu_family_BER)<-BER_family
for(i in 1:21){
  colnames(otu_family_BER)[i]<-paste0('BER_',i)
  aa<-eval(parse(text=paste0('BER_',i)))
  for(j in BER_family){
    if(length(which(aa$Family %in% j))==0){otu_family_BER[j,i]<-0}
    else{otu_family_BER[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_BER$ubiquity<-0
for(i in 1:nrow(otu_family_BER)){
  otu_family_BER[i,'ubiquity']<-length(which(otu_family_BER[i,-ncol(otu_family_BER)]!=0))/(ncol(otu_family_BER)-1)
}

#### BOG
otu_family_BOG<-as.data.frame(matrix(nrow=length(BOG_family),ncol=15))
row.names(otu_family_BOG)<-BOG_family
for(i in 1:ncol(otu_family_BOG)){
  colnames(otu_family_BOG)[i]<-paste0('BOG_',i)
  aa<-eval(parse(text=paste0('BOG_',i)))
  for(j in BOG_family){
    if(length(which(aa$Family %in% j))==0){otu_family_BOG[j,i]<-0}
    else{otu_family_BOG[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_BOG$ubiquity<-0
for(i in 1:nrow(otu_family_BOG)){
  otu_family_BOG[i,'ubiquity']<-length(which(otu_family_BOG[i,-ncol(otu_family_BOG)]!=0))/(ncol(otu_family_BOG)-1)
}

#### HAM
otu_family_HAM<-as.data.frame(matrix(nrow=length(HAM_family),ncol=16))
row.names(otu_family_HAM)<-HAM_family
for(i in 1:ncol(otu_family_HAM)){
  colnames(otu_family_HAM)[i]<-paste0('HAM_',i)
  aa<-eval(parse(text=paste0('HAM_',i)))
  for(j in HAM_family){
    if(length(which(aa$Family %in% j))==0){otu_family_HAM[j,i]<-0}
    else{otu_family_HAM[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_HAM$ubiquity<-0
for(i in 1:nrow(otu_family_HAM)){
  otu_family_HAM[i,'ubiquity']<-length(which(otu_family_HAM[i,-ncol(otu_family_HAM)]!=0))/(ncol(otu_family_HAM)-1)
}

#### HGK
otu_family_HGK<-as.data.frame(matrix(nrow=length(HGK_family),ncol=18))
row.names(otu_family_HGK)<-HGK_family
for(i in c(1:8,10:18)){
  colnames(otu_family_HGK)[i]<-paste0('HGK_',i)
  aa<-eval(parse(text=paste0('HGK_',i)))
  for(j in HGK_family){
    if(length(which(aa$Family %in% j))==0){otu_family_HGK[j,i]<-0}
    else{otu_family_HGK[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_HGK$V9<-NULL
otu_family_HGK$ubiquity<-0
for(i in 1:nrow(otu_family_HGK)){
  otu_family_HGK[i,'ubiquity']<-length(which(otu_family_HGK[i,-ncol(otu_family_HGK)]!=0))/(ncol(otu_family_HGK)-1)
}

#### ILR
otu_family_ILR<-as.data.frame(matrix(nrow=length(ILR_family),ncol=24))
row.names(otu_family_ILR)<-ILR_family
for(i in 1:ncol(otu_family_ILR)){
  colnames(otu_family_ILR)[i]<-paste0('ILR_',i)
  aa<-eval(parse(text=paste0('ILR_',i)))
  for(j in ILR_family){
    if(length(which(aa$Family %in% j))==0){otu_family_ILR[j,i]<-0}
    else{otu_family_ILR[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_ILR$ubiquity<-0
for(i in 1:nrow(otu_family_ILR)){
  otu_family_ILR[i,'ubiquity']<-length(which(otu_family_ILR[i,-ncol(otu_family_ILR)]!=0))/(ncol(otu_family_ILR)-1)
}

#### LON
otu_family_LON<-as.data.frame(matrix(nrow=length(LON_family),ncol=24))
row.names(otu_family_LON)<-LON_family
for(i in c(1,3:6,8:24)){
  colnames(otu_family_LON)[i]<-paste0('LON_',i)
  aa<-eval(parse(text=paste0('LON_',i)))
  for(j in LON_family){
    if(length(which(aa$Family %in% j))==0){otu_family_LON[j,i]<-0}
    else{otu_family_LON[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_LON$V2<-NULL
otu_family_LON$V7<-NULL
otu_family_LON$ubiquity<-0
for(i in 1:nrow(otu_family_LON)){
  otu_family_LON[i,'ubiquity']<-length(which(otu_family_LON[i,-ncol(otu_family_LON)]!=0))/(ncol(otu_family_LON)-1)
}

#### MAR
otu_family_MAR<-as.data.frame(matrix(nrow=length(MAR_family),ncol=10))
row.names(otu_family_MAR)<-MAR_family
for(i in 1:ncol(otu_family_MAR)){
  colnames(otu_family_MAR)[i]<-paste0('MAR_',i)
  aa<-eval(parse(text=paste0('MAR_',i)))
  for(j in MAR_family){
    if(length(which(aa$Family %in% j))==0){otu_family_MAR[j,i]<-0}
    else{otu_family_MAR[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_MAR$ubiquity<-0
for(i in 1:nrow(otu_family_MAR)){
  otu_family_MAR[i,'ubiquity']<-length(which(otu_family_MAR[i,-ncol(otu_family_MAR)]!=0))/(ncol(otu_family_MAR)-1)
}

#### NYC
otu_family_NYC<-as.data.frame(matrix(nrow=length(NYC_family),ncol=26))
row.names(otu_family_NYC)<-NYC_family
for(i in 1:ncol(otu_family_NYC)){
  colnames(otu_family_NYC)[i]<-paste0('NYC_',i)
  aa<-eval(parse(text=paste0('NYC_',i)))
  for(j in NYC_family){
    if(length(which(aa$Family %in% j))==0){otu_family_NYC[j,i]<-0}
    else{otu_family_NYC[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_NYC$ubiquity<-0
for(i in 1:nrow(otu_family_NYC)){
  otu_family_NYC[i,'ubiquity']<-length(which(otu_family_NYC[i,-ncol(otu_family_NYC)]!=0))/(ncol(otu_family_NYC)-1)
}

#### OFA
otu_family_OFA<-as.data.frame(matrix(nrow=length(OFA_family),ncol=20))
row.names(otu_family_OFA)<-OFA_family
for(i in 1:ncol(otu_family_OFA)){
  colnames(otu_family_OFA)[i]<-paste0('OFA_',i)
  aa<-eval(parse(text=paste0('OFA_',i)))
  for(j in OFA_family){
    if(length(which(aa$Family %in% j))==0){otu_family_OFA[j,i]<-0}
    else{otu_family_OFA[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_OFA$ubiquity<-0
for(i in 1:nrow(otu_family_OFA)){
  otu_family_OFA[i,'ubiquity']<-length(which(otu_family_OFA[i,-ncol(otu_family_OFA)]!=0))/(ncol(otu_family_OFA)-1)
}

#### PXO
otu_family_PXO<-as.data.frame(matrix(nrow=length(PXO_family),ncol=20))
row.names(otu_family_PXO)<-PXO_family
for(i in 1:ncol(otu_family_PXO)){
  colnames(otu_family_PXO)[i]<-paste0('PXO_',i)
  aa<-eval(parse(text=paste0('PXO_',i)))
  for(j in PXO_family){
    if(length(which(aa$Family %in% j))==0){otu_family_PXO[j,i]<-0}
    else{otu_family_PXO[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_PXO$ubiquity<-0
for(i in 1:nrow(otu_family_PXO)){
  otu_family_PXO[i,'ubiquity']<-length(which(otu_family_PXO[i,-ncol(otu_family_PXO)]!=0))/(ncol(otu_family_PXO)-1)
}

#### SAC
otu_family_SAC<-as.data.frame(matrix(nrow=length(SAC_family),ncol=18))
row.names(otu_family_SAC)<-SAC_family
for(i in 1:ncol(otu_family_SAC)){
  colnames(otu_family_SAC)[i]<-paste0('SAC_',i)
  aa<-eval(parse(text=paste0('SAC_',i)))
  for(j in SAC_family){
    if(length(which(aa$Family %in% j))==0){otu_family_SAC[j,i]<-0}
    else{otu_family_SAC[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_SAC$ubiquity<-0
for(i in 1:nrow(otu_family_SAC)){
  otu_family_SAC[i,'ubiquity']<-length(which(otu_family_SAC[i,-ncol(otu_family_SAC)]!=0))/(ncol(otu_family_SAC)-1)
}

#### SAO
otu_family_SAO<-as.data.frame(matrix(nrow=length(SAO_family),ncol=24))
row.names(otu_family_SAO)<-SAO_family
for(i in 1:ncol(otu_family_SAO)){
  colnames(otu_family_SAO)[i]<-paste0('SAO_',i)
  aa<-eval(parse(text=paste0('SAO_',i)))
  for(j in SAO_family){
    if(length(which(aa$Family %in% j))==0){otu_family_SAO[j,i]<-0}
    else{otu_family_SAO[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_SAO$ubiquity<-0
for(i in 1:nrow(otu_family_SAO)){
  otu_family_SAO[i,'ubiquity']<-length(which(otu_family_SAO[i,-ncol(otu_family_SAO)]!=0))/(ncol(otu_family_SAO)-1)
}

#### SOF
otu_family_SOF<-as.data.frame(matrix(nrow=length(SOF_family),ncol=10))
row.names(otu_family_SOF)<-SOF_family
for(i in 1:ncol(otu_family_SOF)){
  colnames(otu_family_SOF)[i]<-paste0('SOF_',i)
  aa<-eval(parse(text=paste0('SOF_',i)))
  for(j in SOF_family){
    if(length(which(aa$Family %in% j))==0){otu_family_SOF[j,i]<-0}
    else{otu_family_SOF[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_SOF$ubiquity<-0
for(i in 1:nrow(otu_family_SOF)){
  otu_family_SOF[i,'ubiquity']<-length(which(otu_family_SOF[i,-ncol(otu_family_SOF)]!=0))/(ncol(otu_family_SOF)-1)
}

#### STO
otu_family_STO<-as.data.frame(matrix(nrow=length(STO_family),ncol=20))
row.names(otu_family_STO)<-STO_family
for(i in 1:ncol(otu_family_STO)){
  colnames(otu_family_STO)[i]<-paste0('STO_',i)
  aa<-eval(parse(text=paste0('STO_',i)))
  for(j in STO_family){
    if(length(which(aa$Family %in% j))==0){otu_family_STO[j,i]<-0}
    else{otu_family_STO[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_STO$ubiquity<-0
for(i in 1:nrow(otu_family_STO)){
  otu_family_STO[i,'ubiquity']<-length(which(otu_family_STO[i,-ncol(otu_family_STO)]!=0))/(ncol(otu_family_STO)-1)
}

#### TOK
otu_family_TOK<-as.data.frame(matrix(nrow=length(TOK_family),ncol=25))
row.names(otu_family_TOK)<-TOK_family
for(i in 1:ncol(otu_family_TOK)){
  colnames(otu_family_TOK)[i]<-paste0('TOK_',i)
  aa<-eval(parse(text=paste0('TOK_',i)))
  for(j in TOK_family){
    if(length(which(aa$Family %in% j))==0){otu_family_TOK[j,i]<-0}
    else{otu_family_TOK[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_TOK$ubiquity<-0
for(i in 1:nrow(otu_family_TOK)){
  otu_family_TOK[i,'ubiquity']<-length(which(otu_family_TOK[i,-ncol(otu_family_TOK)]!=0))/(ncol(otu_family_TOK)-1)
}

#### check ubiquity
check_ubiquity(otu_family_AKL,common_family)
check_ubiquity(otu_family_BER,common_family)
check_ubiquity(otu_family_BOG,common_family)
check_ubiquity(otu_family_HAM,common_family)
check_ubiquity(otu_family_HGK,common_family)
check_ubiquity(otu_family_ILR,common_family)
check_ubiquity(otu_family_LON,common_family)
check_ubiquity(otu_family_MAR,common_family)
check_ubiquity(otu_family_NYC,common_family)
check_ubiquity(otu_family_OFA,common_family)
check_ubiquity(otu_family_PXO,common_family)
check_ubiquity(otu_family_SAC,common_family)
check_ubiquity(otu_family_SAO,common_family)
check_ubiquity(otu_family_SOF,common_family)
check_ubiquity(otu_family_STO,common_family)
check_ubiquity(otu_family_TOK,common_family)

rm(aa,i,j)

## ubiquity
otu_family_ubiquity<-as.data.frame(matrix(nrow=length(common_family),ncol=16))
row.names(otu_family_ubiquity)<-common_family
colnames(otu_family_ubiquity)<-c('AKL','BER','BOG','HAM',
                                  'HGK','ILR','LON','MAR',
                                  'NYC','OFA','PXO','SAC',
                                  'SAO','SOF','STO','TOK')
for(i in c('AKL','BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
  eval(parse(text=paste0('otu_family_ubiquity$',i,'<-','check_ubiquity(otu_family_',i,',common_family)$ubiquity')))
}

## common family cbind - all samples
otu_common_family<-otu_family_AKL[which(row.names(otu_family_AKL) %in% common_family),-ncol(otu_family_AKL)]
for(i in c('BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
  aa<-eval(parse(text=paste0('otu_family_',i)))
  bb<-aa[which(row.names(aa) %in% common_family),-ncol(aa)]
  otu_common_family<-cbind(otu_common_family,bb)
}


## all family
otu_all_family<-as.data.frame(matrix(nrow=length(all_family),ncol=302))
row.names(otu_all_family)<-all_family
colnames(otu_all_family)<-c(paste0('AKL','_',c(1:14)),
                             paste0('BER','_',c(1:21)),
                             paste0('BOG','_',c(1:15)),
                             paste0('HAM','_',c(1:16)),
                             paste0('HGK','_',c(1:8,10:18)),
                             paste0('ILR','_',c(1:24)),
                             paste0('LON','_',c(1,3:6,8:24)),
                             paste0('MAR','_',c(1:10)),
                             paste0('NYC','_',c(1:26)),
                             paste0('OFA','_',c(1:20)),
                             paste0('PXO','_',c(1:20)),
                             paste0('SAC','_',c(1:18)),
                             paste0('SAO','_',c(1:24)),
                             paste0('SOF','_',c(1:10)),
                             paste0('STO','_',c(1:20)),
                             paste0('TOK','_',c(1:25)))

for(i in colnames(otu_all_family)){
  aa<-eval(parse(text=i))
  for(j in all_family){
    if(length(which(aa$Family %in% j))==0){otu_all_family[j,i]<-0}
    else{otu_all_family[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}

otu_all_family$ubiquity<-NA
for(i in 1:nrow(otu_all_family)){
  otu_all_family[i,'ubiquity']<-(302-length(which(otu_all_family[i,1:302]==0)))/302
}
otu_all_family[common_family,'ubiquity']


################################################################# order ###########################################################################
##
##
##
##
#### AKL
otu_order_AKL<-as.data.frame(matrix(nrow=length(AKL_order),ncol=14))
row.names(otu_order_AKL)<-AKL_order
for(i in 1:14){
  colnames(otu_order_AKL)[i]<-paste0('AKL_',i)
  aa<-eval(parse(text=paste0('AKL_',i)))
  for(j in AKL_order){
    if(length(which(aa$Order %in% j))==0){otu_order_AKL[j,i]<-0}
    else{otu_order_AKL[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_AKL$ubiquity<-0
for(i in 1:nrow(otu_order_AKL)){
  otu_order_AKL[i,'ubiquity']<-length(which(otu_order_AKL[i,-ncol(otu_order_AKL)]!=0))/(ncol(otu_order_AKL)-1)
}

#### BER
otu_order_BER<-as.data.frame(matrix(nrow=length(BER_order),ncol=21))
row.names(otu_order_BER)<-BER_order
for(i in 1:21){
  colnames(otu_order_BER)[i]<-paste0('BER_',i)
  aa<-eval(parse(text=paste0('BER_',i)))
  for(j in BER_order){
    if(length(which(aa$Order %in% j))==0){otu_order_BER[j,i]<-0}
    else{otu_order_BER[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_BER$ubiquity<-0
for(i in 1:nrow(otu_order_BER)){
  otu_order_BER[i,'ubiquity']<-length(which(otu_order_BER[i,-ncol(otu_order_BER)]!=0))/(ncol(otu_order_BER)-1)
}

#### BOG
otu_order_BOG<-as.data.frame(matrix(nrow=length(BOG_order),ncol=15))
row.names(otu_order_BOG)<-BOG_order
for(i in 1:ncol(otu_order_BOG)){
  colnames(otu_order_BOG)[i]<-paste0('BOG_',i)
  aa<-eval(parse(text=paste0('BOG_',i)))
  for(j in BOG_order){
    if(length(which(aa$Order %in% j))==0){otu_order_BOG[j,i]<-0}
    else{otu_order_BOG[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_BOG$ubiquity<-0
for(i in 1:nrow(otu_order_BOG)){
  otu_order_BOG[i,'ubiquity']<-length(which(otu_order_BOG[i,-ncol(otu_order_BOG)]!=0))/(ncol(otu_order_BOG)-1)
}

#### HAM
otu_order_HAM<-as.data.frame(matrix(nrow=length(HAM_order),ncol=16))
row.names(otu_order_HAM)<-HAM_order
for(i in 1:ncol(otu_order_HAM)){
  colnames(otu_order_HAM)[i]<-paste0('HAM_',i)
  aa<-eval(parse(text=paste0('HAM_',i)))
  for(j in HAM_order){
    if(length(which(aa$Order %in% j))==0){otu_order_HAM[j,i]<-0}
    else{otu_order_HAM[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_HAM$ubiquity<-0
for(i in 1:nrow(otu_order_HAM)){
  otu_order_HAM[i,'ubiquity']<-length(which(otu_order_HAM[i,-ncol(otu_order_HAM)]!=0))/(ncol(otu_order_HAM)-1)
}

#### HGK
otu_order_HGK<-as.data.frame(matrix(nrow=length(HGK_order),ncol=18))
row.names(otu_order_HGK)<-HGK_order
for(i in c(1:8,10:18)){
  colnames(otu_order_HGK)[i]<-paste0('HGK_',i)
  aa<-eval(parse(text=paste0('HGK_',i)))
  for(j in HGK_order){
    if(length(which(aa$Order %in% j))==0){otu_order_HGK[j,i]<-0}
    else{otu_order_HGK[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_HGK$V9<-NULL
otu_order_HGK$ubiquity<-0
for(i in 1:nrow(otu_order_HGK)){
  otu_order_HGK[i,'ubiquity']<-length(which(otu_order_HGK[i,-ncol(otu_order_HGK)]!=0))/(ncol(otu_order_HGK)-1)
}

#### ILR
otu_order_ILR<-as.data.frame(matrix(nrow=length(ILR_order),ncol=24))
row.names(otu_order_ILR)<-ILR_order
for(i in 1:ncol(otu_order_ILR)){
  colnames(otu_order_ILR)[i]<-paste0('ILR_',i)
  aa<-eval(parse(text=paste0('ILR_',i)))
  for(j in ILR_order){
    if(length(which(aa$Order %in% j))==0){otu_order_ILR[j,i]<-0}
    else{otu_order_ILR[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_ILR$ubiquity<-0
for(i in 1:nrow(otu_order_ILR)){
  otu_order_ILR[i,'ubiquity']<-length(which(otu_order_ILR[i,-ncol(otu_order_ILR)]!=0))/(ncol(otu_order_ILR)-1)
}

#### LON
otu_order_LON<-as.data.frame(matrix(nrow=length(LON_order),ncol=24))
row.names(otu_order_LON)<-LON_order
for(i in c(1,3:6,8:24)){
  colnames(otu_order_LON)[i]<-paste0('LON_',i)
  aa<-eval(parse(text=paste0('LON_',i)))
  for(j in LON_order){
    if(length(which(aa$Order %in% j))==0){otu_order_LON[j,i]<-0}
    else{otu_order_LON[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_LON$V2<-NULL
otu_order_LON$V7<-NULL
otu_order_LON$ubiquity<-0
for(i in 1:nrow(otu_order_LON)){
  otu_order_LON[i,'ubiquity']<-length(which(otu_order_LON[i,-ncol(otu_order_LON)]!=0))/(ncol(otu_order_LON)-1)
}

#### MAR
otu_order_MAR<-as.data.frame(matrix(nrow=length(MAR_order),ncol=10))
row.names(otu_order_MAR)<-MAR_order
for(i in 1:ncol(otu_order_MAR)){
  colnames(otu_order_MAR)[i]<-paste0('MAR_',i)
  aa<-eval(parse(text=paste0('MAR_',i)))
  for(j in MAR_order){
    if(length(which(aa$Order %in% j))==0){otu_order_MAR[j,i]<-0}
    else{otu_order_MAR[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_MAR$ubiquity<-0
for(i in 1:nrow(otu_order_MAR)){
  otu_order_MAR[i,'ubiquity']<-length(which(otu_order_MAR[i,-ncol(otu_order_MAR)]!=0))/(ncol(otu_order_MAR)-1)
}

#### NYC
otu_order_NYC<-as.data.frame(matrix(nrow=length(NYC_order),ncol=26))
row.names(otu_order_NYC)<-NYC_order
for(i in 1:ncol(otu_order_NYC)){
  colnames(otu_order_NYC)[i]<-paste0('NYC_',i)
  aa<-eval(parse(text=paste0('NYC_',i)))
  for(j in NYC_order){
    if(length(which(aa$Order %in% j))==0){otu_order_NYC[j,i]<-0}
    else{otu_order_NYC[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_NYC$ubiquity<-0
for(i in 1:nrow(otu_order_NYC)){
  otu_order_NYC[i,'ubiquity']<-length(which(otu_order_NYC[i,-ncol(otu_order_NYC)]!=0))/(ncol(otu_order_NYC)-1)
}

#### OFA
otu_order_OFA<-as.data.frame(matrix(nrow=length(OFA_order),ncol=20))
row.names(otu_order_OFA)<-OFA_order
for(i in 1:ncol(otu_order_OFA)){
  colnames(otu_order_OFA)[i]<-paste0('OFA_',i)
  aa<-eval(parse(text=paste0('OFA_',i)))
  for(j in OFA_order){
    if(length(which(aa$Order %in% j))==0){otu_order_OFA[j,i]<-0}
    else{otu_order_OFA[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_OFA$ubiquity<-0
for(i in 1:nrow(otu_order_OFA)){
  otu_order_OFA[i,'ubiquity']<-length(which(otu_order_OFA[i,-ncol(otu_order_OFA)]!=0))/(ncol(otu_order_OFA)-1)
}

#### PXO
otu_order_PXO<-as.data.frame(matrix(nrow=length(PXO_order),ncol=20))
row.names(otu_order_PXO)<-PXO_order
for(i in 1:ncol(otu_order_PXO)){
  colnames(otu_order_PXO)[i]<-paste0('PXO_',i)
  aa<-eval(parse(text=paste0('PXO_',i)))
  for(j in PXO_order){
    if(length(which(aa$Order %in% j))==0){otu_order_PXO[j,i]<-0}
    else{otu_order_PXO[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_PXO$ubiquity<-0
for(i in 1:nrow(otu_order_PXO)){
  otu_order_PXO[i,'ubiquity']<-length(which(otu_order_PXO[i,-ncol(otu_order_PXO)]!=0))/(ncol(otu_order_PXO)-1)
}

#### SAC
otu_order_SAC<-as.data.frame(matrix(nrow=length(SAC_order),ncol=18))
row.names(otu_order_SAC)<-SAC_order
for(i in 1:ncol(otu_order_SAC)){
  colnames(otu_order_SAC)[i]<-paste0('SAC_',i)
  aa<-eval(parse(text=paste0('SAC_',i)))
  for(j in SAC_order){
    if(length(which(aa$Order %in% j))==0){otu_order_SAC[j,i]<-0}
    else{otu_order_SAC[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_SAC$ubiquity<-0
for(i in 1:nrow(otu_order_SAC)){
  otu_order_SAC[i,'ubiquity']<-length(which(otu_order_SAC[i,-ncol(otu_order_SAC)]!=0))/(ncol(otu_order_SAC)-1)
}

#### SAO
otu_order_SAO<-as.data.frame(matrix(nrow=length(SAO_order),ncol=24))
row.names(otu_order_SAO)<-SAO_order
for(i in 1:ncol(otu_order_SAO)){
  colnames(otu_order_SAO)[i]<-paste0('SAO_',i)
  aa<-eval(parse(text=paste0('SAO_',i)))
  for(j in SAO_order){
    if(length(which(aa$Order %in% j))==0){otu_order_SAO[j,i]<-0}
    else{otu_order_SAO[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_SAO$ubiquity<-0
for(i in 1:nrow(otu_order_SAO)){
  otu_order_SAO[i,'ubiquity']<-length(which(otu_order_SAO[i,-ncol(otu_order_SAO)]!=0))/(ncol(otu_order_SAO)-1)
}

#### SOF
otu_order_SOF<-as.data.frame(matrix(nrow=length(SOF_order),ncol=10))
row.names(otu_order_SOF)<-SOF_order
for(i in 1:ncol(otu_order_SOF)){
  colnames(otu_order_SOF)[i]<-paste0('SOF_',i)
  aa<-eval(parse(text=paste0('SOF_',i)))
  for(j in SOF_order){
    if(length(which(aa$Order %in% j))==0){otu_order_SOF[j,i]<-0}
    else{otu_order_SOF[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_SOF$ubiquity<-0
for(i in 1:nrow(otu_order_SOF)){
  otu_order_SOF[i,'ubiquity']<-length(which(otu_order_SOF[i,-ncol(otu_order_SOF)]!=0))/(ncol(otu_order_SOF)-1)
}

#### STO
otu_order_STO<-as.data.frame(matrix(nrow=length(STO_order),ncol=20))
row.names(otu_order_STO)<-STO_order
for(i in 1:ncol(otu_order_STO)){
  colnames(otu_order_STO)[i]<-paste0('STO_',i)
  aa<-eval(parse(text=paste0('STO_',i)))
  for(j in STO_order){
    if(length(which(aa$Order %in% j))==0){otu_order_STO[j,i]<-0}
    else{otu_order_STO[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_STO$ubiquity<-0
for(i in 1:nrow(otu_order_STO)){
  otu_order_STO[i,'ubiquity']<-length(which(otu_order_STO[i,-ncol(otu_order_STO)]!=0))/(ncol(otu_order_STO)-1)
}

#### TOK
otu_order_TOK<-as.data.frame(matrix(nrow=length(TOK_order),ncol=25))
row.names(otu_order_TOK)<-TOK_order
for(i in 1:ncol(otu_order_TOK)){
  colnames(otu_order_TOK)[i]<-paste0('TOK_',i)
  aa<-eval(parse(text=paste0('TOK_',i)))
  for(j in TOK_order){
    if(length(which(aa$Order %in% j))==0){otu_order_TOK[j,i]<-0}
    else{otu_order_TOK[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_TOK$ubiquity<-0
for(i in 1:nrow(otu_order_TOK)){
  otu_order_TOK[i,'ubiquity']<-length(which(otu_order_TOK[i,-ncol(otu_order_TOK)]!=0))/(ncol(otu_order_TOK)-1)
}

#### check ubiquity
check_ubiquity(otu_order_AKL,common_order)
check_ubiquity(otu_order_BER,common_order)
check_ubiquity(otu_order_BOG,common_order)
check_ubiquity(otu_order_HAM,common_order)
check_ubiquity(otu_order_HGK,common_order)
check_ubiquity(otu_order_ILR,common_order)
check_ubiquity(otu_order_LON,common_order)
check_ubiquity(otu_order_MAR,common_order)
check_ubiquity(otu_order_NYC,common_order)
check_ubiquity(otu_order_OFA,common_order)
check_ubiquity(otu_order_PXO,common_order)
check_ubiquity(otu_order_SAC,common_order)
check_ubiquity(otu_order_SAO,common_order)
check_ubiquity(otu_order_SOF,common_order)
check_ubiquity(otu_order_STO,common_order)
check_ubiquity(otu_order_TOK,common_order)

## ubiquity
otu_order_ubiquity<-as.data.frame(matrix(nrow=length(common_order),ncol=16))
row.names(otu_order_ubiquity)<-common_order
colnames(otu_order_ubiquity)<-c('AKL','BER','BOG','HAM',
                                  'HGK','ILR','LON','MAR',
                                  'NYC','OFA','PXO','SAC',
                                  'SAO','SOF','STO','TOK')
for(i in c('AKL','BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
  eval(parse(text=paste0('otu_order_ubiquity$',i,'<-','check_ubiquity(otu_order_',i,',common_order)$ubiquity')))
}

rm(aa,i,j)

## common order cbind - all samples
otu_common_order<-otu_order_AKL[which(row.names(otu_order_AKL) %in% common_order),-ncol(otu_order_AKL)]
for(i in c('BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
  aa<-eval(parse(text=paste0('otu_order_',i)))
  bb<-aa[which(row.names(aa) %in% common_order),-ncol(aa)]
  otu_common_order<-cbind(otu_common_order,bb)
}

## all order
otu_all_order<-as.data.frame(matrix(nrow=length(all_order),ncol=302))
row.names(otu_all_order)<-all_order
colnames(otu_all_order)<-c(paste0('AKL','_',c(1:14)),
                            paste0('BER','_',c(1:21)),
                            paste0('BOG','_',c(1:15)),
                            paste0('HAM','_',c(1:16)),
                            paste0('HGK','_',c(1:8,10:18)),
                            paste0('ILR','_',c(1:24)),
                            paste0('LON','_',c(1,3:6,8:24)),
                            paste0('MAR','_',c(1:10)),
                            paste0('NYC','_',c(1:26)),
                            paste0('OFA','_',c(1:20)),
                            paste0('PXO','_',c(1:20)),
                            paste0('SAC','_',c(1:18)),
                            paste0('SAO','_',c(1:24)),
                            paste0('SOF','_',c(1:10)),
                            paste0('STO','_',c(1:20)),
                            paste0('TOK','_',c(1:25)))

for(i in colnames(otu_all_order)){
  aa<-eval(parse(text=i))
  for(j in all_order){
    if(length(which(aa$Order %in% j))==0){otu_all_order[j,i]<-0}
    else{otu_all_order[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}

otu_all_order$ubiquity<-NA
for(i in 1:nrow(otu_all_order)){
  otu_all_order[i,'ubiquity']<-(302-length(which(otu_all_order[i,1:302]==0)))/302
}
otu_all_order[common_order,'ubiquity']


####### high ubiquity across cities ###########
####### high ubiquity across cities ###########
####### high ubiquity across cities ###########
####### high ubiquity across cities ###########
####### high ubiquity across cities ###########
############################################################### all species #########################################################################
all_species_count<-as.data.frame(matrix(nrow=length(all_species),ncol=16))
row.names(all_species_count)<-all_species
colnames(all_species_count)<-c('AKL','BER','BOG','HAM',
                                'HGK','ILR','LON','MAR',
                                'NYC','OFA','PXO','SAC',
                                'SAO','SOF','STO','TOK')
all_species_count[,]<-0
for(i in c('AKL','BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
  aa<-eval(parse(text=paste0(i,'_species')))
  all_species_count[which(row.names(all_species_count) %in% aa),i]<-1
}
all_species_count$ubiquity<-apply(all_species_count,1,sum)/16

############################################################### all family #########################################################################
all_family_count<-as.data.frame(matrix(nrow=length(all_family),ncol=16))
row.names(all_family_count)<-all_family
colnames(all_family_count)<-c('AKL','BER','BOG','HAM',
                               'HGK','ILR','LON','MAR',
                               'NYC','OFA','PXO','SAC',
                               'SAO','SOF','STO','TOK')
all_family_count[,]<-0
for(i in c('AKL','BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
  aa<-eval(parse(text=paste0(i,'_family')))
  all_family_count[which(row.names(all_family_count) %in% aa),i]<-1
}
all_family_count$ubiquity<-apply(all_family_count,1,sum)/16

############################################################### all order #########################################################################
all_order_count<-as.data.frame(matrix(nrow=length(all_order),ncol=16))
row.names(all_order_count)<-all_order
colnames(all_order_count)<-c('AKL','BER','BOG','HAM',
                              'HGK','ILR','LON','MAR',
                              'NYC','OFA','PXO','SAC',
                              'SAO','SOF','STO','TOK')
all_order_count[,]<-0
for(i in c('AKL','BER','BOG','HAM',
           'HGK','ILR','LON','MAR',
           'NYC','OFA','PXO','SAC',
           'SAO','SOF','STO','TOK')){
  aa<-eval(parse(text=paste0(i,'_order')))
  all_order_count[which(row.names(all_order_count) %in% aa),i]<-1
}
all_order_count$ubiquity<-apply(all_order_count,1,sum)/16

rm(aa,bb,i,j)

################################################################# Combine ###########################################################################
##
##
##
##
otu_common_combine<-rbind(otu_common_species,otu_common_family,otu_common_order)
save(list=c('otu_common_species',
            'otu_common_family',
            'otu_common_order',
            'otu_common_combine',
            'all_family_count',
            'all_order_count',
            'all_species_count',
            'otu_all_family',
            'otu_all_order',
            'otu_all_species',
            'city_feature'
),
file='Main dataset 2 OTU tables based on common features.RData')















