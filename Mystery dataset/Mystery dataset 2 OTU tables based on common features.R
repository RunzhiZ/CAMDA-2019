load('Mystery dataset 1 OTU tables for all samples.RData')
l_species_mystery<-list(as.vector(Brisbane_species),
                        as.vector(Doha_species),
                        as.vector(Kiev_species),
                        as.vector(Oslo_species),
                        as.vector(Paris_species),
                        as.vector(Rio_de_Janerio_species),
                        as.vector(Santiago_de_Chile_species),
                        as.vector(Vienna_species))
common_species_mystery<-Reduce(intersect,l_species_mystery)

l_family_mystery<-list(as.vector(Brisbane_family),
                       as.vector(Doha_family),
                       as.vector(Kiev_family),
                       as.vector(Oslo_family),
                       as.vector(Paris_family),
                       as.vector(Rio_de_Janerio_family),
                       as.vector(Santiago_de_Chile_family),
                       as.vector(Vienna_family))
common_family_mystery<-Reduce(intersect,l_family_mystery)

l_order_mystery<-list(as.vector(Brisbane_order),
                      as.vector(Doha_order),
                      as.vector(Kiev_order),
                      as.vector(Oslo_order),
                      as.vector(Paris_order),
                      as.vector(Rio_de_Janerio_order),
                      as.vector(Santiago_de_Chile_order),
                      as.vector(Vienna_order))
common_order_mystery<-Reduce(intersect,l_order_mystery)

######## species ########
######## species ########
######## species ########
######## species ########
#### Brisbane
otu_species_Brisbane<-as.data.frame(matrix(nrow=length(Brisbane_species),ncol=7))
row.names(otu_species_Brisbane)<-Brisbane_species
for(i in c(1:4,6:7)){
  colnames(otu_species_Brisbane)[i]<-paste0('Brisbane_',i)
  aa<-eval(parse(text=paste0('Brisbane_',i)))
  for(j in Brisbane_species){
    if(length(which(aa$Species %in% j))==0){otu_species_Brisbane[j,i]<-0}
    else{otu_species_Brisbane[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_Brisbane$V5<-NULL
otu_species_Brisbane$ubiquity<-0
for(i in 1:nrow(otu_species_Brisbane)){
  otu_species_Brisbane[i,'ubiquity']<-length(which(otu_species_Brisbane[i,-ncol(otu_species_Brisbane)]!=0))/(ncol(otu_species_Brisbane)-1)
}
#### Doha
otu_species_Doha<-as.data.frame(matrix(nrow=length(Doha_species),ncol=3))
row.names(otu_species_Doha)<-Doha_species
for(i in 1:3){
  colnames(otu_species_Doha)[i]<-paste0('Doha_',i)
  aa<-eval(parse(text=paste0('Doha_',i)))
  for(j in Doha_species){
    if(length(which(aa$Species %in% j))==0){otu_species_Doha[j,i]<-0}
    else{otu_species_Doha[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_Doha$ubiquity<-0
for(i in 1:nrow(otu_species_Doha)){
  otu_species_Doha[i,'ubiquity']<-length(which(otu_species_Doha[i,-ncol(otu_species_Doha)]!=0))/(ncol(otu_species_Doha)-1)
}
#### Kiev
otu_species_Kiev<-as.data.frame(matrix(nrow=length(Kiev_species),ncol=8))
row.names(otu_species_Kiev)<-Kiev_species
for(i in 2:8){
  colnames(otu_species_Kiev)[i]<-paste0('Kiev_',i)
  aa<-eval(parse(text=paste0('Kiev_',i)))
  for(j in Kiev_species){
    if(length(which(aa$Species %in% j))==0){otu_species_Kiev[j,i]<-0}
    else{otu_species_Kiev[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_Kiev$V1<-NULL
otu_species_Kiev$ubiquity<-0
for(i in 1:nrow(otu_species_Kiev)){
  otu_species_Kiev[i,'ubiquity']<-length(which(otu_species_Kiev[i,-ncol(otu_species_Kiev)]!=0))/(ncol(otu_species_Kiev)-1)
}
#### Oslo
otu_species_Oslo<-as.data.frame(matrix(nrow=length(Oslo_species),ncol=12))
row.names(otu_species_Oslo)<-Oslo_species
for(i in 1:12){
  colnames(otu_species_Oslo)[i]<-paste0('Oslo_',i)
  aa<-eval(parse(text=paste0('Oslo_',i)))
  for(j in Oslo_species){
    if(length(which(aa$Species %in% j))==0){otu_species_Oslo[j,i]<-0}
    else{otu_species_Oslo[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_Oslo$ubiquity<-0
for(i in 1:nrow(otu_species_Oslo)){
  otu_species_Oslo[i,'ubiquity']<-length(which(otu_species_Oslo[i,-ncol(otu_species_Oslo)]!=0))/(ncol(otu_species_Oslo)-1)
}
#### Paris
otu_species_Paris<-as.data.frame(matrix(nrow=length(Paris_species),ncol=8))
row.names(otu_species_Paris)<-Paris_species
for(i in c(1:2,5:8)){
  colnames(otu_species_Paris)[i]<-paste0('Paris_',i)
  aa<-eval(parse(text=paste0('Paris_',i)))
  for(j in Paris_species){
    if(length(which(aa$Species %in% j))==0){otu_species_Paris[j,i]<-0}
    else{otu_species_Paris[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_Paris$V3<-NULL
otu_species_Paris$V4<-NULL
otu_species_Paris$ubiquity<-0
for(i in 1:nrow(otu_species_Paris)){
  otu_species_Paris[i,'ubiquity']<-length(which(otu_species_Paris[i,-ncol(otu_species_Paris)]!=0))/(ncol(otu_species_Paris)-1)
}
#### Rio_de_Janerio
otu_species_Rio_de_Janerio<-as.data.frame(matrix(nrow=length(Rio_de_Janerio_species),ncol=12))
row.names(otu_species_Rio_de_Janerio)<-Rio_de_Janerio_species
for(i in 1:12){
  colnames(otu_species_Rio_de_Janerio)[i]<-paste0('Rio_de_Janerio_',i)
  aa<-eval(parse(text=paste0('Rio_de_Janerio_',i)))
  for(j in Rio_de_Janerio_species){
    if(length(which(aa$Species %in% j))==0){otu_species_Rio_de_Janerio[j,i]<-0}
    else{otu_species_Rio_de_Janerio[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_Rio_de_Janerio$ubiquity<-0
for(i in 1:nrow(otu_species_Rio_de_Janerio)){
  otu_species_Rio_de_Janerio[i,'ubiquity']<-length(which(otu_species_Rio_de_Janerio[i,-ncol(otu_species_Rio_de_Janerio)]!=0))/(ncol(otu_species_Rio_de_Janerio)-1)
}
#### Santiago_de_Chile
otu_species_Santiago_de_Chile<-as.data.frame(matrix(nrow=length(Santiago_de_Chile_species),ncol=6))
row.names(otu_species_Santiago_de_Chile)<-Santiago_de_Chile_species
for(i in 1:6){
  colnames(otu_species_Santiago_de_Chile)[i]<-paste0('Santiago_de_Chile_',i)
  aa<-eval(parse(text=paste0('Santiago_de_Chile_',i)))
  for(j in Santiago_de_Chile_species){
    if(length(which(aa$Species %in% j))==0){otu_species_Santiago_de_Chile[j,i]<-0}
    else{otu_species_Santiago_de_Chile[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_Santiago_de_Chile$ubiquity<-0
for(i in 1:nrow(otu_species_Santiago_de_Chile)){
  otu_species_Santiago_de_Chile[i,'ubiquity']<-length(which(otu_species_Santiago_de_Chile[i,-ncol(otu_species_Santiago_de_Chile)]!=0))/(ncol(otu_species_Santiago_de_Chile)-1)
}
#### Vienna
otu_species_Vienna<-as.data.frame(matrix(nrow=length(Vienna_species),ncol=5))
row.names(otu_species_Vienna)<-Vienna_species
for(i in 1:5){
  colnames(otu_species_Vienna)[i]<-paste0('Vienna_',i)
  aa<-eval(parse(text=paste0('Vienna_',i)))
  for(j in Vienna_species){
    if(length(which(aa$Species %in% j))==0){otu_species_Vienna[j,i]<-0}
    else{otu_species_Vienna[j,i]<-sum(aa[which(aa$Species==j),'count'])}
  }
}
otu_species_Vienna$ubiquity<-0
for(i in 1:nrow(otu_species_Vienna)){
  otu_species_Vienna[i,'ubiquity']<-length(which(otu_species_Vienna[i,-ncol(otu_species_Vienna)]!=0))/(ncol(otu_species_Vienna)-1)
}

######## family ########
######## family ########
######## family ########
######## family ########
#### Brisbane
otu_family_Brisbane<-as.data.frame(matrix(nrow=length(Brisbane_family),ncol=7))
row.names(otu_family_Brisbane)<-Brisbane_family
for(i in c(1:4,6:7)){
  colnames(otu_family_Brisbane)[i]<-paste0('Brisbane_',i)
  aa<-eval(parse(text=paste0('Brisbane_',i)))
  for(j in Brisbane_family){
    if(length(which(aa$Family %in% j))==0){otu_family_Brisbane[j,i]<-0}
    else{otu_family_Brisbane[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_Brisbane$V5<-NULL
otu_family_Brisbane$ubiquity<-0
for(i in 1:nrow(otu_family_Brisbane)){
  otu_family_Brisbane[i,'ubiquity']<-length(which(otu_family_Brisbane[i,-ncol(otu_family_Brisbane)]!=0))/(ncol(otu_family_Brisbane)-1)
}
#### Doha
otu_family_Doha<-as.data.frame(matrix(nrow=length(Doha_family),ncol=3))
row.names(otu_family_Doha)<-Doha_family
for(i in 1:3){
  colnames(otu_family_Doha)[i]<-paste0('Doha_',i)
  aa<-eval(parse(text=paste0('Doha_',i)))
  for(j in Doha_family){
    if(length(which(aa$Family %in% j))==0){otu_family_Doha[j,i]<-0}
    else{otu_family_Doha[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_Doha$ubiquity<-0
for(i in 1:nrow(otu_family_Doha)){
  otu_family_Doha[i,'ubiquity']<-length(which(otu_family_Doha[i,-ncol(otu_family_Doha)]!=0))/(ncol(otu_family_Doha)-1)
}
#### Kiev
otu_family_Kiev<-as.data.frame(matrix(nrow=length(Kiev_family),ncol=8))
row.names(otu_family_Kiev)<-Kiev_family
for(i in 2:8){
  colnames(otu_family_Kiev)[i]<-paste0('Kiev_',i)
  aa<-eval(parse(text=paste0('Kiev_',i)))
  for(j in Kiev_family){
    if(length(which(aa$Family %in% j))==0){otu_family_Kiev[j,i]<-0}
    else{otu_family_Kiev[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_Kiev$V1<-NULL
otu_family_Kiev$ubiquity<-0
for(i in 1:nrow(otu_family_Kiev)){
  otu_family_Kiev[i,'ubiquity']<-length(which(otu_family_Kiev[i,-ncol(otu_family_Kiev)]!=0))/(ncol(otu_family_Kiev)-1)
}
#### Oslo
otu_family_Oslo<-as.data.frame(matrix(nrow=length(Oslo_family),ncol=12))
row.names(otu_family_Oslo)<-Oslo_family
for(i in 1:12){
  colnames(otu_family_Oslo)[i]<-paste0('Oslo_',i)
  aa<-eval(parse(text=paste0('Oslo_',i)))
  for(j in Oslo_family){
    if(length(which(aa$Family %in% j))==0){otu_family_Oslo[j,i]<-0}
    else{otu_family_Oslo[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_Oslo$ubiquity<-0
for(i in 1:nrow(otu_family_Oslo)){
  otu_family_Oslo[i,'ubiquity']<-length(which(otu_family_Oslo[i,-ncol(otu_family_Oslo)]!=0))/(ncol(otu_family_Oslo)-1)
}
#### Paris
otu_family_Paris<-as.data.frame(matrix(nrow=length(Paris_family),ncol=8))
row.names(otu_family_Paris)<-Paris_family
for(i in c(1:2,5:8)){
  colnames(otu_family_Paris)[i]<-paste0('Paris_',i)
  aa<-eval(parse(text=paste0('Paris_',i)))
  for(j in Paris_family){
    if(length(which(aa$Family %in% j))==0){otu_family_Paris[j,i]<-0}
    else{otu_family_Paris[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_Paris$V3<-NULL
otu_family_Paris$V4<-NULL
otu_family_Paris$ubiquity<-0
for(i in 1:nrow(otu_family_Paris)){
  otu_family_Paris[i,'ubiquity']<-length(which(otu_family_Paris[i,-ncol(otu_family_Paris)]!=0))/(ncol(otu_family_Paris)-1)
}
#### Rio_de_Janerio
otu_family_Rio_de_Janerio<-as.data.frame(matrix(nrow=length(Rio_de_Janerio_family),ncol=12))
row.names(otu_family_Rio_de_Janerio)<-Rio_de_Janerio_family
for(i in 1:12){
  colnames(otu_family_Rio_de_Janerio)[i]<-paste0('Rio_de_Janerio_',i)
  aa<-eval(parse(text=paste0('Rio_de_Janerio_',i)))
  for(j in Rio_de_Janerio_family){
    if(length(which(aa$Family %in% j))==0){otu_family_Rio_de_Janerio[j,i]<-0}
    else{otu_family_Rio_de_Janerio[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_Rio_de_Janerio$ubiquity<-0
for(i in 1:nrow(otu_family_Rio_de_Janerio)){
  otu_family_Rio_de_Janerio[i,'ubiquity']<-length(which(otu_family_Rio_de_Janerio[i,-ncol(otu_family_Rio_de_Janerio)]!=0))/(ncol(otu_family_Rio_de_Janerio)-1)
}
#### Santiago_de_Chile
otu_family_Santiago_de_Chile<-as.data.frame(matrix(nrow=length(Santiago_de_Chile_family),ncol=6))
row.names(otu_family_Santiago_de_Chile)<-Santiago_de_Chile_family
for(i in 1:6){
  colnames(otu_family_Santiago_de_Chile)[i]<-paste0('Santiago_de_Chile_',i)
  aa<-eval(parse(text=paste0('Santiago_de_Chile_',i)))
  for(j in Santiago_de_Chile_family){
    if(length(which(aa$Family %in% j))==0){otu_family_Santiago_de_Chile[j,i]<-0}
    else{otu_family_Santiago_de_Chile[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_Santiago_de_Chile$ubiquity<-0
for(i in 1:nrow(otu_family_Santiago_de_Chile)){
  otu_family_Santiago_de_Chile[i,'ubiquity']<-length(which(otu_family_Santiago_de_Chile[i,-ncol(otu_family_Santiago_de_Chile)]!=0))/(ncol(otu_family_Santiago_de_Chile)-1)
}
#### Vienna
otu_family_Vienna<-as.data.frame(matrix(nrow=length(Vienna_family),ncol=5))
row.names(otu_family_Vienna)<-Vienna_family
for(i in 1:5){
  colnames(otu_family_Vienna)[i]<-paste0('Vienna_',i)
  aa<-eval(parse(text=paste0('Vienna_',i)))
  for(j in Vienna_family){
    if(length(which(aa$Family %in% j))==0){otu_family_Vienna[j,i]<-0}
    else{otu_family_Vienna[j,i]<-sum(aa[which(aa$Family==j),'count'])}
  }
}
otu_family_Vienna$ubiquity<-0
for(i in 1:nrow(otu_family_Vienna)){
  otu_family_Vienna[i,'ubiquity']<-length(which(otu_family_Vienna[i,-ncol(otu_family_Vienna)]!=0))/(ncol(otu_family_Vienna)-1)
}

######## order ########
######## order ########
######## order ########
######## order ########
#### Brisbane
otu_order_Brisbane<-as.data.frame(matrix(nrow=length(Brisbane_order),ncol=7))
row.names(otu_order_Brisbane)<-Brisbane_order
for(i in c(1:4,6:7)){
  colnames(otu_order_Brisbane)[i]<-paste0('Brisbane_',i)
  aa<-eval(parse(text=paste0('Brisbane_',i)))
  for(j in Brisbane_order){
    if(length(which(aa$Order %in% j))==0){otu_order_Brisbane[j,i]<-0}
    else{otu_order_Brisbane[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_Brisbane$V5<-NULL
otu_order_Brisbane$ubiquity<-0
for(i in 1:nrow(otu_order_Brisbane)){
  otu_order_Brisbane[i,'ubiquity']<-length(which(otu_order_Brisbane[i,-ncol(otu_order_Brisbane)]!=0))/(ncol(otu_order_Brisbane)-1)
}
#### Doha
otu_order_Doha<-as.data.frame(matrix(nrow=length(Doha_order),ncol=3))
row.names(otu_order_Doha)<-Doha_order
for(i in 1:3){
  colnames(otu_order_Doha)[i]<-paste0('Doha_',i)
  aa<-eval(parse(text=paste0('Doha_',i)))
  for(j in Doha_order){
    if(length(which(aa$Order %in% j))==0){otu_order_Doha[j,i]<-0}
    else{otu_order_Doha[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_Doha$ubiquity<-0
for(i in 1:nrow(otu_order_Doha)){
  otu_order_Doha[i,'ubiquity']<-length(which(otu_order_Doha[i,-ncol(otu_order_Doha)]!=0))/(ncol(otu_order_Doha)-1)
}
#### Kiev
otu_order_Kiev<-as.data.frame(matrix(nrow=length(Kiev_order),ncol=8))
row.names(otu_order_Kiev)<-Kiev_order
for(i in 2:8){
  colnames(otu_order_Kiev)[i]<-paste0('Kiev_',i)
  aa<-eval(parse(text=paste0('Kiev_',i)))
  for(j in Kiev_order){
    if(length(which(aa$Order %in% j))==0){otu_order_Kiev[j,i]<-0}
    else{otu_order_Kiev[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_Kiev$V1<-NULL
otu_order_Kiev$ubiquity<-0
for(i in 1:nrow(otu_order_Kiev)){
  otu_order_Kiev[i,'ubiquity']<-length(which(otu_order_Kiev[i,-ncol(otu_order_Kiev)]!=0))/(ncol(otu_order_Kiev)-1)
}
#### Oslo
otu_order_Oslo<-as.data.frame(matrix(nrow=length(Oslo_order),ncol=12))
row.names(otu_order_Oslo)<-Oslo_order
for(i in 1:12){
  colnames(otu_order_Oslo)[i]<-paste0('Oslo_',i)
  aa<-eval(parse(text=paste0('Oslo_',i)))
  for(j in Oslo_order){
    if(length(which(aa$Order %in% j))==0){otu_order_Oslo[j,i]<-0}
    else{otu_order_Oslo[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_Oslo$ubiquity<-0
for(i in 1:nrow(otu_order_Oslo)){
  otu_order_Oslo[i,'ubiquity']<-length(which(otu_order_Oslo[i,-ncol(otu_order_Oslo)]!=0))/(ncol(otu_order_Oslo)-1)
}
#### Paris
otu_order_Paris<-as.data.frame(matrix(nrow=length(Paris_order),ncol=8))
row.names(otu_order_Paris)<-Paris_order
for(i in c(1:2,5:8)){
  colnames(otu_order_Paris)[i]<-paste0('Paris_',i)
  aa<-eval(parse(text=paste0('Paris_',i)))
  for(j in Paris_order){
    if(length(which(aa$Order %in% j))==0){otu_order_Paris[j,i]<-0}
    else{otu_order_Paris[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_Paris$V3<-NULL
otu_order_Paris$V4<-NULL
otu_order_Paris$ubiquity<-0
for(i in 1:nrow(otu_order_Paris)){
  otu_order_Paris[i,'ubiquity']<-length(which(otu_order_Paris[i,-ncol(otu_order_Paris)]!=0))/(ncol(otu_order_Paris)-1)
}
#### Rio_de_Janerio
otu_order_Rio_de_Janerio<-as.data.frame(matrix(nrow=length(Rio_de_Janerio_order),ncol=12))
row.names(otu_order_Rio_de_Janerio)<-Rio_de_Janerio_order
for(i in 1:12){
  colnames(otu_order_Rio_de_Janerio)[i]<-paste0('Rio_de_Janerio_',i)
  aa<-eval(parse(text=paste0('Rio_de_Janerio_',i)))
  for(j in Rio_de_Janerio_order){
    if(length(which(aa$Order %in% j))==0){otu_order_Rio_de_Janerio[j,i]<-0}
    else{otu_order_Rio_de_Janerio[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_Rio_de_Janerio$ubiquity<-0
for(i in 1:nrow(otu_order_Rio_de_Janerio)){
  otu_order_Rio_de_Janerio[i,'ubiquity']<-length(which(otu_order_Rio_de_Janerio[i,-ncol(otu_order_Rio_de_Janerio)]!=0))/(ncol(otu_order_Rio_de_Janerio)-1)
}
#### Santiago_de_Chile
otu_order_Santiago_de_Chile<-as.data.frame(matrix(nrow=length(Santiago_de_Chile_order),ncol=6))
row.names(otu_order_Santiago_de_Chile)<-Santiago_de_Chile_order
for(i in 1:6){
  colnames(otu_order_Santiago_de_Chile)[i]<-paste0('Santiago_de_Chile_',i)
  aa<-eval(parse(text=paste0('Santiago_de_Chile_',i)))
  for(j in Santiago_de_Chile_order){
    if(length(which(aa$Order %in% j))==0){otu_order_Santiago_de_Chile[j,i]<-0}
    else{otu_order_Santiago_de_Chile[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_Santiago_de_Chile$ubiquity<-0
for(i in 1:nrow(otu_order_Santiago_de_Chile)){
  otu_order_Santiago_de_Chile[i,'ubiquity']<-length(which(otu_order_Santiago_de_Chile[i,-ncol(otu_order_Santiago_de_Chile)]!=0))/(ncol(otu_order_Santiago_de_Chile)-1)
}
#### Vienna
otu_order_Vienna<-as.data.frame(matrix(nrow=length(Vienna_order),ncol=5))
row.names(otu_order_Vienna)<-Vienna_order
for(i in 1:5){
  colnames(otu_order_Vienna)[i]<-paste0('Vienna_',i)
  aa<-eval(parse(text=paste0('Vienna_',i)))
  for(j in Vienna_order){
    if(length(which(aa$Order %in% j))==0){otu_order_Vienna[j,i]<-0}
    else{otu_order_Vienna[j,i]<-sum(aa[which(aa$Order==j),'count'])}
  }
}
otu_order_Vienna$ubiquity<-0
for(i in 1:nrow(otu_order_Vienna)){
  otu_order_Vienna[i,'ubiquity']<-length(which(otu_order_Vienna[i,-ncol(otu_order_Vienna)]!=0))/(ncol(otu_order_Vienna)-1)
}

rm(aa,bb,i,j)
################################################################# Combine ###########################################################################
##
##
##
##
obtain_common<-function(data,common){
  output<-data[row.names(data)%in%common,-ncol(data)]
  return(output)
}
otu_common_species_mystery<-cbind(obtain_common(otu_species_Brisbane,common_species_mystery),
                                  obtain_common(otu_species_Doha,common_species_mystery),
                                  obtain_common(otu_species_Kiev,common_species_mystery),
                                  obtain_common(otu_species_Oslo,common_species_mystery),
                                  obtain_common(otu_species_Paris,common_species_mystery),
                                  obtain_common(otu_species_Rio_de_Janerio,common_species_mystery),
                                  obtain_common(otu_species_Santiago_de_Chile,common_species_mystery),
                                  obtain_common(otu_species_Vienna,common_species_mystery))

otu_common_family_mystery<-cbind(obtain_common(otu_family_Brisbane,common_family_mystery),
                                 obtain_common(otu_family_Doha,common_family_mystery),
                                 obtain_common(otu_family_Kiev,common_family_mystery),
                                 obtain_common(otu_family_Oslo,common_family_mystery),
                                 obtain_common(otu_family_Paris,common_family_mystery),
                                 obtain_common(otu_family_Rio_de_Janerio,common_family_mystery),
                                 obtain_common(otu_family_Santiago_de_Chile,common_family_mystery),
                                 obtain_common(otu_family_Vienna,common_family_mystery))

otu_common_order_mystery<-cbind(obtain_common(otu_order_Brisbane,common_order_mystery),
                                obtain_common(otu_order_Doha,common_order_mystery),
                                obtain_common(otu_order_Kiev,common_order_mystery),
                                obtain_common(otu_order_Oslo,common_order_mystery),
                                obtain_common(otu_order_Paris,common_order_mystery),
                                obtain_common(otu_order_Rio_de_Janerio,common_order_mystery),
                                obtain_common(otu_order_Santiago_de_Chile,common_order_mystery),
                                obtain_common(otu_order_Vienna,common_order_mystery))



save(list=c('otu_common_species_mystery',
            'otu_common_family_mystery',
            'otu_common_order_mystery'
),
file='Mystery dataset 2 OTU tables based on common features.RData')

save(list=c('otu_species_Brisbane',
            'otu_species_Doha',
            'otu_species_Kiev',
            'otu_species_Oslo',
            'otu_species_Paris',
            'otu_species_Rio_de_Janerio',
            'otu_species_Santiago_de_Chile',
            'otu_species_Vienna',
            'otu_family_Brisbane',
            'otu_family_Doha',
            'otu_family_Kiev',
            'otu_family_Oslo',
            'otu_family_Paris',
            'otu_family_Rio_de_Janerio',
            'otu_family_Santiago_de_Chile',
            'otu_family_Vienna',
            'otu_order_Brisbane',
            'otu_order_Doha',
            'otu_order_Kiev',
            'otu_order_Oslo',
            'otu_order_Paris',
            'otu_order_Rio_de_Janerio',
            'otu_order_Santiago_de_Chile',
            'otu_order_Vienna'),
     file='Mystery dataset for prediction.RData')



