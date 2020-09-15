setwd('E:/University of Florida/My research/CAMDA 2019/Data/Data of mystery dataset/Mystery cities')
#### read city ####
#### read city ####
#### read city ####
read.city<-function(city_name){
  eval(parse(text=paste0(city_name,"<-read.table('",city_name,".txt',header=T)")))
  eval(parse(text=paste0('return(',city_name,')')))
}
Brisbane<-read.city('Brisbane')
Doha<-read.city('Doha')
Kiev<-read.city('Kiev')
Oslo<-read.city('Oslo')
Paris<-read.city('Paris')
Rio_de_Janerio<-read.city('Rio_de_Janerio')
Santiago_de_Chile<-read.city('Santiago_de_Chile')
Vienna<-read.city('Vienna')

#### input label ####
#### input label ####
#### input label ####
input.label<-function(city){
  eval(parse(text=paste0(city,"$label<-paste0('",city,"','_',1:nrow(",city,"))")))
  eval(parse(text=paste0('return(',city,')')))
}
Brisbane<-input.label('Brisbane')
Doha<-input.label('Doha')
Kiev<-input.label('Kiev')
Oslo<-input.label('Oslo')
Paris<-input.label('Paris')
Rio_de_Janerio<-input.label('Rio_de_Janerio')
Santiago_de_Chile<-input.label('Santiago_de_Chile')
Vienna<-input.label('Vienna')

#### read data ####
#### read data ####
#### read data ####
setwd('E:/University of Florida/My research/CAMDA 2019/Data/Data of mystery dataset/Mystery samples')
#### Brisbane ####
#### Brisbane ####
for(i in c(1:4,6,7)){
  aa<-read.delim(paste0(Brisbane[i,1],"/wgs_from_biom_final.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0(Brisbane[i,2],'<-aa')))
  bb<-read.delim(paste0(Brisbane[i,1],"/rep_set_aligned_tax_assignments.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0('otu_',i,'<-bb')))
}

for(i in c(1:4,6,7)){
  eval(parse(text=paste0('colnames(Brisbane_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0('Brisbane_',i,'<-merge(Brisbane_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

Brisbane_species<-NULL
for(i in c(1:4,6,7)){
  bb<-eval(parse(text=paste0('unique(Brisbane_',i,'$Species)')))
  Brisbane_species<-union(Brisbane_species,bb)
}
Brisbane_species<-na.omit(Brisbane_species)

Brisbane_family<-NULL
for(i in c(1:4,6,7)){
  bb<-eval(parse(text=paste0('unique(Brisbane_',i,'$Family)')))
  Brisbane_family<-union(Brisbane_family,bb)
}
Brisbane_family<-na.omit(Brisbane_family)

Brisbane_order<-NULL
for(i in c(1:4,6,7)){
  bb<-eval(parse(text=paste0('unique(Brisbane_',i,'$Order)')))
  Brisbane_order<-union(Brisbane_order,bb)
}
Brisbane_order<-na.omit(Brisbane_order)

rm(aa,i,j,bb)

#### Doha ####
#### Doha ####
for(i in 1:nrow(Doha)){
  aa<-read.delim(paste0(Doha[i,1],"/wgs_from_biom_final.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0(Doha[i,2],'<-aa')))
  bb<-read.delim(paste0(Doha[i,1],"/rep_set_aligned_tax_assignments.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0('otu_',i,'<-bb')))
}

for(i in 1:nrow(Doha)){
  eval(parse(text=paste0('colnames(Doha_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0('Doha_',i,'<-merge(Doha_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

Doha_species<-NULL
for(i in 1:nrow(Doha)){
  bb<-eval(parse(text=paste0('unique(Doha_',i,'$Species)')))
  Doha_species<-union(Doha_species,bb)
}
Doha_species<-na.omit(Doha_species)

Doha_family<-NULL
for(i in 1:nrow(Doha)){
  bb<-eval(parse(text=paste0('unique(Doha_',i,'$Family)')))
  Doha_family<-union(Doha_family,bb)
}
Doha_family<-na.omit(Doha_family)

Doha_order<-NULL
for(i in 1:nrow(Doha)){
  bb<-eval(parse(text=paste0('unique(Doha_',i,'$Order)')))
  Doha_order<-union(Doha_order,bb)
}
Doha_order<-na.omit(Doha_order)

rm(aa,i,j,bb)

#### Kiev ####
for(i in 2:nrow(Kiev)){
  aa<-read.delim(paste0(Kiev[i,1],"/wgs_from_biom_final.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0(Kiev[i,2],'<-aa')))
  bb<-read.delim(paste0(Kiev[i,1],"/rep_set_aligned_tax_assignments.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0('otu_',i,'<-bb')))
}

for(i in 2:nrow(Kiev)){
  eval(parse(text=paste0('colnames(Kiev_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0('Kiev_',i,'<-merge(Kiev_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

Kiev_species<-NULL
for(i in 2:nrow(Kiev)){
  bb<-eval(parse(text=paste0('unique(Kiev_',i,'$Species)')))
  Kiev_species<-union(Kiev_species,bb)
}
Kiev_species<-na.omit(Kiev_species)

Kiev_family<-NULL
for(i in 2:nrow(Kiev)){
  bb<-eval(parse(text=paste0('unique(Kiev_',i,'$Family)')))
  Kiev_family<-union(Kiev_family,bb)
}
Kiev_family<-na.omit(Kiev_family)

Kiev_order<-NULL
for(i in 2:nrow(Kiev)){
  bb<-eval(parse(text=paste0('unique(Kiev_',i,'$Order)')))
  Kiev_order<-union(Kiev_order,bb)
}
Kiev_order<-na.omit(Kiev_order)

rm(aa,i,j,bb)

#### Oslo ####
for(i in 1:nrow(Oslo)){
  aa<-read.delim(paste0(Oslo[i,1],"/wgs_from_biom_final.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0(Oslo[i,2],'<-aa')))
  bb<-read.delim(paste0(Oslo[i,1],"/rep_set_aligned_tax_assignments.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0('otu_',i,'<-bb')))
}

for(i in 1:nrow(Oslo)){
  eval(parse(text=paste0('colnames(Oslo_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0('Oslo_',i,'<-merge(Oslo_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

Oslo_species<-NULL
for(i in 1:nrow(Oslo)){
  bb<-eval(parse(text=paste0('unique(Oslo_',i,'$Species)')))
  Oslo_species<-union(Oslo_species,bb)
}
Oslo_species<-na.omit(Oslo_species)

Oslo_family<-NULL
for(i in 1:nrow(Oslo)){
  bb<-eval(parse(text=paste0('unique(Oslo_',i,'$Family)')))
  Oslo_family<-union(Oslo_family,bb)
}
Oslo_family<-na.omit(Oslo_family)

Oslo_order<-NULL
for(i in 1:nrow(Oslo)){
  bb<-eval(parse(text=paste0('unique(Oslo_',i,'$Order)')))
  Oslo_order<-union(Oslo_order,bb)
}
Oslo_order<-na.omit(Oslo_order)

rm(aa,i,j,bb)

#### Paris ####
for(i in c(1,2,5:8)){
  aa<-read.delim(paste0(Paris[i,1],"/wgs_from_biom_final.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0(Paris[i,2],'<-aa')))
  bb<-read.delim(paste0(Paris[i,1],"/rep_set_aligned_tax_assignments.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0('otu_',i,'<-bb')))
}

for(i in c(1,2,5:8)){
  eval(parse(text=paste0('colnames(Paris_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0('Paris_',i,'<-merge(Paris_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

Paris_species<-NULL
for(i in c(1,2,5:8)){
  bb<-eval(parse(text=paste0('unique(Paris_',i,'$Species)')))
  Paris_species<-union(Paris_species,bb)
}
Paris_species<-na.omit(Paris_species)

Paris_family<-NULL
for(i in c(1,2,5:8)){
  bb<-eval(parse(text=paste0('unique(Paris_',i,'$Family)')))
  Paris_family<-union(Paris_family,bb)
}
Paris_family<-na.omit(Paris_family)

Paris_order<-NULL
for(i in c(1,2,5:8)){
  bb<-eval(parse(text=paste0('unique(Paris_',i,'$Order)')))
  Paris_order<-union(Paris_order,bb)
}
Paris_order<-na.omit(Paris_order)

rm(aa,i,j,bb)

#### Rio_de_Janerio ####
for(i in 1:nrow(Rio_de_Janerio)){
  aa<-read.delim(paste0(Rio_de_Janerio[i,1],"/wgs_from_biom_final.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0(Rio_de_Janerio[i,2],'<-aa')))
  bb<-read.delim(paste0(Rio_de_Janerio[i,1],"/rep_set_aligned_tax_assignments.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0('otu_',i,'<-bb')))
}

for(i in 1:nrow(Rio_de_Janerio)){
  eval(parse(text=paste0('colnames(Rio_de_Janerio_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0('Rio_de_Janerio_',i,'<-merge(Rio_de_Janerio_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

Rio_de_Janerio_species<-NULL
for(i in 1:nrow(Rio_de_Janerio)){
  bb<-eval(parse(text=paste0('unique(Rio_de_Janerio_',i,'$Species)')))
  Rio_de_Janerio_species<-union(Rio_de_Janerio_species,bb)
}
Rio_de_Janerio_species<-na.omit(Rio_de_Janerio_species)

Rio_de_Janerio_family<-NULL
for(i in 1:nrow(Rio_de_Janerio)){
  bb<-eval(parse(text=paste0('unique(Rio_de_Janerio_',i,'$Family)')))
  Rio_de_Janerio_family<-union(Rio_de_Janerio_family,bb)
}
Rio_de_Janerio_family<-na.omit(Rio_de_Janerio_family)

Rio_de_Janerio_order<-NULL
for(i in 1:nrow(Rio_de_Janerio)){
  bb<-eval(parse(text=paste0('unique(Rio_de_Janerio_',i,'$Order)')))
  Rio_de_Janerio_order<-union(Rio_de_Janerio_order,bb)
}
Rio_de_Janerio_order<-na.omit(Rio_de_Janerio_order)

rm(aa,i,j,bb)

#### Santiago_de_Chile ####
for(i in 1:nrow(Santiago_de_Chile)){
  aa<-read.delim(paste0(Santiago_de_Chile[i,1],"/wgs_from_biom_final.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0(Santiago_de_Chile[i,2],'<-aa')))
  bb<-read.delim(paste0(Santiago_de_Chile[i,1],"/rep_set_aligned_tax_assignments.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0('otu_',i,'<-bb')))
}

for(i in 1:nrow(Santiago_de_Chile)){
  eval(parse(text=paste0('colnames(Santiago_de_Chile_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0('Santiago_de_Chile_',i,'<-merge(Santiago_de_Chile_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

Santiago_de_Chile_species<-NULL
for(i in 1:nrow(Santiago_de_Chile)){
  bb<-eval(parse(text=paste0('unique(Santiago_de_Chile_',i,'$Species)')))
  Santiago_de_Chile_species<-union(Santiago_de_Chile_species,bb)
}
Santiago_de_Chile_species<-na.omit(Santiago_de_Chile_species)

Santiago_de_Chile_family<-NULL
for(i in 1:nrow(Santiago_de_Chile)){
  bb<-eval(parse(text=paste0('unique(Santiago_de_Chile_',i,'$Family)')))
  Santiago_de_Chile_family<-union(Santiago_de_Chile_family,bb)
}
Santiago_de_Chile_family<-na.omit(Santiago_de_Chile_family)

Santiago_de_Chile_order<-NULL
for(i in 1:nrow(Santiago_de_Chile)){
  bb<-eval(parse(text=paste0('unique(Santiago_de_Chile_',i,'$Order)')))
  Santiago_de_Chile_order<-union(Santiago_de_Chile_order,bb)
}
Santiago_de_Chile_order<-na.omit(Santiago_de_Chile_order)

rm(aa,i,j,bb)

#### Vienna ####
for(i in 1:nrow(Vienna)){
  aa<-read.delim(paste0(Vienna[i,1],"/wgs_from_biom_final.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0(Vienna[i,2],'<-aa')))
  bb<-read.delim(paste0(Vienna[i,1],"/rep_set_aligned_tax_assignments.txt"), header=FALSE, comment.char="#",stringsAsFactors=FALSE)
  eval(parse(text=paste0('otu_',i,'<-bb')))
}

for(i in 1:nrow(Vienna)){
  eval(parse(text=paste0('colnames(Vienna_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0('Vienna_',i,'<-merge(Vienna_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

Vienna_species<-NULL
for(i in 1:nrow(Vienna)){
  bb<-eval(parse(text=paste0('unique(Vienna_',i,'$Species)')))
  Vienna_species<-union(Vienna_species,bb)
}
Vienna_species<-na.omit(Vienna_species)

Vienna_family<-NULL
for(i in 1:nrow(Vienna)){
  bb<-eval(parse(text=paste0('unique(Vienna_',i,'$Family)')))
  Vienna_family<-union(Vienna_family,bb)
}
Vienna_family<-na.omit(Vienna_family)

Vienna_order<-NULL
for(i in 1:nrow(Vienna)){
  bb<-eval(parse(text=paste0('unique(Vienna_',i,'$Order)')))
  Vienna_order<-union(Vienna_order,bb)
}
Vienna_order<-na.omit(Vienna_order)

rm(aa,i,j,bb)

#############################################################################################################################
################################################ save the data ##############################################################
setwd("E:/University of Florida/My research/CAMDA 2019/Data/Data of mystery dataset")
save.image(file='Mystery dataset 1 OTU tables for all samples.RData')


