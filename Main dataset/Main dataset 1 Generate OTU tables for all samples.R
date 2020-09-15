#############################################################################################################################
################################################ AKL: 14 samples#############################################################

AKL_1 <- read.delim("AKL/AKL-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_2 <- read.delim("AKL/AKL-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_3 <- read.delim("AKL/AKL-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_4 <- read.delim("AKL/AKL-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_5 <- read.delim("AKL/AKL-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_6 <- read.delim("AKL/AKL-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_7 <- read.delim("AKL/AKL-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_8 <- read.delim("AKL/AKL-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_9 <- read.delim("AKL/AKL-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_10 <- read.delim("AKL/AKL-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_11 <- read.delim("AKL/AKL-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_12 <- read.delim("AKL/AKL-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_13 <- read.delim("AKL/AKL-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
AKL_14 <- read.delim("AKL/AKL-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)

otu_1 <- read.delim("AKL/AKL-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("AKL/AKL-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("AKL/AKL-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("AKL/AKL-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("AKL/AKL-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("AKL/AKL-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("AKL/AKL-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("AKL/AKL-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("AKL/AKL-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("AKL/AKL-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("AKL/AKL-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("AKL/AKL-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("AKL/AKL-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("AKL/AKL-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)

for(i in 1:14){
  eval(parse(text=paste0('colnames(AKL_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('AKL_',i,'<-merge(AKL_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

AKL_species<-NULL
for(i in 1:14){
  bb<-eval(parse(text=paste0('unique(AKL_',i,'$Species)')))
  AKL_species<-union(AKL_species,bb)
}
AKL_species<-na.omit(AKL_species)

AKL_family<-NULL
for(i in 1:14){
  bb<-eval(parse(text=paste0('unique(AKL_',i,'$Family)')))
  AKL_family<-union(AKL_family,bb)
}
AKL_family<-na.omit(AKL_family)

AKL_order<-NULL
for(i in 1:14){
  bb<-eval(parse(text=paste0('unique(AKL_',i,'$Order)')))
  AKL_order<-union(AKL_order,bb)
}
AKL_order<-na.omit(AKL_order)

rm(aa,i,j,bb)

#############################################################################################################################
################################################ BER: 21 samples#############################################################

BER_1 <- read.delim("BER/BER-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_2 <- read.delim("BER/BER-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_3 <- read.delim("BER/BER-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_4 <- read.delim("BER/BER-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_5 <- read.delim("BER/BER-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_6 <- read.delim("BER/BER-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_7 <- read.delim("BER/BER-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_8 <- read.delim("BER/BER-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_9 <- read.delim("BER/BER-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_10 <- read.delim("BER/BER-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_11 <- read.delim("BER/BER-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_12 <- read.delim("BER/BER-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_13 <- read.delim("BER/BER-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_14 <- read.delim("BER/BER-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_15 <- read.delim("BER/BER-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_16 <- read.delim("BER/BER-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_17 <- read.delim("BER/BER-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_18 <- read.delim("BER/BER-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_19 <- read.delim("BER/BER-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_20 <- read.delim("BER/BER-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BER_21 <- read.delim("BER/BER-021/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)


otu_1 <- read.delim("BER/BER-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("BER/BER-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("BER/BER-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("BER/BER-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("BER/BER-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("BER/BER-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("BER/BER-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("BER/BER-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("BER/BER-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("BER/BER-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("BER/BER-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("BER/BER-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("BER/BER-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("BER/BER-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("BER/BER-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("BER/BER-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("BER/BER-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("BER/BER-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("BER/BER-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("BER/BER-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_21 <- read.delim("BER/BER-021/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)

for(i in 1:21){
  eval(parse(text=paste0('colnames(BER_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('BER_',i,'<-merge(BER_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

BER_species<-NULL
for(i in 1:21){
  bb<-eval(parse(text=paste0('unique(BER_',i,'$Species)')))
  BER_species<-union(BER_species,bb)
}
BER_species<-na.omit(BER_species)

BER_family<-NULL
for(i in 1:21){
  bb<-eval(parse(text=paste0('unique(BER_',i,'$Family)')))
  BER_family<-union(BER_family,bb)
}
BER_family<-na.omit(BER_family)

BER_order<-NULL
for(i in 1:21){
  bb<-eval(parse(text=paste0('unique(BER_',i,'$Order)')))
  BER_order<-union(BER_order,bb)
}
BER_order<-na.omit(BER_order)

rm(i,j,bb,aa)

#############################################################################################################################
################################################ BOG: 15 samples#############################################################

BOG_1 <- read.delim("BOG/BOG-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_2 <- read.delim("BOG/BOG-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_3 <- read.delim("BOG/BOG-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_4 <- read.delim("BOG/BOG-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_5 <- read.delim("BOG/BOG-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_6 <- read.delim("BOG/BOG-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_7 <- read.delim("BOG/BOG-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_8 <- read.delim("BOG/BOG-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_9 <- read.delim("BOG/BOG-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_10 <- read.delim("BOG/BOG-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_11 <- read.delim("BOG/BOG-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_12 <- read.delim("BOG/BOG-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_13 <- read.delim("BOG/BOG-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_14 <- read.delim("BOG/BOG-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
BOG_15 <- read.delim("BOG/BOG-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)

otu_1 <- read.delim("BOG/BOG-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("BOG/BOG-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("BOG/BOG-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("BOG/BOG-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("BOG/BOG-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("BOG/BOG-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("BOG/BOG-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("BOG/BOG-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("BOG/BOG-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("BOG/BOG-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("BOG/BOG-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("BOG/BOG-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("BOG/BOG-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("BOG/BOG-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("BOG/BOG-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)

for(j in 1:15){
  bb<-data.frame()
  otu_table<-eval(parse(text=paste0('BOG_',j,'[c(-1,-2),]')))
  uniq<-unique(otu_table$V2)
  uniq<-uniq[uniq!='']
  for(i in uniq){
    bb[which(uniq==i),1]<-i
    bb[which(uniq==i),2]<-sum(otu_table[which(otu_table$V2==i),'V3'])
  }
  eval(parse(text=paste0('BOG_',j,'<-bb')))
}

for(j in 1:15){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  bb<-otu_info[!duplicated(otu_info$V2),-1]
  eval(parse(text=paste0('otu_',j,'<-bb')))
}
rm(otu_info,otu_table,uniq)

for(i in c(1:15)){
  eval(parse(text=paste0('colnames(BOG_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('BOG_',i,'<-merge(BOG_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

BOG_species<-NULL
for(i in c(1:15)){
  bb<-eval(parse(text=paste0('unique(BOG_',i,'$Species)')))
  BOG_species<-union(BOG_species,bb)
}
BOG_species<-na.omit(BOG_species)

BOG_family<-NULL
for(i in c(1:15)){
  bb<-eval(parse(text=paste0('unique(BOG_',i,'$Family)')))
  BOG_family<-union(BOG_family,bb)
}
BOG_family<-na.omit(BOG_family)

BOG_order<-NULL
for(i in c(1:15)){
  bb<-eval(parse(text=paste0('unique(BOG_',i,'$Order)')))
  BOG_order<-union(BOG_order,bb)
}
BOG_order<-na.omit(BOG_order)

rm(aa,bb,cc,i,j)


#############################################################################################################################
################################################ HAM: 16 samples#############################################################

HAM_1 <- read.delim("HAM/HAM-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_2 <- read.delim("HAM/HAM-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_3 <- read.delim("HAM/HAM-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_4 <- read.delim("HAM/HAM-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_5 <- read.delim("HAM/HAM-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_6 <- read.delim("HAM/HAM-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_7 <- read.delim("HAM/HAM-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_8 <- read.delim("HAM/HAM-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_9 <- read.delim("HAM/HAM-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_10 <- read.delim("HAM/HAM-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_11 <- read.delim("HAM/HAM-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_12 <- read.delim("HAM/HAM-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_13 <- read.delim("HAM/HAM-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_14 <- read.delim("HAM/HAM-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_15 <- read.delim("HAM/HAM-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HAM_16 <- read.delim("HAM/HAM-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)

otu_1 <- read.delim("HAM/HAM-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("HAM/HAM-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("HAM/HAM-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("HAM/HAM-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("HAM/HAM-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("HAM/HAM-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("HAM/HAM-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("HAM/HAM-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("HAM/HAM-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("HAM/HAM-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("HAM/HAM-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("HAM/HAM-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("HAM/HAM-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("HAM/HAM-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("HAM/HAM-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("HAM/HAM-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)


for(i in c(1:16)){
  eval(parse(text=paste0('colnames(HAM_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('HAM_',i,'<-merge(HAM_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

HAM_species<-NULL
for(i in c(1:16)){
  bb<-eval(parse(text=paste0('unique(HAM_',i,'$Species)')))
  HAM_species<-union(HAM_species,bb)
}
HAM_species<-na.omit(HAM_species)

HAM_family<-NULL
for(i in c(1:16)){
  bb<-eval(parse(text=paste0('unique(HAM_',i,'$Family)')))
  HAM_family<-union(HAM_family,bb)
}
HAM_family<-na.omit(HAM_family)

HAM_order<-NULL
for(i in c(1:16)){
  bb<-eval(parse(text=paste0('unique(HAM_',i,'$Order)')))
  HAM_order<-union(HAM_order,bb)
}
HAM_order<-na.omit(HAM_order)

rm(aa,bb,i,j)

#############################################################################################################################
################################################ HGK: 1:8,10:18 samples######################################################
HGK_1 <- read.delim("HGK/HGK-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_2 <- read.delim("HGK/HGK-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_3 <- read.delim("HGK/HGK-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_4 <- read.delim("HGK/HGK-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_5 <- read.delim("HGK/HGK-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_6 <- read.delim("HGK/HGK-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_7 <- read.delim("HGK/HGK-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_8 <- read.delim("HGK/HGK-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
#HGK_9 <- read.delim("HGK/HGK-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_10 <- read.delim("HGK/HGK-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_11 <- read.delim("HGK/HGK-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_12 <- read.delim("HGK/HGK-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_13 <- read.delim("HGK/HGK-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_14 <- read.delim("HGK/HGK-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_15 <- read.delim("HGK/HGK-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_16 <- read.delim("HGK/HGK-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_17 <- read.delim("HGK/HGK-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
HGK_18 <- read.delim("HGK/HGK-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)

otu_1 <- read.delim("HGK/HGK-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("HGK/HGK-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("HGK/HGK-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("HGK/HGK-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("HGK/HGK-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("HGK/HGK-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("HGK/HGK-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("HGK/HGK-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
#otu_9 <- read.delim("HGK/HGK-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("HGK/HGK-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("HGK/HGK-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("HGK/HGK-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("HGK/HGK-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("HGK/HGK-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("HGK/HGK-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("HGK/HGK-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("HGK/HGK-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("HGK/HGK-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)


for(i in c(1:8,10:18)){
  eval(parse(text=paste0('colnames(HGK_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('HGK_',i,'<-merge(HGK_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

HGK_species<-NULL
for(i in c(1:8,10:18)){
  bb<-eval(parse(text=paste0('unique(HGK_',i,'$Species)')))
  HGK_species<-union(HGK_species,bb)
}
HGK_species<-na.omit(HGK_species)

HGK_family<-NULL
for(i in c(1:8,10:18)){
  bb<-eval(parse(text=paste0('unique(HGK_',i,'$Family)')))
  HGK_family<-union(HGK_family,bb)
}
HGK_family<-na.omit(HGK_family)

HGK_order<-NULL
for(i in c(1:8,10:18)){
  bb<-eval(parse(text=paste0('unique(HGK_',i,'$Order)')))
  HGK_order<-union(HGK_order,bb)
}
HGK_order<-na.omit(HGK_order)

rm(aa,bb,i,j)

#############################################################################################################################
################################################ ILR: 24 samples#############################################################

ILR_1 <- read.delim("ILR/ILR-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_2 <- read.delim("ILR/ILR-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_3 <- read.delim("ILR/ILR-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_4 <- read.delim("ILR/ILR-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_5 <- read.delim("ILR/ILR-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_6 <- read.delim("ILR/ILR-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_7 <- read.delim("ILR/ILR-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_8 <- read.delim("ILR/ILR-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_9 <- read.delim("ILR/ILR-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_10 <- read.delim("ILR/ILR-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_11 <- read.delim("ILR/ILR-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_12 <- read.delim("ILR/ILR-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_13 <- read.delim("ILR/ILR-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_14 <- read.delim("ILR/ILR-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_15 <- read.delim("ILR/ILR-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_16 <- read.delim("ILR/ILR-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_17 <- read.delim("ILR/ILR-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_18 <- read.delim("ILR/ILR-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_19 <- read.delim("ILR/ILR-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_20 <- read.delim("ILR/ILR-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_21 <- read.delim("ILR/ILR-021/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_22 <- read.delim("ILR/ILR-022/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_23 <- read.delim("ILR/ILR-023/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
ILR_24 <- read.delim("ILR/ILR-024/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)


otu_1 <- read.delim("ILR/ILR-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("ILR/ILR-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("ILR/ILR-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("ILR/ILR-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("ILR/ILR-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("ILR/ILR-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("ILR/ILR-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("ILR/ILR-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("ILR/ILR-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("ILR/ILR-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("ILR/ILR-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("ILR/ILR-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("ILR/ILR-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("ILR/ILR-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("ILR/ILR-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("ILR/ILR-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("ILR/ILR-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("ILR/ILR-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("ILR/ILR-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("ILR/ILR-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_21 <- read.delim("ILR/ILR-021/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_22 <- read.delim("ILR/ILR-022/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_23 <- read.delim("ILR/ILR-023/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_24 <- read.delim("ILR/ILR-024/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)

for(i in 1:24){
  eval(parse(text=paste0('colnames(ILR_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('ILR_',i,'<-merge(ILR_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

ILR_species<-NULL
for(i in 1:24){
  bb<-eval(parse(text=paste0('unique(ILR_',i,'$Species)')))
  ILR_species<-union(ILR_species,bb)
}
ILR_species<-na.omit(ILR_species)

ILR_family<-NULL
for(i in 1:24){
  bb<-eval(parse(text=paste0('unique(ILR_',i,'$Family)')))
  ILR_family<-union(ILR_family,bb)
}
ILR_family<-na.omit(ILR_family)

ILR_order<-NULL
for(i in 1:24){
  bb<-eval(parse(text=paste0('unique(ILR_',i,'$Order)')))
  ILR_order<-union(ILR_order,bb)
}
ILR_order<-na.omit(ILR_order)

rm(aa,bb,i,j)


#############################################################################################################################
################################################ LON: 1,3:6,8:24 samples#####################################################

LON_1 <- read.delim("LON/LON-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_3 <- read.delim("LON/LON-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_4 <- read.delim("LON/LON-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_5 <- read.delim("LON/LON-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_6 <- read.delim("LON/LON-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_8 <- read.delim("LON/LON-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_9 <- read.delim("LON/LON-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_10 <- read.delim("LON/LON-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_11 <- read.delim("LON/LON-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_12 <- read.delim("LON/LON-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_13 <- read.delim("LON/LON-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_14 <- read.delim("LON/LON-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_15 <- read.delim("LON/LON-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_16 <- read.delim("LON/LON-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_17 <- read.delim("LON/LON-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_18 <- read.delim("LON/LON-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_19 <- read.delim("LON/LON-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_20 <- read.delim("LON/LON-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_21 <- read.delim("LON/LON-021/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_22 <- read.delim("LON/LON-022/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_23 <- read.delim("LON/LON-023/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
LON_24 <- read.delim("LON/LON-024/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)


otu_1 <- read.delim("LON/LON-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("LON/LON-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("LON/LON-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("LON/LON-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("LON/LON-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("LON/LON-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("LON/LON-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("LON/LON-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("LON/LON-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("LON/LON-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("LON/LON-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("LON/LON-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("LON/LON-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("LON/LON-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("LON/LON-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("LON/LON-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("LON/LON-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("LON/LON-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_21 <- read.delim("LON/LON-021/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_22 <- read.delim("LON/LON-022/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_23 <- read.delim("LON/LON-023/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_24 <- read.delim("LON/LON-024/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)

for(i in c(1,3:6,8:24)){
  eval(parse(text=paste0('colnames(LON_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('LON_',i,'<-merge(LON_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

LON_species<-NULL
for(i in c(1,3:6,8:24)){
  bb<-eval(parse(text=paste0('unique(LON_',i,'$Species)')))
  LON_species<-union(LON_species,bb)
}
LON_species<-na.omit(LON_species)

LON_family<-NULL
for(i in c(1,3:6,8:24)){
  bb<-eval(parse(text=paste0('unique(LON_',i,'$Family)')))
  LON_family<-union(LON_family,bb)
}
LON_family<-na.omit(LON_family)

LON_order<-NULL
for(i in c(1,3:6,8:24)){
  bb<-eval(parse(text=paste0('unique(LON_',i,'$Order)')))
  LON_order<-union(LON_order,bb)
}
LON_order<-na.omit(LON_order)

rm(aa,bb,i,j)

#############################################################################################################################
################################################ MAR: 10 samples#############################################################

MAR_1 <- read.delim("MAR/MAR-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_2 <- read.delim("MAR/MAR-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_3 <- read.delim("MAR/MAR-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_4 <- read.delim("MAR/MAR-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_5 <- read.delim("MAR/MAR-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_6 <- read.delim("MAR/MAR-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_7 <- read.delim("MAR/MAR-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_8 <- read.delim("MAR/MAR-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_9 <- read.delim("MAR/MAR-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
MAR_10 <- read.delim("MAR/MAR-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)

otu_1 <- read.delim("MAR/MAR-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("MAR/MAR-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("MAR/MAR-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("MAR/MAR-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("MAR/MAR-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("MAR/MAR-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("MAR/MAR-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("MAR/MAR-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("MAR/MAR-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("MAR/MAR-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)


for(i in 1:10){
  eval(parse(text=paste0('colnames(MAR_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('MAR_',i,'<-merge(MAR_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

MAR_species<-NULL
for(i in 1:10){
  bb<-eval(parse(text=paste0('unique(MAR_',i,'$Species)')))
  MAR_species<-union(MAR_species,bb)
}
MAR_species<-na.omit(MAR_species)

MAR_family<-NULL
for(i in 1:10){
  bb<-eval(parse(text=paste0('unique(MAR_',i,'$Family)')))
  MAR_family<-union(MAR_family,bb)
}
MAR_family<-na.omit(MAR_family)

MAR_order<-NULL
for(i in 1:10){
  bb<-eval(parse(text=paste0('unique(MAR_',i,'$Order)')))
  MAR_order<-union(MAR_order,bb)
}
MAR_order<-na.omit(MAR_order)

rm(aa,bb,i,j)

#############################################################################################################################
################################################ NYC: 26 samples#############################################################

NYC_1 <- read.delim("NYC/NYC-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_2 <- read.delim("NYC/NYC-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_3 <- read.delim("NYC/NYC-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_4 <- read.delim("NYC/NYC-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_5 <- read.delim("NYC/NYC-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_6 <- read.delim("NYC/NYC-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_7 <- read.delim("NYC/NYC-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_8 <- read.delim("NYC/NYC-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_9 <- read.delim("NYC/NYC-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_10 <- read.delim("NYC/NYC-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_11 <- read.delim("NYC/NYC-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_12 <- read.delim("NYC/NYC-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_13 <- read.delim("NYC/NYC-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_14 <- read.delim("NYC/NYC-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_15 <- read.delim("NYC/NYC-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_16 <- read.delim("NYC/NYC-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_17 <- read.delim("NYC/NYC-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_18 <- read.delim("NYC/NYC-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_19 <- read.delim("NYC/NYC-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_20 <- read.delim("NYC/NYC-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_21 <- read.delim("NYC/NYC-021/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_22 <- read.delim("NYC/NYC-022/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_23 <- read.delim("NYC/NYC-023/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_24 <- read.delim("NYC/NYC-024/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_25 <- read.delim("NYC/NYC-025/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
NYC_26 <- read.delim("NYC/NYC-026/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)

otu_1 <- read.delim("NYC/NYC-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("NYC/NYC-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("NYC/NYC-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("NYC/NYC-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("NYC/NYC-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("NYC/NYC-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("NYC/NYC-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("NYC/NYC-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("NYC/NYC-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("NYC/NYC-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("NYC/NYC-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("NYC/NYC-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("NYC/NYC-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("NYC/NYC-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("NYC/NYC-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("NYC/NYC-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("NYC/NYC-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("NYC/NYC-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("NYC/NYC-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("NYC/NYC-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_21 <- read.delim("NYC/NYC-021/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_22 <- read.delim("NYC/NYC-022/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_23 <- read.delim("NYC/NYC-023/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_24 <- read.delim("NYC/NYC-024/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_25 <- read.delim("NYC/NYC-025/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_26 <- read.delim("NYC/NYC-026/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)


for(j in 1:26){
  bb<-data.frame()
  otu_table<-eval(parse(text=paste0('NYC_',j,'[c(-1,-2),]')))
  uniq<-unique(otu_table$V2)
  uniq<-uniq[uniq!='']
  for(i in uniq){
    bb[which(uniq==i),1]<-i
    bb[which(uniq==i),2]<-sum(otu_table[which(otu_table$V2==i),'V3'])
  }
  eval(parse(text=paste0('NYC_',j,'<-bb')))
}

for(j in 1:26){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  bb<-otu_info[!duplicated(otu_info$V2),-1]
  eval(parse(text=paste0('otu_',j,'<-bb')))
}
rm(otu_info,otu_table,uniq)

for(i in c(1:26)){
  eval(parse(text=paste0('colnames(NYC_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('NYC_',i,'<-merge(NYC_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

NYC_species<-NULL
for(i in c(1:26)){
  bb<-eval(parse(text=paste0('unique(NYC_',i,'$Species)')))
  NYC_species<-union(NYC_species,bb)
}
NYC_species<-na.omit(NYC_species)

NYC_family<-NULL
for(i in c(1:26)){
  bb<-eval(parse(text=paste0('unique(NYC_',i,'$Family)')))
  NYC_family<-union(NYC_family,bb)
}
NYC_family<-na.omit(NYC_family)

NYC_order<-NULL
for(i in c(1:26)){
  bb<-eval(parse(text=paste0('unique(NYC_',i,'$Order)')))
  NYC_order<-union(NYC_order,bb)
}
NYC_order<-na.omit(NYC_order)


rm(aa,bb,i,j)

#############################################################################################################################
################################################ OFA: 20 samples#############################################################

OFA_1 <- read.delim("OFA/OFA-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_2 <- read.delim("OFA/OFA-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_3 <- read.delim("OFA/OFA-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_4 <- read.delim("OFA/OFA-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_5 <- read.delim("OFA/OFA-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_6 <- read.delim("OFA/OFA-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_7 <- read.delim("OFA/OFA-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_8 <- read.delim("OFA/OFA-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_9 <- read.delim("OFA/OFA-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_10 <- read.delim("OFA/OFA-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_11 <- read.delim("OFA/OFA-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_12 <- read.delim("OFA/OFA-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_13 <- read.delim("OFA/OFA-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_14 <- read.delim("OFA/OFA-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_15 <- read.delim("OFA/OFA-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_16 <- read.delim("OFA/OFA-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_17 <- read.delim("OFA/OFA-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_18 <- read.delim("OFA/OFA-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_19 <- read.delim("OFA/OFA-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
OFA_20 <- read.delim("OFA/OFA-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)


otu_1 <- read.delim("OFA/OFA-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("OFA/OFA-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("OFA/OFA-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("OFA/OFA-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("OFA/OFA-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("OFA/OFA-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("OFA/OFA-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("OFA/OFA-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("OFA/OFA-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("OFA/OFA-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("OFA/OFA-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("OFA/OFA-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("OFA/OFA-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("OFA/OFA-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("OFA/OFA-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("OFA/OFA-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("OFA/OFA-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("OFA/OFA-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("OFA/OFA-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("OFA/OFA-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)



for(j in 1:20){
  bb<-data.frame()
  otu_table<-eval(parse(text=paste0('OFA_',j,'[c(-1,-2),]')))
  uniq<-unique(otu_table$V2)
  uniq<-uniq[uniq!='']
  for(i in uniq){
    bb[which(uniq==i),1]<-i
    bb[which(uniq==i),2]<-sum(otu_table[which(otu_table$V2==i),'V3'])
  }
  eval(parse(text=paste0('OFA_',j,'<-bb')))
}

for(j in 1:20){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  bb<-otu_info[!duplicated(otu_info$V2),-1]
  eval(parse(text=paste0('otu_',j,'<-bb')))
}
rm(otu_info,otu_table,uniq)

for(i in c(1:20)){
  eval(parse(text=paste0('colnames(OFA_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('OFA_',i,'<-merge(OFA_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

OFA_species<-NULL
for(i in c(1:20)){
  bb<-eval(parse(text=paste0('unique(OFA_',i,'$Species)')))
  OFA_species<-union(OFA_species,bb)
}
OFA_species<-na.omit(OFA_species)

OFA_family<-NULL
for(i in c(1:20)){
  bb<-eval(parse(text=paste0('unique(OFA_',i,'$Family)')))
  OFA_family<-union(OFA_family,bb)
}
OFA_family<-na.omit(OFA_family)

OFA_order<-NULL
for(i in c(1:20)){
  bb<-eval(parse(text=paste0('unique(OFA_',i,'$Order)')))
  OFA_order<-union(OFA_order,bb)
}
OFA_order<-na.omit(OFA_order)


rm(aa,bb,i,j)

#############################################################################################################################
################################################ PXO: 20 samples#############################################################

PXO_1 <- read.delim("PXO/PXO-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_2 <- read.delim("PXO/PXO-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_3 <- read.delim("PXO/PXO-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_4 <- read.delim("PXO/PXO-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_5 <- read.delim("PXO/PXO-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_6 <- read.delim("PXO/PXO-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_7 <- read.delim("PXO/PXO-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_8 <- read.delim("PXO/PXO-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_9 <- read.delim("PXO/PXO-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_10 <- read.delim("PXO/PXO-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_11 <- read.delim("PXO/PXO-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_12 <- read.delim("PXO/PXO-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_13 <- read.delim("PXO/PXO-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_14 <- read.delim("PXO/PXO-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_15 <- read.delim("PXO/PXO-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_16 <- read.delim("PXO/PXO-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_17 <- read.delim("PXO/PXO-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_18 <- read.delim("PXO/PXO-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_19 <- read.delim("PXO/PXO-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
PXO_20 <- read.delim("PXO/PXO-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)

otu_1 <- read.delim("PXO/PXO-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("PXO/PXO-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("PXO/PXO-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("PXO/PXO-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("PXO/PXO-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("PXO/PXO-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("PXO/PXO-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("PXO/PXO-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("PXO/PXO-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("PXO/PXO-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("PXO/PXO-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("PXO/PXO-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("PXO/PXO-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("PXO/PXO-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("PXO/PXO-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("PXO/PXO-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("PXO/PXO-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("PXO/PXO-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("PXO/PXO-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("PXO/PXO-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)


for(i in 1:20){
  eval(parse(text=paste0('colnames(PXO_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('PXO_',i,'<-merge(PXO_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

PXO_species<-NULL
for(i in 1:20){
  bb<-eval(parse(text=paste0('unique(PXO_',i,'$Species)')))
  PXO_species<-union(PXO_species,bb)
}
PXO_species<-na.omit(PXO_species)

PXO_family<-NULL
for(i in 1:20){
  bb<-eval(parse(text=paste0('unique(PXO_',i,'$Family)')))
  PXO_family<-union(PXO_family,bb)
}
PXO_family<-na.omit(PXO_family)

PXO_order<-NULL
for(i in 1:20){
  bb<-eval(parse(text=paste0('unique(PXO_',i,'$Order)')))
  PXO_order<-union(PXO_order,bb)
}
PXO_order<-na.omit(PXO_order)

rm(aa,bb,i,j)

#############################################################################################################################
################################################ SAC: 18 samples#############################################################

SAC_1 <- read.delim("SAC/SAC-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_2 <- read.delim("SAC/SAC-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_3 <- read.delim("SAC/SAC-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_4 <- read.delim("SAC/SAC-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_5 <- read.delim("SAC/SAC-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_6 <- read.delim("SAC/SAC-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_7 <- read.delim("SAC/SAC-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_8 <- read.delim("SAC/SAC-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_9 <- read.delim("SAC/SAC-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_10 <- read.delim("SAC/SAC-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_11 <- read.delim("SAC/SAC-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_12 <- read.delim("SAC/SAC-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_13 <- read.delim("SAC/SAC-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_14 <- read.delim("SAC/SAC-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_15 <- read.delim("SAC/SAC-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_16 <- read.delim("SAC/SAC-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_17 <- read.delim("SAC/SAC-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAC_18 <- read.delim("SAC/SAC-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)


otu_1 <- read.delim("SAC/SAC-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("SAC/SAC-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("SAC/SAC-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("SAC/SAC-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("SAC/SAC-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("SAC/SAC-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("SAC/SAC-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("SAC/SAC-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("SAC/SAC-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("SAC/SAC-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("SAC/SAC-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("SAC/SAC-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("SAC/SAC-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("SAC/SAC-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("SAC/SAC-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("SAC/SAC-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("SAC/SAC-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("SAC/SAC-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)



for(j in 1:18){
  bb<-data.frame()
  otu_table<-eval(parse(text=paste0('SAC_',j,'[c(-1,-2),]')))
  uniq<-unique(otu_table$V2)
  uniq<-uniq[uniq!='']
  for(i in uniq){
    bb[which(uniq==i),1]<-i
    bb[which(uniq==i),2]<-sum(otu_table[which(otu_table$V2==i),'V3'])
  }
  eval(parse(text=paste0('SAC_',j,'<-bb')))
}

for(j in 1:18){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  bb<-otu_info[!duplicated(otu_info$V2),-1]
  eval(parse(text=paste0('otu_',j,'<-bb')))
}
rm(otu_info,otu_table,uniq)

for(i in c(1:18)){
  eval(parse(text=paste0('colnames(SAC_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('SAC_',i,'<-merge(SAC_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

SAC_species<-NULL
for(i in c(1:18)){
  bb<-eval(parse(text=paste0('unique(SAC_',i,'$Species)')))
  SAC_species<-union(SAC_species,bb)
}
SAC_species<-na.omit(SAC_species)

SAC_family<-NULL
for(i in c(1:18)){
  bb<-eval(parse(text=paste0('unique(SAC_',i,'$Family)')))
  SAC_family<-union(SAC_family,bb)
}
SAC_family<-na.omit(SAC_family)

SAC_order<-NULL
for(i in c(1:18)){
  bb<-eval(parse(text=paste0('unique(SAC_',i,'$Order)')))
  SAC_order<-union(SAC_order,bb)
}
SAC_order<-na.omit(SAC_order)


rm(aa,bb,i,j)

#############################################################################################################################
################################################ SAO: 24 samples#############################################################

SAO_1 <- read.delim("SAO/SAO-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_2 <- read.delim("SAO/SAO-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_3 <- read.delim("SAO/SAO-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_4 <- read.delim("SAO/SAO-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_5 <- read.delim("SAO/SAO-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_6 <- read.delim("SAO/SAO-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_7 <- read.delim("SAO/SAO-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_8 <- read.delim("SAO/SAO-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_9 <- read.delim("SAO/SAO-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_10 <- read.delim("SAO/SAO-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_11 <- read.delim("SAO/SAO-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_12 <- read.delim("SAO/SAO-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_13 <- read.delim("SAO/SAO-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_14 <- read.delim("SAO/SAO-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_15 <- read.delim("SAO/SAO-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_16 <- read.delim("SAO/SAO-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_17 <- read.delim("SAO/SAO-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_18 <- read.delim("SAO/SAO-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_19 <- read.delim("SAO/SAO-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_20 <- read.delim("SAO/SAO-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_21 <- read.delim("SAO/SAO-021/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_22 <- read.delim("SAO/SAO-022/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_23 <- read.delim("SAO/SAO-023/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SAO_24 <- read.delim("SAO/SAO-024/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)


otu_1 <- read.delim("SAO/SAO-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("SAO/SAO-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("SAO/SAO-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("SAO/SAO-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("SAO/SAO-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("SAO/SAO-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("SAO/SAO-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("SAO/SAO-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("SAO/SAO-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("SAO/SAO-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("SAO/SAO-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("SAO/SAO-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("SAO/SAO-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("SAO/SAO-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("SAO/SAO-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("SAO/SAO-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("SAO/SAO-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("SAO/SAO-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("SAO/SAO-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("SAO/SAO-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_21 <- read.delim("SAO/SAO-021/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_22 <- read.delim("SAO/SAO-022/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_23 <- read.delim("SAO/SAO-023/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_24 <- read.delim("SAO/SAO-024/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)

for(i in 1:24){
  eval(parse(text=paste0('colnames(SAO_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('SAO_',i,'<-merge(SAO_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

SAO_species<-NULL
for(i in 1:24){
  bb<-eval(parse(text=paste0('unique(SAO_',i,'$Species)')))
  SAO_species<-union(SAO_species,bb)
}
SAO_species<-na.omit(SAO_species)

SAO_family<-NULL
for(i in 1:24){
  bb<-eval(parse(text=paste0('unique(SAO_',i,'$Family)')))
  SAO_family<-union(SAO_family,bb)
}
SAO_family<-na.omit(SAO_family)

SAO_order<-NULL
for(i in 1:24){
  bb<-eval(parse(text=paste0('unique(SAO_',i,'$Order)')))
  SAO_order<-union(SAO_order,bb)
}
SAO_order<-na.omit(SAO_order)

rm(aa,bb,i,j)

#############################################################################################################################
################################################ SOF: 10 samples#############################################################

SOF_1 <- read.delim("SOF/SOF-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_2 <- read.delim("SOF/SOF-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_3 <- read.delim("SOF/SOF-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_4 <- read.delim("SOF/SOF-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_5 <- read.delim("SOF/SOF-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_6 <- read.delim("SOF/SOF-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_7 <- read.delim("SOF/SOF-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_8 <- read.delim("SOF/SOF-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_9 <- read.delim("SOF/SOF-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
SOF_10 <- read.delim("SOF/SOF-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)


otu_1 <- read.delim("SOF/SOF-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("SOF/SOF-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("SOF/SOF-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("SOF/SOF-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("SOF/SOF-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("SOF/SOF-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("SOF/SOF-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("SOF/SOF-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("SOF/SOF-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("SOF/SOF-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)

for(i in 1:10){
  eval(parse(text=paste0('colnames(SOF_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('SOF_',i,'<-merge(SOF_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

SOF_species<-NULL
for(i in 1:10){
  bb<-eval(parse(text=paste0('unique(SOF_',i,'$Species)')))
  SOF_species<-union(SOF_species,bb)
}
SOF_species<-na.omit(SOF_species)

SOF_family<-NULL
for(i in 1:10){
  bb<-eval(parse(text=paste0('unique(SOF_',i,'$Family)')))
  SOF_family<-union(SOF_family,bb)
}
SOF_family<-na.omit(SOF_family)

SOF_order<-NULL
for(i in 1:10){
  bb<-eval(parse(text=paste0('unique(SOF_',i,'$Order)')))
  SOF_order<-union(SOF_order,bb)
}
SOF_order<-na.omit(SOF_order)

rm(aa,bb,i,j)

#############################################################################################################################
################################################ STO: 20 samples#############################################################

STO_1 <- read.delim("STO/STO-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_2 <- read.delim("STO/STO-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_3 <- read.delim("STO/STO-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_4 <- read.delim("STO/STO-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_5 <- read.delim("STO/STO-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_6 <- read.delim("STO/STO-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_7 <- read.delim("STO/STO-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_8 <- read.delim("STO/STO-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_9 <- read.delim("STO/STO-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_10 <- read.delim("STO/STO-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_11 <- read.delim("STO/STO-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_12 <- read.delim("STO/STO-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_13 <- read.delim("STO/STO-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_14 <- read.delim("STO/STO-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_15 <- read.delim("STO/STO-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_16 <- read.delim("STO/STO-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_17 <- read.delim("STO/STO-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_18 <- read.delim("STO/STO-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_19 <- read.delim("STO/STO-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
STO_20 <- read.delim("STO/STO-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)


otu_1 <- read.delim("STO/STO-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("STO/STO-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("STO/STO-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("STO/STO-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("STO/STO-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("STO/STO-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("STO/STO-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("STO/STO-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("STO/STO-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("STO/STO-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("STO/STO-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("STO/STO-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("STO/STO-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("STO/STO-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("STO/STO-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("STO/STO-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("STO/STO-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("STO/STO-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("STO/STO-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("STO/STO-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)

for(i in 1:20){
  eval(parse(text=paste0('colnames(STO_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('STO_',i,'<-merge(STO_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

STO_species<-NULL
for(i in 1:20){
  bb<-eval(parse(text=paste0('unique(STO_',i,'$Species)')))
  STO_species<-union(STO_species,bb)
}
STO_species<-na.omit(STO_species)

STO_family<-NULL
for(i in 1:20){
  bb<-eval(parse(text=paste0('unique(STO_',i,'$Family)')))
  STO_family<-union(STO_family,bb)
}
STO_family<-na.omit(STO_family)

STO_order<-NULL
for(i in 1:20){
  bb<-eval(parse(text=paste0('unique(STO_',i,'$Order)')))
  STO_order<-union(STO_order,bb)
}
STO_order<-na.omit(STO_order)

rm(aa,bb,i,j)

#############################################################################################################################
################################################ TOK: 25 samples#############################################################

TOK_1 <- read.delim("TOK/TOK-001/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_2 <- read.delim("TOK/TOK-002/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_3 <- read.delim("TOK/TOK-003/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_4 <- read.delim("TOK/TOK-004/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_5 <- read.delim("TOK/TOK-005/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_6 <- read.delim("TOK/TOK-006/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_7 <- read.delim("TOK/TOK-007/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_8 <- read.delim("TOK/TOK-008/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_9 <- read.delim("TOK/TOK-009/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_10 <- read.delim("TOK/TOK-010/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_11 <- read.delim("TOK/TOK-011/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_12 <- read.delim("TOK/TOK-012/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_13 <- read.delim("TOK/TOK-013/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_14 <- read.delim("TOK/TOK-014/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_15 <- read.delim("TOK/TOK-015/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_16 <- read.delim("TOK/TOK-016/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_17 <- read.delim("TOK/TOK-017/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_18 <- read.delim("TOK/TOK-018/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_19 <- read.delim("TOK/TOK-019/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_20 <- read.delim("TOK/TOK-020/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_21 <- read.delim("TOK/TOK-021/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_22 <- read.delim("TOK/TOK-022/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_23 <- read.delim("TOK/TOK-023/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_24 <- read.delim("TOK/TOK-024/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)
TOK_25 <- read.delim("TOK/TOK-025/wgs_from_biom_final.txt", header=FALSE, comment.char="#",stringsAsFactors=FALSE)

otu_1 <- read.delim("TOK/TOK-001/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_2 <- read.delim("TOK/TOK-002/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_3 <- read.delim("TOK/TOK-003/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_4 <- read.delim("TOK/TOK-004/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_5 <- read.delim("TOK/TOK-005/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_6 <- read.delim("TOK/TOK-006/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_7 <- read.delim("TOK/TOK-007/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_8 <- read.delim("TOK/TOK-008/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_9 <- read.delim("TOK/TOK-009/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_10 <- read.delim("TOK/TOK-010/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_11 <- read.delim("TOK/TOK-011/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_12 <- read.delim("TOK/TOK-012/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_13 <- read.delim("TOK/TOK-013/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_14 <- read.delim("TOK/TOK-014/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_15 <- read.delim("TOK/TOK-015/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_16 <- read.delim("TOK/TOK-016/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_17 <- read.delim("TOK/TOK-017/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_18 <- read.delim("TOK/TOK-018/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_19 <- read.delim("TOK/TOK-019/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_20 <- read.delim("TOK/TOK-020/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_21 <- read.delim("TOK/TOK-021/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_22 <- read.delim("TOK/TOK-022/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_23 <- read.delim("TOK/TOK-023/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_24 <- read.delim("TOK/TOK-024/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)
otu_25 <- read.delim("TOK/TOK-025/rep_set_aligned_tax_assignments.txt", header=FALSE,stringsAsFactors=FALSE)


for(j in 2:25){
  bb<-data.frame()
  otu_table<-eval(parse(text=paste0('TOK_',j,'[c(-1,-2),]')))
  uniq<-unique(otu_table$V2)
  uniq<-uniq[uniq!='']
  for(i in uniq){
    bb[which(uniq==i),1]<-i
    bb[which(uniq==i),2]<-sum(otu_table[which(otu_table$V2==i),'V3'])
  }
  eval(parse(text=paste0('TOK_',j,'<-bb')))
}

for(j in 2:25){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  bb<-otu_info[!duplicated(otu_info$V2),-1]
  eval(parse(text=paste0('otu_',j,'<-bb')))
}
rm(otu_info,otu_table,uniq)

for(i in c(1:25)){
  eval(parse(text=paste0('colnames(TOK_',i,')[1:2]<-c(\'id\',\'count\')')))
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
  eval(parse(text=paste0('TOK_',i,'<-merge(TOK_',i,',otu_',i,',by=\'id\')')))
  eval(parse(text=paste0('rm(otu_',i,')')))
}

TOK_species<-NULL
for(i in c(1:25)){
  bb<-eval(parse(text=paste0('unique(TOK_',i,'$Species)')))
  TOK_species<-union(TOK_species,bb)
}
TOK_species<-na.omit(TOK_species)

TOK_family<-NULL
for(i in c(1:25)){
  bb<-eval(parse(text=paste0('unique(TOK_',i,'$Family)')))
  TOK_family<-union(TOK_family,bb)
}
TOK_family<-na.omit(TOK_family)

TOK_order<-NULL
for(i in c(1:25)){
  bb<-eval(parse(text=paste0('unique(TOK_',i,'$Order)')))
  TOK_order<-union(TOK_order,bb)
}
TOK_order<-na.omit(TOK_order)

rm(aa,bb,i,j)



#############################################################################################################################
################################################ save the data ##############################################################

save.image(file='Main dataset 1 OTU tables for all samples.RData')














