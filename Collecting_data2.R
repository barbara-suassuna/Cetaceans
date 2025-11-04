#install.packages("taxize")
#load packages
library(taxize)
library(dplyr)
library(rentrez)
library(seqinr)

#Part I: 
#Find the TSN (Taxonomic Serial Number) for Cetacea within ITIS db
tax_cetacea <- get_ids("Cetacea", db="itis")
tax_cetacea$itis 
classification(tax_cetacea$itis)

#Get all species (list form) under Cetacea
cetacea_id <- 180403
parenttaxa <- downstream(cetacea_id, downto = "species", db = "itis") 
parenttaxa$`180403`$taxonname #all species scientific names

#Part II: repeat it for Outgroup data
tax_out <- get_ids("Hippopotamidae", db="itis")
tax_out$itis 
classification(tax_out$itis)
hippo_id <- 624917
outgroup <- downstream(hippo_id, downto = "species", db = "itis") 
outgroup$`624917`$taxonname #all species scientific names

#Merge the 2 lists
species_outgroup <- c(parenttaxa$`180403`$taxonname, outgroup$`624917`$taxonname)

#Part III: Collect the data for each gene

##1 - COX1 GENE
#Create a Table
whale<-matrix(,nrow=(length(species_outgroup)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
cox1_sequence<-character()

#LOOP
for(i in 1: length(unique(species_outgroup))){
  print(i)
  whale[i,1]<-species_outgroup[i]
  seqout<-entrez_search(db="nuccore", term=paste(species_outgroup[i], "[ORGN] AND 1:1600[SLEN] AND COX1[GENE]",
                                                 sep=""), retmax=1)
  if(length(seqout$ids)<1){
    whale[i,2]<-NA
    whale[i,3]<-NA
    whale[i,4]<-NA
  }
  
  else{
    seqout1<-entrez_summary(db="nuccore", id=seqout$ids)
    whale[i,2]<-seqout1$accessionversion
    whale[i,3]<-seqout1$title
    whale[i,4]<-seqout1$slen
    
    seqout<-entrez_fetch(db="nuccore", id=seqout1$accessionversion, rettype="fasta")
    seqout<-sub(">([^\n]*)", paste0(">", seqout1$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    cox1_sequence<-c(cox1_sequence, seqout)
  }
}

#Save FASTA file
write(cox1_sequence, file="cox1_gene.fasta")

#Save TABLE in a file
write.csv(whale, "cox1-info.csv")

##2 - MC1R GENE 
#Create a Table
whale<-matrix(,nrow=(length(parenttaxa$`180403`$taxonname)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
mc_sequence<-character()

for(i in 1: length(unique(species_outgroup))){
  print(i)
  whale[i,1]<-species_outgroup[i]
  seqout<-entrez_search(db="nuccore", term=paste(species_outgroup[i], "[ORGN] AND 1:4500[SLEN] AND 
                                                 MC1R[ALL]",
                                                 sep=""), retmax=1)
  if(length(seqout$ids)<1){
    whale[i,2]<-NA
    whale[i,3]<-NA
    whale[i,4]<-NA
  }
  
  else{
    seqout1<-entrez_summary(db="nuccore", id=seqout$ids)
    whale[i,2]<-seqout1$accessionversion
    whale[i,3]<-seqout1$title
    whale[i,4]<-seqout1$slen
    
    seqout<-entrez_fetch(db="nuccore", id=seqout1$accessionversion, rettype="fasta")
    seqout<-sub(">([^\n]*)", paste0(">", seqout1$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    mc_sequence<-c(mc_sequence, seqout)
  }
}

#Save FASTA file
write(mc_sequence, file="mc_gene.fasta")

#Save TABLE in a file
write.csv(whale, "mc-info.csv")

##3 - ENAM GENE
#Create a Table
whale<-matrix(,nrow=(length(species_outgroup)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
enam_sequence<-character()
i <- 1
#LOOP 
for(i in 1: length(unique(species_outgroup))){
  print(i)
  whale[i,1]<-species_outgroup[i]
  seqout<-entrez_search(db="nuccore", term=paste(species_outgroup[i], "[ORGN], AND ENAM[ALL] AND 1:6000[SLEN]",
                                                 sep=""), retmax=1)
  if(length(seqout$ids)<1){
    whale[i,2]<-NA
    whale[i,3]<-NA
    whale[i,4]<-NA
  }
  
  else{
    seqout1<-entrez_summary(db="nuccore", id=seqout$ids)
    whale[i,2]<-seqout1$accessionversion
    whale[i,3]<-seqout1$title
    whale[i,4]<-seqout1$slen
    
    seqout<-entrez_fetch(db="nuccore", id=seqout1$accessionversion, rettype="fasta")
    seqout<-sub(">([^\n]*)", paste0(">", seqout1$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    enam_sequence<-c(enam_sequence, seqout)
  }
}

#Save FASTA file
write(enam_sequence, file="enam_gene.fasta")

#Save TABLE in a file
write.csv(whale, "enam-info.csv")

##4 - CSN2 GENE
#Create a Table
whale<-matrix(,nrow=(length(species_outgroup)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
csn_sequence<-character()

#LOOP 
for(i in 1: length(unique(species_outgroup))){
  print(i)
  whale[i,1]<-species_outgroup[i]
  seqout<-entrez_search(db="nuccore", term=paste(species_outgroup[i], "[ORGN], AND CSN2[ALL] AND 1:1200[SLEN]",
                                                 sep=""), retmax=1)
  if(length(seqout$ids)<1){
    whale[i,2]<-NA
    whale[i,3]<-NA
    whale[i,4]<-NA
  }
  
  else{
    seqout1<-entrez_summary(db="nuccore", id=seqout$ids)
    whale[i,2]<-seqout1$accessionversion
    whale[i,3]<-seqout1$title
    whale[i,4]<-seqout1$slen
    
    seqout<-entrez_fetch(db="nuccore", id=seqout1$accessionversion, rettype="fasta")
    seqout<-sub(">([^\n]*)", paste0(">", seqout1$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    csn_sequence<-c(csn_sequence, seqout)
  }
}

#Save FASTA file
write(csn_sequence, file="csn_gene.fasta")

#Save TABLE in a file
write.csv(whale, "csn-info.csv")

##5 - PRM1 GENE
#Create a Table
whale<-matrix(,nrow=(length(species_outgroup)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
prm1_sequence<-character()

for(i in 1: length(unique(species_outgroup))){
  print(i)
  whale[i,1]<-species_outgroup[i]
  seqout<-entrez_search(db="nuccore", term=paste(species_outgroup[i], "[ORGN], AND PRM1[ALL]",
                                                 sep=""), retmax=1)
  if(length(seqout$ids)<1){
    whale[i,2]<-NA
    whale[i,3]<-NA
    whale[i,4]<-NA
  }
  
  else{
    seqout1<-entrez_summary(db="nuccore", id=seqout$ids)
    whale[i,2]<-seqout1$accessionversion
    whale[i,3]<-seqout1$title
    whale[i,4]<-seqout1$slen
    
    seqout<-entrez_fetch(db="nuccore", id=seqout1$accessionversion, rettype="fasta")
    seqout<-sub(">([^\n]*)", paste0(">", seqout1$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    prm1_sequence<-c(prm1_sequence, seqout)
  }
}

#Save FASTA file
write(prm1_sequence, file="prm1_gene.fasta")

#Save TABLE in a file
write.csv(whale, "prm1-info.csv")

##6 - MCPH1 GENE 
#Create a Table
whale<-matrix(,nrow=(length(species_outgroup)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
mcph_sequence<-character()

for(i in 1: length(unique(species_outgroup))){
  print(i)
  whale[i,1]<-species_outgroup[i]
  seqout<-entrez_search(db="nuccore", term=paste(species_outgroup[i], "[ORGN] AND 1:4500[SLEN] AND 
                                                 MCPH1[ALL]",
                                                 sep=""), retmax=1)
  if(length(seqout$ids)<1){
    whale[i,2]<-NA
    whale[i,3]<-NA
    whale[i,4]<-NA
  }
  
  else{
    seqout1<-entrez_summary(db="nuccore", id=seqout$ids)
    whale[i,2]<-seqout1$accessionversion
    whale[i,3]<-seqout1$title
    whale[i,4]<-seqout1$slen
    
    seqout<-entrez_fetch(db="nuccore", id=seqout1$accessionversion, rettype="fasta")
    seqout<-sub(">([^\n]*)", paste0(">", seqout1$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    mcph_sequence<-c(mcph_sequence, seqout)
  }
}

#Save FASTA file
write(mcph_sequence, file="mcph_gene.fasta")

#Save TABLE in a file
write.csv(whale, "mcph-info.csv")

##7 - ND2 GENE 
#Create a Table
whale<-matrix(,nrow=(length(species_outgroup)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
nd_sequence<-character()

for(i in 1: length(unique(species_outgroup))){
  print(i)
  whale[i,1]<-species_outgroup[i]
  seqout<-entrez_search(db="nuccore", term=paste(species_outgroup[i], "[ORGN] AND 1:4500[SLEN] AND 
                                                 ND2[ALL]",
                                                 sep=""), retmax=1)
  if(length(seqout$ids)<1){
    whale[i,2]<-NA
    whale[i,3]<-NA
    whale[i,4]<-NA
  }
  
  else{
    seqout1<-entrez_summary(db="nuccore", id=seqout$ids)
    whale[i,2]<-seqout1$accessionversion
    whale[i,3]<-seqout1$title
    whale[i,4]<-seqout1$slen
    
    seqout<-entrez_fetch(db="nuccore", id=seqout1$accessionversion, rettype="fasta")
    seqout<-sub(">([^\n]*)", paste0(">", seqout1$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    nd_sequence<-c(nd_sequence, seqout)
  }
}

#Save FASTA file
write(mcph_sequence, file="nd_gene.fasta")

#Save TABLE in a file
write.csv(whale, "nd-info.csv")





