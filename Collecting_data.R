#install.packages("worrms")
library(worrms)
library(taxize)
library(dplyr)
library(rentrez)
library(seqinr)
library(XML)
library(httr)

#Find the TSN (Taxonomic Serial Number) for Cetacea within ITIS db
tax_cetacea <- get_ids("Cetacea", db="itis")
tax_cetacea$itis 
classification(tax_cetacea$itis)

#Get all species (list form) under Cetacea
cetacea_id <- 180403
parenttaxa <- downstream(cetacea_id, downto = "species", db = "itis") 

parenttaxa$`180403`$rankname #confirms all have a "species" status
parenttaxa$`180403`$taxonname #all species scientific names

#1 - PRM1 GENE
#Create a Table
whale<-matrix(,nrow=(length(parenttaxa$`180403`$taxonname)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
prm1_sequence<-character()

#i<-1

#LOOP 
### We're referring to Scientific Names' List found in ITIS db (parenttaxa$`180403`$taxonname)
### Nuccore db will be used to get accnum, seqlen, and seqname

for(i in 1: length(unique(parenttaxa$`180403`$taxonname))){
  print(i)
  whale[i,1]<-parenttaxa$`180403`$taxonname[i]
  seqout<-entrez_search(db="nuccore", term=paste(parenttaxa$`180403`$taxonname[i], "[ORGN], AND PRM1[ALL]",
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
write.csv(whale, "prm1_gene_seq_info.csv")


#2- CYTB GENE
#Create a Table
whale<-matrix(,nrow=(length(parenttaxa$`180403`$taxonname)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
cytb_sequence<-character()

#LOOP
for(i in 1: length(unique(parenttaxa$`180403`$taxonname))){
  print(i)
  whale[i,1]<-parenttaxa$`180403`$taxonname[i]
  seqout<-entrez_search(db="nuccore", term=paste(parenttaxa$`180403`$taxonname[i], "[ORGN] AND 300:1300[SLEN] AND CYTB[GENE]",
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
    cytb_sequence<-c(cytb_sequence, seqout)
  }
}

#Save FASTA file
write(cytb_sequence, file="cytb_gene.fasta")

#Save TABLE in a file
write.csv(whale, "cytb_gene_seq_info.csv")

#3 - COX1 GENE
#Create a Table
whale<-matrix(,nrow=(length(parenttaxa$`180403`$taxonname)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
cox1_sequence<-character()

#LOOP
for(i in 1: length(unique(parenttaxa$`180403`$taxonname))){
  print(i)
  whale[i,1]<-parenttaxa$`180403`$taxonname[i]
  seqout<-entrez_search(db="nuccore", term=paste(parenttaxa$`180403`$taxonname[i], "[ORGN] AND 1:1600[SLEN] AND COX1[GENE]",
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
write.csv(whale, "cox1_gene_seq_info.csv")


#4 - MC1R GENE 
#Create a Table
whale<-matrix(,nrow=(length(parenttaxa$`180403`$taxonname)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
mc_sequence<-character()

for(i in 1: length(unique(parenttaxa$`180403`$taxonname))){
  print(i)
  whale[i,1]<-parenttaxa$`180403`$taxonname[i]
  seqout<-entrez_search(db="nuccore", term=paste(parenttaxa$`180403`$taxonname[i], "[ORGN] AND 1:4500[SLEN] AND 
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
write.csv(whale, "mc_gene_seq_info.csv")


#5 - MCPH1 GENE 
#Create a Table
whale<-matrix(,nrow=(length(parenttaxa$`180403`$taxonname)), ncol=4)
colnames(whale)<-c("SpeciesName", "AccNum", "SeqName", "SeqLen")
mcph_sequence<-character()

for(i in 1: length(unique(parenttaxa$`180403`$taxonname))){
  print(i)
  whale[i,1]<-parenttaxa$`180403`$taxonname[i]
  seqout<-entrez_search(db="nuccore", term=paste(parenttaxa$`180403`$taxonname[i], "[ORGN] AND 1:4500[SLEN] AND 
                                                 MCPH1[ALL] AND biomol_genomic[PROP]",
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
write.csv(whale, "mcph_gene_seq_info.csv")

