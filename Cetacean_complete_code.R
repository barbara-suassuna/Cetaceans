library(rentrez)
library(XML)
library(dplyr)
library(httr)
library(seqinr)

pull_taxa<-entrez_search(db="taxonomy", term="txid9721[SBTR]", retmax=2000)

#how long our table will be
length(unique(pull_taxa$ids))

#Create a Table
cetacea_data<-matrix(,nrow=length(unique(pull_taxa$ids)), ncol=5)
colnames(cetacea_data)<-c("taxID", "SpeciesName", "AccNum", "SeqName", "SeqLen")

#LOOP

#1 - Prestin (SLC26A5) gene
prestin_slc26a5_sequence<-character()
for(i in 1:length(unique(pull_taxa$ids))){
  print(i)
  cetacea_data[i,1]<-pull_taxa$ids[i]
  record<-entrez_summary(db="taxonomy", id=pull_taxa$ids[i])
  cetacea_data[i,2]<-record$scientificname
  seqout<-entrez_search(db="nuccore", term=paste("txid",pull_taxa$ids[i],"[ORGN] AND 1:4500[SLEN] 
                                                 AND PRESTIN[ALL]", sep=""), 
                        retmax=15)
  
  if(length(seqout$ids)<1){
    cetacea_data[i,3]<-NA
    cetacea_data[i,4]<-NA
    cetacea_data[i,5]<-NA
  }
  else{
    outsum<-entrez_summary(db="nuccore", id=seqout$ids)
    cetacea_data[i,3]<-outsum$accessionversion
    cetacea_data[i,4]<-outsum$title
    cetacea_data[i,5]<-outsum$slen
    
    seqout<-entrez_fetch(db="nuccore", id=outsum$accessionversion, rettype="fasta")
    
    seqout<-sub(">([^\n]*)", paste0(">",outsum$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    prestin_slc26a5_sequence<-c(prestin_slc26a5_sequence, seqout)
  }
  
}
#SAVING FASTA file
write(prestin_slc26a5_sequence, file="slc26a5_file.fasta")

#SAVE the TABLE
write.csv(cetacea_data, "cetacean_slc26a5_seq_info.csv")

#TEST the FILE
testing<-read.fasta("slc26a5_file.fasta")

#2 - CYTB gene
cytb_sequence<-character()
for(i in 1:length(unique(pull_taxa$ids))){
  print(i)
  cetacea_data[i,1]<-pull_taxa$ids[i]
  record<-entrez_summary(db="taxonomy", id=pull_taxa$ids[i])
  cetacea_data[i,2]<-record$scientificname
  seqout<-entrez_search(db="nuccore", term=paste("txid",pull_taxa$ids[i],"[ORGN] AND 1:4500[SLEN] 
                                                 AND CYTB[GENE]", sep=""), 
                        retmax=15)
  
  if(length(seqout$ids)<1){
    cetacea_data[i,3]<-NA
    cetacea_data[i,4]<-NA
    cetacea_data[i,5]<-NA
  }
  else{
    outsum<-entrez_summary(db="nuccore", id=seqout$ids)
    cetacea_data[i,3]<-outsum$accessionversion
    cetacea_data[i,4]<-outsum$title
    cetacea_data[i,5]<-outsum$slen
    
    seqout<-entrez_fetch(db="nuccore", id=outsum$accessionversion, rettype="fasta")
    
    seqout<-sub(">([^\n]*)", paste0(">",outsum$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    cytb_sequence<-c(cytb_sequence, seqout)
  }
}

#SAVING the FASTA file
write(cytb_sequence, file="cytb_file.fasta")

#SAVE the TABLE
write.csv(cetacea_data, "cetacean_cytb_seq_info.csv")

#TEST the file
testing<-read.fasta("cytb_file.fasta")

#3 - COI gene
coi_sequence<-character()
for(i in 1:length(unique(pull_taxa$ids))){
  print(i)
  cetacea_data[i,1]<-pull_taxa$ids[i]
  record<-entrez_summary(db="taxonomy", id=pull_taxa$ids[i])
  cetacea_data[i,2]<-record$scientificname
  seqout<-entrez_search(db="nuccore", term=paste("txid",pull_taxa$ids[i],"[ORGN] AND 1:4500[SLEN] 
                                                 AND COI[GENE]", sep=""), 
                        retmax=15)
  
  if(length(seqout$ids)<1){
    cetacea_data[i,3]<-NA
    cetacea_data[i,4]<-NA
    cetacea_data[i,5]<-NA
  }
  else{
    outsum<-entrez_summary(db="nuccore", id=seqout$ids)
    cetacea_data[i,3]<-outsum$accessionversion
    cetacea_data[i,4]<-outsum$title
    cetacea_data[i,5]<-outsum$slen
    
    seqout<-entrez_fetch(db="nuccore", id=outsum$accessionversion, rettype="fasta")
    
    seqout<-sub(">([^\n]*)", paste0(">",outsum$organism), seqout)
    seqout<-gsub(" ", "_", seqout)
    coi_sequence<-c(coi_sequence, seqout)
  }
}

#SAVING the FASTA file
write(coi_sequence, file="coi_file.fasta")

#SAVE the TABLE
write.csv(cetacea_data, "cetacean_coi_seq_info.csv")

#TEST the file
testing<-read.fasta("coi_file.fasta")
