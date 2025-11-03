install.packages("ape")
library(evobiR)
library(ape)
library(msaR)
library(phangorn)

#1 - COX1 GENE
cox <- read.FASTA("raw_cox1_gene_alignment.fasta")
length(cox)

length(unique(names(cox))) #this takes cares of any duplicate names
# Look at the Multiple Sequence Alignment 
msaR(cox) # This shows the alignment in R
# It is best to clean up the alignment in each gene rather than in the super matrix

#trimming
cox2<- phyDat(cox, type = "DNA")
#this is where we trim out data, 60%
cox3 <- cox2[, colMeans(as.character(cox2)== "-") < 0.6]
cox4 <- as.DNAbin(cox3)
msaR(cox4)
# This writes a new file with the trimmed alignment
write.FASTA(cox4, "cox1_trimmed.fasta")

#2 - CSN2 GENE
csn <- read.FASTA("raw_csn_gene_alignment.fasta")
length(csn)

length(unique(names(csn)))  
msaR(csn) 

#trimming
csn2<- phyDat(csn, type = "DNA")
#this is where we trim out data, 60%
csn3 <- csn2[, colMeans(as.character(csn2)== "-") < 0.6]
csn4 <- as.DNAbin(csn3)
msaR(cytb4)
# This writes a new file with the trimmed alignment
write.FASTA(csn4, "csn_trimmed.fasta")

#3 - PMR1 GENE
pmr <- read.FASTA("raw_pmr1_gene_alignment.fasta")
length(pmr)

length(unique(names(pmr)))  
msaR(pmr) 

#trimming
pmr2<- phyDat(pmr, type = "DNA")
#this is where we trim out data, 60%
pmr3 <- pmr2[, colMeans(as.character(pmr2)== "-") < 0.6]
pmr4 <- as.DNAbin(pmr3)
msaR(pmr4)
# This writes a new file with the trimmed alignment
write.FASTA(pmr4, "pmr1_trimmed.fasta")


#4 MC1R GENE
mc <- read.FASTA("raw_mc_gene_alignment.fasta")
length(mc)

length(unique(names(mc)))  
msaR(mc) 

#trimming
mc2<- phyDat(mc, type = "DNA")
#this is where we trim out data, 60%
mc3 <- mc2[, colMeans(as.character(mc2)== "-") < 0.6]
mc4 <- as.DNAbin(mc3)
msaR(mc4)
# This writes a new file with the trimmed alignment
write.FASTA(mc4, "mc_trimmed.fasta")

#5 MCPH1 GENE
mcph <- read.FASTA("raw_mcph_gene_alignment.fasta")
length(mcph)

length(unique(names(mcph))) 
msaR(mcph)
#trimming
mcph2<- phyDat(mcph, type = "DNA")
#this is where we trim out data, 60%
mcph3 <- mcph2[, colMeans(as.character(mcph2)== "-") < 0.6]
mcph4 <- as.DNAbin(mcph3)
msaR(mcph4)
write.FASTA(mc4, "mcph1_trimmed.fasta")
