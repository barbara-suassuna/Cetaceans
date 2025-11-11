#install.packages("evobiR")
library(evobiR)
library(ape)
library(msaR)
library(phangorn)

#1 - COX1 GENE
cox <- read.FASTA("raw_cox1.fasta")
length(cox)
msaR(cox) # This shows the alignment in R

#trimming
cox2<- phyDat(cox, type = "DNA")
#this is where we trim out data, 60%
cox3 <- cox2[, colMeans(as.character(cox2)== "-") < 0.6]
cox4 <- as.DNAbin(cox3)
msaR(cox4)
# This writes a new file with the trimmed alignment
write.FASTA(cox4, "cox1_trimmed.fasta")

#2 - MC1R GENE
mc <- read.FASTA("raw_mc.fasta")
length(mc)

msaR(mc) # This shows the alignment in R

#trimming
mc2<- phyDat(mc, type = "DNA")
#this is where we trim out data, 60%
mc3 <- mc2[, colMeans(as.character(mc2)== "-") < 0.6]
mc4 <- as.DNAbin(mc3)
msaR(mc4)
# This writes a new file with the trimmed alignment
write.FASTA(mc4, "mc_trimmed.fasta")

#3- ENAM GENE
enam <- read.FASTA("raw_enam.fasta")
length(enam)
msaR(enam) # This shows the alignment in R

#trimming
enam2<- phyDat(enam, type = "DNA")
#this is where we trim out data, 60%
enam3 <- enam2[, colMeans(as.character(enam2)== "-") < 0.6]
enam4 <- as.DNAbin(enam3)
msaR(enam4)
# This writes a new file with the trimmed alignment
write.FASTA(enam4, "enam_trimmed.fasta")

#4- CSN2 GENE
csn <- read.FASTA("raw_csn.fasta")
length(csn)
msaR(csn) # This shows the alignment in R

#trimming
csn2<- phyDat(csn, type = "DNA")
#this is where we trim out data, 60%
csn3 <- csn2[, colMeans(as.character(csn2)== "-") < 0.6]
csn4 <- as.DNAbin(csn3)
msaR(csn4)
# This writes a new file with the trimmed alignment
write.FASTA(csn4, "csn_trimmed.fasta")

#5 - PRM1 GENE
prm1 <- read.FASTA("raw_prm1.fasta")
length(prm1)

msaR(prm1) # This shows the alignment in R

#trimming
prm1.2<- phyDat(prm1, type = "DNA")
#this is where we trim out data, 60%
prm1.3 <- prm1.2[, colMeans(as.character(prm1.2)== "-") < 0.6]
prm1.4 <- as.DNAbin(prm1.3)
msaR(prm1.4)
# This writes a new file with the trimmed alignment
write.FASTA(prm1.4, "prm1_trimmed.fasta")

#6- MCPH1 GENE
mcph <- read.FASTA("raw_mcph.fasta")
length(mcph)

msaR(mcph) #This shows the alignment in R

#trimming
mcph2<- phyDat(mcph, type = "DNA")
#this is where we trim out data, 40%
mcph3 <- mcph2[, colMeans(as.character(mcph2)== "-") < 0.4]
mcph4 <- as.DNAbin(mcph3)
msaR(mcph4)
#This writes a new file with the trimmed alignment
write.FASTA(mcph4, "mcph1_trimmed.fasta")
