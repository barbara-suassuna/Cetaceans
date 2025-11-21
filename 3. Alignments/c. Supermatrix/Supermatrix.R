#Load packages
library(ape)
library(evobiR)

#All aligned fasta files need to be in a single folder in your WD~
SuperMatrix(missing = "-", prefix = "supermatrix_data1", save = T)
data <- read.FASTA("supermatrix_data1.fasta")


