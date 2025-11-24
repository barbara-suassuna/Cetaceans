#load packages
library(phytools)
library(ggtree)
library(viridis)
library(ape)
library(ggplot2)
library(geiger)

#Read in the Discrete Data matrix
discrete.data <- read.csv ("Discrete Data Matrix.csv")

# Replace spaces in Species column to match FASTA files
discrete.data$Species <- gsub(" ", "_", discrete.data$Species)


##Trait 1
dentition <- discrete.data$Species
names(dentition) <- discrete.data$Dentition

##Trait 2
feeding_type <- discrete.data$Species
names(feeding_type) <- discrete.data$Feeding

##Trait 3
diet <- discrete.data$Species
names(diet) <- discrete.data$Diet

?geom_text
#Read the tree file for ML tree
read.ML <- "supermatrix_data1.fasta.treefile"
ML.tree<-read.tree (read.ML)
ML.tree$tip.label
#Reroot the tree
root.tree<- root(ML.tree, outgroup = "Hippopotamus_amphibius", resolve.root = TRUE)

tree1 <-ggtree(root.tree,cex=0.3, layout="rectangular") + 
  geom_tiplab(align=T,fontface="bold.italic", size=2) +
  xlim(0,0.20)
 
          
  



