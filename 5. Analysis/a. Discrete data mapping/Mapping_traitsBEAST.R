#install.packages("phytools")
library(phytools)
library(ggtree)
library(ape)
library(ggplot2)
library(geiger)
library(phangorn)

#-----------------------#
#Read our BEAST2 tree file and Time-calibrate our tree file
#------------------------#
bays_tree<- read.beast("Whales.MCC.tree")
balance<-max(node.depth.edgelength(bays_tree@phylo))
bays_tree@phylo$edge.length<-bays_tree@phylo$edge.length/balance
bays_tree@phylo$edge.length<-bays_tree@phylo$edge.length*55
max(node.depth.edgelength(bays_tree@phylo))

#Plot tree
plotTree(bays_tree@phylo, fsize=0.5)

#Read in the Discrete Data matrix csv file
discrete.data <- read.csv ("DiscreteDataMatrix.csv")


#Read Discrete Data
newdata<-data.frame(discrete.data$Dentition, discrete.data$Feeding) 
rownames(newdata)<-discrete.data$Species 
colnames(newdata)<-c("Dentition", "Feeding") 
newdata$Dentition<-as.factor(newdata$Dentition) 
newdata$Feeding<-as.factor(newdata$Feeding) 

#Plot Bayesian tree with a heatmap for discrete data
dis.tree<-plotTree.datamatrix(bays_tree@phylo, X=newdata, fsize=0.5)
dis.tree$colors
##define your colors based on the function above
cols <- list(
  Dentition = setNames(c("#7FC97F","#BEAED4","#FDC086","#FFFF99"), levels(newdata$Dentition)),
  Feeding   = setNames(c("#1B9E77","#D95F02","#7570B3","#E7298A"), levels(newdata$Feeding))
)

# Legend 1: Dentition
legend("topright",
       legend = names(cols$Dentition),
       fill   = cols$Dentition,
       title  = "Dentition",
       cex    = 0.4,
       inset  = c(0.02, 0.02))

# Legend 2: Feeding
legend("right",
       legend = names(cols$Feeding),
       fill   = cols$Feeding,
       title  = "Feeding",
       cex    = 0.4,
       inset  = c(0.007, 0.02))

##Export you plot as PDF
