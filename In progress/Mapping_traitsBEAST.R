#install.packages("ggnewscale")
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

# #Add scale bar
# add.scale.bar(bays_tree@phylo, length = 55, x = 2,          
#               y = 5)
# add.scale.bar(bays_tree@phylo, x = 0, y = 5, length = 55)

#Read in the Discrete Data matrix
discrete.data <- read.csv ("DiscreteDataMatrix.csv")


#Read Discrete Data
newdata<-data.frame( discrete.data$Echolocation, discrete.data$Dentition, discrete.data$Feeding) 
rownames(newdata)<-discrete.data$Species 
colnames(newdata)<-c( "Echolocation","Dentition", "Feeding") 
newdata$Echolocation <- as.factor(newdata$Echolocation)
newdata$Dentition<-as.factor(newdata$Dentition) 
newdata$Feeding<-as.factor(newdata$Feeding) 

#Plot ML tree with a heatmap for discrete data
dis.tree<-plotTree.datamatrix(bays_tree@phylo, X=newdata, fsize=0.5)
dis.tree$colors
##define your colors based on the function above
cols <- list(
  Echolocation      = setNames(c("#7FC97F","#FDC086"), levels(newdata$Echolocation)),
  Dentition = setNames(c("#1B9E77","#D95F02","#7570B3","#E7298A" ), levels(newdata$Dentition)),
  Feeding   = setNames(c( "#A6CEE3","#1F78B4","#B2DF8A","#33A02C"), levels(newdata$Feeding))
)

# Legend 1: Echolocation
legend("topright",
       legend = names(cols$Echolocation),
       fill   = cols$Echolocation,
       title  = "Echolocation",
       cex    = 0.5,
       inset  = c(0.01, 0.05))

# Legend 2: Dentition
legend("right",
       legend = names(cols$Dentition),
       fill   = cols$Dentition,
       title  = "Dentition",
       cex    = 0.4,
       inset  = c(0.02, 0.02))

# Legend 3: Feeding
legend("bottomright",
       legend = names(cols$Feeding),
       fill   = cols$Feeding,
       title  = "Feeding",
       cex    = 0.4,
       inset  = c(0.007, 0.02))

##Export you plot as PDF
