#install.packages("ggnewscale")
library(phytools)
library(ggtree)
library(ape)
library(ggplot2)
library(geiger)
library(phangorn)

#---------------------------------
#Read the tree file for ML tree
#---------------------------------
ML.tree <- read.tree ("supermatrix_data1.fasta.treefile")
ML.tree$tip.label
#Re-root the tree
root.tree<- root(ML.tree, outgroup = "Hippopotamus_amphibius", resolve.root = TRUE)

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
dis.tree<-plotTree.datamatrix(root.tree, X=newdata, fsize=0.5)
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

#-----------------------------------#
#Stochastic mapping discrete data
#-----------------------------------#
discrete.data <- read.csv ("DiscreteDataMatrix.csv",row.names=1,stringsAsFactors=TRUE)

echo <- setNames(discrete.data$Echolocation, rownames(discrete.data))
levels(echo)

#--Find best fit model
#Fit ER model
fit.ER <- fitMk (root.tree, echo, model="ER")
print(fit.ER)

#Fir ARD model
fit.ARD <- fitMk (root.tree, echo, model="ARD")
print(fit.ARD)

#COMPARE models and choose best one
AIC(fit.ER,fit.ARD)
#####"ER" is the best model with a lower AIC value

#Make simMap (and simulate a 1000 different trees)
smaps.echo<-make.simmap(root.tree,echo,
                             model="ER",nsim=1000,pi="fitzjohn")

#Posterior probability distribution on number of changes
dd<-density(smaps.echo)
dd
plot(dd,colors=c("black","#E7298A"),alpha=0.5)

#Find the number of changes of traits for each sampled tree
ss<-lapply(smaps.echo,summary)
ss
head(ss,5) #this looks only at the first 5 trees

#Get a summary of what changes seem to happen more often
n.changes<-t(sapply(ss,function(x) c(x$Tr[2,1],x$Tr[1,2])))
dimnames(n.changes)<-list(1:length(smaps.echo),
                          c("echo->no_echo","no_echo->echo"))

head(n.changes,60) #looks at first 60 trees' summary of changes

##Set colors
cols<-setNames(c("black","#E7298A"),c("echo","no_echo"))

##Set map
i<-1
plot(smaps.echo[[i]],cols,ftype="i",fsize=0.6,
     outline=TRUE,lwd=0.75,ylim=c(-5,Ntip(root.tree)))
markChanges(smaps.echo[[i]],
            setNames(c("black","black"),names(cols)),cex=0.5,lwd=1)
markChanges(smaps.echo[[i]],cols,cex=0.5,lwd=2)
text(paste("Number of changes:",sum(n.changes[i,])),
     x=mean(par()$usr[1:2]),
     y=-4,cex=1.5)
legend("right",
              legend = names(cols),
              fill   = cols,
              title  = "Echolocation",
              cex    = 1,
              inset  = c(0.075, 0.075))
