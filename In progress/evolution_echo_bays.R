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

discrete.data <- read.csv ("DiscreteDataMatrix.csv",row.names=1,stringsAsFactors=TRUE)

echo <- setNames(discrete.data$Echolocation, rownames(discrete.data))
levels(echo)

#--Find best fit model
#Fit ER model
fit.ER <- fitMk (bays_tree@phylo, echo, model="ER")
print(fit.ER)

#Fit ARD model
fit.ARD <- fitMk (bays_tree@phylo, echo, model="ARD")
print(fit.ARD)

#COMPARE models and choose best one
AIC(fit.ER,fit.ARD)
#####"ER" is the best model with a lower AIC value

#Make simMap (and simulate a 1000 different trees)
smaps.echo<-make.simmap(bays_tree@phylo,echo,
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

head(n.changes,30) #looks at first 30 trees' summary of changes


##Set colors
cols<-setNames(c("black","#E7298A"),c("echo","no_echo"))

##Set map for first most common
i<-11
plot(smaps.echo[[i]],cols,ftype="i",fsize=0.6,
     outline=TRUE,lwd=0.75,ylim=c(-2,Ntip(bays_tree@phylo$tip.label)))
markChanges(smaps.echo[[i]],
            setNames(c("black","black"),names(cols)),cex=0.5,lwd=1)
markChanges(smaps.echo[[i]],cols,cex=0.5,lwd=2)
text(paste("Number of changes:",sum(n.changes[i,])),
     x=mean(par()$usr[1:2]),
     y=-4,cex=1.5)
legend("left",
       legend = names(cols),
       fill   = cols,
       title  = "Echolocation",
       cex    = 1,
       inset  = c(0.075, 0.075))
?plot
##Set map for second most common
i<-2
plot(smaps.echo[[i]],cols,ftype="i",fsize=0.6,
     outline=TRUE,lwd=0.75,ylim=c(-5,Ntip(root.tree)))
markChanges(smaps.echo[[i]],
            setNames(c("black","black"),names(cols)),cex=0.5,lwd=1)
markChanges(smaps.echo[[i]],cols,cex=0.5,lwd=2)
text(paste("Number of changes:",sum(n.changes[i,])),
     x=mean(par()$usr[1:2]),
     y=-4,cex=1.5)
legend("left",
       legend = names(cols),
       fill   = cols,
       title  = "Echolocation",
       cex    = 1,
       inset  = c(0.075, 0.075))
