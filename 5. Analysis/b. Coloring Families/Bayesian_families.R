library(ggtree) 
library(treeio) 
library(ggsci)
library(ggplot2)
library(tidytree)
library(phytools)

#-----Time-calibrate MCC (Maximum Clade Credibility) Tree-----#
##Upload your desired tree in read.tree ("")
bays_tree<- read.beast("Whales.MCC.tree")
balance<-max(node.depth.edgelength(bays_tree@phylo))
bays_tree@phylo$edge.length<-bays_tree@phylo$edge.length/balance
bays_tree@phylo$edge.length<-bays_tree@phylo$edge.length*55
max(node.depth.edgelength(bays_tree@phylo))

##Plot phylogeny with tip labels and time scale at the bottom
p <- ggtree(bays_tree, layout = "rectangular") +
  geom_tiplab(align = TRUE, fontface = "bold.italic", size = 2)+
  theme_tree2()+
  labs(caption= "Million of years")+ #adds scale bar at the bottom
  xlim(-60, 15)
p <- revts(p) #reverts time scale
p

#Add posterior values to the nodes 
p1 <- p +
  geom_nodelab(aes(label = round(posterior, 3)), 
               hjust = 1.5, vjust = -0.5, size = 1.8) + # adds posterior value to internal nodes
  xlim_tree(10) # adjust x axis to show full tip names
  
p1

# --- Define MRCA nodes --- #
mrca_node <- MRCA(p, c("Lagenorhynchus_albirostris", "Cephalorhynchus_commersonii"))
mrca_node2 <- MRCA(p, c("Phocoena_sinus","Neophocaena_phocaenoides"))
mrca_node3 <- MRCA(p, c("Monodon_monoceros", "Delphinapterus_leucas"))
mrca_node4 <- MRCA(p, c("Ziphius_cavirostris", "Berardius_bairdii"))
mrca_node5 <- MRCA(p, c("Megaptera_novaeangliae","Balaenoptera_acutorostrata"))
mrca_node6 <- MRCA(p, c("Eubalaena_glacialis", "Balaena_mysticetus"))
mrca_node7  <- MRCA(p, c("Kogia_sima", "Kogia_breviceps"))

# Build dataframe linking nodes to families
df <- data.frame(
  node = c(mrca_node, mrca_node2, mrca_node3, mrca_node4, mrca_node5, mrca_node6, mrca_node7),
  family = c("Delphinidae", "Phocoenidae", "Monodontidae",
             "Ziphiidae", "Balaenopteridae+Eschrichtiidae", 
             "Balaenidae","Kogiidae")
)

# --- Plot tree with highlights --- #
p2 <- p1 +
  geom_hilight(data = df, aes(node = node, fill = family), alpha = 0.4, extendto = 40) +
  scale_fill_manual(values = c(
    "Delphinidae" = "#E69F00",   # Orange
    "Phocoenidae" = "turquoise",  # Turquoise
    "Monodontidae" = "#009E73", # Bluish Green
    "Ziphiidae" =  "#CC79A7",   # Reddish Purple
    "Balaenopteridae+Eschrichtiidae" = "#0072B2", # Blue
    "Balaenidae" = "hotpink",   # Hotpink
    "Kogiidae" = "#F0E442"   # Yellow   
  )) +
  labs(fill = "Family")+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),legend.position = c(0, 1), # Place legend in top-left corner
  legend.justification = c(0, 1))


# Display 
p2

# Save plot as pdf

