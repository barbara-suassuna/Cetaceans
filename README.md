Cetaceans Phylogenetics

This repository documents the workflow used to reconstruct a Cetacean phylogeny based on six target genes. It is organized into numbered folders outlining the major steps of the project.

Summary of Workflow

1. Data collection: An R script is provided to mine sequence data from NCBI. This script retrieves the target genes and outputs both FASTA files and CSV tables containing the extracted DNA sequence information.
2. Filtered taxa: A filtering script is provided to reduce the initial list of taxa (~104) to include only those with sequence data for a minimum of three genes.
3. Alignments: This folder contains the raw alignments generated with MAFFT, the cleaned alignments produced using an R script, and a script to create a concatenated supermatrix with all the alignments. These alignments serve as the input for phylogenetic inference.
4. Phylogenetic Inference: Both IQ-TREE and BEAST2 were used to reconstruct a phylogeny. Each has its own subfolder with the outputs.
     - IQ-TREE: the concatenated supermatrix was analyzed to create a maximum-likelihood tree.
     - BEAST2: a partitioned NEXUS file was created with the python script provided and imported into BEAUti to set up the XML file. The latter was input into BEAST2. 
6. Analysis: This folder includes scripts to visualize the trees with support values, create a heatmap of discrete traits, and perform stochastic character mapping of echolocation. There is a subfolder for each analysis.
7. Bibliography: References used throughout the project.
8. Final Tree: Final phylogeny from the analyses above.
