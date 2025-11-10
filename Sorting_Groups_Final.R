library(seqinr)

# Extract species names from CSV
extract_species_names <- function(csv_file_path) {
  species_list <- c()
  lines <- readLines(csv_file_path)
  lines <- lines[-1]  # Skip header
  
  for (line in lines) {
    columns <- strsplit(line, ",")[[1]]
    if (length(columns) > 1) {
      species_name <- gsub('"', '', trimws(columns[2]))
      species_list <- c(species_list, species_name)
    }
  }
  return(species_list)
}

# Initialize species_dict with empty named lists
species_names <- extract_species_names("cox1-info.csv") # We can use any of our csv files for this step. 
species_dict <- setNames(vector("list", length(species_names)), species_names)

# Function to clean FASTA names
clean_fasta_names <- function(fasta) {
  gsub("_", " ", names(fasta))
}

# Function to process a FASTA file and update species_dict
process_fasta <- function(fasta_path, species_dict) {
  fasta <- read.fasta(fasta_path)
  fasta_species <- clean_fasta_names(fasta)
  
  for (species in names(species_dict)) {
    match <- species %in% fasta_species
    species_dict[[species]][[fasta_path]] <- match
  }
  
  return(species_dict)
}

# Process multiple FASTA files
fasta_files <- c("cox1_gene.fasta", "csn_gene.fasta", "enam_gene.fasta", 
                 "mc_gene.fasta", "mcph_gene.fasta", "prm1_gene.fasta")
for (file in fasta_files) {
  species_dict <- process_fasta(file, species_dict)
}

# for (species in names(species_dict)) {
#   cat("Species:", species, "\n")
#   print(species_dict[[species]])
# }

# Set your gene threshold
min_genes_required <- 3

# Prepare summary vectors
species_summary <- data.frame(
  species = character(),
  total_matches = integer(),
  threshold_met = logical(),
  stringsAsFactors = FALSE
)

# Build the summary
for (species in names(species_dict)) {
  match_vector <- unlist(species_dict[[species]])
  true_count <- sum(match_vector)
  threshold_status <- true_count >= min_genes_required
  
  species_summary <- rbind(
    species_summary,
    data.frame(
      species = species,
      total_matches = true_count,
      threshold_met = threshold_status,
      stringsAsFactors = FALSE
    )
  )
}

# View the summary
print(species_summary)

#Write the table to a csv file
write.csv(species_summary, "species_gene_threshold_summary.csv", 
          row.names = FALSE)

# Identify species that meet the threshold
passing_species <- species_summary$species[species_summary$threshold_met]

# Function to filter and write FASTA
filter_fasta <- function(fasta_path, output_path, passing_species) {
  fasta <- read.fasta(fasta_path)
  fasta_names <- clean_fasta_names(fasta)
  
  # Keep only sequences whose cleaned name is in passing_species
  keep_indices <- which(fasta_names %in% passing_species)
  filtered_fasta <- fasta[keep_indices]
  
  write.fasta(sequences = filtered_fasta,
              names = names(filtered_fasta),
              file.out = output_path)
}

# Apply filtering to each FASTA file
for (file in fasta_files) {
  output_file <- paste0("filtered_", file)
  filter_fasta(file, output_file, passing_species)
}

species_dict$`Hippopotamus amphibius`
