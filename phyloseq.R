# Aim: Create a phyloseq object for DNA 16S sequencing data 
# Performed in R V4.4.1

BiocManager::install("phyloseq")
BiocManager::install("biomformat")

# libraries
library(phyloseq)
library(biomformat)

# BIOM data loading
biom_file <- "FROGS.biom1"
biom_data <- import_biom(biom_file)

# metadata loading 
sample_metadata <- read.csv("metadata.csv", sep = ";")
sample_data <- sample_data(sample_metadata)

# phyloseq object
phyloseq <- merge_phyloseq(biom_data, sample_data)

# rarefaction 
set.seed(123)
phyloseq_rar <- rarefy_even_depth(phyloseq, replace = FALSE)



