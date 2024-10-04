# Aim: Plot ordinations to explore microbial communities 

# Performed in R V4.4.1

#librairies

library(phyloseq)
library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

## creation of a raw abundance dataframe gathering samples as rows and genera taxa names as columns 
# Used phyloseq_rar object created from 1_phyloseq.R

otu_table2 = as.data.frame(phyloseq_rar@otu_table)
tax_table2 = as.data.frame(phyloseq_rar@tax_table)

# merge otu_table and tax_table from a phyloseq object rarefied
tax_otu_combined = merge(otu_table2, tax_table2, by = 0)

# normalization
genus_samp = tax_otu_combined %>% select(-Kingdom, -Phylum, -Class, -Order, -Family, -Species, -Row.names) 

# sum abundance by row for the same genus
genus_samp = setDT(genus_samp)[, lapply(.SD, sum, na.rm = TRUE), by = Genus]

# formatting of the final data frame
genus_samp = as.data.frame(genus_samp)
rownames(genus_samp) = genus_samp$Genus
genus_samp= genus_samp %>% select (-Genus)
genus_samp_t = as.data.frame(t(genus_samp)) # transposition to obtain samples as rows and genera as columns 

# write file
write.csv(genus_samp_t, "df_genus.csv")

## ordination

# load files
df <- read.csv("df_genus.csv", sep=";")
metadata <- read.csv("metadata.csv", sep=";")

# add useful metadata
df$material = metadata$material
df$time = metadata$time

# normalization 
df.norm <- subset(df, select = -c(X, material, time)) # material and time are examples of added metadata

# trasnformation and creation of dissimilarity matrix
spe.h <- decostand(df.norm, method = "hellinger")
spe.bray <- vegdist(spe.h)

# creation of the PCoA on sites
spe.b.pcoa <- cmdscale(spe.bray, k = (nrow(spe.h) - 1), eig = TRUE)

# species coordinates 
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], spe.h)

