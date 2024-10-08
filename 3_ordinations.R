## Aim: Principal Coordinate Analysis (PCoA) plot to analyze the microbial community composition associated with PHBV/Cellulose composites during anaerobic digestion.
# The plot visualizes how microbial communities shift over time and across different cellulose percentages. 
# Samples are represented as points, with color and shape indicating cellulose content and sample type. 
# Arrows illustrate the influence of specific microbial taxa, helping to identify patterns in community structure relative to the degradation process.

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
df$cellulose = metadata$cellulose
df$time = metadata$time
df$sample = metadata$sample

# Subset the data to remove unnecessary columns (metadata)
df.norm <- subset(df, select = -c(X, material, time, sample)) # Exclude metadata columns

# Transform the data and create a dissimilarity matrix using Hellinger transformation
spe.h <- decostand(df.norm, method = "hellinger")
spe.bray <- vegdist(spe.h)

# Perform Principal Coordinate Analysis (PCoA) on the dissimilarity matrix
spe.b.pcoa <- cmdscale(spe.bray, k = (nrow(spe.h) - 1), eig = TRUE)

# Extract site coordinates from the PCoA result
sites <- as.data.frame(spe.b.pcoa$points)

# Perform weighted averaging to associate species with PCoA axes
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], spe.h)
species <- as.data.frame(spe.wa)

# Calculate the variance explained by the first two PCoA axes
explained_var <- spe.b.pcoa$eig / sum(spe.b.pcoa$eig) * 100
pretty_pe <- format(round(explained_var, digits = 1), nsmall = 1, trim = TRUE)
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
             glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

# Prepare the sites data frame for plotting
colnames(sites) <- c("pcoa1", "pcoa2")
sites <- merge(sites, df.norm, by = "row.names", all = TRUE)

# Identify the top 10 most contributive species for each axis
species$abs_pcoa1 <- abs(species$pcoa1)
sp_axis1_top10 <- species[order(-species$abs_pcoa1), ][1:10, ]

species$abs_pcoa2 <- abs(species$pcoa2)
sp_axis2_top10 <- species[order(-species$abs_pcoa2), ][1:10, ]

# Combine the top 10 species from both axes
top_20 <- merge(sp_axis1_top10, sp_axis2_top10, all = TRUE)

# Create the PCoA plot using ggplot2
pcoa_plot <- ggplot(data = sites, aes(x = pcoa1, y = pcoa2)) +
  geom_point(aes(shape = sample,
                 fill = as.factor(cellulose), 
                 color = as.factor(cellulose), 
                 size = time)) +
  geom_text_repel(aes(label = interaction_samples), segment.color = "black", size = 3) +
  geom_segment(data = top_scores_taxa_plot_axis_1_2, 
               aes(xend = pcoa1, yend = pcoa2, x = 0, y = 0),
               arrow = arrow(type = "closed", length = unit(0.05, "inches")), color = "gray50") +
  geom_text(data = top_scores_taxa_plot_axis_1_2, 
            aes(x = pcoa1, y = pcoa2, label = taxa), 
            size = 3, hjust = -0.1, vjust = 0.1) +
  labs(x = labels[1], y = labels[2]) +
  scale_size_continuous(name = "Time", range = c(2, 15), guide = guide_legend(title = "Time")) +
  scale_color_manual(values = c("0" = "red", "20" = "green", "40" = "blue", "100" = "purple"), 
                     guide = guide_legend(title = "% Cellulose")) +
  scale_shape_manual(values = c("Blank" = 22, "Samples" = 23, "Cellulose" = 24), 
                     guide = guide_legend(title = "Sample Type")) +
  labs(fill = "% Cellulose", shape = "Sample Type") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

# Display the PCoA plot
print(pcoa_plot)
