# Aim: Create barplots with main phyla and the associated genera

# Performed in R V4.4.1

library(phyloseq)
library(fantaxtic)
library (ggplot2)

# Use phyloseq_rar object created from 1_phyloseq.R

## Normalization 
# remove unknown and multi-affiliation taxa level: here at genus level
taxa_to_remove <- tax_table(phyloseq_rar)[, "Genus"] %in% c("Multi-affiliation", "unknown genus")
phyloseq_rar <- subset_taxa(phyloseq_rar, !taxa_to_remove)

# new sample order 
sample_order <- c("PHBV_1", "PHBV2", "PHBV3", ...)

## Nested top 10 most abundant phyla with the top 2 most abundant genera and plot
top_nested <- nested_top_taxa(phyloseq_rar, top_tax_level = "Phylum", nested_tax_level = "Genus", n_top_taxa = 10, n_nested_taxa = 2)
p <- plot_nested_bar(ps_obj = top_nested$ps_obj, top_level = "Rank2", nested_level = "Rank6", sample_order = sample_order)
ggsave("plot_barplots.tiff", p, width = 10, height = 6, units = "in", dpi = 300)

