# Aim: Plot ordinations to explore microbial communities : exmample with a Principal Coordinate Analysis (PCoA)

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

# sites 

sites <- as.data.frame(spe.b.pcoa$points)

# species 
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], spe.h)
species = as.data.frame(spe.wa)

## Example of PCoA with 2 axis
# variance explained by the first two axis 

explained_var <- spe.b.pcoa$eig / sum(spe.b.pcoa$eig) * 100
pretty_pe <- format(round(explained_var, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

# prepare sites data frame for plot

colnames(sites) <- c("pcoa1", "pcoa2")
sites <- merge(sites, df.norm, by = "row.names", all = TRUE)

# if top 10 most contributive species to the first and second axis plotted on the PCoA  

# Top 10 most contributive taxa on the first axis 1

species$abs_pcoa1 <- abs(species$pcoa1)
sp_axis1 <- species[order(-species$abs_pcoa1), ]
sp_axis1_top10 <- sp_axis1[1:10, ]

# Top 10 most contributive taxa on the first axis 2

species$abs_pcoa2 <- abs(species$pcoa2)
sp_axis2 <- species[order(-species$abs_pcoa2), ]
sp_axis2_top10 <- sp_axis2[1:10, ]

# Top 20 most contributive taxa 

top_20 = merge(sp_axis1_top10, sp_axis2_top10, all= TRUE)

# plot on ggplot2

pcoa_plot <- ggplot(data = sites, aes(x = pcoa1, y = pcoa2)) +
  geom_point(aes(color = material,
                 size = time)
             na.rm = TRUE, stroke = 5) +
  geom_text_repel(aes(label = material), segment.color = "black", size = 3, max.overlaps = 50) +
  geom_segment(data = top_20, aes(xend = pcoa1, yend = pcoa2, x = 0, y = 0),
               arrow = arrow(type = "closed", length = unit(0.05, "inches")), color = "gray50") +
  geom_text_repel(data = top_20, aes(label = taxa), size = 5, nudge_x = 0.01, nudge_y = 0.01, max.overlaps = 50) +
  labs(x = labels[1], y = labels[2]) +
  scale_color_manual(values = c("powder" = "cadetblue3", "film" = "orange"), 
                     labels = c("PHBV/Cellulose 80/20 in powder", "PHBV/Cellulose 80/20 in film"),
                     guide = guide_legend(title = "Shape of the material")) +
  scale_size_manual(values = c("TO" = 1, 
                               "Tinter_Round n°1" = 5,  
                               "Tfinal_Round n°1" = 8, 
                               "Tinter_Round n°2" = 12, 
                               "Tfinal_Round n°2" = 16), 
                    labels = c("TO", 
                               expression(T[inter]*"_Round n°1"),
                               expression(T[final]*"_Round n°1"), 
                               expression(T[inter]*"_Round n°2"), 
                               expression(T[final]*"_Round n°2")), 
                    guide = guide_legend(title = "Sampling time")) +
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
  stat_ellipse(aes(linetype = factor(group)), level = 0.95, show.legend = TRUE) +  # Include group as linetype
  scale_linetype_manual(values = c("1" = "solid", "2" = "dashed"), 
                      labels = c("Group 1", "Group 2"), 
                        guide = guide_legend(title = "Group", override.aes = list(color = NA)))

print(pcoa_plot)
