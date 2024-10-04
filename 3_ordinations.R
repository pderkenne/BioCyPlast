# Aim: Plot ordinations to explore microbial communities 

# Performed in R V4.4.1

#librairies

library(phyloseq)
library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)


otu_table2 = as.data.frame(carbom@otu_table)
tax_table2 = as.data.frame(carbom@tax_table)

tax_otu_combined = merge(otu_table2, tax_table2, by = 0)

genus_samp = tax_otu_combined %>% select(-Rank1, -Rank2, -Rank3, -Rank4, -Rank5, -Rank7, -Row.names) 

# sum abundance by row for the same genus
genus_samp = setDT(genus_samp)[, lapply(.SD, sum, na.rm = TRUE), by = Rank6]

genus_samp = as.data.frame(genus_samp)
rownames(genus_samp) = genus_samp$Rank6
genus_samp= genus_samp %>% select (-Rank6)
genus_samp_t = as.data.frame(t(genus_samp))
genus_samp_t <- genus_samp_t %>% select(-c("unknown genus", "Multi-affiliation"))

write.csv(genus_samp_t, "df_genus.csv")
write.csv(genus_samp_t, "df_genus_unknown_added.csv")
