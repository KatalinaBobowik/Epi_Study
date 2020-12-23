
all_taxa <- sapply(strsplit(rownames(tax_table(merged_phylo_counts_zComposition)), "[_.]"), `[`, 4)
all_taxa <- all_taxa %>% gsub("unk", NA, .) %>% gsub("\\bp\\b", NA, .)
all_taxa <- all_taxa[-which(is.na(all_taxa))]
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")

taxa_matrix <- matrix(ncol=5, nrow=length(all_taxa))
rownames(taxa_matrix) <- all_taxa
for (i in 1:5){
  for (taxa in all_taxa){
   logged_phyla_counts <- log10(colSums(otu_table(AllREadsSE_Indo_Counts_physeq)[grep(taxa,tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"])])+1)
   taxa_matrix[taxa,i] <- anova(lm(zComposition.ord$vectors[,i] ~ logged_phyla_counts))[1,5]
 }
}
