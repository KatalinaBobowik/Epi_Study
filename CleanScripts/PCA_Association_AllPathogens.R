
all_taxa <- sapply(strsplit(rownames(tax_table(merged_phylo_counts_zComposition)), "[_.]"), `[`, 4)
all_taxa <- all_taxa %>% gsub("unk", NA, .) %>% gsub("\\bp\\b", NA, .)
all_taxa <- all_taxa[-which(is.na(all_taxa))]
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")
taxa_names(AllREadsSE_Indo_Counts_physeq) <- make.unique(paste(tax_table(AllREadsSE_Indo_Counts_physeq)[,"Superkingdom"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Kingdom"], tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"], sep="_"))


# taxa_matrix <- matrix(ncol=5, nrow=length(all_taxa))
# rownames(taxa_matrix) <- all_taxa
# for (i in 1:5){
#   for (taxa in all_taxa){
#    logged_phyla_counts <- log10(colSums(otu_table(AllREadsSE_Indo_Counts_physeq)[grep(taxa,tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"])])+1)
#    taxa_matrix[taxa,i] <- anova(lm(zComposition.ord$vectors[,i] ~ logged_phyla_counts))[1,5]
#  }
# }

taxa_matrix <- matrix(ncol=5, nrow=length(all_taxa))
rownames(taxa_matrix) <- all_taxa
for (i in 1:5){
  for (taxa in all_taxa){
   logged_phyla_counts <- log10(colSums(otu_table(AllREadsSE_Indo_Counts_physeq)[grep(taxa, rownames(tax_table(AllREadsSE_Indo_Counts_physeq))),])+1)
   taxa_matrix[taxa,i] <- anova(lm(zComposition.ord$vectors[,i] ~ logged_phyla_counts))[1,5]
 }
}

write.table(taxa_matrix, file=paste0(outputdir,"ANOVA_taxa_matrix.txt"), sep="\t")


-----

# alternative

all_taxa <- rownames(tax_table(merged_phylo_counts_zComposition))
zComposition.ord <- ordinate(merged_phylo_counts_zComposition, method = "PCoA", distance = "euclidean")


taxa_matrix <- matrix(ncol=5, nrow=length(all_taxa))
rownames(taxa_matrix) <- all_taxa
for (i in 1:5){
  for (taxa in all_taxa){
   logged_phyla_counts <- log10(otu_table(pop_comparison)[taxa,]+1)
   taxa_matrix[taxa,i] <- anova(lm(zComposition.ord$vectors[,i] ~ matrix(logged_phyla_counts)))[1,5]
 }
}

---

# Perform ANOVA and get matrix of significance levels for all taxa
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








