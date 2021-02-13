# created by KSB, 20.01.19
# Downsample Indonesian and Malian populations to 12, 1000 times

# load packages
library(tidyverse)             
library(breakaway)        
library(DivNet)
library(taxize)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison_TB/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison_TB/"
# outputdir = "/data/cephfs/punim0586/kbobowik/R/"

###########################
# Alpha diversity testing #
###########################

# load in population data
load(paste0(inputdir, "Mali_Counts_physeq.Rda"))

# fill in taxa that CCMetagen couldn't identify with taxize for each of the datasets 
# First Mali
index_Mali=which(tax_table(Mali_Counts_physeq)[,"Phylum"]=="unk_p")
specieslist_Mali <- c(unname(tax_table(Mali_Counts_physeq)[index_Mali,"Family"]))
retrieved_taxa_Mali=tax_name(query = c(specieslist_Mali), get = c("phylum"), db = "ncbi")
retrieved_taxa_Mali[is.na(retrieved_taxa_Mali)]="unk_p"
tax_table(Mali_Counts_physeq)[index_Mali,"Phylum"]=retrieved_taxa_Mali[,"phylum"]

# set up matrix to fill in for Shannon doversity
downsampled_shannon_matrix = matrix(ncol=1, nrow = 1000)
colnames(downsampled_shannon_matrix) = c("Mali")

# set up matrix to fill in for Simpson doversity
downsampled_simpson_matrix = matrix(ncol=1, nrow = 1000)
colnames(downsampled_simpson_matrix) = c("Mali")

# set seed so that results are reproducible
set.seed(12435)

downsampled_diversity_test <- function(subsample, iteration) {

  # subsample data for Mali dataset
  Mali_subSample = sample(colnames(otu_table(Mali_Counts_physeq)), subsample, replace=F)
  SubSampled_Mali_Counts_physeq = prune_samples(Mali_subSample, Mali_Counts_physeq)
  # remove taxa with only 0's in the phyloseq object
  any(taxa_sums(SubSampled_Mali_Counts_physeq ) == 0)
  SubSampled_Mali_Counts_physeq=prune_taxa(taxa_sums(SubSampled_Mali_Counts_physeq) > 0, SubSampled_Mali_Counts_physeq)
  print(Mali_subSample)

  # Remove Viridiplantae and Metazoa
  phylo_counts_withSingletons <- subset_taxa(SubSampled_Mali_Counts_physeq, (Kingdom!="Viridiplantae"))
  phylo_counts_withSingletons <- subset_taxa(phylo_counts_withSingletons, (Kingdom!="Metazoa"))
  phylo_counts_withSingletons <- subset_taxa(phylo_counts_withSingletons, (Superkingdom!="unk_sk"))
  # remove any empty rows
  phylo_counts_withSingletons <- prune_taxa(taxa_sums(phylo_counts_withSingletons) > 0, phylo_counts_withSingletons)

  # compare diversity at the phylum level
  pop_comparison <- phylo_counts_withSingletons %>%
    tax_glom("Phylum")

  # If we don't change the sample names here from hyphens to periods, we'll get an error later
  sample_names(pop_comparison) <- gsub("\\-", ".", sample_names(pop_comparison))

  # Run divnet
  dv_pop_comparison <- divnet(pop_comparison, ncores = 4)

  summary_df_shannon <- as.data.frame(dv_pop_comparison$shannon %>% summary)

  #write.table(summary_df_shannon, file=paste0(outputdir, "Summary_diversity_Shannon_Indonesia_",iteration,".txt"), col.names=T, row.names=F, quote=F, sep="\t")

  # get shannon values
  Mali_estimate_shannon <<- mean(summary_df_shannon$estimate)
 
  # now simpson
  summary_df_simpson <- as.data.frame(dv_pop_comparison$simpson %>% summary)
  #write.table(summary_df_simpson, file=paste0(outputdir, "Summary_diversity_Simpson_Indonesia_",iteration,".txt"), col.names=T, row.names=F, quote=F, sep="\t")

  Mali_estimate_simpson <<- mean(summary_df_simpson$estimate)
}

for(i in 1:1000){
  downsampled_diversity_test(12, i)
  # shannon diversity
  downsampled_shannon_matrix[i,1] <- Mali_estimate_shannon
  # simpson
  downsampled_simpson_matrix[i,1] <- Mali_estimate_simpson
}

write.table(downsampled_shannon_matrix, file=paste0(outputdir, "Shannon_diversity_Mali_1000.txt"), col.names=T, row.names=F, quote=F, sep="\t")
write.table(downsampled_simpson_matrix, file=paste0(outputdir, "Simpson_diversity_Mali_1000.txt"), col.names=T, row.names=F, quote=F, sep="\t")
