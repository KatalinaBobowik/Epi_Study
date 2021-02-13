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

###########################
# Alpha diversity testing #
###########################

# load in Indonesian population data
load(paste0(inputdir, "AllREadsSE_Indo_Counts_physeq.Rda"))

# Indonesia
index_Indo=which(tax_table(AllREadsSE_Indo_Counts_physeq)[,"Phylum"]=="unk_p")
specieslist_Indo <- c(unname(tax_table(AllREadsSE_Indo_Counts_physeq)[index_Indo,"Family"]))
retrieved_taxa_Indo=tax_name(query = c(specieslist_Indo), get = c("phylum"), db = "ncbi")
retrieved_taxa_Indo[is.na(retrieved_taxa_Indo)]="unk_p"
tax_table(AllREadsSE_Indo_Counts_physeq)[index_Indo,"Phylum"]=retrieved_taxa_Indo[,"phylum"]

# load in Malian population data
load(paste0(inputdir, "Mali_Counts_physeq.Rda"))

# fill in taxa that CCMetagen couldn't identify with taxize for each of the datasets 
index_Mali=which(tax_table(Mali_Counts_physeq)[,"Phylum"]=="unk_p")
specieslist_Mali <- c(unname(tax_table(Mali_Counts_physeq)[index_Mali,"Family"]))
retrieved_taxa_Mali=tax_name(query = c(specieslist_Mali), get = c("phylum"), db = "ncbi")
retrieved_taxa_Mali[is.na(retrieved_taxa_Mali)]="unk_p"
tax_table(Mali_Counts_physeq)[index_Mali,"Phylum"]=retrieved_taxa_Mali[,"phylum"]

# UK
load(paste0(inputdir, "European_Counts_physeq.Rda"))
# fill in taxa that CCMetagen couldn't identify with taxize for each of the datasets 
index_UK=which(tax_table(European_Counts_physeq)[,"Phylum"]=="unk_p")
specieslist_UK <- c(unname(tax_table(European_Counts_physeq)[index_UK,"Family"]))
retrieved_taxa_UK=tax_name(query = c(specieslist_UK), get = c("phylum"), db = "ncbi")
retrieved_taxa_UK[is.na(retrieved_taxa_UK)]="unk_p"
tax_table(European_Counts_physeq)[index_UK,"Phylum"]=retrieved_taxa_UK[,"phylum"]

# set seed so that results are reproducible
set.seed(12345)


downsampled_diversity_test <- function(subsample, iteration) {

downsampled_matrix = matrix(nrow=1000,ncol=3)
colnames(downsampled_matrix) = c("Indo_est", "Mali_est", "UK_est")

for (i in 1:1000){
  # subsample data for Indonesian dataset
  Indo_subSample = sample(colnames(otu_table(AllREadsSE_Indo_Counts_physeq)), subsample, replace=F)
  SubSampled_Indo_Counts_physeq = prune_samples(Indo_subSample,AllREadsSE_Indo_Counts_physeq)
  # remove taxa with only 0's in the phyloseq object
  any(taxa_sums(SubSampled_Indo_Counts_physeq) == 0)
  SubSampled_Indo_Counts_physeq=prune_taxa(taxa_sums(SubSampled_Indo_Counts_physeq) > 0, SubSampled_Indo_Counts_physeq)
  print(Indo_subSample)

  # subsample data for Mali dataset
  Mali_subSample = sample(colnames(otu_table(Mali_Counts_physeq)), subsample, replace=F)
  SubSampled_Mali_Counts_physeq = prune_samples(Mali_subSample, Mali_Counts_physeq)
  # remove taxa with only 0's in the phyloseq object
  any(taxa_sums(SubSampled_Mali_Counts_physeq ) == 0)
  SubSampled_Mali_Counts_physeq=prune_taxa(taxa_sums(SubSampled_Mali_Counts_physeq) > 0, SubSampled_Mali_Counts_physeq)
  print(Mali_subSample)

  # Remove Viridiplantae and Metazoa
  merged_phylo_counts_withSingletons <- merge_phyloseq(SubSampled_Indo_Counts_physeq, SubSampled_Mali_Counts_physeq, European_Counts_physeq)
  phylo_counts_withSingletons <- subset_taxa(merged_phylo_counts_withSingletons, (Kingdom!="Viridiplantae"))
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
  dv_pop_comparison <- pop_comparison %>% divnet(X = "SamplePop", ncores = 8)

  summary_df_shannon <- as.data.frame(dv_pop_comparison$shannon %>%
  summary %>%
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

  Indo_est = mean(summary_df_shannon[which(summary_df_shannon$SamplePop == "Indonesia"),"estimate"])
  Mali_est = mean(summary_df_shannon[which(summary_df_shannon$SamplePop == "Mali"),"estimate"])
  UK_est = mean(summary_df_shannon[which(summary_df_shannon$SamplePop == "UK"),"estimate"])

  downsampled_matrix[i, "Indo_est"] = Indo_est
  downsampled_matrix[i, "Mali_est"] = Mali_est
  downsampled_matrix[i, "UK_est"] = UK_est

}

write.table(downsampled_matrix, file=paste0(outputdir, "Shannon_diversity_1000_covariatesInModel.txt"), col.names=T, row.names=F, quote=F, sep="\t")

# now that we have our function, repeat this 1,000 times
iterations = 1000
downsampledSummary = list()
# set the downsampled size to 12 (maximum number of individuals in the UK dataset)
downsampledSummary = replicate(iterations, downsampled_diversity_test(12), simplify=FALSE)

# create function to get Shannon diversity estimates for each population
getEstimate_Shannon = function(population){
  sapply(1:iterations, function(x) mean(filter(as.data.frame(downsampledSummary[[x]]$shannon %>% summary %>% 
  add_column("SamplePop" = c(rep("Indonesia",12),rep("Mali",12),rep("UK",12)))), SamplePop == population)$estimate))
}

# extract shannon estimates for each population
Indonesian_estimate = getEstimate_Shannon("Indonesia")
Malian_estimate = getEstimate_Shannon("Mali")
UK_estimate = getEstimate_Shannon("UK")

# save Shannon diversity table
Shannon_diversity_1000 = data.frame(Indonesian_estimate, Malian_estimate, UK_estimate)
write.table(Shannon_diversity_1000, file=paste0(outputdir, "Shannon_diversity_1000_covariatesInModel.txt"), col.names=T, row.names=F, quote=F, sep="\t")

# create function to get Shannon diversity estimates for each population
getEstimate_Simpson = function(population){
  sapply(1:iterations, function(x) mean(filter(as.data.frame(downsampledSummary[[x]]$simpson %>% summary %>% 
  add_column("SamplePop" = c(rep("Indonesia",12),rep("Mali",12),rep("UK",12)))), SamplePop == population)$estimate))
}

# extract shannon estimates for each population
Indonesian_estimate = getEstimate_Simpson("Indonesia")
Malian_estimate = getEstimate_Simpson("Mali")
UK_estimate = getEstimate_Simpson("UK")

# save Shannon diversity table
Simpson_diversity_1000 = data.frame(Indonesian_estimate, Malian_estimate, UK_estimate)
write.table(Simpson_diversity_1000, file=paste0(outputdir, "Simpson_diversity_1000_covariatesInModel.txt"), col.names=T, row.names=F, quote=F, sep="\t")

