
# created by KSB, 20.01.19
# Downsample Indonesian and Malian populations to 12, 1000 times

# load packages
library(tidyverse)             
library(breakaway)        
library(DivNet)
library(taxize)

# Set paths:
inputdir = "/data/cephfs/punim0586/kbobowik/R/"
outputdir = "/data/cephfs/punim0586/kbobowik/R/"


###########################
# Alpha diversity testing #
###########################

# set seed so that results are reproducible
set.seed(12345)

downsampled_diversity_test <- function(subsample) {
  # load in population data
  load(paste0(inputdir, "Mali_Counts_physeq.Rda"))
  load(paste0(inputdir, "European_Counts_physeq.Rda"))
  load(paste0(inputdir, "AllREadsSE_Indo_Counts_physeq.Rda"))

  # subsample data for Indonesian dataset
  Indo_subSample = sample(colnames(otu_table(AllREadsSE_Indo_Counts_physeq)), subsample, replace=F)
  AllREadsSE_Indo_Counts_physeq = prune_samples(Indo_subSample,AllREadsSE_Indo_Counts_physeq)
  # remove taxa with only 0's in the phyloseq object
  any(taxa_sums(AllREadsSE_Indo_Counts_physeq) == 0)
  AllREadsSE_Indo_Counts_physeq=prune_taxa(taxa_sums(AllREadsSE_Indo_Counts_physeq) > 0, AllREadsSE_Indo_Counts_physeq)

  # subsample data for Malian dataset
  Malian_subSample = sample(colnames(otu_table(Mali_Counts_physeq)), subsample, replace=F)
  Mali_Counts_physeq = prune_samples(Malian_subSample, Mali_Counts_physeq)
  # remove taxa with only 0's in the phyloseq object
  any(taxa_sums(Mali_Counts_physeq) == 0)
  Mali_Counts_physeq=prune_taxa(taxa_sums(Mali_Counts_physeq) > 0, Mali_Counts_physeq)

  # Get phyloseq object with singletons by merging original phyloseq objects
  merged_phylo_counts_withSingletons <- merge_phyloseq(AllREadsSE_Indo_Counts_physeq, Mali_Counts_physeq, European_Counts_physeq)
  # remove any empty rows
  merged_phylo_counts_withSingletons <- prune_taxa(taxa_sums(merged_phylo_counts_withSingletons) > 0, merged_phylo_counts_withSingletons)

  # compare diversity at the phylum level
  pop_comparison <<- merged_phylo_counts_withSingletons %>%
    tax_glom("Phylum")

  # If we don't change the sample names here from hyphens to periods, we'll get an error later
  sample_names(pop_comparison) <- gsub("\\-", ".", sample_names(pop_comparison))

  # Run divnet
  dv_pop_comparison <- divnet(pop_comparison, ncores = 4)
}

# now that we have our function, repeat this 1,000 times
iterations = 1000
downsampledSummary = list()
# set the downsampled size to 12 (maximum number of individuals in the UK dataset)
downsampledSummary = replicate(iterations, downsampled_diversity_test(12), simplify=FALSE)

# create function to get mean Shannon diversity estimates for each population
getEstimate_Shannon = function(population){
  sapply(1:iterations, function(x) mean(filter(as.data.frame(downsampledSummary[[x]]$shannon %>% summary %>% 
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.))), SamplePop == population)$estimate))
}

# extract shannon estimates for each population
Indonesian_estimate = getEstimate_Shannon("Indonesia")
Malian_estimate = getEstimate_Shannon("Mali")
UK_estimate = getEstimate_Shannon("UK")

# save Shannon diversity table
Shannon_diversity_1000 = data.frame(Indonesian_estimate, Malian_estimate, UK_estimate)
write.table(Shannon_diversity_1000, file=paste0(outputdir, "Shannon_diversity_1000.txt"), col.names=T, row.names=F, quote=F, sep="\t")

