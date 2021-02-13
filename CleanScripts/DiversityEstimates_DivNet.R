
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

# load in population data
load(paste0(inputdir, "Mali_Counts_physeq.Rda"))
load(paste0(inputdir, "European_Counts_physeq.Rda"))
load(paste0(inputdir, "AllREadsSE_Indo_Counts_physeq.Rda"))

# Get phyloseq object with singletons by merging original phyloseq objects
merged_phylo_counts_withSingletons <- merge_phyloseq(AllREadsSE_Indo_Counts_physeq, Mali_Counts_physeq, European_Counts_physeq)
# remove any empty rows
merged_phylo_counts_withSingletons <- prune_taxa(taxa_sums(merged_phylo_counts_withSingletons) > 0, merged_phylo_counts_withSingletons)

# Now that we have data with singletons, we now need to merge our data at a specified taxonomic level. DivNet is computationally expensive, and therefore a higher level is much, much faster.
# We'll therefore test how our groups look like at the Phylum level. Then, we'll run DivNet without specifying any hypothesis testing.

# comparing diversity at the phylum level
pop_comparison <- merged_phylo_counts_withSingletons %>%
  tax_glom("Phylum")

# If we don't change the sample names here from hyphens to periods, we'll get an error later
sample_names(pop_comparison) <- gsub("\\-", ".", sample_names(pop_comparison))

# Run divnet without specifying any hypothesis testing
dv_pop_comparison <- divnet(pop_comparison, ncores = 4)

# DivNet will output an object with estimates for multiple different alpha (and beta) diversity measures 
# Let's plot the results of shannon diversity
summary_df_shannon <- as.data.frame(dv_pop_comparison$shannon %>%
  summary %>%
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

# save shannon diversity summary
write.table(summary_df_shannon, file=paste0(outputdir,"summary_df_shannon_noHypothesisTesting.txt"))

# Take out the Shannon diversity metric from DivNet and plot it
ggplot(summary_df_shannon, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_boxplot(width=0.08, outlier.color = NA) + 
  scale_fill_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + ggtitle("Shannon Diversity") +
  ylab("Estimate of Shannon Diversity") + xlab("") + theme_bw(base_size = 16) + theme(legend.position="bottom")

# Beta diversity -------------

# First, let's look at Bray-curtis dissimilarity at the individual sample level
bray_est <- simplifyBeta(dv_pop_comparison, pop_comparison, "bray-curtis", "SamplePop")
write.table(bray_est, file=paste0(outputdir,"bray_diversity_estimate.txt"))

# add in group comparisons and plot
bray_est$group=paste(bray_est$Covar1,bray_est$Covar2,sep="_")
bray_est$Int <- interaction(bray_est$Covar1, bray_est$Covar2)
bray_est$Int=factor(bray_est$Int, levels = c("UK.Mali", "Mali.Indonesia", "Mali.Mali", "Indonesia.Indonesia", "UK.Indonesia", "UK.UK"))

# plot results
ggplot(data= bray_est, aes(x= Int, y=beta_est, fill=group)) + geom_violin(alpha=0.7) + geom_boxplot(width=0.1) + xlab("Population Comparisons") + theme_bw(base_size = 18) +
  theme(legend.position="none") + ggtitle("Bray-Curtis Distance Estimate") + 
  ylab("Bray-Curtis Distance") + scale_fill_manual(values = c("#228833","#EE7733",MaliCol,IndonesiaCol,"#EE3377",UKCol))

