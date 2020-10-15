# load in packages
library(breakaway)
library(DivNet)
library(tidyverse)
library(phyloseq)

# Rationale

# Phyloseq has some built in tools for exploring alpha
# diversity, they underestimate richness, underestimate uncertainty, 
# and don't allow hypothesis testing.
# Breakaway was specifically designed for these tasks

# Divnet accounts for differences in sequencing depth and estimates 
# the number of missing species based on the sequence depth and
# number of rare taxa in the data

##################################
# Divnet Estimates of Diversity #
##################################

# Divnet is a method for estimating within- and between-community diversity in ecosystems where 
# taxa interact via an ecological network

# The idea behind DivNet, is that you don't care about your samples, you care about the population that 
# your samples were drawn from. For example, if we knew the true relative abundances of all of the phyla 
# living in blood, we could calculate the true Shannon diversity and compare it to the Shannon diversity 
# we measured in our blood samples. Therefore, in order to make this comparison, we need to estimate the 
# Shannon diversity of microbes living in blood, which is what DivNet does.

# comparing diversity at the phylum level
pop_comparison <- merged_phylo_counts_withSingletons %>%
  tax_glom("Phylum")

# If we don't change the sample names here from hypehns to periods, we'll get an error later
sample_names(pop_comparison) <- gsub("\\-", ".", sample_names(pop_comparison))

# Run divnet without specifying any hypothesis testing
dv_pop_comparison <- divnet(pop_comparison, ncores = 4)

# DivNet outputs a list of the estimates shannon, simpson (alpha diversity)
# bray-curtis, euclidean (beta diversity)
dv_pop_comparison %>% names

# Now let's plot the results of shannon and Simpson diversity
summary_df_shannon <- as.data.frame(dv_pop_comparison$shannon %>%
  summary %>%
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

ggplot(summary_df_shannon, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_jitter(height = 0, width = .2) + geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + ggtitle("Shannon Diversity") +
  ylab("Estimate of Shannon Diversity")

# Simpson diversity index is a similarity index where the higher the value, the lower in diversity. It measures
# the probability that two individuals randomly selected from a sample will belong to the same species (or some category other than species).
# With this index, 0 represents infinite diversity and 1, no diversity.

summary_df_simpson <- as.data.frame(dv_pop_comparison$simpson %>%
  summary %>%
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>%
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

# Since a larger Simpson index value equates to a lower diversity index, many people find this confusing and not
# very intuitive. Therefore, the inverse Simpsone Index, or 1 - Simpson Index, is also commonly used.
# Let's plot that now. 

# Subtract the Simpson estimate from one
summary_df_simpson$estimate = 1-summary_df_simpson$estimate
# Plot
ggplot(summary_df_simpson, aes(y = estimate, x = SamplePop, fill = SamplePop)) + geom_violin(alpha=0.7) + 
  geom_jitter(height = 0, width = .2) + geom_boxplot(width=0.08, outlier.color = NA) +
  scale_fill_manual(values=c(IndonesiaCol,MaliCol,UKCol)) + ggtitle("Simpson's Diversity Index") +
  ylab("Estimate of Simpson Diversity")

# You can also plot alpha diversity indices for each individual, along with their standard error.
plot(dv_pop_comparison$shannon, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol))
plot(dv_pop_comparison$simpson, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol))

#################################
# Adding covariates into Divnet #
#################################

# Now let's test the hypothesis that the diversity is differnet between islands.
# We are now estimating the diversity of island/population being an ecosystem, so we're
# focusing on the ecosystem, not just the samples. If we want to reproduce the results of our study,
# it's better to focus on the populations that the samples come from, not the samples themselves

# test the hypothesis that the diversity is differnet between islands
dv_pop_comparison_cov <- pop_comparison %>%
  divnet(X = "SamplePop", ncores = 8)

# Plot the results for each individual
plot(dv_pop_comparison_cov$shannon, pop_comparison, col = "SamplePop") + scale_colour_manual(values=c(IndonesiaCol,MaliCol,UKCol))

# test that these populations are actually different
testDiversity(dv_pop_comparison_cov, "shannon")

# We can see that all of the populations are statistically significant from one another and that 
# the Malian population has the highest Shannon diversity. We can also see that the UK population has 0.3X less 
# diversity than the Indonesian population

# We can also do the same thing for Simpson:
testDiversity(dv_pop_comparison_cov, "simpson")
plot(dv_pop_comparison_cov$simpson, pop_comparison, col = "SamplePop")

##########################
# Testing beta diversity #
##########################

# For beta diversity, we'll plot two results - one for each sample and one with the hypothesis that
# islands are different

# First, let's look at Bray-curtis dissimilarity at the individual smple level
bray_est <- simplifyBeta(dv_pop_comparison, pop_comparison, "bray-curtis", "SamplePop")

# add in group comparisons and plot
bray_est$group=paste(bray_est$Covar1,bray_est$Covar2,sep="_")
ggplot(bray_est, aes(x = interaction(Covar1, Covar2), y = beta_est, fill=group)) +
  geom_violin(alpha=0.7) + geom_boxplot(width=0.1) + xlab("Population Comparisons") + 
  theme(legend.position="none") + ggtitle("Bray-Curtis Distance Estimate") +
  ylab("Bray-Curtis Distance")

# Now do this for island-level comaprison

bray_est_island <- simplifyBeta(dv_pop_comparison_cov, pop_comparison, "bray-curtis", "SamplePop")

# Now plot
simplifyBeta(dv_pop_comparison_cov, pop_comparison, "bray-curtis", "SamplePop") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")


#####################
# Running Breakaway #
#####################

ba <- breakaway(pop_comparison)

# Plot estimates
plot(ba, pop_comparison, color = "SamplePop")

# Take the estimates and turn them into a data frame
summary_df <- as.data.frame(summary(ba) %>% 
  add_column("SampleNames" = pop_comparison %>% otu_table %>% sample_names) %>% 
  add_column("SamplePop" =  pop_comparison %>% sample_data %>% .[,"SamplePop"] %>% as.matrix(.) %>% .[,1] %>% unname(.)))

ggplot(summary_df, aes(y = estimate, x = SamplePop, color = SamplePop)) + geom_violin()
 
# Let's test the hypothesis that different populations have the same microbial diversity
# betta() works like a regression model but it accounts for the uncertainty in estimating diversity
bt <- betta(summary(ba)$estimate, summary(ba)$error, make_design_matrix(pop_comparison, "SamplePop"))
bt$table

# betta() estimates that the mean Phylum-level diversity in Indonesians is 10.87 orders.
# It estimates that the diversity in Mali is significantly higher (on average 9.3 orders) while the
# UK population does not have significantly different diversity than the Indonesian population

