# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(vegan)
library(ggfortify)
library(microbiome)
library(microbiomeutilities)
library(viridis)
library(tibble)
library(knitr)
library(data.table)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/100KSamples/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/100KSamples/"

# load in Indonesian and control data

# Indonesian Count
load(paste0(inputdir,"subSampled_100K_Indo_Counts_physeq.Rda"))
# Indonesian RPM
load(paste0(inputdir,"subSampled_100K_Indo_RPM_physeq.Rda"))
# control counts
load(paste0(inputdir,"control_100K_Counts_physeq.Rda"))
# control RPM
load(paste0(inputdir,"control_100K_RPM_physeq.Rda"))

#########################
### Dataset filtering ###
#########################

# merge both phyloseq count objects together
merged_phylo_counts=merge_phyloseq(subSampled_100K_Indo_Counts_physeq, control_100K_Counts_physeq)
merged_phylo_RPM=merge_phyloseq(subSampled_100K_Indo_RPM_physeq, control_100K_RPM_physeq)

tdt = data.table(tax_table(merged_phylo_counts),
                 TotalCounts = taxa_sums(merged_phylo_counts),
                 OTU = taxa_names(merged_phylo_counts), 
                 SamplePop=meta(merged_phylo_counts)$SamplePop)

# save this table
write.table(tdt,file=paste0(outputdir,"OTU_dataTable_Summary.txt"),sep="\t")

# save a histogram of the number of counts
pdf(paste0(outputdir,"histogram_counts_noFiltering.pdf"), width=15)
ggplot(tdt, aes(log10(TotalCounts))) + 
  geom_histogram() + facet_wrap(~SamplePop) +
  ggtitle("Histogram of Total Counts")
dev.off()

pdf(paste0(outputdir,"histogram_log10counts_noFiltering.pdf"), width=15)
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + facet_wrap(~SamplePop) +
  ggtitle("Histogram of Total Counts")
dev.off()

sum(tdt$TotalCounts)
# [1] 463417

# summarise the dataset before filtering
summarize_phyloseq(merged_phylo_counts)
# [1] "1] Min. number of reads = 466"
# [1] "2] Max. number of reads = 30322"
# [1] "3] Total number of reads = 463417"
# [1] "4] Average number of reads = 2694.28488372093"
# [1] "5] Median number of reads = 2278.5"
# [1] "7] Sparsity = 0.963604495709489"
# [1] "6] Any OTU sum to 1 or less? YES"
# [1] "8] Number of singletons = 306"
# [1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)20.4545454545455"
# [1] "10] Number of sample variables are: 2"

# How many OTUs still survive in dutch and Indonesian samples?
dutch_prefilter=as.character(sample_data(merged_phylo_counts)$SampleName[which(sample_data(merged_phylo_counts)$SamplePop=="Netherlands")])
dutch_prefilter=prune_samples(dutch_prefilter,merged_phylo_counts)
dutch_prefilter <- prune_taxa(taxa_sums(dutch_prefilter) > 0, dutch_prefilter)
# otu_table()   OTU Table:         [ 469 taxa and 49 samples ]
# sample_data() Sample Data:       [ 49 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 469 taxa by 8 taxonomic ranks ]
indo_prefilter=as.character(sample_data(merged_phylo_counts)$SampleName[which(sample_data(merged_phylo_counts)$SamplePop=="Indonesia")])
indo_prefilter=prune_samples(indo_prefilter,merged_phylo_counts)
indo_prefilter <- prune_taxa(taxa_sums(indo_prefilter) > 0, indo_prefilter) 
# otu_table()   OTU Table:         [ 1496 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 1496 taxa by 8 taxonomic ranks ]

# make summary of the number of reads for Dutch, Indonesian, and all samples
nOTUs_prefilter=c(nrow(tax_table(merged_phylo_counts)), nrow(tax_table(indo_prefilter)), nrow(tax_table(dutch_prefilter)))
nReads_prefilter=c(sum(colSums(otu_table(merged_phylo_counts))),sum(colSums(otu_table(indo_prefilter))),sum(colSums(otu_table(dutch_prefilter))))

# make a filtering summary that we will append to later for the number of OTUs and reads at each step of filtering
filtering_summary=data.frame(OTUs_Prefiltering=nOTUs_prefilter,Reads_Prefiltering=nReads_prefilter)
rownames(filtering_summary) = c("All","Indonesian","Dutch")

# Now filter out Viridiplantae and metazoa, as they are not of interest
merged_phylo_counts=subset_taxa(merged_phylo_counts, (Kingdom!="Viridiplantae"))
merged_phylo_counts=subset_taxa(merged_phylo_counts, (Kingdom!="Metazoa"))
# Remove Metazoa and Viridiplantae for RPM
merged_phylo_RPM=subset_taxa(merged_phylo_RPM, (Kingdom!="Viridiplantae"))
merged_phylo_RPM=subset_taxa(merged_phylo_RPM, (Kingdom!="Metazoa"))

# save the unfiltered data for downstream use
merged_phylo_counts_unfiltered=merged_phylo_counts

# How many OTUs still survive in dutch and Indonesian samples?
dutch_noPlants_noChordates=as.character(sample_data(merged_phylo_counts)$SampleName[which(sample_data(merged_phylo_counts)$SamplePop=="Netherlands")])
dutch_noPlants_noChordates=prune_samples(dutch_noPlants_noChordates,merged_phylo_counts)
dutch_noPlants_noChordates<- prune_taxa(taxa_sums(dutch_noPlants_noChordates) > 0, dutch_noPlants_noChordates)
# otu_table()   OTU Table:         [ 469 taxa and 49 samples ]
# sample_data() Sample Data:       [ 49 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 469 taxa by 8 taxonomic ranks ]
indo_noPlants_noChordates=as.character(sample_data(merged_phylo_counts)$SampleName[which(sample_data(merged_phylo_counts)$SamplePop=="Indonesia")])
indo_noPlants_noChordates=prune_samples(indo_noPlants_noChordates,merged_phylo_counts)
indo_noPlants_noChordates <- prune_taxa(taxa_sums(indo_noPlants_noChordates) > 0, indo_noPlants_noChordates) 
# otu_table()   OTU Table:         [ 1074 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 1074 taxa by 8 taxonomic ranks ]

# save this information to the filter df
nOTUs_noChord_noMet=c(nrow(tax_table(merged_phylo_counts)), nrow(tax_table(indo_noPlants_noChordates)), nrow(tax_table(dutch_noPlants_noChordates)))
nReads_noChord_noMet=c(sum(colSums(otu_table(merged_phylo_counts))),sum(colSums(otu_table(indo_noPlants_noChordates))),sum(colSums(otu_table(dutch_noPlants_noChordates))))
filtering_summary$OTUs_noChord_noMet=nOTUs_noChord_noMet
filtering_summary$Reads_noChord_noMet=nReads_noChord_noMet

# Filter out singletons
otu_table(merged_phylo_counts)[otu_table(merged_phylo_counts)<=1]<-0
otu_table(merged_phylo_RPM)[otu_table(merged_phylo_RPM)<=1]<-0

# filter out any 0's that remain in the dataframe
any(taxa_sums(merged_phylo_counts) == 0)
# TRUE
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
merged_phylo_counts_dup=merged_phylo_counts
# summarise the phyloseq dataset
summarize_phyloseq(merged_phylo_counts)
# [1] "1] Min. number of reads = 18"
# [1] "2] Max. number of reads = 29600"
# [1] "3] Total number of reads = 212892"
# [1] "4] Average number of reads = 1237.74418604651"
# [1] "5] Median number of reads = 178.5"
# [1] "7] Sparsity = 0.973387949260042"
# [1] "6] Any OTU sum to 1 or less? NO"
# [1] "8] Number of singletons = 0"
# [1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)0"
# [1] "10] Number of sample variables are: 2"
# [1] "SampleName" "SamplePop" 

# Also, summarise phyloseq object after filtering
# otu_table()   OTU Table:         [ 660 taxa and 172 samples ]
# sample_data() Sample Data:       [ 172 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 660 taxa by 8 taxonomic ranks ]

# How many OTUs still survive in dutch and Indonesian samples?
dutch_noSingletons=as.character(sample_data(merged_phylo_counts)$SampleName[which(sample_data(merged_phylo_counts)$SamplePop=="Netherlands")])
dutch_noSingletons=prune_samples(dutch_noSingletons,merged_phylo_counts)
dutch_noSingletons <- prune_taxa(taxa_sums(dutch_noSingletons) > 0, dutch_noSingletons)
# otu_table()   OTU Table:         [ 248 taxa and 49 samples ]
# sample_data() Sample Data:       [ 49 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 248 taxa by 8 taxonomic ranks ] 
indo_noSingletons=as.character(sample_data(merged_phylo_counts)$SampleName[which(sample_data(merged_phylo_counts)$SamplePop=="Indonesia")])
indo_noSingletons=prune_samples(indo_noSingletons,merged_phylo_counts)
indo_noSingletons <- prune_taxa(taxa_sums(indo_noSingletons) > 0, indo_noSingletons) 
# otu_table()   OTU Table:         [ 534 taxa and 123 samples ]
# sample_data() Sample Data:       [ 123 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ 534 taxa by 8 taxonomic ranks ]

# add this information to the filtering table
nOTUs_noSingletons=c(nrow(tax_table(merged_phylo_counts)), nrow(tax_table(indo_noSingletons)), nrow(tax_table(dutch_noSingletons)))
nReads_noSingletons=c(sum(colSums(otu_table(merged_phylo_counts))),sum(colSums(otu_table(indo_noSingletons))),sum(colSums(otu_table(dutch_noSingletons))))
filtering_summary$OTUs_noSingletons=nOTUs_noSingletons
filtering_summary$Reads_noSingletons=nReads_noSingletons
# save table
write.table(filtering_summary,file=paste0(outputdir,"Filtering_Summary.txt"),sep="\t")

# make this into a figure
filtering_summary=t(filtering_summary)
filtering_summary=as.data.frame(filtering_summary)
stage=sapply(strsplit(rownames(filtering_summary), "[_.]"), `[`, 2)
type=sapply(strsplit(rownames(filtering_summary), "[_.]"), `[`, 1)
filtering_summary$Stage=stage
filtering_summary$Type=type
melted_filtering_summary=melt(filtering_summary)
melted_filtering_summary$Stage <- factor(melted_filtering_summary$Stage, levels = c("Prefiltering", "noChord", "noSingletons"))

# p<-ggplot(melted_a, aes(x=stage, y=value, group=variable)) +
#   geom_line(aes(color=variable))+
#   geom_point(aes(color=variable)) + facet_wrap(~Type, scales = "free")

pdf(paste0(outputdir,"filtering_summary.pdf"), width=15)
ggplot(data=melted_filtering_summary, aes(x=stage, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired") + facet_wrap(~Type, scales = "free")
dev.off()

# phylo counts summary table after filtering
tdt = data.table(tax_table(merged_phylo_counts),
                 TotalCounts = taxa_sums(merged_phylo_counts),
                 OTU = taxa_names(merged_phylo_counts), 
                 SamplePop=meta(merged_phylo_counts)$SamplePop)

write.table(tdt,file=paste0(outputdir,"OTU_dataTable_Summary_AfterFilteringSingletons.txt"),sep="\t")

pdf(paste0(outputdir,"histogram_counts_afterFiltering.pdf"), width=15)
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + facet_wrap(~SamplePop) +
  ggtitle("Histogram of Total Counts")
dev.off()

pdf(paste0(outputdir,"histogram_counts_afterFiltering_logged.pdf"), width=15)
ggplot(tdt, aes(log10(TotalCounts))) + 
  geom_histogram() + facet_wrap(~SamplePop) +
  ggtitle("Histogram of Total Counts")
dev.off()

# Check distribution of how many reads/samples?
SeqDepth = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth = SeqDepth

# Make a histogra of the depth distribution for both groups
pdf(paste0(outputdir,"histogram_seqDepth_afterFiltering.pdf"), width=15)
ggplot(meta(merged_phylo_counts)) +
    geom_histogram(aes(x = SeqDepth), alpha= 0.6) + facet_wrap(~SamplePop)
dev.off()

# log-transformed data
pdf(paste0(outputdir,"histogram_seqDepth_afterFiltering_logged.pdf"), width=15)
ggplot(meta(merged_phylo_counts)) +
    geom_histogram(aes(x = log10(SeqDepth)), alpha= 0.6) + facet_wrap(~SamplePop)
dev.off()

# get the minimum and maximum sequencing depth
min(SeqDepth)
# 18
max(SeqDepth)
# 29600

# barplot of library sizes
pdf(paste0(outputdir,"barplot_librarySize_afterFiltering.pdf"), width=15)
ggbarplot(meta(merged_phylo_counts), "SampleName", "SeqDepth", fill = "SamplePop") + rotate_x_text()
dev.off()

# get barplot of total counts per individual
nOTUs = colSums(otu_table(merged_phylo_counts)!=0)
sample_data(merged_phylo_counts)$nOTUs = nOTUs

# barplot of library sizes
pdf(paste0(outputdir,"barplot_nOTUs_afterFiltering.pdf"), width=15)
ggbarplot(meta(merged_phylo_counts), "SampleName", "nOTUs", fill = "SamplePop") + rotate_x_text()
dev.off()

# Any empty samples?
=# character(0)

#########################
### Exploratory plots ###
#########################

mergedlog <- transform_sample_counts(merged_phylo_counts, function(x) log(1 + x))
out.wuf.log <- ordinate(mergedlog, method = "PCoA", distance = "euclidean")
pdf(paste0(outputdir,"SamplePop_PCoA_Euclidean_afterFiltering.pdf"), width=10)
plot_ordination(mergedlog, out.wuf.log, color="SamplePop", axes = 1:2, label="SampleName")
plot_ordination(mergedlog, out.wuf.log, color="SamplePop", axes = 2:3, label="SampleName")
plot_ordination(mergedlog, out.wuf.log, color="SamplePop", axes = 3:4, label="SampleName")

dev.off()

# if you look at the third and fourth dimension, clustering appears by Plasmodium samples
plot_ordination(mergedlog, out.wuf.log, color="SamplePop", axes = 3:4, label="SampleName")

# Let's highlight all taxa by their phylum to see what's going on
all_phylum=unique(tax_table(merged_phylo_counts)[,"Phylum"])

for (phylum in all_phylum){
	logged_phyla_counts = log10(colSums(otu_table(merged_phylo_counts)[grep(phylum,tax_table(merged_phylo_counts)[,"Phylum"])])+1)
	sample_data(merged_phylo_counts)[[phylum]] = logged_phyla_counts
	mergedlog <- transform_sample_counts(merged_phylo_counts, function(x) log(1 + x))
	out.wuf.log <- ordinate(mergedlog, method = "PCoA", distance = "euclidean")
	pdf(paste0(outputdir,phylum,"_PCoA_Euclidean_afterFiltering.pdf"), width=10)
	print(plot_ordination(mergedlog, out.wuf.log, color=phylum, axes = 1:2, label="SampleName"))
	print(plot_ordination(mergedlog, out.wuf.log, color=phylum, axes = 2:3, label="SampleName"))
	print(plot_ordination(mergedlog, out.wuf.log, color=phylum, axes = 3:4, label="SampleName"))
	dev.off()
}

# You can also plot samples by taxa
pdf(paste0(outputdir,"OTU_PCoA_Euclidean.pdf"), width=10)
p1=plot_ordination(mergedlog, out.wuf.log, type="taxa", color="Phylum", title="taxa", axes = 3:4)
dev.off()
# separate this out by using facet_wrap
pdf(paste0(outputdir,"OTU_PCoA_Euclidean.pdf"), width=10)
p1 + facet_wrap(~Phylum, 3)
dev.off()


# Bray-Curtis dissimilarity.
out.pcoa.log <- ordinate(mergedlog,  method = "MDS", distance = "bray")
pdf(paste0(outputdir,"SampleDistance_MDS_BrayDissimilarity.pdf"), width=10)
print(plot_ordination(mergedlog, out.pcoa.log, color="SamplePop", axes = 1:2, label="SampleName"))
print(plot_ordination(mergedlog, out.pcoa.log, color="SamplePop", axes = 2:3, label="SampleName"))
print(plot_ordination(mergedlog, out.pcoa.log, color="SamplePop", axes = 3:4, label="SampleName"))
dev.off()

# We can also look at how prevalent each OTU is. (Note, this doesn't work at the kingdom level but it does work for lower levels)
prev.otu <- plot_taxa_prevalence(merged_phylo_counts, "Phylum")
pdf(paste0(outputdir,"OTU_Prevalence.pdf"), width=10)
print(prev.otu)
dev.off()

#######################
### Alpha Diversity ###
#######################

# Of specific interest in microbial ecology is the diversity of microbial communities.
# First we will calculate eveness. We will use unfiltered data with singletons included.

pAlpha = plot_richness(merged_phylo_counts_unfiltered,
                       color = "SamplePop",
                       measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                       title = "Alpha Diveristy")
pdf(paste0(outputdir,"AlphaDiversity_AllMeasures.pdf"), width=10)
print(pAlpha)
dev.off()


## Species evenness refers to how close in numbers each species in an environment is. We'll calulate the species evnness
# by using the Simpson evenness metric
ps.even <- evenness(merged_phylo_counts, index = "all") 
kable(head(ps.even))
ps1a.meta <- meta(merged_phylo_counts)
ps1a.meta$simpson <- ps.even$simpson 
hist(ps1a.meta$simpson)
shapiro.test(ps1a.meta$simpson)
qqnorm(ps1a.meta$simpson)


p1 <- ggviolin(ps1a.meta, x = "SamplePop", y = "simpson",
 add = "boxplot", fill = "SamplePop", 
 palette = c("#a6cee3", "#b2df8a"),
 legend = "right") 
p1 <- p1 + stat_compare_means(method = "t.test", label.y=1.25, label.x=0.5)

pdf(paste0(outputdir,"SimpsonEvenness.pdf"), width=10)
print(p1)
dev.off()

######################
### Beta Diversity ###
######################

# we will keep only those OTUs that are detected at least once in 5 out of total 14 samples
ps4 <- core(merged_phylo_counts_unfiltered, detection = 1, prevalence = 2/nsamples(merged_phylo_counts))

ps4.rel <- microbiome::transform(ps4, "compositional")
bx.ord_pcoa_bray <- ordinate(ps4.rel, "PCoA", "bray")

# Make an ordination plot using bray's dissimilarity
beta.ps1 <- plot_ordination(ps4.rel, 
                            bx.ord_pcoa_bray, 
                            color="SamplePop", 
                            label = "SampleName") + 
  geom_point(aes(), size= 4) + 
  theme(plot.title = element_text(hjust = 0, size = 12))


# add in an ellipse
beta.ps1 + scale_color_brewer(palette = "Dark2") + stat_ellipse() + theme_bw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# get a test on dispersion
metadf.bx <- data.frame(sample_data(ps4.rel))
bray_ps.bxn <- phyloseq::distance(physeq = ps4.rel, method = "bray")
adonis.test <- adonis(bray_ps.bxn ~ SamplePop, data = metadf.bx)
adonis.test
dist <- vegdist(t(abundances(ps4.rel)))
anova(betadisper(dist, metadf.bx$SamplePop))

# Quantifying group divergence / spread ---------------------

# Divergence of a given sample set can be quantified as the average dissimilarity of each sample 
# from the group mean; the dissimilarity can be quantified by beta diversity, for instance. This was 
# applied in group-level comparisons for instance in Salonen et al. ISME J 2014. 
# They focused on homogeneity using inverse correlation, whereas here we focus on divergence using 
# correlation but the measure is essentially the same. For more information, check Beta diversity and 
# microbiome divergence

#Calculate group divergences within the Biopsy and Stool samples
b.st <- as.data.frame(divergence(subset_samples(ps4, SamplePop == "Indonesia")))
b.bx <- as.data.frame(divergence(subset_samples(ps4, SamplePop == "Netherlands")))

# Plot the divergence
dif.a <- reshape2::melt(b.st)
dif.a$variable="Indonesia"
dif.b <- reshape2::melt(b.bx)
dif.b$variable="Netherlands"
dif.g <- rbind(dif.a,dif.b)
ggpubr::ggboxplot(dif.g, "variable", "value", 
                  ylab = "Divergence", 
                  xlab = "Sample Type", 
                  add = "jitter",
                  fill = "variable",
                  palette = c("#a6cee3", "#b2df8a"))

#############################
### Composition Analysis ###
#############################

# to do: get regular barplots of samples

ps1a.com <- merged_phylo_counts
taxic <- as.data.frame(ps1a.com@tax_table) 

# Add the OTU ids from OTU table into the taxa table at the end.
taxic$OTU <- rownames(taxic) 

# You can see that we now have extra taxonomy levels.
colnames(taxic)

# convert it into a matrix.
taxmat <- as.matrix(taxic)

# convert into phyloseq compaitble file.
new.tax <- tax_table(taxmat)  

# incroporate into phyloseq Object
tax_table(ps1a.com) <- new.tax 
pseq.ph <- aggregate_top_taxa(ps1a.com, "Family", top = 11)
p.phy <- plot_composition(pseq.ph, sample.sort = NULL, otu.sort = NULL,
  x.label = "SamplePop", plot.type = "barplot", verbose = FALSE)

print(p.phy + scale_fill_brewer(palette = "Paired") + theme_bw() + rotate_x_text())
# it would be nice to have the Taxonomic names in italics.
# for that we set this
guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15, 
    face = "italic", colour = "Black", angle = 0)))

pseq.ph.rel <- microbiome::transform(pseq.ph, "compositional")

plot.comp.rel <- plot_composition(pseq.ph.rel, x.label = "SampleType") + 
  theme(legend.position = "bottom") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Relative abundance") + guide_italics + 
  theme(legend.title = element_text(size=18))

plot.comp.rel + scale_fill_brewer( "Family",palette = "Paired")

# Boxplots

# Check what the top ten genera are and how they differ between sample types.
pn <- plot_taxa_boxplot(ps4.rel, "Family", 10, "SamplePop", color = "Set2", "Relative abundance of top 10 genera")
print(pn)




