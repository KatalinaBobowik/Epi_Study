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

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/100KControls/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/IndoVsControls/"
controldir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Controls/"

# compare to Control dataset ----------------------

# read in control dataset
controlDF = read.table(paste0(controldir,"Control_Taxa_DF.txt"))
controlDF$taxa=rownames(controlDF)

# earlier, we took Plasmodium out so reset taxa to include Plasmodium
taxa2 = otu_table(abund_raw, taxa_are_rows = TRUE)
taxa=as.data.frame(taxa)
taxa$taxa=rownames(taxa)

# combine both DFs

a=left_join(taxa, controlDF, by = 'taxa')
a[is.na(a)] <- 0
taxaNames = a$taxa
a$taxa <- NULL
#a$taxa = taxaNames
a=t(a)
b=as.data.frame(a)
b$group = c(rep("Indonesian",123),rep("Dutch",49))

# plot pca
pca <- prcomp(b[1:228], scale. = F, center=T)
autoplot(pca, data = b, colour = 'group', x=1, y=2, label = TRUE)
autoplot(pca, data = b, colour = 'group', x=2, y=3, label = TRUE)
autoplot(pca, data = b, colour = 'group', x=3, y=4, label = TRUE)

# brings out plasmodium
autoplot(pca,  data = a, legend=F, label = TRUE, x=4, y=5)

# should I scale?
# As a rule of thumb, if all your variables are measured on the same scale 
# and have the same unit, it might be a good idea *not* to scale the variables 
# (i.e., PCA based on the covariance matrix). If you want to maximize variation, 
# it is fair to let variables with more variation contribute more.

# for pathogens
pca <- prcomp(a[,1:172], scale. = F, center=T)
autoplot(pca, data = a, x=1, y=2, label = TRUE)
autoplot(pca, data = a, x=2, y=3, label = TRUE)


# Tutorial ----------------------------

# load in Phyloseq objects
load(paste0(inputdir,"CCMeta_physeq_100K_AllTaxaRanks.Rda"))
load(paste0(controldir,"CCMeta_physeq_Controls_AllTaxaRanks.Rda"))

# merge both phyloseq objects together
merged_phylo=merge_phyloseq(CCMeta_physeq_100K, CCMeta_physeq_Controls)
# add in sample information
samples_df=data.frame(SampleName=colnames(otu_table(merged_phylo)), SamplePop=c(rep("Indonesian",123),rep("Dutch",49)))
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(merged_phylo) <- samples

rank_names(merged_phylo)
# [1] "Superkingdom" "Kingdom"      "Phylum"       "Class"        "Order"       
# [6] "Family"       "Genus"        "Species" 

# Check distribution of how many reads/samples?
SeqDepth = colSums(otu_table(merged_phylo))
sample_data(merged_phylo)$SeqDepth = SeqDepth
qplot(log10(SeqDepth), geom = "histogram") + theme_bw()
# do this without the log transform
qplot(SeqDepth, geom = "histogram") + theme_bw()

# Make a histogra of the depth distribution for both groups
ggplot(meta(merged_phylo)) +
    geom_histogram(aes(x = log10(SeqDepth)), alpha= 0.6) + facet_wrap(~SamplePop) + theme_bw()

# without the log transform
ggplot(meta(merged_phylo)) +
    geom_histogram(aes(x = SeqDepth), alpha= 0.6) + facet_wrap(~SamplePop) + theme_bw()

# get the minimum and maximum sequencing depth
min(SeqDepth)
# 402.836
max(SeqDepth)
# 672030.3

# barplot of library sizes
ggbarplot(meta(merged_phylo), "SampleName", "SeqDepth", fill = "SamplePop") + rotate_x_text()

# maybe you can delete this? Basically we're trimming our any singletons in our data
ps1b <- prune_taxa(taxa_sums(merged_phylo) > 1, merged_phylo)
hist(log10(taxa_sums(ps1b)))
# The data is left tailed, which is common for microbiome count data. 

# We can also look at how prevalent each OTU is. (Note, this doesn't work at the kingdom level but it does work for lower levels)
prev.otu <- plot_taxa_prevalence(ps1b, "Phylum")
print(prev.otu)

# Diversity ------------------------------------------------------

# Of specific interest in microbial ecology is the diversity of microbial communities.
# First we will calculate eveness. We will use unfiltered data with singletons.

ps.even <- evenness(merged_phylo, index = "all") 
kable(head(ps.even))
ps1a.meta <- meta(merged_phylo)
ps1a.meta$simpson <- ps.even$simpson 
hist(ps1a.meta$simpson)
shapiro.test(ps1a.meta$simpson)
qqnorm(ps1a.meta$simpson)


p1 <- ggviolin(ps1a.meta, x = "SamplePop", y = "simpson",
 add = "boxplot", fill = "SamplePop", 
 palette = c("#a6cee3", "#b2df8a"),
 legend = "right") 

print(p1)
p1 <- p1 + stat_compare_means(method = "t.test")
print(p1)

# Beta diversity ---------------

# we will keep only those OTUs that are detected at least 5 times in 5 out of total 14 samples
ps4 <- core(merged_phylo, detection = 5, prevalence = 5/nsamples(merged_phylo))
hist(log10(taxa_sums(ps4)))

ps4.rel <- microbiome::transform(ps4, "compositional")
bx.ord_pcoa_bray <- ordinate(ps4.rel, "PCoA", "bray")

#Scree plot
plot_scree(bx.ord_pcoa_bray) + theme_bw()

beta.ps1 <- plot_ordination(ps4.rel, 
                            bx.ord_pcoa_bray, 
                            color="SamplePop", 
                            label = "SampleName") + 
  geom_point(aes(shape = SamplePop), size= 4) + 
  theme(plot.title = element_text(hjust = 0, size = 12))

beta.ps1 <- beta.ps1 + theme_bw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Now we can join the biopsy and stool from same subject
beta.ps2 <- beta.ps1 + geom_line() + scale_color_brewer(palette = "Dark2")
beta.ps2

beta.ps3 <- plot_ordination(ps4.rel, 
                            bx.ord_pcoa_bray, 
                            color="SamplePop", 
                            label = "SampleName") + 
  geom_point(size= 4) + 
  theme(plot.title = element_text(hjust = 0, size = 12))

beta.ps3 <- beta.ps3 + theme_bw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

beta.ps3 + scale_color_brewer(palette = "Dark2") + stat_ellipse()

metadf.bx <- data.frame(sample_data(ps4.rel))
bray_ps.bxn <- phyloseq::distance(physeq = ps4.rel, method = "bray")

adonis.test <- adonis(bray_ps.bxn ~ SamplePop, data = metadf.bx)
adonis.test

dist <- vegdist(t(abundances(ps4.rel)))
anova(betadisper(dist, metadf.bx$SamplePop))

# Quantifying group divergence / spread ---------------------

b.st <- as.data.frame(divergence(subset_samples(ps4, SamplePop == "Indonesian")))
b.bx <- as.data.frame(divergence(subset_samples(ps4, SamplePop == "Dutch")))

# Plot the divergence
dif.a <- reshape2::melt(b.st)
dif.b <- reshape2::melt(b.bx)
dif.g <- rbind(dif.a,dif.b)

ggpubr::ggboxplot(dif.g, "variable", "value", 
                  ylab = "Divergence", 
                  xlab = "Sample Type", 
                  add = "jitter",
                  fill = "variable",
                  palette = c("#a6cee3", "#b2df8a"))

# Composition ---------------

ps1a.com <- merged_phylo

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

# Check which are the top ten genera and how they differ between sample types.

pn <- plot_taxa_boxplot(ps4.rel, "Family", 10, "SamplePop", color = "Set2", "Relative abundance of top 10 genera")
print(pn)

ps1.bx <- subset_samples(merged_phylo, SamplePop == "Indonesian")

ps1.bx.rel <- microbiome::transform(ps1.bx, "compositional")
#your original pseq/relative abundance file
#if
colnames(tax_table(ps1.bx.rel))

## [1] "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

#last column has SVs/OTU ids then you can skip the following five steps and go to aggregate_taxa step.

taxic <- as.data.frame(ps1.bx.rel@tax_table)
taxic$OTU <- row.names(taxic)
#convert it into a matrix.
taxmat <- as.matrix(taxic)

#convert into phyloseq compaitble file.
new.tax <- tax_table(taxmat)

#incroporate into phyloseq Object
tax_table(ps1.bx.rel) <- new.tax

#the presence of NA is an issue.

tax_table(ps1.bx.rel)[,"Family"][is.na(tax_table(ps1.bx.rel)[,"Family"])] <- paste0(tolower(substring("Family", 1, 1)), "__")

#at family level
ps1.bx.gen <- aggregate_taxa(ps1.bx.rel, "Family")

#Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
#(1e-3) = 0.001% abundance; change "-3" to -2 to increase to 0.01%

p <- plot_core(ps1.bx.gen, plot.type = "heatmap", 
               colours = rev(brewer.pal(10, "Spectral")),
               min.prevalence = 0.9, 
               prevalences = prevalences, 
               detections = detections) +
  xlab("Detection Threshold (Relative Abundance (%))")
print(p)

# Heatmap -------------------

ps1.c <- format_to_besthit(merged_phylo)

heat.sample <- plot_taxa_heatmap(merged_phylo, subset.top = 5,
    VariableA = "SamplePop",
    heatcolors = rev(brewer.pal(100, "Blues")),
    transformation = 'clr')




