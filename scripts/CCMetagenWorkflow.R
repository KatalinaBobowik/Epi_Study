# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 
# help from this pipeline: https://f1000research.com/articles/5-1492/v1

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)
library(ggpubr)
library(vegan)
library(BiocStyle)
library(knitr)
library(rmarkdown)
library(ggplot2)
library(gridExtra)
library(dada2)
library(DECIPHER)
library(phangorn)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# Import data and convert to a phyloseq object
# raw_CCMetagen_data <- read.csv(paste0(inputdir,"IndoUnmapped_species_table_noRepeats_RPM.csv"),check.names=FALSE)
raw_CCMetagen_data <- read.csv(paste0(inputdir,"IndoUnmapped_species_table_noRepeats_RPM.csv"),check.names=FALSE)

# Data preprocessing ------------------------------------

# remove plants from the analysis
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Kingdom %in% c("Viridiplantae","Fungi")),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Phylum %in% c("Bacillariophyta","Cnidaria","Echinodermata","Mollusca")),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Order %in% "Albuginales"),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(raw_CCMetagen_data[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
#rownames(taxa_raw) <- taxa_raw[,"Family"]
abund_raw <- as.matrix(raw_CCMetagen_data[,-which(colnames(raw_CCMetagen_data) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
#rownames(abund_raw) <- CCMetagen_data[,which(names(CCMetagen_data)=="Family")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

CCMeta_physeq = phyloseq(taxa, tax)
samplenames <- colnames(otu_table(CCMeta_physeq))
island <- sapply(strsplit(samplenames, "[-.]"), `[`, 1)

samples_df=data.frame(SampleName=colnames(otu_table(CCMeta_physeq)), SamplePop=island)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(CCMeta_physeq) <- samples

# do the same for count data ------
count_data <- read.csv(paste0(inputdir,"Indo_noFiltering_Subsampled_Counts.csv"),check.names=FALSE)

# remove plants from the analysis
count_data=count_data[-which(count_data$Kingdom %in% c("Viridiplantae","Fungi")),]
count_data=count_data[-which(count_data$Phylum %in% c("Bacillariophyta","Cnidaria","Echinodermata","Mollusca")),]
count_data=count_data[-which(count_data$Order %in% "Albuginales"),]
count_data=count_data[-which(count_data$Phylum %in% c("Chordata","Arthropoda")),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(count_data[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
#rownames(taxa_raw) <- taxa_raw[,"Family"]
abund_raw <- as.matrix(count_data[,-which(colnames(count_data) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
#rownames(abund_raw) <- CCMetagen_data[,which(names(CCMetagen_data)=="Family")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

CCMeta_Counts = phyloseq(taxa, tax)
samplenames <- colnames(otu_table(CCMeta_Counts))
island <- sapply(strsplit(samplenames, "[-.]"), `[`, 1)

samples_df=data.frame(SampleName=colnames(otu_table(CCMeta_Counts)), SamplePop=island)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(CCMeta_Counts) <- samples


#############################################

# get summary of each taxonomic rank
plot_taxa_summary(CCMeta_physeq, "Phylum")
summarize_taxa(CCMeta_physeq, "Phylum")

# Define prevalence of each taxa
# (in how many samples did each taxa appear at least once)
prev0 = apply(X = otu_table(CCMeta_physeq),
                MARGIN = ifelse(taxa_are_rows(CCMeta_physeq), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                      TotalAbundance = taxa_sums(CCMeta_physeq),
                      tax_table(CCMeta_physeq))

# QC pipeline
Step 1: Look at taxa prevalence
ggplot(prevdf, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum)

Step 2: Rarefaction curves

Step 3: Filter data if necessary

Step 5: Convert to RPM
It is usually necessary to transform microbiome count data to account for differences in library size, variance, scale, etc.

Step 6: Plot abundance of reads (do this in counts and RPM)
# Check distribution of how many reads/samples?
SeqDepth = colSums(otu_table(CCMeta_Counts))
sample_data(CCMeta_Counts)$SeqDepth = SeqDepth
qplot(log10(SeqDepth), geom = "histogram") + theme_bw()
# do this without the log transform
qplot(SeqDepth, geom = "histogram") + theme_bw()

# Make a histogra of the depth distribution for both groups
ggplot(meta(CCMeta_Counts)) +
    geom_histogram(aes(x = log10(SeqDepth)), alpha= 0.6) + facet_wrap(~SamplePop) + theme_bw()

# without the log transform
ggplot(meta(CCMeta_Counts)) +
    geom_histogram(aes(x = SeqDepth), alpha= 0.6) + facet_wrap(~SamplePop) + theme_bw()

# barplot of library sizes
ggbarplot(meta(CCMeta_Counts), "SampleName", "SeqDepth", fill = "SamplePop") + rotate_x_text()

What about the total reads per sample, and what does the distribution look like?

readsumsdf = data.frame(nreads = sort(taxa_sums(CCMeta_Counts), TRUE), sorted = 1:ntaxa(CCMeta_Counts), 
    type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(CCMeta_Counts), 
    TRUE), sorted = 1:nsamples(CCMeta_Counts), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")


Step 7: look at principal coordinates analysis (PCoA) for outliers
keepSamples=names(which(colSums(otu_table(CCMeta_physeq))>0))
ps=prune_samples(keepSamples, CCMeta_physeq)

sample_data(pslog)$age_binned <- cut(sample_data(pslog)$age,
  				          breaks = c(0, 100, 200, 400))
out.wuf.log <- ordinate(ps, method = "PCoA", distance = "euclidean")
plot_ordination(ps, out.wuf.log)

evals <- out.wuf.log$values$Eigenvalues
plot_ordination(ps, out.wuf.log, color = "SamplePop") +
  labs(col = "SamplePop") +
  coord_fixed(sqrt(evals[2] / evals[1]))

plot_ordination(ps, out.wuf.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))

# alpha diversity
plot_richness(CCMeta_Counts, x="SamplePop", measures=c("Shannon", "Simpson"), color="SamplePop")

### functions

library("phyloseq")
library("data.table")
library("ggplot2")

# from here: https://github.com/joey711/phyloseq/issues/418
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(meanRA = mean(RelativeAbundance),
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

plot_taxa_summary = function(physeq, Rank, GroupBy = NULL){
  # Get taxa summary table 
  dt1 = summarize_taxa(physeq, Rank = Rank, GroupBy = GroupBy)
  # Set factor appropriately for plotting
  RankCol = which(colnames(dt1) == Rank)
  setorder(dt1, -meanRA)
  dt1[, RankFac := factor(dt1[[Rank]], 
                          levels = rev(dt1[[Rank]]))]
  dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
  dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
  # Set zeroes to one-tenth the smallest value
  ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
  ebarMinFloor <- ebarMinFloor / 10
  dt1[(ebarMin == 0), ebarMin := ebarMinFloor]

  pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
    scale_x_log10() +
    xlab("Mean Relative Abundance") +
    ylab(Rank) +
    theme_bw()
  if(!is.null(GroupBy)){
    # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
    pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                size = 5)
  } else {
    # Don't include error bars for faceted version
    pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                        xmin = ebarMin))
  }
  return(pRank)
}

# Test
data("GlobalPatterns")
plot_taxa_summary(GlobalPatterns, "Phylum")
summarize_taxa(CCMeta_physeq_Controls, "Phylum")


