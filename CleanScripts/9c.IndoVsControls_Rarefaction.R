
# title: 4.IndoSamples_EffectOfSequencingDepth.R
# author: katalinabobowik
# date: '2020-09-13'


# Loading packages and colour setup

# The code below will install the packages needed to run the analyses. We're also setting
# up our directories to run locally, but I've provided all of the paths so they 
# can be easily modified if they need to be run on the server.

require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(metacoder)
library(tidyverse)             
library(phyloseq)                   
library(DESeq2)                       
library(microbiome)               
library(vegan)                         
library(picante)                     
library(ALDEx2)                      
library(metagenomeSeq)          
library(HMP)                             
library(dendextend)               
library(selbal)                       
library(rms)
library(breakaway)        
library(microbiomeutilities)
library(mixOmics)
library(SRS)

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
AllReadsdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/All_Reads/"
Indo250Kdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Indo_250K/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/ControlSampleComparison/"
refdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"

# set ggplot colour theme to white
theme_set(theme_bw())

# Set up colour scheme
IndonesiaCol="#4477AA"
NetherlandsCol="#EE6677"

# Read in the data
AllREadsSE_Indo_Counts <- read.csv(paste0(refdir,"Counts_Indo_AllReads_SE.csv"),check.names=FALSE)

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(AllREadsSE_Indo_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(AllREadsSE_Indo_Counts[,-which(colnames(AllREadsSE_Indo_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
AllREadsSE_Indo_Counts_physeq = phyloseq(taxa, tax)
# add in sample information, starting with Island
samplenames <- colnames(otu_table(AllREadsSE_Indo_Counts_physeq))
pop <- rep("Indonesia",ncol(otu_table(AllREadsSE_Indo_Counts_physeq)))
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(AllREadsSE_Indo_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(AllREadsSE_Indo_Counts_physeq) <- samples

# get information on phyloseq objects
AllREadsSE_Indo_Counts_physeq
# get replicates
trimmed_samplenames = gsub("Batch1",'',samplenames) %>% gsub("Batch2",'', .) %>% gsub("Batch3",'', .)
trimmed_samplenames = sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", trimmed_samplenames)
replicate_index = which(duplicated(trimmed_samplenames) | duplicated(trimmed_samplenames, fromLast = TRUE))
replicates = samplenames[replicate_index]

# add sequencing depth information to the Physeq object in order to filter replicates by seqDepth
SeqDepth = colSums(otu_table(AllREadsSE_Indo_Counts_physeq))
sample_data(AllREadsSE_Indo_Counts_physeq)$SeqDepth = SeqDepth

# find out which replicates have the highest sequencing depth
sample_data(AllREadsSE_Indo_Counts_physeq)[replicates,]
replicateDF=as.data.frame(sample_data(AllREadsSE_Indo_Counts_physeq)[replicates,])
replicateDF$SampleName = sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", replicateDF$SampleName)
replicateDF$SampleName = gsub("Batch1",'',replicateDF$SampleName) %>% gsub("Batch2",'', .) %>% gsub("Batch3",'', .)
replicateDF=replicateDF[with(replicateDF, order(-SeqDepth)), ]
removeReplicates=rownames(replicateDF[which(duplicated(replicateDF$SampleName)),])
keepReplicates=rownames(sample_data(AllREadsSE_Indo_Counts_physeq))[-which(rownames(sample_data(AllREadsSE_Indo_Counts_physeq)) %in% removeReplicates)]

# prune these out
AllREadsSE_Indo_Counts_physeq=prune_samples(keepReplicates,AllREadsSE_Indo_Counts_physeq)
AllREadsSE_Indo_Counts_physeq

# Control samples
control_100K_Counts <- read.csv(paste0(refdir,"Controls_NoFiltering_Counts.csv"),check.names=FALSE)
# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(control_100K_Counts[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(control_100K_Counts[,-which(colnames(control_100K_Counts) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
control_100K_Counts_physeq = phyloseq(taxa, tax)
control_100K_Counts_physeq
# text file made from script 'GetDiseaseStatus_ControlSamples.sh'
diseaseStatus=read.table(paste0(refdir,"ControlSamples_DiseaseStatus.txt"))
colnames(diseaseStatus)=c("Samples","diseaseStatus")
controlSamples=as.character(diseaseStatus[grep("CSF", diseaseStatus$diseaseStatus),"Samples"])
# prune out samples we don't want
control_100K_Counts_physeq=prune_samples(controlSamples,control_100K_Counts_physeq)

# add in sample information, i.e., the sample names and population they're from
samplenames <- colnames(otu_table(control_100K_Counts_physeq))
pop <- rep("Netherlands",ncol(otu_table(control_100K_Counts_physeq)))
# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(control_100K_Counts_physeq)), SamplePop=pop)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(control_100K_Counts_physeq) <- samples

control_100K_Counts_physeq

# merge the data together
merged_phylo_counts=merge_phyloseq(AllREadsSE_Indo_Counts_physeq, control_100K_Counts_physeq)
merged_phylo_counts

# summarise the dataset 
summarize_phyloseq(merged_phylo_counts)

# Data processing
SeqDepth = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth = SeqDepth

# log-transformed data
ggplot(meta(merged_phylo_counts)) +
    geom_histogram(aes(x = SeqDepth), alpha= 0.6) + facet_wrap(~SamplePop)

# Separate species' abundances and taxonomy columns
rarecurve_counts <- otu_table(merged_phylo_counts)

# Try with different filtering thresholds:
for (i in c(1,5)){
 	rarecurve_counts[rarecurve_counts<=i]<-0
	rarecurve(t(otu_table(rarecurve_counts, taxa_are_rows = TRUE)), step=20, col=c(rep(IndonesiaCol,sum(sample_data(merged_phylo_counts)[,"SamplePop"]=="Indonesia")),rep(NetherlandsCol,sum(sample_data(merged_phylo_counts)[,"SamplePop"]=="Netherlands"))),label=F, xlab="Counts",ylab="Number Species",main=paste("Removing Reads",i,"And Below",sep=" "),xlim=c(0,50000))
}

# Filter out singletons
merged_phylo_counts_unfiltered=merged_phylo_counts
otu_table(merged_phylo_counts)[otu_table(merged_phylo_counts)<=1]<-0

# filter out any 0's that remain in the dataframe
any(taxa_sums(merged_phylo_counts) == 0)
# TRUE
merged_phylo_counts <- prune_taxa(taxa_sums(merged_phylo_counts) > 0, merged_phylo_counts)
# summarise the phyloseq dataset
summarize_phyloseq(merged_phylo_counts)

# Summary stats after filtering
# Check distribution of how many reads/samples?
SeqDepth = colSums(otu_table(merged_phylo_counts))
sample_data(merged_phylo_counts)$SeqDepth = SeqDepth

# log-transformed data
ggplot(meta(merged_phylo_counts)) +
    geom_histogram(aes(x = log10(SeqDepth)), alpha= 0.6) + facet_wrap(~SamplePop)

# barplot of library sizes
ggplot(meta(merged_phylo_counts), aes(SampleName, SeqDepth)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
scale_fill_manual(values = c(IndonesiaCol,NetherlandsCol)) + rotate_x_text()

# get barplot of total counts per individual
nOTUs = colSums(otu_table(merged_phylo_counts)!=0)
sample_data(merged_phylo_counts)$nOTUs = nOTUs

# barplot of OTUs
ggplot(meta(merged_phylo_counts), aes(SampleName, nOTUs)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
scale_fill_manual(values = c(IndonesiaCol,NetherlandsCol)) + rotate_x_text()

# filter out plants and animals
merged_phylo_counts=subset_taxa(merged_phylo_counts, (Kingdom!="Viridiplantae"))
merged_phylo_counts=subset_taxa(merged_phylo_counts, (Kingdom!="Metazoa"))

# Data normalisation
# CLR
# add an offset of 1 to the counts. This is necessary, since performing a log on 0 values is undefined. We'll then perform a log ration transformation of the data from the mixOmics package.
offset_otu=otu_table(merged_phylo_counts)+1
transform_counts=t(otu_table(offset_otu))
data_clr <- logratio.transfo(as.matrix(transform_counts), logratio = 'CLR', offset = 0) 

# Exploratory plots

# Make a duplicated phyloseq object to use for plotting
class(data_clr)="matrix"
taxa = otu_table(t(data_clr), taxa_are_rows = TRUE)
merged_phylo_counts_clr=merged_phylo_counts
otu_table(merged_phylo_counts_clr)=taxa

out.wuf.log <- ordinate(merged_phylo_counts_clr, method = "PCoA", distance = "euclidean")
plot_ordination(merged_phylo_counts_clr, out.wuf.log, color="SamplePop", axes = 1:2, label="SampleName")
plot_ordination(merged_phylo_counts_clr, out.wuf.log, color="SamplePop", axes = 3:4, label="SampleName")

# Relative abundance of microbiome data

# add a new column containing family names and superkingdom
tax_table(merged_phylo_counts)[,"Superkingdom"] = paste(tax_table(merged_phylo_counts)[,"Superkingdom"], tax_table(merged_phylo_counts)[,"Family"], sep="_")
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])
tax_table(merged_phylo_counts)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(merged_phylo_counts)[,"Superkingdom"])

aggregated_phyloCounts <- aggregate_top_taxa(merged_phylo_counts, "Superkingdom", top = 19)
relative_phyloCounts <- microbiome::transform(aggregated_phyloCounts, "compositional")

# Plot
p=plot_bar(relative_phyloCounts, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))  

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(15)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(2)
PaletteUnk = colorRampPalette(c("black"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteUnk,PaletteVirus)

phyloseq::plot_bar(relative_phyloCounts, fill = "Superkingdom") +
  geom_bar(aes(fill = Superkingdom), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ SamplePop, scales = "free") + scale_fill_manual(values=Merged_Palette) +
  theme(panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Rarefy data
minControl=min(sample_sums(phyloseq::subset_samples(merged_phylo_counts, SamplePop == "Netherlands")))
keep=names(which(sample_sums(merged_phylo_counts)>=minControl))
pruned_data=prune_samples(keep, merged_phylo_counts)
any(taxa_sums(pruned_data) == 0)
# TRUE
pruned_data <- prune_taxa(taxa_sums(pruned_data) > 0, pruned_data)
any(taxa_sums(pruned_data) == 0)
# FALSE
pruned_data_df=as.data.frame(otu_table(pruned_data))
SRS=SRS(pruned_data_df,minControl)
rownames(SRS)=rownames(pruned_data_df)

# integrate back into phyloseq ibject
taxa = otu_table(SRS, taxa_are_rows = TRUE)
otu_table(pruned_data)=taxa
any(taxa_sums(pruned_data) == 0)
# TRUE
pruned_data <- prune_taxa(taxa_sums(pruned_data) > 0, pruned_data)
any(taxa_sums(pruned_data) == 0)

SeqDepthPruned = sample_sums(pruned_data)
sample_data(pruned_data)$SeqDepthPruned = SeqDepthPruned

# barplot of library sizes
ggplot(meta(pruned_data), aes(SampleName, SeqDepthPruned)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
scale_fill_manual(values = c(IndonesiaCol,NetherlandsCol)) + rotate_x_text()

# get barplot of total counts per individual
nOTUs = colSums(otu_table(pruned_data)!=0)
sample_data(pruned_data)$nOTUs = nOTUs

# barplot of OTUs
ggplot(meta(pruned_data), aes(SampleName, nOTUs)) + geom_bar(stat = "identity", aes(fill = SamplePop)) +
scale_fill_manual(values = c(IndonesiaCol,NetherlandsCol)) + rotate_x_text()

# Exploratory plots after rarefaction
mergedlog <- transform_sample_counts(pruned_data, function(x) log(1 + x))
out.wuf.log <- ordinate(mergedlog, method = "PCoA", distance = "euclidean")
plot_ordination(mergedlog, out.wuf.log, color="SamplePop", axes = 1:2, label="SampleName")
plot_ordination(mergedlog, out.wuf.log, color="SamplePop", axes = 3:4, label="SampleName")

# Let's highlight all taxa by their phylum to see what's going on
all_phylum=c("Apicomplexa", "Firmicutes", "unk_p")

for (phylum in all_phylum){
	logged_phyla_counts = log10(colSums(otu_table(merged_phylo_counts)[grep(phylum,tax_table(merged_phylo_counts)[,"Phylum"])])+1)
	sample_data(merged_phylo_counts)[[phylum]] = logged_phyla_counts
	mergedlog <- transform_sample_counts(merged_phylo_counts, function(x) log(1 + x))
	out.wuf.log <- ordinate(mergedlog, method = "PCoA", distance = "euclidean")
	print(plot_ordination(mergedlog, out.wuf.log, color=phylum, axes = 1:2, label="SampleName"))
}

p1=plot_ordination(mergedlog, out.wuf.log, type="taxa", color="Phylum", title="taxa", axes = 3:4)
# separate this out by using facet_wrap
p1 + facet_wrap(~Phylum, 3)

# add a new column containing family names and superkingdom
tax_table(pruned_data)[,"Superkingdom"] = paste(tax_table(pruned_data)[,"Superkingdom"], tax_table(pruned_data)[,"Family"], sep="_")
tax_table(pruned_data)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(pruned_data)[,"Superkingdom"])
tax_table(pruned_data)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(pruned_data)[,"Superkingdom"])
tax_table(pruned_data)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(pruned_data)[,"Superkingdom"])

aggregated_phyloCounts <- aggregate_top_taxa(pruned_data, "Superkingdom", top = 19)
relative_phyloCounts <- microbiome::transform(aggregated_phyloCounts, "compositional")

# Plot
p=plot_bar(relative_phyloCounts, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
   
PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(16)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(2)
PaletteUnk = colorRampPalette(c("black"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(1)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteUnk,PaletteVirus)

pdf(paste0(outputdir,"relative_phyloCounts.pdf"), width=15)
phyloseq::plot_bar(relative_phyloCounts, fill = "Superkingdom") +
  geom_bar(aes(fill = Superkingdom), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ SamplePop, scales = "free") + scale_fill_manual(values=Merged_Palette) +
  theme(panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

# Boxplots
phyloseq::psmelt(pruned_data) %>%
ggplot(data = ., aes(x = SamplePop, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free")

# Hierarchical clustering
ps_rel_abund = phyloseq::transform_sample_counts(pruned_data, function(x){x / sum(x)})
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c(Indonesia = IndonesiaCol, `Netherlands` = NetherlandsCol)
labels_colors(ward) <- colorCode[meta$SamplePop][order.dendrogram(ward)]
#Plot
pdf(paste0(outputdir,"hierarchicalClustering.pdf"), width=15)
plot(ward)
dev.off()

# Differnetial abundance testing
differential_abund=pruned_data
otu_table(differential_abund)=otu_table(differential_abund)+1
taxa_names(differential_abund)=make.unique(tax_table(differential_abund)[,"Phylum"])
diagdds = phyloseq_to_deseq2(differential_abund, ~ SamplePop)
# had help from this: https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564/2
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(differential_abund)[rownames(sigtab), ], "matrix"))
head(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
#sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
ggplot(sigtab, aes(x=log2FoldChange, y=Phylum, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) 

# Alpha diversity

#Generate a data.frame with adiv measures
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(pruned_data, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(pruned_data, measures = "Simpson"),
  "Status" = phyloseq::sample_data(pruned_data)$SamplePop)
head(adiv)

adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon","Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon","Simpson"))) %>%
  ggplot(aes(x = Status, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Status), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none") + scale_colour_manual(values=c(IndonesiaCol,NetherlandsCol))

adiv %>%
  group_by(Status) %>%
  dplyr::summarise(median_observed = median(Observed),
            median_shannon = median(Shannon),
            median_simpson = median(Simpson))

wilcox.test(Observed ~ Status, data = adiv, exact = FALSE, conf.int = TRUE)
wilcox.test(Shannon ~ Status, data = adiv, conf.int = TRUE)              
wilcox.test(Simpson ~ Status, data = adiv, conf.int = TRUE)              

# Beta Diversity
ps_clr <- microbiome::transform(pruned_data, "clr")
# RDA without constraints is PCA        
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(pruned_data, ord_clr, type="samples", color="SamplePop",label = "SampleName") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = SamplePop), linetype = 2) + scale_colour_manual(values=c(IndonesiaCol,NetherlandsCol))

ps4.rel <- microbiome::transform(pruned_data, "compositional")
bx.ord_pcoa_bray <- ordinate(ps4.rel, "PCoA", "bray")

# Make an ordination plot using bray's dissimilarity
beta.ps1 <- plot_ordination(ps4.rel, 
                            bx.ord_pcoa_bray, 
                            color="SamplePop", 
                            label = "SampleName") + 
  geom_point(aes(), size= 4) + 
  theme(plot.title = element_text(hjust = 0, size = 12))

# add in an ellipse
beta.ps1 + stat_ellipse() + theme_bw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_manual(values=c(IndonesiaCol,NetherlandsCol))

