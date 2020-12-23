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
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Indo_250K/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/Indo_250K/"

# load in 250K Indonesian phyloseq counts object
load(paste0(inputdir,"subSampled_250K_Indo_Counts_physeq.Rda"))
# load in 250K Indonesian phyloseq RPM object
load(paste0(inputdir,"subSampled_250K_Indo_RPM_physeq.Rda"))

# Remove Metazoa and Viridiplantae
subSampled_250K_Indo_Counts_physeq=subset_taxa(subSampled_250K_Indo_Counts_physeq, (Kingdom!="Viridiplantae"))
subSampled_250K_Indo_Counts_physeq=subset_taxa(subSampled_250K_Indo_Counts_physeq, (Kingdom!="Metazoa"))
# Remove Metazoa and Viridiplantae for RPM
subSampled_250K_Indo_RPM_physeq=subset_taxa(subSampled_250K_Indo_RPM_physeq, (Kingdom!="Viridiplantae"))
subSampled_250K_Indo_RPM_physeq=subset_taxa(subSampled_250K_Indo_RPM_physeq, (Kingdom!="Metazoa"))

# filter reads 1 and below
otu_table(subSampled_250K_Indo_Counts_physeq)[otu_table(subSampled_250K_Indo_Counts_physeq)<=1]<-0
otu_table(subSampled_250K_Indo_RPM_physeq)[otu_table(subSampled_250K_Indo_Counts_physeq)<=1]<-0

# Data preprocessing ------------------------------------

# Check distribution of how many reads/samples?
SeqDepth = colSums(otu_table(subSampled_250K_Indo_Counts_physeq))
sample_data(subSampled_250K_Indo_Counts_physeq)$SeqDepth = SeqDepth

# make a histogran of sequencing depth
pdf(paste0(outputdir,"HistogramSeqDepth_250KIndonesian.pdf"), width=15)
ggplot(meta(subSampled_250K_Indo_Counts_physeq)) +
    geom_histogram(aes(x = SeqDepth), alpha= 0.6)
dev.off()

# get the minimum and maximum sequencing depth
min(SeqDepth)
# 0
max(SeqDepth)
# 3951

table(colSums(otu_table(subSampled_250K_Indo_Counts_physeq)))
 #   0    2    3    4    6    8   11   12   13   15   17   22   25   37   40   60 
 #  84   15    4    1    1    1    1    1    1    1    2    1    1    1    1    1 
 
 # 128  354  532  941 3470 3951 
 #   1    1    1    1    1    1 

 # You can see that 84 samples have 0 values, and 39 have values

# which samples actually have data?
which(colSums(otu_table(subSampled_250K_Indo_Counts_physeq))>0)
#  [1] "MPI-025"           "MPI-061"           "MPI-065"          
#  [4] "MPI-296"           "MPI-334"           "MPI-335"          
#  [7] "MPI-345"           "MPI-376"           "MPI-389"          
# [10] "MTW-MDB-024"       "MTW-MDB-034"       "MTW-TLL-006"      
# [13] "MTW-TLL-007"       "MTW-TLL-012"       "MTW-TLL-016"      
# [16] "MTW-TLL-017"       "MTW-TLL-023"       "MTW-TLL-027"      
# [19] "MTW-TLL-032"       "MTW-TLL-034"       "MTW-TLL-035"      
# [22] "MTW-TLL-040"       "SMB-ANK-009"       "SMB-ANK-011"      
# [25] "SMB-ANK-013"       "SMB-ANK-016Batch2" "SMB-ANK-027Batch1"
# [28] "SMB-ANK-027Batch2" "SMB-ANK-028"       "SMB-PDT001"       
# [31] "SMB-PTB028"        "SMB-RIN009"        "SMB-RIN016"       
# [34] "SMB-WNG-001"       "SMB-WNG-003"       "SMB-WNG-004"      
# [37] "SMB-WNG-014"       "SMB-WNG-018"       "SMB-WNG-023" 

# Strangely, one replicate samples, SMB-ANK-016 (in batch two) has data, but its replicate does not
# one other replicate, SMB-ANK-027, has two samples for batch one and batch two

# barplot of library sizes
pdf(paste0(outputdir,"BarplotSeqDepth_counts_250KIndonesian.pdf"), width=15)
ggbarplot(meta(subSampled_250K_Indo_Counts_physeq), "SampleName", "SeqDepth", fill = "SamplePop") + rotate_x_text()
dev.off()

# Do this for RPM
SeqDepth_RPM = colSums(otu_table(subSampled_250K_Indo_RPM_physeq))
sample_data(subSampled_250K_Indo_RPM_physeq)$SeqDepth = SeqDepth_RPM

pdf(paste0(outputdir,"BarplotSeqDepth_RPM_250KIndonesian.pdf"), width=15)
ggbarplot(meta(subSampled_250K_Indo_RPM_physeq), "SampleName", "SeqDepth", fill = "SamplePop") + rotate_x_text()
dev.off()

# save OTU table
scaledOTUs = sapply(1:123, function(x) otu_table(subSampled_250K_Indo_RPM_physeq)[,x]/sum(otu_table(subSampled_250K_Indo_RPM_physeq)[,x]))
colnames(scaledOTUs) = colnames(otu_table(subSampled_250K_Indo_RPM_physeq))
rownames(scaledOTUs) = rownames(otu_table(subSampled_250K_Indo_RPM_physeq))
write.table(scaledOTUs, file=paste0(outputdir,"scaledOTUs.txt"))
# save unscaled OTU table
write.table(otu_table(subSampled_250K_Indo_RPM_physeq), file=paste0(outputdir,"OTUtable.txt"))


##############################
### Pathogen Investigation ###
##############################


for (rank in c("Phylum","Class","Order","Family")) {
	p=plot_bar(subSampled_250K_Indo_RPM_physeq, fill = rank)
	fig <- p + geom_bar(stat="identity", position="stack")
	pdf(paste0(outputdir,"AllReads_",rank,".pdf"), width=15)
	fig
	dev.off()
}

# add a new column containing family names and superkingdom
tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"] = paste(tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"], tax_table(subSampled_250K_Indo_RPM_physeq)[,"Family"], sep="_")
tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"] <- gsub("Bacteria_$", "Bacteria_unclassified", tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"])
tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"] <- gsub("Eukaryota_$", "Eukaryota_unclassified", tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"])
tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"] <- gsub("Viruses_$", "Viruses_unclassified", tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"])

# Plot
p=plot_bar(subSampled_250K_Indo_RPM_physeq, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
 # Bacteria Eukaryota   unk		Viruses 
 #    32         12      1  		3 

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(32)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(12)
PaletteUnk = colorRampPalette(c("burlywood"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(3)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteUnk,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Superkingdom), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Barplot_Filtered_RPM.pdf"), width=15)
fig
dev.off()

# without Plasmodium
plasmodium_rows=grep("Plasmodiidae",tax_table(subSampled_250K_Indo_RPM_physeq)[,"Superkingdom"])
plasmodium_otu=rownames(tax_table(subSampled_250K_Indo_RPM_physeq)[plasmodium_rows,])
allTaxa = taxa_names(subSampled_250K_Indo_RPM_physeq)
allTaxa <- allTaxa[!(allTaxa %in% plasmodium_otu)]
subSampled_250K_Indo_RPM_physeq_noPlasmo=prune_taxa(allTaxa, subSampled_250K_Indo_RPM_physeq)

# Plot
p=plot_bar(subSampled_250K_Indo_RPM_physeq_noPlasmo, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))
 # Bacteria Eukaryota   unk		Viruses 
 #    32         11      1  		3 

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(32)
PaletteEukaryote = colorRampPalette(c("#fd8d3c","#800026"))(11)
PaletteUnk = colorRampPalette(c("burlywood"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(3)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteUnk,PaletteVirus)

fig <- p + scale_fill_manual(values=Merged_Palette) +
  geom_bar(aes(fill=Superkingdom), stat="identity", position="stack") +
  guides(fill=guide_legend(ncol=2))

pdf(paste0(outputdir,"Barplot_Filtered_RPM_NoPlasmo.pdf"), width=15)
fig
dev.off()

# PCA
mergedlog <- transform_sample_counts(subSampled_250K_Indo_Counts_physeq, function(x) log(1 + x))
out.wuf.log <- ordinate(mergedlog, method = "PCoA", distance = "euclidean")
pdf(paste0(outputdir,"PCA_Plasmo.pdf"), width=10)
plot_ordination(mergedlog, out.wuf.log, color="SamplePop", axes = 1:2, label="SampleName")
dev.off()

pdf(paste0(outputdir,"PCA_Plasmo.pdf"), width=10)
plot_ordination(mergedlog, out.wuf.log, color="SamplePop", axes = 1:2) +
geom_text(aes(label = sample_data(subSampled_250K_Indo_Counts_physeq)$SampleName), size = 5, vjust = 1.5) + 
geom_point(size = 4) +
theme_bw()
dev.off()

