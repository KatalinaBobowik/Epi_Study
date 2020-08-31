# Script for creating barplot of number of reads at all different processing stages
# Code developed by Katalina Bobowik, 10.12.2018

# load library and colour palette
library(viridis)
library(Rcmdr)
library(reshape2)
library(ggplot2)

# set working directory
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/QC/STAR/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/QC/"

# set ggplot colour theme to white
theme_set(theme_bw())

# read in summary file and tidy up 
a=read.table(paste0(inputdir,"mqc_star_alignment_plot_1.txt"), header=T, sep="\t")
# remove data from first pass of star
a=a[-grep("STARpass1", a$Sample),]
  
pdf(paste0(outputdir,"QCofReads_STAR.pdf"), width=12)
ggplot(data, aes(fill=forcats::fct_rev(variable), y=value, x=Sample)) + 
    geom_bar(position="fill", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank()) + 
    scale_fill_manual(values=c("#7F0000", "#B1084C", "#F7A35C", "#7CB5EC", "#437BB1")) +
    ylab("Percentage of Reads")
dev.off()