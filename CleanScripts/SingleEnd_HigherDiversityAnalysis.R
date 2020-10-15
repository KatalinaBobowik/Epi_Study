# Investigate why single ended data is more taxonomically rich than paired end data
# Note- this is on subsampled reads

library(ggplot2)
library(ggpubr)
library(phyloseq)
library(reshape2)

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"
batchInfodir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"
# load in data
forward_strand <- read.csv(paste0(inputdir,"Counts_Indo_AllReads_SE.csv"),check.names=FALSE)
reverse_strand <- read.csv(paste0(inputdir,"Counts_Indo_AllReads_SE_ReverseRead.csv"),check.names=FALSE)

#############################
# Analysis of dissimilarity #
#############################

# The first thing we want to do is test the dissimilarity between forward and reverse reads.
# We'll do this by making a dataframe containing the number of similar taxa for each taxonomic rank

# make a taxonomic rank variable
taxonomicRanks=c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Superkingdom")
samplenames=colnames(reverse_strand)[1:123]

# make 'not in' own function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Now creat the df
dissimilarityDF=data.frame(Samples=as.character(),taxonomicRank=as.character(),nDissimilar=as.numeric(),percDissimilar=as.numeric(),stringsAsFactors=FALSE)
counter=0
for (rank in taxonomicRanks){
  for (sample in samplenames){
    counter=counter+1
    a=forward_strand[,c(sample,rank)]
    a=a[which(a[,1]>0),]
    b=reverse_strand[,c(sample,rank)]
    b=b[which(b[,1]>0),]
    nDisagree=length(which(a[,rank]%!in%b[,rank]))
    percDissimilar=nDisagree/nrow(a)
    dissimilarityDF[counter,]=c(sample,rank,nDisagree,percDissimilar)
  }
}

# convert last two columns to numeric
dissimilarityDF$nDissimilar=as.numeric(dissimilarityDF$nDissimilar)
dissimilarityDF$percDissimilar=as.numeric(dissimilarityDF$percDissimilar)
# order factors
dissimilarityDF$taxonomicRank <- factor(dissimilarityDF$taxonomicRank, levels = taxonomicRanks)

# Make violin plot
pdf(paste0(outputdir,"ForwardStrand_Vs_RevereseStrand_Similarity.pdf"), height = 8, width = 10)
ggboxplot(dissimilarityDF, x = "taxonomicRank", y = "percDissimilar", fill="taxonomicRank",add=c("boxplot"))
dev.off()

# Using phyloseq

# Forward strand
taxa_raw <- as.matrix(forward_strand[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(forward_strand[,-which(colnames(forward_strand) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
forward = phyloseq(taxa, tax)
# add in population info
pop <- rep("forward",ncol(otu_table(forward)))
# get batch info
load(paste0(batchInfodir, "indoRNA.read_counts.TMM.filtered.Rda"))
batch_df = data.frame(rownames(y$samples),y$samples$batch)
colnames(batch_df)=c("Sample","Batch")
batch_df$Sample=gsub("_firstBatch","",batch_df$Sample) %>% gsub("_secondBatch","",.) %>% gsub("_thirdBatch","",.)
samplenames_forward=gsub("Batch1","",colnames(otu_table(forward))) %>% gsub("Batch2","",.) %>% gsub("Batch3","",.)
batch=batch_df[match(samplenames_forward, batch_df$Sample),"Batch"]

# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(forward)), SamplePop=pop, batch=batch)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(forward) <- samples
taxa_names(forward) <- paste(tax_table(forward)[,"Superkingdom"], tax_table(forward)[,"Kingdom"], tax_table(forward)[,"Phylum"], tax_table(forward)[,"Class"], tax_table(forward)[,"Order"], tax_table(forward)[,"Family"], tax_table(forward)[,"Genus"], tax_table(forward)[,"Species"], sep="_")


# reverse strand
taxa_raw <- as.matrix(reverse_strand[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
abund_raw <- as.matrix(reverse_strand[,-which(colnames(reverse_strand) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
colnames(abund_raw)=paste(colnames(abund_raw),"rev",sep="_")
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)
reverse = phyloseq(taxa, tax)
# add in population info
pop <- rep("reverse",ncol(otu_table(reverse)))
samplenames_reverse=gsub("_rev","",colnames(otu_table(reverse))) %>% gsub("Batch1","",.) %>% gsub("Batch2","",.) %>% gsub("Batch3","",.)
batch=batch_df[match(samplenames_reverse, batch_df$Sample),"Batch"]

# make this into a df and add to the Phloseq object
samples_df=data.frame(SampleName=colnames(otu_table(reverse)), SamplePop=pop, batch=batch)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(reverse) <- samples
taxa_names(reverse) <- paste(tax_table(reverse)[,"Superkingdom"], tax_table(reverse)[,"Kingdom"], tax_table(reverse)[,"Phylum"], tax_table(reverse)[,"Class"], tax_table(reverse)[,"Order"], tax_table(reverse)[,"Family"], tax_table(reverse)[,"Genus"], tax_table(reverse)[,"Species"], sep="_")

# merge
merged_phylo_counts=merge_phyloseq(forward, reverse)

# Make correlation matrix for each taxonomic rank
taxa_summary=NULL
taxaNames=c("Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species")
for (taxa in taxaNames){
  print(taxa)
  Species=tax_glom(merged_phylo_counts,taxa)
  a=as.data.frame(otu_table(Species))
  correlation=cor(a,method="spearman")
  same=diag(correlation[c(1:123),c(124:ncol(correlation))])
  diff=correlation[c(1:123),c(124:ncol(correlation))]
  diag(diff)=NA
  b=melt(diff)
  b$Variability=NA
  b[which(is.na(b$value)),"Variability"]="within"
  b[which(!is.na(b$value)),"Variability"]="between"
  b[which(is.na(b$value)),"value"]=same
  b$taxa=taxa
  taxa_summary=rbind(taxa_summary, b)
}

taxa_summary$taxa <- factor(taxa_summary$taxa, levels = taxaNames)
pdf(paste0(outputdir,"ForwardStrand_Vs_RevereseStrand_Correlation.pdf"), height = 8, width = 10)
ggplot(taxa_summary, aes(x=taxa, y=value, fill=Variability)) +
  geom_boxplot() + theme_bw() + ylab("Spearman pairwise correlation") +
  theme(axis.title.x=element_blank()) + scale_fill_brewer(palette="BrBg")
dev.off()

# Look at beta diversity plots

# First, regular PCA plot
ps4.rel <- microbiome::transform(merged_phylo_counts, "clr")
# RDA without constraints is PCA        
bx.ord_pcoa_bray <- ordinate(ps4.rel, "RDA")

# Make an ordination plot using bray's dissimilarity
beta.ps1 <- plot_ordination(ps4.rel, 
                            bx.ord_pcoa_bray, 
                            color="batch", 
                            label = "SampleName") + 
  geom_point(aes(), size= 4) + 
  theme(plot.title = element_text(hjust = 0, size = 12))

# add in an ellipse
pdf(paste0(outputdir,"ForwardStrand_Vs_BetaDiversity.pdf"), height = 8, width = 10)
beta.ps1 + stat_ellipse() + theme_bw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# Or by bray distance dissimilarity
ps4.rel <- microbiome::transform(merged_phylo_counts, "compositional")
bx.ord_pcoa_bray <- ordinate(ps4.rel, "PCoA", "bray")
# Make an ordination plot using bray's dissimilarity
beta.ps1 <- plot_ordination(ps4.rel, 
                            bx.ord_pcoa_bray, 
                            color="batch", 
                            label = "SampleName") + 
  geom_point(aes(), size= 4) + 
  theme(plot.title = element_text(hjust = 0, size = 12))

# add in an ellipse
pdf(paste0(outputdir,"ForwardStrand_Vs_BetaDiversity.pdf"), height = 8, width = 10)
beta.ps1 + stat_ellipse() + theme_bw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
