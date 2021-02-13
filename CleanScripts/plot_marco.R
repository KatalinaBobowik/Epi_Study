# script created by Heini Natri and files given by Heini Natri
# file paths edited by KSB
# 12.02.2021

# load libraries
library(ggplot2)

# Set up colour schemes
KorowaiCol="#F21A00"
MentawaiCol="#3B9AB2"
SumbaCol="#EBCC2A"

# set paths
refdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/SEIndonesianSamples/"

# Plotting ENSG00000019169, MARCO, rs13425622
gt <- as.data.frame(t(read.table(paste0(refdir,"MARCO_gt.tsv"))))
colnames(gt) <- c("sample", "gt")
exp <- as.data.frame(t(read.table(paste0(refdir,"MARCO_exp.tsv"))))
colnames(exp) <- c("sample", "exp")
head(exp)
exp$exp <- as.numeric(as.character(exp$exp))

gt_exp <- merge(gt, exp, by="sample")

boxplot <- ggplot(gt_exp, aes(x=gt, y=exp))+
  geom_boxplot()

# Boxplot with jitter
y_title <- expression(paste("Normalized ", italic("MARCO"), " expression"))

# set up colour scheme for each population
pop <- sapply(strsplit(as.character(gt_exp$sample), "[-.]"), `[`, 1)
pop = gsub("MPI","KOR", pop)
gt_exp$Population = pop

boxplot <- ggplot(gt_exp, aes(x=gt, y=exp)) +
  geom_boxplot(outlier.shape = NA, color="black") + 
  geom_jitter(shape=19, size=3, width = 0.25, height = 0.25, aes(colour=Population), alpha = 0.7) + scale_color_manual(values=c(KorowaiCol,MentawaiCol,SumbaCol)) +
  #scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
  theme_bw() +
  scale_x_discrete(name="# of minor alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
  scale_y_continuous(name=y_title) +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "top")

boxplot

pdf(paste0(outputdir,"MARCO_eQTL_effect.pdf"))
boxplot
dev.off()

