# Script to tidy up the output I made with the Snakemake workflow 'Snakefile_ReadSummary'

# load library
library(tidyverse) 

# set paths
inputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/QC/"
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/QC/"

# 75 BP --------------------------

# Indonesian data
Indonesian_Reads_Summary <- read.table(paste0(inputdir,"Indonesian_ReadDataSummary_75BP.txt"))
grouped_samplenames_Indo <- Indonesian_Reads_Summary %>% group_by(V1)
Indo_Summary <- grouped_samplenames_Indo %>% tidyr::spread(V2,V3)
write.table(as.data.frame(Indo_Summary), file=paste0(outputdir,"Indonesian_Reads_Summary_75BP.txt"), quote = F, sep = "\t", row.names = F)

# Malian data
Malian_Reads_Summary <- read.table(paste0(inputdir,"MalianReadsSummary_75BP.txt"))
grouped_samplenames_Mali <- Malian_Reads_Summary %>% group_by(V1)
Mali_Summary <- grouped_samplenames_Mali %>% tidyr::spread(V2,V3)
write.table(as.data.frame(Mali_Summary), file=paste0(outputdir,"Malian_Reads_Summary.txt"), quote = F, sep = "\t", row.names = F)

# UK data
UK_Reads_Summary <- read.table(paste0(inputdir,"TB_ReadDataSummary_75BP.txt"))
grouped_samplenames_UK <- UK_Reads_Summary %>% group_by(V1)
UK_Summary <- grouped_samplenames_UK %>% tidyr::spread(V2,V3)
write.table(as.data.frame(UK_Summary), file=paste0(outputdir,"UK_Reads_Summary.txt"), quote = F, sep = "\t", row.names = F)

# 101 BP --------------------------

# Indonesian data
Indonesian_Reads_Summary <- read.table(paste0(inputdir,"Indonesian_ReadDataSummary_101BP.txt"))
grouped_samplenames_Indo <- Indonesian_Reads_Summary %>% group_by(V1)
Indo_Summary <- grouped_samplenames_Indo %>% tidyr::spread(V2,V3)
write.table(as.data.frame(Indo_Summary), file=paste0(outputdir,"Indonesian_Reads_Summary_101BP.txt"), quote = F, sep = "\t", row.names = F)

# alternative, non tidyr version
# for (name in as.character(unique(Indonesian_Reads_Summary[,1]))){
# 	print(Indonesian_Reads_Summary[which(Indonesian_Reads_Summary[,1]==name),3])
# }
