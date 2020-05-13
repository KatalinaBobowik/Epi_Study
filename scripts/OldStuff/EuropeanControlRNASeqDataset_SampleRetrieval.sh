# Created 14.02.19
# Find out what unmapped sequences from STAR bam files are mapping to
 
# first, laod the SRA-Toolkit module
module load SRA-Toolkit
module load web_proxy

# set outdir
outdir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data

# The first six samples we're interested in are from the following study: "Variation in RNA-Seq Transcriptome Profiles of Peripheral Whole Blood from Healthy Individuals with and without Globin Depletion".
Globindepletion_paper="SRR1060753 SRR1060754 SRR1060755 SRR1060756 SRR1060757 SRR1060758"

# we can also get samples from Irene's paper: https://www.ncbi.nlm.nih.gov/pubmed/24885439
RIN_paper="SRR1300821 SRR1300811 SRR1300801 SRR1300791"

# combine both
healthy_samples="$Globindepletion_paper $RIN_paper"


# command for one sample:
for i in $healthy_samples; do
	echo fastq-dump --gzip --outdir $outdir $i
done > unmappedReads_Controls_Array.txt
