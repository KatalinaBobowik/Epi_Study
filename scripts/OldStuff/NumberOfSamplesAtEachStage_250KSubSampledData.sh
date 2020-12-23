#!/usr/bin/env bash

# Indonesian data --------------
QCDir=/data/scratch/projects/punim0586/kat/QC
arrayDir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files

# raw data 
MaliDir=/scratch/punim0586/kat/ControlSamples_Mali
TBDir=/scratch/punim0586/kat/ControlSamples_TB
for file in `find $MaliDir $TBDir -name '*_1.fastq.gz'`; do
	sample=`basename $file _1.fastq.gz`
	echo "echo $sample \`expr \$(zcat $file|wc -l)/4|bc\` >> ${QCDir}/totalRawReads.txt"

done > ${arrayDir}/rawReadsCounter_ControlsMaliUK.txt


# Trimmed reads with adapters removed
trimDir=/data/scratch/projects/punim0586/kat/101BPIndonesianSamplesTrimmed

for file in ${trimDir}/*_1.fastq.gz; do
	sample=`basename $file _1.fastq.gz`
	echo "echo $sample \`expr \$(zcat $file|wc -l)/4|bc\` >> ${QCDir}/trimmedNoAdapterIndo.txt"
done > ${arrayDir}/trimmedReadsCounter.txt

# mapped and unmapped reads
unmappedDir=/data/scratch/projects/punim0586/kat/101BPIndonesianSamplesMapped
for file in ${unmappedDir}/*_Unmapped.out.mate1; do
	sample=`basename $file _Unmapped.out.mate1`
	echo $sample `expr $(cat $file|wc -l)/4|bc` >> ${QCDir}/unmappedReadTotal.txt
done

# human contaminant reads
kneaddataDir=/data/scratch/projects/punim0586/kat/101BPIndoKneaddata
for file in ${kneaddataDir}/*kneaddata_paired_Homo_sapiens_assembly19.fasta_bmtagger_contam_1.fastq; do
	sample=`basename $file _1_kneaddata_paired_Homo_sapiens_assembly19.fasta_bmtagger_contam_1.fastq`
	echo $sample `expr $(cat $file|wc -l)/4|bc` >> ${QCDir}/noHumanReads.txt
done

# repeats removed
for file in ${kneaddataDir}/*kneaddata.repeats.removed.1.fastq; do
	sample=`basename $file _1_kneaddata.repeats.removed.1.fastq`
	echo $sample `expr $(cat $file|wc -l)/4|bc` >> ${QCDir}/noRepeats.txt
done

# sort everything 
sort -k1 ${QCDir}/totalRawReads.txt > ${QCDir}/totalRawReads_Sorted.txt
sort -k1 ${QCDir}/trimmedNoAdapterIndo.txt > ${QCDir}/trimmedNoAdapterIndo_Sorted.txt
sort -k1 ${QCDir}/unmappedReadTotal.txt > ${QCDir}/unmappedReadTotal_Sorted.txt
sort -k1 ${QCDir}/noHumanReads.txt > ${QCDir}/noHumanReads_Sorted.txt
sort -k1 ${QCDir}/noRepeats.txt > ${QCDir}/noRepeats_Sorted.txt

paste ${QCDir}/totalRawReads_Sorted.txt ${QCDir}/trimmedNoAdapterIndo_Sorted.txt ${QCDir}/unmappedReadTotal_Sorted.txt ${QCDir}/noHumanReads_Sorted.txt ${QCDir}/noRepeats_Sorted.txt > ${QCDir}/allReads_AllStages.txt
cat allReads_AllStages.txt | awk -F" " '{print $1, $2, $4, $6, $8, $10}' > allReads_AllStages_Ordered.txt

cat * | sort -k1
