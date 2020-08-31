#!/usr/bin/env bash

# subSample data to 250,000 reads
# subsampling code taken from: Tutorial: BASH one-liners for subsampling reads
outdir=/data/scratch/projects/punim0586/kat/50BPControl_SubSampled
for sample in /data/scratch/projects/punim0586/kat/Bipolar_controlSamples_EpiStudy_MappedTwoPass/*_Unmapped.out.mate1; do
	sampleID=`basename $sample _Unmapped.out.mate1`
	path=`dirname $sample`
	f1_output=${outdir}/${sampleID}_1.fastq
	f2_output=${outdir}/${sampleID}_2.fastq
	paste ${path}/${sampleID}_Unmapped.out.mate1 ${path}/${sampleID}_Unmapped.out.mate2 | \
    awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' | \
    awk -v k=250000 \
    'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | \
    awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 > "'$f1_output'";print$2"\n"$4"\n"$6"\n"$8 > "'$f2_output'"}'
done

# the ending '0:N:  00' in the name line of the fastq file (the first of the four lines) needs to be removed
# otherwise, KneadData will throw an error
inputdir=/data/scratch/projects/punim0586/kat/50BPIndo_SubSampled
outdir=/data/scratch/projects/punim0586/kat/50BPIndo_SubSampled_EndingRemoved
for sample in ${inputdir}/*_1.fastq; do
	sampleID=`basename $sample _1.fastq`
	sed 's/\s.*$//' ${inputdir}/${sampleID}_1.fastq > ${outdir}/${sampleID}_1.fastq
	sed 's/\s.*$//' ${inputdir}/${sampleID}_2.fastq > ${outdir}/${sampleID}_2.fastq
done

