#!/usr/bin/env bash

# the ending '0:N:  00' in the name line of the fastq file (the first of the four lines) needs to be removed
# otherwise, KneadData will throw an error
inputdir=/data/scratch/projects/punim0586/kat/101BPIndonesianSamplesMapped
outdir=/data/scratch/projects/punim0586/kat/101BPIndo_EndingRemoved
for sample in ${inputdir}/*_Unmapped.out.mate1; do
	sampleID=`basename $sample _Unmapped.out.mate1`
	sed 's/\s.*$//' ${inputdir}/${sampleID}_Unmapped.out.mate1 > ${outdir}/${sampleID}_1.fastq
	sed 's/\s.*$//' ${inputdir}/${sampleID}_Unmapped.out.mate2 > ${outdir}/${sampleID}_2.fastq
done
