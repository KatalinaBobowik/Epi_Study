#!/usr/bin/env bash

# subSample data, taken from: Tutorial: BASH one-liners for subsampling reads
outdir=/data/scratch/projects/punim0586/kat/101BPIndo_SubSampled
for sample in /data/scratch/projects/punim0586/kat/101BPIndonesianSamplesMapped/*_Unmapped.out.mate1; do
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

inputdir=/data/scratch/projects/punim0586/kat/101BPIndo_SubSampled
outdir=/data/scratch/projects/punim0586/kat/101BPIndo_SubSampled_EndingRemoved
for sample in ${inputdir}/*_1.fastq; do
	sampleID=`basename $sample _1.fastq`
	sed 's/\s.*$//' ${inputdir}/${sampleID}_1.fastq > ${outdir}/${sampleID}_1.fastq
	sed 's/\s.*$//' ${inputdir}/${sampleID}_2.fastq > ${outdir}/${sampleID}_2.fastq
done

# for control data
# subSample data, taken from: Tutorial: BASH one-liners for subsampling reads
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

# for 50BP trimmed data
outdir=/data/scratch/projects/punim0586/kat/50BPIndo_SubSampled
for sample in /data/scratch/projects/punim0586/kat/IndonesianSamplesMapped/*_Unmapped.out.mate1; do
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

inputdir=/data/scratch/projects/punim0586/kat/50BPIndo_SubSampled
outdir=/data/scratch/projects/punim0586/kat/50BPIndo_SubSampled_EndingRemoved
for sample in ${inputdir}/*_1.fastq; do
	sampleID=`basename $sample _1.fastq`
	sed 's/\s.*$//' ${inputdir}/${sampleID}_1.fastq > ${outdir}/${sampleID}_1.fastq
	sed 's/\s.*$//' ${inputdir}/${sampleID}_2.fastq > ${outdir}/${sampleID}_2.fastq
done



paste /data/scratch/projects/punim0586/kat/101BPIndonesianSamplesMapped/SMB-PTB028_Unmapped.out.mate1 /data/scratch/projects/punim0586/kat/101BPIndonesianSamplesMapped/SMB-PTB028_Unmapped.out.mate2 | awk '{ printf("%s",$0); n++;
if(n%4==0) { printf("\n");} else { printf("\t");} }' |
awk -v k=100000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |
awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 > "/data/scratch/projects/punim0586/kat/101BPIndo_SubSampled/SMB-PTB028_1.fastq";print
$2"\n"$4"\n"$6"\n"$8 > "/data/scratch/projects/punim0586/kat/101BPIndo_SubSampled/SMB-PTB028_2.fastq"}'

# Run kneaddata
kneaddata --input /data/scratch/projects/punim0586/kat/101BPIndo_SubSampled/MTW-TLL-019_1.fastq \
--input /data/scratch/projects/punim0586/kat/101BPIndo_SubSampled/MTW-TLL-019_2.fastq \
--output kneaddata_output5 --bypass-trim --run-trf -db SILVA_128_LSUParc_SSUParc_ribosomal_RNA

# run KMA
/data/cephfs/punim0586/kbobowik/CCmetagen/kma/kma \
-ipe /data/scratch/projects/punim0586/kat/kneaddataGenome/Homo_sapiens_BMTagger_v0.1/kneaddata_output13/SRR8367969_Unmapped.out_kneaddata.repeats.removed.1.fastq /data/scratch/projects/punim0586/kat/kneaddataGenome/Homo_sapiens_BMTagger_v0.1/kneaddata_output13/SRR8367969_Unmapped.out_kneaddata.repeats.removed.2.fastq \
-o testKMA3 -ef -t_db \
/scratch/punim0586/kat/ncbi_nt_no_env_11jun2019 -t 12 -1t1 -mem_mode -apm f -and

# run CCMetagen
source activate py38
python /data/cephfs/punim0586/kbobowik/CCmetagen/CCMetagen/CCMetagen.py \
-i testKMA3.res \
-o SRR8367969 \
--depth_unit rpm --mapstat testKMA3.mapstat \
--depth 1

# merge data
CCMetagen_merge.py --input_fp /data/scratch/projects/punim0586/kat/50BPControlSamplesCCMetagen --keep_or_remove r --filtering_tax_level Phylum --taxa_list Chordata,Arthropoda --output_fp ControlUnmapped_species_table_noRepeats_RPM


