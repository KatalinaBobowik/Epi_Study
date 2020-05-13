# 02/03/2020

# Run KMA
rootDir=/data/cephfs/punim0586/kbobowik
input_dir=/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass
output_dir=/data/cephfs/punim0586/kbobowik/CCmetagen/IndoSample
array_dir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files
mkdir $output_dir
nt_db=/scratch/punim0586/kat/ncbi_nt_no_env_11jun2019

# Run KMA on all
for r1 in $input_dir/{sample_output,second_batch/sample_output,third_batch/sample_output}/R1_unmapped_*.fastq; do
	sample=`basename $r1`
	sampleID_part1=$(basename $r1 Aligned.sortedByCoord.out.fastq)
	sampleID=`echo ${sampleID_part1##*_}`
	dir=`dirname $r1`
	r2_filename=`echo "$sample" | sed -r 's/R1/R2/g'`
	r2=`echo ${dir}/${r2_filename}`
	echo ${rootDir}/CCmetagen/kma/kma -ipe $r1 $r2 -o ${output_dir}/${sampleID} -ef -t_db $nt_db -t 4 -1t1 -mem_mode -and -apm f
done > ${array_dir}/KMA_AllIndoFiles.txt

# saved as sbatch ${rootDir}/EpidemiologicalSurvey/scripts/shell_scripts/KMA_AllIndoSamples.sh 
