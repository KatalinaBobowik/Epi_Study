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
	# flag descriptions: -ipe Inputfile(s), paired end, -o output; -ef used to calculate abundance in reads per million (RPM); -t_db database; -t threads
	# -1t1 One read to one template, no splicing performed; -mem_mode *.index and *.seq are not loaded into memory, which enables one to map against larger databases. Templates are chosen using k-mer counting
	# -apm Paired end method, “p”: Reward if pairing the reads, “u”: unite best template matches in each read if possible, “f” force paired reads to pair.
	# Both mrs and p_value thresholds has to reached to in order to report a template hit.
done > ${array_dir}/KMA_AllIndoFiles.txt

# saved as sbatch ${rootDir}/EpidemiologicalSurvey/scripts/shell_scripts/KMA_AllIndoSamples.sh 
