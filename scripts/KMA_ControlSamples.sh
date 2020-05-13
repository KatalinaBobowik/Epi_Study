# 25/03/2020

# Run KMA
rootDir=/data/cephfs/punim0586/kbobowik
input_dir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/trimmedControls
output_dir=/data/cephfs/punim0586/kbobowik/CCmetagen/ControlSamples
array_dir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files
nt_db=/scratch/punim0586/kat/ncbi_nt_no_env_11jun2019

# Run KMA on all
for r1 in $input_dir/paired_trimmedOutput_R1_*.fastq; do
	sample=`basename $r1`
	sampleID_part1=$(basename $r1 _100K.fastq)
	sampleID=`echo ${sampleID_part1#"paired_trimmedOutput_R1_"}`
	dir=`dirname $r1`
	r2_filename=`echo "$sample" | sed -r 's/R1/R2/g'`
	r2=`echo ${dir}/${r2_filename}`
	# flag descriptions: -ipe Inputfile(s), paired end, -o output; -ef used to calculate abundance in reads per million (RPM); -t_db database; -t threads
	# -1t1 One read to one template, no splicing performed; -mem_mode *.index and *.seq are not loaded into memory, which enables one to map against larger databases. Templates are chosen using k-mer counting
	# -apm Paired end method, “p”: Reward if pairing the reads, “u”: unite best template matches in each read if possible, “f” force paired reads to pair.
	# Both mrs and p_value thresholds has to reached to in order to report a template hit.
	echo ${rootDir}/CCmetagen/kma/kma -ipe $r1 $r2 -o ${output_dir}/${sampleID} -ef -t_db $nt_db -t 4 -1t1 -mem_mode -and -apm f
done > ${array_dir}/KMA_AllControlFiles.txt

# saved as sbatch ${rootDir}/EpidemiologicalSurvey/scripts/shell_scripts/KMA_AllIndoSamples.sh 


---------------------------------------------
## Test KMA on control samples without thresholds
# Run KMA on all
for r1 in $input_dir/paired_trimmedOutput_R1_*.fastq; do
	sample=`basename $r1`
	sampleID_part1=$(basename $r1 _100K.fastq)
	sampleID=`echo ${sampleID_part1#"paired_trimmedOutput_R1_"}`
	dir=`dirname $r1`
	r2_filename=`echo "$sample" | sed -r 's/R1/R2/g'`
	r2=`echo ${dir}/${r2_filename}`
	# flag descriptions: -ipe Inputfile(s), paired end, -o output; -ef used to calculate abundance in reads per million (RPM); -t_db database; -t threads
	# -1t1 One read to one template, no splicing performed; -mem_mode *.index and *.seq are not loaded into memory, which enables one to map against larger databases. Templates are chosen using k-mer counting
	# -apm Paired end method, “p”: Reward if pairing the reads, “u”: unite best template matches in each read if possible, “f” force paired reads to pair.
	# -and Both mrs and p_value thresholds has to reached to in order to report a template hit.
	echo ${rootDir}/CCmetagen/kma/kma -ipe $r1 $r2 -o ${output_dir}/${sampleID} -ef -t_db $nt_db -t 4 -1t1 -mem_mode -and -apm f
done > ${array_dir}/KMA_AllControlFiles.txt


${rootDir}/CCmetagen/kma/kma -ipe $r1 $r2 -o /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/${sampleID} -t_db $nt_db -t 4 -1t1 -mem_mode -and -apm p
${rootDir}/CCmetagen/kma/kma -ipe $r1 $r2 -o /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/${sampleID} -ef -t_db $nt_db -t 4 -1t1 -mem_mode -and -apm f

${rootDir}/CCmetagen/kma/kma -i /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/GSM2139514_unmapped_108000B_100K.fastq -o /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/test -t_db $nt_db -t 4 -1t1 -mem_mode -and


