# Created 28.02.20
# Run CCMetagen aligner on unmapped reads of Indonesian samples

# ask for a lot of memory in sinteractive 
sinteractive --ntasks=1 --cpus-per-task=32 --mem=500000 --time=04:00:00 -p mig

# Run conda session

# conda create -n py37 python=3.7
# conda install -c anaconda pandas
# conda install -c bioconda java-jdk
# conda install -c anaconda numpy
# conda install -c etetoolkit ete3
# conda install -c bioconda krona

source activate py37

# Download and install KMA
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make

# download CCMetagen and add it to your path
git clone https://github.com/vrmarcelino/CCMetagen

# add the CCMetagen python scripts to path
vi ~/.bashrc
export PATH="$PATH:/data/cephfs/punim0586/kbobowik/CCmetagen/CCMetagen"
source ~/.bashrc

# get sequences for database without environmental or artificia sequences
wget "https://cloudstor.aarnet.edu.au/plus/s/vfKH9S8c5FVGBjV/download?path=%2F&files=ncbi_nt_no_env_11jun2019.zip&downloadStartSecret=m36d6mrye9s"

# Map with KMA 

# For paired-end files:
# where -t_db is the reference database 
# -t is the flag for the number of threads 
# -ipe is the flag for mate1 and mat2 of a paired-end metagenome/metatranscriptome sample (fastq or fasta) 
# $SAMPLE is the path to a single-end metagenome/metatranscriptome file (reads or contigs)

rootDir="/data/cephfs/punim0586/kbobowik"

input_dir=/data/cephfs/punim0586/kbobowik/CCmetagen/CCMetagen/tutorial/figs_tutorial/Subset_used
output_dir=01_KMA_res
mkdir $output_dir
nt_db=/scratch/punim0586/kat/ncbi_nt_no_env_11jun2019

for r1 in $input_dir/*R1.fastq; do
	r2=${r1/R1.fastq/R2.fastq}
	o_part1=$output_dir/${r1/$input_dir\//''}
	o=${o_part1/.R*/}
	echo $o
	/usr/bin/time -v ${rootDir}/CCmetagen/kma/kma -ipe $r1 $r2 -o /scratch/punim0586/kat/kat -t_db $nt_db -t 4 -1t1 -mem_mode -and -apm f
done

# output of usr/bin/time
	# User time (seconds): 122.24
	# System time (seconds): 74.14
	# Percent of CPU this job got: 101%
	# Elapsed (wall clock) time (h:mm:ss or m:ss): 3:12.55
	# Average shared text size (kbytes): 0
	# Average unshared data size (kbytes): 0
	# Average stack size (kbytes): 0
	# Average total size (kbytes): 0
	# Maximum resident set size (kbytes): 135654768
	# Average resident set size (kbytes): 0
	# Major (requiring I/O) page faults: 0
	# Minor (reclaiming a frame) page faults: 90395
	# Voluntary context switches: 2344188
	# Involuntary context switches: 3150
	# Swaps: 0
	# File system inputs: 0
	# File system outputs: 17936
	# Socket messages sent: 0
	# Socket messages received: 0
	# Signals delivered: 0
	# Page size (bytes): 4096
	# Exit status: 0

# run CCMetagen
CCMetagen.py -i /scratch/punim0586/kat.res -o ~/katresults












