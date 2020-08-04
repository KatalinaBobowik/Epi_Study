# 02/03/2020

# Run conda session (if running sinteractive session)
source activate py37

# download the RNA database
kneaddata_database --download ribosomal_RNA bowtie2 $DIR

INPUT_DIR=/scratch/punim0586/kat/101BPIndonesianSamplesMapped
DATABASE=/data/cephfs/punim0586/kbobowik/genome/RNA
# $DATABASE = the index of the KneadData database
# $OUTPUT_DIR = the output directory

kneaddata --input paired_trimmedOutput_SMB-WNG-023_1.fastq --input paired_trimmedOutput_SMB-WNG-023_2.fastq \
-db demo_db --output kneaddata_output  --fastqc /usr/local/easybuild/software/fastqc/0.11.8/


--trimmomatic \
 --run-trf --run-fastqc-start --run-fastqc-end --trimmomatic-options="MINLEN:90" --cut-adapters 

 ##

# note I had to rename my fasta file to have a .fastq extension.
# I also had to rename the trf074b.inx file to just be called 'trf'
 kneaddata --input forward_sub.fastq \
 --input /scratch/punim0586/kat/raw_data_copied/SMB-WNG-028_1.fastq.gz \
 -o /data/cephfs/punim0586/kbobowik/mySample \
 --fastqc /usr/local/easybuild/software/fastqc/0.11.8/
 
 --trimmomatic /data/cephfs/punim0586/kbobowik/bin/Trimmomatic-0.36/ \
 --trf /data/cephfs/punim0586/kbobowik/bin/TRF/ \
 --fastqc /usr/local/easybuild/software/fastqc/0.11.8/




 /scratch/punim0586/kat/raw_data_copied/SMB-WNG-028_1.fastq.gz