#!/bin/bash
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END
# Use this email address:
#SBATCH --mail-user=katalina.bobowik@gmail.com
# Export path variables.
#SBATCH --export=ALL
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=2:00:00
# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi
# Run the job from the directory where it was launched (default)

# Run the simulations:
source activate py37

rootDir="/data/cephfs/punim0586/kbobowik"
input_dir=/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/sample_output
output_dir=/data/cephfs/punim0586/kbobowik/CCmetagen/IndoSample

mkdir $output_dir
nt_db=/scratch/punim0586/kat/ncbi_nt_no_env_11jun2019
r1=R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_MPI-048Aligned.sortedByCoord.out.fastq
r2=R2_unmapped_STAR_Hg38_second_paired_trimmedOutput_MPI-048Aligned.sortedByCoord.out.fastq

${rootDir}/CCmetagen/kma/kma -ipe $r1 $r2 -o $output_dir -t_db $nt_db -t 4 -1t1 -mem_mode -and -apm f

# run CCMetagen
CCMetagen.py -i /data/cephfs/punim0586/kbobowik/CCmetagen/IndoSample.res -o ~/CCMetagenOutput/katresults
