# Created 22.09.20

# The first step is to get one sample from the following study: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP134018&o=acc_s%3Aa
# We'll just choose a random control sample for now

# download with fasterq-dump
cd /data/scratch/projects/punim0586/kat/TuberculosisSA_Controls
fasterq-dump SRR6809895

# First map the reads and get all of the unmapped reads

# subSample data to 100,000 reads
# subsampling code taken from: Tutorial: BASH one-liners for subsampling reads
sampleID=`basename SRR6809895_1.fastq _1.fastq`
f1_output=SRR6809895_100K_1.fastq
f2_output=SRR6809895_100K_1.fastq
paste SRR6809895_1.fastq SRR6809895_2.fastq | \
awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' | \
awk -v k=100000 \
'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | \
awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 > "'$f1_output'";print$2"\n"$4"\n"$6"\n"$8 > "'$f2_output'"}'

# Subsample 100K reads randomly
kneaddata --input SRR11994658.fastq --output control2 \
--trimmomatic-options="LEADING:20 TRAILING:20 MINLEN:90" --run-trf \
-db /data/scratch/projects/punim0586/kat/kneaddataGenome/Homo_sapiens_BMTagger_v0.1/Homo_sapiens_assembly19.fasta \
--run-bmtagger

# run KMA
/data/cephfs/punim0586/kbobowik/CCmetagen/kma/kma \
-i SRR11994658_kneaddata.repeats.removed.fastq -o /data/scratch/projects/punim0586/kat/100KBPControls/kma/SRR11994658 \
-ef -t_db /data/scratch/projects/punim0586/kat/ncbi_nt_no_env_11jun2019 \
-t 12 -1t1 -mem_mode -and

# CCMetagen
python /data/cephfs/punim0586/kbobowik/CCmetagen/CCMetagen/CCMetagen.py -i {input.resultFile} \
-o /data/scratch/projects/punim0586/kat/100KBPControls/CCMetagen/{wildcards.sample} \
--depth_unit rpm --mapstat {input.mapstatFile} --depth 1

