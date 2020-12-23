# Execute snakemake pipeline to get read summary of all datasets

# execute snakemake command
snakemake -s Snakefile_ReadSummary -j 30 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} -c {threads} -p {cluster.partition}" --use-conda

# combine all files and write output to a file
cat /scratch/punim0586/kat/QC/ReadSummary/*.txt > Indonesian_ReadDataSummary_101BP.txt

# move over to local computer
scp kbobowik@spartan.hpc.unimelb.edu.au:/scratch/punim0586/kat/QC/Indonesian_ReadDataSummary_75BP.txt /Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Epi_Study_SuppTables

# execute 'ConvertReadsSummaryToTidyMatrix.R'