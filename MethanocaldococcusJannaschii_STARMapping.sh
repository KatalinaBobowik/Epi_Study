# created by KSB, 29.01.2021
# The aim of this script is to get all Methanogen sequences from the output of KMA fasta files
# Once those are retrieved, get the Methanogen whole geneome and map the fasta files to the Methanogen genome with STAR
# Finally, convert these to BigWig format in order to view with IGV

# load modules
module load star
module load deeptools


MaliDir=/data/cephfs/punim0586/kbobowik/MalianData
genomeDir=/data/cephfs/punim0586/kbobowik/genome/MethanococcusGenome/
outDir=/data/cephfs/punim0586/kbobowik/MalianData/STAR_MappedReads
cleanedReads_outdir=/data/cephfs/punim0586/kbobowik/MalianData/STAR_Clean_MappedReads
unmappedDir=/scratch/punim0586/kat/MaliSamplesMapped
cleanDir=/scratch/punim0586/kat/MaliSamplesKneadData_75BP

### Download Methanocaldococcus jannaschii DSM 2661 genome from NCBI ----------------------

# url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/
# Download both the genome fasta file and the gtf file

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/GCA_000091665.1_ASM9166v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/GCA_000091665.1_ASM9166v1_genomic.gtf.gz 


# Generate genome indices
# genomeSAindexNbases is set to 10 for small genomes (otherwise, thrpws bus error)
STAR --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFiles /data/cephfs/punim0586/kbobowik/genome/GCA_000091665.1_ASM9166v1_genomic.fna --sjdbGTFfile /data/cephfs/punim0586/kbobowik/genome/GCA_000091665.1_ASM9166v1_genomic.gtf --runThreadN 12 --genomeSAindexNbases 10

# map all unmapped reads from Malian samples
for file in ${unmappedDir}/*.mate1; do
	f1=`basename $file`
	f2="${f1%.*}.mate2"
	ID=${f1%_*}
	echo STAR --genomeDir ${genomeDir} --readFilesIn ${unmappedDir}/${f1} ${unmappedDir}/${f2} --runThreadN 6 --outFileNamePrefix ${outDir}/${ID}_unmapped_ --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 999 --genomeSAindexNbases 10
	# bamCoverage -b $file -o ${file}.bw
done > /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files/MalianMappedMethanogens_Array.txt

# index files
for file in ${outDir}/*_unmapped_Aligned.sortedByCoord.out.bam; do
	#filename=`basename $file _MethanogensAligned.sortedByCoord.out.bam`
	samtools index $file
done

# Do this for cleaned reads with human contamination and repeats removed
# map all unmapped reads from Malian samples
for file in ${cleanDir}/*kneaddata.repeats.removed.fastq; do
	f1=`basename $file`
	ID=`basename $f1 _Unmapped.out_kneaddata.repeats.removed.fastq`
	echo STAR --genomeDir ${genomeDir} --readFilesIn $file --runThreadN 6 --outFileNamePrefix ${cleanedReads_outdir}/${ID}_unmapped_ --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 999 --genomeSAindexNbases 10
	# bamCoverage -b $file -o ${file}.bw
done > /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files/MalianMappedMethanogens_CleanData_Array.txt

for file in ${cleanedReads_outdir}/*_unmapped_Aligned.sortedByCoord.out.bam; do
	#filename=`basename $file _MethanogensAligned.sortedByCoord.out.bam`
	samtools index $file
done



