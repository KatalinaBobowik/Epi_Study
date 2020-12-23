# 3/12/2020

# Get Flavivirus sequences for Pai

# set up directories
CCMetagen_dir=/data/scratch/projects/punim0586/kat/101BPIndonesianSamplesCCMetagen_SE
KMA_dir=/data/scratch/projects/punim0586/kat/101BPIndonesianSamplesKMA_SE
Pai_dir=/data/scratch/projects/punim0586/kat/FlavivirusFiles_ForPai

for file in ${CCMetagen_dir}/*csv; do
	ID=`basename $file .csv`
	# get CCMetagen files that have been assigned as Flavivirus and extract their ncbi IDs
	grep Flaviviridae $file | cut -f2 -d"|" | cut -f1 -d" " > ${CCMetagen_dir}/temp.txt
	# Extract Flavivirus sequences from KMA frag files by using the ncbi IDs (obtained above)
	zcat ${KMA_dir}/${ID}.frag.gz | grep -f ${CCMetagen_dir}/temp.txt > ${Pai_dir}/${ID}_Flavivirus.frag
done

# delete empty files
find $Pai_dir -size 0 -print -delete

# compress files
tar cvzf Flavivirus.tar.gz ${Pai_dir}/*.frag
