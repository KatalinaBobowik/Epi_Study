# Created 17.01.19
# Find out what unmapped sequences from STAR bam files are mapping to
 
# Load modules -------------------

# load samtools 1.9
module load SAMtools
# load diamond version 0.9.10
module load diamond


# Run Diamond ------------------------

# set up directory name for Diamond
diamond="/data/cephfs/punim0586/kbobowik/Diamond"

# make an array script to speed up the process
for file in /data/cephfs/punim0586/kbobowik/Diamond/unmapped_reads/*_unmapped; do
	sample=`basename $file`
	sampleID=${sample%_*}
	# here, we're specifying the database with the '-d' flag. This is set to nr.dmnd. The '-q' flag specifies the query sequence, which is our filtered unmapped read files form the previous step. The '-o' flag specifies the output, which we're putting in the Diamond folder under 'matches_output' sub-folder. Finally, we'll specify '--top 1' which will report alignments within 1% of the top alignment score for each query
	echo diamond blastx -d /data/cephfs/punim0586/kbobowik/Diamond/nr.dmnd -q $file -o /data/cephfs/punim0586/kbobowik/Diamond/matches_output/matches_oneHit_${sampleID}.m8 -k 1 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore stitle
done > ${diamond}/diamond_array_script_OneHit.txt

# Get a summary of how often each accession number appears for the top 1000 hits -----------------------

# we can do this for all samples in an array
for file in ${diamond}/matches_output/matches_oneHit*.m8; do
	filename=`basename $file ".m8"`
	# get 11th column, which is the accession number and organism name
	cut -f 11 $file | sort | uniq -c | sort -rn | head -n 1000 > ${diamond}/matches_output/${filename}_top1000hits_summary.txt
done


# Get organism names ------------------------------------


# Now all we have to do is tidy up our file so that all we have are the accession numbers and organism names
for file in ${diamond}/matches_output/matches_oneHit_*_top1000hits_summary.txt; do
	filename=`basename $file "_top1000hits_summary.txt"`
	path=`dirname $file`
	# we will be extracting the number of times each accession number occurred from this file
	while IFS= read -r line; do
		# get the organism name from each line
		org=`echo $line | awk -F'[][]' '{print $2}'`
		# get the accession number from each line and get rid of the '>' sign in front
		freq=`echo $line | awk '{print $1}'`
		acc=`echo $line | awk '{print $2}'`
		if [[ -z "$org" ]]; then
			# If there is no scientfic name inbetween the brackets, just return the accession number
			echo -e "$acc\t$acc\t$freq"
		else
			# else, give the accession number and organism name
			echo -e "$acc\t$org\t$freq"
		fi
	# read output into a temporary file 
	done<$file>${filename}_AccOrgCount.txt
done




