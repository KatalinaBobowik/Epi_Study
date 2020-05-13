# Created 17.01.19
# Find out what unmapped sequences from STAR bam files are mapping to
 
# Load modules -------------------

# load samtools 1.9
module load SAMtools
# load diamond version 0.9.10
module load diamond

# Filter reads for only unmapped reads ---------------

# the samtools -f flag outputs alignmemnts only matching this flag. The '4' flag (or 00100) corresponds to all unmapped reads (as per the Samtools manual: http://www.htslib.org/doc/samtools-1.2.html and the Picard 'explain samtools flags' website: http://broadinstitute.github.io/picard/explain-flags.html).
# This is the command for one file
samtools view -f 4 alignment.bam > unmapped.bam

# We need to execute this for all of our files. It takes a while to exceutew this for one file (around 3 minutes) so we can make an array script for each file to speed up the process. 
for file in /data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/sample_output/STAR*.bam; do
	filename=`basename $file`
	path=`dirname $file`
	echo samtools view -f 4 $file '>' ${path}/unmapped_${filename}
done > unmapped_array.txt


# Make all unmapped read files into the appropriate Diamond format (i.e., query template name and sequence) --------------

# set up directory name for star files
star="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass"

# extract just the sequence ID and the read via awk, which prints the following components : '>' , then the name of the read, then a name extension of the paired reads (/1 or /2), then a Carriage return, then the sequence and a final carriage return (as per: https://www.biostars.org/p/157434/)
for file in ${star}/{sample_output,second_batch/sample_output,third_batch/sample_output}/unmapped*.bam; do
	# get batch information ftom the 8th field seperator in the sequence file
	batch=`echo $file | cut -d/ -f 8`
	# the first batch of sequences isn't labelled as 'first batch' (there is onlu 'sample_output' in teh 8th field separator) so we need to add this information in. 
	if [[ $batch == "sample_output" ]]; then
		batch="first_batch"
	fi
	shortenedFile=`basename $file Aligned.sortedByCoord.out.bam`
	# ${variable##pattern} is like $variable, minus the longest matching pattern from front-end
	sampleID=${shortenedFile##*_}
	batchPlusID="${batch}_${sampleID}"
	samtools view $file | awk '{printf(">%s/%s\n%s\n",$1,(and(int($2),0x40)?1:2),$10);}' > /data/cephfs/punim0586/kbobowik/Diamond/unmapped_reads/${batchPlusID}_unmapped
	echo ${sampleID} done
done

# check to make sure there are 123 samples in all
ls * | wc -l
# 123 files are all present so we can now run Diamond

# Run Diamond ------------------------

# set up directory name for Diamond
diamond="/data/cephfs/punim0586/kbobowik/Diamond"
# initial diamon setup - first set up database
# download the database as a single file nr.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast//db/FASTA/nr.gz

# create the database
diamond makedb --in nr.gz -d nr

# make an array script to speed up the process
for file in /data/cephfs/punim0586/kbobowik/Diamond/unmapped_reads/*_unmapped; do
	sample=`basename $file`
	sampleID=${sample%_*}
	# here, we're specifying the database with the '-d' flag. This is set to nr.dmnd. The '-q' flag specifies the query sequence, which is our filtered unmapped read files form the previous step. The '-o' flag specifies the output, which we're putting in the Diamond folder under 'matches_output' sub-folder. Finally, we'll specify '--top 1' which will report alignments within 1% of the top alignment score for each query
	echo diamond blastx -d /data/cephfs/punim0586/kbobowik/Diamond/nr.dmnd -q $file -o /data/cephfs/punim0586/kbobowik/Diamond/matches_output/matches_${sampleID}.m8 --top 1
done > ${diamond}/diamond_array_script.txt

# Get a summary of how often each accession number appears for the top 1000 hits -----------------------

# Now that we've run Diamond and have sequence matches for our query sequences, we want to extract the accession number and find out how many times it occurs in each file. we'll do this for the top 1000 hits, since any more will start to get messy.
# Get the second columns (which is the seqence ID) and count how many times each accession number appears.
cut -f 2 matches_first_batch_MPI-025.m8 | sort  | uniq -c | sort -rn | head -n 1000 > summary.txt

# we can do this for all samples in an array
for file in ${diamond}/matches_output/*.m8; do
	filename=`basename $file ".m8"`
	echo cut -f 2 $file '|' sort '|' uniq -c '|' sort -rn '|' head -n 1000 '>' ${diamond}/matches_output/${filename}_top1000hits_summary.txt
done > ${diamond}/accessionNumber_summaryArray.txt

# for some reason, the delimeter is a space so we'll change this to tab-delimited to make things easier downstream
for file in ${diamond}/matches_output/*_top1000hits_summary.txt; do
	filename=`basename $file ".txt"`
	awk -v OFS="\t" '$1=$1' $file > ${diamond}/matches_output/${filename}_tabDelimited.txt
done

# Get organism names ------------------------------------

# from the top 1000 hits accession file, we'll get the second field (the accession number) and use the file database file 'nr.oneline.out' to match accession numbers to organism names
for file in ${diamond}/matches_output/*_top1000hits_summary.txt; do
	filename=`basename $file "_top1000hits_summary.txt"`
	echo awk \'{print \$\2}\' $file '|' grep -w -f - /data/cephfs/punim0586/kbobowik/Diamond/nr/nr.oneline.out '>' ${diamond}/matches_output/${filename}_accessionNumber_OrganismName.txt
done > ${diamond}/MatchedOrganismName_Array.txt

# Now all we have to do is tidy up our file so that all we have are the accession numbers and organism names
for file in ${diamond}/matches_output/*_accessionNumber_OrganismName.txt; do
	filename=`basename $file "_accessionNumber_OrganismName.txt"`
	path=`dirname $file`
	# we will be extracting the number of times each accession number occurred from this file
	summaryFile=${path}/${filename}_top1000hits_summary_tabDelimited.txt
	while IFS= read -r line; do
		# get the organism name from each line
		org=`echo $line | awk 'NR>1{print $1}' RS='[' FS=']'`
		# get the accession number from each line and get rid of the '>' sign in front
		acc=`echo $line | grep -Eo '^[^ ]+' | sed 's/^[ >]*//'` 
		if [[ -z "$org" ]]; then
			# If there is no scientfic name inbetween the brackets, just return the accession number
			echo -e "$acc\t$acc"
		else
			# else, give the accession number and organism name
			echo -e "$acc\t$org" 
		fi
	# read output into a temporary file 
	done<$file>organisedAccOutput.txt
	# Get the number of accession number occurrences from the summary file and merge with the organism name + accession number file. Save this as a new file
	awk -v FILE_A=$summaryFile -v OFS="\t" 'BEGIN { while ( ( getline < FILE_A ) > 0 ) { VAL = $0 ; sub( /^[^ ]+ /, "", VAL ) ; DICT[ $2 ] = VAL } } { print $0, DICT[ $1 ] }' organisedAccOutput.txt > ${path}/${filename}_AccOrgCount.txt
done



