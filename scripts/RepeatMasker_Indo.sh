# 02/03/2020

# Run RepeatMasker on all samples in order to get reads which may be mapping to other organisms
# due to sequence repeats

repeatMasker=/data/cephfs/punim0586/kbobowik/bin/RepeatMasker
ccmetagenDir=/home/kbobowik/CCMetagenoutputNorRNA
repeatMaskerOutdir=/data/cephfs/punim0586/kbobowik/CCmetagen/IndoSample
arraydir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files

module load Perl
# get .fsa files which are the output from the KMA algorithm
for file in ${repeatMaskerOutdir}/*.fsa; do
	# run repeatmasker on the KMA fasta files
	echo ${repeatMasker}/RepeatMasker $file
done > ${arraydir}/RepeatMasker_IndoSamples.txt

# for each sample, get accession numbers coming from organisms that had repeats 
for file in ${repeatMaskerOutdir}/*.fsa.out; do
	ID=`basename $file .fsa.out`
	accession=`cat $file | cut -d '|' -f 2 |  cut -d ' ' -f 1 | sort | uniq`
	echo $accession | sed 's/ /\\|/g' > ${repeatMaskerOutdir}/${ID}.txt
	# grep `cat ${ID}.txt` ${ccmetagenDir}/${ID}_noRNA.csv | cut -d ' ' -f 2,3 | sort | uniq | grep -v "Human" | grep -v "Homo" | grep -v "Pan" > ${ccmetagenDir}/${ID}_OrganismsWithRepeats.txt
done

# remove these accession numbers from the CCmetagen no RNA files
for file in ${repeatMaskerOutdir}/*.txt; do
	ID=`basename $file .txt`
	grep -vf $file ${ccmetagenDir}/${ID}_noRNA.csv > /home/kbobowik/CCMetagenoutputNoRepeats/${ID}_OrganismsNoRepeats.csv
done






