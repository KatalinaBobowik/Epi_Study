# 02/03/2020

# Run conda session (if running sinteractive session)
source activate py37

# directories
Control50BP=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/controlSamples_SnakemakeWorkflow/output/50BPControlSamplesCCMetagen_TE
Indo101BP=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/101BPTrimmedIndonesianSamples/output/101BPIndonesianSamplesCCMetagen_TE
Indo50BP=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/50BPTrimmedIndonesianSamples/output/CCMetagenOutput_TE

# remove all ribosomal RNA hits --------------------------

# make no RNA directory for all three datasets
mkdir ${Control50BP}_noRNA/ ${Indo101BP}_noRNA/ ${Indo50BP}_noRNA/
# copy folde into new directories
cp ${Control50BP}/* ${Control50BP}_noRNA/
cp ${Indo101BP}/* ${Indo101BP}_noRNA/
cp ${Indo50BP}/* ${Indo50BP}_noRNA/

# remove ribosomal hits from CCMetagen csv file
export folders=(${Control50BP}_noRNA ${Indo101BP}_noRNA ${Indo50BP}_noRNA)

for folder in ${folders[@]}; do
	for sample in ${folder}/*.csv; do
		ID=`basename $sample .csv`
		grep -iv Ribosomal $sample | grep -iv rRNA > ${folder}/${ID}_noRNA.csv
	done
	# remove orginal .csv files
	cd $folder | ls | grep '.csv$' | grep -v 'noRNA.csv$' | xargs rm
	# Produce summary table
	CCMetagen_merge.py --input_fp ${folder}/ --keep_or_remove r --filtering_tax_level Phylum --taxa_list Chordata,Arthropoda --output_fp ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_RPM
done

# get total number of accession hits and total number of samples for each species
IFS=$'\n'
for folder in ${folders[@]}; do
	speciesIndex=`head -n1 ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_RPM.csv |tr ',' '\n'|wc -l`
	for species in `cat ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_RPM.csv| cut -f$speciesIndex -d"," | sort | uniq`; do
		allInstances=`grep -i $species ${folder}/*_noRNA.csv | cut -f2 -d"." | sort | uniq | printf "\`\echo $species\` \`wc -l\`\n"`
		totalSamples=`grep -i $species ${folder}/*_noRNA.csv | cut -f1 -d":" | sort | uniq | wc -l`
		printf "`echo $allInstances` `echo $totalSamples` \n"
	done > ${folder}/species_accessionHits_totalSamples.txt
done



