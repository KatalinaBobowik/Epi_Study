#!/bin/bash

# Run conda session (if running sinteractive session)
eval "$(conda shell.bash hook)"
conda activate py37

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
	find $folder -name '*.csv' -and -not -name '*noRNA.csv' -type f -exec rm '{}' \;
	# Produce summary table
	CCMetagen_merge.py --input_fp ${folder}/ --keep_or_remove r --filtering_tax_level Phylum --taxa_list Chordata,Arthropoda --output_fp ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_RPM
done

# get total number of accession hits and total number of samples for each species an order
IFS=$'\n'
for folder in ${folders[@]}; do
	speciesIndex=`head -n1 ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_RPM.csv |tr ',' '\n'|wc -l`
	for species in `cat ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_RPM.csv| cut -f$speciesIndex -d"," | sort | uniq`; do
		allInstances=`grep -i $species ${folder}/*_noRNA.csv | cut -f2 -d"|" | cut -f1 -d" " | sort | uniq | wc -l`
		totalSamples=`grep -i $species ${folder}/*_noRNA.csv | cut -f1 -d":" | sort | uniq | wc -l`
		maxAcessionNumberLength=`grep -i $species ${folder}/*_noRNA.csv | cut -f3 -d"\"" | cut -f4 -d "," | sort -nr | head -n1`
		printf "`echo $species`\t`echo $allInstances`\t`echo $totalSamples`\t`echo $maxAcessionNumberLength`\t\n"
	done > ${folder}/species_accessionHits_totalSamples.txt
	orderIndex=`echo "$(($speciesIndex-3))"`
	for order in `cat ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_RPM.csv| cut -f$orderIndex -d"," | sort | uniq`; do
		allInstances=`grep -i $order ${folder}/*_noRNA.csv | cut -f2 -d"|" | cut -f1 -d" " | sort | uniq | wc -l`
		totalSamples=`grep -i $order ${folder}/*_noRNA.csv | cut -f1 -d":" | sort | uniq | wc -l`
		maxAcessionNumberLength=`grep -i $family ${folder}/*_noRNA.csv | cut -f3 -d"\"" | cut -f4 -d "," | sort -nr | head -n1`
		printf "`echo $order`\t`echo $allInstances`\t`echo $totalSamples`\t`echo $maxAcessionNumberLength`\t\n"
	done > ${folder}/order_accessionHits_totalSamples.txt
done

# Remove pathogens from CCMetagen merged file where reads only map to one accession number and are less than 1000BP
for folder in ${folders[@]}; do
	for sample in ${folder}/*_noRNA.csv; do
		ID=`basename $sample _noRNA.csv`
		grep -vf <(awk -F "\t" '{if($2==1)print $0}' ${folder}/order_accessionHits_totalSamples.txt | cut -f1) ${folder}/${ID}_noRNA.csv > ${folder}/${ID}_lowAccessionRemoval_order.csv
	done
	# remove orginal .csv files
	find $folder -name '*noRNA.csv' -and -not -name '*lowAccessionRemoval_order.csv' -type f -exec rm '{}' \;
	# Remove original csv file made by CCMetagen_merge
	find $folder -name 'ControlUnmapped_species_table_noRepeats_noRNA_RPM.csv' -type f -exec rm '{}' \;
	# Produce summary table
	CCMetagen_merge.py --input_fp ${folder}/ --keep_or_remove r --filtering_tax_level Phylum --taxa_list Chordata,Arthropoda --output_fp ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_noLowAccession_RPM
done

