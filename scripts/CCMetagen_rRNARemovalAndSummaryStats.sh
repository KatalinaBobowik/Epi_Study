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
	for i in $(seq 0 6); do
		index=`expr $speciesIndex - $i`
		for taxonomicRank in `cat ${folder}/ControlUnmapped_species_table_noRepeats_noRNA_RPM.csv| cut -f$index -d"," | sort | uniq`; do
			allInstances=`grep -i $taxonomicRank ${folder}/*_noRNA.csv | cut -f2 -d"|" | cut -f1 -d" " | sort | uniq | wc -l`
			totalSamples=`grep -i $taxonomicRank ${folder}/*_noRNA.csv | cut -f1 -d":" | sort | uniq | wc -l`
			maxAcessionNumberLength=`grep -i $taxonomicRank ${folder}/*_noRNA.csv | cut -f3 -d"\"" | cut -f4 -d "," | sort -nr | head -n1`
			fullTaxonomicRank=`grep -i $taxonomicRank ${folder}/*_noRNA.csv | awk -F ',' '{i = 8; for (--i; i >= 0; i--){ printf "%s,",$(NF-i)} print ""}' | head -n1`
			printf "`echo $fullTaxonomicRank`\t`echo $taxonomicRank`\t`echo $allInstances`\t`echo $totalSamples`\t`echo $maxAcessionNumberLength`\t\n"
		done > ${folder}/${i}_accessionHits_totalSamples.txt
	done
done

# do the same thing as abive but comine everything into one file --------

# get only accession number and taconomic ranks
for folder in ${folders[@]}; do
	for file in ${folder}/*_noRNA.csv; do
		noHeaderfile=`tail -n +2 $file`
		accession=`echo "$noHeaderfile" | cut -f2 -d"|" | cut -f1 -d"."`
		taxonomicRanks=`echo "$noHeaderfile" | cut -f3 -d"\"" | cut -f13-21 -d","`
		paste <(echo "$accession") <(echo "$taxonomicRanks") --delimiters '\t'
	done > ${folder}/accessionAndTaxonomicRanks.csv
	# get only unique instances of taxonomic ranks
	cat ${folder}/accessionAndTaxonomicRanks.csv | cut -f2 | sort | uniq > ${folder}/uniqueTaxonomicRanks.txt
	# filter out chordata and arthropoda
	grep -v 'Chordata\|Arthropoda' ${folder}/uniqueTaxonomicRanks.txt > ${folder}/uniqueTaxonomicRanks_Filtered.txt
	while read line; do
		kingdom=`echo $line | cut -f2 -d","`
		kingdom_instance=`[[ ! -z "$kingdom" ]] && awk -v pat="$kingdom" -F, '$2==pat' ${folder}/accessionAndTaxonomicRanks.csv | cut -f1 | sort | uniq | wc -l || echo 0`
		phylum=`echo $line | cut -f3 -d","`
		phylum_instance=`[[ ! -z "$phylum" ]] && awk -v pat="$phylum" -F, '$3==pat' ${folder}/accessionAndTaxonomicRanks.csv | cut -f1 | sort | uniq | wc -l || echo 0`
		class=`echo $line | cut -f4 -d","`
		class_instance=`[[ ! -z "$class" ]] && awk -v pat="$class" -F, '$4==pat' ${folder}/accessionAndTaxonomicRanks.csv | cut -f1 | sort | uniq | wc -l || echo 0`
		order=`echo $line | cut -f5 -d","`
		order_instance=`[[ ! -z "$order" ]] && awk -v pat="$order" -F, '$5==pat' ${folder}/accessionAndTaxonomicRanks.csv | cut -f1 | sort | uniq | wc -l || echo 0`
		family=`echo $line | cut -f6 -d","`
		family_instance=`[[ ! -z "$family" ]] && awk -v pat="$family" -F, '$6==pat' ${folder}/accessionAndTaxonomicRanks.csv | cut -f1 | sort | uniq | wc -l || echo 0`
		genus=`echo $line | cut -f7 -d","`
		genus_instance=`[[ ! -z "$genus" ]] && awk -v pat="$genus" -F, '$7==pat' ${folder}/accessionAndTaxonomicRanks.csv | cut -f1 | sort | uniq | wc -l || echo 0`
		species=`echo $line | cut -f8 -d","`
		species_instance=`[[ ! -z "$species" ]] && awk -v pat="$species" -F, '$8==pat' ${folder}/accessionAndTaxonomicRanks.csv | cut -f1 | sort | uniq | wc -l || echo 0`
		printf "`echo $kingdom`\t`echo $kingdom_instance`\t`echo $phylum`\t`echo $phylum_instance`\t`echo $class`\t`echo $class_instance`\t`echo $order`\t`echo $order_instance`\t`echo $family`\t`echo $family_instance`\t`echo $genus`\t`echo $genus_instance`\t`echo $species`\t`echo $species_instance`\n"
	done < ${folder}/uniqueTaxonomicRanks_Filtered.txt > ${folder}/summaryAccessionTable.txt
done

in `cat ${folder}/uniqueTaxonomicRanks_Filtered.txt`; 


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

