# 02/03/2020

# Run conda session (if running sinteractive session)
source activate py37

# run CCMetagen
inputdir=/data/cephfs/punim0586/kbobowik/CCmetagen/IndoSample
array_dir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files

for sample in ${inputdir}/*.res; do
	sampleID=$(basename $sample .res)
	#echo CCMetagen.py -i $sample -o ~/CCMetagenOutput/${sampleID}
	echo CCMetagen.py -i $sample -o ~/CCMetagenOutput/${sampleID} --depth_unit rpm --mapstat ${sample/.res}.mapstat --depth 1
done > ${array_dir}/CCMetagen_AllIndoFiles.txt

# Produce summary table
CCMetagen_merge.py --input_fp ~/CCMetagenOutput/ --keep_or_remove r --filtering_tax_level Phylum --taxa_list Chordata,Arthropoda --output_fp IndoUnmapped_species_table_Helminths_RPM



