# 02/03/2020

# Run conda session (if running sinteractive session)
source activate py37

# run CCMetagen
inputdir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/controlSamples_SnakemakeWorkflow/output/50BPControlSamplesCCMetagen_TE
array_dir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files
noRNA_dir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/101BPTrimmedIndonesianSamples/output/101BPIndonesianSamplesCCMetagen_TE_noRNA

for sample in ${inputdir}/*.res; do
	sampleID=$(basename $sample .res)
	#echo CCMetagen.py -i $sample -o ~/CCMetagenOutput/${sampleID}
	echo CCMetagen.py -i $sample -o ~/CCMetagenOutput/${sampleID} --depth_unit rpm --mapstat ${sample/.res}.mapstat --depth 1
done > ${array_dir}/CCMetagen_AllIndoFiles.txt