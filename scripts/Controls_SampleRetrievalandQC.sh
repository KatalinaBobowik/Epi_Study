# created by KSB, 24.03.2020

# load modules
module load fastqc
module load MultiQC

# specify directory
inputdir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data
refdir=/data/cephfs/punim0586/kbobowik/Sumba/ReferenceFiles/EpiStudy
bindir=/data/cephfs/punim0586/kbobowik/bin
fastQCdir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/FastQC/
trimdir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/trimmomatic
arraydir=/data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/scripts/array_files

# All samples were first downloaded from GEO 
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80974/suppl/GSE80974_RAW.tar'

# untar the file
tar -xvf GSE80974_RAW.tar

# samples which are not neededed are removed
cd ${inputdir}/Controls
for i in *.fastq.gz; do
	ID=${i%.fastq.gz}
	ID2=`echo $ID | cut -d'_' -f1`
    # control samples file from Ruby script GetDiseaseStatus_ControlSamples.sh
    if ! grep -qxFe "$ID2" ${refdir}/ControlSampleIDs.txt; then
        echo "Deleting: ${ID2}_*"
        rm ${ID2}_*
    fi
done
                                                           
# Split reads into two separate files
# use script deinterleave_fastq.sh, taken from  nathanhaigh - https://gist.github.com/nathanhaigh/3521724
for file in ${inputdir}/Controls/*.fastq.gz; do
	filename=`basename $file .gz`
	path=`dirname $file`
	gzip -dc $file | ./deinterleave_fastq.rb ${path}/R1_${filename} ${path}/R2_${filename}
done

# run fastqc
for file in ${inputdir}/Controls/*.fastq; do
	echo $file
	fastqc $file --outdir=$fastQCdir
done

# run multiQC
cd $fastQCdir
multiQC .

# run trimmomatic
for file in ${inputdir}/Controls/R1*.fastq; do
  path=`dirname $file`
  f1=`basename $file`
  f2=`echo "$f1" | sed -r 's/R1/R2/g'`
  echo java -jar ${bindir}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 -trimlog ${inputdir}/trimmedControls/${f1%_*}.fastq.trim.log $path/${f1} $path/${f2} ${inputdir}/trimmedControls/paired_trimmedOutput_${f1} ${inputdir}/trimmedControls/unpaired_trimmedOutput_${f1} ${inputdir}/trimmedControls/paired_trimmedOutput_${f2} ${inputdir}/trimmedControls/unpaired_trimmedOutput_${f2} LEADING:20 TRAILING:20 MINLEN:90 '&>' ${inputdir}/trimmedControls/${f1%_*}.txt
done > ${arraydir}/trimmomaticControls_metagenome.txt

java -jar ${bindir}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 -trimlog /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/${f1%_*}.fastq.trim.log $path/${f1} $path/${f2} /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/paired_trimmedOutput_${f1} /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/unpaired_trimmedOutput_${f1} /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/paired_trimmedOutput_${f2} /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/unpaired_trimmedOutput_${f2} LEADING:20 TRAILING:20 MINLEN:90 &> /data/cephfs/punim0586/kbobowik/EpidemiologicalSurvey/data/test/${f1%_*}.txt
