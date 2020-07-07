configfile: "config.yaml"

rule all:
	input:
		#expand("qc/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
		#expand("qc/fastqc/{sample}_2_fastqc.html", sample=config["samples"]),
		#expand("qc/fastqc/{sample}_1_fastqc.zip", sample=config["samples"]),
		#expand("qc/fastqc/{sample}_2_fastqc.zip", sample=config["samples"]),
		#expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.fastq.gz", sample=config["samples"]),
		#expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.fastq.gz", sample=config["samples"]),
		#expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.unpaired.fastq.gz", sample=config["samples"]),
		#expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.unpaired.fastq.gz", sample=config["samples"]),
		directory("genome"),
		expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/{sample}_Aligned.sortedByCoord.out.bam", sample=config["samples"]),
		expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/{sample}_Unmapped.out.mate1", sample=config["samples"]),
		expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/{sample}_Unmapped.out.mate1", sample=config["samples"]),
		expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_KMA/{sample}.res", sample=config["samples"]),
		expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_KMA/{sample}.mapstat", sample=config["samples"]),
		expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_CCMetagen/{sample}.csv", sample=config["samples"]),
		expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_CCMetagen/{sample}.html", sample=config["samples"]),
		expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_CCMetagen/{sample}.tsv", sample=config["samples"])

rule trimmomatic_pe:
	input:
		r1=expand("/data/cephfs/punim0586/shared/raw_data/{dir}/{sample}.txt", dir=["indoRNA", "indoRNA_second_batch", "indoRNA_third_batch"], sample=config["samples"]),
		r2=expand("/data/cephfs/punim0586/shared/raw_data/{dir}/{sample}.txt", dir=["indoRNA", "indoRNA_second_batch", "indoRNA_third_batch"], sample=config["samples"]),
	output:
		r1="/scratch/punim0586/kat/IndonesianSamplesTrimmed/{sample}_1.fastq.gz",
		r2="/scratch/punim0586/kat/IndonesianSamplesTrimmed/{sample}_2.fastq.gz",
		# reads where trimming entirely removed the mate
		r1_unpaired="/scratch/punim0586/kat/IndonesianSamplesTrimmed/{sample}_1.unpaired.fastq.gz",
		r2_unpaired="/scratch/punim0586/kat/IndonesianSamplesTrimmed/{sample}_2.unpaired.fastq.gz"
	log:
		"logs/trimmomatic/{sample}.log"
	shell:
		"module load web_proxy; module load Java; java -jar /data/cephfs/punim0586/kbobowik/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE " 
		"-threads 12 -phred33 {input.r1} {input.r2} {output.r1} "
		"{output.r1_unpaired} {output.r2} {output.r2_unpaired} CROP:50 LEADING:20 TRAILING:20 MINLEN:45 -trimlog {log}"

rule fastqc:
	input:
		["/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_1.fastq.gz", "/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_2.fastq.gz"]
	output:
        	"qc/fastqc/{sample}_1_fastqc.html",
		"qc/fastqc/{sample}_2_fastqc.html",
		"qc/fastqc/{sample}_1_fastqc.zip",
		"qc/fastqc/{sample}_2_fastqc.zip"
	conda:
		"envs/QC.yaml"
	log:
		"logs/fastqc/{sample}.log"
	shell:
        	"module load web_proxy; fastqc {input} -t 12 --extract --outdir=qc/fastqc/ &> {log}"

rule star_index:
	input:
		fasta="/data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa",
		gtf="/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf"
	output:
		directory("genome")
	log:
		"logs/star_index/genome.log"
	shell:
		"module load web_proxy; module load STAR; mkdir genome; STAR --runMode genomeGenerate "
		"--genomeDir genome "
		"--genomeFastaFiles {input.fasta} "
		"--sjdbGTFfile {input.gtf} "
		"--sjdbOverhang 49 "
		"--runThreadN 12 &> {log}"

rule star_pe_firstmap:
	input:
		f1="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.fastq.gz",
		f2="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.fastq.gz",
		dir=directory("genome")
	output:
		"/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/{sample}_Aligned.sortedByCoord.out.bam",
		"/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/{sample}_Unmapped.out.mate1",
		"/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/{sample}_Unmapped.out.mate2"
	log:
		"logs/star/{sample}.log"
	shell:
		# unfortunately I have no better way of doing it- define sample name then use that as dirname
		"module load web_proxy; module load STAR; x=`basename {input.f1} _1.fastq.gz`; STAR --genomeDir {input.dir} "
		"--readFilesIn {input.f1} {input.f2} "
		"--runThreadN 12 --sjdbOverhang 49 "
		"--readFilesCommand zcat "
		"--outFileNamePrefix /scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/${{x}}_ "
		"--outSAMunmapped Within "
		"--outReadsUnmapped Fastx "
		"--outSAMtype BAM SortedByCoordinate &> {log}"

rule KMA:
	input:
		f1="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/{sample}_Unmapped.out.mate1",
		f2="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Mapped/{sample}_Unmapped.out.mate2"
	output:
		"/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_KMA/{sample}.res"
		"/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_KMA/{sample}.mapstat"
	shell:
		"/data/cephfs/punim0586/kbobowik/CCmetagen/kma/kma "
		"-ipe {input.f1} {input.f2} " # -ipe Inputfile(s), paired end
		"-o /scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_KMA/{sample} " # -o output
		"-ef " # -ef used to calculate abundance in reads per million (RPM)
		"-t_db /scratch/punim0586/kat/ncbi_nt_no_env_11jun2019 " # -t_db database; 
		"-t 12 -1t1 " # -t threads; -1t1 One read to one template, no splicing performed
		"-mem_mode " # -mem_mode *.index and *.seq are not loaded into memory, which enables one to map against larger databases
		"-apm f " # -apm Paired end method, “f” force paired reads to pair
		"-and" # -and Both mrs and p_value thresholds has to reached to in order to report a template hit

rule CCMetagen:
	input:
		resultFile="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_KMA/{sample}.res",
		mapstatFile="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_KMA/{sample}.mapstat"
	output:
		prefix="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_CCMetagen/{sample}",
		csv="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_CCMetagen/{sample}.csv",
		html="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_CCMetagen/{sample}.html",
		tsv="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_CCMetagen/{sample}.tsv"
	shell:
		"source activate py37 "
		"CCMetagen.py -i {input.resultFile} -o {output.prefix} "
		"--depth_unit rpm " # depth is in reads per million
		"--mapstat {input.mapstatFile} " # if results are reported in reads per million, the mapstat file generated by KMA needs to be supplied
		"--depth 1" # filter out matches with less than one read per million

