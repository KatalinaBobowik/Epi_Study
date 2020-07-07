configfile: "config.yaml"

rule all:
        input:
                expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_1.fastq", sample=config["samples"]),
                expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_2.fastq", sample=config["samples"]),
                expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}.fastq", sample=config["samples"]),
                expand("qc/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
                expand("qc/fastqc/{sample}_2_fastqc.html", sample=config["samples"]),
                expand("qc/fastqc/{sample}_fastqc.html", sample=config["samples"]),
                expand("qc/fastqc/{sample}_1_fastqc.zip", sample=config["samples"]),
                expand("qc/fastqc/{sample}_2_fastqc.zip", sample=config["samples"]),
                expand("qc/fastqc/{sample}_fastqc.zip", sample=config["samples"]),
                expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.fastq", sample=config["samples"]),
                expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.fastq", sample=config["samples"]),
                expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.unpaired.fastq", sample=config["samples"]),
                expand("/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.unpaired.fastq", sample=config["samples"]),
                directory("genome"),
                expand("star/firstpass/{sample}_Aligned.sortedByCoord.out.bam", sample=config["samples"]),
                "star/STAR_Hg38_allFilesSJ.out.tab",
                directory("genome_secondPass"),
                expand("star/secondPass/{sample}_Aligned.sortedByCoord.out.bam", sample=config["samples"])

rule download_files:
        output:
                "/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_1.fastq",
                "/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_2.fastq",
                "/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}.fastq"
        conda:
                "envs/download.yaml"
        shell:
                "fasterq-dump {wildcards.sample} -O /scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{wildcards.sample}/"

rule fastqc:
        input:
                ["/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_1.fastq", "/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_2.fastq", "/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}.fastq"]
        output:
                "qc/fastqc/{sample}_1_fastqc.html",
                "qc/fastqc/{sample}_2_fastqc.html",
                "qc/fastqc/{sample}_fastqc.html",
                "qc/fastqc/{sample}_1_fastqc.zip",
                "qc/fastqc/{sample}_2_fastqc.zip",
                "qc/fastqc/{sample}_fastqc.zip"
        conda:
                "envs/QC.yaml"
        shell:
                "module load web_proxy; fastqc {input} --extract --outdir=qc/fastqc/"

rule trimmomatic_pe:
        input:
                r1="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_1.fastq",
                r2="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_2.fastq"
        output:
                r1="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.fastq",
                r2="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.fastq",
                # reads where trimming entirely removed the mate
                r1_unpaired="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.unpaired.fastq",
                r2_unpaired="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.unpaired.fastq"
        log:
                "logs/trimmomatic/{sample}.log"
        conda:
                "envs/trimmomatic.yaml"
        shell:
                "module load web_proxy; module load Java; java -jar /data/cephfs/punim0586/kbobowik/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE "
                "-threads 6 -phred33 {input.r1} {input.r2} {output.r1} "
                "{output.r1_unpaired} {output.r2} {output.r2_unpaired} LEADING:20 TRAILING:20 MINLEN:45 -trimlog {log}"
rule star_index:
        input:
                fasta="/data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa",
                gtf="/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf"
        output:
                directory("genome")
        conda:
                "envs/STAR.yaml"
        log:
                "logs/star_index/genome.log"
        shell:
                "module load web_proxy; mkdir genome; STAR --runMode genomeGenerate "
                "--genomeDir genome "
                "--genomeFastaFiles {input.fasta} "
                "--sjdbGTFfile {input.gtf} "
                "--sjdbOverhang 49 "
                "--runThreadN 12 {log}"

rule star_pe_firstmap:
        input:
                f1="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.fastq",
                f2="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.fastq"
        output:
                "star/firstpass/{sample}_Aligned.sortedByCoord.out.bam"
        log:
                "logs/star/{sample}.log"
        conda:
                "envs/STAR.yaml"
        shell:
                # unfortunately I have no better way of doing it- define sample name then use that as dirname
                "module load web_proxy; x=`basename {input.f1} _1.fastq`; STAR --genomeDir genome "
                "--readFilesIn {input.f1} {input.f2} "
                "--runThreadN 8 --sjdbOverhang 49 "
                "--outFileNamePrefix star/firstpass/${{x}}_ "
                "--outSAMtype BAM SortedByCoordinate &> {log}"

rule combine_splice_junctions:
        # combine all of the SJ.out.tab files from from all files created in first pass into one file
        input:
                ""
        output:
                "star/STAR_Hg38_allFilesSJ.out.tab"
        shell:
                "cat star/firstpass/*_SJ.out.tab > {output}"

rule star_index_second:
        input:
                fasta="/data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa",
                gtf="/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf"
        output:
                directory("genome_secondPass")
                                                                                         116,3-17      75%
        conda:
                "envs/STAR.yaml"
        shell:
                "module load web_proxy; mkdir genome_secondPass; STAR --runMode genomeGenerate "
                "--genomeDir genome_secondPass " # Path to output
                "--genomeFastaFiles {input.fasta} " # Path to the fasta file with the genome sequences
                "--sjdbGTFfile {input.gtf} " # path to the GTF file with annotations
                "--sjdbOverhang 49 " # this length should be equal to the ReadLength-1
                "--sjdbFileChrStartEnd star/STAR_Hg38_allFilesSJ.out.tab " #this is the file we used to get all of the combined splice junctions
                "--runThreadN 12"

rule star_pe_secondmap:
        input:
                f1="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_1.fastq",
                f2="/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy_Trimmed/{sample}_2.fastq"
        output:
                "star/secondPass/{sample}_Aligned.sortedByCoord.out.bam"
        log:
                "logs/star/{sample}.log"
        conda:
                "envs/STAR.yaml"
        shell:
                # unfortunately I have no better way of doing it- define sample name then use that as dirname
               "module load web_proxy; x=`basename {input.f1} _1.fastq`; STAR --genomeDir genome_secondPass "
                "--readFilesIn {input.f1} {input.f2} "
                "--runThreadN 8 --sjdbOverhang 49 "
                "--outFileNamePrefix star/secondPass/${{x}}_ "
                "--outSAMunmapped Within "
                "--outReadsUnmapped Fastx "
                "--outSAMtype BAM SortedByCoordinate"
                                                                                         146,3-17      Bot
