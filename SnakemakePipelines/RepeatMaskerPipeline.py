configfile: "config.yaml"

DIR = ["/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/sample_output", "/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/second_batch/sample_output", "/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/third_batch/sample_output"]
sample=config["samples"]

rule all:
        input:
            	expand("/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta", sample=config["samples"]),
            	expand("/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta", sample=config["samples"]),
            	expand("/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta.masked", sample=config["samples"]),
				expand("/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R2_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta.masked", sample=config["samples"])

rule fastq_to_fasta:
        input:
                r1=lambda wildcards: [os.path.join(DIR[i], 'R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_' + x + 'Aligned.sortedByCoord.out.fastq') for i,x in enumerate(sample) if x == wildcards.sample],
                r2=lambda wildcards: [os.path.join(DIR[i], 'R2_unmapped_STAR_Hg38_second_paired_trimmedOutput_' + x + 'Aligned.sortedByCoord.out.fastq') for i,x in enumerate(sample) if x == wildcards.sample],
        output:
                r1="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta",
                r2="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R2_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta"
        shell:
                "sed -n '1~4s/^@/>/p;2~4p' {input.r1} > {output.r1} ; "
                "sed -n '1~4s/^@/>/p;2~4p' {input.r2} > {output.r2}"


rule RepeatMasker:
		input:
				r1="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta",
                r2="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R2_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta"
		output:
				"/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta.masked",
				"/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R2_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta.masked"

		envs:
				"repeatMasker.yaml"
		shell:
				"/data/cephfs/punim0586/kbobowik/bin/RepeatMasker/RepeatMasker -species human {input.r1} ; "
				"/data/cephfs/punim0586/kbobowik/bin/RepeatMasker/RepeatMasker -species human {input.r2}"

rule KMA:
		input:
				f1="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/allSamples_FastaFiles/R1_unmapped_STAR_Hg38_second_paired_trimmedOutput_{sample}Aligned.sortedByCoord.out.fasta.masked",
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

