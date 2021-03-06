# 25/09/2020

# get control datasets from SRA

# The first study is BioProject PRJNA638819, SRA study SRP266904, GEO study GSE152255
# SRA numbers:
dengueSRA=SRR11994658,SRR11994667,SRR11994676,SRR11994685,SRR11994694,SRR11994703,SRR11994712,SRR11994721,SRR11994730,SRR11994739,SRR11994748
MaliHealthy=SRR1026870
SRR1026873
SRR1026875
SRR1026879
SRR1026881
SRR1026885
SRR1026889
SRR1026891
SRR1026895
SRR1026898
SRR1026900
SRR1026902
SRR1026906
SRR1026910
SRR1026914
SRR1026918
SRR1026922
SRR1026924
SRR1026926
SRR1026928
SRR1026930
SRR1026932
SRR1026936
SRR1026940
SRR1026942
SRR1026944
SRR1026946
SRR1026948
SRR1026950
SRR1026952
SRR1026954
SRR1026956
SRR1026958
SRR1026960
SRR1026962
SRR1026964
SRR1026968
SRR1026970
SRR1026974
SRR1026978
SRR1026980
SRR1026982
SRR1026986
SRR1026988
SRR1026990
SRR1026992
SRR1026994
SRR1026996
SRR1027000
SRR1027002
SRR1027004
SRR1027006
SRR1027008
SRR1027010
rule download_files:
        output:
                "/scratch/punim0586/kat/Bipolar_controlSamples_EpiStudy/{sample}/{sample}_1.fastq.gz"
        log:
                "logs/enaDownload/{sample}.log"
        shell:
                "fasterq-dump {wildcards.sample}"