#Get all the samples (.fastq) from the folder 
allsamples= glob_wildcards("000.fastq/{sample}.fastq").sample

#Check if the samples are detected
print("Detected samples:", allsamples) 

#To run without argument, specifies the final targets
rule all:
    input:
        expand("010.fastqc/{sample}_fastqc.html", sample=allsamples),
        expand("010.fastqc/{sample}_fastqc.zip", sample=allsamples)
    
#Fastqc rule
rule fastqc:
    input:
        "000.fastq/{sample}.fastq"
    output:
        "010.fastqc/{sample}_fastqc.html",
        "010.fastqc/{sample}_fastqc.zip"
    shell:
        "fastqc {input} --outdir=010.fastqc"
