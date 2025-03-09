#Get all the samples (.fastq) from the folder 
allsamples= glob_wildcards("000.fastq/{sample}.fastq").sample

print("Detected samples:", allsamples) 

rule all:
    input:
        expand("010.fastqc/{sample}_fastqc.html", sample=allsamples),
        expand("010.fastqc/{sample}_fastqc.zip", sample=allsamples)

rule fastqc:
    input:
        "000.fastq/{sample}.fastq"
    output:
        "010.fastqc/{sample}_fastqc.html",
        "010.fastqc/{sample}_fastqc.zip"
    shell:
        "fastqc {input} --outdir=010.fastqc"
