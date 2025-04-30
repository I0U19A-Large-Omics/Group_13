# Setup paths and constants
genome_db = "/staging/leuven/stg_00079/teaching/hg38_9/chr9.fa"
snpeff_jar = ("/lustre1/project/stg_00079/teaching/I0U19a_conda_2025/share/snpeff-5.2-1/snpEff.jar")
snpeff_db_folder = "/staging/leuven/stg_00079/teaching/snpeff_db/"

# Change this to match your r-number's last digit
GROUP_NUMBER = 6 

# Define your user ID for iRODS upload folder
USER_ID = "r1017376" 

# Define iRODS paths
IRODS_SOURCE_PATH = f"/gbiomed/home/large_omics_course/fastq/group_{GROUP_NUMBER}"
IRODS_OUTPUT_PATH = f"/gbiomed/home/large_omics_course/output/{USER_ID}"

# Add this function after your checkpoint rule
def get_sample_names(wildcards):
    checkpoint_output = checkpoints.download_from_irods.get(**wildcards).output[0]
    sample_names, = glob_wildcards(f"000.fastq/group_{GROUP_NUMBER}/{{name}}.fq.gz")
    return sample_names

# Function to get BAM files for variant calling
def get_bam_files(wildcards):
    sample_names = get_sample_names(wildcards)
    return expand("020.bwa/{sample}.bam", sample=sample_names)

# This rule now uses checkpoints to handle downloading first, then processing
rule all:
    input:
        download_complete="000.fastq/download_complete.flag",
        vcf="050.snpeff/snps.annotated.vcf",
        database="060.database/variants.db",
        upload_complete="upload_complete.flag"

# Rule to download data from iRODS
checkpoint download_from_irods:
    output:
        flag="000.fastq/download_complete.flag"
    shell:
        """
        # Create directory for FASTQ files
        mkdir -p 000.fastq
        
        # Download all FASTQ files from iRODS
        iget -r "{IRODS_SOURCE_PATH}" 000.fastq/

        # Create flag file indicating downloads are complete
        touch {output.flag}
        """


rule fastqc:
    input:
        fq=f"000.fastq/group_{GROUP_NUMBER}/{{name}}.fq.gz",
        download_complete="000.fastq/download_complete.flag"

    output:
        zip="010.fastqc/{name}_fastqc.zip",
        html="010.fastqc/{name}_fastqc.html",
        summarytxt="010.fastqc/{name}_fastqc/summary.txt",
        basequal=report("010.fastqc/{name}_fastqc/Images/per_base_quality.png",
                        category='FastQC',
                        subcategory='Per base quality',
                        labels={"sample": "{name}"})
    shell:
        """
        echo "Input Fastq: {input.fq} "
        fastqc -o 010.fastqc {input.fq}

        echo "Unzip the output"
        ( cd 010.fastqc ; unzip -o {wildcards.name}_fastqc.zip )
        """

rule fastqc_report_image:
    input:
        summarytxt="010.fastqc/{name}_fastqc/summary.txt"
    output:
        statuspng=report("010.fastqc/{name}_fastqc/summary.png",
                         category='FastQC',
                         subcategory='Status',
                         labels={"sample": "{name}"})
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        #load data
        data = pd.read_csv(input.summarytxt, sep="\t", header=None)
        data.columns = ['status', 'test', 'sample']

        #assign dummy x value for scatterplot
        data['x'] = 1

        #create image
        fig = plt.figure(figsize=(4,5))
        ax = plt.gca()
        sns.scatterplot(data, x='x', y='test', hue='status', s=200, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.title(wildcards.name)
        plt.savefig(output.statuspng)

rule bwa:
    input:
        fq=f"000.fastq/group_{GROUP_NUMBER}/{{sample}}.fq.gz",
        download_complete="000.fastq/download_complete.flag"

    output:
        bam="020.bwa/{sample}.bam",
        bai="020.bwa/{sample}.bam.bai",
    params:
        db=genome_db,
    shell:
        """
        bwa mem {params.db} {input.fq} \
            | samtools sort - \
            > {output.bam}
        samtools index {output.bam}
        """

rule variant_calling:
    input:
        db=genome_db,
        bams=get_bam_files,
    output:
        vcf="030.samtools/snps.vcf",
    shell:
        """
        bcftools mpileup -Ou -f {input.db} {input.bams} \
             | bcftools call -mv -Ov -o {output.vcf}
        """

rule variant_cleanup:
    input:
        db=genome_db,
        vcf="030.samtools/snps.vcf"
    output:
        vcf="040.cleaned/snps.cleaned.vcf"
    shell:
        """
        ( cat {input.vcf} \
           | vt decompose - \
           | vt normalize -n -r {input.db} - \
           | vt uniq - \
           | vt view -f "QUAL>20" -h - \
           > {output.vcf} )
        """

rule snpeff:
    input:
        vcf = "040.cleaned/snps.cleaned.vcf",
    params:
        snpeff_db_folder = snpeff_db_folder,
        snpeff_jar = snpeff_jar,
    log:
        err="050.snpeff/snakemake.err",
    output:
        vcf = "050.snpeff/snps.annotated.vcf",
        html = "050.snpeff/snpEff_summary.html",
        genetxt = "050.snpeff/snpEff_genes.txt",
    shell:
        """

        mkdir -p 050.snpeff

        java -Xmx4096m -jar \
            {params.snpeff_jar} eff hg19 \
            -dataDir {params.snpeff_db_folder} \
            {input.vcf} > {output.vcf}

        # move output files to the snpeff output folder
        mv snpEff_genes.txt snpEff_summary.html 050.snpeff

        """

rule vcf_to_sqlite:
    input:
        vcf="050.snpeff/snps.annotated.vcf"
    output:
        db="060.database/variants.db"
    run:
        import sqlite3
        import os
        import pandas as pd
        import vcfpy  # Using vcfpy like in ParseVCF
        
        # Create output directory if it doesn't exist
        os.makedirs("060.database", exist_ok=True)
        
        # Connect to SQLite database
        conn = sqlite3.connect(output.db)
        
        # Empty lists to store the data (like in ParseVCF)
        snp_records = []
        effect_records = []
        call_records = []
        
        # These are the columns from a SNPeff record
        effect_rec_names = """snp allele effect impact gene gene_id feature_type feature_id 
                          biotype rank hgvs.c hgvs.p cdna_pos cds_pos prot_pos distance_to_feature
                          messages""".split()
        
        # Open the vcf iterator using vcfpy
        reader = vcfpy.Reader.from_path(input.vcf)
        
        # Process records
        for record in reader:
            # We expect no multiallelic SNPs (just like ParseVCF)
            assert len(record.ALT) == 1
            
            # Get the ALT allele
            alt = record.ALT[0]
            
            # Create SNP name for joins
            snp_name = f"{record.CHROM}:{record.POS}:{record.REF}:{alt.value}"
            
            # Store SNP record
            snp_records.append(
                dict(snp=snp_name,
                     chrom=record.CHROM,
                     pos=record.POS,
                     quality=record.QUAL,
                     ref=record.REF,
                     type=alt.type,
                     alt=alt.value))
            
            # Process calls (samples)
            for call_record in record.calls:
                # Get sample name
                sample = os.path.basename(call_record.sample)
                sample = sample.replace('.bam', '')
                
                # Calculate genotype_simple (count of ALT alleles)
                genotype_simple = call_record.data['GT'].count('1')
                
                # Store call record
                call_records.append(
                    dict(snp=snp_name,
                         sample=sample,
                         genotype=call_record.data['GT'],
                         genotype_simple=genotype_simple))
            
            # Process SnpEff annotations
            if 'ANN' in record.INFO:
                for ann in record.INFO['ANN']:
                    ann_fields = ann.split('|')
                    # Create a dictionary of all fields
                    eff_record = dict(zip(effect_rec_names, [snp_name] + ann_fields))
                    
                    # Convert distance to feature to integer (with error handling)
                    try:
                        if 'distance_to_feature' in eff_record and eff_record['distance_to_feature']:
                            eff_record['distance_to_feature'] = int(eff_record['distance_to_feature'])
                    except:
                        eff_record['distance_to_feature'] = -1
                        
                    effect_records.append(eff_record)
        
        # Convert lists to DataFrames
        snp_df = pd.DataFrame.from_records(snp_records)
        call_df = pd.DataFrame.from_records(call_records)
        effect_df = pd.DataFrame.from_records(effect_records)
        
        # Save to database
        snp_df.to_sql('snp', conn, if_exists='replace', index=False)
        call_df.to_sql('snp_call', conn, if_exists='replace', index=False)
        effect_df.to_sql('snp_effect', conn, if_exists='replace', index=False)
        
        conn.close()

# Rule to upload results to iRODS
rule upload_to_irods:
    input:
        vcf="050.snpeff/snps.annotated.vcf",
        db="060.database/variants.db",
        snakefile="snakefile"  # Also upload the Snakefile itself
    output:
        flag="upload_complete.flag"
    shell:
        """
        # Create the output directory in iRODS if it doesn't exist
        imkdir -p {IRODS_OUTPUT_PATH}
        
        # Upload the VCF file
        iput -f {input.vcf} {IRODS_OUTPUT_PATH}/
        
        # Upload the database
        iput -f {input.db} {IRODS_OUTPUT_PATH}/
        
        # Upload the Snakefile
        iput -f {input.snakefile} {IRODS_OUTPUT_PATH}/
        
        # Create flag file indicating upload is complete
        touch {output.flag}
        """

#Does work anymore, but kept just in case it's needed later

#rule differential_snps:
#    input:
#        db="060.database/variants.db"
#    output:
#        tsv="070.reports/high_impact_differential_snps.tsv",
#    run:
#        import sqlite3
#        import os
#        import pandas as pd
#        
#        # Create output directory if it doesn't exist
#        os.makedirs("070.reports", exist_ok=True)
#        
#        # Connect to SQLite database
#        conn = sqlite3.connect(input.db)
#    
#        
#        # Original query for high impact differential SNPs
#        query = """
#        WITH tumor_calls AS (
#            SELECT sc.snp, sc.sample as tumor_sample, sc.genotype as tumor_genotype
#            FROM snp_call sc
#            WHERE sc.sample LIKE '%_T'
#        ),
#        normal_calls AS (
#            SELECT sc.snp, sc.sample as normal_sample, sc.genotype as normal_genotype
#            FROM snp_call sc
#            WHERE sc.sample LIKE '%_N'
#        )
#        SELECT 
#            s.chrom, s.pos, s.ref, s.alt, 
#            se.gene, se.effect, se.impact, se.feature_id as transcript,
#            tc.tumor_sample, tc.tumor_genotype,
#            nc.normal_sample, nc.normal_genotype
#        FROM snp s
#        JOIN snp_effect se ON s.snp = se.snp
#        JOIN tumor_calls tc ON s.snp = tc.snp
#        JOIN normal_calls nc ON s.snp = nc.snp
#        WHERE se.impact IN ('HIGH')
#        AND (tc.tumor_genotype != nc.normal_genotype OR 
#             (tc.tumor_genotype = '1/1' AND nc.normal_genotype = '0/1'))
#        ORDER BY s.chrom, s.pos
#        """
#        
#        # Execute query and get results
#        df = pd.read_sql_query(query, conn)
#        
#        # Save to TSV file (even if empty)
#        df.to_csv(output.tsv, sep='\t', index=False)
#        
#        conn.close()

