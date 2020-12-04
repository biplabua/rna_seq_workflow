
import pandas as pd
configfile: "config.yaml"
samples_df = pd.read_csv('meta_data.csv').set_index("sample", drop=False)
SAMPLES = list(samples_df['sample'])

# download yeast rna-seq data from Schurch et al, 2016 study
rule all:
    input:
        expand("raw_data/{sample}.fq.gz", sample=SAMPLES),
        #"Saccharomyces_cerevisiae.R64-1-1.101.gtf.gz",
        "rna_coding.fasta.gz",
        "rna_coding_index",
        expand("out_salmon/{sample}", sample=SAMPLES),
        "log2fc.csv"

        #"quants/{sample}_quant"

# rule to download each individual file specified in samples_df
rule download_reads:
    output: "raw_data/{sample}.fq.gz" 
    params:
        # dynamically grab the download link from the "dl_link" column in the samples data frame
        download_link = lambda wildcards: samples_df.loc[wildcards.sample, "id"]
    shell:
        "fastq-dump -Z --gzip {params.download_link} > {output}"

rule download_transcriptome:
    output:
        "rna_coding.fasta.gz"
    params:
        transcriptome=config["transcript_sequence"]
    shell:
        "wget {params.transcriptome} -O {output}"

rule make_salmon_index:
    input:
        "rna_coding.fasta.gz"
    output:
        directory("rna_coding_index")
    conda:
        "environment.yaml"
    shell:
        "salmon index -t {input} -i {output}"

rule salmon_count:
    input:
        fq="raw_data/{sample}.fq.gz",
        index="rna_coding_index"
    output:
        directory("out_salmon/{sample}")
    conda:
        "environment.yaml"
    shell:
        "salmon quant -i {input.index} -l A -r {input.fq} --validateMappings -o {output}"
rule dge_deseq2:
    input:
        meta_data="meta_data.csv",
        file_path=expand("out_salmon/{sample}", sample=SAMPLES)
    output:
        out_file="log2fc.csv"
    conda:
        "dge.yaml"
    script:
        "deseq2.R"
