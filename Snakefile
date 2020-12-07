import pandas as pd
configfile: "config.yaml"
samples_df = pd.read_csv('meta_data.csv').set_index("sample", drop=False)
SAMPLES = list(samples_df['sample'])

myoutput = list()

if config["paired"] == True:
    myoutput.append(expand("raw_data/{sample}_1.fastq", sample=SAMPLES)),
    myoutput.append(expand("raw_data/{sample}_2.fastq", sample=SAMPLES))
if config["paired"]== False:
    myoutput.append(expand("raw_data/{sample}.fq.gz", sample=SAMPLES))
if config["paired"] == True:
    ruleorder: salmon_count_pe > salmon_count_se
if config["paired"]== False:
    ruleorder: salmon_count_se > salmon_count_pe
rule all:
    input:
        myoutput,
        "rna_coding.fasta.gz",
        "rna_coding_index",
        expand("out_salmon/{sample}", sample=SAMPLES),
        "results/log2fc.csv",
        "results/volcano_plot.png"
# Rule to download fastq file from SRA
rule download_reads:
    output: "raw_data/{sample}.fq.gz" 
    params:
        # dynamically grab the SRA ID from the "id" column in the samples_df data frame
        download_link = lambda wildcards: samples_df.loc[wildcards.sample, "id"]
    conda:
        "environment.yaml"
    shell:
        "fastq-dump -Z --gzip {params.download_link} > {output}"

rule download_reads_pe:
    output: "raw_data/{sample}_1.fastq", "raw_data/{sample}_2.fastq"
    params:
        # dynamically grab the SRA ID from the "id" column in the samples_df data frame
        download_link = lambda wildcards: samples_df.loc[wildcards.sample, "id"],
        output_dir = directory("raw_data")
    conda:
        "environment.yaml"
    shell:
        "fastq-dump --split-3 {params.download_link} -O {params.output_dir}"
# Rule to download transcritome sequence fro the link in config.yaml file
rule download_transcriptome:
    output:
        "rna_coding.fasta.gz"
    params:
        transcriptome=config["transcript_sequence"]
    shell:
        "wget {params.transcriptome} -O {output}"
# Rule to make salmon index of transcript sequence
rule make_salmon_index:
    input:
        "rna_coding.fasta.gz"
    output:
        directory("rna_coding_index")
    conda:
        "environment.yaml"
    shell:
        "salmon index -t {input} -i {output}"
# Rule to get count from fastq file using Salmon 
rule salmon_count_se:
    input:
        fq="raw_data/{sample}.fq.gz",
        index="rna_coding_index"
    output:
        directory("out_salmon/{sample}")
    conda:
        "environment.yaml"
    shell:
        "salmon quant -i {input.index} -l A -r {input.fq} --validateMappings -o {output}"

rule salmon_count_pe:
    input:
        fq1="raw_data/{sample}_1.fastq",
        fq2="raw_data/{sample}_2.fastq",
        index="rna_coding_index"
    output:
        directory("out_salmon/{sample}")
    conda:
        "environment.yaml"
    shell:
        "salmon quant -i {input.index} -l A -1 {input.fq1} -2 {input.fq2} --validateMappings -o {output}"
# Rule for differential expression between two sample
rule dge_deseq2:
    input:
        meta_data="meta_data.csv",
        file_path=expand("out_salmon/{sample}", sample=SAMPLES)
    output:
        out_file="results/log2fc.csv"
    conda:
        "dge.yaml"
    script:
        "scripts/deseq2.R"
rule volcano_plot:
    input: 
        log2fc="results/log2fc.csv"
    output: 
        volcano_plot="results/volcano_plot.png"
    conda:
        "volcano_plot.yaml"
    script:
        "scripts/volcano_plot.R"







