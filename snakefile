import glob
import os
import subprocess

# Capture available cores dynamically
def get_threads():
    result = subprocess.run(
        "lscpu | grep '^Processeur(s)' | awk '{print $2-1}'", 
        shell=True, 
        capture_output=True, 
        text=True
    )
    return int(result.stdout.strip()) if result.stdout.strip() else 1

THREADS = get_threads()  # Using one fewer than total cores to avoid overloading
print("Threads available:", THREADS)


# Automatically extract sample names from raw data directory
raw_data_dir = "data/raw_fastq"
SAMPLES = [os.path.basename(f).split("_1.fastq")[0] for f in glob.glob(f"{raw_data_dir}/*_1.fastq")]


print("Samples detected :")
print(SAMPLES)

adapter_file_path = "/home/utilisateur/Bureau/snakemake_differential_gene_expression_pipeline/data/adapters/TruSeq3-PE.fa"

rule all:
    input:
        expand("data/trimmed/{sample}_trimmed_1.fastq", sample=SAMPLES),
        expand("data/fastqc/{sample}_trimmed_1_fastqc.html", sample=SAMPLES),
        "multiqc_report.html",
        "data/reference/transcriptome_index/refseq.bin",
        expand("data/quants/{sample}_quant/quant.sf", sample=SAMPLES)

rule trimmomatic_pe_fq:
    input:
        r1="data/raw_fastq/{sample}_1.fastq",
        r2="data/raw_fastq/{sample}_2.fastq",
        adapter= adapter_file_path
    output:
        r1="data/trimmed/{sample}_trimmed_1.fastq",
        r2="data/trimmed/{sample}_trimmed_2.fastq",
        r1_unpaired="data/trimmed/{sample}_trimmed_1.unpaired.fastq",
        r2_unpaired="data/trimmed/{sample}_trimmed_2.unpaired.fastq"
    params:
        extra="LEADING:3 TRAILING:3 MINLEN:36",
    threads: THREADS
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "mkdir -p data/trimmed"
        "trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} "
        "ILLUMINACLIP:{input.adapter}:2:30:10:2:keepBothReads {params.extra}"

rule fastqc_trimmed:
    input:
        "data/trimmed/{sample}_trimmed_1.fastq",
        "data/trimmed/{sample}_trimmed_2.fastq"
    output:
        "data/fastqc/{sample}_trimmed_1_fastqc.html",
        "data/fastqc/{sample}_trimmed_2_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir="data/fastqc"
    threads: THREADS
    shell:
        "mkdir -p {params.outdir}"
        "fastqc -o {params.outdir} -t {threads} {input}"

rule multiqc_trimmed:
    input:
        expand("data/fastqc/{sample}_trimmed_1_fastqc.html", sample=SAMPLES),
        expand("data/fastqc/{sample}_trimmed_2_fastqc.html", sample=SAMPLES)
    output:
        "multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc data/fastqc/ -o . -f"

rule salmon_index:
    input:
        "data/reference/M_chelonae_transcripts.fasta"
    output:
        directory("data/reference/transcriptome_index")
    threads: THREADS
    conda:
        "envs/salmon.yaml"
    params:
        decoys="data/reference/decoys.txt"
    shell:
        "salmon index -t {input} -i {output} -d {params.decoys} --threads {threads}"

rule salmon_quant_reads:
    input:
        r1="data/trimmed/{sample}_trimmed_1.fastq",
        r2="data/trimmed/{sample}_trimmed_2.fastq",
        index="data/reference/transcriptome_index"
    output:
        "data/quants/{sample}_quant/quant.sf"
    params:
        libtype="A",
        extra="",
        q1="data/quants/{sample}_quant/"
    threads: THREADS
    conda:
        "envs/salmon.yaml"
    shell:
        " mkdir -p data/reference/transcriptome_index"
        "salmon quant -i {input.index} -l SR -r {input.r1} -p {threads} --validateMappings -o {params.q1}"
