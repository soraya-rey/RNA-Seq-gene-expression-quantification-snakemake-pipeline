# RNA-Seq Differential Expression Analysis Pipeline

## Overview

This pipeline is based on the [RNA-Seq Differential Expression Analysis tutorial](https://github.com/NIGMS/RNA-Seq-Differential-Expression-Analysis/blob/master/Tutorial_2_Snakemake.ipynb) and has been adapted to be less specific to example datasets, making it applicable to a wider range of RNA-Seq data.

This pipeline performs the following steps:

1. **Read Trimming**: Removes low-quality bases and adapter sequences from raw sequencing reads.
2. **Read Quality Control (QC)**: Assesses the quality of trimmed reads using FastQC.
3. **Read Mapping**: Aligns the cleaned reads to a reference transcriptome using Salmon.
4. **Counting Mapped Reads**: Quantifies gene expression by counting the number of reads mapped to each gene.

This pipeline does not yet include differential expression analysis but is designed to handle various RNA-Seq datasets.

## Required Input Data

To run this pipeline, the following input data are required:

1. **Paired-end sequencing data**: 
   - FastQ files containing the sequencing reads for each sample. Each sample should have two FastQ files, typically named with suffixes `_1.fastq` and `_2.fastq` (e.g., `sample_1.fastq`, `sample_2.fastq`).

2. **Adapter sequences**: 
   - Adapter files in FASTA format for trimming the reads. These files should contain the adapter sequences that need to be removed from the reads.

3. **Decoy sequences**: 
   - A FASTA file containing decoy sequences for indexing and mapping. This file can help improve the accuracy of alignment by accounting for potential off-target mappings.

4. **Reference genome**: 
   - A FASTA file of the reference genome or transcriptome to which the reads will be aligned.


## Steps and Packages Used

The pipeline consists of several steps, each utilizing specific bioinformatics tools and packages:

1. **Read Trimming**:
   - **Tool**: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
   - **Purpose**: Removes adapter sequences and low-quality bases from the FastQ files. It enhances the quality of the data for subsequent analysis.

2. **Quality Control**:
   - **Tool**: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
   - **Purpose**: Generates a report that provides an overview of the quality of the sequencing data. It helps identify potential issues with the FastQ files.

3. **Indexing**:
   - **Tool**: [Salmon](https://salmon.readthedocs.io/en/latest/)
   - **Purpose**: Creates a transcriptome index from the reference transcriptome FASTA file. This index is used for efficient mapping of reads to transcripts.

4. **Quantification**:
   - **Tool**: [Salmon](https://salmon.readthedocs.io/en/latest/)
   - **Purpose**: Quantifies transcript abundance from the trimmed FastQ files using the created index. The output is an expression matrix that can be used for differential expression analysis.

5. **MultiQC**:
   - **Tool**: [MultiQC](https://multiqc.info/)
   - **Purpose**: Compiles the output from FastQC and presents it in a single report. This helps in assessing the overall quality of the sequencing runs.



## How to Use the Pipeline

To run the RNA-Seq differential expression analysis pipeline, navigate to the directory containing your `Snakefile` and execute the following command:

```bash
snakemake all --use-conda -p
```