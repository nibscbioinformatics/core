# Small Genome Alignment and Typing Pipeline

## Summary

This is a Nextflow pipeline developed by the NIBSC NGS Core Bioinformatics team to perform alignment to a reference and call variants for input sequencing data. It requires FASTQ format input data and a genome reference, and outputs variant calls and coverage data.

## Software

This pipeline uses the tool CutAdapt to perform trimming to remove adapter sequences given an adapter FASTA reference file. This tool is also set to remove low quality bases from the ends of reads. The pipeline uses BWA for alignment to the reference. The tool LoFreq is then used to call variants, allowing for calling at low minor allele read frequency. We generate CSV output of these variant calls for each sample, with information of the number of reads supporting calls.

## Required inputs

In order to launch the analysis, you will need the following input files

- The set of gzipped fastq files to analyse
- The folder and name of the adapter file
- The genome reference, which must be indexed before use

All paths should be specified as **absolute paths**.

## Running the pipeline

In order to run the pipeline on an ordinary dataset, you should submit a job which will execute a Nextflow run command.

Within an interactive session on a HPC compute node, this can be completed with the command:

module load NextFlow/latest
nextflow run ~/CODE/core/workflows/typing/cutlo.nf --filepattern '<your data folder>/*{_R1_001,_R2_001}.fastq.gz' --outdir <desired output folder> --referencefolder <folder in which reference genome and index is generated> --referencefile <genome fasta> -c ~/CODE/core/workflows/typing/nextflow.config



