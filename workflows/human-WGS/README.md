# Human WGS Alignment and GATK Variant Calling Workflow

## Summary

This is a Nextflow pipeline developed by the NIBSC NGS Core Bioinformatics team to follow GATK "best practice" in processing human whole genome sequencing data. It performs alignment and then applies germline and somatic variant calling methods.

## Software

This pipeline performs alignment to the reference human genome hg19 using BWA. It processes alignments using SAMTools and GATK utilities. It uses GATK Base Recalibration to recalibrate base quality scores using reference indel locations. The default reference files for use with this pipeline are located on the HPC for hg19 and will be used automatically unless specified. Germline variant calling is performed using GATK HaplotypeCaller. Somatic unpaired sample variant calling is performed using GATK Mutect2 which references a panel of normals and germline resource. We use GATK VariantFiltration on the respective variant outputs to provide a high confidence call set. Possible effects of variants are predicted using the tool SNPEff.

## Required inputs

In order to launch the analysis, you will need the following input files

- The set of gzipped fastq files of human WGS samples
- If a reference other than hg19 is to be used, this needs to be prepared together with appropriate panel files

All paths should be specified as **absolute paths**.

## Running the pipeline

In order to run the pipeline on a ordinary dataset, you should submit a job which will execute a Nextflow run command.

Within an interactive session on a HPC compute node, the command can be run with:

module load NextFlow/latest

nextflow run ~/CODE/core/workflows/human-WGS/gatk.nf --outdir /folder/for/output --inputfolder /folder/with/input/files

This assumes that there are fastq.gz files within the input folder. Alternatively, it is possible to run this pipeline using BAM input data. For this, give the additional flag --frombam when running.
