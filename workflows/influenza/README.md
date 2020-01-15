# Influenza Classifier

## Scope

This pipeline addresses the needs of HGR process in selecting influenza serotypes from co-infection of parental viral strains.
It analyses the sequencing reads of each sample against a database of parental genomes and reports an abundance table by each of the 8 influenza genes.

## Expected Inputs

The pipeline will expect:

- A FASTA file containing the parental genome segments, with a header formatted as ">GENE-STRAIN"
- A FASTQ file containing the reads resulting from the sequencing

## Running the pipeline

In order to launch the Nextflow script, the software needs to be loaded as environmental module in Slurm, with the following command:

```
module load NextFlow/latest
```

We assume that our *core* Github pipeline has been cloned under the user home, in the folder CODE and therefore is available at

```
$HOME/CODE/core/
```

The script expects 4 mandatory arguments:

1. the FASTA file of parental genomes
2. the path to a folder containing all FASTQ files (without subfolders) of the samples
3. a name to be assigned to the BLAST database
4. the path to an existing folder where the results will be copied

The workflow can then be launched in the following way

```
nextflow \
-c ~/CODE/core/workflows/nextflow.config \
run ~/CODE/core/workflows/influenza/influenza_classifier.nf \
-with-trace \
-with-timeline nf_influenza_timeline.htm \
-with-report nf_influenza_report.htm \
--reads /path/to/project/raw_data \
--origin path/to/project/raw_data/parental_genomes.fasta \
--db "parentalDB" \
--output_dir /path/to/project/analysis/run_identifier
```

## Expected Outputs

The workflow will create a subfolder for each sample, containing:

- A FASTA file of the original reads, converted from the FASTQ: this can be used to further analyse specific reads of a sample matching a parental gene of interest
- A blast results file, formatted in custom format for the analysis

A final HTML report, as well as pipeline monitoring reports, will be copied on the parent folder.


## Methods

The workflow creates a blast database out of the FASTA file of parental genomes and uses FASTQ reads as query each run through *blastn*.
*blastn* is configured to report only the top result, with the following command:

```
blastn \
-query ${sampleId}.fa \
-db ${dbName} \
-max_target_seqs 1 \
-num_threads ${task.cpus} \
-outfmt '6 qseqid sseqid sgi qstart qend sstart send pident mismatch nident evalue' \
| sort -k 1,1 -k11,11g > "${sampleId}_blast_results.txt"
```
