# Metagenomics Analyses

The code included in this folder performs different metagenomics analyses, as described in the paragraphs below.

## Taxonomic Profiling with Kraken

### Nextflow Pipeline





## Taxonomic and Functional Profiling with Metaphlan and Humann2

The following paragraphs will briefly explain how to run Humann2, using the scripts we have prepared.


### Nextflow Pipeline

Our Nextflow pipeline launches Humann2 on all the fastq read pairs present in a specified input folder: the pipeline performs a three tiered analysis as described in Humann2 paper.

Expected inputs:

- an input folder where all read pairs for all samples are located
- an output folder where new sample folders will be created and joined tables will be copied

The outputs will include:

- Metaphlan taxonomic profiling for each sample
- Abundance tables for genes and pathways for each sample
- Joined gene and pathways abundance tables
- Joined gene and pathways abundance tables, re-normalised in cpm

In order to run the Nextflow pipeline, the user will submit a batch job running the master Nextflow script, as follows:

```
sbatch \
-D `pwd` \
-o "`pwd`/nextflow_metafunctional_job_%j.out" \
~/CODE/core/workflows/metagenomics/run_nextflow_humann2.job \
/path/to/project/inputs_folder \
/path/to/project/outputs_folder
```



### Single Sample - Slurm submission

If you are running on NIBSC cluster, and you wish to analyse one or few samples only, you can also submit a batch jobs through our Slurm scheduler.

**NOTE:**
Currently *Humann2* does not take into account of paired end reads information: **before running Humann2, you need to concatenate the 2 fastq files into one**, and use the concatenated file as input of this script.

In order to run this script you need 3 elements:

- the concatenated FASTQ input file
- the folder where you want to write the project data
- the name of the folder where you will write the sample results

For the first 2 elements, **always use the absolute path**, while the last element (i.e. name of the output) you should **only indicate the name** of the new folder to be created, under the project data folder.

You can now launch a job running the analysis by entering in the command line:

```
sbatch \
-D `pwd` \
-o "`pwd`/humann2_job_%j.out" \
submit_humann2_analysis.job \
/path/to/your/input_file_concatenated.fastq \
/path/to/the/chosen/project_folder \
sample_name

```

If you have a small number of samples (concatenated files) in the same folder, you can launch them all with a simple bash loop

```
for sample in `ls /path/to/folder/*_concatenated.fastq`
do
sampleFile=`echo $sample | awk -F/ '{print $NF}'`
sampleName=${sampleFile%.fastq}

echo "#### launching analysis on sample ${sampleName}"

sbatch \
-D `pwd` \
-o "`pwd`/humann2_job_%j.out" \
submit_humann2_analysis.job \
${sample} \
/path/to/the/chosen/project_folder \
${sampleName}

done
```
