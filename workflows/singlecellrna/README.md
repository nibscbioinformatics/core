# scRNA Sequencing Pipeline

## Software

This Nexflow pipeline is using **10X CellRanger** to perform the UMI handling and counting of the transcripts, and then uses **Seurat** to perform the data aggregation, normalisation and further downstream analyses.

The pipeline can be run anywhere using a Singularity container setup with all necessary tools. The Singularity container is specified in the config file.

Note: the config file also makes sure the container mounts the path */usr/share* where the data are located in our system. If you run it somewhere else you might want to change this to an appropriate path.

## Required inputs

In order to launch the analysis, you will need the following input files

- A metadata file, containing three columns: the sample name, the name of the file (i.e. the identifier of the FASTQ file without the lane information), and the path to the FASQT file
- 10X CellRanger reference: the path to downloaded and unzipped references from 10X you would like to use
- Path to an output folder, where the results will be stored

All paths should be specified as **absolute paths**.


## Testing the pipeline

We assume you have cloned the repository in your home folder, under the path ```$HOME/CODE``` and therefore the path to the code will be ```$HOME/CODE/core```.

**if** *(and only if)* you want to test the pipeline and have a small sample set, the workflow can be run by opening first an interactive session on a node:

```
srun --mem=3g --pty /bin/bash
```

and then using the following command

```
module load NextFlow/latest

nextflow \
-c ~/CODE/core/workflows/singlecellrna/scRNA.config \
run ~/CODE/core/workflows/singlecellrna/scrna_short.nf \
-with-trace \
-with-timeline nf_scRNA_timeline.htm \
-with-report nf_scRNA_report.htm \
--metadata /path/to/metadata/file.txt \
--reference /path/to/GRCh38/cellranger/refdata-cellranger-GRCh38-3.0.0 \
--output_dir /path/to/your/destination
```
**Note**:
the **singularity image** has been created and stored into our system, and the path is specified in the config file.
You can create the image by using the image definition available under this folder **scRNA_image.def**.
The config file also specifies a mount point as */usr/share*: you will need to change that, if you want to mount other locations in your cluster.

## Running the pipeline

**In all other cases**, in order to run the pipeline on a ordinary dataset, you should instead submit a job which will execute the above command.

The job script is also available in this folder and can be submitted by passing the same arguments required by nextflow, as follows:

```
sbatch \
-D `pwd` \
-o "`pwd`/nextflow_job_%j.out" \
~/CODE/core/workflows/singlecellrna/run_nextflow_scRNA.job \
/path/to/metadata/file.txt \
/path/to/GRCh38/cellranger/refdata-cellranger-GRCh38-3.0.0 \
/path/to/your/destination
```

You can also adjust the requested wall time (currently maxing to one day by default) by adding the option `-t 5:00:00` (in this example would change it to 5 hours).
