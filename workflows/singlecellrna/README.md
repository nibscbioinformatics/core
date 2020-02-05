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


## Running the pipeline

We assume you have cloned the repository in your home folder, under the path ```$HOME/CODE``` and therefore the path to the code will be ```$HOME/CODE/core```.

The pipeline can be run with the following command

```
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
Note that the **singularity image** has been created and stored into our system, and the path is specified in the config file.
You can create the image by using the image definition available under this folder **scRNA_image.def**. 
