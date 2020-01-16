This documents explains, with examples, how to construct a file to demultiplex different sample data.


Our FASTQ files should be names in this format:
	Sample_S1_L00X_R1_001.fastq.gz
	
	Example:
	673168_13_S6_L002_I1_001.fastq.gz


In case we have multiple lanes sequenced for the same sample, one line of our document should look like this:
	cellranger count --id=name_of_sample --sample=fastq_id_1,fastq_id_2 --fastq=/path/to/FASTQ/files_lane_1,/path/to/FASTQ/files_lane_2 --transcriptome=/path/to/reference_genome

	Example:
	cellranger count --id=bone_marrow --sample=673168_13,676445_13 --fastq=~/MSC/bsg-ftp.well.ox.ac.uk/190416_A00711_0029_BHJ5GMDMXX/FASTQ/673168_13,~/MSC/bsg-ftp.well.ox.ac.uk/190426_A00711_0030_AHJMGJDMXX/FASTQ/676445_13/ --transcriptome=~/cell_ranger/human_reference/refdata-cellranger-GRCh38-3.0.0.tar/refdata-cellranger-GRCh38-3.0.0/


In case we have 1 lane sequenced for 1 sample, one line of our document should look like this:
	cellranger count --id=name_of_sample --sample=fastq_id --fastq=/path/to/FASTQ/files_lane --transcriptome=/path/to/reference_genome

	Example:
	cellranger count --id=bone_marrow_mb --sample=673168_14 --fastq=~/MSC/bsg-ftp.well.ox.ac.uk/190416_A00711_0029_BHJ5GMDMXX/FASTQ/673168_14/ --transcriptome=~/cell_ranger/human_reference/refdata-cellranger-GRCh38-3.0.0.tar/refdata-cellranger-GRCh38-3.0.0/
	
	
Notes:
	- fastq_id_1 and fastq_id_2 should be the Sample in our original FASTQ file, in our case 673168_13
	- /path/to/FASTQ/files_lane folder should have FASTQ files in the correct format
	- Example of file available, with both approaches, named 'MSC_count.sh'
	- More info: See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct
	