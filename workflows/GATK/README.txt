WHAT IS IT?

This is a pipeline based on the GATK workflow to detect germline CNV variants in a WGS or WES data.
It is constituted by 2 parts, and so it is divided in 2 different nextflow files: GATK_CNV_COHORT.nf and GATK_CNV_CASE.nf. The first file refers to the construction of a model to be used as a reference (a model is optimally made out from as many 'normal' samples as possible to statistically determine what constitutes as a CNV variant or not in the tested sample), and the second one to the detection itself of CNV variants in a sample. A diagram of these processes is found in the schematic.png file.


WHAT ARE THE PROCESSES?

Summarized information from each process is commented inside the pipeline. To get more insight, you can find more infrmation here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants


HOW TO RUN?

It is needed a GATK container to run the pipeline, and this can be download by using Singularity with the following command:
singularity build --sandbox gatk_4.1.3.0/ docker://broadinstitute/gatk
Every additional required file is mentioned where to get it in the SOP XXXX.

To run each pipeline is recommended to be attached to a config file which determines input/output files path and other metadata. Here it is described what processes are done through te container.
