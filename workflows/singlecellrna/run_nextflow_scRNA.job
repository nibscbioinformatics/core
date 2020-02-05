#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5g
#SBATCH -t 1-00:00:00
#SBATCH --job-name NFmaster

metadata=$1
reference=$2
output=$3

module load NextFlow/19.10.0

nextflow \
-c ~/CODE/core/workflows/singlecellrna/scRNA.config \
run ~/CODE/core/workflows/singlecellrna/scrna_short.nf \
-with-trace \
-with-timeline nf_scRNA_timeline.htm \
-with-report nf_scRNA_report.htm \
--metadata ${metadata} \
--reference ${reference} \
--output_dir ${output}
