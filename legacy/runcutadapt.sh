#!/bin/bash
adapter=$1 #/sequencing/projects/193/TruSeq-adapters-recommended.txt
forwardfile=$2 #/sequencing/nextseq/processed/180817/fastq/raw/193_YW_10_S4_L001_R1_001.fastq.gz
reversefile=$3 #/sequencing/nextseq/processed/180817/fastq/raw/193_YW_10_S4_L001_R2_001.fastq.gz
outfor=$4 #/home/AD/tbleazar/test/trimmed.for.fastq.gz
outrev=$5 #/home/AD/tbleazar/test/trimmed.rev.fastq.gz

cutadapt -a file:$adapter -A file:$adapter -g file:$adapter -G file:$adapter -o $outfor -p $outrev $forwardfile $reversefile -q 30,30 --minimum-length 50 --times 40
