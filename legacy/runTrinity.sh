#!/bin/bash

forwardfastq=$1 #/sequencing/projects/193/trimmed/193_YW_10_S4_L001.p1.fastq.gz
reversefastq=$2 #/sequencing/projects/193/trimmed/193_YW_10_S4_L001.p2.fastq.gz
outputdir=$3 #/sequencing/projects/193/processed/193_YW_10_S4_trinity

Trinity --seqType fq --max_memory 200G --left $forwardfastq --right $reversefastq --SS_lib_type FR --CPU 32 --output $outputdir
