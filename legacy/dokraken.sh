#!/bin/bash

forward=$1 #/home/AD/tbleazar/165/input/180430/fastq/30/165-ChildC-D1day37-17-18_S3_L001.s1.fastq.gz
reverse=$2 #the reverse fastq.gz
outdir=$3 #/home/AD/tbleazar/165/testkraken

export PATH="/usr/local/share/Bracken-master:$PATH"

mkdir -p $outdir
cd $outdir

echo Launching Kraken for input fastq $forward $reverse
/usr/local/bin/kraken --threads 4 --paired --preload --db /raid/kraken/standard_kraken_db/ --classified-out $outdir/kraken.classified --unclassified-out $outdir/kraken.unclassified --out $outdir/kraken.output $forward $reverse
/usr/local/bin/kraken-mpa-report --db /raid/kraken/standard_kraken_db/ $outdir/kraken.output >  $outdir/kraken.mpa.report
kraken-report --db /raid/kraken/standard_kraken_db/ $outdir/kraken.output > $outdir/kraken.report
bracken -d /raid/kraken/standard_kraken_db -i $outdir/kraken.report -o $outdir/bracken.report -r 140 -t 5
