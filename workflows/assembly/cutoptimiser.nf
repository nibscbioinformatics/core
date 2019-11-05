//This is a nextflow script to trim and then use velvetoptimiser to assemble

params.filepattern = "/usr/share/sequencing/miseq/output/161007_M01745_0131_000000000-ATN6R/Data/Intensities/BaseCalls/049*{_L001_R1_001,_L001_R2_001}.fastq.gz"
params.adapterfile = "/usr/share/sequencing/references/adapters//NexteraPE-PE.fa"
params.outdir = "/home/AD/tbleazar/test"
params.cpus = "32"

Channel
  .fromFilePairs(params.filepattern)
  .set { readpairs }

process docutadapt {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  set ( sampleprefix, file(samples) ) from readpairs

  output:
  set ( sampleprefix, file("${sampleprefix}_L001_R1_001.trimmed.fastq.gz"), file("${sampleprefix}_L001_R2_001.trimmed.fastq.gz") ) into (trimmingoutput1, trimmingoutput2)
  file("${sampleprefix}.trim.out") into trimouts

  script:
  """
  module load CUTAdapt/latest
  cutadapt -a file:${params.adapterfile} -A file:${params.adapterfile} -g file:${params.adapterfile} -G file:${params.adapterfile} -o ${sampleprefix}_L001_R1_001.trimmed.fastq.gz -p ${sampleprefix}_L001_R2_001.trimmed.fastq.gz ${samples[0]} ${samples[1]} -q 30,30 --minimum-length 50 --times 40 -e 0.1 --max-n 0 > ${sampleprefix}.trim.out 2> ${sampleprefix}.trim.err
  """
}

process dotrimlog {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  file "logdir/*" from trimouts.toSortedList()

  output:
  file("read-summary.csv") into trimlogend

  script:
  """
  module load anaconda/Py2/python2
  python /home/AD/tbleazar/pipelines/trimlogger.py logdir read-summary.csv
  """
}

process dovelvet {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 32
  queue 'WORK'
  time '24h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(forwardfile), file(reversefile) ) from trimmingoutput1

  output:
  set ( sampleprefix, file("${sampleprefix}_velvet"), file("${sampleprefix}_velvetlog.txt") ) into velvetoutput

  """
  module load velvetoptimiser/2.2.6
  VelvetOptimiser.pl -s 25 -e 145 -t ${params.cpus} -p vorun -m "-0.001" -f '-fastq.gz -shortPaired -separate $forwardfile $reversefile'
  mv `ls | grep vorun_data` ${sampleprefix}_velvet
  mv vorun_logfile.txt ${sampleprefix}_velvetlog.txt
  """
}

process dospades {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 32
  queue 'WORK'
  time '24h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(forwardfile), file(reversefile) ) from trimmingoutput2

  output:
  set ( sampleprefix, file("${sampleprefix}_spades") ) into spadesoutput

  """
  module load BWA/latest
  python /home/AD/tbleazar/spades/SPAdes-3.13.1-Linux/bin/spades.py -o ${sampleprefix}_spades -1 $forwardfile -2 $reversefile -t ${params.cpus} -m 240 --cov-cutoff 10.0
  """
}

