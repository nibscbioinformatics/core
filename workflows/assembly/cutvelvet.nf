//This is a nextflow script to trim and then use velveth and velvetg with hardcoded parameters based on an experiment to assemble

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

process dovelveth {
  cpus 32
  queue 'WORK'
  time '24h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(forwardfile), file(reversefile) ) from trimmingoutput1

  output:
  set ( sampleprefix, file("${sampleprefix}_assembly") ) into (velvethoutput1, velvethoutput2)

  """
  module load velvetoptimiser/2.2.6
  velveth ${sampleprefix}_assembly 75 -fastq.gz -shortPaired -separate $forwardfile $reversefile
  """
}

process dovelvetglong {
  publishDir "$params.outdir/analysis", mode: "copy"
  queue 'WORK'
  time '24h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(assemblydir) ) from velvethoutput1

  output:
  set ( sampleprefix, file("${sampleprefix}_long/contigs.fa") ) into longoutput

  """
  module load velvetoptimiser/2.2.6
  mkdir ${sampleprefix}_long
  cp ${assemblydir}/* ${sampleprefix}_long
  velvetg ${sampleprefix}_long -exp_cov 15 -cov_cutoff 0.5 -ins_length 555 -min_contig_lgth 5000
  """
}

process dovelvetgkb {
  publishDir "$params.outdir/analysis", mode: "copy"
  queue 'WORK'
  time '24h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(assemblydir) ) from velvethoutput2

  output:
  set ( sampleprefix, file("${sampleprefix}_assembled") ) into kboutput

  """
  module load velvetoptimiser/2.2.6
  mkdir ${sampleprefix}_assembled
  cp ${assemblydir}/* ${sampleprefix}_assembled
  velvetg ${sampleprefix}_assembled -exp_cov 15 -cov_cutoff 0.5 -ins_length 555 -min_contig_lgth 1000
  """
}

