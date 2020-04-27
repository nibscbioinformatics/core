//This is a nextflow script to be run with
//nextflow run trimonly.nf
//Performs basic cutadapt trimming

params.reads = null
params.adapterfile = "/usr/share/sequencing/references/adapters/NexteraPE-PE.fa"
params.trimoutdir = null

Channel
    .fromFilePairs("$params.reads/*_{R1,R2}*.fastq.gz")
    .ifEmpty { error "Cannot find any reads matching ${params.reads}"}
    .set { readpairs }

process docutadapt {
  publishDir "$params.trimoutdir", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  set ( sampleprefix, file(samples) ) from readpairs

  output:
  set ( sampleprefix, file("${sampleprefix}_L001_R1_001.trimmed.fastq.gz"), file("${sampleprefix}_L001_R2_001.trimmed.fastq.gz") ) into trimmingoutput
  file("${sampleprefix}.trim.out") into trimouts

  script:
  """
  module load CUTAdapt/latest
  cutadapt -a file:${params.adapterfile} -A file:${params.adapterfile} -g file:${params.adapterfile} -G file:${params.adapterfile} -o ${sampleprefix}_L001_R1_001.trimmed.fastq.gz -p ${sampleprefix}_L001_R2_001.trimmed.fastq.gz ${samples[0]} ${samples[1]} -q 30,30 --minimum-length 50 --times 40 -e 0.1 --max-n 0 > ${sampleprefix}.trim.out 2> ${sampleprefix}.trim.err
  """
}

process dotrimlog {
  publishDir "$params.trimoutdir", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  file "logdir/*" from trimouts.toSortedList()

  output:
  file("trimming-summary.csv") into trimlogend

  script:
  """
  module load anaconda/Py2/python2
  python $HOME/CODE/core/utilities/logger.py logdir trimming-summary.csv cutadapt
  """
}

