//Performs basic cutadapt trimming, then alignment with bwa, then variant calling with lofreq

params.filepattern = "/usr/share/sequencing/projects/EGFR_CRISPR/raw_data/*_L001_R1_001.fastq.gz"
params.adapterfile = "/home/AD/tbleazar/usrlocalmark/NexteraPE-PE.fa"
params.outdir = "/usr/share/sequencing/projects/EGFR_CRISPR"
params.referencefolder = "/usr/share/sequencing/references/homo_sapiens/hg19/gatk_references" //and it is required that this have an indexed genome already
params.referencefile = "hg19_v0_Homo_sapiens_assembly19.fasta"
params.cpus = "32"
params.depthchart = false

forwardreads = Channel
  .fromPath(params.filepattern)

references = Channel
  .fromPath(params.referencefolder)

references.into {
  ref1
  ref2
  ref3
  ref4
}

process docutadapt {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  file(forwardfile) from forwardreads

  output:
  set ( val(sampleprefix), file("${sampleprefix}_L001_R1_001.trimmed.fastq.gz") ) into (trimmingoutput1, trimmingoutput2)
  file("${sampleprefix}.trim.out") into trimouts

  script:
  sampleprefix = (forwardfile.name).replace("_L001_R1_001.fastq.gz","")
  """
  module load CUTAdapt/latest
  cutadapt -a file:${params.adapterfile} -g file:${params.adapterfile} -o ${sampleprefix}_L001_R1_001.trimmed.fastq.gz $forwardfile -q 30,30 --minimum-length 50 --times 40 -e 0.1 --max-n 0 > ${sampleprefix}.trim.out 2> ${sampleprefix}.trim.err
  """
}

process doalignment {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '12h'
  memory '10 GB'

  input:
  set (sampleprefix, file(forwardtrimmed)) from trimmingoutput1
  file refs from ref1.first()

  output:
  set (sampleprefix, file("${sampleprefix}.unsorted.sam") ) into samfile

  script:
  """
  module load BWA/latest
  bwa mem -t ${params.cpus} -R '@RG\\tID:${sampleprefix}\\tSM:${sampleprefix}\\tPL:Illumina' ${refs}/${params.referencefile} ${forwardtrimmed} > ${sampleprefix}.unsorted.sam
  """
}

process sorttobam {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '24 GB'

  input:
  set ( sampleprefix, file(unsortedsam) ) from samfile

  output:
  set ( sampleprefix, file("${sampleprefix}.sorted.bam") ) into sortedbam

  """
  module load SAMTools/latest
  samtools sort -o ${sampleprefix}.sorted.bam -O BAM -@ ${params.cpus} ${unsortedsam}
  """
}

process indelqual {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '24 GB'

  input:
  set ( sampleprefix, file(markedbamfile) ) from sortedbam
  file refs from ref3.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.indelqual.bam") ) into (indelqualforindex, indelqualforcall)

  """
  module load LoFREQ/latest
  lofreq indelqual --dindel -f ${refs}/${params.referencefile} -o ${sampleprefix}.indelqual.bam $markedbamfile
  """
}

process samtoolsindex {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '24 GB'

  input:
  set ( sampleprefix, file(indelqualfile) ) from indelqualforindex

  output:
  set ( sampleprefix, file("${indelqualfile}.bai") ) into samindex

  """
  module load SAMTools/latest
  samtools index $indelqualfile
  """
}

forcall = indelqualforcall.join(samindex)
forcall.into {
  forcall1
  forcall2
}

process varcall {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '24 GB'

  input:
  set ( sampleprefix, file(indelqualfile), file(samindexfile) ) from forcall1
  file refs from ref4.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.lofreq.vcf") ) into finishedcalls

  """
  module load LoFREQ/latest
  lofreq call -f ${refs}/${params.referencefile} --no-default-filter --min-cov 10 -o ${sampleprefix}.lofreq.vcf --call-indels $indelqualfile
  """
}

process dodepth {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '50 GB'

  when:
  params.depthchart

  input:
  set ( sampleprefix, file(indelqualfile), file(samindexfile) ) from forcall2

  output:
  set ( sampleprefix, file("${sampleprefix}.samtools.depth") ) into samdepthout

  """
  module load SAMTools/latest
  samtools depth -aa $indelqualfile > ${sampleprefix}.samtools.depth
  """
}
