#!/opt/software/conda2/envs/NextFlow/bin/nextflow

//This is a nextflow script to perform GATK4 best practice processing and then apply LoFreq for trimmed WGS fastq data

params.referencefolder = "/home/AD/tbleazar/228redo/gatk_references/" //we require there to already be indexed references here
params.outdir = "/home/AD/tbleazar/228redo/test"
params.filepattern = "/home/AD/tbleazar/228redo/testinput/228_HT1080_test{_R1_001,_R2_001}.fastq.gz"
//params.filepattern = "/usr/share/sequencing/projects/228/trimmed/228_HT1080_WT_P15_S1{_R1_001,_R2_001}.fastq.gz"
params.cpus = "32"

Channel
  .fromFilePairs(params.filepattern)
  .set { readpairs }

references = Channel
  .fromPath(params.referencefolder)

references.into {
  ref1
  ref2
  ref3
}

process doalignment {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '100 GB'

  input:
  set ( sampleprefix, file(samples) ) from readpairs
  file refs from ref1.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.unsorted.sam") ) into samfile

  """
  module load BWA/latest
  bwa mem -t ${params.cpus} -M -R '@RG\\tID:${sampleprefix}\\tSM:${sampleprefix}\\tPL:Illumina' ${refs}/hg19_v0_Homo_sapiens_assembly19.fasta ${samples[0]} ${samples[1]} > ${sampleprefix}.unsorted.sam
  """
}

process sorttobam {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(unsortedsam) ) from samfile

  output:
  set ( sampleprefix, file("${sampleprefix}.sorted.bam") ) into sortedbam

  """
  module load SAMTools/latest
  samtools sort -o ${sampleprefix}.sorted.bam -O BAM -@ ${params.cpus} ${unsortedsam}
  """
}

process markduplicates {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(sortedbamfile) ) from sortedbam

  output:
  set ( sampleprefix, file("${sampleprefix}.marked.bam") ) into (markedbamfortable, markedbamforapply)

  """
  module load GATK/4.1.3.0
  gatk MarkDuplicates -I $sortedbamfile -M ${sampleprefix}.metrics.txt -O ${sampleprefix}.marked.bam
  """
}

process baserecalibrationtable {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(markedbamfile) ) from markedbamfortable
  file refs from ref2.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.recal_data.table") ) into recaltable

  """
  module load GATK/4.1.3.0
  gatk BaseRecalibrator -I $markedbamfile --known-sites ${refs}/b37_dbsnp_138.b37.vcf --known-sites ${refs}/b37_Mills_and_1000G_gold_standard.indels.b37.vcf -O ${sampleprefix}.recal_data.table -R ${refs}/hg19_v0_Homo_sapiens_assembly19.fasta
  """
}

forrecal = recaltable.join(markedbamforapply)

process applybaserecalibration {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(recalibrationtable), file(markedbamfile) ) from forrecal

  output:
  set ( sampleprefix, file("${sampleprefix}.bqsr.bam") ) into (recalibratedforindex, recalibratedforcaller)

  """
  module load GATK/4.1.3.0
  gatk ApplyBQSR -I $markedbamfile -bqsr $recalibrationtable -O ${sampleprefix}.bqsr.bam
  """
}

process indexrecalibrated {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(bqsrfile) ) from recalibratedforindex

  output:
  set ( sampleprefix, file("${bqsrfile}.bai") ) into indexedbam

  """
  module load SAMTools/latest
  samtools index $bqsrfile
  """
}

forcaller = recalibratedforcaller.join(indexedbam)

process haplotypecall {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(bamfile), file(baifile) ) from forcaller
  file refs from ref3.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.hapcalled.vcf") ) into calledhaps

  """
  module load GATK/4.1.3.0
  gatk HaplotypeCaller -R ${refs}/hg19_v0_Homo_sapiens_assembly19.fasta -O ${sampleprefix}.hapcalled.vcf -I $bamfile --native-pair-hmm-threads ${params.cpus}
  """
}
