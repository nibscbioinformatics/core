#!/opt/software/conda2/envs/NextFlow/bin/nextflow

//This is a nextflow script to perform GATK4 best practice processing and apply haplotypecaller and mutect2

params.referencefolder = "/home/AD/tbleazar/228redo/gatk_references/" //we require there to already be indexed references here
params.outdir = "/home/AD/tbleazar/228redo/output"
//params.filepattern = "/home/AD/tbleazar/228redo/testinput/228_HT1080_test{_R1_001,_R2_001}.fastq.gz"
params.filepattern = "/usr/share/sequencing/projects/228/trimmed/228_HT1080_WT_P15_S1{_R1_001,_R2_001}.fastq.gz"
params.cpus = "32"
params.dbsnp = "b37_dbsnp_138.b37.vcf"
params.goldindels = "b37_Mills_and_1000G_gold_standard.indels.b37.vcf"
params.genomefasta = "hg19_v0_Homo_sapiens_assembly19.fasta"
params.normpanel = "somatic-b37_Mutect2-WGS-panel-b37.vcf"

Channel
  .fromFilePairs(params.filepattern)
  .set { readpairs }

references = Channel
  .fromPath(params.referencefolder)

references.into {
  ref1
  ref2
  ref3
  ref4
  ref5
  ref6
}

process doalignment {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '20 GB'

  input:
  set ( sampleprefix, file(samples) ) from readpairs
  file refs from ref1.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.unsorted.sam") ) into samfile

  """
  module load BWA/latest
  bwa mem -t ${params.cpus} -M -R '@RG\\tID:${sampleprefix}\\tSM:${sampleprefix}\\tPL:Illumina' ${refs}/${params.genomefasta} ${samples[0]} ${samples[1]} > ${sampleprefix}.unsorted.sam
  """
}

process sorttobam {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

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
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

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
  cpus 8
  queue 'WORK'
  time '16h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(markedbamfile) ) from markedbamfortable
  file refs from ref2.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.recal_data.table") ) into recaltable

  """
  module load GATK/4.1.3.0
  gatk BaseRecalibrator -I $markedbamfile --known-sites ${refs}/${params.dbsnp} --known-sites ${refs}/${params.goldindels} -O ${sampleprefix}.recal_data.table -R ${refs}/${params.genomefasta}
  """
}

forrecal = recaltable.join(markedbamforapply)

process applybaserecalibration {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 8
  queue 'WORK'
  time '16h'
  memory '40 GB'

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
  cpus 8
  queue 'WORK'
  time '16h'
  memory '40 GB'

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
forcaller.into {
  forcaller1
  forcaller2
}

process haplotypecall {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(bamfile), file(baifile) ) from forcaller1
  file refs from ref3.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.hapcalled.vcf") ) into calledhaps

  """
  module load GATK/4.1.3.0
  gatk HaplotypeCaller -R ${refs}/${params.genomefasta} -O ${sampleprefix}.hapcalled.vcf -I $bamfile --native-pair-hmm-threads ${params.cpus} --dbsnp ${refs}/${params.dbsnp}
  """
}

process mutectcall {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 32
  queue 'WORK'
  time '48h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(bamfile), file(baifile) ) from forcaller2
  file refs from ref4.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.mutcalled.vcf") ) into calledmuts

  """
  module load GATK/4.1.3.0
  gatk Mutect2 -R ${refs}/${params.genomefasta} -O ${sampleprefix}.mutcalled.vcf -I $bamfile --native-pair-hmm-threads ${params.cpus} --panel-of-normals ${refs}/${params.normpanel}
  """
}

rawvars = calledhaps.join(calledmuts)

process snpindelsplit {
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(hapfile), file(mutfile) ) from rawvars

  output:
  set ( sampleprefix, file("${sampleprefix}.hapcalled.snp.vcf"), file("${sampleprefix}.hapcalled.indel.vcf"), file("${sampleprefix}.mutcalled.snp.vcf"), file("${sampleprefix}.mutcalled.indel.vcf") ) into splitupvars

  """
  module load GATK/4.1.3.0
  gatk SelectVariants -V $hapfile -O ${sampleprefix}.hapcalled.snp.vcf -select-type SNP
  gatk SelectVariants -V $hapfile -O ${sampleprefix}.hapcalled.indel.vcf -select-type INDEL
  gatk SelectVariants -V $mutfile -O ${sampleprefix}.mutcalled.snp.vcf -select-type SNP
  gatk SelectVariants -V $mutfile -O ${sampleprefix}.mutcalled.indel.vcf -select-type INDEL
  """
}

process hardfilter {
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(hapsnp), file(hapindel), file(mutsnp), file(mutindel) ) from splitupvars
  file refs from ref5.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.filtered.snp.vcf"), file("${sampleprefix}.germline.filtered.indel.vcf"), file("${sampleprefix}.somatic.filtered.snp.vcf"), file("${sampleprefix}.somatic.filtered.indel.vcf") ) into filteredvars

  """
  module load GATK/4.1.3.0
  gatk VariantFiltration -O ${sampleprefix}.germline.filtered.snp.vcf -V $hapsnp -R ${refs}/${params.genomefasta} --filter-name snpfilter --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
  gatk VariantFiltration -O ${sampleprefix}.germline.filtered.indel.vcf -V $hapindel -R ${refs}/${params.genomefasta} --filter-name indelfilter --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0"
  gatk VariantFiltration -O ${sampleprefix}.somatic.filtered.snp.vcf -V $mutsnp -R ${refs}/${params.genomefasta} --filter-name snpfilter --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
  gatk VariantFiltration -O ${sampleprefix}.somatic.filtered.indel.vcf -V $mutindel -R ${refs}/${params.genomefasta} --filter-name indelfilter --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0"
  """
}

process remergevars {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(germlinesnp), file(germlineindel), file(somaticsnp), file(somaticindel) ) from filteredvars

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.vcf"), file("${sampleprefix}.germline.vcf.idx"), file("${sampleprefix}.somatic.vcf"), file("${sampleprefix}.somatic.vcf.idx") ) into (germsomvars1, germsomvars2)

  """
  module load GATK/4.1.3.0
  gatk MergeVcfs -I $germlinesnp -I $germlineindel -O ${sampleprefix}.germline.vcf
  gatk MergeVcfs -I $somaticsnp -I $somaticindel -O ${sampleprefix}.somatic.vcf
  """
}

process variantevaluation {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(germline), file(germlineindex), file(somatic), file(somaticindex) ) from germsomvars1
  file refs from ref6.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.eval.grp"), file("${sampleprefix}.somatic.eval.grp") ) into variantevaluations

  """
  module load GATK/4.1.3.0
  gatk VariantEval -eval $germline -O ${sampleprefix}.germline.eval.grp -R ${refs}/${params.genomefasta} -D ${refs}/${params.dbsnp}
  gatk VariantEval -eval $somatic -O ${sampleprefix}.somatic.eval.grp -R ${refs}/${params.genomefasta} -D ${refs}/${params.dbsnp}
  """
}

process effectprediction {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(germline), file(germlineindex), file(somatic), file(somaticindex) ) from germsomvars2

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.annotated.vcf"), file("${sampleprefix}.somatic.annotated.vcf") ) into annotatedvars

  """
  module load snpeff/4.3.1t
  snpEff -Xmx8g hg19 $germline > ${sampleprefix}.germline.annotated.vcf
  snpEff -Xmx8g hg19 $somatic > ${sampleprefix}.somatic.annotated.vcf
  """
}


