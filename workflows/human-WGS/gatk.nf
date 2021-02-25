#!/opt/software/conda2/envs/NextFlow/bin/nextflow

//This is a nextflow script to perform GATK4 best practice processing and apply haplotypecaller and mutect2

references = Channel
  .fromPath(params.referencefolder)

references.into {
  ref1
  ref2
  ref3
  ref4
  ref5
  ref6
  ref7
  ref8
}

if (params.fastqs) {
  Channel
      .fromFilePairs("$params.fastqs/*_{R1,R2}*.fastq.gz")
      .ifEmpty { error "Cannot find any reads matching ${params.fastqs}"}
      .set { samples_ch }
}
else if (params.bams) {
  if (params.realign){
    Channel
        .fromPath("$params.bams/*.bam")
        .ifEmpty { error "Cannot find any file matching ${params.bams}"}
        .map { file -> tuple(file.baseName, file) }
        .set { bams_ch }
    println("## NB: you have chosen to realign the bam files \n")
  }
  else {
    Channel
        .fromPath("$params.bams/*.bam")
        .ifEmpty { error "Cannot find any file matching ${params.bams}"}
        .map { file -> tuple(file.baseName, file) }
        .set { aligned_bams_ch }
    println("## NB: you have chosen NOT to realign the bam files \n")
  }
}
else {
  error "you have not specified any input parameters"
}


// In case we don't already have alignments, we can just select the samples and
// align the fastq to be used for the creation of the PON

if (params.fastq) {

  process AlignSamples {

    tag "BWA alignment"
    cpus 8
    queue 'WORK'
    time '24h'
    memory '24 GB'

    publishDir "${params.output_dir}/${sampleId}", mode: 'copy'

    input:
    set sampleId, file(reads) from samples_ch

    output:
    tuple val(sampleId), file("${sampleId}_sorted.bam") into aligned_bams_ch


    script:
    """
    module load BWA/latest
    module load SAMTools/1.10

    bwa mem -t ${task.cpus} -M -R \"@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:Illumina\" \
    ${refs}/${params.genomefasta} ${reads} | \
    samtools sort -@ ${task.cpus} -o ${sampleId}_sorted.bam -
    """
  }

}

if (params.bams && params.realign) {

  process realignBams {

    tag "BWA alignment"
    cpus 8
    queue 'WORK'
    time '24h'
    memory '24 GB'

    // no publish dir, meaning all files are temporary
    // and will be stored in work directory until deletion
    // because in this case the aligned bam file is still UMI tagged

    input:
    set sampleId, file(umappedBam) from bams_ch
	file refs from ref1.first()

    output:
    tuple val(sampleId), file("${sampleId}_sorted.bam") into aligned_bams_ch


    script:
    """
    module load BWA/latest
    module load SAMTools/latest

    samtools bam2fq -T RX ${umappedBam} | \
    bwa mem -p -t ${task.cpus} -C -M -R \"@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:Illumina\" \
    ${refs}/${params.genomefasta} - | \
    samtools sort -@ ${task.cpus} -o ${sampleId}_sorted.bam -
    """

  }

}

process markduplicates {
  cpus 8
  queue 'WORK'
  time '14h'
  memory '180 GB'

  input:
  set ( sampleprefix, file(sortedbamfile) ) from aligned_bams_ch

  output:
  set ( sampleprefix, file("${sampleprefix}.marked.bam") ) into (markedbamfortable, markedbamforapply)

  """
  module load GATK/4.1.3.0
  gatk MarkDuplicates -I $sortedbamfile -M ${sampleprefix}.metrics.txt -O ${sampleprefix}.marked.bam
  """
}

process baserecalibrationtable {
  cpus 8
  queue 'WORK'
  time '16h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(markedbamfile) ) from markedbamfortable
  file refs from ref2.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.recal_data.table") ) into recaltable
  
  script:
  
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
  cpus 32
  queue 'WORK'
  time '48h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(bamfile), file(baifile) ) from forcaller1
  file refs from ref3.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.hapcalled.vcf") ) into calledhaps

  script:
 
  """
  module load GATK/4.1.3.0
  gatk HaplotypeCaller -R ${refs}/${params.genomefasta} -O ${sampleprefix}.hapcalled.vcf -I $bamfile --native-pair-hmm-threads ${params.cpus} --dbsnp ${refs}/${params.dbsnp}
  """
}

process mutectcall {
  cpus 32
  queue 'WORK'
  time '48h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(bamfile), file(baifile) ) from forcaller2
  file refs from ref4.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.mutcalled.vcf"), file("${sampleprefix}.mutcalled.vcf.stats") ) into calledmuts

  script:
  
  """
  module load GATK/4.1.3.0
  gatk Mutect2 -R ${refs}/${params.genomefasta} -O ${sampleprefix}.mutcalled.vcf -I $bamfile --native-pair-hmm-threads ${params.cpus} --panel-of-normals /home/AD/praposo/WGS/GIAB/somatic-b37_Mutect2-WGS-panel-b37.vcf --germline-resource ${refs}/${params.gnomad}
  """
}

process mutectfilter {
  cpus 32
  queue 'WORK'
  time '24h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(mutvcf), file(mutstats) ) from calledmuts
  file refs from ref7.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.mutcalled.filtered.vcf") ) into filteredmuts
  
  script:
  
  """
  module load GATK/4.1.3.0
  gatk FilterMutectCalls -R ${refs}/${params.genomefasta} -V $mutvcf -O ${sampleprefix}.mutcalled.filtered.vcf
  """
}

rawvars = calledhaps.join(filteredmuts)

process snpindelsplit {
  cpus 8
  queue 'WORK'
  time '8h'
  memory '120 GB'

  input:
  set ( sampleprefix, file(hapfile), file(mutfile) ) from rawvars

  output:
  set ( sampleprefix, file("${sampleprefix}.hapcalled.snp.vcf"), file("${sampleprefix}.hapcalled.indel.vcf"), file("${sampleprefix}.mutcalled.snp.vcf"), file("${sampleprefix}.mutcalled.indel.vcf") ) into splitupvars

  """
  module load GATK/4.1.3.0
  gatk --java-options "-Xmx120G" SelectVariants -V $hapfile -O ${sampleprefix}.hapcalled.snp.vcf -select-type SNP
  gatk --java-options "-Xmx120G" SelectVariants -V $hapfile -O ${sampleprefix}.hapcalled.indel.vcf -select-type INDEL
  gatk --java-options "-Xmx120G" SelectVariants -V $mutfile -O ${sampleprefix}.mutcalled.snp.vcf -select-type SNP
  gatk --java-options "-Xmx120G" SelectVariants -V $mutfile -O ${sampleprefix}.mutcalled.indel.vcf -select-type INDEL
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
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(germlinesnp), file(germlineindel), file(somaticsnp), file(somaticindel) ) from filteredvars

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.vcf"), file("${sampleprefix}.germline.vcf.idx"), file("${sampleprefix}.somatic.vcf"), file("${sampleprefix}.somatic.vcf.idx") ) into (germsomvars1, germsomvars2, germsomvars3)

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
  set ( sampleprefix, file("${sampleprefix}.germline.annotated.vcf"), file("${sampleprefix}.somatic.annotated.vcf") ) into effectpredicted

  """
  module load snpeff/4.3.1t
  snpEff -Xmx8g hg19 $germline > ${sampleprefix}.germline.annotated.vcf
  snpEff -Xmx8g hg19 $somatic > ${sampleprefix}.somatic.annotated.vcf
  """
}

process annotation {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(germline), file(germlineindex), file(somatic), file(somaticindex) ) from germsomvars3
  file refs from ref8.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.funcotated.maf"), file("${sampleprefix}.somatic.funcotated.maf") ) into annotatedvars

  """
  module load GATK/4.1.3.0
  module load snpsift/4.3.1t
  
  gatk Funcotator -R ${refs}/${params.genomefasta} -V $germline -O ${sampleprefix}.germline.funcotated.maf --output-file-format MAF --data-sources-path ${params.funcotator} --ref-version hg19 --annotation-default tumor_barcode:${sampleprefix} --remove-filtered-variants true --transcript-selection-mode BEST_EFFECT
  gatk Funcotator -R ${refs}/${params.genomefasta} -V $somatic -O ${sampleprefix}.somatic.funcotated.maf --output-file-format MAF --data-sources-path ${params.funcotator} --ref-version hg19 --annotation-default tumor_barcode:${sampleprefix} --remove-filtered-variants true --transcript-selection-mode BEST_EFFECT
  """
}
