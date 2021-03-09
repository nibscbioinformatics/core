#!/opt/software/conda2/envs/NextFlow/bin/nextflow

// Copyright (C) 2019 NIBSC/MHRA
// Author: Tom Bleazard (tom.bleazard@nibsc.org)
// Author: Pedro Raposo (pedro.raposo@nibsc.org)

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Initialisation of parameters before passed by command line or config file

params.bams = null
params.realign = null
params.fastqs = null
params.output_dir = null
params.reference = null
params.intervals = null
params.scale = null
params.cpus = 32
params.dbsnp = null
params.goldindels = null
params.pon = null
params.gnomad = null
params.funcotator = null
params.version = null
params.help = null


log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  Detect SNP and indels using GATK Toolkit  "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) NIBSC/MHRA"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""
log.info """\
        PARAMETERS RECEIVED:
        --------------------------------
        READS FOLDER: ${params.fastqs}
        BAMS FOLDER: ${params.bams}
        REALIGN THE BAMS? ${params.realign}
        REFERENCE: ${params.reference}
        INTERVALS: ${params.intervals}
        SCALE: ${params.scale}
        CPUS: ${params.cpus}
        DBSNP FILE: ${params.dbsnp}
        GOLDINDEL FILE: ${params.goldindels}
        PANEL OF NORMALS FILE: ${params.pon}
        GNOMAD FILE: ${params.gnomad}
        FUNCOTATOR FOLDER: ${params.funcotator}
        GENOME VERSION: ${params.version}
		OUTPUT FOLDER: ${params.output_dir}
        """
        .stripIndent()

		
if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  Detect germline CNVs using GATK Toolkit  "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run nibsbioinformatics/core/workflows/GATK/build_pon.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--fastqs                       FASTQ folder                 Folder where paired end fastq reads files are located"
    log.info "--bams                         BAMS folder                  Folder where bam files are located"
    log.info "--realign                      BOOLEAN                      Use if you need to realign a bamfile input"
    log.info "--reference                    FASTA REFERENCE file         Fasta reference file"
    log.info "--intervals                    INTERVALS list               Intervals used in analysis"
	log.info "--scale                        SCALE string                 Scale of analysis ('wgs' or 'exome')"
	log.info "--cpus                         CPUS integer                 CPUs to be used in the analysis"
	log.info "--dbsnp                        DBSNP file                   dbSNP file"
	log.info "--goldindels                   GOLDINDEL file               GoldIndel file"
	log.info "--pon                          PANEL OF NORMALS file        Panel of normals made from 'built_pon.nf'"
	log.info "--gnomad                       GNOMAD file                  Gnomad file"
	log.info "--funcotator                   FUNCOTATOR folder            Downloaded FUNCotator folder"
	log.info "--version                      Genome VERSION string        Genome version used ('hg19' or 'hg38')"
    log.info "--output_dir                   OUTPUT FOLDER                Folder where output reports and data will be copied"	
    exit 1
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
    module load SAMTools/latest

    bwa mem -t ${task.cpus} -M -R \"@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:Illumina\" \
    ${params.reference} ${reads} | \
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

    output:
    tuple val(sampleId), file("${sampleId}_sorted.bam") into aligned_bams_ch


    script:
    """
    module load BWA/latest
    module load SAMTools/latest

    samtools bam2fq -T RX ${umappedBam} | \
    bwa mem -p -t ${task.cpus} -C -M -R \"@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:Illumina\" \
    ${params.reference} - | \
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
  module load SAMTools/latest
  gatk MarkDuplicates -I $sortedbamfile -M ${sampleprefix}.metrics.txt -O ${sampleprefix}.marked.bam
  samtools index ${sampleprefix}.marked.bam
  """
}

process baserecalibrationtable {
  cpus 8
  queue 'WORK'
  time '16h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(markedbamfile) ) from markedbamfortable

  output:
  set ( sampleprefix, file("${sampleprefix}.recal_data.table") ) into recaltable
  
  script:
  optinalCommand = ""
  if (params.scale == "exome") {
    optionalCommand = "-L ${params.intervals}"
  }
  
  """
  module load GATK/4.1.3.0
  gatk BaseRecalibrator -I $markedbamfile --known-sites ${params.dbsnp} --known-sites ${params.goldindels} -O ${sampleprefix}.recal_data.table -R ${params.reference} ${optinalCommand}
  """
}

forrecal = recaltable.join(markedbamforapply)

process applybaserecalibration {
  publishDir "$params.output_dir/alignments", mode: "copy"
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
  publishDir "$params.output_dir/alignments", mode: "copy"
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

  output:
  set ( sampleprefix, file("${sampleprefix}.hapcalled.vcf") ) into calledhaps

  script:
  optinalCommand = ""
  if (params.scale == "exome") {
    optionalCommand = "-L ${params.intervals}"
  }
  
  """
  module load GATK/4.1.3.0
  gatk HaplotypeCaller -R ${params.reference} -O ${sampleprefix}.hapcalled.vcf -I $bamfile --native-pair-hmm-threads ${params.cpus} --dbsnp ${params.dbsnp} ${optionalCommand}
  """
}

process mutectcall {
  cpus 32
  queue 'WORK'
  time '48h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(bamfile), file(baifile) ) from forcaller2

  output:
  set ( sampleprefix, file("${sampleprefix}.mutcalled.vcf"), file("${sampleprefix}.mutcalled.vcf.stats") ) into calledmuts

  script:
  optinalCommand = ""
  if (params.scale == "exome") {
    optionalCommand = "-L ${params.intervals}"
  }
  
  """
  module load GATK/4.1.3.0
  gatk Mutect2 -R ${params.reference} -O ${sampleprefix}.mutcalled.vcf -I $bamfile --native-pair-hmm-threads ${params.cpus} --panel-of-normals ${params.pon} --germline-resource ${params.gnomad} ${optionalCommand}
  """
}

process mutectfilter {
  cpus 32
  queue 'WORK'
  time '24h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(mutvcf), file(mutstats) ) from calledmuts

  output:
  set ( sampleprefix, file("${sampleprefix}.mutcalled.filtered.vcf") ) into filteredmuts
  
  script:
  optinalCommand = ""
  if (params.scale == "exome") {
    optionalCommand = "-L ${params.intervals}"
  }
  
  """
  module load GATK/4.1.3.0
  gatk IndexFeatureFile -F $mutvcf
  gatk FilterMutectCalls -R ${params.reference} -V $mutvcf -O ${sampleprefix}.mutcalled.filtered.vcf ${optionalCommand}
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

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.filtered.snp.vcf"), file("${sampleprefix}.germline.filtered.indel.vcf"), file("${sampleprefix}.somatic.filtered.snp.vcf"), file("${sampleprefix}.somatic.filtered.indel.vcf") ) into filteredvars

  """
  module load GATK/4.1.3.0
  gatk VariantFiltration -O ${sampleprefix}.germline.filtered.snp.vcf -V $hapsnp -R ${params.reference} --filter-name snpfilter --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
  gatk VariantFiltration -O ${sampleprefix}.germline.filtered.indel.vcf -V $hapindel -R ${params.reference} --filter-name indelfilter --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0"
  gatk VariantFiltration -O ${sampleprefix}.somatic.filtered.snp.vcf -V $mutsnp -R ${params.reference} --filter-name snpfilter --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
  gatk VariantFiltration -O ${sampleprefix}.somatic.filtered.indel.vcf -V $mutindel -R ${params.reference} --filter-name indelfilter --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0"
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
  publishDir "$params.output_dir/analysis", mode: "copy"
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(germline), file(germlineindex), file(somatic), file(somaticindex) ) from germsomvars1

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.eval.grp"), file("${sampleprefix}.somatic.eval.grp") ) into variantevaluations

  """
  module load GATK/4.1.3.0
  gatk VariantEval -eval $germline -O ${sampleprefix}.germline.eval.grp -R ${params.reference} -D ${params.dbsnp}
  gatk VariantEval -eval $somatic -O ${sampleprefix}.somatic.eval.grp -R ${params.reference} -D ${params.dbsnp}
  """
}

process effectprediction {
  publishDir "$params.output_dir/analysis", mode: "copy"
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
  publishDir "$params.output_dir/analysis", mode: "copy"
  cpus 8
  queue 'WORK'
  time '8h'
  memory '40 GB'

  input:
  set ( sampleprefix, file(germline), file(germlineindex), file(somatic), file(somaticindex) ) from germsomvars3

  output:
  set ( sampleprefix, file("${sampleprefix}.germline.funcotated.maf"), file("${sampleprefix}.somatic.funcotated.maf") ) into annotatedvars

  """
  module load GATK/4.1.3.0
  
  gatk Funcotator -R ${params.reference} -V $germline -O ${sampleprefix}.germline.funcotated.maf --output-file-format MAF --data-sources-path ${params.funcotator} --ref-version ${params.version} --annotation-default tumor_barcode:${sampleprefix} --remove-filtered-variants true --transcript-selection-mode BEST_EFFECT
  gatk Funcotator -R ${params.reference} -V $somatic -O ${sampleprefix}.somatic.funcotated.maf --output-file-format MAF --data-sources-path ${params.funcotator} --ref-version ${params.version} --annotation-default tumor_barcode:${sampleprefix} --remove-filtered-variants true --transcript-selection-mode BEST_EFFECT
  """
}
