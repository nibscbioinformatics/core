#!/opt/software/conda2/envs/NextFlow/bin/nextflow

// Copyright (C) 2019 NIBSC/MHRA
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

// initialisation of parameters before passed by command line or config file
params.picard = null
params.bams = null
params.realign = null
params.fastqs = null
params.output_dir = null
params.reference = null
params.dictionary = null
params.intervals = null
params.cohort_ploidy_model = null
cohort_cnv_caller_model = null
params.help = null


log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  Detect germline CNVs using GATK Toolkit    "
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
        DICTIONARY: ${params.dictionary}
        INTERVALS: ${params.intervals}
        COHORT PLOIDY MODEL: ${params.cohort_ploidy_model}
        COHORT CNV CALLER MODEL: ${params.cohort_cnv_caller_model}
        OUTPUT FOLDER ${params.output_dir}
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
    log.info "--fastqs                       FASTQ FOLDER                 Folder where paired end fastq reads files are located"
    log.info "--bams                         BAMS FOLDER                  Folder where bam files are located"
    log.info "--realign                      BOOLEAN                      Use if you need to realign a bamfile input"
    log.info "--reference                    FASTA reference              Fasta reference file"
    log.info "--dictionary                   DICT file                    Dict file from reference file"
    log.info "--intervals                    INTERVALS list               Filtered intervals used in cohort mode"
	log.info "--contig_ploidy                CONTIG PLOIDY list           TSV file with contig ploidy priors"
	log.info "--cohort_ploidy_model          COHORT PLOIDY model          Cohort ploidy model folder"
	log.info "--cohort_cnv_caller_model      COHORT CNV CALLER model      Cohort CNV caller model folder"
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

if (params.fastqs) {

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
	module load SAMTools
	module load BWA
	
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
    set sampleId, file("${sampleId}_sorted.bam") into aligned_bams_ch

    script:
    """
	module load SAMTools
	module load BWA
	
	java -jar ${params.picard} AddOrReplaceReadGroups I=${umappedBam} O=${sampleId}_group.bam RGID=${sampleId} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${sampleId}
    
	samtools bam2fq -T RX ${sampleId}_group.bam | \
    bwa mem -p -t ${task.cpus} -C -M -R \"@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:Illumina\" \
    ${params.reference} - | \
    samtools sort -@ ${task.cpus} -o ${sampleId}_sorted.bam -
    """
  }
}

process CollectReadCounts {

  tag "Collect Read Counts"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
  set sampleId, file(bamfile) from aligned_bams_ch

  output:
  set sampleId, file("${sampleId}.counts.hdf5") into hdf5_ch
  file("${bamfile}.bai")

  script:
  """
  samtools index $bamfile
  
  gatk CollectReadCounts \
  -I ${bamfile} \
  -L ${params.intervals} \
  --interval-merging-rule OVERLAPPING_ONLY \
  -O "${sampleId}.counts.hdf5" \
  --tmp-dir /home/AD/praposo/Temp 
  """
}

process DetermineGermlineContigPloidy {

  tag "Determine Contig Ploidy"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}", mode: 'copy'

  input:
  set sampleId, file(hdf5) from hdf5_ch

  output:
  set file("ploidy_case_${sampleId}"), sampleId, file("${hdf5}") into germline_contig_ploidy_ch
  
  script:
  """
  gatk DetermineGermlineContigPloidy \
  -I ${hdf5} \
  --model ${params.cohort_ploidy_model} \
  --output-prefix case \
  --tmp-dir /home/AD/praposo/Temp \
  --output "ploidy_case_${sampleId}" \
  --verbosity DEBUG
  """
}

process GermlineCNVCaller {

  tag "CNV Caller"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}", mode: 'copy'

  input:
  set file(ploidy_case), sampleId, file(hdf5) from germline_contig_ploidy_ch

  output:
  set file("${ploidy_case}"), sampleId, file("cnv_caller_case_${sampleId}") into cnv_caller_ch

  script:
  """
  gatk GermlineCNVCaller \
  --run-mode CASE \
  -I ${hdf5} \
  --model ${params.cohort_cnv_caller_model} \
  --contig-ploidy-calls ${ploidy_case}/case-calls \
  --output "cnv_caller_case_${sampleId}" \
  --output-prefix case \
  --verbosity DEBUG \
  --tmp-dir /home/AD/praposo/Temp
  """
}

process PostprocessGermlineCNVCalls {

  tag "CNV Caller"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}", mode: 'copy'

  input:
  set file(ploidy_case), sampleId, file(cnv_caller) from cnv_caller_ch

  output:
  file("genotyped_intervals_case_${sampleId}.vcf.gz")
  file("genotyped_segments_case_${sampleId}.vcf.gz")
  file("genotyped_copy_ratios_case_${sampleId}.vcf.gz")

  script:
  """
  gatk PostprocessGermlineCNVCalls \
  --model-shard-path ${params.cohort_cnv_caller_model} \
  --calls-shard-path ${cnv_caller}/case-calls \
  --contig-ploidy-calls ${ploidy_case}/case-calls \
  --allosomal-contig X --allosomal-contig Y \
  --output-genotyped-intervals "genotyped_intervals_case_${sampleId}.vcf.gz" \
  --output-genotyped-segments "genotyped_segments_case_${sampleId}.vcf.gz" \
  --output-denoised-copy-ratios "genotyped_copy_ratios_case_${sampleId}.vcf.gz" \
  --sequence-dictionary ${params.dictionary} \
  --tmp-dir /home/AD/praposo/Temp \
  --verbosity ERROR
  """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed\n" : "Oops .. something went wrong\n" )
}