#!/opt/software/conda2/envs/NextFlow/bin/nextflow

// Copyright (C) 2020 NIBSC/MHRA
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
params.dictionary = null
params.intervals = null
params.mappability_track = null
params.segmental_dup_track = null
params.lcr_intervals = null
params.par_intervals = null
params.contig_ploidy = null
params.scale = null
params.tmp_dir = null
params.help = null


log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  Build a Cohort Model for germline CNV detection using GATK Toolkit  "
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
        MAPPABILITY TRACK: ${params.mappability_track}
        SEGMENTAL DUP TRACK: ${params.segmental_dup_track}
        LCR INTERVALS: ${params.lcr_intervals}
        PAR INTERVALS: ${params.par_intervals}
        CONTIG PLOIDY: ${params.contig_ploidy}
		SCALE: ${params.scale}
		TEMPORARY DIRECTORY: ${params.tmp_dir}
        OUTPUT FOLDER ${params.output_dir}
        """
        .stripIndent()


if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  Build model for the detection of germline CNVs using GATK Toolkit  "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run nibsbioinformatics/core/workflows/GATK/GATK_CNV_COHORT.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--fastqs                       FASTQ FOLDER                 Folder where paired end fastq reads files are located"
    log.info "--bams                         BAMS FOLDER                  Folder where bam files are located"
    log.info "--realign                      BOOLEAN                      Use if you need to realign a bamfile input"
    log.info "--reference                    FASTA reference              Fasta reference file"
    log.info "--dictionary                   DICT file                    Dict file from reference file"	
	log.info "--intervals                    INTERVALS list               BED file in case of capture (exome or panel)"
    log.info "--mappability_track            MAPPABILITY track            Umap single-read mappability track BED file (exome or panel)"
    log.info "--segmental_dup_track          SEGMENTAL DUP track          Segmental-duplication track BED file (exome or panel)"	
	log.info "--lcr_intervals                LCR INTERVAL list            BED file with low-complexity regions to exclude (exome or panel)"
	log.info "--par_intervals                PAR INTERVAL list            BED file with pseudoautossomal regions to exclude (exome or panel)"
	log.info "--scale                        STRING                       'wgs' for whole-genome sequencing or 'exome' for exome or panel"
	log.info "--contig_ploidy                CONTIG PLOIDY list           TSV file with contig ploidy priors"
	log.info "--tmp_dir                      TEMPORARY folder             Folder where the temporary data is stored"	
    log.info "--output_dir                   OUTPUT FOLDER                Folder where output reports and data will be copied"
    exit 1
}

if (params.dictionary) {
  Channel
      .fromPath("$params.dictionary")
      .ifEmpty { error "Cannot find any dictonary file matching ${params.dictionary}"}
      .set { dictionary_ch }
}
if (params.mappability_track) {
  Channel
      .fromPath("$params.mappability_track")
      .ifEmpty { error "Cannot find any mappability track matching ${params.mappability_track}"}
      .set { map_track_ch }
}
if (params.segmental_dup_track) {
  Channel
      .fromPath("$params.segmental_dup_track")
      .ifEmpty { error "Cannot find any segmental-duplication track matching ${params.segmental_dup_track}"}
      .set { seg_dup_track_ch }
}
if (params.contig_ploidy) {
  Channel
      .fromPath("$params.contig_ploidy")
      .ifEmpty { error "Cannot find any contig ploidy prior matching ${params.contig_ploidy}"}
      .set { contig_ploidy_ch }
}
if (params.lcr_intervals) {
  Channel
      .fromPath("$params.lcr_intervals")
      .ifEmpty { error "Cannot find any LCR BED file matching ${params.lcr_intervals}"}
      .set { lcr_ch }
}
if (params.par_intervals) {
  Channel
      .fromPath("$params.par_intervals")
      .ifEmpty { error "Cannot find any PAR BED file matching ${params.par_intervals}"}
      .set { par_ch }
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
// align the fastq to be used for the creation of the model

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


// If we have alignments but we want to realign with BWA MEM

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
    set val(sampleId), file("${sampleId}.chr1toY.bam") into aligned_bams_ch
	file("${sampleId}.chr1toY.bam.bai")

    script:
    """
	module load SAMTools/latest
	module load BWA/latest
	
    samtools bam2fq -T RX ${umappedBam} | \
    bwa mem -p -t ${task.cpus} -C -M -R \"@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:Illumina\" \
    ${params.reference} - | \
    samtools sort -@ ${task.cpus} -o ${sampleId}_sorted.bam -
	
	samtools view -L ${params.chromosomes} -o ${sampleId}.chr1toY.bam ${sampleId}_sorted.bam
	samtools index ${sampleId}.chr1toY.bam
    """
  }
}


// For WES, this steps padds the intervals and for the WGS it divides
// the genome into bins, to be used as segments in the CNV detection

process PreprocessIntervals {

  tag "Preprocess Intervals"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}", mode: 'copy'
  
  output:
  file("preprocessed.interval_list") into intervals_ch

  script:
  scaleCommand = ""
  if ( params.scale == "exome" ) {scaleCommand = "-L ${params.intervals} --bin-length 0 --padding 240"} 
  else if ( params.scale == "wgs" ) {scaleCommand = "--padding 0 --bin-length 1000"}
  
  """
  gatk PreprocessIntervals \
  ${scaleCommand} \
  -R ${params.reference} \
  --interval-merging-rule OVERLAPPING_ONLY \
  -O preprocessed.interval_list \
  --tmp-dir ${params.tmp_dir}
  """
}

intervals_ch.into {
  intervals_branch_1
  intervals_branch_2
  intervals_branch_3
  intervals_branch_4
}


// It is counted the number of reads in the sample for each interval

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
  file(preprocessed_intervals) from intervals_branch_1

  output:
  file("${sampleId}.counts.hdf5") into hdf5_ch

  script:
  """
  gatk CollectReadCounts \
  -I ${bamfile} \
  -L ${preprocessed_intervals} \
  --interval-merging-rule OVERLAPPING_ONLY \
  -O "${sampleId}.counts.hdf5" \
  --tmp-dir ${params.tmp_dir} 
  """
}

hdf5_ch.into {
  hdf5_branch_1
  hdf5_branch_2
  hdf5_branch_3
  hdf5_branch_4
  }

  
// Annotate intervals to remove problematic genomic regions defined by CG content,
// or explicitily mappability and segmentation duplication
  
process AnnotateIntervals {

  tag "Annotate Intervals"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}", mode: 'copy'

  input:
  file(mappability_track) from map_track_ch
  file(segmental_dup_track) from seg_dup_track_ch
  file(preprocessed_intervals) from intervals_branch_2

  output:
  file("annotated_intervals.tsv") into annotated_intervals_ch

  script:
  optinalCommand = ""
  if (params.segmental_dup_track || params.mappability_track) {
    if (!params.segmental_dup_track) {
      optionalCommand = "--mappability-track ${params.mappability_track} "
	} else if (!params.mappability_track) {
      optionalCommand = "--segmental-duplication-track ${params.segmental_dup_track} "
	} else {
	  optionalCommand = "--mappability-track ${params.mappability_track} --segmental-duplication-track ${params.segmental_dup_track} "
	}
  }

  """
  gatk AnnotateIntervals \
  -R ${params.reference} \
  -L ${preprocessed_intervals} ${optionalCommand} \
  --interval-merging-rule OVERLAPPING_ONLY \
  -O "annotated_intervals.tsv" \
  --tmp-dir ${params.tmp_dir} 
  """
}

annotated_intervals_ch.into {
  annotated_intervals_branch_1
  annotated_intervals_branch_2
}


// Exclude counted reads based on the previously annotated intervals,
// genomic low complexity regions, or pseudoautossomal regions 

process FilterIntervals {

  tag "Filter Intervals"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}", mode: 'copy'

  input:
  file(annotated_intervals) from annotated_intervals_branch_1
  file(preprocessed_intervals) from intervals_branch_3
  file(hdf5) from hdf5_branch_1.collect()
  file(lcr_intervals) from lcr_ch
  file(par_intervals) from par_ch

  output:
  file("filtered_intervals.intervals") into filtered_intervals_ch

  script:
  hdf5List = []
  hdf5.each() { a -> hdf5List.add("-I " + a) }
  inputHdf5 = hdf5List.join(" ")

  """
  gatk FilterIntervals \
  -L ${preprocessed_intervals} \
  -XL ${lcr_intervals} \
  -XL ${par_intervals} \
  ${inputHdf5} \
  -imr OVERLAPPING_ONLY \
  -O "filtered_intervals.intervals" \
  --tmp-dir ${params.tmp_dir} 
  """
}

filtered_intervals_ch.into {
  filtered_intervals_branch_1
  filtered_intervals_branch_2
  filtered_intervals_branch_3
}


// It determines baseline contig ploidies using the total read count per contig
// for autosomal and allosomal chromosomes

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
  file(hdf5) from hdf5_branch_2.collect()
  file(contig_ploidy) from contig_ploidy_ch
  file(intervals) from filtered_intervals_branch_2

  output:
  file("ploidy_cohort") into germline_contig_ploidy_ch

  script:
  hdf5List = []
  hdf5.each() { a -> hdf5List.add("-I " + a) }
  inputHdf5 = hdf5List.join(" ")

  """
  gatk DetermineGermlineContigPloidy \
  -imr OVERLAPPING_ONLY \
  --contig-ploidy-priors ${contig_ploidy} \
  --output-prefix cohort \
  -L ${intervals} \
  ${inputHdf5} \
  --tmp-dir ${params.tmp_dir} \
  --output ploidy_cohort \
  --verbosity DEBUG
  """
}

germline_contig_ploidy_ch.into {
  germline_contig_ploidy_branch_1
  germline_contig_ploidy_branch_2
}


// Detects the CNV variantions within the sample, based on its contig ploidy.
// Extra parameters are included for WGS to increase sensitivity

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
  file(hdf5) from hdf5_branch_4.collect()
  file(contig_ploidy) from germline_contig_ploidy_branch_1
  file(annotated_intervals) from annotated_intervals_branch_2
  file(intervals) from filtered_intervals_branch_1

  output:
  file("cnv_caller_cohort") into cnv_caller_ch

  script:
  hdf5List = []
  hdf5.each() { a -> hdf5List.add("-I " + a) }
  inputHdf5 = hdf5List.join(" ")
  
  optionalCommand = ""
  if ( params.scale == "wgs") {optionalCommand = "--class-coherence-length 1000.0 --cnv-coherence-length 1000.0 --enable-bias-factors false --interval-psi-scale 1.0E-6 --log-mean-bias-standard-deviation 0.01 --sample-psi-scale 1.0E-6"}

  """
  gatk GermlineCNVCaller \
  --run-mode COHORT \
  -L ${intervals} \
  --annotated-intervals ${annotated_intervals} \
  ${inputHdf5} \
  -imr OVERLAPPING_ONLY \
  --contig-ploidy-calls ${contig_ploidy}/cohort-calls \
  --output cnv_caller_cohort \
  --output-prefix cohort \
  --verbosity DEBUG \
  ${optionalCommand} \
  --tmp-dir ${params.tmp_dir}
  """
}


// Consolidates the CNV variant results into VCF files 

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
  file(cnv_caller) from cnv_caller_ch
  file(germline_contig_ploidy) from germline_contig_ploidy_branch_2
  file(dictionary) from dictionary_ch

  output:
  file("genotyped_intervals_cohort.vcf.gz")
  file("genotyped_segments_cohort.vcf.gz")
  file("genotyped_copy_ratios_cohort.vcf.gz")

  script:
  """
  gatk PostprocessGermlineCNVCalls \
  --model-shard-path ${cnv_caller}/cohort-model \
  --calls-shard-path ${cnv_caller}/cohort-calls \
  --contig-ploidy-calls ${germline_contig_ploidy}/cohort-calls \
  --allosomal-contig X --allosomal-contig Y \
  --output-genotyped-intervals genotyped_intervals_cohort.vcf.gz \
  --output-genotyped-segments genotyped_segments_cohort.vcf.gz \
  --output-denoised-copy-ratios genotyped_copy_ratios_cohort.vcf.gz \
  --sequence-dictionary ${dictionary} \
  --tmp-dir ${params.tmp_dir} \
  --verbosity ERROR
  """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed\n" : "Oops .. something went wrong\n" )
}