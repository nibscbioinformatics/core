#!/opt/software/conda2/envs/NextFlow/bin/nextflow

// Copyright (C) 2019 NIBSC/MHRA
// Author: Francesco Lescai francesco.lescai@nibsc.org
// Author: Thomas Bleazard thomas.bleazard@nibsc.org

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
params.bams = null
params.realign = null
params.fastqs = null
params.output_dir     = "."
params.reference = null
params.germline_resource = null
params.intervals = null
params.help = null

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  Build a Panel of Normals for Use with Mutect2    "
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
        GERMLINE RES: ${params.germline_resource}
        INTERVALS: ${params.intervals}
        OUTPUT FOLDER ${params.output_dir}
        """
        .stripIndent()

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run nibsbioinformatics/core/workflows/somatic/build_pon.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--fastqs                       FASTQ FOLDER              Folder where paired end fastq reads files are located"
    log.info "--bams                         BAMS FOLDER               Folder where bam files are located"
    log.info "--realign                      BOOLEAN                   Use if you need to realign a bamfile input"
    log.info "--reference                   FASTA reference           Fasta reference file"
    log.info "--germline_resource            GERMLINE resource         Gnomad file AF-only VCF to be used as germline resource"
    log.info "--intervals                    INTERVALS list            BED file in case of capture (exome or panel)"
    log.info "--output_dir                   OUTPUT FOLDER             Folder where output reports and data will be copied"
    exit 1
}


if (params.fastqs) {
  Channel
      .fromFilePairs("$params.fastqs/*_{R1,R2}*.fastq.gz,")
      .ifEmpty { error "Cannot find any reads matching ${params.fastqs}"}
      .set { samples_ch }
}
if (params.bams) {
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
    module load SAMTools/1.10

    samtools bam2fq -T RX ${umappedBam} | \
    bwa mem -p -t ${task.cpus} -C -M -R \"@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:Illumina\" \
    ${params.reference} - | \
    samtools sort -@ ${task.cpus} -o ${sampleId}_sorted.bam -
    """

  }

}


// in the following process we call the VCF file of each sample
// used to create the panel of normal
// we use the bam files as channel, as it is created either by an alignment process
// or by the path in the beginning
// in both cases is a tuple with sample name so we can use it for the name
// of the VCF file

process CreatTumorOnlyCalls {

  tag "Tumor Calls"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}/${sampleName}", mode: 'copy'

  input:
  set sampleName, file(bamfile) from aligned_bams_ch

  output:
  file("${sampleName}_normal.vcf.gz") into vcf_ch

  script:
  """
  module load GATK/4.1.3.0
  module load SAMTools/1.10

  samtools index $bamfile

  gatk Mutect2 \
  -R $params.reference \
  -I $bamfile \
  --max-mnp-distance 0 \
  -O "${sampleName}_normal.vcf.gz"

  """
}

// Then we need to create a GenomicsDB
// the main issue here is to use a variable number of VCF files
// each introduced with -V into the GATK command line
// so I'm going to use an array and then a join

process GenomicsDB {

  tag "GenomicsDB"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}", mode: 'copy'

  input:
  file(vcfs) from vcf_ch.collect()

  output:
  file("pon_db") into pondb_ch

  script:
  vcfList = []
  vcfs.each() { a -> vcfList.add("-V " + a) }
  inputVcfs = vcfList.join(" ")
  intervalsCommand = ""
  if (params.intervals) {
    intervalsCommand = "-L ${params.intervals} "
  }

  """
  module load GATK/4.1.3.0

  gatk GenomicsDBImport \
  -R reference.fasta ${intervalsCommand}\
  --genomicsdb-workspace-path pon_db \
  ${inputVcfs}
  """

}


// once we have all the files we can now create the panel of normals


process CreatePanelOfNormals {

  tag "createPON"
  cpus 1
  queue 'WORK'
  memory { 12.GB * task.attempt }
  time { 24.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}", mode: 'copy'

  input:
  file("pon_db") from pondb_ch

  output:
  file("pon.vcf.gz")

  script:
  """
  module load GATK/4.1.3.0

  gatk CreateSomaticPanelOfNormals \
  -R $params.reference \
  --germline-resource $params.germline_resource \
  -V gendb://pon_db \
  -O pon.vcf.gz
  """
}



workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed\n" : "Oops .. something went wrong\n" )
}
