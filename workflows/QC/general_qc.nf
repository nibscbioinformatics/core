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

params.help = null

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  Generic FASTQC and MultiQC report    "
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
        OUTPUT FOLDER ${params.output_dir}
        """
        .stripIndent()

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run nibsbioinformatics/core/workflows/QC/general_qc.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--fastqs                       FASTQ FOLDER              Folder where paired end fastq reads files are located"
    log.info "--bams                         BAMS FOLDER               Folder where bam files are located"
    log.info "--output_dir                   OUTPUT FOLDER             Folder where output reports and data will be copied"
    exit 1
}

// initialisation of parameters before passed by command line or config file

params.output_dir     = "."

if (params.fastqs) {
  Channel
      .fromFilePairs("$params.fastqs/*_{R1,R2}*.fastq.gz,")
      .ifEmpty { error "Cannot find any reads matching ${params.fastqs}"}
      .set { samples_ch }
  params.bams = null
}
if (params.bams) {
  Channel
      .fromFilePairs("$params.bams/*.bam")
      .ifEmpty { error "Cannot find any file matching ${params.bams}"}
      .set { samples_ch }
  params.fastqs = null
}
else {
  error "you have not specified any input parameters"
}


process dofastqc {
  publishDir "$params.qcdir", mode: "copy"
  cpus 1
  queue 'WORK'
  time '8h'
  memory '10 GB'

  input:
  set val(sampleId), file(reads) from samples_ch

  output:
  file "*_fastqc.{zip,html}" into fastqcout_ch

  """
  module load FastQC/latest
  fastqc $reads
  """
}


process domultiqc {
  publishDir "$params.qcdir", mode: "copy"
  cpus 8
  queue 'WORK'
  time '8h'
  memory '50 GB'

  input:
  file "qccollected/*" from fastqcout_ch.collect()

  output:
  file "multiqc_report.html" into multiqc_report
  file "multiqc_data"

  """
  module load multiqc/1.7
  multiqc --interactive qccollected
  """
}
