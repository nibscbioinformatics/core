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
params.fastqs = null
params.output_dir     = "."
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
        REFERENCE: ${params.reference}
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
    log.info "--references                   FASTA reference           Fasta reference file"
    log.info "--output_dir                   OUTPUT FOLDER             Folder where output reports and data will be copied"
    exit 1
}



if (params.fastqs) {
  Channel
      .fromFilePairs("$params.fastqs/*_{R1,R2}*.fastq.gz,")
      .ifEmpty { error "Cannot find any reads matching ${params.fastqs}"}
      .set { samples_ch }
  mode = 'fastq'
}
if (params.bams) {
  Channel
      .fromPath("$params.bams/*.bam")
      .ifEmpty { error "Cannot find any file matching ${params.bams}"}
      .set { samples_ch }
  mode = 'bam'
}
else {
  error "you have not specified any input parameters"
}
