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
log.info "  Functional Metagenomics with Humann2    "
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
        READS FOLDER: ${params.reads}
        OUTPUT FOLDER ${params.output_dir}
        """
        .stripIndent()

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run nibsbioinformatics/core/workflows/metagenomics/humann2_functional.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--reads                         READS FOLDER              Folder where paired end fastq reads are located"
    log.info "--output_dir                    OUTPUT FOLDER             Folder where output sub-folders and results will be copied"
    exit 1
}

// initialisation of parameters before passed by command line or config file

params.reads          = null
params.output_dir     = "."

Channel
    .fromFilePairs("$params.reads/*_{R1,R2}*.fastq.gz")
    .ifEmpty { error "Cannot find any reads matching ${params.reads}"}
    .set { samples_ch }

process characteriseReads {

  tag "humann2 $sampleId"
  cpus 8
  queue 'WORK'
  time '12h'
  memory '32 GB'
  containerOptions = "-B ${params.reads} -B ${params.output_dir}"

  publishDir "${params.output_dir}", mode: 'copy'

  input:
  set sampleId, file(reads) from samples_ch

  output:
  file("${sampleId}/*.tsv")
  file("${sampleId}_concat_humann2_temp/${sampleId}_concat_metaphlan_bowtie2.txt")
  file("${sampleId}_concat_humann2_temp/${sampleId}_concat_metaphlan_bugs_list.tsv")

  script:

  """
  cat $reads >${sampleId}_concat.fastq.gz

  humann2 \
  --input ${sampleId}_concat.fastq.gz \
  --output ${sampleId} \
  --threads ${task.cpus}
  """

}
