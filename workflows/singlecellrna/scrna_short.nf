#!/opt/software/conda2/envs/NextFlow/bin/nextflow

// Copyright (C) 2019 NIBSC/MHRA
// Author: Francesco Lescai francesco.lescai@nibsc.org
// Author: Leo Perfect leo.perfect@nibsc.org
// Author: Pedro Raposo pedro.raposo@nibsc.org

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
log.info "  Single Cell RNA Sequencing Workflow    "
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
        PROJECT FOLDER: ${params.project}
        OUTPUT FOLDER ${params.output_dir}
        REFERENCE ${params.reference}
        """
        .stripIndent()

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run nibsbioinformatics/core/workflows/singlecellrna/scrna_short.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--project                       PROJECT FOLDER            Folder containing one sub folder per sample (each with reads)"
    log.info "--output_dir                    OUTPUT FOLDER             Folder where output sub-folders and results will be copied"
    log.info "--reference                     REFERENCE TRANSCRIPTOME   File referering to the reference transcriptome to be used"
    exit 1
}

// initialisation of parameters before passed by command line or config file

params.project          = null
params.output_dir     = "."

samples_ch = Channel.fromPath( "$params.project/*", type: 'dir' )


// This first process uses the CellRanger suite in order to process the reads per sample
// collapse the UMIs and identify cell barcodes
// creating the genome alingments as well as the expression counts

process CellRangerCount {

  tag "counting"
  cpus 32
  queue 'WORK'
  time '2d'
  memory '250 GB'

  publishDir "$params.output_dir/counting/$sampleId", mode: 'copy'

  input:
  set sampleId, val(sampleFolder) from samples_ch

  output:
  file(XXXX) into counts_ch
  file(YYYY) into alignments_ch
  file(ZZZZ) into reports_ch

  script:
  """
  cellranger count XXXXXXXXXXX
  """

}


// Next we use the Seurat package in order to aggregage the previously generated counts


process Aggregate {

  tag "aggregate"
  cpus 32
  queue 'WORK'
  time '2d'
  memory '250 GB'

  publishDir "$params.output_dir/aggregated", mode: 'copy'

  input:
  file(counts) from counts_ch.collect()

  output:
  file(XXXX) into aggregate_filtered_ch, aggregate_unfiltered_ch

  script:
  """
  Rscript $HOME/CODE/core/workflows/singlecellrna/aggregate.R
  """

}


// A minimum set of exploratory analyses are then run on unfiltered and filtered data

process ExploreUnfiltered {

  tag "exploreUnfiltered"
  cpus 12
  queue 'WORK'
  time '24h'
  memory '250 GB'

  publishDir "$params.output_dir/reports", mode: 'copy'

  input:
  file(robject) from aggregate_unfiltered_ch

  output:
  file('unfiltered_report.html') into unfiltered_report_ch

  script:
  """
  Rscript $HOME/CODE/core/workflows/singlecellrna/unfiltered_report.R
  """
}


process ExploreFiltered {

  tag "exploreFiltered"
  cpus 12
  queue 'WORK'
  time '24h'
  memory '250 GB'

  publishDir "$params.output_dir/reports", mode: 'copy'

  input:
  file(robject) from aggregate_filtered_ch

  output:
  file('filtered_report.html') into filtered_report_ch

  script:
  """
  Rscript $HOME/CODE/core/workflows/singlecellrna/filtered_report.R
  """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed\n" : "Oops .. something went wrong\n" )
}
