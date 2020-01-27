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
        METADATA: ${params.metadata}
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
    log.info "--metadata                      PROJECT METADATA          Metadata file with 3 columns linking sample name and fastqs"
    log.info "--output_dir                    OUTPUT FOLDER             Folder where output sub-folders and results will be copied"
    log.info "--reference                     TRANSCRIPTOME FOLDER      Folder where the reference transcriptome for cellranger is located"
    exit 1
}

// initialisation of parameters before passed by command line or config file

params.metadata          = null
params.output_dir     = "."


// samples metadata need to be specified in tabular format with the following information (per column)
// SAMPLE NAME (which goes into cell ranger)
// SAMPLE ID(s) again in cell ranger specification from fastq files for merging
// FASTQ FOLDER(s) where files with a name formatted to begin with the provided sample ID are present
// the file has to be tab separated because a coma is used to separate ids and folders

Channel
      .fromPath("${params.metadata}")
      .splitCsv(header: ['sampleID', 'fastqIDs', 'fastqLocs'], sep: '\t')
      .map{ row-> tuple(row.sampleID, row.fastqIDs, row.fastqLocs) }
      .set { metadata_ch }


// This first process uses the CellRanger suite in order to process the reads per sample
// collapse the UMIs and identify cell barcodes
// creating the genome alingments as well as the expression counts



process CellRangerCount {

  tag "counting"
  cpus 32
  queue 'WORK'
  time '2d'
  memory '250 GB'

  publishDir "$params.output_dir/counting/$sampleName", mode: 'copy'

  input:
  set sampleName, fastqIDs, fastqLocs from metadata_ch


  output:
  file("$sampleName/outs/metrics_summary.csv") into cellranger_summary_ch
  file("$sampleName/outs/filtered_feature_bc_matrix/*.gz") into count_files_ch
  tuple val("$sampleName"), file("$sampleName/outs/filtered_feature_bc_matrix/") into count_folders_ch
  tuple file("$sampleName/outs/possorted_genome_bam.bam"), file("$sampleName/outs/possorted_genome_bam.bam.bai") into alignments_ch

  script:

  """
  cellranger count \
  --id=${sampleName} \
  --sample=${fastqIDs} \
  --fastqs=${fastqLocs} \
  --transcriptome=$params.reference
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
  set sampleData from count_folders_ch.collect()

  output:
  file('aggregated_object.RData') into (aggregate_filtered_ch, aggregate_unfiltered_ch)

  script:
  sampleData.each() { k,v -> sampleNamesList << k; countFoldersList << v}
  sampleNames = sampleNamesList.join(",")
  countFolders = countFoldersList.join(",")

  """
  Rscript -e "workdir<-getwd()
  rmarkdown::render('$HOME/CODE/core/workflows/singlecellrna/seurat_scripts/aggregate.Rmd',
    params = list(
      sample_paths = \"$sampleNames\",
      sample_names = \"$countFolders\",
      output_path = workdir),
      knit_root_dir=workdir,
      output_dir=workdir)"
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
  file(aggregatedObj) from aggregate_unfiltered_ch

  output:
  file('analyse_unfiltered.html') into unfiltered_report_ch
  file('aggregated_object_analyzed_unfiltered.RData') into unfiltered_object_ch

  script:
  """
  Rscript -e "workdir<-getwd()
    rmarkdown::render('$HOME/CODE/core/workflows/singlecellrna/seurat_scripts/analyse_unfiltered.Rmd',
    params = list(input_path = $aggregatedObj),
    knit_root_dir=workdir,
    output_dir=workdir)"
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
  file(aggregatedObj) from aggregate_filtered_ch

  output:
  file('analyse_filtered.html') into filtered_report_ch
  file('aggregated_object_analyzed_filtered.RData') into filtered_object_ch

  script:
  """
  Rscript -e "workdir<-getwd()
    rmarkdown::render('$HOME/CODE/core/workflows/singlecellrna/seurat_scripts/analyse_filtered.Rmd',
    params = list(input_path = $aggregatedObj),
    knit_root_dir=workdir,
    output_dir=workdir))"
  """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed\n" : "Oops .. something went wrong\n" )
}
