#!/opt/software/conda2/envs/NextFlow/bin/nextflow

// Copyright (C) 2019 NIBSC/MHRA
// Author: Francesco Lescai francesco.lescai@nibsc.org

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
log.info "  influenza_classifier: Classify composition of influenza HGR samples    "
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
        DATABASE FASTA FILE: ${params.origin}
        DESTINATION FOLDER: ${params.output_dir}
        """
        .stripIndent()

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run nibsbioinformatics/core/workflows/influenza/influenza_classifier.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--reads                         READS FOLDER              Folder with sample reads"
    log.info "--origin                        PARENTS FASTA FILE        Fasta file with sequences of viral parents"
    log.info "--output_dir                    OUTPUT FOLDER             Output for classification results"
    exit 1
}

// initialisation of parameters before passed by command line or config file

params.reads          = null
params.output_dir     = "."
params.origin         = null

Channel
    .fromFilePairs("$params.output_dir/*_{R1,R2}*.fq.gz")
    .ifEmpty { error "Cannot find any reads matching ${params.reads}"}
    .set { samples_ch }

database_fasta_ch = Channel.fromPath(params.origin)


process createBlastDatabase {

  tag "${dbFasta.baseName}"
  // note:
  // the line below presumes you have cloned our github repository
  // under your home directory, in a folder called CODE
  // this way, the pipeline remains portable on any platform by any
  // user, as long as you have cloned our main repository in this way

  conda "$HOME/CODE/core/workflows/influenza/influenza_conda.yml"
  publishDir "/usr/share/sequencing/references/influenzaDBs", mode: 'copy'

  input:
  file dbFasta from database_fasta_ch

  output:
  file "${dbFasta.baseName}.*" into blast_database_ch

  """
  makeblastdb -in $dbFasta -out ${dbFasta.baseName} -parse_seqids -dbtype nucl
  """

}


process blastSearch {

  tag "$sampleId"
  // note:
  // the line below presumes you have cloned our github repository
  // under your home directory, in a folder called CODE
  // this way, the pipeline remains portable on any platform by any
  // user, as long as you have cloned our main repository in this way

  conda "$HOME/CODE/core/workflows/influenza/influenza_conda.yml"
  publishDir "${params.output_dir}/${sampleId}", mode: 'copy'

  input:
  set sampleId, file(reads) from samples_ch
  file dbBlastPrefix from blast_database_ch

  output:
  file("${sampleId}_blast_results.txt") into blast_results

  """
  zcat $reads | seqkit fq2fa -o ${sampleId}.fa

  blastn \
  -query ${sampleId}.fa \
  -db $dbBlastPrefix \
  -max_target_seqs 1 \
  -num_threads 24 \
  -outfmt '6 qseqid sseqid sgi qstart qend sstart send pident mismatch nident evalue' \
  | sort -k 1,1 -k11,11g > "${sampleId}_blast_results.txt"

  """

}




workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed\n" : "Oops .. something went wrong\n" )
}
