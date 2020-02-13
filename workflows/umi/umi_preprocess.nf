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
log.info "  Pre-process UMIs for downstream workflows    "
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
        REFERENCE ${params.reference}
        BARCODE1 STRUCTURE ${params.read_structure1}
        BARCODE2 STRUCTURE ${params.read_structure2}
        """
        .stripIndent()

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run nibsbioinformatics/core/workflows/umi/umi_preprocess.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--reads                         READS FOLDER              Folder containing all reads of the project"
    log.info "--output_dir                    OUTPUT FOLDER             Folder where output sub-folders and results will be copied"
    log.info "--reference                     REFERENCE                 File referering to the reference to be used"
    log.info "--read_structure1               BARCODE1 STRUCTURE         Structure of Reads R1 Barcode to be parsed (follow fgbio/Picard format)"
    log.info "--read_structure2               BARCODE2 STRUCTURE         Structure of Reads R2 Barcode to be parsed (follow fgbio/Picard format)"
exit 1
}

// initialisation of parameters before passed by command line or config file

params.reads          = null
params.output_dir     = "."

Channel
    .fromFilePairs("$params.reads/*_{R1,R2}*.fastq.gz")
    .ifEmpty { error "Cannot find any reads matching ${params.reads}"}
    .set { samples_ch }


process FastqToBAM {

  tag "ConvertToBAM"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '6 GB'

  // no publish dir, meaning all files are temporary
  // and will be stored in work directory until deletion

  input:
  set sampleId, file(reads) from samples_ch

  output:
  tuple val(sampleId), file("${sampleId}_converted.bam") into converted_bams_ch


  script:
  """
  module load fgbio/1.1.0
  mkdir tmpFolder

  fgbio --tmp-dir=${PWD}/tmpFolder \
  FastqToBam \
  -i $reads \
  -o "${sampleId}_converted.bam" \
  --read-structures $params.read_structure1 $params.read_structure2 \
  --sample $sampleId \
  --library $sampleId
  """

}

// now we need to see how we can do that because the reads are now
// in a BAM format, unmapped
// however, if the read structure is specified correctly,
// ideally the actual sequence should only contain the template portion
// making the trimming step unnecessary
// see blog http://nilshomer.com/2017/07/05/single-strand-umi-somatic-variant-calling/

// process TrimReads {
//
//   tag "trim reads"
//   cpus 1
//   queue 'WORK'
//   time '2h'
//   memory '6 GB'
//
//
// }


// this step is necessary to allow the following read grouping
// because mapping information and mapping quality are essential
// inputs to next step
// bwa mem doesn't use BAM as input, while it's whole-package bwa kit does
// consider using bwa-kit instead

process AlignBamFile {

  tag "BWA alignment"
  cpus 8
  queue 'WORK'
  time '24h'
  memory '24 GB'

  // no publish dir, meaning all files are temporary
  // and will be stored in work directory until deletion
  // because in this case the aligned bam file is still UMI tagged

  input:
  set sampleId, file(convertedBam) from converted_bams_ch

  output:
  tuple val(sampleId), file("${sampleId}_unsorted.bam") into aligned_bams_ch


  script:
  """
  module load BWA/latest
  module load SAMTools/1.10

  samtools bam2fq -T RX ${convertedBam} | \
  bwa mem -t ${task.cpus} -C -M -R \"@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:Illumina\" \
  ${params.reference} - | \
  samtools view -bS - > ${sampleId}_unsorted.bam
  """

}

// now that the reads are aligned, we can group them by UMI as important preliminary step
// to create the consensus
// one thing is still missing though, which is required by GroupReadsByUmi
// and is the mandatory presence of the MQ tag, informing of the mapping quality of the read mate pair
// this information is essential for the aggressive filtering the tool does, i.e. making sure that
// the UMI reads actually come from the same molecule of origin

// We have chose the Adjacency method, following the nice paper and blog explanation integrated in both
// UMItools and FGBIO
// https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/

process GroupReadsByUmi {

  tag "group UMIs"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '8 GB'

  // no publish dir, meaning all files are temporary
  // and will be stored in work directory until deletion

  input:
  set sampleId, file(alignedBam) from aligned_bams_ch

  output:
  file("${sampleId}_umi_histogram.txt") into umi_histogram_ch
  tuple val(sampleId), file("${sampleId}_umi-grouped.bam") into umi_grouped_bams_ch

  script:
  """
  module load fgbio/1.1.0
  module load samblaster/0.1.24
  module load SAMTools/1.10

  mkdir tmpFolder

  samtools view -h $alignedBam | \
  samblaster -M --addMateTags | \
  samtools view -Sb - >${sampleId}_unsorted_tagged.bam

  fgbio --tmp-dir=${PWD}/tmpFolder \
  GroupReadsByUmi \
  -s Adjacency \
  -i ${sampleId}_unsorted_tagged.bam \
  -o ${sampleId}_umi-grouped.bam \
  -f ${sampleId}_umi_histogram.txt
  """

}

// Now that the reads are organised by UMI groups

process CallMolecularConsensusReads {

  tag "XXXXXXXX"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '8 GB'

  publishDir "$params.output_dir", mode: 'copy'

  input:
  set sampleId, file(groupedBamFile) from umi_grouped_bams_ch

  output:
  file("${sampleId}_umi-consensus.bam") into consensus_bam_ch

  script:
  """
  module load fgbio/1.1.0
  mkdir tmpFolder

  fgbio --tmp-dir=${PWD}/tmpFolder \
  CallMolecularConsensusReads \
  -i $groupedBamFile \
  -o ${sampleId}_umi-consensus.bam \
  -M 1 -S Coordinate
  """

}


workflow.onComplete {

  if( workflow.success ) {
    log.info("\nDone! Workflow completed\n")
    log.info("Removing all intermediate files now\n")
    log.info("Removing ${workflow.workDir}\n")
    deleteWork = workflow.workDir.deleteDir()
    log.info("Removing ${workflow.launchDir}/.nextflow/\n")
    mycache = file("${workflow.launchDir}/.nextflow")
    deleteCache = mycache.deleteDir()
  }
  else {
    log.info("Oops .. something went wrong\n")
    log.info("Pipeline execution stopped with the following message: ${workflow.errorMessage}")
  }
}
