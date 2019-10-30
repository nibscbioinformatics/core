#!/opt/software/conda2/envs/NextFlow/bin/nextflow

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.index = "/usr/share/sequencing/references/rnaseq/Homo_sapiens/ENSEMBL/gencode.v32.transcripts.index"
params.reads = "$baseDir/reads/*_{1,2}.fastq.gz"
params.transcriptome = "/usr/share/sequencing/references/rnaseq/Homo_sapiens/ENSEMBL/gencode.v32.transcripts.fa.gz"
params.outdir = "results"


/*
* Use the following to print an initial message
*/

log.info """\
         SALMON ANALYSIS OF RNASEQ DATA
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

transcriptome_file = file(params.transcriptome)
transcriptome_index = file(params.index)

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch }
