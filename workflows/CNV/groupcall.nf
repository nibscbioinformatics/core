#!/opt/software/conda2/envs/NextFlow/bin/nextflow

//This is a nextflow script to perform CNV calling based on a collection of tumour and a collection of normal samples

params.referencefolder = "/usr/share/sequencing/references/homo_sapiens/hg19/gatk_references/" //we require there to already be indexed references here
params.outdir = "/home/AD/tbleazar/test/output"
params.tumourdir = "/usr/share/sequencing/projects/228/trimmed/"
params.normaldir = "/usr/share/sequencing/projects/228/trimmed/"
params.cpus = "32"
params.genomefasta = "hg19_v0_Homo_sapiens_assembly19.fasta"

tumours = Channel
  .fromPath(params.tumourdir)
  
normals = Channel
  .fromPath(params.normaldir)

references = Channel
  .fromPath(params.referencefolder)

process runcnvkit {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 32
  queue 'WORK'
  time '72h'
  memory '40 GB'

  input:
  file refs from references
  file tums from tumours
  file norms from normals
  
  output:
  file(${"results"}) into cnvkitoutput

  """
  module load cnvkit/0.9.6
  cnvkit.py batch ${tums}/*Tumor.bam --normal ${norms}/*Normal.bam \
    --targets my_baits.bed \
    --fasta ${refs}/${params.genomefasta} \
    --output-reference my_reference.cnn --output-dir results/ \
    --diagram --scatter 
  """
}
