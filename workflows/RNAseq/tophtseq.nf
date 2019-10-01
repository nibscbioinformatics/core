//This is a nextflow script to process RNA-seq fastq data through trimming to tophat to htseq-count

params.indir = "/usr/share/sequencing/projects/272/raw_data/"
params.adapterfile = "/home/AD/tbleazar/usrlocalmark/TruSeq-adapters-recommended.fa"
params.outdir = "/usr/share/sequencing/projects/272/"
params.referencefolder = "/usr/share/sequencing/projects/132/reference/" //and it is required that this have an indexed genome already
params.referencegenome = "GRCh38.fa"
params.genegtf = "hsa.ensembl.92.chr.gtf"
params.cpus = "32"

inputdirectory = file(params.indir)

inputdirectory.into {
  inputdirectory1
  inputdirectory2
}

references = Channel
  .fromPath(params.referencefolder)

references.into {
  ref1
  ref2
  ref3
  ref4
}

process mergeforwardlanes {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  file(inputdir) from inputdirectory1

  output:
  file('*.forward.fastq.gz') into forwardfiles

  """
  for samplename in `ls $inputdir | grep _L001_R1_001.fastq.gz | sed 's/_L001_R1_001.fastq.gz//g'`
  do
  cat $inputdir/\${samplename}_L001_R1_001.fastq.gz $inputdir/\${samplename}_L002_R1_001.fastq.gz $inputdir/\${samplename}_L003_R1_001.fastq.gz $inputdir/\${samplename}_L004_R1_001.fastq.gz > \${samplename}.forward.fastq.gz
  done
  """
}

process forwardsets {
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB

  input:
  file(forwardfile) from forwardfiles

  output:
  set ( val(stdout), file($forwardfile) ) into forwardsets

  """
  echo $forwardfile | sed 's/.forward.fastq.gz//g'
  """
}

process mergereverselanes {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  file(inputdir) from inputdirectory2

  output:
  file('*.reverse.fastq.gz') into reversefiles

  """
  for samplename in `ls $inputdir | grep _L001_R2_001.fastq.gz | sed 's/_L001_R2_001.fastq.gz//g'`
  do
  cat $inputdir/\${samplename}_L001_R2_001.fastq.gz $inputdir/\${samplename}_L002_R2_001.fastq.gz $inputdir/\${samplename}_L003_R2_001.fastq.gz $inputdir/\${samplename}_L004_R2_001.fastq.gz > \${samplename}.reverse.fastq.gz
  done
  """
}

process reversesets {
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB

  input:
  file(reversefile) from reversefiles

  output:
  set ( val(stdout), file($reversefile) ) into reversesets

  """
  echo $reversefile | sed 's/.reverse.fastq.gz//g'
  """
}

fortrim = forwardsets.join(reversesets)

process docutadapt {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '8 GB'

  input:
  set ( sampleprefix, file(forwardraw), file(reverseraw) ) from fortrim

  output:
  set ( sampleprefix, file("${sampleprefix}_L001_R1_001.trimmed.fastq.gz"), file("${sampleprefix}_L001_R2_001.trimmed.fastq.gz") ) into (trimmingoutput1, trimmingoutput2)
  file("${sampleprefix}.trim.out") into trimouts

  """
  module load CUTAdapt/latest
  cutadapt -a file:${params.adapterfile} -A file:${params.adapterfile} -g file:${params.adapterfile} -G file:${params.adapterfile} -o ${sampleprefix}_L001_R1_001.trimmed.fastq.gz -p ${sampleprefix}_L001_R2_001.trimmed.fastq.gz ${forwardraw} ${reverseraw} -q 30,30 --minimum-length 50 --times 40 -e 0.1 --max-n 0 > ${sampleprefix}.trim.out 2> ${sampleprefix}.trim.err
  """
}

process dotrimlog {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  file "logdir/*" from trimouts.toSortedList()

  output:
  file("trim-summary.csv") into trimlogend

  """
  module load anaconda/Py2/python2
  python /home/AD/tbleazar/pipelines/trimlogger.py logdir trim-summary.csv
  """
}

process dotophat {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '24h'
  memory '100 GB'

  input:
  set ( sampleprefix, file(forwardtrimmed), file(reversetrimmed) ) from trimmingoutput1
  file refs from ref1.first()

  output:
  set ( sampleprefix, file("tophat_out/accepted_hits.bam") ) in bamfiles

  """
  module load TopHat/2.1.1
  module load Bowtie/latest
  module load SAMTools/latest
  referencebase=`echo ${params.referencegenome} | sed 's/.fasta//g' | sed 's/.fa//g'`
  gtfindexed=${params.genegtf}.indexed
  tophat2 -p ${params.cpus} --library-type fr-firststrand --mate-inner-dist 50 -G ${refs}/${params.genegtf} --transcriptome-index ${refs}/\${gtfindexed} \${referencebase} $forwardtrimmed $reversetrimmed
  """
}

process addreadgroups {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '12h'
  memory '8 GB'

  input:
  set ( sampleprefix, file(bamfile) ) from bamfiles

  output:
  set ( sampleprefix, file("${sampleprefix}.rg.bam") ) into readgrouped

  """
  module load GATK/4.1.3.0
  gatk AddOrReplaceReadGroups -I $bamfile -O ${sampleprefix}.rg.bam -LB LIBRARY -PL ILLUMINA -PU PLATFORMUNIT -SM ${sampleprefix}
  """
}

process indexbam {
  publishDir "$params.outdir/alignments", mode: "copy"
  cpus 32
  queue 'WORK'
  time '12h'
  memory '8 GB'

  input:
  set ( sampleprefix, file(rgfile) ) from readgrouped

  output:
  set ( sampleprefix, file("${rgfile}"), file("${rgfile}.bai") ) into bamindexed

  """
  module load SAMTools/latest
  samtools index $rgfile
  """
}

process dohtseq {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 32
  queue 'WORK'
  time '12h'
  memory '20 GB'

  input:
  set ( sampleprefix, file(bamfile), file(bamindex) ) from bamindexed
  file refs from ref2.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.htseq.out") ) into htseqcounted

  """
  module load HTSeq-count/0.11.2
  htseq-count -f bam -r pos -a 5 -s reverse -t exon -i gene_name -m union --nonunique=all $bamfile ${refs}/${params.genegtf} > ${sampleprefix}.htseq.out
  """
}



