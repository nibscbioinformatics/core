//This is a nextflow script to collect all read files from a directory then merge if needed
//Then perform trimming using cutadapt
//Then follow dokraken.sh commands to run kraken and bracken on each sample

params.mergelanes = false
params.indir = "/usr/share/sequencing/projects/272/raw_data"
params.outdir = "/usr/share/sequencing/projects/272"
params.adapterfile = "/usr/share/sequencing/references/adapters/NexteraPE-PE.fa"
params.cpus = "8"
params.krakendb = "/usr/share/sequencing/references/metagenomes/standard_kraken_db/"

krakendbdir = Channel
  .fromPath(params.krakendb)

krakendbdir.into {
  dbdir1
  dbdir2
}

//merging four lanes only when params.mergelanes = true
if (params.mergelanes) {
  inputdirectory1 = file(params.indir)
  inputdirectory2 = file(params.indir)
} else {
  inputdirectory1 = Channel.empty()
  inputdirectory2 = Channel.empty()
}

process mergeforwardlanes {
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  file(inputdir) from inputdirectory1

  output:
  file('*.forward.fastq.gz') into forwardfiles

  when:
  params.mergelanes

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
  memory '4 GB'

  input:
  file(forwardfile) from forwardfiles.flatten()

  output:
  set ( val(sampleprefix), file("${sampleprefix}_L001_R1_001.fastq.gz") ) into forwardsets

  when:
  params.mergelanes

  script:
  sampleprefix = (forwardfile.name).replace(".forward.fastq.gz","")
  """
  basename=`echo $forwardfile | sed 's/.forward.fastq.gz//g'`
  mv $forwardfile \${basename}_L001_R1_001.fastq.gz
  """
}

process mergereverselanes {
  cpus 1
  queue 'WORK'
  time '12h'
  memory '4 GB'

  input:
  file(inputdir) from inputdirectory2

  output:
  file('*.reverse.fastq.gz') into reversefiles

  when:
  params.mergelanes

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
  memory '4 GB'

  input:
  file(reversefile) from reversefiles.flatten()

  output:
  set ( val(sampleprefix), file("${sampleprefix}_L001_R2_001.fastq.gz") ) into reversesets

  when:
  params.mergelanes

  script:
  sampleprefix = (reversefile.name).replace(".reverse.fastq.gz","")
  """
  basename=`echo $reversefile | sed 's/.reverse.fastq.gz//g'`
  mv $reversefile \${basename}_L001_R2_001.fastq.gz
  """
}
//end of merging section

//conditionally readschannel is either a join from the forward and reverse sets, or it is just the channel from filepairs
if (params.mergelanes) {
  readschannel = forwardsets.join(reversesets)
} else {
  filepattern = params.indir + "/*{_L001_R1_001,_L001_R2_001}.fastq.gz"
  Channel.fromFilePairs(filepattern, flat: true).set{readschannel}
}

process docutadapt {
  cpus 1
  queue 'WORK'
  time '12h'
  memory '8 GB'

  input:
  set ( sampleprefix, file(forwardraw), file(reverseraw) ) from readschannel

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

process dokraken {
  cpus 8
  queue 'WORK'
  time '24h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(forward), file(reverse) ) from trimmingoutput1
  file(krakendatabase) from dbdir1.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.kraken.mpa.report"), file("${sampleprefix}.kraken.report") ) into krakenoutput

  """
  module load kraken/1.1.1
  kraken --threads ${params.cpus} --paired --preload --db $krakendatabase --classified-out ${sampleprefix}.kraken.classified --unclassified-out ${sampleprefix}.kraken.unclassified --output ${sampleprefix}.kraken.output $forward $reverse
  kraken-mpa-report --db $krakendatabase ${sampleprefix}.kraken.output > ${sampleprefix}.kraken.mpa.report
  kraken-report --db $krakendatabase ${sampleprefix}.kraken.output > ${sampleprefix}.kraken.report
  """
}

process dobracken {
  publishDir "$params.outdir/analysis", mode: "copy"
  cpus 8
  queue 'WORK'
  time '24h'
  memory '240 GB'

  input:
  set ( sampleprefix, file(mpareport), file(krakenreport) ) from krakenoutput
  file(krakendatabase) from dbdir2.first()

  output:
  set ( sampleprefix, file("${sampleprefix}.bracken.report") ) into brackenoutput

  """
  module load bracken/2.5
  bracken -d $krakendatabase -i $krakenreport -o ${sampleprefix}.bracken.report -r 140 -t 5
  """
}
