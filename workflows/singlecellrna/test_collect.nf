#!/opt/software/conda2/envs/NextFlow/bin/nextflow


samples_ch = Channel.value([["sample1", "/path/to/sample1"],["sample2", "/path/tooe/sample2"]])


process testTuples {

  input:
  tuple sample, path from samples_ch.collect()

  exec:
  println "$sample"
  println "$path"

}
