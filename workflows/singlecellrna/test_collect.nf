#!/opt/software/conda2/envs/NextFlow/bin/nextflow


samples_ch = Channel.value([["sample1", "/path/to/sample1"],["sample2", "/path/tooe/sample2"]])


process testTuples {

  input:
  val(testInput) from samples_ch.collect()

  exec:
  println "now printing the sample"
  println "......."
  println "-------------------"
  print "$testInput"
  println "-------------------"
  println "......."
  println "......."

}
