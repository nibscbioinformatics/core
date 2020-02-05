#!/opt/software/conda2/envs/NextFlow/bin/nextflow


samples_ch = Channel.value([["sample1", "/path/to/sample1"],["sample2", "/path/tooe/sample2"]])


process testTuples {

  input:
  val(testInput) from samples_ch.collect()

  exec:
  ids = []
  paths = []
  testInput.each() { a,b -> ids.add(a); paths.add(b) }
  idall = ids.join(",")
  pathall = paths.join(",")
  println "now printing the sample"
  println ".......\n"
  println "-------------------\n"
  print "$ids\n"
  println "-------------------\n"
  print "$paths\n"
  println "-------------------\n"
  println "collated ids \n"
  println "$idall\n"
  println "collated paths \n"
  println "$pathall\n"
  println ".......\n"
  println ".......\n"

}
