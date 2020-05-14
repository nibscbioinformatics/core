#This is a script to filter a lofreq vcf file based on the alternative allele read frequency

import sys

vcffile = open(sys.argv[1])
fileout = open(sys.argv[2], "w")

freqcutoff = 0.5
depthcutoff = 100

for line in vcffile:
  if line[0] == "#":
    fileout.write(line)
    continue
  collect = line.rstrip().split("\t")
  infofield = collect[7].split(";")
  depth = int(infofield[0].split("=")[1])
  frequency = float(infofield[1].split("=")[1])
  if depth >= depthcutoff and frequency >= freqcutoff:
    fileout.write(line)
 
vcffile.close()
fileout.close()  
