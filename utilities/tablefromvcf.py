#This is a script to read a LoFreq VCF file and extract a CSV table with Sample,Chromosome,Position,Ref,Alt,Ref_Reads,Alt_Reads,Proportion,Basic_Pass

import sys

filein = open(sys.argv[1])
fileout = open(sys.argv[2], "w")

#VCF files have lines like:
###CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#7       55259288        .       TC      T       152     PASS    DP=32327;AF=0.000588;SB=0;DP4=32308,0,19,0;INDEL;HRUN=1
#7       55259290        .       GC      G       239     PASS    DP=32327;AF=0.000804;SB=0;DP4=32300,0,26,0;INDEL;HRUN=2
#7       55259293        .       A       AG      58      PASS    DP=32327;AF=0.000309;SB=0;DP4=32315,0,10,0;INDEL;HRUN=1
#7       55259294        .       GC      G       463     PASS    DP=32327;AF=0.001299;SB=0;DP4=32284,0,42,0;INDEL;HRUN=2
#7       55259294        .       G       T       12517   PASS    DP=32327;AF=0.045040;SB=0;DP4=30838,0,1456,0

fileout.write("Chromosome,Position,Ref,Alt,Ref_Reads,Alt_Reads,Proportion,Basic_Pass\n")
for line in filein:
  if line[0] == "#":
    continue
  collect = line.rstrip().split("\t")
  chromosome = collect[0]
  position = collect[1]
  ref = collect[3]
  alt = collect[4]
  infofield = collect[7].split(";")
  proportion = infofield[1].split("=")[1]
  refreads = int(infofield[3].split("=")[1].split(",")[0]) + int(infofield[3].split("=")[1].split(",")[1])
  altreads = int(infofield[3].split("=")[1].split(",")[2]) + int(infofield[3].split("=")[1].split(",")[3])
  basicpass = altreads > 100 and float(proportion) > 0.01
  fileout.write(chromosome+","+position+","+ref+","+alt+","+str(refreads)+","+str(altreads)+","+proportion+","+str(basicpass)+"\n")
filein.close()
fileout.close()
  
  
