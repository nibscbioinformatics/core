#This is a utility script for use with Nextflow pipelines to capture log files in a given folder
#Has two modes for now, cutadapt and flagstat


import os
import sys

trimsdir = sys.argv[1]
fileout = open(sys.argv[2], "w")
filemode = sys.argv[3]

basesin = {} #cutadapt stats
basesout = {}
pairsin = {}
pairsout = {}

rawfilenames = os.listdir(trimsdir)

if filemode == "cutadapt":
  #First the cutadapt log collection
  filenames = []
  for filename in rawfilenames:
    if ".trim.out" in filename:
      filenames.append(filename)
  filenames = sorted(filenames)
  print(filenames)
  for filename in filenames:
      basesin[filename] = 0
      basesout[filename] = 0
      pairsin[filename] = 0
      pairsout[filename] = 0
      filein = open(trimsdir+"/"+filename)
      for line in filein:
          if "Total read pairs processed:" in line:
              collect = line.rstrip().split()
              pairsin[filename] += int(collect[4].replace(",",""))
          if "Pairs written (passing filters):" in line:
              collect = line.rstrip().split()
              pairsout[filename] += int(collect[4].replace(",",""))
          if "Total basepairs processed:" in line:
              collect = line.rstrip().split()
              basesin[filename] += int(collect[3].replace(",",""))
          if "Total written (filtered):" in line:
              collect = line.rstrip().split()
              basesout[filename] += int(collect[3].replace(",",""))
      filein.close()
  #Now writing out this information in a table
  fileout.write("Sample-Name,Total-Read-Pairs-Sequenced,Total-Bases-Sequenced,Trimming-Survivor-Read-Pairs,Trimming-Survivor-Bases,Proportion-Reads-Surviving\n")
  for filename in filenames:
      fileout.write(filename.replace(".trim.out","")+","+str(pairsin[filename])+","+str(basesin[filename])+","+str(pairsout[filename])+","+str(basesout[filename])+","+str(float(pairsout[filename])/float(pairsin[filename]))+"\n")
  fileout.close()




