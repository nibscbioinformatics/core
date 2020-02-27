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

#The flagstat files are like:
#514693 + 0 in total (QC-passed reads + QC-failed reads)   TOTAL IN READS
#0 + 0 secondary                                           MULTIPLE ALIGNMENTS AFTER THE FIRST FOR A READ
#35 + 0 supplementary                                      CHIMERIC SPLIT ALIGNMENTS AFTER THE FIRST REPRESENTATIVE
#27597 + 0 duplicates                                      MARKED AS DUPLICATES
#77753 + 0 mapped (15.11% : N/A)                           TOTAL = MATEMAPPEDREADS + SINGLETONS + SECONDARY + SUPPLEMENTARY
#514658 + 0 paired in sequencing                           ACTUAL READ PAIRS AFTER TRIMMING?
#257329 + 0 read1
#257329 + 0 read2
#69396 + 0 properly paired (13.48% : N/A)                  READ IN A PAIR ALIGNED PROPERLY TOGETHER
#69400 + 0 with itself and mate mapped                     READ THAT HAS ITSELF AND MATE MAPPED
#8318 + 0 singletons (1.62% : N/A)                         READ WITH ITS MATE UNMAPPED
#0 + 0 with mate mapped to a different chr
#0 + 0 with mate mapped to a different chr (mapQ>=5)
if filemode == "flagstat":
  mappedinpair = {}
  mappedsingletons = {}
  readsin = {}
  markedduplicates = {}
  secondary = {}
  supplementary = {}
  #Collect all the flagstat.out files and make a table from them
  filenames = []
  for filename in rawfilenames:
    if ".flagstat.out" in filename:
      filenames.append(filename)
  filenames = sorted(filenames)
  print(filenames)
  for filename in filenames:
    filein = open(trimsdir+"/"+filename)
    for line in filein:
      if "with itself and mate mapped" in line:
        mappedinpair[filename] = int(line.rstrip().split()[0])
      if "singletons" in line:
        mappedsingletons[filename] = int(line.rstrip().split()[0])
      if "in total (QC-passed reads + QC-failed reads)" in line:
        readsin[filename] = int(line.rstrip().split()[0])
      if "duplicates" in line:
        markedduplicates[filename] = int(line.rstrip().split()[0])
      if "secondary" in line:
        secondary[filename] = int(line.rstrip().split()[0])
      if "supplementary" in line:
        supplementary[filename] = int(line.rstrip().split()[0])
    filein.close()
  #Now writing out this information
  fileout.write("Sample_Name,Reads_Processed,Aligned_Reads,Secondary_Alignments,Supplementary_Alignments,Singleton_Alignments,Marked_Duplicates\n")
  for filename in filenames:
    fileout.write(",".join([filename.replace(".flagstat.out",""), str(readsin[filename]), str(mappedinpair[filename] + mappedsingletons[filename]), str(secondary[filename]), str(supplementary[filename]), str(mappedsingletons[filename]), str(markedduplicates[filename])]))
    fileout.write("\n")
  fileout.close()








