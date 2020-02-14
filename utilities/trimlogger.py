#For now just works on the TopHat2 pipeline
#This is a script to go through all of the cutadapt and all of the tophat2 logs to grab trimming and alignment stats
#Log files for cutadapt are under the log folder in files like 220_NSLAT_lentivirus_transduced_1_S10_L001_R1_001.fastq.gz.cutadapt.log
#Log files for tophap2 are under the aligned folder in align_summary.txt files

import os
import sys

trimsdir = sys.argv[1]
fileout = open(sys.argv[2], "w")

#leftin = {} #tophat2 stats
#leftmapped = {}
#rightin = {}
#rightmapped = {}
#pairsmapped = {}
#pairsmultiply = {}
basesin = {} #cutadapt stats
basesout = {}
pairsin = {}
pairsout = {}

filenames = os.listdir(trimsdir)
filenames = sorted(filenames)
print(filenames)

#First the cutadapt logs
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

#print(pairsin)
#print(pairsout)
#print(basesin)
#print(basesout)

#now the tophat2 logs
#for samplename in samplenames:
#    filein = open(working+"/aligned/"+samplename+"/align_summary.txt")
#    left = filein.readline()
#    leftin[samplename] = filein.readline().rstrip().split()[2]
#    leftmapped[samplename] = filein.readline().rstrip().split()[2]
#    junk = filein.readline()
#    right = filein.readline()
#    rightin[samplename] = filein.readline().rstrip().split()[2]
#    rightmapped[samplename] = filein.readline().rstrip().split()[2]
#    junk = filein.readline()
#    junk = filein.readline()
#    junk = filein.readline()
#    pairsmapped[samplename] = filein.readline().rstrip().split()[2]
#    pairsmultiply[samplename] = filein.readline().rstrip().split()[2]
#    filein.close()

#print(leftin)
#print(pairsmultiply)

#Now writing out this information in a table
fileout.write("Sample-Name,Total-Read-Pairs-Sequenced,Total-Bases-Sequenced,Trimming-Survivor-Read-Pairs,Trimming-Survivor-Bases,Proportion-Reads-Surviving\n")
for filename in filenames:
    fileout.write(filename.replace(".trim.out","")+","+str(pairsin[filename])+","+str(basesin[filename])+","+str(pairsout[filename])+","+str(basesout[filename])+","+str(float(pairsout[filename])/float(pairsin[filename]))+"\n")
fileout.close()




