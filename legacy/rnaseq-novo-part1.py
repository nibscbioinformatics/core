#This is a script to take RNA-seq data and perform processing with it to generate de novo transcriptome assembly

adapter = "/sequencing/projects/193/TruSeq-adapters-recommended.txt"
threads = "32"
runcutadapt = "/sequencing/projects/193/runcutadapt.sh"
runTrinity = "/sequencing/projects/193/runTrinity.sh"

import os
import sys

projectdir = sys.argv[1] #/sequencing/projects/193/
projectnumber = sys.argv[2] #193
rawdir = sys.argv[3] #/sequencing/nextseq/processed/180817/fastq/raw/

#Start writing to stdout
print("#!/bin/bash")
print("cd "+projectdir)

trimmeddir = projectdir+"/trimmed/"
processeddir = projectdir+"/processed/"
logdir = projectdir+"/log/"
print("mkdir -p "+trimmeddir)
print("mkdir -p "+processeddir)
print("mkdir -p "+logdir)

#First get a definitively ordered collection of the raw fastq files and also of their lane-merged names
forwardfiles = []
reversefiles = []
cutnames = []
for filename in os.listdir(rawdir):
  if filename[:len(projectnumber)] != projectnumber:
    continue #we don't want to handle any files that don't match our project number
  if "R1_001.fastq.gz" in filename:
    forwardfiles.append(filename)
  if "R2_001.fastq.gz" in filename:
    reversefiles.append(filename)
  if "_L004_R2_001.fastq.gz" in filename:
    cutnames.append(filename.replace("_L004_R2_001.fastq.gz",""))
forwardfiles.sort()
reversefiles.sort()
cutnames.sort()

#Now write trimming qsub commands with cutadapt
for i in range(len(forwardfiles)):
  forwardfile = rawdir+"/"+forwardfiles[i]
  reversefile = forwardfile.replace("R1_001.fastq.gz","R2_001.fastq.gz")
  pone = trimmeddir+forwardfiles[i].replace("_R1_001.fastq.gz",".p1.fastq.gz")
  ptwo = trimmeddir+forwardfiles[i].replace("_R1_001.fastq.gz",".p2.fastq.gz")
  logout = logdir+forwardfiles[i].replace("_R1_001.fastq.gz",".trim.out")
  logerr = logdir+forwardfiles[i].replace("_R1_001.fastq.gz",".trim.err")
  print("qsub -N qjob_"+projectnumber+"_trims -pe parallel 1 -cwd -S /bin/bash -o "+logout+" -e "+logerr+" "+runcutadapt+" "+adapter+" "+forwardfile+" "+reversefile+" "+pone+" "+ptwo)


