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

#Having run the trimming and Trinity commands, now generating summary stats for the assemblies

#Now run stats on the trinity assemblies
for i in range(len(cutnames)):
  cutname = cutnames[i]
  contigs = processeddir+"/"+cutname+"_trinity/Trinity.fasta"
  print("/usr/local/bin/trinityrnaseq-Trinity-v2.8.3/util/TrinityStats.pl "+contigs+" > "+processeddir+"/"+cutname+"_trinity/assemblystats.txt")
  


