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

#After having run the trimming, now write the following commands

#Now concatenate the lanes and make Trinity folders for them, and delete the single unmerged lanes
for i in range(len(cutnames)):
  cutname = cutnames[i]
  print("mkdir -p "+processeddir+"/"+cutname+"_trinity")
  formerge = ["_L001.p1.fastq.gz","_L002.p1.fastq.gz","_L003.p1.fastq.gz","_L004.p1.fastq.gz"]
  revmerge = ["_L001.p2.fastq.gz","_L002.p2.fastq.gz","_L003.p2.fastq.gz","_L004.p2.fastq.gz"]
  for j in range(len(formerge)):
    formerge[j] = trimmeddir+cutname+formerge[j]
  for j in range(len(revmerge)):
    revmerge[j] = trimmeddir+cutname+revmerge[j]
  print("cat "+formerge[0]+" "+formerge[1]+" "+formerge[2]+" "+formerge[3]+" > "+trimmeddir+cutname+".p1.fastq.gz")
  print("cat "+revmerge[0]+" "+revmerge[1]+" "+revmerge[2]+" "+revmerge[3]+" > "+trimmeddir+cutname+".p2.fastq.gz")

#Now run the actual Trinity commands with qsub commands
for i in range(len(cutnames)):
  cutname = cutnames[i]
  logout = logdir+cutname+".trinity.out"
  logerr = logdir+cutname+".trinity.err"
  print("qsub -N qjob_"+projectnumber+"_trinity -pe parallel 32 -cwd -S /bin/bash -o "+logout+" -e "+logerr+" "+runTrinity+" "+trimmeddir+cutname+".p1.fastq.gz "+trimmeddir+cutname+".p2.fastq.gz "+processeddir+"/"+cutname+"_trinity")


  


