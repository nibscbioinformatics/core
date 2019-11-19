#This is a python script to go through a fastq.gz file and look at the read names, then throw out those from offending tiles
#Read names are of the form (for tile 1101 here):
#@M01745:221:000000000-C5LDC:1:1101:16763:1014 1:N:0:1

import gzip
import os

def clearfile(rawfile, outfile):
    goodtile = True
    filein = gzip.open(rawfile, "rt")
    fileout = gzip.open(outfile, "wt")
    linecount = 0
    for line in filein:
        linecount += 1
        if linecount % 4 == 1:
            tile = line.rstrip().split(":")[4]
            #print(line)
            #print(tile)
            if int(tile) == 2105 or int(tile) == 2111:
                goodtile = False
            else:
                goodtile = True
        if goodtile:
            fileout.write(line)

rawdir = "/usr/share/sequencing/miseq/output/191115_M01745_0254_000000000-CP2V4/Data/Intensities/BaseCalls/"
outdir = "/usr/share/sequencing/projects/293/raw_data/"

forwardfiles = []
for filename in os.listdir(rawdir):
    if "293-" in filename and "L001_R1_001.fastq.gz" in filename:
        forwardfiles.append(filename)

for filename in forwardfiles:
    print(filename)
    clearfile(rawdir+filename, outdir+filename)
    print(filename.replace("L001_R1_001.fastq.gz","L001_R2_001.fastq.gz"))
    clearfile(rawdir+(filename.replace("L001_R1_001.fastq.gz","L001_R2_001.fastq.gz")), outdir+(filename.replace("L001_R1_001.fastq.gz","L001_R2_001.fastq.gz")))
