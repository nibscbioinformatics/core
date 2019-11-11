#This is a script to write a bash script with  commands to loop running Kraken over all fastq.gz files in a given directory
#Usage:
#python loopkraken.py /path/to/directory/which/has/fastq.gz/files/inside /path/to/desired/results/directory > nameofbashscripttowrite.sh
#Output:
#The script will write a script with qsub commands to write kraken reports to subdirectories inside the specified results directory
#After running loopkraken.py as above, logon to hpc-head, then:
#chmod u+x nameofbashscripttowrite.sh
#./nameofbashscripttowrite.sh
#After the submitted jobs finish, a report file will be found at merged.kraken.report.tsv in the desired results directory

import sys
import os

print("#!/bin/bash")

fastqdir = sys.argv[1]
resultsdir = sys.argv[2]

forwardfiles = []
for filename in os.listdir(fastqdir):
  if "_R1_001.fastq.gz" in filename:
    forwardfiles.append(filename)

for filename in forwardfiles:
  targetdir = resultsdir+"/"+(filename.replace("_R1_001.fastq.gz",""))
  reversefile = filename.replace("_R1_001.fastq.gz","_R2_001.fastq.gz")
  print("mkdir -p "+targetdir)
  print('/usr/local/mark/dokraken.sh '+fastqdir+'/'+filename+' '+fastqdir+'/'+reversefile+' '+targetdir)

#Now write a script to run kraken-mpa-report commands in the resultsdir
if not os.path.exists(resultsdir):
  os.makedirs(resultsdir)
fileout = open(resultsdir+"/reportingscript.sh","w")
fileout.write("#!/bin/bash\n")
fileout.write("/usr/local/bin/kraken-mpa-report --header-line --db /raid/kraken/standard_kraken_db/")
for filename in forwardfiles:
  targetdir = resultsdir+"/"+(filename.replace("_R1_001.fastq.gz",""))
  fileout.write(" "+targetdir+"/kraken.output")
fileout.write(" > "+resultsdir+"/merged.kraken.report.tsv\n")
fileout.close()

#holding this job until all the kraken runs finish
print("chmod u+x "+resultsdir+"/reportingscript.sh")
print(resultsdir+'/reportingscript.sh')
