#!/bin/bash -l
#SBATCH -p WORK # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH -c 16 # number of cpus per task
#SBATCH --mem 100 # memory pool for all cores in Gb
#SBATCH -t 23:30:30 # max time (HH:MM:SS)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

echo "### START - the job is starting at"
date
starttime=`date +"%s"`
echo
echo "the job is running on the node $SLURM_NODELIST"
echo "job number $SLURM_JOB_ID"
echo "STAT:jobName:$SLURM_JOB_ID\.out"
echo "STAT:exechosts:$SLURM_NODELIST"
echo

cd $PWD

sequencedate=$1 #200122
diroutput=$2 # /usr/share/sequencing/nextseq/output/200122_NB501506_0046_AH3CYJAFX2/

mkdir -p /usr/share/sequencing/nextseq/processed/${sequencedate}/InterOp
mkdir -p /usr/share/sequencing/nextseq/processed/${sequencedate}/log
mkdir -p /usr/share/sequencing/nextseq/processed/${sequencedate}/bcl2fastq
mkdir -p /usr/share/sequencing/nextseq/processed/${sequencedate}/fastq

module load bcl2fastq/2.20

#Make sure that the appropriate sample sheet is at:
#/usr/share/sequencing/nextseq/processed/samplesheets/${sequencedate}_SampleSheet.csv

bcl2fastq --barcode-mismatches 0 -R $diroutput --interop-dir /usr/share/sequencing/nextseq/processed/${sequencedate}/InterOp --sample-sheet $diroutput/SampleSheet.csv -o /usr/share/sequencing/nextseq/processed/${sequencedate}/bcl2fastq --no-lane-splitting > /usr/share/sequencing/nextseq/processed/${sequencedate}/log/bcl2fastq.out 2> /usr/share/sequencing/nextseq/processed/${sequencedate}/log/bcl2fastq.err

#Make symlinks to the fastq files into the folder fastq (so they are all in the same folder together for ease of use)
cd /usr/share/sequencing/nextseq/processed/${sequencedate}/fastq
for s in `find ../bcl2fastq -name "*.fastq.gz"` ; do ln -s $s ; done

chown -R tbleazar.hpc_sequencing_writers /usr/share/sequencing/nextseq/processed/${sequencedate}
chmod -R 755 /usr/share/sequencing/nextseq/processed/${sequencedate}

echo "####END job finished"
endtime=`date +"%s"`
duration=$((endtime - starttime))
echo "STAT:startTime:$starttime"
echo "STAT:doneTime:$endtime"
echo "STAT:runtime:$duration"
