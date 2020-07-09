#!/bin/bash
#SBATCH -p WORK # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH -c 8 # number of cpus per task
#SBATCH --mem 40 # memory pool for all cores in Gb
#SBATCH -t 90:30:30 # max time (HH:MM:SS)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

echo "### START - the job is starting at"
date
starttime=`date +"%s"`
echo
echo "the job is running on the node $SLURM_NODELIST"
echo "job number $SLURM_JOB_ID"
echo "STAT:jobName:$SLURM_JOB_ID"
echo "STAT:exechosts:$SLURM_NODELIST"
echo

cd $PWD
module load Bowtie/latest
module load BWA/latest
module load SAMTools/latest
module load GATK/4.1.3.0

reffasta=$1
refdir=$2

cd $refdir
refroot=`echo $reffasta | sed 's/.fasta$//g' | sed 's/.fa$//g'`
bowtie2-build $reffasta $reffasta
bwa index $reffasta
samtools faidx $reffasta
gatk --java-options "-Xmx20g" CreateSequenceDictionary --REFERENCE $reffasta --OUTPUT ${refroot}.dict

echo "####END job finished"
endtime=`date +"%s"`
duration=$((endtime - starttime))
echo "STAT:startTime:$starttime"
echo "STAT:doneTime:$endtime"
echo "STAT:runtime:$duration"
