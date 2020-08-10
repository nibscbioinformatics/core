cd /usr/share/sequencing/nextseq/output
for nextseqfolder in `ls`
do
seqdate=`echo $nextseqfolder | cut -c1-6`
if test -f "$nextseqfolder/RunCompletionStatus.xml"; then
if ! test -e "/usr/share/sequencing/nextseq/processed/$seqdate"; then
echo "RUNNING FOR FOLDER $nextseqfolder"
mkdir -p /usr/share/sequencing/nextseq/processed/$seqdate
sbatch --output="/usr/share/sequencing/nextseq/processed/$seqdate/bcltrigger.out" --error="/usr/share/sequencing/nextseq/processed/$seqdate/bcltrigger.err" /home/AD/tbleazar/CODE/core/utilities/dobcl2fastq.sh $seqdate /usr/share/sequencing/nextseq/output/$nextseqfolder
fi
fi
done
