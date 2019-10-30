fastq=$1 ## assumes its fastq.gz
fasta=${fastq%fastq.gz}.fasta

echo "processing $fastq into $fasta"
paste - - - - < $fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $fasta
echo "completed"
echo
