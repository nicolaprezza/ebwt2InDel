# usage: snp2vcf.sh calls.snp ref.fasta output.vcf

# calls.snp: calls generated by ebwt2snp
# ref.fasta: reference fasta file
# output.vcf the output vcf file

# creates BWA index in the reference folder if needed
# converts calls.snp into a vcf file by aligning the calls on the reference genome using BWA MEM

# requires: seqtk, bwa

calls=$1
ref=$2
out=$3

echo "input calls: "$calls
echo "input reference: "$ref
echo "output: "$out

seqtk seq -F 'h' $calls > $out.fastq

if [ ! -f $ref.bwt ]; then
	bwa index $ref
fi

bwa mem $ref $out.fastq -o $out.sam

rm $out.fastq

sam2vcf -f $ref -s $out.sam -v $out

rm $out.sam

