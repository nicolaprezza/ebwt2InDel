#Required: 
#HARC (https://github.com/shubhamchandak94/HARC)
#BCR_LCP_GSA (https://github.com/giovannarosone/BCR_LCP_GSA)


#run the ebwt2InDel pipeline in parallel
#usage: pebwt2InDel.sh input_fasta threads RAM read_len output_directory/ harc_folder/ 

#NOTE: use absolute paths ending with /

#breaks the reads into pieces of length "read_len" and runs the pipeline on "threads" parallel threads
#NOTE: "read_len" must be shorter than or equal to the length of most reads; remainders are lost!

#RAM is the RAM used by BCR (in MB)

START=$(date +%s.%N)

input_fasta=$1
p=$2
RAM=$3
len=$4
output_directory=$5
harc_folder=$6

if [[ "$p" -le 1 ]]; then p=2; fi #p must be >1


echo "Output directory = "$output_directory

#b
echo "Preparing input for HARC"
cat $input_fasta | grep -v ^'>' | fold -w $len | awk -v L="$len" 'length()==L {print}' | tr 'N' 'A' | awk '{printf "@a\n";print;printf "+\nH\n"}' > ${output_directory}input.fq 


#step2: sort reads using HARC
echo "Running HARC"
cd $harc_folder
/usr/bin/time -v ./harc -c ${output_directory}input.fq -t $p > ${output_directory}harc.compress.stdout 2>&1
/usr/bin/time -v ./harc -d ${output_directory}input.fq.harc -t $p > ${output_directory}harc.decompress.stdout 2>&1

#rm ${output_directory}input.fq.harc

sorted_file=${output_directory}input.fq.dna.d

echo "Split fasta"
#step3: split into p (or p-1) pieces

cd ${output_directory}
lines=$((`cat $sorted_file | wc -l`/(p-1)))
split -d -l $lines $sorted_file

#number of pieces (=number of threads)
N=`ls -l ${output_directory}x*  | wc -l`

echo "Processing in parallel all pieces"

foo () {

i=$1
if [[ "$i" -le 9 ]]; then i=0$i; fi
echo "processing piece" ${output_directory}x$i

mkdir ${output_directory}BCRx$i
mv ${output_directory}x$i ${output_directory}BCRx$i

#convert to fasta
cat ${output_directory}BCRx$i/x$i | awk '{printf ">a\n";print;}' > ${output_directory}BCRx$i/in.fa

rm ${output_directory}BCRx$i/x$i

#run BCR
cd ${output_directory}BCRx$i
echo "running BCR_LCP_GSA "${output_directory}BCRx$i/in.fa ${output_directory}BCRx$i/in $RAM
/usr/bin/time -v BCR_LCP_GSA ${output_directory}BCRx$i/in.fa ${output_directory}BCRx$i/in $RAM > ${output_directory}BCRx$i/BCR.stdout 2>&1

#run ebwt2InDel
/usr/bin/time -v ebwt2InDel -1 ${output_directory}BCRx$i/in.ebwt -o ${output_directory}BCRx$i/out.snp -m 3 > ${output_directory}BCRx$i/ebwt2InDel.stdout 2>&1

}

N=$((N-1))

for t in $(seq 0 $N); do foo $t & done; wait

echo "Concatenating outputs."

cat ${output_directory}BCRx*/out.snp > ${output_directory}/variants.snp

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)

for t in $(seq 0 $N); do  

 if [[ "$t" -le 9 ]]; then t=0$t; fi
 rm -rf ${output_directory}BCRx$t

done;

rm -rf ${output_directory}input*
rm -rf ${output_directory}harc*

echo "Done. Time (seconds): "$DIFF


