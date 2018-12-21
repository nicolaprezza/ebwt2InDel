# Usage: pipeline.sh <reads1.fastq> <reads2.fastq> <reference.fasta>
#
# behaviour (the steps are skipped if the file they will produce already exists)
# 1.  Converts input reads to fasta -> reads1.fasta reads2.fasta 
# 2.  Adds the reverse complements to the reads and concatenates the two read files -> reads1.reads2.frc.fasta
# 3.  Builds EGSA -> reads1.reads2.frc.fasta.gesa
# 4.  Run ebwt2clust -> reads1.reads2.frc.fasta.clusters 
# 5.  Run clust2snp -> reads1.reads2.frc.snp
# 6.  Run snp2fastq -> reads1.reads2.frc.snp.fastq
# 7.  Builds BWA MEM index of reference.fasta -> reference.fasta.{amb,ann,bwt,fai,pac,sa} files
# 8.  Create reference of reads1.fasta using BWA MEM + bcftools + vcfconsensus -> reads1.reference.fasta
# 9.  Builds BWA MEM index of reads1.reference.fasta -> reads1.reference.fasta.{amb,ann,bwt,fai,pac,sa} files
# 10. Aligns reads1.reads2.frc.snp.fastq on reads1.reference.fasta -> reads1.reads2.frc.snp.sam
# 11. Generates VCF (ebwt2snp calls) using sam2vcf -> reads1.reads2.frc.snp.sam.vcf. Filter keeping only variants testified by at least <c> reads -> reads1.reads2.frc.<c>.snp.sam.vcf
# 12. Generates VCF (bcftools calls) using BWA MEM + bcftools -> reads1.reads2.bcftools.vcf
# 13. Generate report containing precision/recall of ebwt2snp pipeline (using BWA+bcftools pipeline as ground truth)

# Requires the following executables to be globally visible (in addition to the executables of ebwt2snp): 
# - fastq2fasta.sh (https://github.com/nicolaprezza/bioinfo-tools)
# - vcf2fasta.sh (https://github.com/nicolaprezza/bioinfo-tools)
# - seqtk
# - egsa (https://github.com/felipelouza/egsa)
# - bwa
# - samtools
# - bcftools
# - bgzip
# - vcftools

#Byte size of EGSA fields
GSAtext=4
GSAsuff=1
LCP=1

WD=$(dirname $(readlink -f "$1"))

#extract full path of working directory without trailing /
echo "Working directory: "${WD}

#extract file name and extension of reads1
READS1=`basename $1`
EXT1=${READS1##*.}
READS1="${READS1%.*}"

echo "input 1: "${READS1}"."${EXT1}

#extract file name and extension of reads2
READS2=`basename $2`
EXT2=${READS2##*.}
READS2="${READS2%.*}"

#extract file name of reference
REF=`basename $3`

#parameter -c in clust2snp (max reads in consensus. Ideally big, but bigger slows down computation)
#C=$4

#parameter -m in clust2snp (we keep this fixed: at least 4 reads per individual in cluster)
M=4

TIME_EGSA=${WD}/egsa.time
TIME_EBWTCLUST=${WD}/ebwt2clust.time
TIME_CLUST2SNP=${WD}/clust2snp.time
TIME_BWAMEM=${WD}/bwamem.time
TIME_BWAIDX=${WD}/bwaindex.time
TIME_BCFTOOLS=${WD}/bcftools.time
REPORT_OUT=${WD}/report

samtool_version=`samtools 2>&1 >/dev/null | grep Version | cut -d' ' -f 2 | cut -d. -f 1`

echo "input 2: "${READS2}"."${EXT2}

# 1.  Converts input reads to fasta, if not already in fasta -> reads1.fasta reads2.fasta 

if [ ! -f ${WD}/${READS1}.fasta ]; then
	echo ${WD}/${READS1}.fasta" not found. Converting fastq to fasta ..."
	fastq2fasta.sh ${WD}/${READS1}.fastq > ${WD}/${READS1}.fasta
fi

if [ ! -f ${WD}/${READS2}.fasta ]; then
	echo ${WD}/${READS2}.fasta" not found. Converting fastq to fasta ..."
	fastq2fasta.sh ${WD}/${READS2}.fastq > ${WD}/${READS2}.fasta
fi

#number of reads in individual 1
N=`cat ${WD}/${READS1}.fasta | grep '^>' | wc -l`
N=$((N*2))

# 2.  Adds the reverse complements to the reads and concatenates the two read files -> reads1.reads2.frc.fasta

if [ ! -f ${WD}/${READS1}.${READS2}.frc.fasta ]; then
	echo "Adding reverse complement to the reads and building main fasta "${WD}/${READS1}.${READS2}.frc.fasta" ..."
	seqtk seq -r ${WD}/${READS1}.fasta > ${WD}/${READS1}.rc.fasta
	seqtk seq -r ${WD}/${READS2}.fasta > ${WD}/${READS2}.rc.fasta
	cat ${WD}/${READS1}.fasta ${WD}/${READS1}.rc.fasta ${WD}/${READS2}.fasta ${WD}/${READS2}.rc.fasta > ${WD}/${READS1}.${READS2}.frc.fasta
	rm ${WD}/${READS1}.rc.fasta ${WD}/${READS2}.rc.fasta
fi

# 3.  If EGSA/BCR files do not exist (both), then builds EGSA -> reads1.reads2.frc.fasta.gesa
if [ ! -f ${WD}/${READS1}.${READS2}.frc.fasta.gesa ]; then
	if [ ! -f ${WD}/${READS1}.${READS2}.frc.fasta.out ]; then
		echo "building EGSA ..."
		/usr/bin/time -v egsa -vvv ${WD}/${READS1}.${READS2}.frc.fasta 0 > ${TIME_EGSA} 2>&1
		#/usr/bin/time -v gsufsort ${WD}/${READS1}.${READS2}.frc.fasta --gesa ${GSAtext} ${GSAsuff} ${LCP} > ${TIME_EGSA} 2>&1
		mv ${WD}/${READS1}.${READS2}.frc.fasta.${GSAtext}.${GSAsuff}.${LCP}.1.gesa ${WD}/${READS1}.${READS2}.frc.fasta.gesa
		
		rm -rf ${WD}/tmp
		rm -rf ${WD}/partition
	fi
fi

# 4.  Run ebwt2clust -> reads1.reads2.frc.fasta.clusters 

if [ ! -f ${WD}/${READS1}.${READS2}.frc.fasta.clusters ]; then
	echo "running ebwt2clust ..."
	/usr/bin/time -v ebwt2clust -i ${WD}/${READS1}.${READS2}.frc.fasta -m $((M*2)) -x ${LCP} -y ${GSAtext} -z ${GSAsuff} > ${TIME_EBWTCLUST} 2>&1
fi

# 5.  Run clust2snp -> reads1.reads2.frc.snp

if [ ! -f ${WD}/${READS1}.${READS2}.frc.snp ]; then
	echo "running clust2snp ..."
	/usr/bin/time -v clust2snp -i ${WD}/${READS1}.${READS2}.frc.fasta -n $N -x ${LCP} -y ${GSAtext} -z ${GSAsuff} > ${TIME_CLUST2SNP} 2>&1
fi

# 6.  Run snp2fastq -> reads1.reads2.frc.snp.fastq

#if [ ! -f ${WD}/${READS1}.${READS2}.frc.snp.fastq ]; then
#	echo "converting .snp to .fastq ..."
#	snp2fastq ${WD}/${READS1}.${READS2}.frc.snp
#fi


# 7.  Builds BWA MEM index of reference.fasta -> reference.fasta.{amb,ann,bwt,fai,pac,sa} files

if [ ! -f ${WD}/${REF}.amb ]; then
	echo "Indexing "${WD}/${REF}" (BWA) ..."
	bwa index ${WD}/${REF}
fi

# 8.  Create reference of reads1.fasta using BWA MEM + bcftools + vcfconsensus

if [ "$samtool_version" -eq "1" ]; then
	if [ ! -f ${WD}/${READS1}.reference.fasta ]; then
		echo "Computing reference of "${READS1}". Storing it to file "${WD}/${READS1}.reference.fasta" ..."
		bwa mem ${WD}/${REF} ${WD}/${READS1}.fastq -o ${WD}/alignment.tmp.sam
		samtools view -b -S ${WD}/alignment.tmp.sam > ${WD}/alignment.tmp.bam
		samtools sort ${WD}/alignment.tmp.bam > ${WD}/alignment.tmp.sorted.bam
		samtools index ${WD}/alignment.tmp.sorted.bam
		bcftools mpileup -f ${WD}/${REF} ${WD}/alignment.tmp.sorted.bam | bcftools call -mv -o ${WD}/calls.tmp.vcf

		#filter VCF
		cat ${WD}/calls.tmp.vcf | grep '^#' > ${WD}/H
		cat ${WD}/calls.tmp.vcf | grep -v ^# |  awk '$6>14' > ${WD}/B
		cat ${WD}/H ${WD}/B > ${WD}/calls_filtered.tmp.vcf
		mv ${WD}/calls.tmp.vcf ${WD}/calls_used_to_build_reference_unfiltered.vcf
		bgzip ${WD}/calls_filtered.tmp.vcf
		rm ${WD}/H ${WD}/B

		vcf2fasta.sh ${WD}/calls_filtered.tmp.vcf.gz ${WD}/${REF} > ${WD}/${READS1}.reference.fasta
		mv ${WD}/calls_filtered.tmp.vcf.gz ${WD}/calls_used_to_build_reference.vcf.gz
		rm ${WD}/*.tmp*
	fi
else
	if [ ! -f ${WD}/${READS1}.reference.fasta ]; then
		echo "Computing reference of "${READS1}". Storing it to file "${WD}/${READS1}.reference.fasta" ..."
		bwa mem ${WD}/${REF} ${WD}/${READS1}.fastq -o ${WD}/alignment.tmp.sam
		samtools view -b -S ${WD}/alignment.tmp.sam > ${WD}/alignment.tmp.bam
		samtools sort ${WD}/alignment.tmp.bam ${WD}/alignment.tmp.sorted
		samtools index ${WD}/alignment.tmp.sorted.bam

		#OLD SAMTOOLS/BCFTOOLS
		samtools mpileup -uD -f ${WD}/${REF} ${WD}/alignment.tmp.sorted.bam > ${WD}/mpileup.tmp
		bcftools view -vc ${WD}/mpileup.tmp > ${WD}/calls.tmp.vcf
		
		#filter VCF
		cat ${WD}/calls.tmp.vcf | grep '^#' > ${WD}/H
		cat ${WD}/calls.tmp.vcf | grep -v ^# |  awk '$6>14' > ${WD}/B
		cat ${WD}/H ${WD}/B > ${WD}/calls_filtered.tmp.vcf
		mv ${WD}/calls.tmp.vcf ${WD}/calls_used_to_build_reference_unfiltered.vcf
		bgzip ${WD}/calls_filtered.tmp.vcf
		rm ${WD}/H ${WD}/B

		vcf2fasta.sh ${WD}/calls_filtered.tmp.vcf.gz ${WD}/${REF} > ${WD}/${READS1}.reference.fasta
		mv ${WD}/calls_filtered.tmp.vcf.gz ${WD}/calls_used_to_build_reference.vcf.gz
		rm ${WD}/*.tmp*
	fi
fi

# 9.  Builds BWA MEM index of reads1.reference.fasta -> reads1.reference.fasta.{amb,ann,bwt,fai,pac,sa} files

if [ ! -f ${WD}/${READS1}.reference.fasta.amb ]; then
	echo "Indexing "${WD}/${READS1}.reference.fasta" (BWA) ..."
	/usr/bin/time -v bwa index ${WD}/${READS1}.reference.fasta > ${TIME_BWAIDX} 2>&1
fi

# 10.  Aligns reads1.reads2.frc.snp.fastq on reads1.reference.fasta -> reads1.reads2.frc.snp.sam

#if [ ! -f ${WD}/${READS1}.${READS2}.frc.snp.sam ]; then
#	echo "Aligning "${WD}/${READS1}.${READS2}.frc.snp.fastq" on "${WD}/${READS1}.reference.fasta" ..."
#	bwa mem ${WD}/${READS1}.reference.fasta ${WD}/${READS1}.${READS2}.frc.snp.fastq -o ${WD}/${READS1}.${READS2}.frc.snp.sam 
#fi

# 11. Generates VCF (ebwt2snp calls) using sam2vcf -> reads1.reads2.frc.snp.sam.vcf

#if [ ! -f ${WD}/${READS1}.${READS2}.frc.${C}.snp.sam.vcf ]; then
#	echo "Generating ebwt2snp's VCF in "${WD}/${READS1}.${READS2}.frc.snp.sam.vcf" ..."
#	sam2vcf -s ${WD}/${READS1}.${READS2}.frc.snp.sam

	#Filter the VCF: keep only variations testified by at least C reads (in both REF and ALT)
#	cat ${WD}/${READS1}.${READS2}.frc.snp.sam.vcf | awk -v c="$C" '$8>=c && $9>=c' > ${WD}/${READS1}.${READS2}.frc.${C}.snp.sam.vcf
#fi

# 12. Generates VCF (bcftools calls) using BWA MEM + bcftools -> reads1.reads2.bcftools.vcf

if [ "$samtool_version" -eq "1" ]; then
	if [ ! -f ${WD}/${READS1}.${READS2}.bcftools.vcf ]; then
		echo "Calling SNPs using BWA + bcftools. Storing SNPs to file "${WD}/${READS1}.${READS2}.bcftools.vcf" ..."
		/usr/bin/time -v bwa mem ${WD}/${READS1}.reference.fasta ${WD}/${READS2}.fastq -o ${WD}/alignment.tmp.sam 2> ${TIME_BWAMEM}
		/usr/bin/time -v samtools view -b -S ${WD}/alignment.tmp.sam > ${WD}/alignment.tmp.bam 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v samtools sort ${WD}/alignment.tmp.bam > ${WD}/alignment.tmp.sorted.bam 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v samtools index ${WD}/alignment.tmp.sorted.bam 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v bcftools mpileup -f ${WD}/${READS1}.reference.fasta ${WD}/alignment.tmp.sorted.bam | bcftools call -mv -o ${WD}/${READS1}.${READS2}.bcftools_unfilt.vcf 2>> ${TIME_BCFTOOLS}
		
		#filter VCF
		cat ${WD}/${READS1}.${READS2}.bcftools_unfilt.vcf | grep '^#' > ${WD}/H
		cat ${WD}/${READS1}.${READS2}.bcftools_unfilt.vcf | grep -v ^# |  awk '$6>14' > ${WD}/B
		cat ${WD}/H ${WD}/B > ${WD}/${READS1}.${READS2}.bcftools.vcf
		rm ${WD}/H ${WD}/B

		rm ${WD}/*.tmp*
	fi
else
	if [ ! -f ${WD}/${READS1}.${READS2}.bcftools.vcf ]; then
		echo "Calling SNPs using BWA + bcftools. Storing SNPs to file "${WD}/${READS1}.${READS2}.bcftools.vcf" ..."
		/usr/bin/time -v bwa mem ${WD}/${READS1}.reference.fasta ${WD}/${READS2}.fastq -o ${WD}/alignment.tmp.sam 2> ${TIME_BWAMEM}
		/usr/bin/time -v samtools view -b -S ${WD}/alignment.tmp.sam > ${WD}/alignment.tmp.bam 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v samtools sort ${WD}/alignment.tmp.bam ${WD}/alignment.tmp.sorted 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v samtools index ${WD}/alignment.tmp.sorted.bam 2>> ${TIME_BCFTOOLS}

		#OLD SAMTOOLS/BCFTOOLS
		/usr/bin/time -v samtools mpileup -uD -f ${WD}/${READS1}.reference.fasta ${WD}/alignment.tmp.sorted.bam > ${WD}/mpileup.tmp 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v bcftools view -vc ${WD}/mpileup.tmp > ${WD}/${READS1}.${READS2}.bcftools_unfilt.vcf

		#filter VCF
		cat ${WD}/${READS1}.${READS2}.bcftools_unfilt.vcf | grep '^#' > ${WD}/H
		cat ${WD}/${READS1}.${READS2}.bcftools_unfilt.vcf | grep -v ^# |  awk '$6>14' > ${WD}/B
		cat ${WD}/H ${WD}/B > ${WD}/${READS1}.${READS2}.bcftools.vcf
		rm ${WD}/H ${WD}/B

		rm ${WD}/*.tmp*
	fi
fi


# 13. Compute precision/sensitivity of ebwt2snp pipeline using different thresholds for minimum coverage

if [ ! -f ${WD}/${READS1}.${READS2}.report_3 ]; then

	for i in $(seq 3 10); do
		filter_snp ${WD}/${READS1}.${READS2}.frc.snp $i > ${WD}/${READS1}.${READS2}.frc.cov_${i}.snp
		snp_vs_vcf -v ${WD}/${READS1}.${READS2}.bcftools.vcf -c ${WD}/${READS1}.${READS2}.frc.cov_${i}.snp -f ${WD}/${READS1}.reference.fasta > ${WD}/${READS1}.${READS2}.report_${i}
	done
fi



