# ebwt2snp - version 2

### Overview

New version of ebwt2snp. This version works directly on the input BWTs and does not require building the LCP and GSA arrays. This is faster and much more space-efficient. 

The *ebwt2snp* suite can be used to discover SNPs/indels between two sets of reads (fasta/fastq) *without* aligning them to  a reference genome (alignment-free, reference-free) by just comparing the Burrows-Wheeler Transforms of the two datasets. The output is a fasta file (in KisSNP2 format) where sequences are the contexts surrounding the identified SNPs/indels.  
The suite finds its main use in applications where no reference genome is known (alignment-free, reference-free variation discovery). The following modules for SNP/indel discovery are available:

- **ebwt2snp** Takes as input the BWTs of two read collections and detects SNPs/indels between the two datasets. Output: a ".snp" file (this is actually a fasta file in KisSNP++ format containing pairs of reads testifying the variations).
- **filter_snp** Filters the .snp file generated by ebwt2snp so to keep only SNP/indels covered by at least m reads. Increasing m greatly improves precision (but decreases sensitivity), and it is suggested when the input dataset is high-covered (e.g. a value m = 5 is suggested for coverages >= 25x).

To build the BWT of the reads sets, you can use https://github.com/felipelouza/egsa, https://github.com/giovannarosone/BCR_LCP_GSA, or https://github.com/felipelouza/egap. Note also that the **ebwt2snp** pipeline finds many SNPs/indels twice: one time on the forward strand and one on the reverse complement strand.

If a ground-truth VCF file (of first against second individual) and the reference of the first individual are available, one can validate the .snp file generated by **ebwt2snp** using the executable **snp_vs_vcf** (this works on any file in KisSNP++ format). **For now, validation is available for SNPs only (no indels).**

If a reference (of the first individual) is available, one can extend the above pipeline to produce a vcf file. For this, one can use the tools

- **snp2fastq** converts the ".snp" file produced by the **ebwt2snp** pipeline into a ".fastq" file (with fake base qualities) ready to be aligned (e.g. using BWA-MEM) against the reference of the fist individual.
- **sam2vcf** converts the ".sam" file produced by aligning the above ".fastq" into a ".vcf" file containing the variations. 

We call **snp2vcf** the pipeline **snp2fastq -> bwa-mem -> sam2vcf**. Note that bwa-mem requires the reference of the first individual to be available. This reference can be computed, for example, using a standard bwa-mem+{sam,bcf,vcf}tools pipeline. 

To conclude, one can validate the VCF against a ground-truth VCF generated using a standard pipeline (for example, bwa-mem + {sam,bcf,vcf}tools) by using the tool **vcf_vs_vcf**. This is equivalent to validating the .snp file against a ground-truth VCF using **snp_vs_vcf**. 

The paper describing the theory behind the tool (eBWT positional clustering) has been published in:

---

*Nicola Prezza, Nadia Pisanti, Marinella Sciortino and Giovanna Rosone: Detecting mutations by eBWT. WABI 2018. Leibniz International Proceedings in Informatics, LIPIcs , 2018, 113, art. no. 3, Schloss Dagstuhl--Leibniz-Zentrum für Informatik.*

---

A pre-print version can be found here: https://arxiv.org/abs/1805.01876. 


### Install

~~~~
#download ebwt2snp, eGap, and BCR
git clone https://github.com/nicolaprezza/ebwt2snp-v2
git clone https://github.com/felipelouza/egap
git clone https://github.com/giovannarosone/BCR\_LCP\_GSA

#build ebwt2snp
cd ebwt2snp-v2
mkdir build
cd build
cmake ..
make

#build egsa
cd ../../egap
make
~~~~

### Run

Enter the folder with the two fasta files _reads1.fasta_  and _reads2.fasta_ (i.e. the reads of the two samples). We assume that executables 'eGap', 'ebwt2snp' are global. 

~~~~
#Step 1: optional, but considerably increases sensitivity of the tool. Insert in reads1.fasta also the reverse-complement of the reads. Repeat with reads2.fasta

#Build the BWT of the sets of reads
eGap reads1.fasta
eGap reads2.fasta

#Call SNPs (do this in the same folder containing all other files)
ebwt2snp -1 reads1.fasta.bwt -2 reads1.fasta.bwt -o output.snp -t 0

#File ALL.snp.fasta now contains identified SNPs/indels. Note: the third field between "|" in the read-names of this file indicates the number of times the variant is observed (maximum value specified with option -m in ebwt2snp). You can further filter this file according to this field in order to improve accuracy. For this, use executable filter_snp:

#Filter only events supported by at least 5
filter_snp output.snp 5 > output.5.snp

~~~~
