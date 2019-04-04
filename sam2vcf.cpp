// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * First: obtain the KissReads2 file (.snp) with ebwt2snp. Convert it to fastq with seqtk:
 *
 * seqtk seq -F 'h' input.snp > output.fastq
 *
 * Then, align it with BWA-MEM against the reference genome. Finally, this piece of code converts all differences seen in the sam file
 * into variants in a VCF file.
 *
 *
 */

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <set>
#include <cstring>
#include "include.hpp"
#include <algorithm>
#include <map>

using namespace std;

int indel_deduplicate_def = 5;//if 2 indels are within this number of bases, keep only one of the two.
int indel_deduplicate = 0;

bool non_isolated = true;

bool only_exact = false;

int tot_events=0;
int used_events=0;
int good_events=0;

void help(){

	cout << "sam2vcf [OPTIONS]" << endl << endl <<
	"Turns all mismatches/indels seen in a sam file into VCF entries." << endl <<
	"Options:" << endl <<
		"-h          Print this help." << endl <<
		"-f <arg>    Reference fasta file. REQUIRED." << endl <<
		"-s <arg>    Input SAM file. REQUIRED" << endl <<
		"-v <arg>    Output vcf file. REQUIRED." << endl <<
		"-m <arg>    Maximum number of differences, i.e. mismatches + indels (default:5)" << endl;
	exit(0);
}

void parse_cigar(string cigar, int& S, int& M1, int& DI, int& M2, char& ty){

	S=0;
	M1=0;
	DI=0;
	M2=0;

	vector<char> types;
	vector<int> len;

	int L=0;

	for(char c:cigar){

		if('0' <= c and c <= '9'){

			L = L*10 + (int(c)-int('0'));

		}else{

			types.push_back(c);
			len.push_back(L);
			L=0;

		}

	}

	if(types.size()==0) return;

	int k = 0;

	if(types[k]=='S'){

		S = len[k];
		k++;

	}

	if(types[k]=='M'){

		M1 = len[k];
		k++;
		ty = 'M';

	}else return;

	if(types.size()>k){

		if(types[k]=='D' or types[k]=='I'){

			DI = len[k];
			ty = types[k];
			k++;

		}else return;

		if(types.size()>k and types[k]=='M'){

			M2 = len[k];

		}

	}

}

int main(int argc, char** argv){

	if(argc < 2) help();

	string input_fasta;
	string input_sam;
	string output;
	int max_mism = 5;

	int opt;
	while ((opt = getopt(argc, argv, "hs:f:v:m:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'x':
				non_isolated = false;
			break;
			case 'e':
				only_exact = true;
			break;
			case 'd':
				indel_deduplicate = atoi(optarg);
			break;
			case 'f':
				input_fasta = string(optarg);
			break;
			case 's':
				input_sam = string(optarg);
			break;
			case 'v':
				output = string(optarg);
			break;
			case 'm':
				max_mism = atoi(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	if(input_fasta.compare("")==0) help();
	if(input_sam.compare("")==0) help();
	if(output.compare("")==0) help();

	ifstream ifasta(input_fasta);
	ifstream isam(input_sam);
	ofstream out(output);

	//reference indexed by contig: e.g. string chr1 = ref["chr1"]
	std::map<string,string> ref;

	cout << "Loading reference ... " << flush;

	/*
	 *
	 * READ REFERENCE, LOAD CONTIGS IN MEMORY
	 *
	 */

	string line;
	string contig;

	vector<string> contigs;

	int nonisolated_snps=0;

	while(not ifasta.eof()){

		getline(ifasta, line);

		if(line[0] == '>'){//new contig name

			contig = line.substr(1);
			contigs.push_back(contig);

			ref[contig] = string();

		}else{

			for (auto & c: line) c = toupper(c);
			ref[contig].append(line);

		}

	}

	cout << "done." << endl;

	cout << "Contig\tlength" << endl;
	for(auto c : contigs){

		cout << c << "\t" << ref[c].length() << endl;

	}

	ifasta.close();

	string str;
	out << "#CHROM\tPOS\tID\tREF\tALT\tTYPE" << endl;

	int ID = 1; //event ID

	while(getline(isam, str)){

		if(str[0]!='@' and str[0]!='['){//skip header

			tot_events++;

			string name;//entry name
			string flag;
			string chr;
			string pos;
			string mapq;
			string cigar;
			string rnext;
			string pnext;
			string tlen;
			string seq;
			string qual;
			string NM;

			std::istringstream iss(str);
			getline(iss, name, '\t');
			getline(iss, flag, '\t');
			getline(iss, chr, '\t');
			getline(iss, pos, '\t');
			getline(iss, mapq, '\t');
			getline(iss, cigar, '\t');
			getline(iss, rnext, '\t');
			getline(iss, pnext, '\t');
			getline(iss, tlen, '\t');
			getline(iss, seq, '\t');
			getline(iss, qual, '\t');
			getline(iss, NM, '\t');

			uint64_t pos_int = atoi(pos.c_str());

			int n_mism = 0;

			//get number of mismatches (includes indels and mismatches)
			{

				std::istringstream stream(NM);
				string nm_s;
				getline(stream, nm_s, ':');
				getline(stream, nm_s, ':');
				getline(stream, nm_s, ':');

				n_mism = atoi(nm_s.c_str());

			}

			//if mismatch/indel is present
			if(n_mism>0 and n_mism <= max_mism and ref[chr].size()>0){

				good_events++;

				//only admitted cigar formats: wS xM or xMy{I,D}zM, where w,x,y,z are integers

				int S=0,M1=0,DI=0,M2=0;
				char type = 0;

				parse_cigar(cigar,S,M1,DI,M2,type);

				//starting positions in read and genome (skip S chars)
				int k = S;

				if(M1>0){

					used_events++;

					//scan the first x bases in REF and ALT
					for(int j=0;j<M1;++j){

						if(ref[chr][pos_int + j - 1] != seq[k]){

							out << chr << "\t" << (pos_int + j) << "\t" << (ID++) << "\t" << ref[chr][pos_int + j - 1] << "\t" << seq[k] << "\tSNP" << endl;

						}

						k++;

					}

					//if there's an indel
					if(DI>0){

						//insert of length y in reference
						if(type == 'I'){

							string REF = ref[chr].substr(pos_int + M1 - 2,1);

							string ALT = seq.substr(k-1,DI+1);

							k += DI;

							out << chr << "\t" << (pos_int + (M1-1)) << "\t" << (ID++) << "\t" << REF << "\t" << ALT << "\tINDEL" << endl;

						}else if (type == 'D'){//deletion in reference

							string REF = ref[chr].substr(pos_int + M1 - 2,DI+1);
							string ALT = seq.substr(k-1,1);
							out << chr << "\t" << (pos_int + (M1-1)) << "\t" << (ID++) << "\t" << REF << "\t" << ALT << "\tINDEL" << endl;

						}

					}

					if(M2>0){//scan the last z bases and look for further SNPs

						for(int j=0;j<M2;++j){

							uint64_t start = pos_int + M1 + (type=='D'?DI:0) -1;//start is 0-based

							if(ref[chr][start+j] != seq[k]){

								out << chr << "\t" << (start + j +1) << "\t" << (ID++) << "\t" << ref[chr][start+j] << "\t" << seq[k] << "\tSNP" << endl;

							}

							k++;

						}

					}

				}

			}

		}

	}

	cout << "Done." << endl;
	cout << "Alignments used to generate VCF entries/alignments with 0 < errors <= max_diff: " << used_events << "/" << good_events << endl;

}
