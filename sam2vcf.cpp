// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

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

using namespace std;

int indel_deduplicate_def = 5;//if 2 indels are within this number of bases, keep only one of the two.
int indel_deduplicate = 0;

bool non_isolated = true;

bool only_exact = false;

void help(){

	cout << "sam2vcf [OPTIONS]" << endl << endl <<
	"Converts the aligned calls (with bwa-mem) 'calls.sam' of clust2snp into a vcf file 'calls.sam.vcf'." << endl <<
	"Options:" << endl <<
		"-h          Print this help." << endl <<
		"-x          Disable non-isolated SNPs (default: enabled)." << endl <<
		"-s <arg>    Input SAM file. REQUIRED" << endl <<
		"-d <arg>    Keep only one indel in pairs within <arg> bases. Default: " <<  indel_deduplicate_def << "." << endl <<
		"-e          Keep only exact alignments." << endl;
	exit(0);
}

struct vcf_entry{

	string chr;
	uint64_t pos;

	string REF;
	string ALT;

	bool indel;

	bool exact;

	uint64_t cov_ref;
	uint64_t cov_alt;

	void print(){

		cout << endl;
		cout << chr << " " << pos << " " << REF << " " << ALT << " " << indel << " " << exact << " " << cov_ref << " " << cov_alt << endl;

	}

	bool operator<(const vcf_entry & a) const{

		if(chr.compare(a.chr) < 0 ) return true;
		else if(chr.compare(a.chr) == 0 ){

			if(pos == a.pos){

				string snp1 = REF;
				snp1.append(ALT);

				string snp2 = a.REF;
				snp2.append(a.ALT);

				return snp1.compare(snp2) < 0;

			}

			return pos < a.pos;

		}

		return false;

	}

	bool operator==(const vcf_entry & a) const{

		return chr.compare(a.chr)==0 and pos == a.pos and REF.compare(a.REF) == 0 and ALT.compare(a.ALT) == 0;

	}

};

int main(int argc, char** argv){

	if(argc < 2) help();

	string infile;

	int opt;
	while ((opt = getopt(argc, argv, "xhd:s:e")) != -1){
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
			case 's':
				infile = string(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	indel_deduplicate = indel_deduplicate==0 ? indel_deduplicate_def : indel_deduplicate;

	if(infile.compare("")==0) help();

	string outfile = infile;
	outfile.append(".vcf");

	ifstream is(infile);
	ofstream of(outfile);

	string str;

	of << "#CHROM\tPOS\tID\tREF\tALT\tTYPE\tEXACT\tCOV_REF\tCOV_ALT" << endl;
	//cout << "#CHROM\tPOS\tID\tREF\tALT\tINFO" << endl;

	vector<vcf_entry> VCF;

	while(getline(is, str)){

		if(str[0]!='@' and str[0]!='['){//skip header

			string name;//entry name
			std::istringstream iss(str);
			getline(iss, name, '\t');

			string event_nr;//INDEL or SNP
			string type;//INDEL or SNP
			string REF;//reference allele
			string ALT;//alternative allele
			string REF_dna;//reference dna
			string ALT_dna;//alternative dna
			uint64_t pos;//alignment position
			uint64_t COV_REF;//reads supporting variation on individual 1
			uint64_t COV_ALT;//reads supporting variation on individual 2
			int snp_pos;
			string chr;
			string cigar;
			string mismatches;

			string token;

			//First individual = reference = read DNA
			//Second individual = ALT = DNA in header

			{

			std::istringstream iss1(name);

			getline(iss1, type, '_');
			getline(iss1, event_nr, '_');
			getline(iss1, token, '_');
			snp_pos = atoi(token.c_str());

			string REFALT;
			getline(iss1, REFALT, '_');

			std::istringstream iss2(REFALT);
			getline(iss2, REF, '/');
			ALT=string();
			getline(iss2, ALT, '/');

			getline(iss1, token, '_');
			COV_REF = atoi(token.c_str());
			getline(iss1, token, '_');
			COV_ALT = atoi(token.c_str());

			getline(iss1, ALT_dna, '_');

			}

			getline(iss, token, '\t');//flag
			unsigned int f = atoi(token.c_str());
			getline(iss, chr, '\t');
			getline(iss, token, '\t');
			pos = atoi(token.c_str());
			getline(iss, token, '\t');
			getline(iss, cigar, '\t');

			getline(iss, token, '\t');
			getline(iss, token, '\t');
			getline(iss, token, '\t');
			getline(iss, REF_dna, '\t');

			getline(iss, token, '\t');//quality

			getline(iss, mismatches, '\t');//NM:i:x

			//exact alignment: 0 mismatches and 0 skips/indels
			bool exact = 	std::count(cigar.begin(), cigar.end(), 'S') == 0 and
							std::count(cigar.begin(), cigar.end(), 'I') == 0 and
							std::count(cigar.begin(), cigar.end(), 'D') == 0 and
							mismatches.compare("NM:i:0")==0;

			bool reversed = (f & (unsigned int)16) != 0;
			bool indel = type.compare("INDEL")==0;

			if(reversed){//then REF_DNA has been reverse-complemented by the aligner. apply reverse-complement also ALT_DNA and the variant

				ALT_dna = RC(ALT_dna);
				REF = RC(REF);
				ALT = RC(ALT);

				if(indel){

					snp_pos--;

				}

			}


			/*cout << "event_nr " << event_nr << endl;
			cout << "type " << type<< endl;//INDEL or SNP
			cout <<  "REF " <<  REF<< endl;//reference allele
			cout <<  "ALT " <<  ALT<< endl;//alternative allele
			cout <<  "REF_dna " <<  REF_dna<< endl;//reference dna
			cout <<  "ALT_dna " <<  ALT_dna<< endl;//alternative dna
			cout <<  "pos " <<  pos<< endl;//alignment position
			cout <<  "COV_REF " <<  COV_REF<< endl;//reads supporting variation on individual 1
			cout <<  "COV_ALT " <<  COV_ALT<< endl;//reads supporting variation on individual 2
			cout <<  "snp_pos " <<  snp_pos<< endl;
			cout <<  "chr " <<  chr<< endl;
			cout <<  "cigar " <<  cigar<< endl;
			cout <<  "mismatches " <<  mismatches<< endl;
			cout <<  "flag " <<  f<< endl;
			cout <<  "indel? " <<  indel<< endl;
			cout <<  "reversed? " <<  reversed<< endl<< endl;*/


			//adjust snp_pos in the case we are on FW strand
			if(not reversed){

				if(indel){

					if(REF.length()>0){

						//insert in REF
						int indel_len = REF.length();
						snp_pos = ((REF_dna.length() - snp_pos) - indel_len) -1;

					}else{

						//insert in ALT
						snp_pos = (REF_dna.length() - snp_pos) - 1;

					}

				}else{

					snp_pos = (REF_dna.length() - snp_pos)-1;

				}

			}

			if((f==0 or f==16) and snp_pos >= 0 and chr.compare("*") != 0){

				vcf_entry v;

				if(indel){

					if(REF.length()>0){

						v = {
										chr,
										pos + snp_pos,
										REF_dna.substr(snp_pos,REF.length()+1),
										ALT_dna.substr(snp_pos,1),
										true,
										exact,
										COV_REF,
										COV_ALT
						};

					}else{

						v = {
										chr,
										pos + snp_pos,
										REF_dna.substr(snp_pos,1),
										ALT_dna.substr(snp_pos,ALT.length()+1),
										true,
										exact,
										COV_REF,
										COV_ALT
						};

					}


				}else{

					//cout << "SNP:\n" << REF_dna << endl << ALT_dna << endl << REF_dna.substr(snp_pos,1) << " " << ALT_dna.substr(snp_pos,1)  << endl;

					v = {
									chr,
									pos + snp_pos,
									REF,
									ALT,
									false,
									exact,
									COV_REF,
									COV_ALT
					};

				}

				if((not only_exact) or exact){

					VCF.push_back(v);

					/*
					 * find non-isolated SNPs
					 */

					if(non_isolated){

						if(indel){

							if(reversed){

								//non-isolated SNPs are on the right of the end position of indel

								int indel_length = REF.length() > 0 ? REF.length() : ALT.length();

								int L = std::max(REF_dna.length(), ALT_dna.length());//length of fragment containing the insert

								int len_right = L - (snp_pos+indel_length) -1; //length of right part in common (with potential SNPs)

								int snp_pos_ref = REF.length() > 0 ? snp_pos + indel_length +1 : snp_pos +1;
								int snp_pos_alt = REF.length() > 0 ? snp_pos +1 : snp_pos + indel_length +1;

								for(int i=0;i<len_right;++i){

									if(REF_dna[snp_pos_ref+i] != ALT_dna[snp_pos_alt+i]){

										vcf_entry v  = {
															chr,
															pos + i,
															REF_dna.substr(snp_pos_ref+i,1),
															ALT_dna.substr(snp_pos_alt+i,1),
															false,
															exact,
															COV_REF,
															COV_ALT
											};

										VCF.push_back(v);

									}

								}

							}else{

								//non-isolated SNPs are on the left of the start position of indel

								for(int i=0;i<snp_pos;++i){

									if(REF_dna[i] != ALT_dna[i]){

										vcf_entry v  = {
															chr,
															pos + i,
															REF_dna.substr(i,1),
															ALT_dna.substr(i,1),
															false,
															exact,
															COV_REF,
															COV_ALT
											};

										VCF.push_back(v);

									}

								}

							}

						}else{

							if(reversed){

								//non-isolated SNPs are on the right

								for(int i=snp_pos+1;i<REF_dna.length();++i){

									if(REF_dna[i] != ALT_dna[i]){

										//cout << "SNP (non IS):\n" << REF_dna << endl << ALT_dna << endl << REF_dna.substr(i,1) << " " << ALT_dna.substr(i,1)  << endl;

										vcf_entry v  = {
															chr,
															pos + i,
															REF,
															ALT,
															false,
															exact,
															COV_REF,
															COV_ALT
											};

										VCF.push_back(v);

									}

								}

							}else{

								//non-isolated SNPs are on the left

								for(int i=0;i<snp_pos;++i){

									if(REF_dna[i] != ALT_dna[i]){

										//cout << "SNP (non IS):\n" << REF_dna << endl << ALT_dna << endl << REF_dna.substr(i,1) << " " << ALT_dna.substr(i,1)  << endl;

										vcf_entry v  = {
															chr,
															pos + i,
															REF,
															ALT,
															false,
															exact,
															COV_REF,
															COV_ALT
											};

										VCF.push_back(v);

									}

								}

							}

						}

					}

				}

			}

		}

	}

	//remove duplicates

	uint64_t old_size = VCF.size();

	std::sort(VCF.begin(),VCF.end());

	/*for(int i=0;i<VCF.size()-1;++i){

		if(VCF[i]==VCF[i+1]){

			cout << VCF[i].chr << "\t" << VCF[i].pos << "\t" << ".\t" << VCF[i].REF << "\t" << VCF[i].ALT << "\t" << (VCF[i].indel?"INDEL":"SNP") << endl;
			cout << VCF[i+1].chr << "\t" << VCF[i+1].pos << "\t" << ".\t" << VCF[i+1].REF << "\t" << VCF[i+1].ALT << "\t" << (VCF[i+1].indel?"INDEL":"SNP") << endl;

			cout << endl;

		}

	}*/

	VCF.erase( unique( VCF.begin(), VCF.end() ), VCF.end() );

	//remove close indels

	vector<vcf_entry> VCF_filt;
	for(int i=0;i<VCF.size()-1;++i){

		if(		VCF[i].indel and
				VCF[i+1].indel and
				VCF[i].chr.compare(VCF[i+1].chr)==0 and
				std::abs( int(VCF[i].pos)-int(VCF[i+1].pos) ) <= indel_deduplicate ){

			VCF_filt.push_back(VCF[i]);//keep first indel
			i++;//skip second indel

		}else{

			VCF_filt.push_back(VCF[i]);

		}

	}

	VCF_filt.push_back(VCF[VCF.size()-1]);

	uint64_t new_size = VCF_filt.size();

	cout << (old_size-new_size) << " duplicates found. Saving remaining " << new_size << " unique SNPS/indels." << endl;

	uint64_t n_indels = 0;
	uint64_t n_snps = 0;

	for(auto v:VCF_filt){

		of	<< v.chr << "\t"
			<< v.pos << "\t" << ".\t"
			<< v.REF << "\t"
			<< v.ALT << "\t"
			<< (v.indel?"INDEL":"SNP") << "\t"
			<< v.exact << "\t"
			<< v.cov_ref << "\t"
			<< v.cov_alt << endl;

		n_indels += v.indel;
		n_snps += (not v.indel);

	}

	cout << "Number of SNPs found: " << n_snps << endl;
	cout << "Number of indels found: " << n_indels << endl;

	is.close();
	of.close();

}
