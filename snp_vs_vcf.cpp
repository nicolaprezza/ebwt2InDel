// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include "include.hpp"
#include <unistd.h>
#include <map>
#include <sstream>
#include <set>

using namespace std;

string vcf_path;
string calls_path;
string ref_path;

int k_nonis_def = 31;
int k_nonis=0;
int rlength_def = 100;
int rlength = 0;

void help(){

	cout << "snp_vs_vcf [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-v <arg>    VCF file with the ground-truth SNPs (REQUIRED)" << endl <<
	"-c <arg>    Calls in KisSNP2 format (REQUIRED)" << endl <<
	"-f <arg>    Reference fasta file of first sample (REQUIRED)" << endl <<
	"-k <arg>    Value to define non-isolated SNPs (default: " << k_nonis_def << ")" << endl <<
	"-l <arg>    Max read length (default: " << rlength_def << ")" << endl;
	exit(0);
}

struct call{

	string right_context;//context following SNP (SNP excluded)
	string left_context;//context preceding SNP (SNP excluded)
	char REF;
	char ALT;
	uint64_t ID;
	bool isolated;
	int pos;

};

char RC(char c){

	switch (c){

		case 'A': return 'T'; break;
		case 'C': return 'G'; break;
		case 'G': return 'C'; break;
		case 'T': return 'A'; break;

	}

	return c;

}

string RC(string read){

	string rev(read.rbegin(), read.rend());

	for(auto & c : rev) c = RC(c);

	return rev;

}

string REV(string read){

	string rev(read.rbegin(), read.rend());

	return rev;

}

//is a prefix of b?
bool is_prefix(string &a, string &b){

	bool res = true;

	if(a.length()>b.length()) return false;

	for(int i=0;i<a.length();++i){

		res = res and a[i]==b[i];

	}

	return res;

}

bool comp(const call &a, const call &b){

	return a.right_context.compare(b.right_context) < 0;

}

int main(int argc, char** argv){

	if(argc < 4) help();

	int opt;
	while ((opt = getopt(argc, argv, "hv:c:f:l:k:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'v':
				vcf_path = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			case 'c':
				calls_path = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			case 'f':
				ref_path = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			case 'l':
				rlength = atoi(optarg);
				//cout << "input = " << input << "\n";
			break;
			case 'k':
				k_nonis = atoi(optarg);
				//cout << "input = " << input << "\n";
			break;
			default:
				help();
			return -1;
		}
	}

	rlength = rlength == 0 ? rlength_def : rlength;
	k_nonis = k_nonis == 0 ? k_nonis_def : k_nonis;

	if(vcf_path.compare("")==0 or calls_path.compare("")==0 or ref_path.compare("")==0) help();

	//reference indexed by contig: e.g. string chr1 = ref["chr1"]
	std::map<string,string> ref;

	cout << "Loading reference ... " << flush;

	/*
	 *
	 * READ REFERENCE, LOAD CONTIGS IN MEMORY
	 *
	 */
	ifstream ref_file;
	ref_file.open(ref_path, ios::in);

	string line;
	string contig;

	vector<string> contigs;

	int nonisolated_snps=0;

	while(not ref_file.eof()){

		getline(ref_file, line);

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

	uint64_t N = 0;

	cout << "Contig\tlength" << endl;
	for(auto c : contigs){

		cout << c << "\t" << ref[c].length() << endl;

		N+=ref[c].length();

	}
	//cout << endl;

	ref_file.close();


	/*
	 * LOAD SNPS FROM VCF
	 */

	cout << "Loading VCF ... " << flush;

	uint64_t n_snps=0;//number of SNPs in the VCF

	//read vcf
	ifstream vcf_file;
	vcf_file.open(vcf_path, ios::in);

	vector<call> calls_vcf;
	uint64_t ID = 0;

	while(not vcf_file.eof()){

		getline(vcf_file, line);

		if(line.size()>0 and line[0] != '#'){

			std::istringstream is( line );

			string chr;
			int pos;
			string id;
			string REF;
			string ALT;

			is >> chr >> pos >> id >> REF >> ALT;

			pos--;//coordinates are 1-based in the vcf file

			//keep just SNPs
			if(		(REF.compare("A")==0 or REF.compare("C")==0 or REF.compare("T")==0 or REF.compare("G")==0) and
					(ALT.compare("A")==0 or ALT.compare("C")==0 or ALT.compare("T")==0 or ALT.compare("G")==0)
			){

				if(ref[chr].compare("")!=0){//if chromosome exists in the reference file

					n_snps++;

					//insert forward call

					if(pos >= ref[chr].size()){

						cout << "WARNING: position " << pos << " larger than chromosome " << chr << "'s length " << ref[chr].size() << endl;

					}

					if(pos >= rlength && pos+rlength < ref[chr].size()){

						assert(pos+1<ref[chr].size());
						assert(pos+1+rlength<=ref[chr].size());
						assert(pos>=rlength);
						assert(pos-rlength<ref[chr].size());
						assert(pos<=ref[chr].size());

						string right_context = ref[chr].substr(pos+1,rlength);
						string left_context = REV(ref[chr].substr(pos-rlength,rlength));

						assert(REF.size()>0);
						assert(ALT.size()>0);

						calls_vcf.push_back(call {right_context, left_context, REF[0], ALT[0], ID, true, pos});

						//insert RC call

						left_context = REV(RC(ref[chr].substr(pos+1,rlength)));
						right_context = RC(ref[chr].substr(pos-rlength,rlength));

						calls_vcf.push_back(call {right_context, left_context, RC(REF[0]), RC(ALT[0]), ID, true, pos});

					}

					++ID;

				}else{

					cout << "WARNING: chromosome " << chr << " not found. " << endl;

				}

			}

		}

	}


	if(calls_vcf.size()==0){

		cout << "WARNING: no variants found. Check that chromosome names are the same in the fasta and vcf files. " << endl;

	}

	for(int i=2;calls_vcf.size() > 1 && i<calls_vcf.size()-2;++i){

		if(i%2==0){

			calls_vcf[i].isolated = calls_vcf[i].pos - calls_vcf[i-2].pos >= k_nonis and calls_vcf[i+2].pos - calls_vcf[i].pos >= k_nonis;
			calls_vcf[i+1].isolated = calls_vcf[i].isolated;

			if(not calls_vcf[i].isolated){

				//cout << calls_vcf[i].pos - calls_vcf[i-2].pos << " " << calls_vcf[i+2].pos - calls_vcf[i].pos << endl;
				nonisolated_snps++;

			}

		}

	}

	cout << "done." << endl;

	cout << "Sorting VCF by context ... " << flush;
	std::sort(calls_vcf.begin(), calls_vcf.end(), comp);
	cout << "done." << endl;

	vcf_file.close();

	//vector<call> calls;

	/*
	 *
	 * LOAD CALLS FROM FASTA GENERATED BY SNP-CALLING TOOL, SEARCH THEM IN THE RANGE STRUCTURE OF
	 * CONTAINING VCF CALLS
	 *
	 */

	uint64_t n_calls = 0;//number of SNPs in input fasta file (including hidden SNPs inside the sequences)

	cout << "Checking calls ... " << flush;

	uint64_t FP = 0;
	uint64_t FN = 0;
	uint64_t TP = 0;
	uint64_t TN = 0;

	auto assigned = vector<int>(calls_vcf.size(),0);

	//read calls
	ifstream calls_file;
	calls_file.open(calls_path, ios::in);

	//left/right contexts per sample
	vector<string> left_contexts = vector<string>(2);
	string right_context;

	getline(calls_file, line);

	int last_id = 0;
	int non_isolated_SNPs = 0;

	while(not calls_file.eof()){

		//header of first read

		string h1=line;

		std::istringstream is( line );
		string pos;
		getline(is,pos,'|');

		if(pos.substr(0,4).compare(">SNP")==0){//filter out indels

			getline(is,pos,'|');

			is = istringstream(pos);
			getline(is,pos,':');
			getline(is,pos,':');

			is = istringstream(pos);
			getline(is,pos,'_');

			int ipos = atoi(pos.c_str());//position of leftmost SNP
			ipos = 0;

			string DNA1, DNA2;
			getline(calls_file, DNA1);//DNA of first read

			getline(calls_file, DNA2);//skip header of second read
			string h2=DNA2;

			getline(calls_file, DNA2);//DNA of second read

			//some consistency checks on the file

			/*if(DNA1.substr(DNA1.size() - ipos).compare(DNA2.substr(DNA2.size() - ipos)) != 0){

				cout << "Error: malformed SNP file. Last " << ipos << " bases do not coincide\n";
				cout << h1 << endl << DNA1 << endl << h2 << endl << DNA2 << endl;
				exit(1);

			}*/

			if(DNA1.length()!=DNA2.length()){

				cout << "Error: malformed SNP file. Two reads with different length in a SNP:\n";
				cout << h1 << endl << DNA1 << endl << h2 << endl << DNA2 << endl;
				exit(1);

			}

			/*if(DNA1[DNA1.size()-ipos-1] == DNA2[DNA2.size()-ipos-1]){

				cout << "Error: malformed SNP file. No SNP at specified position:\n";
				cout << h1 << endl << DNA1 << endl << h2 << endl << DNA2 << endl;
				exit(1);

			}*/

			//search all SNPs from the leftmost backwards
			while(ipos<DNA1.size()){

				//if SNP found
				if(DNA1[DNA1.size()-ipos-1] != DNA2[DNA2.size()-ipos-1]){

					n_calls++;

					string right_context1 = DNA1.substr(DNA1.size()-ipos);
					string right_context2 = DNA2.substr(DNA2.size()-ipos);

					string left_context1 = REV(DNA1.substr(0, DNA1.size()-ipos-1));
					string left_context2 = REV(DNA2.substr(0, DNA2.size()-ipos-1));

					char REF = DNA1[DNA1.size()-ipos-1];
					char ALT = DNA2[DNA2.size()-ipos-1];

					call c1 = {right_context1, left_context1, REF, ALT, 0, false};
					call c2 = {right_context2, left_context2, REF, ALT, 0, false};

					non_isolated_SNPs++;

					//Search c1 in the range structure

					auto it = lower_bound(calls_vcf.begin(), calls_vcf.end(),c1, comp);

					uint64_t idx = std::distance(calls_vcf.begin(), it);

					bool found = false;

					/*
					 * while right context of c1 prefixes right context in range of VCF calls
					 */

					while(idx < calls_vcf.size() and is_prefix(right_context1, calls_vcf[idx].right_context)){

						if( ((calls_vcf[idx].ALT == ALT and calls_vcf[idx].REF == REF) or (calls_vcf[idx].ALT == REF and calls_vcf[idx].REF == ALT)) and
							(is_prefix(left_context1, calls_vcf[idx].left_context))	){

							found = true;
							assigned[idx] = 1;

						}

						++idx;

					}

					if(not found){

						//Search c2 in the range structure

						it = lower_bound(calls_vcf.begin(), calls_vcf.end(),c2, comp);

						idx = std::distance(calls_vcf.begin(), it);

						/*
						 * while right context of c2 prefixes right context in range of VCF calls
						 */
						while(idx < calls_vcf.size() and is_prefix(right_context2, calls_vcf[idx].right_context)){

							if( ((calls_vcf[idx].ALT == ALT and calls_vcf[idx].REF == REF) or (calls_vcf[idx].ALT == REF and calls_vcf[idx].REF == ALT)) and
								(is_prefix(left_context2, calls_vcf[idx].left_context))	){

								found = true;
								assigned[idx] = 1;

							}

							++idx;

						}

					}

					if(not found) FP++;

				}

				++ipos;

			}

			non_isolated_SNPs--;

		}else{

			//skip event

			getline(calls_file, line);//DNA of first read
			getline(calls_file, line);//header of second read
			getline(calls_file, line);//DNA of second read

		}

		getline(calls_file, line);//header of first read

	}

	cout << "done." << endl;

	calls_file.close();



	vector<call> found_calls_vcf;//count how many SNPs in the VCF have been found

	std::set<int> found;

	for (uint64_t i = 0; i< calls_vcf.size();++i){

		if(assigned[i]){

			found.insert(calls_vcf[i].ID);

		}

	}

	//nonisolated SNPs found
	std::set<int> found_nonisolated;

	for (uint64_t i = 0; i< calls_vcf.size();++i){

		if(assigned[i] and not calls_vcf[i].isolated){

			found_nonisolated.insert(calls_vcf[i].ID);

		}

	}

	TP = found.size();

	//number of true SNPs that are not found
	FN = n_snps - TP;

	//this is a lower bound because we may call multiple times the same SNP
	TN = (N - n_calls) -  FN;

	cout << endl << "Non-isolated SNPs detected: " << found_nonisolated.size()  << "/" << nonisolated_snps << endl;

	cout << endl;
	cout << "TP = " << TP << endl;
	cout << "TN = " << TN << endl;
	cout << "FP = " << FP << endl;
	cout << "FN = " << FN << endl;

	cout << "sensitivity = TP/(TP+FN) = " << 100*double(TP)/(double(TP)+double(FN)) << "%" << endl;
	cout << "specificity = TN/(TN+FP) = " << 100*double(TN)/(double(TN)+double(FP)) << "%" << endl;
	cout << "precision   = TP/(TP+FP) = " << 100*double(TP)/(double(TP)+double(FP)) << "%" << endl;



}
