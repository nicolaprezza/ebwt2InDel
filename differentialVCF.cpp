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

string vcf_path1;
string vcf_path2;
string ref_path;
string out_folder;


void help(){

	cout << "relativeVCF [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-f <argf>    Reference fasta file used to generate the VCFs. REQUIRED" << endl <<
	"-1 <arg1>    VCF file with SNPs of individual 1 w.r.t. reference. REQUIRED." << endl <<
	"-2 <arg2>    VCF file with SNPs of individual 2 w.r.t. reference. REQUIRED" << endl <<
	"-o <argo>    output directory (ending with slash, e.g. /home/). REQUIRED." << endl << endl <<
	"creates two files in <argo>: a new reference applying the SNPs <arg1> to " << endl <<
	"<argf>, and a new differential VCF file containing the SNPs of individual " << endl <<
	"2 relative to the new reference (therefore to individual 1). Note: only " << endl <<
	"SNPs are admitted in the input VCF files." << endl;
	exit(0);
}

struct call{

	string contig;
	int pos;
	char REF;
	char ALT;

	bool operator<(call other) const{

		if(contig.compare(other.contig) < 0){

			return true;

		}else if(contig.compare(other.contig) == 0){

			return pos < other.pos;

		}

		return false;

	}

};

int main(int argc, char** argv){

	if(argc < 5) help();

	int opt;
	while ((opt = getopt(argc, argv, "h1:2:f:o:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case '1':
				vcf_path1 = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			case '2':
				vcf_path2 = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			case 'f':
				ref_path = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			case 'o':
				out_folder = string(optarg);
				//cout << "input = " << input << "\n";
			break;
			default:
				help();
			return -1;
		}
	}

	if(		vcf_path1.compare("")==0 or
			vcf_path2.compare("")==0 or
			out_folder.compare("")==0 or
			ref_path.compare("")==0 or
			out_folder.compare("")==0 or
			out_folder[out_folder.length()-1] != '/'
	) help();

	string new_ref_file = out_folder;
	new_ref_file.append("new_reference.fa");

	string differential_vcf_file = out_folder;
	differential_vcf_file.append("differential.vcf");


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

	cout << "Loading VCFs ... " << endl;

	uint64_t n_snps=0;//number of SNPs in the VCF

	//read vcf
	ifstream vcf_file1;
	vcf_file1.open(vcf_path1, ios::in);

	ifstream vcf_file2;
	vcf_file2.open(vcf_path2, ios::in);

	set<call> calls_vcf1;
	set<call> calls_vcf2;
	set<call> calls_vcfOut;

	//load VCF1

	int tot_calls = 0;
	int discarded_calls = 0;//keep track of calls that disagree with reference and throw a warning
	int max_warnings = 50; //max number of warnings to show

	while(not vcf_file1.eof()){

		getline(vcf_file1, line);

		if(line[0] != '#' and line.compare("")!=0){

			std::istringstream is( line );

			string chr;
			int pos;
			string id;
			string REF;
			string ALT;

			is >> chr >> pos >> id >> REF >> ALT;

			pos--;//coordinates are 1-based in the vcf file

			if(ref[chr].compare("")!=0 && pos < ref[chr].length() && ref[chr][pos] == REF[0]){

				calls_vcf1.insert(call {chr, pos, REF[0], ALT[0]});

			}else{

				if(discarded_calls<max_warnings){

					cout << "WARNING: call \"" << chr << " " << (pos+1) << " " << REF << " " << ALT << "\" of file " << vcf_path1 << " does not match the reference." << endl;
					cout << "   problem: " << flush;

					if(ref[chr].compare("")==0){

						cout << "contig \"" <<  chr << "\" does not exist." << endl;

					}else if(pos >= ref[chr].length()){

						cout << "VCF position exceeds the contig's length " << ref[chr].length() << endl;

					}else{

						cout << "Reference base " << ref[chr][pos] << " does not match VCF" << endl;

					}

				}

				discarded_calls++;

			}

			tot_calls++;

		}

	}

	//load VCF2

	while(not vcf_file2.eof()){

		getline(vcf_file2, line);

		if(line[0] != '#' and line.compare("")!=0){

			std::istringstream is( line );

			string chr;
			int pos;
			string id;
			string REF;
			string ALT;

			is >> chr >> pos >> id >> REF >> ALT;


			pos--;//coordinates are 1-based in the vcf file

			if(ref[chr].compare("")!=0 && pos < ref[chr].length() && ref[chr][pos] == REF[0]){

				calls_vcf2.insert(call {chr, pos, REF[0], ALT[0]});

			}else{

				if(discarded_calls<max_warnings){

					cout << "WARNING: call \"" << chr << " " << (pos+1) << " " << REF << " " << ALT << "\" of file " << vcf_path2 << " does not match the reference." << endl;
					cout << "   problem: " << flush;

					if(ref[chr].compare("")==0){

						cout << "contig \"" <<  chr << "\" does not exist." << endl;

					}else if(pos >= ref[chr].length()){

						cout << "VCF position exceeds the contig's length " << ref[chr].length() << endl;

					}else{

						cout << "Reference base " << ref[chr][pos] << " does not match VCF" << endl;

					}

				}

				discarded_calls++;

			}

			tot_calls++;

		}

	}

	cout << "done." << endl;

	if(discarded_calls>=max_warnings){

		cout << "(" << (discarded_calls-max_warnings) << " more WARNINGS not shown)" << endl;

	}

	if(discarded_calls>0){

		cout << "WARNING: " << discarded_calls <<  "/" << tot_calls << " VCF entries disagreed with reference and were discarded." << endl;

	}else{

		cout << "All VCF calls agreed with reference." << endl;

	}

	cout << "Computing differential VCF ... " << flush;

	//compute differential VCF

	for(auto ind2 : calls_vcf2){

		assert(ind2.REF == ref[ind2.contig][ind2.pos]);

		auto ind1 = calls_vcf1.find(ind2);

		if(ind1 != calls_vcf1.end()){

			//if call found

			assert(ind2.REF == ind1->REF);//same position, therefore reference base must match

			if(ind2.ALT != ind1->ALT){

				//if the REF is mutated into 2 different bases in the two idividuals

				//same coordinates, but individual 1 becomes REF and individual 2 becomes ALT
				calls_vcfOut.insert(call {ind2.contig, ind2.pos, ind1->ALT, ind2.ALT});

			}//else: same mutation in the two individuals, therefore nothing to report in the differential VCF

			//in both cases we need to mutate the reference

			if(ref[ind2.contig].compare("")!=0){//if chromosome exists in the reference file

				assert(ind2.pos<ref[ind2.contig].size());
				ref[ind2.contig][ind2.pos] = ind1->ALT;

			}

		}else{

			//call not found: individual 1 has the reference base
			//insert the same call in the relative VCF

			calls_vcfOut.insert(ind2);

		}

	}


	for(auto ind1 : calls_vcf1){

		auto ind2 = calls_vcf2.find(ind1);

		//note: take care only of the case where the call is in VCF 1 but not in VCF 2.
		//the other case has already been taken into account in the previous for loop.
		if(ind2 == calls_vcf2.end()){

			assert(ind1.REF == ref[ind1.contig][ind1.pos]);

			if(ref[ind1.contig].compare("")!=0){

				//insert call
				assert(ind1.pos<ref[ind1.contig].size());
				calls_vcfOut.insert(call {ind1.contig, ind1.pos, ind1.ALT, ref[ind1.contig][ind1.pos]});

				//mutate reference
				ref[ind1.contig][ind1.pos] = ind1.ALT;

			}

		}

	}

	ofstream new_ref;
	new_ref.open(new_ref_file);

	int line_length = 60;//line length in fasta

	for(string c : contigs){

		new_ref << ">" << c << endl;

		for(unsigned long i = 0;i<ref[c].length()/line_length + (ref[c].length()%line_length > 0);++i){

			new_ref << ref[c].substr(i*line_length,line_length) << endl;

		}


	}

	new_ref.close();

	ofstream differential_vcf;
	differential_vcf.open(differential_vcf_file);

	differential_vcf << "#CHROM\tPOS\tID\tREF\tALT" << endl;

	for(auto c : calls_vcfOut){

		differential_vcf << c.contig << "\t" << (c.pos+1) << "\t.\t" << c.REF << "\t" << c.ALT << endl;

	}

	cout << "done." << endl;

}
