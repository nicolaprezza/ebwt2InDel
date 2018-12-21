/*
 * merge_bwt.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */


#include <iostream>
#include <fstream>
#include "internal/dna_bwt.hpp"
#include "internal/bwt.hpp"
#include "internal/bwt_merger.hpp"

using namespace std;

string input_bwt1;
string input_bwt2;
string output_file;

bool out_da = false;
uint8_t lcp_size = 0;

void help(){

	cout << "merge_bwt [options]" << endl <<
	"Merges eBWTs of two read sets. Only A,C,G,T, and terminator # are allowed in the input eBWTs." <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-1 <arg>    Input BWT index 1 (REQUIRED)" << endl <<
	"-2 <arg>    Input BWT index 2 (REQUIRED)" << endl <<
	"-o <arg>    Output prefix (REQUIRED)" << endl <<
	"-d          Output document array as an ASCII file of 0/1. Default: do not output." << endl <<
	"-l <arg>    Output LCP of the merged BWT using <arg>=0,1,2,4,8 Bytes" << endl <<
	"            per integer. If arg=0, LCP is not computed (faster). Default: 0." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 4) help();

	int opt;
	while ((opt = getopt(argc, argv, "h1:2:o:l:d")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case '1':
				input_bwt1 = string(optarg);
			break;
			case '2':
				input_bwt2 = string(optarg);
			break;
			case 'o':
				output_file = string(optarg);
			break;
			case 'l':
				lcp_size = atoi(optarg);
			break;
			case 'd':
				out_da = true;
			break;
			default:
				help();
			return -1;
		}
	}

	if(input_bwt1.size()==0) help();
	if(input_bwt2.size()==0) help();
	if(output_file.size()==0) help();

	cout << "Input bwt index file 1: " << input_bwt1 << endl;
	cout << "Input bwt index file 2: " << input_bwt2 << endl;
	cout << "Output prefix: " << output_file << endl;

	cout << "Loading BWTs ... " << endl;

	dna_bwt_t BWT1;
	dna_bwt_t BWT2;

	BWT1.load_from_file(input_bwt1);
	BWT2.load_from_file(input_bwt2);

	cout << "Done. Size of BWTs: " << BWT1.size() << " and " << BWT2.size() << endl;

	switch(lcp_size){

		case 0: { bwt_merger<dna_bwt_t,dna_bwt_t, uint8_t> M0(&BWT1, &BWT2, false, out_da);
				cout << "Storing output to file ... " << endl;
				M0.save_to_file(output_file);
				break; }
		case 1: { bwt_merger<dna_bwt_t,dna_bwt_t, uint8_t> M1(&BWT1, &BWT2, true, out_da);
				cout << "Storing output to file ... " << endl;
				M1.save_to_file(output_file);
				break;}
		case 2: {bwt_merger<dna_bwt_t,dna_bwt_t, uint16_t> M2(&BWT1, &BWT2, true, out_da);
				cout << "Storing output to file ... " << endl;
				M2.save_to_file(output_file);
				break;}
		case 4: {bwt_merger<dna_bwt_t,dna_bwt_t, uint32_t> M4(&BWT1, &BWT2, true, out_da);
				cout << "Storing output to file ... " << endl;
				M4.save_to_file(output_file);
				break;}
		case 8: {bwt_merger<dna_bwt_t,dna_bwt_t, uint64_t> M8(&BWT1, &BWT2, true, out_da);
				cout << "Storing output to file ... " << endl;
				M8.save_to_file(output_file);
				break;}
		default:break;

	}

	cout << "Done. " << endl;

}

