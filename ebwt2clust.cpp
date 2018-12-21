// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include "include.hpp"
#include <unistd.h>

using namespace std;

/*
 * parameters
 *
 */
int min_def = 2; //Discard clusters smaller than this value
int K_def = 16; //require an LCP of at least k inside clusters
int k = 0;
string input;

int lcp_def = 1;
int da_def = 4;
int pos_def = 1;

int lcp = 0;
int da = 0;
int pos = 0;


int min_len=0;

void help(){

	cout << "ebwt2clust [options]" << endl <<
	"Options:" << endl <<
	"-h         Print this help" << endl <<
	"-i <arg>   Input fasta file (REQUIRED)" << endl <<
	"-k <arg>   Minimum LCP required in clusters (default: " << K_def << ")" << endl <<
	"-m <arg>   Discard clusters smaller than this value (default: " << min_def << ")" << endl <<
	"-x <arg>   Byte size of LCP integers in input EGSA/BCR file (default: " << lcp_def <<  ")." << endl <<
	"-y <arg>   Byte size of DA integers (read number) in input EGSA/BCR file (default: " << da_def <<  ")." << endl <<
	"-z <arg>   Byte size of pos integers (position in read) in input EGSA/BCR file (default: " << pos_def <<  ")." << endl << endl <<

	"\nTo run ebwt2clust, you must  first build the Enhanced Generalized  Suffix Array of the input" << endl <<
	"sequences. The EGSA must be stored in the input file's folder adding extension .gesa to the" << endl <<
	"name of the input file (github.com/felipelouza/egsa), or in three files with extensions" << endl <<
	".out, .out.lcp, .out.pairSA computed using the BCR algorithm " << endl <<
	"(https://github.com/giovannarosone/BCR_LCP_GSA). Output is stored in reads.fasta.clusters." << endl;
	 exit(0);
}

void append_entry(ofstream & out, uint64_t start, uint16_t length){

	if(length >=min_len){

		out.write((char*)&start, sizeof(uint64_t));
		out.write((char*)&length, sizeof(uint16_t));

	}

}

/*
 * clusters = regions between local LCP minima (excluding tails where LCP < k)
 */
void cluster_lm(egsa_stream & EGSA,ofstream & out){

	uint64_t null = ~uint64_t(0);

	uint64_t start = null;		//start position of cluster
	uint64_t i = 0;			//current BWT position
	uint64_t off = 0;		//current offset in cluster

	uint64_t prev_lcp = null;//previous LCP value

	//find local minima in the LCP
	t_GSA e1 = EGSA.read_el();
	t_GSA e2 = EGSA.read_el();
	t_GSA e3;

	start = e1.lcp >= k ? 0 :
			e2.lcp >= k ? 1 : null;

	i = 1;//index of e2

	unsigned int n_clust_out = 0;//number of clusters in output

	while(not EGSA.eof()){

		//read next value
		e3 = EGSA.read_el();

		//cout << e3.text << " " << e3.suff << " " << e3.lcp << " " << char(e3.bwt) << endl;

		//e2 is start of a (possibly flat) local minima
		if(		start != null and
					(		(e1.lcp > e2.lcp and e2.lcp <= e3.lcp) or
							e3.lcp < k
					)
					){

			uint16_t length = (i - start) + 1; //this cluster ends in e2
			append_entry(out, start, length);
			n_clust_out++;

			start = null;

		}

		e1 = e2;
		e2 = e3;
		++i;

		if(start==null and e2.lcp >= k){

			start = i;

		}

	}

	/*
	 * check if last cluster was closed. If not, output it
	 */
	if(start != null){

		uint16_t length = (i - start) + 1; //this cluster ends in e2
		append_entry(out, start, length);
		n_clust_out++;

		start = null;

	}

	cout << "Done. " << n_clust_out << " clusters saved to output file." << endl;

}

int main(int argc, char** argv){

	if(argc < 2) help();

	int opt;
	while ((opt = getopt(argc, argv, "hk:i:m:x:y:z:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'k':
				k = atoi(optarg);
			break;
			case 'm':
				min_len = atoi(optarg);
			break;
			case 'i':
				input = string(optarg);
			break;
			case 'x':
				lcp = atoi(optarg);
			break;
			case 'y':
				da = atoi(optarg);
			break;
			case 'z':
				pos = atoi(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	lcp = lcp==0?lcp_def:lcp;
	da = da==0?da_def:da;
	pos = pos==0?pos_def:pos;

	k = k==0?K_def:k;
	min_len = min_len==0?min_def:min_len;

	if(input.compare("")==0) help();

	egsa_stream EGSA(input);
	EGSA.set_bytesizes(lcp,da,pos);

	cout << "This is ebwt2clust. Input file: " << input << endl;

	string filename_out = input;
	filename_out.append(".clusters");
	ofstream out;
	out.open(filename_out, ios::out | ios::binary);

	cluster_lm(EGSA,out);

	out.close();

}
