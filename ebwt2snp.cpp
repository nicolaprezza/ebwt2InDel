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
#include <math.h>
#include <iomanip>
#include "dna_bwt.hpp"
#include "dna_string.hpp"
#include <stack>

using namespace std;

int k_left_def = 31;//extract k_left nucleotides from left of suffix array range, for each entry in the cluster (includes the SNP/indel)
int k_left = 0;//extract k_left nucleotides from left of suffix array range, for each entry in the cluster

int k_right_def = 30;//extract k_right nucleotides from right of suffix array range, only for entry with max LCP
int k_right = 0;//extract k_right nucleotides from right of suffix array range, for each entry in the cluster

double pval_def = 0.99; //choose cluster length so that this fraction of bases are captured (the others are discarded because inside a too small/big cluster)
double pval = 0;

int max_snvs_def = 3;//maximum number of SNVs allowed in left contexts (excluded main SNV).
int max_snvs = 0;//maximum number of SNVs allowed in left contexts

int mcov_out_def = 5;//minimum coverage required in the output events
int mcov_out = 0;//if a SNV is testified at least this number of times, then it is considered as a relevant event

int max_clust_length_def = 150;
int max_clust_length = 0;

//max indel length. If 0, indels are disabled.
int max_gap = 0;
int max_gap_def = 10;

int consensus_reads = 0;
int consensus_reads_def = 20;

//max tolerated errors in left-contexts while building consensus
int max_err_def = 2;
int max_err = 0;

string input1;
string input2;

uint64_t nr_reads1 = 0;//TODO remove

bool bcr = false;

bool discoSNP=true;

int lcp_def = 1;
int da_def = 4;
int pos_def = 1;

int lcp = 0;
int da = 0;
int pos = 0;

int K_def = 16;
int K = 0;

char TERM = '#';

uint64_t n_clust = 0; //number of clusters
uint64_t n_bases = 0; //number of bases in clusters

void help(){

	cout << "clust2snp [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help." << endl <<
	"-1 <arg>    Input eBWT file (A,C,G,T,#) of first individual (REQUIRED)." << endl <<
	"-2 <arg>    Input eBWT file (A,C,G,T,#) of second individual (REQUIRED)." << endl <<
	"-L <arg>    Length of left-context, SNP included (default: " << k_left_def << ")." << endl <<
	"-R <arg>    Length of right context, SNP excluded (default: " << k_right_def << ")." << endl <<
	"-k <arg>    Minimum LCP required in clusters (default: " << K_def << ")" << endl <<
	"-g <arg>    Maximum allowed gap length in indel (default: " << max_gap_def << "). If 0, indels are disabled."<< endl <<
	"-v <arg>    Maximum number of non-isolated SNPs in left-contexts. The central SNP/indel is excluded from this count (default: " << max_snvs_def << ")."<< endl <<
	"-c <arg>    Extract this maximum number of reads per individual to compute consensus of left-context (default: " << consensus_reads_def << ")."<< endl <<
	"-e <arg>    Mismatches allowed between DNA fragments forming consensus of left-context (default: " << max_err_def << ")."<< endl <<
	"-m <arg>    Minimum cluster length per individual (default: " << mcov_out_def << "). The minimum cluster length (for the 2 individuals) is 2*<arg>." <<  endl <<
	"-p <arg>    Automatically choose max cluster length so that this fraction of bases is analyzed (default: " << endl <<
	"            " << pval_def << "). In any case, the maximum cluster length will not exceed the value specified with -M."<< endl <<
	"-M <arg>    Maximum cluster length. Read the description of option -p." << endl <<
	"-x <arg>    Byte size of LCP integers in input EGSA/BCR file (default: " << lcp_def <<  ")." << endl <<
	"-y <arg>    Byte size of DA integers (read number) in input EGSA/BCR file (default: " << da_def <<  ")." << endl <<//TODO remove
	"-z <arg>    Byte size of pos integers (position in read) in input EGSA/BCR file (default: " << pos_def <<  ")." << endl << //TODO remove
	"-t <arg>    ASCII value of terminator character. Default: " << int('#') << " (#)." << endl << endl <<



	"\nTo run clust2snp, you must first build (1) the Enhanced Generalized Suffix Array of the input sequences" << endl <<
	"and the  cluster file built with ebwt2snp. Output is stored in reads.snp (this  is actually a fasta" << endl <<
	"file), where reads.fasta is the input fasta file." << endl << endl <<

	"Output:  SNPs are output in KisSNP2 format as a fasta file. IMPORTANT: in many cases, each SNP/indel is" << endl <<
	"reported twice: one time on the forward strand and one on the reverse strand. " << endl;

	exit(0);
}


/*
 * a pair of DNA segments (encoded as coordinates on reads) containing a potential variant between the
 * two individuals
 */
struct candidate_variant{

	//left contexts, of length k_left. Note: the variant is located at the end of the left context.
	//in the case of SNP, it is in the last letter of the left context.

	vector<uint64_t> left_context_idx_0; //indexes of the read containing the left context. Individual 0
	vector<uint64_t> left_context_pos_0; //starting positions of the left context in the read. Individual 0

	vector<uint64_t> left_context_idx_1; //index of the read containing the left context. Individual 1
	vector<uint64_t> left_context_pos_1; //starting position of the left context in the read. Individual 1

	//right contexdts, of length k_right

	uint64_t right_context_idx; //index of the read containing the left context (the same for both individuals)
	uint64_t right_context_pos; //starting position of the left context in the read

};


/*
 * a pair of DNA segment (this time encoded as strings) containing a potential variant between the
 * two individuals.
 */
struct variant_t{

	//left contexts are of length k_left. Note: the variant is located at the end of the left context.
	//in the case of SNP, it is in the last letter of the left context.

	string left_context_0;
	string left_context_1;

	string right_context;

	int support_0;
	int support_1;

};


/*
 * get a set of reads from file given their rank.
 *
 * reads must be sorted by rank!
 *
 * output: reads and their IDs
 */
void get_reads(string fasta_path, vector<uint64_t> & read_ranks, vector<string> & out_DNA){

	ifstream fasta;
	fasta.open(fasta_path);

	uint64_t j = 0;//current read in file

	string ID;
	string DNA;
	string line;

	//get read ID
	getline(fasta,ID);

	//read DNA
	getline(fasta,line);
	DNA = line;

	//skip newlines, append DNA segments
	while(not fasta.eof() && line[0] != '>'){

		getline(fasta,line);
		if(line[0] != '>') DNA.append(line);

	}

	int perc = 0;
	int last_perc = 0;

	cout << "(2/4) Extracting reads from fasta file ..." << endl;

	for(uint64_t i = 0; i<read_ranks.size();++i){

		while(j < read_ranks[i]){

			ID = line;

			//read DNA
			getline(fasta,line);
			DNA = line;
			while(not fasta.eof() && line[0] != '>'){

				getline(fasta,line);
				if(line[0] != '>') DNA.append(line);

			}

			++j;

		}

		out_DNA.push_back(DNA);

		perc = (i*100)/(read_ranks.size()-1);
		if(perc >= last_perc+10){

			last_perc=perc;
			cout << " " << perc << "% done." << endl;

		}

	}

	fasta.close();

}


/*
 * Hamming distance on strings. If length is different, align them on the right and discard extra chars on left.
 */
int dH(string & a, string & b){

	int len = a.size() < b.size() ? a.size() : b.size();

	int d=0;

	for(int i=0;i<len;++i){

		d += a[a.size()-i-1]!=b[b.size()-i-1];

	}

	return d;

}

/*
 * given two strings a, b of same length:
 *
 * 	1. find if there is an indel at the rightmost end (of max length max_gap)
 * 	2. skip the indel and count number of mismatches in the remaining part
 *
 *
 * 	returns a pair <D,L>, where:
 *
 * 	- D is the number of mismatches before the indel, and
 * 	- L is the indel length. This value is positive if the insertion is on a, and is negative if it is on b
 *
 *	best alignment is the one with minimum distance D + |L|
 *
 * 	examples:
 *
 * 	- distance(ACCTACTG, TTACTTAC) = <1,2>
 * 	- distance(TTACTTAC, ACCTACTG) = <1,-2>
 *
 */
pair<int,int> distance(string & a, string & b){

	assert(a.length()==b.length());

	auto dist_ab = vector<int>(max_gap,0);//insert in a
	auto dist_ba = vector<int>(max_gap,0);//insert in b
	auto dist_no_indel = dH(a,b);

	if(max_gap==0) return {dist_no_indel,0};

	assert(max_gap<=a.length());

	//try insert in a: remove characters from the right of a
	for(int i = 1; i<max_gap+1;++i){

		string a1 = a.substr(0,a.length()-i);
		dist_ab[i-1] = dH(a1,b) + i;

	}

	//try insert in b: remove characters from the right of b
	for(int i = 1; i<max_gap+1;++i){

		string b1 = b.substr(0,b.length()-i);
		dist_ba[i-1] = dH(a,b1) + i;

	}

	int min_ab_idx = std::distance(dist_ab.begin(), std::min_element(dist_ab.begin(),dist_ab.end()));
	int min_ba_idx = std::distance(dist_ba.begin(), std::min_element(dist_ba.begin(),dist_ba.end()));

	if(dist_no_indel < dist_ab[min_ab_idx] and  dist_no_indel < dist_ba[min_ba_idx]){

		//no indels
		return {dist_no_indel,0};

	}else if(dist_ab[min_ab_idx] < dist_ba[min_ba_idx]){

		//insert of length min_ab_idx+1 in a
		return {dist_ab[min_ab_idx] - (min_ab_idx+1), min_ab_idx+1};

	}

	//else

	//insert of length min_ba_idx+1 in b
	return {dist_ba[min_ba_idx] - (min_ba_idx+1), -(min_ba_idx+1)};

}

/*
 * return consensus letter among all i-th letters of strings in S.
 */
unsigned char consensus(vector<string> & S, uint64_t i){

	unsigned char c;

	auto nr = vector<uint64_t>(4,0);

	for(auto s:S){

		assert(i<s.length());
		nr[ base_to_int(s[i]) ]++;

	}

	return int_to_base( std::distance(nr.begin(), std::max_element(nr.begin(),nr.end())) );

}

/*
 * input: vector of strings with the same length
 *
 * output: consensus string. At each position, choose the most represented letter. In case of ties, choose arbitrarily.
 *
 */
string consensus_(vector<string> & S){

	if(S.size()==0) return "";
	if(S.size()==1) return S[0];

	string C = S[0];

	for(uint64_t i = 0;i<C.length();++i) C[i] = consensus(S,i);

	return C;

}

/*
 * input: vector of strings with the same length
 *
 * output: more advanced consensus. First, compute the consensus C with function consensus_. Then, count how many (x) strings lie within distance max_err
 * from C, and return <C,x>.
 *
 * return: consensus and a value indicating the number of strings from the original set S "supporting" the returned consensus.
 *
 */
pair<string, int> consensus(vector<string> & S){

	string C = consensus_(S);

	int support = 0;

	//compute how many strings in S have small Hamming distance from the computed consensus (i.e. they "support" the consensus)
	for(auto s:S)
		if(dH(C,s) <= max_err)
			support++;

	return {C, support};

}

vector<candidate_variant> find_variants(vector<t_GSA> & gsa_cluster){

	vector<candidate_variant>  out;

	auto counts = vector<vector<unsigned int> >(2,vector<unsigned int>(4,0));

	uint64_t max_lcp_val = 0;//value of max LCP in cluster
	uint64_t max_lcp_read_idx = 0;//index of read with max LCP in cluster
	uint64_t max_lcp_read_pos = 0;//position in read where max LCP starts

	for(uint64_t i=0;i<gsa_cluster.size();++i){

		auto e = gsa_cluster[i];

		//find read with max LCP
		if(e.lcp > max_lcp_val){

			max_lcp_val = e.lcp;
			max_lcp_read_idx = e.text;
			max_lcp_read_pos = e.suff;

		}

		bool sample = e.text < nr_reads1 ? 0 : 1;
		counts[sample?1:0][base_to_int(e.bwt)]++;

	}

	//discard cluster if max LCP is less than k_right
	if(max_lcp_val < k_right) return out;

	//compute the lists of frequent characters in indiv 0 and 1
	vector<unsigned char> frequent_char_0;
	vector<unsigned char> frequent_char_1;

	for(int c=0;c<4;++c){

		if(counts[0][c] >= mcov_out) frequent_char_0.push_back(int_to_base(c));
		if(counts[1][c] >= mcov_out) frequent_char_1.push_back(int_to_base(c));

	}

	std::sort(frequent_char_0.begin(), frequent_char_0.end());
	std::sort(frequent_char_1.begin(), frequent_char_1.end());

	//all variations observed in cluster
	auto all_chars = frequent_char_0;
	all_chars.insert(all_chars.begin(), frequent_char_1.begin(), frequent_char_1.end());
	std::sort( all_chars.begin(), all_chars.end() );
	all_chars.erase(std::unique( all_chars.begin(), all_chars.end() ), all_chars.end());

	//filter: remove clusters that cannot reflect a variation
	if(	frequent_char_0.size()==0 or // not covered enough
		frequent_char_1.size()==0 or // not covered enough
		frequent_char_0.size()>2 or // we require at most 2 alleles per individual
		frequent_char_1.size()>2 or // we require  at most 2 alleles per individual
		frequent_char_0 == frequent_char_1 or // same alleles: probably both heterozigous / multiple region (and no variants)
		all_chars.size() > 3	//4 or more distinct frequent characters in the cluster (probably multiple region)
	){

		return out;

	}

	for(auto c0 : frequent_char_0){

		for(auto c1 : frequent_char_1){

			if(c0 != c1){

				//compute max length of left context in indiv. 0 and 1, on the reads whose left
				//contexts end with c0 and c1, respectively.

				vector<uint64_t> left_pos_0;
				vector<uint64_t> left_idx_0;

				vector<uint64_t> left_pos_1;
				vector<uint64_t> left_idx_1;

				for(uint64_t i=0;i<gsa_cluster.size();++i){

					auto e = gsa_cluster[i];
					bool sample = e.text < nr_reads1 ? 0 : 1;
					uint64_t prefix_len = e.suff;
					unsigned char ch = e.bwt;
					uint64_t lcp = e.lcp;

					if(	prefix_len >= k_left and
						ch == c0 and
						sample == 0 and
						lcp >= k_right and //TODO test
						left_idx_0.size()<consensus_reads){

						left_idx_0.push_back(e.text);
						left_pos_0.push_back(e.suff-k_left);

					}

					if(	prefix_len >= k_left and
						ch == c1 and
						sample == 1 and
						lcp >= k_right and //TODO test
						left_idx_1.size()<consensus_reads){

						left_idx_1.push_back(e.text);
						left_pos_1.push_back(e.suff-k_left);

					}

				}

				if(left_idx_0.size()>0 and left_idx_1.size()>0){

					out.push_back(
						{

							left_idx_0, left_pos_0,
							left_idx_1, left_pos_1,
							max_lcp_read_idx, max_lcp_read_pos

						}
					);

				}

			}

		}

	}

	return out;

}

/*
 * extracts from the fasta file the DNA surrounding the variants
 */
vector<variant_t> extract_variants(vector<candidate_variant> & candidate_variants, string fasta_path){

	vector<variant_t> out;

	vector<uint64_t> read_ranks;

	//extract the ranks of all reads we need to process
	for(auto v : candidate_variants){

		read_ranks.insert(read_ranks.end(), v.left_context_idx_0.begin(), v.left_context_idx_0.end());
		read_ranks.insert(read_ranks.end(), v.left_context_idx_1.begin(), v.left_context_idx_1.end());
		read_ranks.push_back(v.right_context_idx);

	}

	//sort and remove duplicates
	std::sort( read_ranks.begin(), read_ranks.end() );
	auto last = std::unique( read_ranks.begin(), read_ranks.end() );
	read_ranks.erase(last, read_ranks.end());

	vector<string> reads;

	//get the reads as strings
	get_reads(fasta_path, read_ranks, reads);

	//invert read_ranks for fast access
	auto read_ranks_inv = vector<uint64_t>(read_ranks[read_ranks.size()-1]+1);

	for(uint64_t i=0;i<read_ranks.size();++i)
		read_ranks_inv[read_ranks[i]] = i;

	cout << "(3/4) Filtering " << candidate_variants.size() <<  " candidates and computing consensus of left-contexts ... " << endl;

	uint64_t idx=0;
	int perc = 0, last_perc=0;

	for(auto v:candidate_variants){

		//left 0
		cons left0(k_left);
		for(int j=0; j<v.left_context_idx_0.size(); ++j){

			uint64_t idx = read_ranks_inv[v.left_context_idx_0[j]];

			for(int i=0; i<k_left;++i)
				left0.increment(i,reads[idx][v.left_context_pos_0[j]+i]);

		}

		int supp0=0;

		for(int j=0; j<v.left_context_idx_0.size(); ++j){

			uint64_t idx = read_ranks_inv[v.left_context_idx_0[j]];

			//compute d_H
			int d_H=0;
			for(int i=0; i<k_left;++i)
				d_H += left0[i] != reads[idx][v.left_context_pos_0[j]+i];

			if(d_H <= max_err) supp0++;

		}

		//left 1
		cons left1(k_left);
		for(int j=0; j<v.left_context_idx_1.size(); ++j){

			uint64_t idx = read_ranks_inv[v.left_context_idx_1[j]];

			for(int i=0; i<k_left;++i)
				left1.increment(i,reads[idx][v.left_context_pos_1[j]+i]);

		}

		int supp1=0;

		for(int j=0; j<v.left_context_idx_1.size(); ++j){

			uint64_t idx = read_ranks_inv[v.left_context_idx_1[j]];

			//compute d_H
			int d_H=0;
			for(int i=0; i<k_left;++i)
				d_H += left1[i] != reads[idx][v.left_context_pos_1[j]+i];

			if(d_H <= max_err) supp1++;

		}

		if(supp0 > 0 and supp1 > 0){

			uint64_t r_idx = read_ranks_inv[v.right_context_idx];

			out.push_back(

				{
					left0.to_string(),
					left1.to_string(),
					reads[r_idx].substr(v.right_context_pos,k_right),
					supp0,
					supp1
				}

			);

		}

		++idx;

		perc = (idx*100)/candidate_variants.size();
		if(perc >= last_perc+1){

			last_perc=perc;
			cout << " " << perc << "% done." << endl;

		}


	}

	return out;

}

/*
 * detect the type of variant (SNP/indel/discard if none) and, if not discarded, output to file the two reads per variant testifying it.
 */
void to_file(vector<variant_t> & output_variants, string & out_path){

	ofstream out_file = ofstream(out_path);

	uint64_t id_nr = 1;
	uint64_t idx = 0;

	int perc = 0;
	int last_perc = 0;

	cout << "(4/4) Computing edit distances and saving SNPs/indels to file ... " << endl;
	for(auto v:output_variants){

		auto d = distance(v.left_context_0,v.left_context_1);

		if(d.first <= max_snvs_def){

			/*
			 * sample 1
			 */

			string ID;

			if(d.second != 0){

				ID =  ">INDEL_higher_path_";

			}else{

				ID =  ">SNP_higher_path_";

			}

			ID.append(to_string(id_nr));
			ID.append("|P_1:");
			ID.append(to_string(v.right_context.size()));
			ID.append("_");

			string snv_type;

			if(d.second==0){

				snv_type += v.left_context_0[v.left_context_0.size()-1];
				snv_type.append("/");
				snv_type += v.left_context_1[v.left_context_1.size()-1];

			}else if(d.second>0){//insert of length d.second in v.left_context_0

				snv_type.append(v.left_context_0.substr(v.left_context_0.size()-d.second));
				snv_type.append("/");

			}else{//insert of length -d.second in v.left_context_1

				snv_type.append("/");
				snv_type.append(v.left_context_1.substr(v.left_context_1.size() - (-d.second)));

			}

			ID.append(snv_type);
			ID.append("|");
			ID.append(std::to_string(v.support_0));//we write the number of reads supporting this variant
			//ID.append("high");
			ID.append("|nb_pol_1");

			out_file << ID << endl;

			string DNA;

			if(d.second==0){

				DNA = v.left_context_0;

			}else if(d.second>0){//insert of length d.second in v.left_context_0

				DNA = v.left_context_0;

			}else{//insert of length -d.second in v.left_context_1

				DNA = v.left_context_0.substr(-d.second);

			}

			DNA.append(v.right_context);

			out_file << DNA << endl;

			/*
			 * sample 2
			 */

			if(d.second != 0){

				ID =  ">INDEL_lower_path_";

			}else{

				ID =  ">SNP_lower_path_";

			}

			ID.append(to_string(id_nr));
			ID.append("|P_1:");
			ID.append(to_string(v.right_context.size()));
			ID.append("_");
			ID.append(snv_type);
			ID.append("|");
			ID.append(std::to_string(v.support_1));
			//ID.append("high");
			ID.append("|nb_pol_1");

			out_file << ID << endl;

			DNA = "";

			if(d.second==0){

				DNA = v.left_context_1;

			}else if(d.second>0){//insert of length d.second in v.left_context_0

				DNA = v.left_context_1.substr(d.second);

			}else{//insert of length -d.second in v.left_context_1

				DNA = v.left_context_1;

			}

			DNA.append(v.right_context);
			out_file << DNA << endl;

			id_nr++;

		}

		idx++;

		perc = (idx*100)/output_variants.size();
		if(perc >= last_perc+10){

			last_perc=perc;
			cout << " " << perc << "% done." << endl;

		}

	}

}



/*
 * scans EGSA, clusters and finds interesting clusters. In chunks, extracts the reads
 * from interesting clusters and aligns them.
 */
void find_events(egsa_stream & EGSA, string & clusters_path, string fasta_path, string out_path){

	ifstream clusters;
	clusters.open(clusters_path, ios::in | ios::binary);

	uint64_t i = 0;//position on suffix array

	//read first egsa entry
	t_GSA e = EGSA.read_el();

	vector<candidate_variant> candidate_variants;

	cout << "(1/4) Filtering relevant clusters ... " << endl;

	uint64_t cl = 0;
	int perc=0;
	int last_perc=0;

	while(not clusters.eof()){

		//1. EXTRACT EGSA CLUSTER

		uint64_t start;
		uint16_t length;

		clusters.read((char*)&start, sizeof(uint64_t));
		clusters.read((char*)&length, sizeof(uint16_t));

		if(length >= mcov_out*2 and length <= max_clust_length){

			while(i < start){

				e = EGSA.read_el();
				++i;

			}

			vector<t_GSA> gsa_cluster;

			while(i < start+length){

				gsa_cluster.push_back(e);
				e = EGSA.read_el();
				++i;

			}

			//now gsa_cluster contains a cluster in the egsa

			//2. EXTRACT EVENTS FROM EGSA CLUSTER

			//find potential variants
			auto v = find_variants(gsa_cluster);

			//append them to the vector of all candidate variants
			candidate_variants.insert(candidate_variants.end(), v.begin(), v.end());

		}

		cl++;

		perc = (cl*100)/(n_clust-1);
		if(perc >= last_perc+10){

			last_perc=perc;
			cout << " " << perc << "% done." << endl;

		}

	}

	cout << "Done. "  << candidate_variants.size() << " potential variants detected (some might be detected twice: on fw and rev strands)" << endl;

	//3. EXTRACT READ SEGMENTS FROM FILE
	//extract from file the interesting parts of the reads and form the variants to be outputted

	vector<variant_t> output_variants = extract_variants(candidate_variants, fasta_path);

	//4. SAVE TO OUTPUT FILE THE VARIANTS

	to_file(output_variants, out_path);

	clusters.close();

}

/*
 * compute coverage statistics, auto-compute max cluster length
 */
void statistics(string & clusters_path){

	ifstream clusters;
	clusters.open(clusters_path, ios::in | ios::binary);

	uint64_t MAX_C_LEN = max_clust_length;

	//init with max cluster length 200
	auto clust_len_freq = vector<uint64_t>(MAX_C_LEN+1,0);

	uint64_t max_len = 0;

	while(not clusters.eof()){

		//read entry in clusters

		uint64_t start;
		uint16_t length;

		clusters.read((char*)&start, sizeof(uint64_t));
		clusters.read((char*)&length, sizeof(uint16_t));

		if(length <= MAX_C_LEN){

			clust_len_freq[length]++;
			max_len = length>max_len ? length : max_len;

		}

		n_clust++;
		n_bases+=length;

	}

	uint64_t max = 0;
	for(int i=1;i<=MAX_C_LEN;++i) max = clust_len_freq[i]*i > max ? clust_len_freq[i]*i : max;

	uint64_t cumulative = 0;

	cout << "\nDistribution of base coverage: "<< endl;
	cout << "\ncluster length\t# bases in a cluster with this length\t cumulative fraction (from 2m = " << 2*mcov_out << ")" << endl;
	for(int i=0;i<=max_len;++i){

		cout << i << "\t" << flush;
		for(uint64_t j=0;j<(100*clust_len_freq[i]*i)/max;++j) cout << "-" << flush;
		cout << "\t" << clust_len_freq[i]*i << flush;

		if(i>=2*mcov_out){

			cumulative += clust_len_freq[i]*i;
			cout << "\t" << double(cumulative)/double(n_bases);

		}

		cout << endl;


	}

	//auto-detect max cluster length

	max_clust_length = 2*mcov_out;//start from the minimum cluster length
	cumulative = clust_len_freq[max_clust_length]*max_clust_length;//cumulative number of bases

	while( double(cumulative)/double(n_bases) < pval and  max_clust_length < MAX_C_LEN){

		max_clust_length++;
		cumulative += clust_len_freq[max_clust_length]*max_clust_length;

	}

/*	max = 0;
	for(int i=1;i<=MAX_C_LEN;++i) max = clust_len_freq[i] > max ? clust_len_freq[i] : max;

	cout << "\nDistribution of cluster length: "<< endl;
	cout << "\ncluster length\t# clusters with this length" << endl;
	for(int i=0;i<=max_len;++i){

		cout << i << "\t" << flush;
		for(uint64_t j=0;j<(100*clust_len_freq[i])/max;++j) cout << "-" << flush;
		cout << "   " << clust_len_freq[i] << endl;


	}*/

	cout << "\nCluster sizes allowed: [" << mcov_out*2 << "," << max_clust_length << "]" << endl;

	clusters.close();

}

void update_DA(sa_leaf L1,sa_leaf L2, vector<bool> & DA, vector<bool> & LCP_threshold, uint64_t & lcp_values, uint64_t & m){

	uint64_t start1 = L1.rn.first + L2.rn.first;//start position of first interval in merged intervals
	uint64_t start2 = L2.rn.first + L1.rn.second;//start position of second interval in merged intervals
	uint64_t end = L1.rn.second + L2.rn.second;//end position (excluded) of merged intervals

	assert(end>start1);

	for(uint64_t i = start1; i<start2; ++i){
		DA[i] = 0;
		m++;
	}

	for(uint64_t i = start2; i<end; ++i){
		DA[i] = 1;
		m++;
	}

	assert(L1.depth==L2.depth);

	for(uint64_t i = start1+1; i<end; ++i){

		LCP_threshold[i] = (L1.depth >= K);

		lcp_values++;

	}

}

void update_DA(sa_leaf L1,sa_leaf L2, vector<bool> & DA, uint64_t & m){

	uint64_t start1 = L1.rn.first + L2.rn.first;//start position of first interval in merged intervals
	uint64_t start2 = L2.rn.first + L1.rn.second;//start position of second interval in merged intervals
	uint64_t end = L1.rn.second + L2.rn.second;//end position (excluded) of merged intervals

	assert(end>start1);

	for(uint64_t i = start1; i<start2; ++i){
		DA[i] = 0;
		m++;
	}

	for(uint64_t i = start2; i<end; ++i){
		DA[i] = 1;
		m++;
	}

	assert(L1.depth==L2.depth);

}

void next_leaves(dna_bwt_t & bwt1, dna_bwt_t & bwt2, sa_leaf & L1, sa_leaf & L2, vector<pair<sa_leaf, sa_leaf> > & TMP_LEAVES, int & t){

	p_range ext_1 = bwt1.LF(L1.rn);
	p_range ext_2 = bwt2.LF(L2.rn);

	//push non-empty leaves on stack in decreasing size order

	t = 0;
	int min_size = 2;

	if(range_length(ext_1.A) + range_length(ext_2.A) >= min_size) TMP_LEAVES[t++] = {{ext_1.A, L1.depth+1},{ext_2.A, L2.depth+1}};
	if(range_length(ext_1.C) + range_length(ext_2.C) >= min_size) TMP_LEAVES[t++] = {{ext_1.C, L1.depth+1},{ext_2.C, L2.depth+1}};
	if(range_length(ext_1.G) + range_length(ext_2.G) >= min_size) TMP_LEAVES[t++] = {{ext_1.G, L1.depth+1},{ext_2.G, L2.depth+1}};
	if(range_length(ext_1.T) + range_length(ext_2.T) >= min_size) TMP_LEAVES[t++] = {{ext_1.T, L1.depth+1},{ext_2.T, L2.depth+1}};

	std::sort( TMP_LEAVES.begin(), TMP_LEAVES.begin()+t, [ ]( const pair<sa_leaf, sa_leaf>& lhs, const pair<sa_leaf, sa_leaf>& rhs )
	{
		return leaf_size(lhs) < leaf_size(rhs);
	});

}


void find_leaves(sa_node x1, sa_node x2, vector<bool> & DA, uint64_t & m){

	//find leaves that were ignored in the first pass
	if(range_length(child_TERM(x1))+range_length(child_TERM(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_TERM(x1),0};
		sa_leaf L2 = {child_TERM(x2),0};

		update_DA(L1,L2,DA,m);

	}

	if(range_length(child_A(x1))+range_length(child_A(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_A(x1),0};
		sa_leaf L2 = {child_A(x2),0};

		update_DA(L1,L2,DA,m);

	}

	if(range_length(child_C(x1))+range_length(child_C(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_C(x1),0};
		sa_leaf L2 = {child_C(x2),0};

		update_DA(L1,L2,DA,m);

	}

	if(range_length(child_G(x1))+range_length(child_G(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_G(x1),0};
		sa_leaf L2 = {child_G(x2),0};

		update_DA(L1,L2,DA,m);

	}

	if(range_length(child_T(x1))+range_length(child_T(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_T(x1),0};
		sa_leaf L2 = {child_T(x2),0};

		update_DA(L1,L2,DA,m);

	}

}

void next_nodes(dna_bwt_t & bwt1, dna_bwt_t & bwt2, sa_node & N1, sa_node & N2, vector<pair<sa_node, sa_node> > & TMP_NODES, int & t){

	p_node left_exts1 = bwt1.LF(N1);
	p_node left_exts2 = bwt2.LF(N2);

	pair<sa_node, sa_node> A = {left_exts1.A, left_exts2.A};
	pair<sa_node, sa_node> C = {left_exts1.C, left_exts2.C};
	pair<sa_node, sa_node> G = {left_exts1.G, left_exts2.G};
	pair<sa_node, sa_node> T = {left_exts1.T, left_exts2.T};

	t = 0;

	if(number_of_children(A) >= 2) TMP_NODES[t++] = A;
	if(number_of_children(C) >= 2) TMP_NODES[t++] = C;
	if(number_of_children(G) >= 2) TMP_NODES[t++] = G;
	if(number_of_children(T) >= 2) TMP_NODES[t++] = T;

	//push right-maximal nodes on stack in decreasing size (i.e. interval length) order

	std::sort( TMP_NODES.begin(), TMP_NODES.begin()+t, [ ]( const pair<sa_node, sa_node>& lhs, const pair<sa_node, sa_node>& rhs )
	{
		return node_size(lhs) < node_size(rhs);
	});

}

void update_lcp_minima(sa_node x, vector<bool> & LCP_minima, uint64_t & n_min){

/*
 * we have a minima at the beginning of each child of size at least 2 of the input node
 */

	//x.first_A is the end position (excluded) of terminator.

	if(x.first_G - x.first_C >= 2 and x.first_C > x.first_A){
		LCP_minima[x.first_C] = true;
		n_min++;
	}

	if(x.first_T - x.first_G >= 2 and x.first_G > x.first_A){
		LCP_minima[x.first_G] = true;
		n_min++;
	}

	if(x.last - x.first_T >= 2 and x.first_T > x.first_A){
		LCP_minima[x.first_T] = true;
		n_min++;
	}

}

int main(int argc, char** argv){

	srand(time(NULL));

	if(argc < 3) help();

	int opt;
	while ((opt = getopt(argc, argv, "h1:2:p:v:L:R:m:g:c:x:y:z:e:k:t:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'b':
				bcr = true;
			break;
			case '1':
				input1 = string(optarg);
			break;
			case '2':
				input2 = string(optarg);
			break;
			case 'm':
				mcov_out = atoi(optarg);
			break;
			case 'k':
				K = atoi(optarg);
			break;
			case 'g':
				max_gap = atoi(optarg);
			break;
			case 'M':
				max_clust_length = atoi(optarg);
			break;
			case 'L':
				k_left = atoi(optarg);
			break;
			case 'c':
				consensus_reads = atoi(optarg);
			break;
			case 'R':
				k_right = atoi(optarg);
			break;
			case 'p':
				pval = atof(optarg);
			break;
			case 'v':
				max_snvs = atoi(optarg);
			break;
			case 'e':
				max_err = atof(optarg);
			break;
			case 'x':
				lcp = atoi(optarg);
			break;
			case 't':
				TERM = atoi(optarg);
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

	K = K == 0 ? K_def : K;
	lcp = lcp==0?lcp_def:lcp;
	da = da==0?da_def:da;
	pos = pos==0?pos_def:pos;
	max_err = max_err==0?max_err_def:max_err;
	consensus_reads = consensus_reads==0?consensus_reads_def:consensus_reads;
	max_gap = max_gap==0?max_gap_def:max_gap;
	max_clust_length = max_clust_length==0?max_clust_length_def:max_clust_length;
	k_left = k_left==0?k_left_def:k_left;
	k_right = k_right==0?k_right_def:k_right;
	pval = pval==0?pval_def:pval;
	max_snvs = max_snvs==0?max_snvs_def:max_snvs;
	mcov_out = mcov_out==0?mcov_out_def:mcov_out;

	if(input1.compare("")==0 or input2.compare("")==0) help();

	cout << "This is clust2snp, version 2." << endl <<
			"Input eBWT files : " << input1 << " and " << input2 << endl <<
			"Left-extending eBWT ranges by " << k_left << " bases." << endl <<
			"Right context length: " << k_right << " bases." << endl << endl;

	cout << "Phase 1/4: loading and indexing eBWTs ... " << flush;

	dna_bwt_t bwt1 = dna_bwt_t(input1,TERM);
	dna_bwt_t bwt2 = dna_bwt_t(input2,TERM);

	cout << "done." << endl;

	uint64_t n1 = bwt1.size();
	uint64_t n2 = bwt2.size();
	uint64_t n = n1+n2;

	cout << "\nPhase 2/4: merging eBWTs." << endl;

	vector<bool> DA = vector<bool>(n,false); //document array
	vector<bool> LCP_threshold = vector<bool>(n,false);//1 where LCP >= K

	uint64_t da_values = 0;//number computed DA values
	uint64_t leaves = 0;//number of visited leaves
	uint64_t max_stack = 0;
	uint64_t lcp_values = 1;//number of computed LCP values

	{

		auto TMP_LEAVES = vector<pair<sa_leaf, sa_leaf> >(5);

		stack<pair<sa_leaf, sa_leaf> > S;
		S.push({bwt1.first_leaf(), bwt2.first_leaf()});

		int last_perc_lcp = -1;
		int last_perc_da = -1;
		int perc_lcp = 0;
		int perc_da = 0;

		while(not S.empty()){

			pair<sa_leaf, sa_leaf> L = S.top();
			S.pop();
			leaves++;

			assert(leaf_size(L)>0);
			max_stack = S.size() > max_stack ? S.size() : max_stack;

			sa_leaf L1 = L.first;
			sa_leaf L2 = L.second;

			update_DA(L1, L2, DA, LCP_threshold, lcp_values, da_values);

			//insert leaf in stack iff size(L1) + size(L2) >= min_size
			//optimization: if we are computing LCP and if size(L1) + size(L2) = 1,
			//then we will find that leaf during the internal nodes traversal (no need to visit leaf here)
			int t = 0;//number of children leaves
			next_leaves(bwt1, bwt2, L1, L2, TMP_LEAVES, t);
			for(int i=t-1;i>=0;--i) S.push(TMP_LEAVES[i]);

			perc_lcp = (100*lcp_values)/n;
			perc_da = (100*da_values)/n;

			if(perc_da > last_perc_da){

				cout << "DA: " << perc_da << "%. ";
				cout << "LCP: " << perc_lcp << "%.";
				cout << endl;

				last_perc_lcp = perc_lcp;
				last_perc_da = perc_da;

			}

		}
	}

	cout << "Computed " << da_values << "/" << n << " DA values." << endl;
	cout << "Computed " << lcp_values << "/" << n << " LCP threshold values." << endl;

	cout << "Max stack depth = " << max_stack << endl;
	cout << "Processed " << leaves << " suffix-tree leaves." << endl << endl;

	cout << "Phase 3/4: computing LCP minima." << endl;

	vector<bool> LCP_minima = vector<bool>(n,false);

	auto TMP_NODES = vector<pair<sa_node, sa_node> >(4);

	uint64_t nodes = 0;//visited ST nodes
	max_stack = 0;

	stack<pair<sa_node, sa_node> > S;
	S.push({bwt1.root(), bwt2.root()});

	int last_perc_lcp = -1;
	int perc_lcp = 0;
	int last_perc_da = -1;
	int perc_da = 0;
	uint64_t n_min = 0;//number of LCP minima

	while(not S.empty()){

		max_stack = S.size() > max_stack ? S.size() : max_stack;

		pair<sa_node, sa_node> N = S.top();
		S.pop();
		nodes++;

		sa_node N1 = N.first;
		sa_node N2 = N.second;
		sa_node merged = merge_nodes(N1, N2);

		//find leaves in the children of N1 and N2 that were
		//skipped in the first pass, and update DA accordingly
		find_leaves(N1, N2, DA, da_values);

		//compute LCP values at the borders of merged's children
		update_lcp_threshold(merged, LCP_threshold, lcp_values, K);

		update_lcp_minima(merged, LCP_minima, n_min);

		//follow Weiner links
		int t = 0;
		next_nodes(bwt1, bwt2, N1, N2, TMP_NODES, t);
		for(int i=t-1;i>=0;--i) S.push(TMP_NODES[i]);

		perc_lcp = (100*lcp_values)/n;
		perc_da = (100*da_values)/n;

		if(perc_da > last_perc_da){

			cout << "DA: " << perc_da << "%. ";
			cout << "LCP: " << perc_lcp << "%.";
			cout << endl;

			last_perc_lcp = perc_lcp;
			last_perc_da = perc_da;

		}

	}

	cout << "Computed " << da_values << "/" << n << " DA values." << endl;
	cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;
	cout << "Found " << n_min << " LCP minima." << endl;
	cout << "Max stack depth = " << max_stack << endl;
	cout << "Processed " << nodes << " suffix-tree nodes." << endl << endl;


	cout << "Phase 4/4: detecting SNPs and indels." << endl;


}
