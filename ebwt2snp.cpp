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

int max_snvs_def = 2;//maximum number of SNVs allowed in left contexts (excluded main SNV).
int max_snvs = 0;//maximum number of SNVs allowed in left contexts

int mcov_out_def = 6;//minimum cluster length / coverage of output events
int mcov_out = 0;

//max indel length. If 0, indels are disabled.
int max_gap = 0;
int max_gap_def = 10;

string input1;
string input2;
string input_da;//document array
string output;

bool bcr = false;

bool diploid = false;

bool discoSNP=true;

int K_def = 16;
int K = 0;

char TERM = '#';

uint64_t cluster_nr = 1;//cluster number while outputting to file

const uint64_t MAX_READ_LEN = 100000;

vector<bool> LCP_minima;
vector<bool> LCP_threshold;
vector<bool> DA;

//parameter -c
int complexity_def = k_right_def-10 < 0 ? 0 : k_right_def-10;
int complexity = 0;

//vector<uint8_t> LCP;

uint64_t n_bases = 0; //number of bases in clusters
uint64_t events = 0;//number of outputted events

uint64_t id_nr = 1;

void help(){

	cout << "ebwt2snp [options]" << endl <<
	"Options:" << endl <<
	"-h          Print this help." << endl <<
	"-1 <arg>    Input eBWT file (A,C,G,T,#) of first reads set (REQUIRED)." << endl <<
	"-2 <arg>    Input eBWT file (A,C,G,T,#) of second reads set. If not specified, perform genotyping of first reads set." << endl <<
	"            If specified, find differences (SNPs/indels) between the two reads sets." << endl <<
	"-d <arg>    Input Document Array. If option -2 is not specified, this file specifies which characters from the input bwt" << endl <<
	"            belong to the first (0) and which from the second (1) individual. Format: ASCII file filled with '0' and '1'." << endl <<
	"-o <arg>    Output .snp file (REQUIRED)." << endl <<
	"-L <arg>    Length of left-context, SNP included (default: " << k_left_def << ")." << endl <<
	"-R <arg>    Length of right context, SNP excluded (default: " << k_right_def << ")." << endl <<
	"-k <arg>    Minimum LCP required in clusters (default: " << K_def << ")" << endl <<
	"-g <arg>    Maximum allowed gap length in indel (default: " << max_gap_def << "). If 0, indels are disabled."<< endl <<
	"-v <arg>    Maximum number of non-isolated SNPs in left-contexts. The central SNP/indel is excluded from this count (default: " << max_snvs_def << ")."<< endl <<
	"-m <arg>    Minimum coverage of output events (default: " << mcov_out_def << ")." <<  endl <<
	"-c <arg>    Discard events with low-complexity right-context. Here, low-complexity means that the context starts with a run of <arg> equal characters." << endl <<
	"            Default: length of right context (-R), minus 10." << endl <<
	"-D          Samples are diploid (default: haploid). Not effective during genotyping (i.e. only -1 specified)." <<  endl <<
	"-t <arg>    ASCII value of terminator character. Default: " << int('#') << " (#)." << endl << endl <<

	"\nTo run ebwt2snp, you must first build the extended Burrows-Wheeler Transform of the input sequences." << endl << endl <<

	"Output format: A fasta file with DNA fragments containing the variations." << endl;

	exit(0);
}

bool file_exists(string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

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

struct variant_single_t{

	//left contexts are of length k_left. Note: the variant is located at the end of the left context.
	//in the case of SNP, it is in the last letter of the left context.

	string left_context;

	string right_context;

	int support;

};

//true iff s starts with a run of >= k equal characters
bool has_run(string & s, int k){

	if(k>s.length()) return false;

	for(int i=1;i<k;++i) if(s[i]!=s[i-1]) return false;

	return true;

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

//find most frequent letter
pair<char,range_t> consensus_letter(p_range pr){

	vector<pair<char,range_t> > freqs;

	freqs.push_back({'A',pr.A});
	freqs.push_back({'C',pr.C});
	freqs.push_back({'G',pr.G});
	freqs.push_back({'T',pr.T});

	sort( freqs.begin( ), freqs.end( ), [ ]( const pair<char,range_t>& lhs, const pair<char,range_t>& rhs )
	{
	   return range_length(lhs.second) > range_length(rhs.second);
	});

	if(range_length(freqs[0].second) == 0) return {0,{0,0}};

	return freqs[0];

}

//extract len characters using backward search from range. At each step, most frequent character is chosen.
//note: the reversed string is returned
void extract_consensus(dna_bwt_t & bwt, range_t R,string & ctx, int & freq, int len){

	if(len==0) return;

	p_range pr = bwt.LF(R);
	auto c = consensus_letter(pr);

	freq = 0;

	if(c.first!=0){

		ctx += c.first;
		freq = range_length(c.second);
		extract_consensus(bwt, c.second, ctx,freq, len-1);

	}

}

//extract len characters using backward search, starting from position in range containing character c
/*void extract_consensus(dna_bwt_t & bwt, range_t range, char c, vector<pair<string, int> > & left_contexts, int len){

	string ctx;
	ctx += c;

	range_t R = bwt.LF(range,c);

	int freq = 0;

	extract_consensus(bwt, R,ctx,freq,len-1);

	reverse(ctx.begin(), ctx.end());

	if(ctx.size()==len) left_contexts.push_back({ctx,freq});

}*/

//extract len characters using backward search, starting from position in range containing character c
void extract_consensus(dna_bwt_t & bwt, range_t range, char c, vector<pair<string, int> > & left_contexts, int len){

	string ctx;
	ctx += c;

	range_t R = bwt.LF(range,c);

	int freq = range_length(R);
	int covg = 0;

	extract_consensus(bwt, R,ctx,covg,len-1);

	reverse(ctx.begin(), ctx.end());

	if(ctx.size()==len) left_contexts.push_back({ctx,freq});

}

/*
 * extract the DNA string of length len starting from SA[i].
 * If # is found while extracting, extraction is interrupted (# is not returned in the result string).
 */
string extract_dna(dna_bwt_t & bwt, uint64_t i, int64_t len){

	string res;

	char c = bwt.F(i);

	while(c != TERM and len > 0){

		res += c;
		i = bwt.FL(i);
		c = bwt.F(i);
		len--;

	}

	return res;

}

/*
 * input: bwts and ranges of cluster (right-excluded)
 *
 * output: the variants found in the cluster between the two individuals
 */
vector<variant_t> find_variants(dna_bwt_t & bwt1, dna_bwt_t & bwt2, range_t range1, range_t range2){

	vector<variant_t>  out;

	auto counts = vector<vector<unsigned int> >(2,vector<unsigned int>(4,0));

	for(uint64_t i=range1.first;i<range1.second;++i) counts[0][base_to_int(bwt1[i])]++;
	for(uint64_t i=range2.first;i<range2.second;++i) counts[1][base_to_int(bwt2[i])]++;

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
	/*auto all_chars = frequent_char_0;
	all_chars.insert(all_chars.begin(), frequent_char_1.begin(), frequent_char_1.end());
	std::sort( all_chars.begin(), all_chars.end() );
	all_chars.erase(std::unique( all_chars.begin(), all_chars.end() ), all_chars.end());*/

	//filter: remove clusters that cannot reflect a variation
	if(	frequent_char_0.size()==0 or // not covered enough
		frequent_char_1.size()==0 or // not covered enough
		frequent_char_0.size() > (diploid ? 2 : 1) or 	// we require at most 2/1 alleles per individual (diploid/haploid)
		frequent_char_1.size() > (diploid ? 2 : 1) or 	// we require  at most 2/1 alleles per individual (diploid/haploid)
		frequent_char_0 == frequent_char_1   			// same alleles: probably both heterozigous / multiple region (and no variants)
		//all_chars.size() > 2 							// too many distinct frequent characters in the cluster (probably multiple region)
	){

		return out;

	}

	//for each frequent char in each of the two individuals, find the associated
	//left-context and corresponding coverage using backward search

	vector<pair<string, int> > left_contexts_0;
	vector<pair<string, int> > left_contexts_1;

	for(auto c : frequent_char_0) extract_consensus(bwt1, range1, c, left_contexts_0, k_left);
	for(auto c : frequent_char_1) extract_consensus(bwt2, range2, c, left_contexts_1, k_left);

	//find right-context

	string right_ctx = "";

	//first position in ranges and merged range (i)
	uint64_t i0 = range1.first;
	uint64_t i1 = range2.first;
	uint64_t i = i0 + i1;

	//find i such that LCP[i] >= k_right
	while(i<range1.second + range2.second and (not LCP_threshold[2*i+1])){

		i0 += DA[i]==0;
		i1 += DA[i]==1;
		++i;
	}

	//found good right context
	if(i<range1.second + range2.second and LCP_threshold[2*i+1]){

		right_ctx = extract_dna(DA[i] ? bwt2 : bwt1, DA[i] ? i1 : i0, k_right);
		assert(right_ctx.length() == k_right);

		//store all variants found
		for(auto L0 : left_contexts_0){

			for(auto L1 : left_contexts_1){

				if(L0.first.size()>0 and L1.first.size()>0){

					if(L0.first[L0.first.size()-1] != L1.first[L1.first.size()-1])
						out.push_back({L0.first,L1.first,right_ctx,L0.second, L1.second});

				}

			}

		}

	}

	return out;

}

/*
 * input: bwt and range of cluster (right-excluded)
 *
 * output: the variants found in the cluster
 */
vector<variant_single_t> find_variants(dna_bwt_t & bwt, range_t range){

	vector<variant_single_t>  out;

	auto counts = vector<unsigned int>(4,0);

	for(uint64_t i=range.first;i<range.second;++i) counts[base_to_int(bwt[i])]++;

	//compute the lists of frequent characters
	vector<unsigned char> frequent_char;

	for(int c=0;c<4;++c){

		if(counts[c] >= mcov_out) frequent_char.push_back(int_to_base(c));

	}

	//all variations observed in cluster

	//filter: remove clusters that cannot reflect a variation
	if(	frequent_char.size()<2 ) return out;

	//for each frequent char find the associated
	//left-context and corresponding coverage using backward search

	vector<pair<string, int> > left_contexts;

	for(auto c : frequent_char) extract_consensus(bwt, range, c, left_contexts, k_left);

	//find right-context

	string right_ctx = "";

	uint64_t i = range.first;

	//find i such that LCP[i] >= k_right
	while(i<range.second and (not LCP_threshold[2*i+1])) ++i;

	//found good right context
	if(i<range.second and LCP_threshold[2*i+1]){

		right_ctx = extract_dna(bwt, i, k_right);
		assert(right_ctx.length() == k_right);

		//store all variants found
		for(auto L : left_contexts){

			if(L.first.size()>0){

				out.push_back({L.first,right_ctx,L.second});

			}

		}

	}

	return out;

}


/*
 * input: bwt and range of cluster (right-excluded)
 *
 * output: the variants found in the cluster
 */
vector<variant_t> find_variants(dna_bwt_t & bwt, vector<bool> & DA, range_t range){

	vector<variant_t>  out;

	auto counts = vector<vector<unsigned int> >(2,vector<unsigned int>(4,0));

	for(uint64_t i=range.first;i<range.second;++i) counts[DA[i]][base_to_int(bwt[i])]++;

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
	/*auto all_chars = frequent_char_0;
	all_chars.insert(all_chars.begin(), frequent_char_1.begin(), frequent_char_1.end());
	std::sort( all_chars.begin(), all_chars.end() );
	all_chars.erase(std::unique( all_chars.begin(), all_chars.end() ), all_chars.end());*/

	//filter: remove clusters that cannot reflect a variation
	if(	frequent_char_0.size()==0 or // not covered enough
		frequent_char_1.size()==0 or // not covered enough
		frequent_char_0.size() > (diploid ? 2 : 1) or 	// we require at most 2/1 alleles per individual (diploid/haploid)
		frequent_char_1.size() > (diploid ? 2 : 1) or 	// we require  at most 2/1 alleles per individual (diploid/haploid)
		frequent_char_0 == frequent_char_1  			// same alleles: probably both heterozigous / multiple region (and no variants)
		//all_chars.size() > 2 							// too many distinct frequent characters in the cluster (probably multiple region)
	){

		return out;

	}

	vector<pair<string, int> > left_contexts_0;
	vector<pair<string, int> > left_contexts_1;

	for(auto c : frequent_char_0) extract_consensus(bwt, range, c, left_contexts_0, k_left);
	for(auto c : frequent_char_1) extract_consensus(bwt, range, c, left_contexts_1, k_left);

	//find right-context

	string right_ctx = "";

	//first position in ranges and merged range (i)
	uint64_t i = range.first;

	//find i such that LCP[i] >= k_right
	while(i<range.second and (not LCP_threshold[2*i+1])) ++i;

	//found good right context
	if(i<range.second and LCP_threshold[2*i+1]){

		right_ctx = extract_dna(bwt, i, k_right);
		assert(right_ctx.length() == k_right);

		//store all variants found
		for(auto L0 : left_contexts_0){

			for(auto L1 : left_contexts_1){

				if(L0.first.size()>0 and L1.first.size()>0){

					if(L0.first[L0.first.size()-1] != L1.first[L1.first.size()-1])
						out.push_back({L0.first,L1.first,right_ctx,L0.second, L1.second});

				}

			}

		}

	}

	return out;

}



/*
 * detect the type of variant (SNP/indel/discard if none) and, if not discarded, output to file the two reads per variant testifying it.
 */
void to_file(vector<variant_t> & output_variants, ofstream & out_file){

	bool found = false;

	id_nr = 1;//ID of the event inside the cluster. Here, events are pairs of reads; each read in the pair has the same ID

	for(auto v:output_variants){

		auto d = distance(v.left_context_0,v.left_context_1);

		if(		(not has_run(v.right_context,complexity)) and
				d.first <= max_snvs and
				v.support_0 >= mcov_out and v.support_1 >= mcov_out){

			found=true;

			//first individual

			string ID = ">";

			ID.append("cluster:");
			ID.append(to_string(cluster_nr));
			ID.append("_id:");
			ID.append(to_string(id_nr));
			ID.append("_right:");
			ID.append(to_string(v.right_context.size()));
			ID.append("_cov:");
			ID.append(std::to_string(v.support_0));//we write the number of reads supporting this variant
			ID.append("_type:");

			if(d.second != 0){

				ID.append("_INDEL_event:");

			}else{

				ID.append("_SNP_event:");

			}

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


			//second individual

			ID = ">";

			ID.append("cluster:");
			ID.append(to_string(cluster_nr));
			ID.append("_id:");
			ID.append(to_string(id_nr));
			ID.append("_right:");
			ID.append(to_string(v.right_context.size()));
			ID.append("_cov:");
			ID.append(std::to_string(v.support_1));//we write the number of reads supporting this variant
			ID.append("_type:");

			if(d.second != 0){

				ID.append("_INDEL_event:");

			}else{

				ID.append("_SNP_event:");

			}

			ID.append(snv_type);

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

	}

	cluster_nr += found;

}

void to_file(vector<variant_single_t> & output_variants, ofstream & out_file){

	if(output_variants.size()<2) return;

	int max_dist = 0;
	int nr_covered = 0;//number of covered-enough events in cluster

	for(int i=0;i<output_variants.size()-1;++i){

		auto d = distance(output_variants[i].left_context,output_variants[i+1].left_context);

		max_dist = d.first>max_dist?d.first:max_dist;

		nr_covered += output_variants[i].support >= mcov_out;

	}

	nr_covered += output_variants[output_variants.size()-1].support >= mcov_out;

	if(max_dist <= max_snvs and nr_covered >= 2){

		id_nr = 1;//ID of event inside cluster. Here, an event is a single read. Each read in the cluster has a different ID.

		for(auto v:output_variants){

			if((not has_run(v.right_context,complexity)) and v.support >= mcov_out){

				string ID = ">";

				ID.append("cluster:");
				ID.append(to_string(cluster_nr));
				ID.append("_id:");
				ID.append(to_string(id_nr));
				ID.append("_right:");
				ID.append(to_string(v.right_context.size()));
				ID.append("_cov:");
				ID.append(std::to_string(v.support));//we write the number of reads supporting this variant

				out_file << ID << endl;

				string DNA;
				DNA = v.left_context;
				DNA.append(v.right_context);

				out_file << DNA << endl;

				id_nr++;
				events++;

			}

		}

	}

	cluster_nr ++;

}


void update_DA(sa_leaf L1,sa_leaf L2, uint64_t & lcp_values, uint64_t & m){

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

		LCP_threshold[2*i] = (L1.depth >= K);
		LCP_threshold[2*i+1] = (L1.depth >= k_right);

		//LCP[i] = L1.depth;

		lcp_values++;

	}

}

void update_LCP_leaf(sa_leaf L, uint64_t & lcp_values){

	for(uint64_t i = L.rn.first+1; i<L.rn.second; ++i){

		LCP_threshold[2*i] = (L.depth >= K);
		LCP_threshold[2*i+1] = (L.depth >= k_right);

		lcp_values++;

	}

}

void update_DA(sa_leaf L1,sa_leaf L2, uint64_t & m){

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

void find_leaves(sa_node x1, sa_node x2, uint64_t & m){

	//find leaves that were ignored in the first pass
	if(range_length(child_TERM(x1))+range_length(child_TERM(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_TERM(x1),0};
		sa_leaf L2 = {child_TERM(x2),0};

		update_DA(L1,L2,m);

	}

	if(range_length(child_A(x1))+range_length(child_A(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_A(x1),0};
		sa_leaf L2 = {child_A(x2),0};

		update_DA(L1,L2,m);

	}

	if(range_length(child_C(x1))+range_length(child_C(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_C(x1),0};
		sa_leaf L2 = {child_C(x2),0};

		update_DA(L1,L2,m);

	}

	if(range_length(child_G(x1))+range_length(child_G(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_G(x1),0};
		sa_leaf L2 = {child_G(x2),0};

		update_DA(L1,L2,m);

	}

	if(range_length(child_T(x1))+range_length(child_T(x2)) == 1){

		//symbolic depth = 0. It will not be used in update_DA
		sa_leaf L1 = {child_T(x1),0};
		sa_leaf L2 = {child_T(x2),0};

		update_DA(L1,L2,m);

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

void update_lcp_minima(sa_node x, uint64_t & n_min){

/*
 * we have a minimum after the end of each child (that is different than #) of size at least 2 of the input node x, except
 * if the candidate minimum position is the last or exceeds the interval of x
 */

	if( x.first_C - x.first_A >= 2 and 	// there are at least 2 'A'
		x.first_C < x.last-1	 		// candidate min in x.first_C is not >= last position
		){

		LCP_minima[x.first_C] = true;
		n_min++;

	}

	if( x.first_G - x.first_C >= 2 and 	// there are at least 2 'C'
		x.first_G < x.last-1	 		// candidate min in x.first_G is not >= last position
		){

		LCP_minima[x.first_G] = true;
		n_min++;

	}

	if( x.first_T - x.first_G >= 2 and 	// there are at least 2 'G'
		x.first_T < x.last-1	 		// candidate min in x.first_T is not >= last position
		){

		LCP_minima[x.first_T] = true;
		n_min++;

	}

}

/*
 * input: two BWTs
 */
void run_two_datasets(){

	cout << "Phase 1/4: loading and indexing eBWTs ... " << flush;

	dna_bwt_t bwt1 = dna_bwt_t(input1,TERM);
	dna_bwt_t bwt2 = dna_bwt_t(input2,TERM);

	cout << "done." << endl;

	uint64_t n1 = bwt1.size();
	uint64_t n2 = bwt2.size();
	uint64_t n = n1+n2;

	cout << "\nPhase 2/4: merging eBWTs." << endl;

	DA = vector<bool>(n,false); //document array

	/*
	 * LCP_threshold[2*i] == 1 iff LCP[i] >= K
	 * LCP_threshold[2*i+1] == 1 iff LCP[i] >= k_right
	 */
	LCP_threshold = vector<bool>(2*n,false);
	//LCP = vector<uint8_t>(n,0);

	uint64_t da_values = 0;//number computed DA values
	uint64_t leaves = 0;//number of visited leaves
	uint64_t max_stack = 0;
	uint64_t lcp_values = 1;//number of computed LCP values

	{

		auto TMP_LEAVES = vector<pair<sa_leaf, sa_leaf> >(4);

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

			update_DA(L1, L2, lcp_values, da_values);

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

	LCP_minima = vector<bool>(n,false);

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
		find_leaves(N1, N2, da_values);

		//compute LCP values at the borders of merged's children
		update_lcp_threshold(merged, LCP_threshold, lcp_values, K, k_right);
		//update_lcp_threshold(merged, LCP_threshold, lcp_values, K, k_right, LCP);

		update_lcp_minima(merged, n_min);

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

	/*cout << "start checking LCP minima " << endl;
	for(uint64_t i = 1; i<n-1;++i){

		if(LCP_minima[i]){

			if(not (LCP[i-1] > LCP[i] and LCP[i+1] >= LCP[i])){
				cout << "Error: LCP_minima=1 but LCP = " << (uint32_t)LCP[i-1] << " " << (uint32_t)LCP[i] << " " <<  (uint32_t)LCP[i+1] << endl;
			}

		}else{

			if(LCP[i-1] > LCP[i] and LCP[i+1] >= LCP[i]){
				cout << "Error: LCP_minima=0 but LCP = " << (uint32_t)LCP[i-1] << " " << (uint32_t)LCP[i] << " " <<  (uint32_t)LCP[i+1] << endl;
			}

		}

	}
	cout << "Done checking LCP minima " << endl;*/

	cout << "Phase 4/4: detecting SNPs and indels." << endl;

	cout << "Output events will be stored in " << output << endl;
	ofstream out_file = ofstream(output);

	uint64_t begin0 = 0;//begin position on indiv. 0
	uint64_t begin1 = 0;//begin position on indiv. 1

	uint64_t i0 = 0;//current position on indiv. 0
	uint64_t i1 = 0;//current position on indiv. 1

	uint64_t clust_len=0;
	bool cluster_open=false;

	int perc = -1;
	int last_perc = -1;

	uint64_t n_clusters = 0; //number of clusters analyzed
	uint64_t clust_size = 0; //cumulative cluster size

	uint64_t MAX_CLUST_LEN = 200;
	auto CLUST_SIZES = vector<uint64_t>(MAX_CLUST_LEN+1,0);

	for(uint64_t i=0;i<n;++i){

		if(LCP_threshold[2*i] and not LCP_minima[i]){

			if(cluster_open){//extend current cluster
				clust_len++;
			}else{//open new cluster
				cluster_open=true;
				clust_len=1;
				begin0=i0;
				begin1=i1;
			}

		}else{

			if(cluster_open){//close current cluster

				clust_size += clust_len;

				if(clust_len<=MAX_CLUST_LEN) CLUST_SIZES[clust_len]+=clust_len;

				if(clust_len>=2*mcov_out){

					n_clusters++;
					vector<variant_t> var = find_variants(bwt1, bwt2, {begin0,i0}, {begin1,i1});//CLUSTER: two datasets (two BWTs)
					to_file(var,out_file);

				}

			}

			cluster_open=false;
			clust_len = 0;

		}

		i0 += (DA[i]==0);
		i1 += (DA[i]==1);

		perc = (100*i)/n;

		if(perc > last_perc){

			cout << perc << "%. ";
			cout << endl;

			last_perc = perc;

		}

	}

	cout 	<< endl << "Done." << endl <<
			"Analyzed " << n_clusters << " clusters." << endl <<
			"Average cluster length: " << double(clust_size)/n_clusters << "." << endl << endl <<
			"Distribution of bases inside clusters (cluster length / number of bases inside clusters of that length): " << endl << endl;



	uint64_t scale = *max_element(CLUST_SIZES.begin(), CLUST_SIZES.end());

	for(int i=0;i<=MAX_CLUST_LEN;++i){

		cout << i << ( i < 10 ? "   " : (i<100 ? "  " : " "));
		for(uint64_t j = 0; j < (100*CLUST_SIZES[i])/scale; ++j) cout << "-";
		cout << " " << CLUST_SIZES[i] << endl;

	}


}


/*
 * input: BWT and DA
 */
void run_two_datasets_da(){

	cout << "Phase 1/4: loading and indexing eBWT ... " << flush;

	dna_bwt_t bwt = dna_bwt_t(input1,TERM);

	cout << "done." << endl;

	uint64_t n = bwt.size();

	cout << "\nPhase 2/4: navigating suffix tree leaves." << endl;

	/*
	 * LCP_threshold[2*i] == 1 iff LCP[i] >= K
	 * LCP_threshold[2*i+1] == 1 iff LCP[i] >= k_right
	 */
	LCP_threshold = vector<bool>(2*n,false);

	uint64_t leaves = 0;//number of visited leaves
	uint64_t max_stack = 0;
	uint64_t lcp_values = 1;//number of computed LCP values

	{

		auto TMP_LEAVES = vector<sa_leaf>(4);

		stack<sa_leaf> S;
		S.push(bwt.first_leaf());

		int last_perc_lcp = -1;
		int perc_lcp = 0;

		while(not S.empty()){

			sa_leaf L = S.top();
			S.pop();
			leaves++;

			assert(leaf_size(L)>0);
			max_stack = S.size() > max_stack ? S.size() : max_stack;

			update_LCP_leaf(L,lcp_values);

			int t = 0;//number of children leaves
			bwt.next_leaves(L, TMP_LEAVES, t, 2);

			for(int i=t-1;i>=0;--i) S.push(TMP_LEAVES[i]);

			perc_lcp = (100*lcp_values)/n;

			if(perc_lcp > last_perc_lcp){

				cout << "LCP: " << perc_lcp << "%.";
				cout << endl;

				last_perc_lcp = perc_lcp;

			}

		}
	}

	cout << "Computed " << lcp_values << "/" << n << " LCP threshold values." << endl;

	cout << "Max stack depth = " << max_stack << endl;
	cout << "Processed " << leaves << " suffix-tree leaves." << endl << endl;

	cout << "Phase 3/4: computing LCP minima." << endl;

	LCP_minima = vector<bool>(n,false);

	auto TMP_NODES = vector<sa_node>(4);

	uint64_t nodes = 0;//visited ST nodes
	max_stack = 0;

	stack<sa_node> S;
	S.push(bwt.root());

	int last_perc_lcp = -1;
	int perc_lcp = 0;
	uint64_t n_min = 0;//number of LCP minima

	while(not S.empty()){

		max_stack = S.size() > max_stack ? S.size() : max_stack;

		sa_node N = S.top();
		S.pop();
		nodes++;

		//compute LCP values at the borders of N's children
		update_lcp_threshold(N, LCP_threshold, lcp_values, K, k_right);

		update_lcp_minima(N, n_min);

		//follow Weiner links
		int t = 0;
		bwt.next_nodes(N, TMP_NODES, t);

		for(int i=t-1;i>=0;--i) S.push(TMP_NODES[i]);

		perc_lcp = (100*lcp_values)/n;

		if(perc_lcp > last_perc_lcp){

			cout << "LCP: " << perc_lcp << "%.";
			cout << endl;

			last_perc_lcp = perc_lcp;

		}

	}

	cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;
	cout << "Found " << n_min << " LCP minima." << endl;
	cout << "Max stack depth = " << max_stack << endl;
	cout << "Processed " << nodes << " suffix-tree nodes." << endl << endl;

	cout << "Phase 3/4: detecting SNPs and indels." << endl;

	cout << "Output events will be stored in " << output << endl;
	ofstream out_file = ofstream(output);

	uint64_t begin = 0;//begin position

	uint64_t clust_len=0;
	bool cluster_open=false;

	int perc = -1;
	int last_perc = -1;

	uint64_t n_clusters = 0; //number of clusters analyzed
	uint64_t clust_size = 0; //cumulative cluster size

	uint64_t MAX_CLUST_LEN = 200;

	vector<bool> DA(n,false);

	ifstream da_file(input_da);
	char x;

	auto CLUST_SIZES = vector<uint64_t>(MAX_CLUST_LEN+1,0);

	//read DA
	for(uint64_t i=0;i<n;++i){

		da_file.read((char*)&x,1);
		DA[i] = (x=='1');

	}

	for(uint64_t i=0;i<n;++i){

		if(LCP_threshold[2*i] and not LCP_minima[i]){

			if(cluster_open){//extend current cluster

				clust_len++;

			}else{//open new cluster

				cluster_open=true;
				clust_len=1;
				begin=i;

			}

		}else{

			if(cluster_open){//close current cluster

				clust_size += clust_len;

				if(clust_len<=MAX_CLUST_LEN) CLUST_SIZES[clust_len]+=clust_len;

				if(clust_len>=2*mcov_out){

					n_clusters++;
					vector<variant_t> var = find_variants(bwt, DA, {begin,i});//CLUSTER: two datasets (one BWT and DA)
					to_file(var,out_file);

				}

			}

			cluster_open=false;
			clust_len = 0;

		}

		perc = (100*i)/n;

		if(perc > last_perc){

			cout << perc << "%. ";
			cout << endl;

			last_perc = perc;

		}

	}

	cout 	<< endl << "Done." << endl <<
			"Analyzed " << n_clusters << " clusters." << endl <<
			"Average cluster length: " << double(clust_size)/n_clusters << "." << endl << endl <<
			"Distribution of bases inside clusters (cluster length / number of bases inside clusters of that length): " << endl << endl;

	uint64_t scale = *max_element(CLUST_SIZES.begin(), CLUST_SIZES.end());

	for(int i=0;i<=MAX_CLUST_LEN;++i){

		cout << i << ( i < 10 ? "   " : (i<100 ? "  " : " "));
		for(uint64_t j = 0; j < (100*CLUST_SIZES[i])/scale; ++j) cout << "-";
		cout << " " << CLUST_SIZES[i] << endl;

	}

	cout << "\nStored to file " << events << " sequences clustered in " << (cluster_nr-1) << " clusters." << endl;

}

/*
 * input: BWT
 */
void run_one_dataset(){

	cout << "Phase 1/4: loading and indexing eBWT ... " << flush;

	dna_bwt_t bwt = dna_bwt_t(input1,TERM);

	cout << "done." << endl;

	uint64_t n = bwt.size();

	cout << "\nPhase 2/4: navigating suffix tree leaves." << endl;

	/*
	 * LCP_threshold[2*i] == 1 iff LCP[i] >= K
	 * LCP_threshold[2*i+1] == 1 iff LCP[i] >= k_right
	 */
	LCP_threshold = vector<bool>(2*n,false);

	uint64_t leaves = 0;//number of visited leaves
	uint64_t max_stack = 0;
	uint64_t lcp_values = 1;//number of computed LCP values

	{

		auto TMP_LEAVES = vector<sa_leaf>(4);

		stack<sa_leaf> S;
		S.push(bwt.first_leaf());

		int last_perc_lcp = -1;
		int perc_lcp = 0;

		while(not S.empty()){

			sa_leaf L = S.top();
			S.pop();
			leaves++;

			assert(leaf_size(L)>0);
			max_stack = S.size() > max_stack ? S.size() : max_stack;

			update_LCP_leaf(L,lcp_values);

			int t = 0;//number of children leaves
			bwt.next_leaves(L, TMP_LEAVES, t, 2);

			for(int i=t-1;i>=0;--i) S.push(TMP_LEAVES[i]);

			perc_lcp = (100*lcp_values)/n;

			if(perc_lcp > last_perc_lcp){

				cout << "LCP: " << perc_lcp << "%.";
				cout << endl;

				last_perc_lcp = perc_lcp;

			}

		}
	}

	cout << "Computed " << lcp_values << "/" << n << " LCP threshold values." << endl;

	cout << "Max stack depth = " << max_stack << endl;
	cout << "Processed " << leaves << " suffix-tree leaves." << endl << endl;

	cout << "Phase 3/4: computing LCP minima." << endl;

	LCP_minima = vector<bool>(n,false);

	auto TMP_NODES = vector<sa_node>(4);

	uint64_t nodes = 0;//visited ST nodes
	max_stack = 0;

	stack<sa_node> S;
	S.push(bwt.root());

	int last_perc_lcp = -1;
	int perc_lcp = 0;
	uint64_t n_min = 0;//number of LCP minima

	while(not S.empty()){

		max_stack = S.size() > max_stack ? S.size() : max_stack;

		sa_node N = S.top();
		S.pop();
		nodes++;

		//compute LCP values at the borders of N's children
		update_lcp_threshold(N, LCP_threshold, lcp_values, K, k_right);

		update_lcp_minima(N, n_min);

		//follow Weiner links
		int t = 0;
		bwt.next_nodes(N, TMP_NODES, t);

		for(int i=t-1;i>=0;--i) S.push(TMP_NODES[i]);

		perc_lcp = (100*lcp_values)/n;

		if(perc_lcp > last_perc_lcp){

			cout << "LCP: " << perc_lcp << "%.";
			cout << endl;

			last_perc_lcp = perc_lcp;

		}

	}

	cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;
	cout << "Found " << n_min << " LCP minima." << endl;
	cout << "Max stack depth = " << max_stack << endl;
	cout << "Processed " << nodes << " suffix-tree nodes." << endl << endl;

	cout << "Phase 4/4: detecting SNPs and indels." << endl;

	cout << "Output events will be stored in " << output << endl;
	ofstream out_file = ofstream(output);

	uint64_t begin = 0;//begin position

	uint64_t clust_len=0;
	bool cluster_open=false;

	int perc = -1;
	int last_perc = -1;

	uint64_t n_clusters = 0; //number of clusters analyzed
	uint64_t clust_size = 0; //cumulative cluster size

	uint64_t MAX_CLUST_LEN = 200;
	auto CLUST_SIZES = vector<uint64_t>(MAX_CLUST_LEN+1,0);

	for(uint64_t i=0;i<n;++i){

		if(LCP_threshold[2*i] and not LCP_minima[i]){

			if(cluster_open){//extend current cluster
				clust_len++;
			}else{//open new cluster
				cluster_open=true;
				clust_len=1;
				begin=i;
			}

		}else{

			if(cluster_open){//close current cluster

				clust_size += clust_len;

				if(clust_len<=MAX_CLUST_LEN) CLUST_SIZES[clust_len]+=clust_len;

				if(clust_len>=2*mcov_out){

					n_clusters++;
					vector<variant_single_t> var = find_variants(bwt, {begin,i});//CLUSTER: one dataset (one BWT)
					to_file(var,out_file);

				}

			}

			cluster_open=false;
			clust_len = 0;

		}

		perc = (100*i)/n;

		if(perc > last_perc){

			cout << perc << "%. ";
			cout << endl;

			last_perc = perc;

		}

	}

	cout 	<< endl << "Done." << endl <<
			"Analyzed " << n_clusters << " clusters." << endl <<
			"Average cluster length: " << double(clust_size)/n_clusters << "." << endl << endl <<
			"Stored to file " << events << " events clustered in " << (cluster_nr-1) << " clusters." << endl << endl <<
			"Distribution of bases inside clusters (cluster length / number of bases inside clusters of that length): " << endl;


	uint64_t scale = *max_element(CLUST_SIZES.begin(), CLUST_SIZES.end());

	for(int i=0;i<=MAX_CLUST_LEN;++i){

		cout << i << ( i < 10 ? "   " : (i<100 ? "  " : " "));
		for(uint64_t j = 0; j < (100*CLUST_SIZES[i])/scale; ++j) cout << "-";
		cout << " " << CLUST_SIZES[i] << endl;

	}

}


int main(int argc, char** argv){

	srand(time(NULL));

	if(argc < 3) help();

	int opt;
	while ((opt = getopt(argc, argv, "h1:2:v:L:R:m:g:k:t:o:Dd:c:")) != -1){
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
			case 'o':
				output = string(optarg);
			break;
			case '2':
				input2 = string(optarg);
			break;
			case 'd':
				input_da = string(optarg);
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
			case 'L':
				k_left = atoi(optarg);
			break;
			case 'R':
				k_right = atoi(optarg);
			break;
			case 'v':
				max_snvs = atoi(optarg);
			break;
			case 't':
				TERM = atoi(optarg);
			break;
			case 'c':
				complexity = atoi(optarg);
			break;
			case 'D':
				diploid = true;
			break;
			default:
				help();
			return -1;
		}
	}

	complexity = complexity == 0 ? complexity_def : complexity;
	K = K == 0 ? K_def : K;
	max_gap = max_gap==0?max_gap_def:max_gap;
	k_left = k_left==0?k_left_def:k_left;
	k_right = k_right==0?k_right_def:k_right;
	max_snvs = max_snvs==0?max_snvs_def:max_snvs;
	mcov_out = mcov_out==0?mcov_out_def:mcov_out;

	if(input1.compare("")==0 or output.compare("")==0) help();

	if(not file_exists(input1)){
		cout << "Error: could not find file " << input1 << endl << endl;
		help();
	}

	if(input2.length() > 0 and not file_exists(input2)){
		cout << "Error: could not find file " << input2 << endl << endl;
		help();
	}

	if(input2.compare("")!=0 and input_da.compare("")!=0){

		cout << "Error: Document array (-d) can only be used with one input BWT file (-1)" << endl << endl;
		help();

	}

	cout << "This is ebwt2snp, version 2." << endl;

	if(input2.length()>0){

		cout << "Running on two samples. Input eBWT files : " << input1 << " and " << input2 << endl;

	}else{

		if(input_da.compare("")!=0){

			cout << "Running on one sample with input Document array. Input eBWT/DA files : " << input1 << " and " << input_da << endl;

		}else{

			cout << "Running on one sample (genotyping). Input eBWT file : " << input1 << endl;

		}

	}

	cout << "Left-extending eBWT ranges by " << k_left << " bases." << endl <<
	"Right context length: " << k_right << " bases." << endl <<
	"Complexity filter: " << complexity << endl <<
	"Storing output events to file " << output << endl <<
	"Minimum coverage of output events: " << mcov_out << endl;

	if(input2.compare("")!=0 or input_da.compare("")!=0){
		if(diploid){

			cout << "Input: Diploid." << endl;

		}else{

			cout << "Input: Haploid." << endl;

		}
	}

	cout << endl;

	if(input2.length()>0){

		run_two_datasets();

	}else{

		if(input_da.length()>0){

			run_two_datasets_da();

		}else{

			run_one_dataset();

		}

	}

}
