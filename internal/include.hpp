// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#ifndef INCLUDE_HPP_
#define INCLUDE_HPP_

#include <fstream>
#include <vector>
#include <cassert>
#include <algorithm>

using namespace std;

typedef uint32_t int_text;
typedef uint32_t int_suff;
typedef uint32_t int_lcp;
typedef uint8_t int8;
typedef uint8_t dataTypelenSeq;	//length of the sequences (in biologic case 100)
typedef uint32_t dataTypeNSeq;	//number of sequences
typedef __int128 uint128_t;

typedef pair<uint64_t,uint32_t> coordinate;//suffix array coordinate (text, suff)

typedef pair<uint64_t,uint64_t> range_t;

/*
 * EGSA
 */
typedef struct{

	int_text	text; //read nr
	int_suff	suff; //starting position of the suffix in the read
	int_lcp 	lcp;
	int8		bwt;

} t_GSA;

/*
 * this class abstracts the EGSA type and allows reading from different formats (EGSA/BCR)
 */
class egsa_stream{

public:

	/*
	 * input: reads.fasta
	 *
	 * automatically detects the index files and format
	 *
	 */
	egsa_stream(string & input_path){

		string egsa_path = input_path;
		egsa_path.append(".gesa");

		EGSA.open(egsa_path, ios::in | ios::binary);

		if(EGSA.is_open()){

			egsa = true;

		}else{//else try BCR

			string LCP_path = input_path;
			LCP_path.append(".out.lcp");

			string BWT_path = input_path;
			BWT_path.append(".out");

			string GSA_path = input_path;
			GSA_path.append(".out.pairSA");

			LCP.open(LCP_path, ios::in | ios::binary);
			BWT.open(BWT_path, ios::in | ios::binary);
			GSA.open(GSA_path, ios::in | ios::binary);

			if(LCP.is_open() and BWT.is_open() and GSA.is_open()){

				bcr = true;

			}else{

				cout << "Error: missing index files." << endl;
				exit(1);

			}

		}

	}

	/*
	 * returns true iff index files exist in the input folder
	 */
	bool index_exists(){

		return egsa or bcr;

	}

	bool eof(){

		if(egsa){

			return EGSA.eof();

		}else if (bcr){

			return LCP.eof();

		}

		cout << "Error: missing index files." << endl;
		exit(1);

	}

	/*
	 * overwrite default type byte-sizes (LCP, document array, suffix in read)
	 */
	void set_bytesizes(int lcp_size, int da_size, int suff_size){

		this->lcp_size = lcp_size;
		this->da_size = da_size;
		this->suff_size = suff_size;

	}

	t_GSA read_el(){

		t_GSA e;

		if(egsa){

			switch(da_size){

				case 1 : uint8_t x8; EGSA.read((char*)&x8, 1); e.text = x8;  break;
				case 2 : uint16_t x16; EGSA.read((char*)&x16, 2); e.text = x16; break;
				case 4 : uint32_t x32; EGSA.read((char*)&x32, 4); e.text = x32; break;
				case 8 : uint64_t x64; EGSA.read((char*)&x64, 8); e.text = x64; break;

			}

			switch(suff_size){

				case 1 : uint8_t x8; EGSA.read((char*)&x8, 1); e.suff = x8;  break;
				case 2 : uint16_t x16; EGSA.read((char*)&x16, 2); e.suff = x16; break;
				case 4 : uint32_t x32; EGSA.read((char*)&x32, 4); e.suff = x32; break;
				case 8 : uint64_t x64; EGSA.read((char*)&x64, 8); e.suff = x64; break;

			}

			switch(lcp_size){

				case 1 : uint8_t x8; EGSA.read((char*)&x8, 1); e.lcp = x8;  break;
				case 2 : uint16_t x16; EGSA.read((char*)&x16, 2); e.lcp = x16; break;
				case 4 : uint32_t x32; EGSA.read((char*)&x32, 4); e.lcp = x32; break;
				case 8 : uint64_t x64; EGSA.read((char*)&x64, 8); e.lcp = x64; break;

			}

			uint8_t x8;
			EGSA.read((char*)&x8, 1);
			e.bwt = x8;

		}else if(bcr){

			switch(suff_size){

				case 1 : uint8_t x8; GSA.read((char*)&x8, 1); e.suff = x8;  break;
				case 2 : uint16_t x16; GSA.read((char*)&x16, 2); e.suff = x16; break;
				case 4 : uint32_t x32; GSA.read((char*)&x32, 4); e.suff = x32; break;
				case 8 : uint64_t x64; GSA.read((char*)&x64, 8); e.suff = x64; break;

			}

			switch(da_size){

				case 1 : uint8_t x8; GSA.read((char*)&x8, 1); e.text = x8;  break;
				case 2 : uint16_t x16; GSA.read((char*)&x16, 2); e.text = x16; break;
				case 4 : uint32_t x32; GSA.read((char*)&x32, 4); e.text = x32; break;
				case 8 : uint64_t x64; GSA.read((char*)&x64, 8); e.text = x64; break;

			}

			switch(lcp_size){

				case 1 : uint8_t x8; LCP.read((char*)&x8, 1); e.lcp = x8;  break;
				case 2 : uint16_t x16; LCP.read((char*)&x16, 2); e.lcp = x16; break;
				case 4 : uint32_t x32; LCP.read((char*)&x32, 4); e.lcp = x32; break;
				case 8 : uint64_t x64; LCP.read((char*)&x64, 8); e.lcp = x64; break;

			}

			uint8_t x8;
			BWT.read((char*)&x8, 1);
			e.bwt = x8;

		}else{

			cout << "Error: missing index files." << endl;
			exit(1);

		}

		return e;

	}

private:

	bool egsa = false;
	bool bcr = false;

	//byte size of components
	int lcp_size = 1; //LCP values
	int da_size = 4; //document array (read number)
	int suff_size = 1; //position inside read

	//the EGSA index
	ifstream EGSA;

	//the BCR index
	ifstream LCP;
	ifstream BWT;
	ifstream GSA;//pairs

};

t_GSA read_el(ifstream & egsa, bool bcr){

	t_GSA e;

	if(bcr){

        dataTypeNSeq txt;
        uint8_t suf;
        uint8_t lcp;

        egsa.read((char*)&txt, sizeof(dataTypeNSeq));
        egsa.read((char*)&suf, sizeof(dataTypelenSeq));
        egsa.read((char*)&lcp, sizeof(dataTypelenSeq));
        egsa.read((char*)&e.bwt, sizeof(uint8_t));

        e.suff = suf;
        e.text = txt;
        e.lcp = lcp;

	}else{

		egsa.read((char*)&e, sizeof(int_text)+sizeof(int_suff)+sizeof(int_lcp)+sizeof(int8));

	}

	return e;

}

unsigned char int_to_base(int i){

	switch(i){

		case 0: return 'A'; break;
		case 1: return 'C'; break;
		case 2: return 'G'; break;
		case 3: return 'T'; break;

	}

	return 'A';

}

int base_to_int(unsigned char c){

	switch(c){

		case 'A': case 'a': return 0; break;
		case 'C': case 'c': return 1; break;
		case 'G': case 'g': return 2; break;
		case 'T': case 't': return 3; break;
		case 'N': case 'n': return rand()%4; break;

	}

	return 0;

}

unsigned char RC(unsigned char c){

	switch(c){

		case 'A': case 'a': return 'T'; break;
		case 'C': case 'c': return 'G'; break;
		case 'G': case 'g': return 'C'; break;
		case 'T': case 't': return 'A'; break;
		default: break;

	}

	return 'N';

}

string RC(string & s){

	if(s.length()==0) return string();

	string rc = s;

	for(int i=0;i<s.length();++i) rc[rc.length()-i-1] = RC(s[i]);

	return rc;

}


string rev(string & a){

	string reversed(a.rbegin(), a.rend());

	return reversed;

}


inline uint16_t clz_u128 (uint128_t u) {

	uint64_t hi = u>>64;
	uint64_t lo = u;

	return hi!= 0 ? __builtin_clzll(hi) : 64 + __builtin_clzll(lo);

}

inline int popcount128(uint128_t u){

	return __builtin_popcountl(uint64_t(u)) + __builtin_popcountl(uint64_t(u>>64));

}

class cons{

public:

	cons(int size){

		counts = vector<vector<int>>(size,vector<int>(4,0));
		C = string(size,'A');

	}

	unsigned char operator[](int i){
		return C[i];
	}

	void increment(int i, unsigned char b){

		int b_i = base_to_int(b);

		counts[i][b_i]++;

		if( counts[i][b_i] > counts[i][ base_to_int(C[i]) ])
			C[i] = b;

	}

	string to_string(){

		return C;

	}

private:

	string C;//the current consensus
	vector<vector<int>> counts;//base counts

};


//const char TERM = '#';

std::ifstream::pos_type filesize(string filename){
    std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

/*
 * representation of a right-maximal substring (SA node) as a list of BWT intervals
 */
struct sa_node{

	//right-maximal substring: string W such that Wa_1, ..., Wa_k occur in the text for
	//at least k>=2 characters a_1, ..., a_k

	uint64_t first_TERM;
	uint64_t first_A;
	uint64_t first_C;
	uint64_t first_G;
	uint64_t first_T;
	uint64_t last;

	//depth = |W|
	uint64_t depth;

	uint64_t key(){
		return first_TERM;
	}

};

struct sa_node_n{

	//right-maximal substring: string W such that Wa_1, ..., Wa_k occur in the text for
	//at least k>=2 characters a_1, ..., a_k

	uint64_t first_TERM;
	uint64_t first_A;
	uint64_t first_C;
	uint64_t first_G;
	uint64_t first_N;
	uint64_t first_T;
	uint64_t last;

	//depth = |W|
	uint64_t depth;

	uint64_t key(){
		return first_TERM;
	}

};

uint64_t node_size(sa_node s){
	return s.last - s.first_TERM;
}

uint64_t node_size(pair<sa_node, sa_node> p){
	return node_size(p.first) + node_size(p.second);
}

uint64_t node_size(sa_node_n s){
	return s.last - s.first_TERM;
}

uint64_t node_size(pair<sa_node_n, sa_node_n> p){
	return node_size(p.first) + node_size(p.second);
}

void print_node(sa_node n){

	cout << "[" << 	n.first_TERM << ", " <<
					n.first_A << ", " <<
					n.first_C << ", " <<
					n.first_G << ", " <<
					n.first_T << ", " <<
					n.last << "]" << endl;

}

void print_node(sa_node_n n){

	cout << "[" << 	n.first_TERM << ", " <<
					n.first_A << ", " <<
					n.first_C << ", " <<
					n.first_G << ", " <<
					n.first_N << ", " <<
					n.first_T << ", " <<
					n.last << "]" << endl;

}

sa_node merge_nodes(sa_node a, sa_node b){

	assert(a.depth == b.depth);

	return {
		a.first_TERM + b.first_TERM,
		a.first_A + b.first_A,
		a.first_C + b.first_C,
		a.first_G + b.first_G,
		a.first_T + b.first_T,
		a.last + b.last,
		a.depth
	};

}

sa_node_n merge_nodes(sa_node_n a, sa_node_n b){

	assert(a.depth == b.depth);

	return {
		a.first_TERM + b.first_TERM,
		a.first_A + b.first_A,
		a.first_C + b.first_C,
		a.first_G + b.first_G,
		a.first_N + b.first_N,
		a.first_T + b.first_T,
		a.last + b.last,
		a.depth
	};

}

/*
 * suffix array leaf = BWT range (inclusive) of W.TERM, for some string W.
 *
 */
struct sa_leaf{

	//rn.first = first position of range. Equivalently, number of suffixes smaller than W.TERM (valid also if W.TERM does not occur)
	//rn.second = last position (excluded) of interval.  Equivalently, number of suffixes smaller than W.TERM + number of occurrences of W.TERM
	//if last == first, then W.TERM does not occur (however, 'first' is in any case number of suffixes smaller than W.TERM)
	range_t rn;

	//depth = |W.TERM|
	uint64_t depth;

	uint64_t key(){
		return rn.first;
	}

};

inline uint64_t range_length(range_t r){
	assert(r.second >= r.first);
	return r.second - r.first;
}

inline uint64_t leaf_size(sa_leaf L){
	return range_length(L.rn);
}

inline uint64_t leaf_size(pair<sa_leaf, sa_leaf> P){
	return leaf_size(P.first) + leaf_size(P.second);
}


struct p_range{

	range_t A;
	range_t C;
	range_t G;
	range_t T;

};

struct p_node{

	sa_node A;
	sa_node C;
	sa_node G;
	sa_node T;

};

struct p_range_n{

	range_t A;
	range_t C;
	range_t G;
	range_t N;
	range_t T;

};

struct p_node_n{

	sa_node_n A;
	sa_node_n C;
	sa_node_n G;
	sa_node_n N;
	sa_node_n T;

};

void print_nodes(p_node p){

	print_node(p.A);
	print_node(p.C);
	print_node(p.G);
	print_node(p.T);

}

struct p_rank{

public:

	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t T;

	p_rank operator+(const p_rank& a) const{

		return {
			a.A + A,
			a.C + C,
			a.G + G,
			a.T + T
		};

	}

	bool operator==(const p_rank& a) const{

		return a.A == A and a.C == C and a.G == G and a.T == T;

	}

	bool operator!=(const p_rank& a) const{

		return a.A != A or a.C != C or a.G != G or a.T != T;

	}

	bool operator<=(const p_rank& a) const{

		return A <= a.A and C <= a.C and G <= a.G and T <= a.T;

	}

};

struct p_rank_n{

public:

	uint64_t A;
	uint64_t C;
	uint64_t G;
	uint64_t N;
	uint64_t T;

	p_rank_n operator+(const p_rank_n& a) const{

		return {
			a.A + A,
			a.C + C,
			a.G + G,
			a.N + N,
			a.T + T
		};

	}

	bool operator==(const p_rank_n& a) const{

		return a.A == A and a.C == C and a.G == G and a.N == N and a.T == T;

	}

	bool operator!=(const p_rank_n& a) const{

		return a.A != A or a.C != C or a.G != G or a.N != N or a.T != T;

	}

	bool operator<=(const p_rank_n& a) const{

		return A <= a.A and C <= a.C and G <= a.G and N <= a.N and T <= a.T;

	}

};

inline p_range fold_ranks(p_rank &a, p_rank &b){

	return {{a.A, b.A},{a.C, b.C},{a.G, b.G},{a.T, b.T}};

}

inline p_range_n fold_ranks(p_rank_n &a, p_rank_n &b){

	return {{a.A, b.A},{a.C, b.C},{a.G, b.G},{a.N, b.N},{a.T, b.T}};

}

inline uint64_t popcount128(__uint128_t x){

	return __builtin_popcountll(uint64_t(x>>64)) + __builtin_popcountll( x & 0xFFFFFFFFFFFFFFFF );

}

range_t child_TERM(sa_node x){
	return {x.first_TERM, x.first_A};
}
range_t child_A(sa_node x){
	return {x.first_A, x.first_C};
}
range_t child_C(sa_node x){
	return {x.first_C, x.first_G};
}
range_t child_G(sa_node x){
	return {x.first_G, x.first_T};
}
range_t child_T(sa_node x){
	return {x.first_T, x.last};
}

range_t child_TERM(sa_node_n x){
	return {x.first_TERM, x.first_A};
}
range_t child_A(sa_node_n x){
	return {x.first_A, x.first_C};
}
range_t child_C(sa_node_n x){
	return {x.first_C, x.first_G};
}
range_t child_G(sa_node_n x){
	return {x.first_G, x.first_N};
}
range_t child_N(sa_node_n x){
	return {x.first_N, x.first_T};
}
range_t child_T(sa_node_n x){
	return {x.first_T, x.last};
}

inline bool has_child_TERM(sa_node N){
	return N.first_A > N.first_TERM;
}
inline bool has_child_A(sa_node N){
	return N.first_C > N.first_A;
}
inline bool has_child_C(sa_node N){
	return N.first_G > N.first_C;
}
inline bool has_child_G(sa_node N){
	return N.first_T > N.first_G;
}
inline bool has_child_T(sa_node N){
	return N.last > N.first_T;
}

inline bool has_child_TERM(sa_node_n N){
	return N.first_A > N.first_TERM;
}
inline bool has_child_A(sa_node_n N){
	return N.first_C > N.first_A;
}
inline bool has_child_C(sa_node_n N){
	return N.first_G > N.first_C;
}
inline bool has_child_G(sa_node_n N){
	return N.first_N > N.first_G;
}
inline bool has_child_N(sa_node_n N){
	return N.first_T > N.first_N;
}
inline bool has_child_T(sa_node_n N){
	return N.last > N.first_T;
}

uint8_t number_of_children(sa_node N){

	return 	uint8_t(N.last>N.first_T) +
			uint8_t(N.first_T>N.first_G) +
			uint8_t(N.first_G>N.first_C) +
			uint8_t(N.first_C>N.first_A) +
			uint8_t(N.first_A>N.first_TERM);

}

uint8_t number_of_children(sa_node_n N){

	return 	uint8_t(N.last>N.first_T) +
			uint8_t(N.first_T>N.first_N) +
			uint8_t(N.first_N>N.first_G) +
			uint8_t(N.first_G>N.first_C) +
			uint8_t(N.first_C>N.first_A) +
			uint8_t(N.first_A>N.first_TERM);

}

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(sa_node N1, sa_node N2){

	return 	uint8_t((N1.last>N1.first_T) or (N2.last>N2.first_T)) +
			uint8_t((N1.first_T>N1.first_G) or (N2.first_T>N2.first_G)) +
			uint8_t((N1.first_G>N1.first_C) or (N2.first_G>N2.first_C)) +
			uint8_t((N1.first_C>N1.first_A) or (N2.first_C>N2.first_A)) +
			uint8_t((N1.first_A>N1.first_TERM) or (N2.first_A>N2.first_TERM));

}

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(sa_node_n N1, sa_node_n N2){

	return 	uint8_t((N1.last>N1.first_T) or (N2.last>N2.first_T)) +
			uint8_t((N1.first_T>N1.first_N) or (N2.first_T>N2.first_N)) +
			uint8_t((N1.first_N>N1.first_G) or (N2.first_N>N2.first_G)) +
			uint8_t((N1.first_G>N1.first_C) or (N2.first_G>N2.first_C)) +
			uint8_t((N1.first_C>N1.first_A) or (N2.first_C>N2.first_A)) +
			uint8_t((N1.first_A>N1.first_TERM) or (N2.first_A>N2.first_TERM));

}

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(pair<sa_node,sa_node> P){

	return number_of_children(P.first, P.second);

}

/*
 * number of children in the union of the two nodes
 */
uint8_t number_of_children(pair<sa_node_n,sa_node_n> P){

	return number_of_children(P.first, P.second);

}

void update_lcp_threshold(sa_node & x, vector<bool> & LCP_threshold, uint64_t & lcp_values, int K, int k_right){

	assert(x.first_A >= x.first_TERM);
	assert(x.first_C >= x.first_A);
	assert(x.first_G >= x.first_C);
	assert(x.first_T >= x.first_G);

	assert(number_of_children(x) >= 2);

	if(has_child_TERM(x) and x.first_A != x.last){
		LCP_threshold[2*x.first_A] = (x.depth>=K);
		LCP_threshold[2*x.first_A+1] = (x.depth>=k_right);
		//LCP[x.first_A] = x.depth;
		lcp_values++;
	}
	if(has_child_A(x) and x.first_C != x.last){
		LCP_threshold[2*x.first_C] = (x.depth>=K);
		LCP_threshold[2*x.first_C+1] = (x.depth>=k_right);
		//LCP[x.first_C] = x.depth;
		lcp_values++;
	}
	if(has_child_C(x) and x.first_G != x.last){
		LCP_threshold[2*x.first_G] = (x.depth>=K);
		LCP_threshold[2*x.first_G+1] = (x.depth>=k_right);
		//LCP[x.first_G] = x.depth;
		lcp_values++;
	}
	if(has_child_G(x) and x.first_T != x.last){
		LCP_threshold[2*x.first_T] = (x.depth>=K);
		LCP_threshold[2*x.first_T+1] = (x.depth>=k_right);
		//LCP[x.first_T] = x.depth;
		lcp_values++;
	}

}

#endif /* INCLUDE_HPP_ */

