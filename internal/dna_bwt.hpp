// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * dna_bwt.hpp
 *
 *  Created on: Dec 3, 2018
 *      Author: nico
 *
 *  BWT for big files: built in chunks. The class is a template on the string type (could be run-length, Huffman...)
 *
 */

#include "include.hpp"
#include "dna_string.hpp"
#include <algorithm>

using namespace std;

#ifndef INTERNAL_DNA_BWT_HPP_
#define INTERNAL_DNA_BWT_HPP_

template<class str_type>
class dna_bwt{

public:

	using sa_node_t = sa_node;

	dna_bwt(){};

	/*
	 * constructor path of a BWT file containing the BWT in ASCII format
	 */
	dna_bwt(string path, char TERM = '#'){

		this->TERM = TERM;

		n = uint64_t(filesize(path));

		BWT = str_type(path, TERM);

		uint64_t n_term = 0;

		//build F column
		for(uint64_t i=0;i<n;++i){

			assert(BWT[i] < 256);

			F_A += (BWT[i]==TERM);
			F_C += (BWT[i]=='A');
			F_G += (BWT[i]=='C');
			F_T += (BWT[i]=='G');

		}

		F_C += F_A;
		F_G += F_C;
		F_T += F_G;

	}

	/*
	 * get full BWT range
	 */
	range_t full_range(){

		//right-exclusive range
		return {0,size()};

	}

	/*
	 * left-extend range by all 4 nucleotides
	 */
	p_range LF(range_t rn){

		assert(rn.second >= rn.first);

		//number of A,C,T,G before start of interval
		p_rank start = BWT.parallel_rank(rn.first);

		//number of A,C,T,G before end of interval (last position of interval included)
		p_rank end;

		if(rn.second>rn.first)
			end	= BWT.parallel_rank(rn.second);
		else
			end = start;

		assert(start <= end);

		p_rank f = {F_A,F_C,F_G,F_T};
		p_rank l = f + start;
		p_rank r = f + end;

		assert(r.A <= l.C);
		assert(r.C <= l.G);
		assert(r.G <= l.T);
		assert(l.T <= n);

		return fold_ranks(l,r);

	}

	char operator[](uint64_t i){

		return BWT[i];

	}

	/*
	 * number of c before position i excluded
	 */
	uint64_t rank(uint64_t i, uint8_t c){

		assert(i<=n);

		return BWT.rank(i,c);

	}

	/*
	 * return number of occurrences of A,C,T,G in the prefix of length i of the text. At most 1 cache miss!
	 */
	p_rank parallel_rank(uint64_t i){

		return BWT.parallel_rank(i);

	}

	uint64_t size(){

		assert(n == BWT.size());
		return n;

	}

	uint64_t serialize(std::ostream& out){

		uint64_t w_bytes = 0;

		out.write((char*)&n,sizeof(n));
		out.write((char*)&F_A,sizeof(uint64_t));
		out.write((char*)&F_C,sizeof(uint64_t));
		out.write((char*)&F_G,sizeof(uint64_t));
		out.write((char*)&F_T,sizeof(uint64_t));

		w_bytes += sizeof(n) + sizeof(uint64_t)*4;

		w_bytes += BWT.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&n,sizeof(n));

		in.read((char*)F_A,sizeof(uint64_t));
		in.read((char*)F_C,sizeof(uint64_t));
		in.read((char*)F_G,sizeof(uint64_t));
		in.read((char*)F_T,sizeof(uint64_t));

		BWT.load(in);

	}

	void save_to_file(string path){

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	/*
	 * path = path of an index file
	 */
	void load_from_file(string path){

		std::ifstream in(path);
		load(in);
		in.close();

	}


	/*
	 * functions for suffix tree navigation
	 */

	sa_node root(){

		return {
			0,
			F_A,
			F_C,
			F_G,
			F_T,
			n,
			0
		};

	}

	/*
	 * depth = LCP inside the leaf.
	 */
	sa_leaf first_leaf(){

		return {{0, F_A}, 0};

	}

	/*
	 * Input: suffix tree node N.
	 * Output: 4 suffix tree nodes (explicit, implicit, or empty) reached applying LF for A,C,G,T from N
	 */
	p_node LF(sa_node N){

		p_rank before_TERM;
		p_rank before_A;
		p_rank before_C;
		p_rank before_G;
		p_rank before_T;
		p_rank before_end;

		before_TERM = parallel_rank(N.first_TERM);

		if(N.first_A == N.first_TERM) before_A = before_TERM;
		else before_A = parallel_rank(N.first_A);

		if(N.first_C == N.first_A) before_C = before_A;
		else before_C = parallel_rank(N.first_C);

		if(N.first_G == N.first_C) before_G = before_C;
		else before_G = parallel_rank(N.first_G);

		if(N.first_T == N.first_G) before_T = before_G;
		else before_T = parallel_rank(N.first_T);

		if(N.last == N.first_T) before_end = before_T;
		else before_end = parallel_rank(N.last);

		return {
			{F_A + before_TERM.A, F_A + before_A.A, F_A + before_C.A, F_A + before_G.A, F_A + before_T.A, F_A + before_end.A, N.depth+1},
			{F_C + before_TERM.C, F_C + before_A.C, F_C + before_C.C, F_C + before_G.C, F_C + before_T.C, F_C + before_end.C, N.depth+1},
			{F_G + before_TERM.G, F_G + before_A.G, F_G + before_C.G, F_G + before_G.G, F_G + before_T.G, F_G + before_end.G, N.depth+1},
			{F_T + before_TERM.T, F_T + before_A.T, F_T + before_C.T, F_T + before_G.T, F_T + before_T.T, F_T + before_end.T, N.depth+1}
		};

	}

	void next_leaves(sa_leaf & L, vector<sa_leaf> & TMP_LEAVES, int & t, int min_n_children){

		t = 0;

		p_range ext = LF(L.rn);

		sa_leaf A = {ext.A, L.depth+1};
		sa_leaf C = {ext.C, L.depth+1};
		sa_leaf G = {ext.G, L.depth+1};
		sa_leaf T = {ext.T, L.depth+1};

		if(leaf_size(A)>=min_n_children) TMP_LEAVES[t++] = A;
		if(leaf_size(C)>=min_n_children) TMP_LEAVES[t++] = C;
		if(leaf_size(G)>=min_n_children) TMP_LEAVES[t++] = G;
		if(leaf_size(T)>=min_n_children) TMP_LEAVES[t++] = T;

		std::sort( TMP_LEAVES.begin(), TMP_LEAVES.begin()+t, [ ]( const sa_leaf& lhs, const sa_leaf& rhs )
		{
			return leaf_size(lhs) < leaf_size(rhs);
		});

	}

	void next_nodes(sa_node & N, vector<sa_node> & TMP_NODES, int & t){

		p_node left_exts = LF(N);

		sa_node A = left_exts.A;
		sa_node C = left_exts.C;
		sa_node G = left_exts.G;
		sa_node T = left_exts.T;

		t = 0;

		if(number_of_children(A) >= 2) TMP_NODES[t++] = A;
		if(number_of_children(C) >= 2) TMP_NODES[t++] = C;
		if(number_of_children(G) >= 2) TMP_NODES[t++] = G;
		if(number_of_children(T) >= 2) TMP_NODES[t++] = T;

		//push right-maximal nodes on stack in decreasing size (i.e. interval length) order

		std::sort( TMP_NODES.begin(), TMP_NODES.begin()+t, [ ]( const sa_node& lhs, const sa_node& rhs )
		{
			return node_size(lhs) < node_size(rhs);
		});

	}

private:

	char TERM = '#';

	uint64_t n = 0;//BWT length

	uint64_t F_A=0; //F array
	uint64_t F_C=0; //F array
	uint64_t F_G=0; //F array
	uint64_t F_T=0; //F array

	//vector<uint64_t> F;
	str_type BWT;

};

typedef dna_bwt<dna_string> dna_bwt_t;

#endif /* INTERNAL_DNA_BWT_HPP_ */
