// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * dna_string.hpp
 *
 *  Created on: Dec 6, 2018
 *      Author: nico
 *
 *  Optimized string with rank on DNA alphabet: {A,C,G,T,TERM}, where the terminator TERM is defined in include.hpp
 *
 *  One access or a parallel rank for all 4 letters A,C,G,T causes only 1 cache miss in the worst case
 *
 *  Max string length: 2^64
 *
 *  Supports very efficient (1 cache miss) parallel rank for (A,C,G,T), and (1 cache miss) single rank for TERM
 *
 *  Data is stored and cache-aligned in blocks of 512 bits (64 bytes)
 *
 *  Size of the string: 4n bits, where n = string length
 *
 *  memory layout:
 *
 *  | 32-bit rank A | 32-bit rank C | 32-bit rank G | 32-bit rank T|
 *  | 128-bit 1st bits | 128-bit 2nd bits | 128-bit 3rd bits |
 *
 *  careful: __uint128_t inverts the first and last 64 bits of the number!
 *
 *
 */

#ifndef INTERNAL_DNA_STRING_HPP_
#define INTERNAL_DNA_STRING_HPP_

#define SUPERBLOCK_SIZE 0x100000000 	//number of characters in a superblock = 2^32 characters
#define BLOCKS_PER_SUPERBLOCK 33554432	//blocks in a superblock
#define BYTES_PER_SUPERBLOCK 2147483648	//bytes in a superblock
#define BLOCK_SIZE 128 					//number of characters inside a block
#define BYTES_PER_BLOCK 64				//bytes in a block of 512 bits
#define ALN 64							//alignment

#include "include.hpp"
#include <cassert>

class dna_string{

public:

	dna_string(){}

	/*
	 * constructor from ASCII file
	 */
	dna_string(string path, char TERM = '#'){

		this->TERM = TERM;

		n = uint64_t(filesize(path));

		n_superblocks = (n+1)/SUPERBLOCK_SIZE + ((n+1)%SUPERBLOCK_SIZE != 0);
		n_blocks = (n+1)/BLOCK_SIZE + ((n+1)%BLOCK_SIZE != 0);
		nbytes = (n_blocks * BYTES_PER_BLOCK);//number of bytes effectively filled with data

		superblock_ranks = vector<p_rank>(n_superblocks);

		/*
		 * this block of code ensures that data is aligned by 64 bytes = 512 bits
		 */
		memory = vector<uint8_t>(nbytes+ALN,0);
		data = memory.data();
		while(uint64_t(data) % ALN != 0) data++;

		//cout << "alignment of data: " << (void*)data << endl;

		{

			ifstream ifs(path);

			string BUF(BLOCK_SIZE,'A');

			for(uint64_t i = 0; i<n; ++i){

				char c;

				ifs.read((char*)&c, sizeof(char));

				BUF[i%BLOCK_SIZE] = c;

				if(c!='A' and c!='C' and c!='G' and c!='T' and c!=TERM){
					cout << "Error while reading file: read forbidden character '" <<  c << "' (ASCII code " << int(c) << ")." << endl <<
					"Only A,C,G,T, and " << TERM << " are admitted in the input BWT!" << endl <<
					"If the unknown character is the terminator, you can solve the problem by adding option \"-t " << int(c) << "\"." << endl;

					exit(1);
				}

				//buffer is full
				if((i%BLOCK_SIZE) == BLOCK_SIZE-1) set(i/BLOCK_SIZE, BUF);

			}

			if(n % BLOCK_SIZE != 0) set(n/BLOCK_SIZE, BUF);

		}

		assert(check_content(path));
		build_rank_support();

	}

	//return i-th character
	char operator[](uint64_t i){

		assert(i<n);

		uint64_t superblock_number = i / SUPERBLOCK_SIZE;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE;
		uint64_t block_number = superblock_off / BLOCK_SIZE;
		uint64_t block_off = superblock_off % BLOCK_SIZE;

		//chars[0..3] contains 1st, 2nd, 3rd least significant bits of the BLOCK_SIZE characters
		__uint128_t* chars = (__uint128_t*)(data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK);

		uint64_t b =	((chars[0]>>(BLOCK_SIZE-(block_off+1)))&0x1) +
						(((chars[1]>>(BLOCK_SIZE-(block_off+1)))&0x1)<<1) +
						(((chars[2]>>(BLOCK_SIZE-(block_off+1)))&0x1)<<2);

		return 	(b == 0)*'A' +
				(b == 1)*'C' +
				(b == 2)*'G' +
				(b == 3)*'T' +
				(b == 4)*TERM;

	}

	/*
	 * Parallel rank of (A,C,T,G) at position i.
	 */
	p_rank parallel_rank(uint64_t i){

		uint64_t superblock_number = i / SUPERBLOCK_SIZE;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE;
		uint64_t block_number = superblock_off / BLOCK_SIZE;
		uint64_t block_off = superblock_off % BLOCK_SIZE;

		p_rank superblock_r = superblock_ranks[superblock_number];
		p_rank block_r = get_counters(superblock_number,block_number);

		return superblock_r + block_r + block_rank(superblock_number, block_number, block_off);

	}

	/*
	 * standard rank. c can be A,C,G,T, or TERM
	 */
	uint64_t rank(uint64_t i, uint8_t c){

		assert(i<=n);

		p_rank pr = parallel_rank(i);

		if(c==TERM) return rank_non_dna(i);

		switch(c){
			case 'A' : return pr.A; break;
			case 'C' : return pr.C; break;
			case 'G' : return pr.G; break;
			case 'T' : return pr.T; break;
		}

		return 0;

	}

	/*
	 * select. Implemented with binary search + rank (much more inefficient than rank!)
	 *
	 * indexes start from 0
	 *
	 */
	uint64_t select(uint64_t i, uint8_t c){

		assert(i<rank(n,c));

		return select(i,0,n,c);

	}


	/*
	 * return number of non-dna symbols in the prefix of length i of the text. At most 1 cache miss!
	 */
	uint64_t rank_non_dna(uint64_t i){

		assert(i<=n);
		auto r = parallel_rank(i);

		assert(r.A + r.C + r.G + r.T <= i);

		return i - (r.A + r.C + r.G + r.T);

	}

	uint64_t serialize(std::ostream& out){

		uint64_t w_bytes = 0;

		out.write((char*)&n,sizeof(n));
		out.write((char*)&nbytes,sizeof(nbytes));
		out.write((char*)&n_superblocks,sizeof(n_superblocks));
		out.write((char*)&n_blocks,sizeof(n_blocks));

		w_bytes += sizeof(n) + sizeof(nbytes) + sizeof(n_superblocks) + sizeof(n_blocks);

		out.write((char*)superblock_ranks.data(),n_superblocks*sizeof(p_rank));
		w_bytes += n_superblocks*sizeof(p_rank);

		out.write((char*)data,nbytes*sizeof(uint8_t));
		w_bytes += nbytes*sizeof(uint8_t);

		return w_bytes;

	}

	void load(std::istream& in) {

		in.read((char*)&n,sizeof(n));
		in.read((char*)&nbytes,sizeof(nbytes));
		in.read((char*)&n_superblocks,sizeof(n_superblocks));
		in.read((char*)&n_blocks,sizeof(n_blocks));

		superblock_ranks = vector<p_rank>(n_superblocks);
		in.read((char*)superblock_ranks.data(),n_superblocks*sizeof(p_rank));

		memory = vector<uint8_t>(nbytes+ALN,0);
		data = memory.data();
		while(uint64_t(data) % ALN != 0) data++;
		in.read((char*)data,nbytes*sizeof(uint8_t));

		assert(check_rank());

	}

	uint64_t size(){
		return n;
	}

private:

	/*
	 * find the i-th c. Assumption: the i-th c is inside range [0,n).
	 */
	uint64_t select(uint64_t i, uint64_t begin, uint64_t end, uint8_t c){

		assert(end>begin);

		if(end==begin+1){

			assert(rank(end,c) == rank(begin,c)+1);
			assert(rank(begin,c)==i);
			return begin;

		}

		uint64_t m = (begin+end)/2;
		uint64_t r = rank(m,c);

		if(r > i) 	return select(i,begin,m,c);
		else		return select(i,m,end,c);

	}


	void build_rank_support(){

		p_rank superblock_r = {};
		p_rank block_r = {};

		for(uint64_t bl = 0; bl < n_blocks-1; ++bl){

			uint64_t superblock_number = bl/BLOCKS_PER_SUPERBLOCK;
			uint64_t block_number = bl%BLOCKS_PER_SUPERBLOCK;

			if(block_number == 0){

				superblock_ranks[superblock_number]=superblock_r;
				block_r = {};

			}

			set_counters(superblock_number, block_number,block_r);

			p_rank local_rank = block_rank(superblock_number, block_number);

			block_r = block_r + local_rank;
			superblock_r = superblock_r + local_rank;

		}

		uint64_t superblock_number = (n_blocks-1)/BLOCKS_PER_SUPERBLOCK;
		uint64_t block_number = (n_blocks-1)%BLOCKS_PER_SUPERBLOCK;

		if(block_number == 0){

			superblock_ranks[superblock_number]=superblock_r;
			block_r = {};

		}

		set_counters(superblock_number, block_number,block_r);

		assert(check_rank());

	}

	/*
	 * set i-th block to s. Assumption: s.length() == BLOCK_SIZE
	 */
	void set(uint64_t i, string & s){

		assert(s.length()==BLOCK_SIZE);
		assert(i<n_blocks);

		uint64_t superblock_number = i / BLOCKS_PER_SUPERBLOCK;
		uint64_t block_number = i % BLOCKS_PER_SUPERBLOCK;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the BLOCK_SIZE characters
		__uint128_t* chars = (__uint128_t*)(data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK);

		chars[2] = 0;
		chars[1] = 0;
		chars[0] = 0;

		for(auto c : s){

			if(c==TERM){

				chars[2] = (chars[2] << 1) | __uint128_t(0x1);
				chars[1] = chars[1] << 1;
				chars[0] = chars[0] << 1;

			}else{

				switch(c){

				case 'A' : 	chars[2] = chars[2] << 1;
							chars[1] = chars[1] << 1;
							chars[0] = chars[0] << 1;break;

				case 'C' : 	chars[2] = chars[2] << 1;
							chars[1] = chars[1] << 1;
							chars[0] = (chars[0] << 1) | __uint128_t(0x1);break;

				case 'G' : 	chars[2] = chars[2] << 1;
							chars[1] = (chars[1] << 1) | __uint128_t(0x1);
							chars[0] = chars[0] << 1;break;

				case 'T' : 	chars[2] = chars[2] << 1;
							chars[1] = (chars[1] << 1) | __uint128_t(0x1);
							chars[0] = (chars[0] << 1) | __uint128_t(0x1);break;

				}

			}

		}

	}


	/*
	 * rank in block given as coordinates: superblock, block, offset in block
	 */
	inline p_rank block_rank(uint64_t superblock_number, uint64_t block_number, uint64_t block_off){

		assert(block_off<BLOCK_SIZE);

		if(block_off < 64) return block_rank64(superblock_number, block_number, block_off);

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the BLOCK_SIZE characters
		__uint128_t* chars = (__uint128_t*)(start);

		__uint128_t PAD = ((~__uint128_t(0))>>block_off);

		__uint128_t b = ~(chars[2] | PAD); //most significant bit, padded and negated

		__uint128_t b0 = b & (~chars[1]);
		__uint128_t b1 = b & chars[1];

		return {

			popcount128(b0 & (~chars[0])),
			popcount128(b0 & (chars[0])),
			popcount128(b1 & (~chars[0])),
			popcount128(b1 & (chars[0]))

		};

	}

	/*
	 * rank in block given as coordinates: superblock, block, offset in block
	 */
	inline p_rank block_rank64(uint64_t superblock_number, uint64_t block_number, uint64_t block_off){

		assert(block_off < 64);

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;

		//chars[0,2,4] contains 1st, 2nd, 3rd bits of the BLOCK_SIZE characters
		uint64_t* chars = (uint64_t*)(start);

		uint64_t PAD = ((~uint64_t(0))>>block_off);

		uint64_t b = ~(chars[5] | PAD); //most significant bit, padded and negated

		uint64_t b0 = b & (~chars[3]);
		uint64_t b1 = b & chars[3];

		return {

			(uint64_t)__builtin_popcountll(b0 & (~chars[1])),
			(uint64_t)__builtin_popcountll(b0 & (chars[1])),
			(uint64_t)__builtin_popcountll(b1 & (~chars[1])),
			(uint64_t)__builtin_popcountll(b1 & (chars[1]))

		};

	}

	/*
	 * rank in whole block given as coordinates: superblock, block
	 */
	inline p_rank block_rank(uint64_t superblock_number, uint64_t block_number){

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the BLOCK_SIZE characters
		__uint128_t* chars = (__uint128_t*)(start);

		__uint128_t b2 = ~chars[2]; //most significant bit

		p_rank res = {

			popcount128(b2 & (~chars[1]) & (~chars[0])),
			popcount128(b2 & (~chars[1]) & (chars[0])),
			popcount128(b2 & (chars[1]) & (~chars[0])),
			popcount128(b2 & (chars[1]) & (chars[0]))

		};

		//assert(check_rank_local(res, superblock_number*SUPERBLOCK_SIZE + block_number*BLOCK_SIZE, block_off));

		return res;

	}

	bool check_rank(){

		p_rank p = {};

		bool res = true;

		for(int i=0;i<size();++i){

			auto r = parallel_rank(i);

			if(p != r){

				res = false;
				/*cout << "Error in local rank at position " << n << " "
				<< p.A << "/" << r.A << " "
				<< p.C << "/" << r.C << " "
				<< p.G << "/" << r.G << " "
				<< p.T << "/" << r.T << endl;*/

			}

			p.A += (operator[](i)=='A');
			p.C += (operator[](i)=='C');
			p.G += (operator[](i)=='G');
			p.T += (operator[](i)=='T');

		}

		auto r = parallel_rank(n);

		if(p != r){

			res = false;
			/*cout << "Error in local rank at position " << n << " "
			<< p.A << "/" << r.A << " "
			<< p.C << "/" << r.C << " "
			<< p.G << "/" << r.G << " "
			<< p.T << "/" << r.T << endl;*/

		}

		if(res){

			cout << "rank is correct" << endl;

		}else{

			cout << "rank is not correct" << endl;

		}

		return res;

	}

	/*
	 * check that the string contains exactly the same characters as the file in path
	 */
	bool check_content(string path){

		ifstream ifs(path);

		bool res = true;

		for(uint64_t i=0;i<n;++i){

			char c;
			ifs.read((char*)&c, sizeof(char));

			if(operator[](i) != c) res = false;

		}

		if(res){

			cout << "string content is valid" << endl;

		}else{

			cout << "string content is not valid" << endl;

		}

		return res;

	}

	/*
	 * set counters in the i-th block to r
	 */
	void set_counters(uint64_t superblock_number, uint64_t block_number, p_rank r){

		//block start
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + block_number*BYTES_PER_BLOCK;
		uint32_t * block_ranks = (uint32_t*)(start+48);

		block_ranks[0] = r.A;
		block_ranks[1] = r.C;
		block_ranks[2] = r.G;
		block_ranks[3] = r.T;

		assert(get_counters(superblock_number,block_number) == r);

	}

	/*
	 * get counters of the i-th block
	 */
	inline p_rank get_counters(uint64_t superblock_number, uint64_t superblock_off){

		//block start
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK + superblock_off*BYTES_PER_BLOCK;
		uint32_t * block_ranks = (uint32_t*)(start+48);

		return {
			block_ranks[0],
			block_ranks[1],
			block_ranks[2],
			block_ranks[3]
		};

	}

	char TERM = '#';

	uint64_t n_superblocks = 0;
	uint64_t n_blocks = 0;

	vector<uint8_t> memory; //allocated memory

	//data aligned with blocks of 64 bytes = 512 bits
	uint8_t * data = NULL;

	vector<p_rank> superblock_ranks;

	uint64_t nbytes = 0; //bytes used in data
	uint64_t n = 0;

};


#endif /* INTERNAL_DNA_STRING_HPP_ */
