#include <fstream>

using namespace std;

typedef uint32_t int_text;
typedef uint32_t int_suff;
typedef uint32_t int_lcp;
typedef uint8_t int8;
typedef uint8_t dataTypelenSeq;	//length of the sequences (in biologic case 100)
typedef uint32_t dataTypeNSeq;	//number of sequences
typedef __int128 uint128_t;

typedef pair<uint64_t,uint32_t> coordinate;//suffix array coordinate (text, suff)



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

