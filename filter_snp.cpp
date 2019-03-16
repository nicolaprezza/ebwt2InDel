// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <set>
#include <cstring>

using namespace std;

void help(){

	cout << "filter_snp calls.snp m [M]" << endl << endl <<
	"Input: a .snp file. Keep only reads with at least coverage m and at most M. Output to stdout." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc != 3 and argc != 4) help();

	string infile = argv[1];
	int m = atoi(argv[2]);
	int M = 0;

	if(argc==4) M = atoi(argv[3]);

	ifstream is(infile);

	string str;
	unsigned int idx=0;

	int cov=0;

	string header;

	while(getline(is, str)){

		if(idx%2==0){//header

			header = str;

			std::istringstream iss1(str);
			std::string token;
			getline(iss1, token, '_');//cluster
			getline(iss1, token, '_');//id
			getline(iss1, token, '_');//right
			getline(iss1, token, '_');//cov

			{
				std::istringstream iss2(token);
				getline(iss2, token, ':');//cov
				getline(iss2, token, ':');//value
				cov = atoi(token.c_str());
			}

		}else{

			if(cov>=m and (M==0 or cov <= M)){

				cout << header << endl << str << endl;

			}

			header="";
			cov=0;

		}

		idx++;

	}

	is.close();

}
