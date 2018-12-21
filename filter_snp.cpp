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

	cout << "filter_snp calls.snp M" << endl << endl <<
	"Input: a .snp file. Filters out only pairs with at least coverage M (in both variants). Output to stdout." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc != 3) help();

	string infile = argv[1];
	int M = atoi(argv[2]);

	ifstream is(infile);

	string str;
	unsigned int idx=0;

	string event_type;
	string event_number;
	string snp_pos;
	string cov0;
	string cov1;
	string event;

	string line1;
	string line2;
	string line3;
	string line4;

	while(getline(is, str)){

		if(idx%4==0){//first line of call

			line1 = str;

			std::istringstream iss_bar(str);
			std::string token;
			getline(iss_bar, token, '|');
			token = token.substr(1);

			{
				std::istringstream iss_underscore(token);
				getline(iss_underscore, event_type, '_');
				getline(iss_underscore, token, '_');
				getline(iss_underscore, token, '_');
				getline(iss_underscore, event_number, '_');
			}

			getline(iss_bar, token, '|');

			{
				std::istringstream iss_dots(token);
				getline(iss_dots, token, ':');
				getline(iss_dots, token, ':');
				std::istringstream iss_underscore(token);
				getline(iss_underscore, snp_pos, '_');
				getline(iss_underscore, event, '_');

			}

			getline(iss_bar, cov0, '|');

		}

		if(idx%4==1){//DNA of first individual (the reference)

			line2 = str;

		}

		if(idx%4==2){//header second individual

			line3 = str;

			std::istringstream iss_bar(str);
			getline(iss_bar, cov1, '|');
			getline(iss_bar, cov1, '|');
			getline(iss_bar, cov1, '|');

		}

		if(idx%4==3){//DNA of second individual (ALT)

			line4 = str;

			if(atoi(cov0.c_str())>=M and atoi(cov1.c_str())>=M){

				cout << line1 << endl << line2 << endl << line3 << endl << line4 << endl;

			}


		}

		idx++;

	}

	is.close();

}
