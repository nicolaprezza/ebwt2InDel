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

bool switch_ = false;

void help(){

	cout << "snp2fastq calls.snp [-i]" << endl << endl <<
	"Converts clust2snp's calls 'calls.snp' into a fastq  file 'calls.snp.fastq'. The  output contains  one " << endl <<
	"read per call, where we put the second individual's DNA in the read's name, and the first individual's " << endl <<
	"DNA in the read DNA. Base qualities are fake (all maximum). If option -i is specified, then individuals" << endl <<
	"are switched." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc != 2 and argc != 3) help();

	if(argc == 3){

		if(string(argv[2]).compare("-i")==0){
			switch_=true;
		}else{
			help();
		}

	}

	string infile = argv[1];

	string outfile = infile;
	outfile.append(".fastq");

	ifstream is(infile);
	ofstream of(outfile);

	string str;
	unsigned int idx=0;

	string header;
	string dna;

	string event_type;
	string event_number;
	string snp_pos;
	string cov0;
	string cov1;
	string event;

	while(getline(is, str)){

		if( ((switch_) and idx%4==2) or ((not switch_) and idx%4==0)  ){//header

			//if switch is false, this is activated on line number 0 (header first indiv)
			//if switch is true, this is activated on line number 2 (header second indiv)

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

			header = event_type;
			header += "_" + event_number + "_" + snp_pos + "_" + event + "_" + cov0;


		}

		if( ((switch_) and idx%4==3) or ((not switch_) and idx%4==1)  ){//DNA

			//if switch is false, this is activated on line number 1 (DNA first indiv)
			//if switch is true, this is activated on line number 3 (DNA second indiv)

			dna = str;

		}

		if( ((switch_) and idx%4==0) or ((not switch_) and idx%4==2)  ){//header

			//if switch is false, this is activated on line number 2 (header second indiv)
			//if switch is true, this is activated on line number 0 (header first indiv)

			std::istringstream iss_bar(str);
			getline(iss_bar, cov1, '|');
			getline(iss_bar, cov1, '|');
			getline(iss_bar, cov1, '|');

			header += "_" + cov1 + "_";

		}

		if( ((switch_) and idx%4==1) or ((not switch_) and idx%4==3)  ){//DNA

			//if switch is false, this is activated on line number 3 (DNA second indiv)
			//if switch is true, this is activated on line number 1 (DNA first indiv)

			header += str;

			//now output fastq entry
			of 	<< "@" << header << endl <<
			dna << endl <<
			"+" << endl <<
			string(dna.length(),'I') << endl;



		}

		idx++;

	}

	is.close();
	of.close();

}
