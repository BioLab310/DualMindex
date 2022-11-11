/*
 * Main.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: bio
 */
#include "basic.h"
#include "build.h"
#include "read.h"
#include "ExactMatch.h"
#include "Bit_sseInstruction.h"
#include <ctime>
#include "BitVec_index.h"
#include "print.h"
#include <bitset>

int main(int argc, char** argv)
{
	struct build_para para;
	const char *seqtest;
	for(int32_t i=1;i<argc;i=i+2)
	{
		if(argv[i][0]=='-'&&argv[i][1]=='h')	//path of reference
		{
			para.ref_path = argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='L') // the length of the reference for testing the index
		{
			para.kmer_len = atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='s')//sa_gap
		{
			para.sa_gap = atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='o')//occ_gap
		{
			para.occ_gap = atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='l')//level or the length of string in testing the index
		{
			para.level = atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='t')//thread_num
		{
			para.thread_num = atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='m')//max_len of cmp
		{
			para.max_len = atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='S')//max_len of cmp
		{
			seqtest = argv[i+1];
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='M')//max_len of cmp
		{
			para.m = argv[i+1][0];
		}
	}
	struct timeval tvs,tve;
	double span;

	switch(para.m)
	{
		//B : build the index
		//T : test the index
		//I : build the compressed
		case 'B':
		{
			char *ndir = ".";
			if(access(ndir,0) == -1)
			{
				printf("%s not exist.\n",ndir);
				mkdir(ndir, 0755);
			}
			build(para, ndir);
			build_occA(para, ndir);
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"the time cost of building the index is: "<<span<<endl;
			cout << "occ_gap:" << para.occ_gap << "sa_gap:" << para.sa_gap << "'s FMindex is done!" << endl;
			break;
		}
		case 'T':
		{
			char *path=(char*)malloc(sizeof(char)*20);
			path=".";
			uint32_t *bet;

			sFMindex FMidx;
			read_bfile2index(path,FMidx,0);
			cout << FMidx.occ_gap << endl;
			cout << FMidx.sa_gap << endl;
			cout << "read_bfile2index done..." << endl;
			gettimeofday(&tvs,NULL);

			char * line;
			line=(char*)malloc(sizeof(char)*(para.level+1));
			line[para.level]='\0';
			char * seq;
			uint32_t seq_length;
			ReadSeq(&seq,&seq_length,para.ref_path);
		//	kmer_len=l;

			for(uint32_t i=0;i<para.kmer_len;i++)
			{
				strncpy(line,seq+i,para.level);
				uint32_t label=1;
				for(uint32_t j=0;j<para.level;j++)
				{
					if(line[j]!='A'&&line[j]!='C'&&line[j]!='G'&&line[j]!='T')
					{
						label=0;
						break;
					}
				}
				if(label==0)
				{
					continue;
				}
				bet=calc_SArangeSeq(FMidx,line);
				if(bet[0] < bet[1])
				{
					for(uint32_t j=bet[0];j<=bet[1];j++)
					{
						uint32_t sa;
						sa=calc_SA(FMidx,j);

		//				if(i!=sa)
		//				{
						cout <<i<<":"<< sa << endl;
		//				}
					}
				}
			}
			gettimeofday(&tve,NULL);
			span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"time of ExactMatch is: "<<span<<endl;
			break;
		}
		case 'I':
		{
			char *ndir = ".";
			if(access(ndir,0) == -1)
			{
				printf("%s not exist.\n",ndir);
				mkdir(ndir, 0755);
			}
			build_C("C","cC",para.sa_gap,para.occ_gap);
			build_occA(para,ndir);
			build_CompactedSA_from_SA("SA","cSA",para.sa_gap);

			cout << "occ_gap:" << para.occ_gap << "sa_gap:" << para.sa_gap << "'s FMindex is done!" << endl;

			gettimeofday(&tve,NULL);
			span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"time cost of build the compressed index is: "<<span<<endl;
			break;
		}
	}
}


