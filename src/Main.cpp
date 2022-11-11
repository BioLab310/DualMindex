/*
 * Main.cpp
 *
 *  Created on: Mar 26, 2021
 *      Author: bio
 */
#include "basic.h"
#include "inputReads.h"
#include "seeding.h"
#include "countMap.h"
#include "inputRef.h"
#include "mappers.h"
#include "DMindex.h"
#include "bitVec.h"
#include "posList.h"
#include "Locate.h"
//#include "outputSam.h"
uint32_t num_read;

int main(int argc, char** argv)
{
	cout << sizeof(struct pos) << endl;
	char *p_read_file;
	char *p_ref_file;
	uint32_t thread_num;
	uint32_t e;
	uint32_t task_size;
	uint32_t read_len;
	uint32_t kmer_len;
	char method;
	char* p_vb_file;
	char * p_Nb_file;
	char * p_Ps_file;
	char * p_Pi_file;
	uint32_t v_bit_gap;
	char *output_file_path;

	struct seed_message seed_test;
	cout << "sizeof(seed_message)" << sizeof(seed_test) << endl;

	for(int32_t i=1;i<argc;i=i+2)
	{
		if(argv[i][0]=='-'&&argv[i][1]=='M')//Method
		{
			method=argv[i+1][0];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='R')//
		{
			p_ref_file = argv[i+1];
		}
		if(argv[i][0]=='-'&&argv[i][1]=='r')//
		{
			p_read_file = argv[i+1];
		}
		if (argv[i][0]=='-'&&argv[i][1]=='r'&&argv[i][2]=='L')
		{
			read_len = atoi(argv[i+1]);
		}
		if (argv[i][0]=='-'&&argv[i][1]=='k'&&argv[i][2]=='L')
		{
			kmer_len = atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='e')//
		{
			e = atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='T')//
		{
			task_size = atoi(argv[i+1]);
		}
		if(argv[i][0]=='-'&&argv[i][1]=='t')//
		{
			thread_num = atoi(argv[i+1]);
		}
		if (argv[i][0]=='-'&&argv[i][1]=='o')
		{
			output_file_path = argv[i + 1];
		}
		if (argv[i][0]=='-'&&argv[i][1]=='V'&&argv[i][2]=='g')
		{
			v_bit_gap = atoi(argv[i + 1]);
		}
		if (argv[i][0]=='-'&&argv[i][1]=='V'&&argv[i][2]=='b')
		{
			p_vb_file = argv[i + 1];
		}
		if (argv[i][0]=='-'&&argv[i][1]=='N'&&argv[i][2]=='b')
		{
			p_Nb_file = argv[i + 1];
		}
		if (argv[i][0]=='-'&&argv[i][1]=='P'&&argv[i][2]=='s')
		{
			p_Ps_file = argv[i + 1];
		}
		if (argv[i][0]=='-'&&argv[i][1]=='P'&&argv[i][2]=='i')
		{
			p_Pi_file = argv[i + 1];
		}
	}

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	cout <<"start..."<<endl;

	switch(method)
	{
		//O : generate the mnrpath,such as aa13
		//V : generate the vb13 and Nb13 files
		//L : generate the position list and pos index
		//F : locate the edge or path
		/******************default setting********************************
			CDBG_para_tmp.v_bit_gap=64;
			CDBG_para_tmp.thread_num=1;
			CDBG_para_tmp.kmer_len=0;
			CDBG_para_tmp.p_ref=NULL;
			CDBG_para_tmp.method='N';
			CDBG_para_tmp.kmer_path=NULL;
		 ****************************************************************/
		case 'O':
		{
			char* nrpath;
			uint64_t nrpath_len;
			cout << "the k-mer length is :" << kmer_len << endl;

			Total_Gen_MNRpath(p_ref_file,kmer_len,&nrpath,&nrpath_len);
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"the time cost of generating the orla path is: "<<span<<endl;
			break;
		}
		case 'o':
		{
			char* nrpath;
			uint64_t nrpath_len;
			cout << "the k-mer length is :" << kmer_len << endl;

			Total_Gen_MNRpath_from_bigger_one(p_ref_file,kmer_len,&nrpath,&nrpath_len);
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"the time cost of generating the orla path is: "<<span<<endl;
			break;
		}
		case 'V':
		{
			genBitVec_total(p_ref_file,kmer_len,thread_num,v_bit_gap);
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"the time cost of generating the bit vector is: "<<span<<endl;
			break;
		}
		case 'L':
		{
			gen_list(p_ref_file,p_vb_file,v_bit_gap,p_Nb_file,kmer_len,thread_num);
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"the time cost of generating the position list is: "<<span<<endl;
			break;
		}
		case 'C':
		{
			uint32_t r;
			r=count_minimal_edges(p_vb_file);
			cout << "the number of minimal edges is:" << r << endl;
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"the time cost of generating the position list is: "<<span<<endl;
			break;
		}
		case 'F':
		{
			sFMindex FMidx;
			char *path=(char*)malloc(sizeof(char)*20);
			path=".";
			read_bfile2index_DM(path,FMidx,0);

			struct DMindex pDM;
			pDM.kmer_len=kmer_len;
			pDM.v_bit_gap=v_bit_gap;
			pDM.FMidx=FMidx;
			read_Pi(&(pDM.p_Pi),p_Pi_file);
			read_Ps(&(pDM.p_Ps),p_Ps_file);
			read_Nb(&(pDM.p_Nb),p_Nb_file);
			read_vb(&(pDM.p_bitVec),p_vb_file);

			struct RefFilePath p_ref_path;
			getRefFilePathes(p_ref_file, &p_ref_path);

			char *seq;
			uint64_t seq_length;

			uint32_t x;
			x=pDM.kmer_len+1;
			struct position_list r;
			char* h;
			h=(char*)malloc(sizeof(char)*(x+1));
			h[x]='\0';

			for(uint32_t ref_i=0;ref_i<1;ref_i++)
			{
				ReadSeq_ref(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
				char *temp_seq;
				temp_seq=seq;

				for(uint32_t i=0;i<100;i++)
				{
					for(uint32_t j=0;j<x;j++)
					{
						h[j]=temp_seq[i+j];
					}
//					Locate_path(pDM,h,&r);
					Locate_edge_whole(pDM,h,&r);
					cout << "i:" << i << endl;
					cout << h << endl;
					cout << r.p_pos_list_len << endl;
//					for(uint32_t j=0;j<r.p_pos_list_len;j++)
//					{
//						cout << (int)r.p_pos_list[j].ref_id << ":" <<r.p_pos_list[j].pos << endl;
//					}
				}
			}
			gettimeofday(&tve,NULL);
			double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
			cout <<"the time cost of locating edge is: "<<span<<endl;
			break;
		}
	}


	//the example of calling "input_read" method
//	char *p_reads;
//	uint32_t start_pos=0;
//	input_reads(&p_reads,&num_read,p_read_path,start_pos,task_size,len_read);
//
//	char a[]="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
//	cout << p_reads << endl;
//	sse3_convert2bit1(a+16, read_vec0_t, read_vec1_t);
//
//	char * p_ref;
//	uint64_t ref_length;
//	ReadSeq(&p_ref,&ref_length,p_ref_path);
//
//	//test for the "seed_indexing" method
//	struct NodeBit** p_BplusTrees;
//	vector <struct seed_message> VS_message;
//	uint32_t *read_sav_seed;
//	read_sav_seed=read_seed_build(num_read,e);
//	p_BplusTrees=seed_indexing(p_reads,num_read,len_read,e,thread_num,VS_message,read_sav_seed);

//	cout <<VS_message.size() << endl;
//	for(uint32_t i=0;i<task_size*(e+1);i++)
//	{
//		cout <<read_sav_seed[i] << endl;
//	}


//	struct naiveMapper_para p;
//	p.VS_message=VS_message;
//	p.e=e;
//	p.len_read=len_read;
//	p.p_BplusTrees=p_BplusTrees;
//	p.p_reads=p_reads;
//	p.p_ref=p_ref;
//	p.read_sav_seed=read_sav_seed;
//	p.ref_length=ref_length;
//	p.task_size=task_size;
//	p.thread_num=thread_num;
//	naiveMapper(p);

//	countMapForSAI(VS_message);
//	ofstream output("seed_sa.txt");
//	for(uint32_t i=0;i<VS_message.size();i++)
//	{
//		output<<VS_message[i].seed_sa<<endl;
//	}
//
//	FILE* outUnipath_c;
//	outUnipath_c=fopen("unipath_counter","wb+");
//	fwrite(read_sav_seed,sizeof(uint32_t),num_read*(e+1),outUnipath_c);

//	// 输出 Sam
//	struct OutputQueue output_queue;
//	// Load reference
//	SequenceBatch reference_sequence_batch;
//	initialize_sequence_batch(&reference_sequence_batch);
//	initialize_sequence_batch_loading(p_ref_path, &reference_sequence_batch);
//	load_all_sequences_into_sequence_batch(&reference_sequence_batch);
//	// initialize output_queue
//	uint32_t output_queue_max_size = 100000;
//	initialize_output_queue(output_file_path, &reference_sequence_batch, thread_num, output_queue_max_size, &output_queue);
//	//	void *p;  //指向指针函数
//	void *output_queue_v;
//	output_queue_v = (void*)malloc(sizeof(struct OutputQueue));
//	output_queue_v = &output_queue;
//	fprintf(stderr, "output: %s\n", output_file_path);
//	output_sam(output_queue_v);

	cout << "end..."<< endl;
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"total time is: "<<span<<endl;
}

