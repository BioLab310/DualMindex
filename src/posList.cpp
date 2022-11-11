/*
 * posList.cpp
 *
 *  Created on: Mar 23, 2022
 *      Author: bio
 */

#include "posList.h"

struct pos **p_pos;
uint32_t *p_pos_len;
uint32_t p_pos_total_len=0;

void * parallel_gen_pos_list(void *arg)
{
	struct Para_gen_list *p=(struct Para_gen_list *)arg;

	sFMindex FMidx=p->FMidx;
	char * seq=p->seq;
	uint64_t seq_length=p->seq_length;
	char *p_bitVec=p->p_bitVec;
	uint32_t v_bit_gap=p->v_bit_gap;
	uint32_t* p_Nb=p->p_Nb;
	uint32_t kmer_len=p->kmer_len;
	uint32_t thread_id=p->thread_id;
	uint32_t thread_num=p->thread_num;
	uint8_t ref_i=p->ref_i;

	char *temp_seq;
	temp_seq=seq;

	for(uint32_t i=0;i<seq_length-kmer_len;i++)
	{
		uint32_t id;
		id=cal_edge_id_new(FMidx,temp_seq+i,kmer_len+1);

		int64_t mini_id=get_minimal_edge_id(p_bitVec,id,v_bit_gap,p_Nb);
		if(mini_id!=-1&&mini_id%thread_num==thread_id)
		{
			if(mini_id<p_pos_total_len)
			{
				if(p_pos_len[mini_id]==0)
				{
					p_pos[mini_id]=(struct pos *)malloc(sizeof(struct pos)*1);
					p_pos_len[mini_id]=1;
				}
				else
				{
					p_pos[mini_id]=(struct pos *)realloc(p_pos[mini_id],sizeof(struct pos)*(p_pos_len[mini_id]+1));
					p_pos_len[mini_id]++;
				}
				p_pos[mini_id][p_pos_len[mini_id]-1].pos=i;//cyy
				p_pos[mini_id][p_pos_len[mini_id]-1].ref_id=ref_i;
			}
			else
			{
//				if(p_pos_total_len==0)
//				{
//					p_pos=(struct pos **)malloc(sizeof(struct pos *)*(mini_id+1));
//					p_pos_len=(uint32_t*)malloc(sizeof(uint32_t)*(mini_id+1));
//					for(uint32_t j=0;j<=mini_id;j++)
//					{
//						p_pos[j]=NULL;
//						p_pos_len[j]=0;
//					}
//
//				}
//				else
//				{
//					p_pos=(struct pos **)realloc(p_pos,sizeof(struct pos *)*(mini_id+1));
//					p_pos_len=(uint32_t*)realloc(p_pos_len,sizeof(uint32_t)*(mini_id+1));
//					for(uint32_t j=p_pos_total_len;j<=mini_id;j++)
//					{
//						p_pos[j]=NULL;
//						p_pos_len[j]=0;
//					}
//
//				}
//				p_pos_total_len=mini_id+1;
//
//				p_pos[mini_id]=(struct pos *)malloc(sizeof(struct pos)*1);
//				p_pos_len[mini_id]=1;
//
//				p_pos[mini_id][p_pos_len[mini_id]-1].pos=i;
//				p_pos[mini_id][p_pos_len[mini_id]-1].ref_id=ref_i;
				cout << "the id of the minimal edge is larger than the total number of the minimal edges!" << endl;
				cout << mini_id << ":"<< p_pos_total_len;

			}
		}
	}
}
void read_Nb(uint32_t** p_Nb,char* p_Nb_file)
{
	uint32_t *p;
	FILE * file = NULL;
	file = fopen(p_Nb_file,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		uint64_t len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		p = (uint32_t *)malloc(sizeof(uint32_t)*(len/4));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况

		fread(p, sizeof(uint32_t),len/4,  file);
		fclose(file);
	}
	*p_Nb=p;
}
void read_vb(char** p_bitVec,char* p_bitVec_file)
{
	char *p;
	FILE * file = NULL;
	file = fopen(p_bitVec_file,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		uint64_t len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		p = (char *)malloc(sizeof(char)*len);//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
		cout << "len:" << len << endl;

		fread(p, sizeof(char),len,  file);
		fclose(file);
	}
	*p_bitVec=p;
}
void read_Pi(uint32_t** p_Pi,char* p_Pi_file)
{
	uint32_t * p;
	FILE * file = NULL;
	file = fopen(p_Pi_file,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		uint64_t len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		p = (uint32_t *)malloc(sizeof(uint32_t)*(len/4));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况

		fread(p, sizeof(uint32_t),len/4,  file);
		fclose(file);
	}
	*p_Pi=p;
}
void read_Ps(struct pos** p_Ps,char* p_Ps_file)
{
	struct pos *p;
	FILE * file = NULL;
	file = fopen(p_Ps_file,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		uint64_t len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		p = (struct pos *)malloc(sizeof(struct pos)*(len/sizeof(struct pos)));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况

		fread(p, sizeof(struct pos),len/sizeof(struct pos), file);
		fclose(file);
	}
	*p_Ps=p;
}
void generate_pos_list(sFMindex &FMidx,char * ref_name,char *p_bitVec,uint32_t v_bit_gap,uint32_t* p_Nb,uint32_t kmer_len,uint32_t thread_num_log,uint32_t minimal_edge_num)
{
	uint32_t thread_num=pow(2,thread_num_log);
	struct RefFilePath p_ref_path;
	getRefFilePathes(ref_name, &p_ref_path);
	char *seq;
	uint64_t seq_length;

	p_pos_total_len=minimal_edge_num;
	p_pos=(struct pos **)malloc(sizeof(struct pos *)*(minimal_edge_num));
	p_pos_len=(uint32_t*)malloc(sizeof(uint32_t)*(minimal_edge_num));
	for(uint32_t j=0;j<p_pos_total_len;j++)
	{
		p_pos[j]=NULL;
		p_pos_len[j]=0;
	}

	uint32_t kk=0;
	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		cout << ref_i+1 << ":"<< p_ref_path.NumberOfPathes <<endl;
//		if(ref_i<43)
//		{
//			continue;
//		}
		ReadSeq_ref(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		if(thread_num==1)
		{
			char *temp_seq;
			temp_seq=seq;

			char tmp;
			for(uint32_t i=0;i<seq_length-kmer_len;i++)
			{
				uint32_t id;
				id=cal_edge_id_new(FMidx,temp_seq+i,kmer_len+1);

				int64_t mini_id=get_minimal_edge_id(p_bitVec,id,v_bit_gap,p_Nb);
				if(mini_id!=-1)
				{
					if(mini_id<p_pos_total_len)
					{
						if(p_pos_len[mini_id]==0)
						{
							p_pos[mini_id]=(struct pos *)malloc(sizeof(struct pos)*1);
							p_pos_len[mini_id]=1;
						}
						else
						{
							p_pos[mini_id]=(struct pos *)realloc(p_pos[mini_id],sizeof(struct pos)*(p_pos_len[mini_id]+1));
							p_pos_len[mini_id]++;
						}
						p_pos[mini_id][p_pos_len[mini_id]-1].pos=i;//cyy
						p_pos[mini_id][p_pos_len[mini_id]-1].ref_id=ref_i;
					}
					else
					{
//						if(p_pos_total_len==0)
//						{
//							p_pos=(struct pos **)malloc(sizeof(struct pos *)*(mini_id+1));
//							p_pos_len=(uint32_t*)malloc(sizeof(uint32_t)*(mini_id+1));
//							for(uint32_t j=0;j<=mini_id;j++)
//							{
//								p_pos[j]=NULL;
//								p_pos_len[j]=0;
//							}
//
//						}
//						else
//						{
//							p_pos=(struct pos **)realloc(p_pos,sizeof(struct pos *)*(mini_id+1));
//							p_pos_len=(uint32_t*)realloc(p_pos_len,sizeof(uint32_t)*(mini_id+1));
//							for(uint32_t j=p_pos_total_len;j<=mini_id;j++)
//							{
//								p_pos[j]=NULL;
//								p_pos_len[j]=0;
//							}
//
//						}
//						p_pos_total_len=mini_id+1;
//
//						p_pos[mini_id]=(struct pos *)malloc(sizeof(struct pos)*1);
//						p_pos_len[mini_id]=1;
//
//						p_pos[mini_id][p_pos_len[mini_id]-1].pos=i;
//						p_pos[mini_id][p_pos_len[mini_id]-1].ref_id=ref_i;
						cout << "error: the id of minimal edge is larger than the total number of minimal edges!" << endl;
						cout << mini_id << ":" <<p_pos_total_len << endl;
						cout << "id" << ":" <<id << endl;

					}
				}
			}
		}
		else
		{
			pthread_t *t;
			t=(pthread_t *)malloc(sizeof(pthread_t)*thread_num);

			struct Para_gen_list *a;
			a=(struct Para_gen_list*)malloc(sizeof(struct Para_gen_list)*thread_num);

			for(uint32_t i=0;i<thread_num;i++)
			{
				a[i].FMidx=FMidx;
				a[i].kmer_len=kmer_len;
				a[i].p_Nb=p_Nb;
				a[i].p_bitVec=p_bitVec;
				a[i].seq=seq;
				a[i].seq_length=seq_length;
				a[i].thread_id=i;
				a[i].thread_num=thread_num;
				a[i].v_bit_gap=v_bit_gap;
				a[i].ref_i=ref_i;
			}
			for(uint32_t i=0;i<thread_num;i++)
			{
				if(pthread_create(t+i, NULL, parallel_gen_pos_list, (void*)(a+i))!=0)
				{
					cout << "error!" << endl;
				}
			}
			for(uint32_t i=0;i<thread_num;i++)
			{
				pthread_join(t[i], NULL);
			}
		}
	}


	char swap[4];
	char *file2;
	file2=(char*)malloc(sizeof(char)*6);
	file2[0]='P';
	file2[1]='s';
	sprintf(swap,"%d",kmer_len);
	strcpy(file2+2,swap);
	FILE *fp1;
	fp1=fopen(file2,"wb+");
	for(uint32_t i=0;i<p_pos_total_len;i++)
	{
		fwrite(p_pos[i],sizeof(struct pos),p_pos_len[i],fp1);
	}
	fclose(fp1);
	free(file2);

	char *file3;
	file3=(char*)malloc(sizeof(char)*6);
	file3[0]='P';
	file3[1]='i';
	sprintf(swap,"%d",kmer_len);
	strcpy(file2+2,swap);
	FILE *fp2;
	fp2=fopen(file3,"wb+");
	uint32_t *p_index;
	p_index=(uint32_t*)malloc(sizeof(uint32_t)*(p_pos_total_len+1));
	p_index[0]=0;

	for(uint32_t i=1;i<=p_pos_total_len;i++)
	{
		p_index[i]=p_index[i-1]+p_pos_len[i-1];
	}
	fwrite(p_index,sizeof(uint32_t),p_pos_total_len+1,fp2);

	fclose(fp2);
	free(file3);
	free(p_index);

	for(uint32_t i=0;i<p_pos_total_len;i++)
	{
		free(p_pos[i]);
	}
	free(p_pos);
	free(p_pos_len);
}
void gen_list(char * ref_name,char *p_bitVec_file,uint32_t v_bit_gap,char* p_Nb_file,uint32_t kmer_len,uint32_t thread_num)
{
	sFMindex FMidx;
	char *path=(char*)malloc(sizeof(char)*20);
	path=".";
	read_bfile2index(path,FMidx,0);

//	char*a="CCATGATTCCATAA";
//	cal_edge_id(FMidx,a);

	uint32_t minimal_edge_num;
	minimal_edge_num=count_minimal_edges(p_bitVec_file);

	cout << "minimal_edge_num:" << minimal_edge_num << endl;

	uint32_t * p_Nb;
	read_Nb(&p_Nb,p_Nb_file);

	char* p_bitVec;
	read_vb(&p_bitVec,p_bitVec_file);

	generate_pos_list(FMidx,ref_name,p_bitVec,v_bit_gap,p_Nb,kmer_len,thread_num,minimal_edge_num);
}

uint32_t count_minimal_edges(char *p_bitVec_file)
{
	int64_t len;
	uint32_t r;
	r=0;
	char *p;
	FILE * file = NULL;
	file = fopen(p_bitVec_file,"rb");
	if(file==NULL)
	{
		cout <<"file can not be open!" << endl;
		return 0;
	}
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		p = (char *)malloc(sizeof(char)*len);//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况

		fread(p, sizeof(char),len,  file);
		fclose(file);
	}

	cout << len << endl;

	for(uint32_t i=0;i<len;i++)
	{
//		uint32_t tmp;
//		tmp=cal_ones(p[i]);
		r=r+cal_ones(p[i]);
//		if(i==4856658)
//		{
//			cout << "r:"<< r <<endl;
//		}
//		if(i==4856659)
//		{
//			cout << "r:"<< r <<endl;
//			cout << "p[i]:" << (uint32_t)p[i]<< endl;
//			cout << "tmp:" << tmp << endl;
//		}
	}

	free(p);

	return r;
}
