/*
 * bitVec.cpp
 *
 *  Created on: Feb 27, 2022
 *      Author: bio
 */

#include "bitVec.h"
uint32_t cal_ones(char a)
{
	uint32_t r=0;
	char x=1;
	char y;
	y=1;
	if((a&y)==y)
	{
		r++;
	}
	for(uint32_t i=1;i<=7;i++)
	{
		y=x<<i;
		if((a&y)==y)
		{
			r++;
		}
	}
	return r;
}
void * parallel_gen_vec(void *arg)
{
	struct Para_gen_vec *p=(struct Para_gen_vec *)arg;
	char* v_bit=p->v_bit;
	char* p_mnrPath=p->p_mnrPath;
	uint32_t kmer_len=p->kmer_len;
	uint32_t start=p->start;
	uint32_t end=p->end;
	struct sFMindex FMidx=p->FMidx;

	char * p_e_tmp;
	p_e_tmp=(char*)malloc(sizeof(char)*(kmer_len+2));
	p_e_tmp[kmer_len+1]='\0';

	for(uint64_t i=start;i<=end;i++)
	{
//		cout << i << endl;
		uint8_t label=1;
		for(uint32_t j=0;j<kmer_len+1;j++)
		{
			if(p_mnrPath[i+j]=='?')
			{
				label=0;
				break;
			}
		}
		if(label==1)
		{
			for(uint32_t k=0;k<kmer_len+1;k++)
			{
				p_e_tmp[k]=p_mnrPath[i+k];
			}
			if(check_Minimal_Edge(p_e_tmp,FMidx))
			{
				uint32_t a=cal_edge_id(FMidx,p_e_tmp);
				uint32_t x=a>>3;
				uint32_t y=a&7;
				uint8_t z=1;
				z=z<<y;
				v_bit[x]=v_bit[x]|z;
			}
		}
	}
	free(p_e_tmp);
}
uint32_t check_Minimal_Edge(char* e,sFMindex &FMidx)
{
	//1, if e is a minimal edge
	//0, otherwise
	uint32_t L;
	L=strlen(e);
	char* p;
	p=(char*)malloc(sizeof(char)*L);
	p[L-1]='\0';
	for(uint32_t i=0;i<L-1;i++)
	{
		p[i]=e[i];
	}

	uint32_t out_Degree,in_Degree;
	char adout,adin;
	out_Degree=cal_Node_outdegree(FMidx,p,&adout);
	if(out_Degree<=1)
	{
		free(p);
		return 0;
	}
	else
	{
		for(uint32_t i=0;i<L-1;i++)
		{
			p[i]=e[i+1];
		}
		in_Degree=cal_Node_indegree(FMidx,p,&adin);
		out_Degree=cal_Node_outdegree(FMidx,p,&adout);
		if(in_Degree>1)
		{
			free(p);
			return 1;
		}
		else if(out_Degree>1)
		{
			free(p);
			return 0;
		}
		else
		{
			while(in_Degree==1&&out_Degree==1)
			{
				for(uint32_t i=0;i<L-1;i++)
				{
					p[i]=p[i+1];
				}
				char cc[4];
				cc[0]=1;
				cc[1]=2;
				cc[2]=4;
				cc[3]=8;
				if((cc[0]&adout)==cc[0])
				{
					p[L-2]='A';
				}
				else if((cc[1]&adout)==cc[1])
				{
					p[L-2]='C';
				}
				else if((cc[2]&adout)==cc[2])
				{
					p[L-2]='G';
				}
				else if((cc[3]&adout)==cc[3])
				{
					p[L-2]='T';
				}
				in_Degree=cal_Node_indegree(FMidx,p,&adin);
				out_Degree=cal_Node_outdegree(FMidx,p,&adout);
			}
			if(in_Degree>1)
			{
				free(p);
				return 1;
			}
			else
			{
				free(p);
				return 0;
			}
		}
	}
}
void genBitVec(char *MNRpath, char **p_bitVec, uint64_t *p_bitVec_len,sFMindex &FMidx,uint32_t kmer_len,uint32_t thread_num)
{
	char * p_mnrPath;
	p_mnrPath=MNRpath;

	uint64_t N_dolla=0;
	uint64_t N_total=1;

	uint64_t ii=0;
	while(p_mnrPath[ii]!='\0')
	{
		if(p_mnrPath[ii]=='?')
		{
			N_dolla++;
		}
		N_total++;
		ii++;
	}

	char * p_e_tmp;
	p_e_tmp=(char*)malloc(sizeof(char)*(kmer_len+2));
	p_e_tmp[kmer_len+1]='\0';

	uint64_t N_bit;
	N_bit=N_total-N_dolla;

	*p_bitVec_len=N_bit;
	cout << "N_bit:" << N_bit << endl;

	uint64_t N_Byte;
	N_Byte=N_bit/8;
	if(N_bit%8!=0)
	{
		N_Byte++;
	}

	char* v_bit;
	v_bit=(char*)malloc(sizeof(char)*N_Byte);
	for(uint64_t i=0;i<N_Byte;i++)
	{
		v_bit[i]=0;
	}

//	for(uint64_t i=0;i<N_total-kmer_len;i++)
//	{
////		cout << i << endl;
//		uint8_t label=1;
//		for(uint32_t j=0;j<kmer_len+1;j++)
//		{
//			if(p_mnrPath[i+j]=='?')
//			{
//				label=0;
//				break;
//			}
//		}
//		if(label==1)
//		{
//			for(uint32_t k=0;k<kmer_len+1;k++)
//			{
//				p_e_tmp[k]=p_mnrPath[i+k];
//			}
//			uint64_t a=cal_edge_id(FMidx,p_e_tmp);
//		}
//	}

	if(thread_num==1)
	{
		//check whether the edge is a minimal edge
		for(uint64_t i=0;i<N_total-kmer_len;i++)
		{
			uint8_t label=1;
			for(uint32_t j=0;j<kmer_len+1;j++)
			{
				if(p_mnrPath[i+j]=='?')
				{
					label=0;
					break;
				}
			}
			if(label==1)
			{
				for(uint32_t k=0;k<kmer_len+1;k++)
				{
					p_e_tmp[k]=p_mnrPath[i+k];
				}
				if(check_Minimal_Edge(p_e_tmp,FMidx))
				{
					uint32_t a=cal_edge_id(FMidx,p_e_tmp);
					if(a>=N_bit)
					{
						cout << "error:the id is larger than the total number of edges!" << endl;
					}

					uint32_t x=a>>3;
					uint32_t y=a&7;
					uint8_t z=1;
					z=z<<y;
					v_bit[x]=v_bit[x]|z;
				}
			}
		}
	}
	else
	{
		pthread_t *t;
		t=(pthread_t *)malloc(sizeof(pthread_t)*thread_num);

		struct Para_gen_vec *a;
		a=(struct Para_gen_vec*)malloc(sizeof(struct Para_gen_vec)*thread_num);

		uint32_t unit=N_total/thread_num;
		for(uint32_t i=0;i<thread_num;i++)
		{
			a[i].FMidx=FMidx;
			a[i].kmer_len=kmer_len;
			a[i].v_bit=v_bit;
			a[i].p_mnrPath=p_mnrPath;
			if(i==0)
			{
				a[i].start=0;
			}
			else
			{
				a[i].start=a[i-1].end+1;
			}
			if(i==thread_num-1)
			{
				a[i].end=N_total-kmer_len-1;
			}
			else
			{
				a[i].end=(i+1)*unit-(((i+1)*unit)%8)-1;
			}
		}
		for(uint32_t i=0;i<thread_num;i++)
		{
			if(pthread_create(t+i, NULL, parallel_gen_vec, (void*)(a+i))!=0)
			{
				cout << "error!" << endl;
			}
		}
		for(uint32_t i=0;i<thread_num;i++)
		{
			pthread_join(t[i], NULL);
		}
	}
	*p_bitVec=v_bit;
//	*p_bitVec_len=N_Byte;
}
void genBitVec_total(char* p_file,uint32_t kmer_len,uint32_t thread_num,uint32_t v_bit_gap)
{
	char *MNRpath;

	FILE * file = NULL;
	file = fopen(p_file,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		uint64_t len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(file);//将文件内部的指针重新指向一个流的开头
		MNRpath = (char *)malloc(len+1*sizeof(char));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
		memset(MNRpath,0,len+1);//将内存空间都赋值为‘\0’

		fread(MNRpath, len, 1, file);
		for(uint32_t i=0;i<len;i++)
		{
			if(MNRpath[i]!='A'&&MNRpath[i]!='C'&&MNRpath[i]!='G'&&MNRpath[i]!='T'&&MNRpath[i]!='?')
			{
				MNRpath[i]='A';
			}
		}
		fclose(file);
	}

	char *p_bitVec;
	uint64_t p_bitVec_len;

	sFMindex FMidx;
	char *path=(char*)malloc(sizeof(char)*20);
	path=".";
	read_bfile2index_DM(path,FMidx,0);

	genBitVec(MNRpath, &p_bitVec, &p_bitVec_len,FMidx,kmer_len,thread_num);

	char *file1;
	file1=(char*)malloc(sizeof(char)*6);
	file1[0]='v';
	file1[1]='b';
	char swap[4];
	sprintf(swap,"%d",kmer_len);
	strcpy(file1+2,swap);

	uint32_t c_num=0;
	c_num=p_bitVec_len>>3;
	if(p_bitVec_len&7!=0)
	{
		c_num++;
	}

	FILE *fp;
	fp=fopen(file1,"wb+");
	fwrite(p_bitVec,sizeof(char),c_num,fp);
	fclose(fp);
	free(file1);

	//calculate the Nb array of the bit vector!
	if((v_bit_gap&7)!=0)
	{
		cout << "error: the v_bit_gap should be 8x!" << endl;
	}


	uint32_t unit_bits=v_bit_gap>>3;
	uint32_t * bit_ones;
	bit_ones=(uint32_t*)malloc(sizeof(uint32_t)*(c_num/unit_bits+1));
	uint32_t bit_ones_size=0;

	uint32_t tmp_num_ones=0;
	for(uint32_t i=0;i<c_num;i++)
	{
		tmp_num_ones+=cal_ones(p_bitVec[i]);
		if((i%unit_bits)==0)
		{
			bit_ones[bit_ones_size]=tmp_num_ones;
			bit_ones_size++;
		}
	}
	if((v_bit_gap&7)!=0)
	{
		cout << "error: the v_bit_gap should be 8x!" << endl;
	}

	char *file2;
	file2=(char*)malloc(sizeof(char)*6);
	file2[0]='N';
	file2[1]='b';
	sprintf(swap,"%d",kmer_len);
	strcpy(file2+2,swap);
	FILE *fp1;
	fp1=fopen(file2,"wb+");
	fwrite(bit_ones,sizeof(uint32_t),c_num/unit_bits+1,fp1);
	fclose(fp1);
	free(file2);

	free(p_bitVec);
}

uint32_t is_minimal_edge(char *p_bitVec,uint32_t id)
{
	//0->0,0
	//1->0,1
	//...
	//6->0,6
	//7->0,7
	//...
	uint32_t x=id>>3;
	uint32_t y=id&7;

//	cout << x << ":" << y << endl;
//	cout << p_bitVec[x] << endl;

	char unit=1;
	unit = unit<<y;
	if((p_bitVec[x]&unit)==unit)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int64_t get_minimal_edge_id(char *p_bitVec,uint32_t id,uint32_t v_bit_gap,uint32_t * p_Nb)
{
	//1->0,0
	//2->0,1
	//...
	//7->0,6
	//8->0,7
	//...
	if(is_minimal_edge(p_bitVec,id)==0)
	{
//		cout << "the edge is not a minimal edge!" << endl;
		return -1;
	}
	uint32_t x=id>>3;
	uint32_t y=id&7;

//	if(id==38853275)
//	{
//		cout << "x:y:" << x << ":"<< y<<endl;
//	}
//
//	int64_t r;
//	r=0;
//
//	char one=1;
//	char one_one=1;
//
//	for(uint32_t j=0;j<x;j++)
//	{
//		r=r+cal_ones(p_bitVec[j]);
//	}
//
//	if(id==38853275)
//	{
//		cout << "r:" << r <<endl;
//	}
//
//	for(uint32_t j=0;j<=y;j++)
//	{
//		one_one=one<<j;
//		if((one_one&p_bitVec[x])==one_one)
//		{
//			r++;
//		}
//	}
//	if(id==38853275)
//	{
//		cout << "r:" << r <<endl;
//	}


//	cout << x << ":" << y << endl;
	uint32_t unit=v_bit_gap>>3;

	uint32_t x1=x/unit;
	uint32_t y1=x%unit;

	int64_t r;
	r=p_Nb[x1];

	uint32_t i=1;
	while(x1*unit+i<x)
	{
		r+=cal_ones(p_bitVec[x1*unit+i]);
		i++;
	}

	char one=1;
	char one_one=1;

	if(y1!=0)
	{
		for(uint32_t j=0;j<=y;j++)
		{
			one_one=one<<j;
			if((one_one&p_bitVec[x])==one_one)
			{
				r++;
			}
		}
	}
	else
	{
		for(uint32_t j=y+1;j<=7;j++)
		{
			one_one=one<<j;
			if((one_one&p_bitVec[x])==one_one)
			{
				r--;
			}
		}
	}
	return r-1;
}
