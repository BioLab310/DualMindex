/*
 * BitVec_build.cpp
 *
 *  Created on: Jan 4, 2020
 *      Author: pluto
 */
#include "BitVec_index.h"
#include "Bit_sseInstruction.h"
#include "build.h"
#include "read.h"

void build_bVec(sFMindexv *pindxv)
{
	char * b = NULL;
	uint32_t len = 0;

	FILE * B_file = NULL;
	B_file = fopen("./index/B","r");

	if(B_file)//打开文件一定要判断是否成功
	{
		fseek(B_file,0,SEEK_END);//将文件内部的指针指向文件末尾
		len = ftell(B_file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(B_file);//将文件内部的指针重新指向一个流的开头
		pindxv->bseqlen = len;
		b = new char[len+1];//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
		memset(b,'\0',(len+1)*sizeof(char));//将内存空间都赋值为‘\0’

		fread(b,sizeof(char),len,B_file);
		fclose(B_file);
	}

	uint32_t unitnum;
	unitnum = len / 64;
	if(len % 64)
	{
		unitnum++;
	}
	if(unitnum % 2)
	{
		unitnum++;
	}
	uint64_t *bSeqvh = new uint64_t[unitnum]();
	uint64_t *bSeqvl = new uint64_t[unitnum]();
	char *b_tmp = b;
	uint32_t unit_index = 0;
	uint32_t loop = 0;
	uint32_t dol_pos = 0;
	while(*b_tmp != '\0')
	{
		switch(*b_tmp)
		{
			case '$':
				dol_pos = b_tmp - b;
				pindxv->dolpos = dol_pos;
				break;
			case 'A':
				break;
			case 'C':
				bSeqvl[unit_index] |= U64HSET1;
				break;
			case 'G':
				bSeqvh[unit_index] |= U64HSET1;
				break;
			case 'T':
				bSeqvl[unit_index] |= U64HSET1;
				bSeqvh[unit_index] |= U64HSET1;
				break;
			default:
				break;
		}
		if(loop == 63)
		{
			unit_index++;
			loop = 0;
			b_tmp++;
			continue;
		}
		else
		{
			bSeqvl[unit_index] >>= 1;
			bSeqvh[unit_index] >>= 1;
		}
		loop++;
		b_tmp++;
	}
	if(loop != 0)
	{
		bSeqvl[unit_index] >>= 63-loop;
		bSeqvh[unit_index] >>= 63-loop;
	}
	FILE * B_fileVec = NULL;
	B_fileVec = fopen("./index/BVec","w");
	fwrite(bSeqvh,sizeof(uint64_t),unitnum,B_fileVec);
	fwrite(bSeqvl,sizeof(uint64_t),unitnum,B_fileVec);
	fclose(B_fileVec);
	pindxv->bvh = (__m128i*)bSeqvh;
	pindxv->bvl = (__m128i*)bSeqvl;
	pindxv->u64num = unit_index;
//	delete []bVech;
//	delete []bVecl;
	delete b;
}

void build_occAVec(sFMindexv *pindxv, struct build_para para)
{
	uint64_t *bVech = NULL;
	uint64_t *bVecl = NULL;
	uint32_t len = 0;

	FILE * BVec_file = NULL;
	BVec_file = fopen("./index/BVec","r");
	uint32_t unitnum;
	uint32_t bVecnum;
	if(BVec_file)//打开文件一定要判断是否成功
	{
		fseek(BVec_file,0,SEEK_END);//将文件内部的指针指向文件末尾
		len = ftell(BVec_file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		rewind(BVec_file);//将文件内部的指针重新指向一个流的开头
		unitnum = len / sizeof(uint64_t);
		unitnum /= 2;
		if(unitnum % 2)
		{
			bVecnum = unitnum + 1;
		}
		else
		{
			bVecnum = unitnum;
		}

		bVech = new uint64_t[bVecnum]();
		bVecl = new uint64_t[bVecnum]();

		fread(bVech,sizeof(uint64_t),unitnum,BVec_file);
		fread(bVecl,sizeof(uint64_t),unitnum,BVec_file);
		fclose(BVec_file);
	}
	uint32_t *a;
	a = new uint32_t[4]();
	__m128i regh, regl;
	FILE * out_OCCAV = NULL;
	out_OCCAV = fopen("./index/OCCAV","wb+");
	uint32_t pos = 0;
	uint32_t remain = 0;
	uint32_t dol_find = 0;
	uint32_t mloop = 0;
	fwrite(a,sizeof(uint32_t)*4,1,out_OCCAV);
	for(uint32_t i = 0; i < bVecnum/2; i++)
	{
		regh = *((__m128i *)&bVech[2*i]);
		regl = *((__m128i *)&bVecl[2*i]);
		if(i == bVecnum/2-1)
		{
			remain = len % 64;
			if(unitnum % 2)
			{
				a[0] = calcVec_Acut(regh, regl, remain, CUTVECL);
			}
			else
			{
				a[0] = calcVec_Acut(regh, regl, remain+64, CUTVECL);
			}
		}
		else
		{
			a[0] += calcVec_Anum(regh, regl);
		}
		a[1] += calcVec_Cnum(regh, regl);
		a[2] += calcVec_Gnum(regh, regl);
		a[3] += calcVec_Tnum(regh, regl);
		pos += 128;
		mloop++;
		if(pos > pindxv->dolpos && dol_find == 0)
		{
			a[0]--;
			dol_find = 1;
		}
		if(remain == 0 && mloop == para.occ_gap / 128)
		{
			fwrite(a,sizeof(uint32_t)*4,1,out_OCCAV);
			mloop = 0;
		}
	}
	fclose(out_OCCAV);
	delete []bVech;
	delete []bVecl;

}

uint32_t calc_OCCV(sFMindexv indxv,struct build_para para,char c,uint32_t pos)
{
	uint32_t offset;
	uint32_t (*p)(__m128i, __m128i);
	uint32_t (*q)(__m128i, __m128i , uint32_t, uint64_t *);
	switch(c)
	{
		case 'A' :
				offset = 0;
				p = calcVec_Anum;
				q = calcVec_Acut;
				break;
		case 'C' :
				offset = 1;
				p = calcVec_Cnum;
				q = calcVec_Ccut;
				break;
		case 'G' :
				offset = 2;
				p = calcVec_Gnum;
				q = calcVec_Gcut;
				break;
		case 'T' :
				offset = 3;
				p = calcVec_Tnum;
				q = calcVec_Tcut;
				break;
		default:
				offset = -1;
				break;
	}
	if(pos % para.occ_gap == 0)
	{
		return *(indxv.occa + (pos/para.occ_gap)*4 + offset);
	}
	else
	{
		uint32_t rslt = *(indxv.occa + (pos/para.occ_gap)*4 + offset);
		uint32_t remain = pos % para.occ_gap;
		uint32_t loop = remain/128 + 1;
		uint32_t index = (pos/para.occ_gap) * (para.occ_gap/128);
		for(uint32_t i = 0; i < loop; i++)
		{
			if(i == loop-1)
			{
				rslt += (*q)(indxv.bvh[index+i],indxv.bvl[index+i],remain,CUTVECL);
			}
			else
			{
				rslt += (*p)(indxv.bvh[index+i],indxv.bvl[index+i]);
			}
			remain -= 128;
		}
		if(offset == 0 && pos > indxv.dolpos && pos < (indxv.dolpos/para.occ_gap+1)*para.occ_gap)
		{
			rslt--;
		}
		return rslt;
	}
}

inline uint32_t calc_OCCVu(sFMindexv indxv,struct build_para para,char c,uint32_t pos, uint32_t methodsel)
{
	if(methodsel == 1)
	{

	}
	else if(methodsel == 2)
	{
		return calc_OCCV(indxv,para,c,pos);
	}
}
int mainbi(int argc, char** argv)
{
	init_CUTVECL();
	init_CUTVECH();
	struct build_para para;
	char c;
	for(int32_t i=1;i<argc;i=i+2)
	{
		if(argv[i][0]=='-'&&argv[i][1]=='h')//Method
		{
			para.ref_path=argv[i+1];
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='s')//sa_gap
		{
			para.sa_gap=atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='o')//sa_gap
		{
			para.occ_gap=atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='l')//level
		{
			para.level=atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='t')//thread_num
		{
			para.thread_num=atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='c')//max_len
		{
			c=*argv[i+1];
		}
	}
	char *dir = "./index";
	if(access(dir,0) == -1)
	{
		printf("%s not exist.\n",dir);
		mkdir(dir, 0755);
	}
	struct timeval tvs,tve;
//	build(para);
//
//	build_occA(para);
//	char index[] = "./index";

	sFMindexv sidxv;
	build_bVec(&sidxv);
	build_occAVec(&sidxv,para);
	uint32_t *occav = NULL;
	occav = (uint32_t *)read_binfile("./index/OCCAV");

	sidxv.occa = occav;
	cout << "c = " << c << endl;
	uint32_t ret, retv;
	double span;
	gettimeofday(&tvs,NULL);
	for(uint32_t i = 0; i < sidxv.bseqlen; i++)
	{
		retv = calc_OCCV(sidxv,para,c,i);
	}
	gettimeofday(&tve,NULL);
	span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"time of construct calc_OCCV is: "<<span<<endl;
	return 0;
}
