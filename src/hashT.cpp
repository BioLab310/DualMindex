/*
 * hashT.cpp
 *
 *  Created on: Feb 16, 2022
 *      Author: bio
 */

#include "hashT.h"

int cmp(uint64_t *a,uint64_t *b,uint32_t l)
{
	int label=1;
	for(uint32_t i=0;i<l;i++)
	{
		if(a[i]!=b[i])
		{
			label=0;
			break;
		}

	}
	return label;
}

uint32_t hash_function(uint64_t *a,struct bit256KmerPara para)
{
	uint64_t tmp=0;
	if(para.kmer64Len>1)
	{
		tmp=tmp|(a[para.kmer64Len-2]<<para.remainer1to64);
		tmp=tmp|a[para.kmer64Len-1];
	}
	else
	{
		tmp=a[para.kmer64Len-1];
	}
	return tmp&(hash_table_size-1);
}

void hash_table_initial(vector< vector<struct ikmer> > &a)
{
//	a.resize(hash_table_size);
	vector<struct ikmer> x;
	x.clear();
//	cout << hash_table.size() << endl;
	for(unsigned i=0;i<hash_table_size;i++)
	{
		a.push_back(x);
	}
}

void hash_table_insert(vector< vector<struct ikmer> > aa,uint64_t *node,uint64_t i,bit256KmerPara para)
{
	struct ikmer a;
	a.i=i;
	a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
	for(uint32_t x=0;x<para.kmer64Len;x++)
	{
		a.kmer[x]=node[x];
	}
	aa[hash_function(node,para)].push_back(a);
}

void hash_table_delete(vector< vector<struct ikmer> > aa,uint64_t *node,uint64_t res_index,bit256KmerPara para)
{
	struct ikmer a;
	a.i=res_index;
	a.kmer=node;

	uint32_t value;
	value=hash_function(node,para);

	if(aa[value].size()==0)
	{
		return;
	}
	else
	{
		if(aa[value].back().i==res_index)
		{
			free(aa[value].back().kmer);
			aa[value].pop_back();
		}
		else
		{
			unsigned i;
			for(i=0;i<aa[value].size();i++)
			{
				if(aa[value][i].i==res_index)
				{
					break;
				}
			}
			free(aa[value][i].kmer);
			aa[value][i]=aa[value].back();
			aa[value].pop_back();
		}
	}
}

int hash_table_query(vector< vector<struct ikmer> > aa,uint64_t *node,uint64_t *res_index,vector <string> b,bit256KmerPara para)
{
	uint32_t value;
	value=hash_function(node,para);

//	if(aa[value].size()==0)
//	{
//		return 0;
//	}
//	else
//	{
//		unsigned i;
//		for(i=0;i<aa[value].size();i++)
//		{
//			if(cmp(aa[value][i].kmer,node,para.kmer64Len)==1)
//			{
//				if(!b[aa[value][i].i].empty())
//				{
//					*res_index=aa[value][i].i;
//					return 1;
//				}
//			}
//		}
//	}
	return 0;
}


void hash_table_destroy(vector< vector<struct ikmer> > a)
{
	if(a.size()==0)
	{
		return;
	}
	else
	{
		for(uint32_t i=0;i<a.size();i++)
		{
			a[i].clear();
		}
		a.clear();
	}
}

