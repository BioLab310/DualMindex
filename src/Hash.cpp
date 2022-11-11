/*
 * Hash.cpp
 *
 *  Created on: Jan 24, 2019
 *      Author: bio
 */
#include "Hash.h"

uint32_t bit256hashFFunction(uint64_t *a,struct bit256KmerPara para)
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
	return tmp&(0xFFFFFF);
}
struct NodeBit** bit256initialHashFTable()
{
	struct NodeBit** pk;
	pk=(struct NodeBit**)malloc(sizeof(struct NodeBit*)*HashFSize);
	for(uint32_t i=0;i<HashFSize;i++)
	{
		pk[i]=NULL;
	}
	return pk;
}
struct NodeBit * bit256initialSingleHashFTable()
{
	struct NodeBit* p;
	p=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p);
	return p;
}
void bit256insertHashFTable(struct NodeBit** pk,struct nodeBit c_tmp_hashtable,struct bit256KmerPara para)
{
	uint32_t tmp=bit256hashFFunction(c_tmp_hashtable.hashValue,para);
	if(pk[tmp]==NULL)
	{
		pk[tmp]=bit256initialSingleHashFTable();
	}
	Insert_Value_bit(pk+tmp,c_tmp_hashtable,para);
}
int64_t getHashFTableValue(struct NodeBit** pk,uint64_t *a,struct bit256KmerPara para)
{
	uint32_t tmp=bit256hashFFunction(a,para);
	if(pk[tmp]!=NULL)
	{
		return MappingHashValueToID_bit(pk[tmp],a,para);
	}
	return -1;
}
void bit256freeHashFTable(struct NodeBit **ph,struct bit256KmerPara para)
{
	for(uint32_t i=0;i<HashFSize;i++)
	{
		if(ph[i]!=NULL)
		{
			destory_tree_bit(ph[i],para);
		}
	}
	free(ph);
}

/***************************hashFor256BitCharKmer******************************/
uint32_t bit256HashFunction(uint64_t* a,struct bit256KmerPara para)
{
	uint64_t tmp=0;
	if(para.kmer64Len>1&&para.remainer1to64!=0)
	{
		tmp=a[para.kmer64Len-2]<<para.remainer1to64;
		tmp=tmp|a[para.kmer64Len-1];
	}
	else
	{
		tmp=a[para.kmer64Len-1];
	}
	return tmp&(HashSize-1);
}
uint32_t bit256HashFunction_left(uint64_t* a,struct bit256KmerPara para)
{
	uint64_t tmp=a[0];
	return tmp&(0xFFFFFFF);
}

struct bit256Hash* initial256BitHashTable()
{
	struct bit256Hash *ph;
	ph=(struct bit256Hash*)malloc(sizeof(struct bit256Hash)*1);
	uint64_t** pk;
	pk=(uint64_t**)malloc(sizeof(uint64_t*)*HashSize);
	uint64_t** ppo;
	ppo=(uint64_t**)malloc(sizeof(uint64_t*)*HashSize);

	uint64_t *pp;
	pp=(uint64_t *)malloc(sizeof(uint64_t)*HashSize);
	memset(pp,0,sizeof(uint64_t)*HashSize);
	uint64_t *pmp;
	pmp=(uint64_t *)malloc(sizeof(uint64_t)*HashSize);
	memset(pmp,0,sizeof(uint64_t)*HashSize);
	ph->p_kmer=pk;
	ph->p_pos=ppo;
	ph->len=pp;
	ph->maxlen=pmp;

	return ph;
}
uint64_t PosSearch256BitTable(uint64_t* p,uint64_t len,uint64_t* y,bit256KmerPara para)
{
	if(len==0)
	{
		return 0;
	}
	uint32_t wid=para.kmer64Len;
	uint64_t s=0;
	uint64_t e=len-1;
	uint64_t m=0;
	uint64_t r=0;
	if(cmp256BitKmer(y,p,wid)==0||cmp256BitKmer(y,p,wid)==2)
	{
		r=0;
	}
	else if(cmp256BitKmer((p+wid*(len-1)),y,wid)==0||cmp256BitKmer((p+wid*(len-1)),y,wid)==2)
	{
		r=len;
	}
	else
	{
		if(len<MinSizeForSearch)
		{
			for(uint64_t i=1;i<len;i++)
			{
				uint32_t x=cmp256BitKmer(y,p+wid*i,wid);
				if(x==0||x==2)
				{
					r=i;
					break;
				}
			}
		}
		else
		{
			while(e-s>1)
			{
				m=(s+e)/2;
				uint32_t x=cmp256BitKmer(y,p+wid*m,wid);
				if(x==0)
				{
					e=m;
				}
				else if(x==1)
				{
					s=m;
				}
				else
				{
					r=m;
					return r;
				}
			}
			r=e;
		}
	}
	return r;
}
void mvtoNext(uint64_t * p, uint64_t * pp,uint64_t pos,uint64_t len, uint32_t wid)
{
	uint64_t fix=len-1+pos;
	for(uint64_t i=pos;i<len;i++)
	{
		uint64_t j=fix-i;//len-1-(i - pos);
		uint64_t l1=(j+1)*wid;
		uint64_t l2=j*wid;
		for(uint64_t k=0;k<wid;k++)
		{
			p[l1+k]=p[l2+k];
		}
		pp[j+1]=pp[j];
	}
}

void insert256BitHashTable(struct bit256Hash *ph,uint64_t* kmer,uint64_t pos,bit256KmerPara para)
{
	uint32_t shift=bit256HashFunction(kmer,para);
	if(shift>HashSize)
	{
		cout << "error!" << endl;
	}
	if(ph->maxlen[shift]==0)
	{
		ph->p_kmer[shift]=(uint64_t*)malloc(sizeof(uint64_t)*\
				(ph->maxlen[shift]+AddSize)*para.kmer64Len);
		ph->p_pos[shift]=(uint64_t*)malloc(sizeof(uint64_t)*\
				(ph->maxlen[shift]+AddSize));
		ph->maxlen[shift]=ph->maxlen[shift]+AddSize;

		uint64_t current=PosSearch256BitTable(ph->p_kmer[shift],ph->len[shift],kmer,para);

		if(current!=ph->len[shift])
		{
			mvtoNext(ph->p_kmer[shift], ph->p_pos[shift],current,ph->len[shift],para.kmer64Len);
		}
		uint64_t x=current*para.kmer64Len;
		for(uint32_t i=0;i<para.kmer64Len;i++)
		{
			ph->p_kmer[shift][x+i]=kmer[i];
		}
		ph->p_pos[shift][current]=pos;

		ph->len[shift]++;
	}
	else if(ph->len[shift]>=ph->maxlen[shift])
	{
		ph->p_kmer[shift]=(uint64_t*)realloc(ph->p_kmer[shift],\
				sizeof(uint64_t)*(ph->maxlen[shift]+AddSize)*para.kmer64Len);
		ph->p_pos[shift]=(uint64_t*)realloc(ph->p_pos[shift],\
				sizeof(uint64_t)*(ph->maxlen[shift]+AddSize));
		ph->maxlen[shift]=ph->maxlen[shift]+AddSize;

		uint64_t current=PosSearch256BitTable(ph->p_kmer[shift],ph->len[shift],kmer,para);

		if(current!=ph->len[shift])
		{
			mvtoNext(ph->p_kmer[shift], ph->p_pos[shift],current,ph->len[shift],para.kmer64Len);
		}
		uint64_t x=current*para.kmer64Len;
		for(uint32_t i=0;i<para.kmer64Len;i++)
		{
			ph->p_kmer[shift][x+i]=kmer[i];
		}
		ph->p_pos[shift][current]=pos;

		ph->len[shift]++;
	}
	else
	{
		uint64_t current=PosSearch256BitTable(ph->p_kmer[shift],ph->len[shift],kmer,para);

		if(current!=ph->len[shift])
		{
			mvtoNext(ph->p_kmer[shift], ph->p_pos[shift],current,ph->len[shift],para.kmer64Len);
		}
		uint64_t x=current*para.kmer64Len;
		for(uint32_t i=0;i<para.kmer64Len;i++)
		{
			ph->p_kmer[shift][x+i]=kmer[i];
		}
		ph->p_pos[shift][current]=pos;

		ph->len[shift]++;
	}
}

void insert256BitHashTable_parallel(struct bit256Hash *ph,\
		uint64_t* kmer,\
		uint64_t pos,\
		bit256KmerPara para,\
		uint32_t thread_num,\
		uint32_t thread_id)
{
	uint32_t shift=bit256HashFunction(kmer,para);
	if(shift>0x10000000)
	{
		cout << "error!" << endl;
	}
	if((shift%thread_num)==thread_id)
	{
		if(ph->maxlen[shift]==0)
		{
			ph->p_kmer[shift]=(uint64_t*)malloc(sizeof(uint64_t)*\
					(ph->maxlen[shift]+AddSize)*para.kmer64Len);
			ph->p_pos[shift]=(uint64_t*)malloc(sizeof(uint64_t)*\
					(ph->maxlen[shift]+AddSize));
			ph->maxlen[shift]=ph->maxlen[shift]+AddSize;

			uint64_t current=PosSearch256BitTable(ph->p_kmer[shift],ph->len[shift],kmer,para);

			if(current!=ph->len[shift])
			{
				mvtoNext(ph->p_kmer[shift], ph->p_pos[shift],current,ph->len[shift],para.kmer64Len);
			}
			uint64_t x=current*para.kmer64Len;
			for(uint32_t i=0;i<para.kmer64Len;i++)
			{
				ph->p_kmer[shift][x+i]=kmer[i];
			}
			ph->p_pos[shift][current]=pos;

			ph->len[shift]++;
		}
		else if(ph->len[shift]>=ph->maxlen[shift])
		{
			ph->p_kmer[shift]=(uint64_t*)realloc(ph->p_kmer[shift],\
					sizeof(uint64_t)*(ph->maxlen[shift]+AddSize)*para.kmer64Len);
			ph->p_pos[shift]=(uint64_t*)realloc(ph->p_pos[shift],\
					sizeof(uint64_t)*(ph->maxlen[shift]+AddSize));
			ph->maxlen[shift]=ph->maxlen[shift]+AddSize;

			uint64_t current=PosSearch256BitTable(ph->p_kmer[shift],ph->len[shift],kmer,para);

			if(current!=ph->len[shift])
			{
				mvtoNext(ph->p_kmer[shift], ph->p_pos[shift],current,ph->len[shift],para.kmer64Len);
			}
			uint64_t x=current*para.kmer64Len;
			for(uint32_t i=0;i<para.kmer64Len;i++)
			{
				ph->p_kmer[shift][x+i]=kmer[i];
			}
			ph->p_pos[shift][current]=pos;

			ph->len[shift]++;
		}
		else
		{
			uint64_t current=PosSearch256BitTable(ph->p_kmer[shift],ph->len[shift],kmer,para);

			if(current!=ph->len[shift])
			{
				mvtoNext(ph->p_kmer[shift], ph->p_pos[shift],current,ph->len[shift],para.kmer64Len);
			}
			uint64_t x=current*para.kmer64Len;
			for(uint32_t i=0;i<para.kmer64Len;i++)
			{
				ph->p_kmer[shift][x+i]=kmer[i];
			}
			ph->p_pos[shift][current]=pos;

			ph->len[shift]++;
		}

	}
}
uint64_t binary256BitSearch(uint64_t * ppo,uint64_t* p,uint64_t len,uint64_t* y,bit256KmerPara para)
{
	if(len==0)
	{
		return 0;
	}
	uint32_t wid=para.kmer64Len;
	uint64_t s=0;
	uint64_t e=len-1;
	uint64_t m;
	if(cmp256BitKmer(p,y,wid)==2)
	{
		return ppo[0];
	}
	else if(cmp256BitKmer((p+wid*(len-1)),y,wid)==2)
	{
		return ppo[len-1];
	}
	else
	{
		while(e-s>1)
		{
			m=(s+e)/2;
			if(cmp256BitKmer(p+wid*m,y,wid)==2)
			{
				return ppo[m];
			}
			else
			{
				if(cmp256BitKmer(p+wid*m,y,wid)==1)
				{
					e=m;
				}
				else
				{
					s=m;
				}
			}
		}
		return 0;
	}
}
uint64_t search256BitHashTable(struct bit256Hash *ph,uint64_t* y,bit256KmerPara para)
{
	uint32_t wid=para.kmer64Len;
	uint32_t shift=bit256HashFunction(y,para);
	uint32_t label=1;
	if(ph->len[shift]<MinSizeForSearch)
	{
		for(uint32_t i=0;i<ph->len[shift];i++)
		{
			if(cmp256BitKmer(ph->p_kmer[shift]+wid*i,y,wid)==2)
			{
				return ph->p_pos[shift][i];
			}
		}
	}
	else
	{
		return binary256BitSearch(ph->p_pos[shift],ph->p_kmer[shift],ph->len[shift],y,para);
	}

	if(label)
	{
		return 0;
	}
}
void freeHash256BitTable(struct bit256Hash *ph)
{
	if(ph!=NULL)
	{
		for(uint32_t i=0;i<HashSize;i++)
		{
			if(ph->maxlen[i]!=0&&ph->p_kmer[i]!=NULL&&ph->p_pos[i]!=NULL)
			{
				free(ph->p_kmer[i]);
				free(ph->p_pos[i]);
			}
		}
		if(ph->p_kmer!=NULL)
		{
			free(ph->p_kmer);
		}
		if(ph->p_pos!=NULL)
		{
			free(ph->p_pos);
		}
		if(ph->len!=NULL)
		{
			free(ph->len);
		}
		if(ph->maxlen!=NULL)
		{
			free(ph->maxlen);
		}
		free(ph);

	}
}
/***************************hashFor256BitCharKmer******************************/

