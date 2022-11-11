/*
 * DMindex_Oper.cpp
 *
 *  Created on: Feb 28, 2022
 *      Author: bio
 */

#include "DMindex_Oper.h"

uint32_t cal_edge_id(sFMindex &FMidx,char* e)
{
	//the edge e should be end with '\0'
	uint32_t sa_id[2];
	calc_SArangeSeq(FMidx,e,&sa_id[0],&sa_id[1]);

	if(sa_id[0]!=sa_id[1])
	{
		cout << "error: the id of an edge is wrong!" << endl;
		cout << e << endl;
		cout << sa_id[0] << endl;
		cout << sa_id[1] << endl;
	}

	return sa_id[0];
}

uint32_t cal_edge_id_new(sFMindex &FMidx,char* e,uint32_t l)
{
	//the edge e should be end with '\0'
	uint32_t sa_id[2];
	calc_SArangeSeq_new(FMidx,e,&sa_id[0],&sa_id[1],l);

	if(sa_id[0]!=sa_id[1])
	{
		cout << "error: the id of an edge is wrong!" << endl;
		cout << e << endl;
		cout << sa_id[0] << endl;
		cout << sa_id[1] << endl;
	}

	return sa_id[0];
}


void cal_SA_oneLeftShift(sFMindex &FMidx,uint32_t l,uint32_t h,char c,uint32_t *x,uint32_t *y)
{
	*x = -1;
	*y = 0;
	*x = char_to_uint8_table[c] != 0x04 ? LF_Mapping_l(FMidx,c,l) : -1;
	*y = LF_Mapping_h(FMidx,c,h);
}

uint32_t cal_Node_indegree(sFMindex &FMidx,char* n,char* ad)
{
	uint32_t p_sa_n[2];
	uint32_t p_sa_n_tmp[2];
	calc_SArangeSeq(FMidx,n,&p_sa_n[0],&p_sa_n[1]);
	if(p_sa_n[0]>p_sa_n[1])
	{
		return 0;
	}

	uint32_t ad_tmp=0;
	char xigma[4];
	xigma[0]='A';
	xigma[1]='C';
	xigma[2]='G';
	xigma[3]='T';

	char c[4];
	c[0]=1;
	c[1]=2;
	c[2]=4;
	c[3]=8;

	uint32_t nn=0;
	for(uint32_t i=0;i<4;i++)
	{
		cal_SA_oneLeftShift(FMidx,p_sa_n[0],p_sa_n[1],xigma[i],&p_sa_n_tmp[0],&p_sa_n_tmp[1]);
		if(p_sa_n_tmp[0]<=p_sa_n_tmp[1])
		{
			ad_tmp=ad_tmp+c[i];
			nn++;
		}
	}

	*ad=ad_tmp;
	return nn;
}

uint32_t cal_Node_outdegree(sFMindex &FMidx,char* n,char* ad)
{
	uint32_t p_sa_n_tmp[2];

	calc_SArangeSeq(FMidx,n,&p_sa_n_tmp[0],&p_sa_n_tmp[1]);
	if(p_sa_n_tmp[0]>p_sa_n_tmp[1])
	{
		cerr<<"error: from cal_Node_outdegree()!" << endl;
		return 0;
	}

	uint32_t L=strlen(n);
	L=L+2;
	char* p_n_tmp;
	p_n_tmp=(char*)malloc(sizeof(char)*L);
	for(uint32_t i=0;i<L-2;i++)
	{
		p_n_tmp[i]=n[i];
	}
	p_n_tmp[L-1]='\0';

	uint32_t ad_tmp=0;
	char xigma[4];
	xigma[0]='A';
	xigma[1]='C';
	xigma[2]='G';
	xigma[3]='T';

	char c[4];
	c[0]=1;
	c[1]=2;
	c[2]=4;
	c[3]=8;

	uint32_t nn=0;
	for(uint32_t i=0;i<4;i++)
	{
		p_n_tmp[L-2]=xigma[i];
		calc_SArangeSeq(FMidx,p_n_tmp,&p_sa_n_tmp[0],&p_sa_n_tmp[1]);
		if(p_sa_n_tmp[0]==p_sa_n_tmp[1])
		{
			ad_tmp=ad_tmp+c[i];
			nn++;
		}
	}

	*ad=ad_tmp;
	free(p_n_tmp);
	return nn;
}

uint32_t is_Node_or_Edge_the_member_of_graph(sFMindex &FMidx,char* n,uint32_t* l,uint32_t* h)
{
	uint32_t p_sa_n[2];
	calc_SArangeSeq(FMidx,n,&p_sa_n[0],&p_sa_n[1]);
	if(p_sa_n[0]<=p_sa_n[1])
	{
		*l=p_sa_n[0];
		*h=p_sa_n[1];
		return 1;
	}
	else
	{
		return 0;
	}
}

