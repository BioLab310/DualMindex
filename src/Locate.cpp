/*
 * Locate.cpp
 *
 *  Created on: Mar 28, 2022
 *      Author: bio
 */

#include "Locate.h"


void left_expand_end(struct DMindex &pDM,char *e,struct shift_edge* p_r,int32_t shift_value)
{
	char* p;
	p=(char*)malloc(sizeof(char)*(pDM.kmer_len+2));
	p[pDM.kmer_len+1]='\0';
	p[pDM.kmer_len]='\0';
	for(uint32_t i=0;i<pDM.kmer_len;i++)
	{
		p[i]=e[i];
	}
	char c;
	c=e[pDM.kmer_len];

	uint32_t l=0;
	uint32_t h=0;

	calc_SArangeSeq(pDM.FMidx,p,&l,&h);

	char xigma[4];
	xigma[0]='A';
	xigma[1]='C';
	xigma[2]='G';
	xigma[3]='T';

	p_r->shift_value=shift_value;
	while(l==h)
	{
		uint32_t l1=1;
		uint32_t h1=0;
		for(uint32_t k=0;k<4;k++)
		{
			cal_SA_oneLeftShift(pDM.FMidx,l,h,xigma[k],&l1,&h1);
			if(l1<=h1)
			{
				c=p[pDM.kmer_len-1];
				for(uint32_t i=0;i<pDM.kmer_len-1;i++)
				{
					p[pDM.kmer_len-1-i]=p[pDM.kmer_len-1-i-1];
				}
				p[0]=xigma[k];
				break;
			}
		}
		calc_SArangeSeq(pDM.FMidx,p,&l,&h);
		p_r->shift_value++;
	}
	p[pDM.kmer_len]=c;
	p[pDM.kmer_len+1]='\0';
	p_r->e=p;
}

void right_expand_end(struct DMindex &pDM,char *e,struct shift_edge* p_r,int32_t shift_value)
{
	char* p;
	p=(char*)malloc(sizeof(char)*(pDM.kmer_len+2));
	p[pDM.kmer_len+1]='\0';
	p[pDM.kmer_len]='\0';
	for(uint32_t i=0;i<pDM.kmer_len;i++)
	{
		p[i]=e[i+1];
	}
	char c;
	c=e[0];

	uint32_t l=0;
	uint32_t h=0;

	calc_SArangeSeq(pDM.FMidx,p,&l,&h);

	char xigma[4];
	xigma[0]='A';
	xigma[1]='C';
	xigma[2]='G';
	xigma[3]='T';

	p_r->shift_value=shift_value;
	while(l==h)
	{
		uint32_t l1=1;
		uint32_t h1=0;
		for(uint32_t k=0;k<4;k++)
		{
			p[pDM.kmer_len]=xigma[k];
			calc_SArangeSeq(pDM.FMidx,p,&l1,&h1);
			if(l1<=h1)
			{
				c=p[0];
				for(uint32_t i=0;i<=pDM.kmer_len-1;i++)
				{
					p[i]=p[i+1];
				}
				p[pDM.kmer_len]='\0';
				break;
			}
		}
		calc_SArangeSeq(pDM.FMidx,p,&l,&h);
		p_r->shift_value--;
	}
	for(uint32_t i=0;i<=pDM.kmer_len-1;i++)
	{
		p[pDM.kmer_len-i]=p[pDM.kmer_len-i-1];
	}
	p[0]=c;
	p[pDM.kmer_len+1]='\0';
	p_r->e=p;
}

int32_t Locate_edge(struct DMindex &pDM,char *e,uint32_t label,struct shift_miniID** p_r,uint32_t *p_r_len,int32_t shift_value)
{
	//label==0,no direction
	//label==1,left
	//label==2,right
//	cout << "label:" << label << endl;
//	cout << "e:" << e << endl;
//	cout << "shift_value:" << shift_value << endl;
	if(label==0)
	{
		uint32_t id;
		uint32_t sa_id[2];
		calc_SArangeSeq(pDM.FMidx,e,&sa_id[0],&sa_id[1]);

		if(sa_id[0]!=sa_id[1])
		{
			cout << "error: this edge is not a edge of the graph!" << endl;
			return -1;
		}
		else
		{
			id=sa_id[0];
			if(is_minimal_edge(pDM.p_bitVec,id))
			{
				//calculate the mini_id
				int64_t miniID;
				miniID=get_minimal_edge_id(pDM.p_bitVec,id,pDM.v_bit_gap,pDM.p_Nb);
				//return the interval
				//(pDM.p_Ps,pDM.p_Pi[miniID],pDM.p_Pi[miniID+1]-1)
				(*p_r)[*p_r_len].mini_edge_id=miniID;
				(*p_r)[*p_r_len].shift_value=shift_value;
				(*p_r_len)++;
				(*p_r)=(struct shift_miniID*)realloc((*p_r),sizeof(struct shift_miniID)*(*p_r_len+1));
			}
			else
			{
				struct shift_edge p_e;
				left_expand_end(pDM,e,&p_e,shift_value);
				id=cal_edge_id(pDM.FMidx,p_e.e);
				if(is_minimal_edge(pDM.p_bitVec,id))
				{
					int64_t miniID;
					miniID=get_minimal_edge_id(pDM.p_bitVec,id,pDM.v_bit_gap,pDM.p_Nb);
					//return the interval
					//(pDM.p_Ps,pDM.p_Pi[miniID],pDM.p_Pi[miniID+1]-1)
					(*p_r)[*p_r_len].mini_edge_id=miniID;
					(*p_r)[*p_r_len].shift_value=p_e.shift_value;
					(*p_r_len)++;
					*p_r=(struct shift_miniID*)realloc((*p_r),sizeof(struct shift_miniID)*(*p_r_len+1));
					free(p_e.e);
				}
				else
				{
					char c_tmp;
					c_tmp=p_e.e[pDM.kmer_len];
					p_e.e[pDM.kmer_len]='\0';

					char ad_tmp;
					uint32_t out_degree_tmp=cal_Node_outdegree(pDM.FMidx,p_e.e,&ad_tmp);
					if(out_degree_tmp==1)
					{
						p_e.e[pDM.kmer_len]=c_tmp;
						Locate_edge(pDM,p_e.e,1,p_r,p_r_len,p_e.shift_value);
						free(p_e.e);
					}
					else
					{
						free(p_e.e);
						right_expand_end(pDM,e,&p_e,shift_value);
						Locate_edge(pDM,p_e.e,2,p_r,p_r_len,p_e.shift_value);
						free(p_e.e);
					}
				}
			}
		}
	}
	else if(label==1)
	{
		uint32_t id;
		char c_tmp=e[pDM.kmer_len];
		e[pDM.kmer_len]='\0';
		char ad_tmp;
		uint32_t in_degree_tmp=cal_Node_indegree(pDM.FMidx,e,&ad_tmp);
//		cout << in_degree_tmp << endl;
//		cout << (uint32_t)ad_tmp << endl;

		if(in_degree_tmp<=1)
		{
			cout << "error: from edge locating, step 1:left expand!" << endl;
			cout << "in_degree:" << in_degree_tmp << endl;
			cout << "ad_tmp:" << (uint32_t)ad_tmp << endl;
		}
		else
		{
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

			for(uint32_t j=0;j<pDM.kmer_len;j++)
			{
				e[pDM.kmer_len-j]=e[pDM.kmer_len-j-1];
			}

			for(uint32_t i=0;i<4;i++)
			{
				if((ad_tmp&c[i])==c[i])
				{
					e[0]=xigma[i];

					struct shift_edge p_e;
					left_expand_end(pDM,e,&p_e,shift_value+1);

					id=cal_edge_id(pDM.FMidx,p_e.e);
					if(is_minimal_edge(pDM.p_bitVec,id))
					{
//						cout << p_e.e << endl;
						int64_t miniID;
						miniID=get_minimal_edge_id(pDM.p_bitVec,id,pDM.v_bit_gap,pDM.p_Nb);
						//return the interval
						//(pDM.p_Ps,pDM.p_Pi[miniID],pDM.p_Pi[miniID+1]-1)
						(*p_r)[*p_r_len].mini_edge_id=miniID;
						(*p_r)[*p_r_len].shift_value=p_e.shift_value;
						(*p_r_len)++;
						*p_r=(struct shift_miniID*)realloc((*p_r),sizeof(struct shift_miniID)*(*p_r_len+1));
						free(p_e.e);
					}
					else
					{
						Locate_edge(pDM,p_e.e,1,p_r,p_r_len,p_e.shift_value);
						free(p_e.e);
					}
				}
			}
		}

	}
	else if(label==2)
	{
		uint32_t id;
		char c_tmp=e[0];
		for(uint32_t j=0;j<pDM.kmer_len;j++)
		{
			e[j]=e[j+1];
		}
		e[pDM.kmer_len]='\0';
		char ad_tmp;
		uint32_t out_degree_tmp=cal_Node_outdegree(pDM.FMidx,e,&ad_tmp);

		if(out_degree_tmp<=1)
		{
			cout << "error: from edge locating, step 2:right expand!" << endl;
		}
		else
		{
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


			for(uint32_t i=0;i<4;i++)
			{
				if((ad_tmp&c[i])==c[i])
				{
					e[pDM.kmer_len]=xigma[i];

					id=cal_edge_id(pDM.FMidx,e);
					if(is_minimal_edge(pDM.p_bitVec,id))
					{
						int64_t miniID;
						miniID=get_minimal_edge_id(pDM.p_bitVec,id,pDM.v_bit_gap,pDM.p_Nb);
						//return the interval
						//(pDM.p_Ps,pDM.p_Pi[miniID],pDM.p_Pi[miniID+1]-1)
						(*p_r)[*p_r_len].mini_edge_id=miniID;
						(*p_r)[*p_r_len].shift_value=shift_value-1;
						(*p_r_len)++;
						*p_r=(struct shift_miniID*)realloc((*p_r),sizeof(struct shift_miniID)*(*p_r_len+1));
					}
					else
					{
						struct shift_edge p_e;
						right_expand_end(pDM,e,&p_e,shift_value-1);
						Locate_edge(pDM,p_e.e,2,p_r,p_r_len,p_e.shift_value);
						free(p_e.e);
					}
				}
			}
		}
	}
	return 0;
}

int32_t Locate_edge_whole(struct DMindex &pDM,char *e,struct position_list *p)
{
//	uint32_t sa_id[2];
//	char pp[]="CGGCAATATGTCT";
//
//	char ad_tmp;
//	uint32_t out_degree_tmp=cal_Node_outdegree(pDM.FMidx,pp,&ad_tmp);
//	cout << out_degree_tmp << ":" << (uint32_t)ad_tmp << endl;
//	uint32_t in_degree_tmp=cal_Node_indegree(pDM.FMidx,pp,&ad_tmp);
//	cout << in_degree_tmp << ":" << (uint32_t)ad_tmp << endl;
//
//	calc_SArangeSeq(pDM.FMidx,pp,&sa_id[0],&sa_id[1]);
//	cout << "said:" << endl;
//	cout << sa_id[0] << ":" << sa_id[1] << endl;

	uint32_t label=0;
	struct shift_miniID* p_r;
	uint32_t p_r_len=0;
	int32_t shift_value=0;
	p_r=(struct shift_miniID*)malloc(sizeof(struct shift_miniID)*1);

	int32_t r=Locate_edge(pDM,e,label,&p_r,&p_r_len,shift_value);
	cout << "r:" << r<< endl;
	cout << "the length of p_r is :" << p_r_len << endl;
	if(r==-1)
	{
		return -1;
	}
	else
	{
		vector<struct position_list> v_tmp;
		struct position_list tmp;
		for(uint32_t i=0;i<p_r_len;i++)
		{
			tmp.p_pos_list=pDM.p_Ps+pDM.p_Pi[p_r[i].mini_edge_id];
			tmp.p_pos_list_len=pDM.p_Pi[p_r[i].mini_edge_id+1]-pDM.p_Pi[p_r[i].mini_edge_id];
			tmp.offset=p_r[i].shift_value;
			v_tmp.push_back(tmp);
		}
		combn_pos_list_multi(v_tmp, p);
	}
	free(p_r);
	return 0;
}

int32_t Locate_path(struct DMindex &pDM,char *h,struct position_list *r)
{
	//-1,the position list is empty.
//	struct shift_miniID* p_r;
//	p_r=(struct shift_miniID*)malloc(sizeof(struct shift_miniID)*1);
//	uint32_t p_r_len=0;

	struct position_list tmp;
	vector<struct position_list> v_tmp;
	v_tmp.clear();
	uint32_t n_tmp;

	uint32_t h_len=strlen(h);
	int32_t min_pos;
	int32_t max_pos;

	min_pos=h_len+1;
	max_pos=-1;

	for(int32_t i=0;i<=h_len-pDM.kmer_len;i++)
	{
		char c;
		c=h[i+pDM.kmer_len];
		h[i+pDM.kmer_len]='\0';

		uint32_t p_sa_n[2];
		calc_SArangeSeq(pDM.FMidx,h+i,&p_sa_n[0],&p_sa_n[1]);
		if(p_sa_n[0]<p_sa_n[1])
		{
			if(min_pos>i)
			{
				min_pos=i;
			}
			if(max_pos<i)
			{
				max_pos=i;
			}
		}
		else if(p_sa_n[0]>p_sa_n[1])
		{
			return -1;
		}
		h[i+pDM.kmer_len]=c;
	}

	cout <<"min:"<< min_pos << endl;
	cout <<"max:"<< max_pos << endl;

	if(min_pos==h_len+1)
	{
		char c;
		c=h[pDM.kmer_len+1];
		h[1+pDM.kmer_len]='\0';

		uint32_t p_sa_n[2];
		calc_SArangeSeq(pDM.FMidx,h,&p_sa_n[0],&p_sa_n[1]);
		if(p_sa_n[0]==p_sa_n[1])
		{
			Locate_edge_whole(pDM,h,&tmp);
//			tmp.offset=tmp.offset-0;
			v_tmp.push_back(tmp);
		}
		else
		{
			return -1;
		}
		h[pDM.kmer_len+1]=c;
	}
	else if(min_pos==max_pos)
	{
		if(min_pos==0)
		{
			char c;
			c=h[pDM.kmer_len+1];
			h[1+pDM.kmer_len]='\0';

			uint32_t p_sa_n[2];
			calc_SArangeSeq(pDM.FMidx,h,&p_sa_n[0],&p_sa_n[1]);
			if(p_sa_n[0]==p_sa_n[1])
			{
//				Locate_edge(pDM,h,0,p_r,p_r_len,0);
				Locate_edge_whole(pDM,h,&tmp);
				v_tmp.push_back(tmp);

			}
			else
			{
				return -1;
			}
			h[pDM.kmer_len+1]=c;
		}
		else if(min_pos==h_len-pDM.kmer_len)
		{
			char c;
			c=h[h_len-pDM.kmer_len-1+pDM.kmer_len+1];
			h[h_len-pDM.kmer_len-1+pDM.kmer_len+1]='\0';

			uint32_t p_sa_n[2];
			calc_SArangeSeq(pDM.FMidx,h+h_len-pDM.kmer_len-1,&p_sa_n[0],&p_sa_n[1]);
			h[h_len-pDM.kmer_len-1+pDM.kmer_len+1]=c;
			if(p_sa_n[0]==p_sa_n[1])
			{
//				Locate_edge(pDM,h+h_len-pDM.kmer_len-1,0,p_r,p_r_len,0);
				c=h[pDM.kmer_len+1];
				h[1+pDM.kmer_len]='\0';

				Locate_edge_whole(pDM,h,&tmp);
				v_tmp.push_back(tmp);

				h[pDM.kmer_len+1]=c;
			}
			else
			{
				return -1;
			}
		}
		else
		{
			char c;
			c=h[min_pos+pDM.kmer_len];
			h[min_pos+pDM.kmer_len]='\0';

			uint32_t out_degree;
			uint32_t in_degree;
			char ad;

			in_degree=cal_Node_indegree(pDM.FMidx,h+min_pos,&ad);
			out_degree=cal_Node_outdegree(pDM.FMidx,h+min_pos,&ad);

			h[min_pos+pDM.kmer_len]=c;

			if(out_degree>1)
			{
				char c;
				c=h[min_pos+pDM.kmer_len+1];
				h[min_pos+1+pDM.kmer_len]='\0';

				uint32_t p_sa_n[2];
				calc_SArangeSeq(pDM.FMidx,h+min_pos,&p_sa_n[0],&p_sa_n[1]);
				if(p_sa_n[0]==p_sa_n[1])
				{
//					Locate_edge(pDM,h+min_pos,0,p_r,p_r_len,0);
					Locate_edge_whole(pDM,h+min_pos,&tmp);
					tmp.offset=tmp.offset-min_pos;
					v_tmp.push_back(tmp);
				}
				else
				{
					return -1;
				}
				h[pDM.kmer_len+1+min_pos]=c;
			}
			if(in_degree>1)
			{
				char c;
				c=h[min_pos-1+pDM.kmer_len+1];
				h[min_pos-1+1+pDM.kmer_len]='\0';

				uint32_t p_sa_n[2];
				calc_SArangeSeq(pDM.FMidx,h+min_pos-1,&p_sa_n[0],&p_sa_n[1]);
				h[min_pos-1+1+pDM.kmer_len]=c;
				if(p_sa_n[0]==p_sa_n[1])
				{
					c=h[pDM.kmer_len+1];
					h[1+pDM.kmer_len]='\0';
//					Locate_edge(pDM,h+min_pos-1,0,p_r,p_r_len,0);
					Locate_edge_whole(pDM,h,&tmp);
					v_tmp.push_back(tmp);

					h[1+pDM.kmer_len]=c;
				}
				else
				{
					return -1;
				}
			}

		}
	}
	else
	{
		char ad;

		if(min_pos>0)
		{
			char c;
			c=h[min_pos+pDM.kmer_len];
			h[min_pos+pDM.kmer_len]='\0';
			uint32_t in_degree;
			in_degree=cal_Node_indegree(pDM.FMidx,h+min_pos,&ad);
			h[min_pos+pDM.kmer_len]=c;

			if(in_degree>1)
			{
				char c;
				c=h[min_pos-1+pDM.kmer_len+1];
				h[min_pos-1+1+pDM.kmer_len]='\0';

				uint32_t p_sa_n[2];
				calc_SArangeSeq(pDM.FMidx,h+min_pos-1,&p_sa_n[0],&p_sa_n[1]);
				h[min_pos-1+1+pDM.kmer_len]=c;
				if(p_sa_n[0]==p_sa_n[1])
				{
//					Locate_edge(pDM,h+min_pos-1,0,p_r,p_r_len,0);
					c=h[min_pos-1+pDM.kmer_len+1];
					h[min_pos-1+1+pDM.kmer_len]='\0';
					Locate_edge_whole(pDM,h+min_pos-1,&tmp);
					tmp.offset=tmp.offset-(min_pos-1);
					v_tmp.push_back(tmp);

					h[min_pos-1+1+pDM.kmer_len]=c;
				}
				else
				{
					return -1;
				}
			}
		}

		if(max_pos<h_len-pDM.kmer_len)
		{
			char c1;
			c1=h[max_pos+pDM.kmer_len];
			h[max_pos+pDM.kmer_len]='\0';
			uint32_t out_degree;
			out_degree=cal_Node_outdegree(pDM.FMidx,h+max_pos,&ad);
			cout <<"out:"<< out_degree << endl;
			h[max_pos+pDM.kmer_len]=c1;

			if(out_degree>1)
			{
				char c;
				c=h[max_pos+pDM.kmer_len+1];
				h[max_pos+1+pDM.kmer_len]='\0';

				uint32_t p_sa_n[2];
				calc_SArangeSeq(pDM.FMidx,h+max_pos,&p_sa_n[0],&p_sa_n[1]);
				if(p_sa_n[0]==p_sa_n[1])
				{
//					Locate_edge(pDM,h+max_pos,0,p_r,p_r_len,0);
					Locate_edge_whole(pDM,h+max_pos,&tmp);
					tmp.offset=tmp.offset-max_pos;
					v_tmp.push_back(tmp);

					h[pDM.kmer_len+1+max_pos]=c;
				}
				else
				{
					return -1;
				}
			}
		}

		n_tmp=v_tmp.size();
		cout << "n_tmp:" << n_tmp << endl;

		for(uint32_t i=min_pos;i<max_pos;i++)
		{
			char c;
			c=h[i+pDM.kmer_len+1];
			h[i+pDM.kmer_len+1]='\0';

			uint32_t sa_id[2];
			calc_SArangeSeq(pDM.FMidx,h+i,&sa_id[0],&sa_id[1]);

			if(sa_id[0]!=sa_id[1])
			{
				return -1;
			}
			else
			{
				int64_t miniID;
				miniID=get_minimal_edge_id(pDM.p_bitVec,sa_id[0],pDM.v_bit_gap,pDM.p_Nb);
				if(miniID!=-1)
				{
					tmp.p_pos_list=pDM.p_Ps+pDM.p_Pi[miniID];
					tmp.p_pos_list_len=pDM.p_Pi[miniID+1]-pDM.p_Pi[miniID];
					tmp.offset=-i;
					v_tmp.push_back(tmp);
					cout << "miniID:"<< miniID<< endl;
					cout << "tmp.p_pos_list_len:" <<tmp.p_pos_list_len<<endl;
					cout << "i:"<< i<< endl;


//					(*p_r)[*p_r_len].mini_edge_id=miniID;
//					(*p_r)[*p_r_len].shift_value=-i;
//					(*p_r_len)++;
//					(*p_r)=(struct shift_miniID*)realloc((*p_r),sizeof(struct shift_miniID)*(*p_r_len+1));
				}
			}
			h[i+pDM.kmer_len+1]=c;
		}
	}

	int32_t x=0;
	cout << "v_tmp:" << v_tmp.size() << endl;
	x=mergy_pos_list_multi(v_tmp, r);
	for(uint32_t i=0;i<n_tmp;i++)
	{
		free(v_tmp[i].p_pos_list);
	}
	return x;
}

bool cmp(struct position_list a,struct position_list b)
{
    return a.p_pos_list_len<b.p_pos_list_len;
}
void mergy_pos_list(struct position_list p1,struct position_list p2,struct position_list *p3)
{
	struct pos* r;
	struct pos* q1;
	struct pos* q2;
    uint32_t s1,lr;
    uint32_t I1,I2;
    int32_t O1,O2;

    if(p1.p_pos_list_len>=p2.p_pos_list_len)
    {
        q1=p1.p_pos_list;
        I1=p1.p_pos_list_len;
        O1=p1.offset;
        q2=p2.p_pos_list;
        I2=p2.p_pos_list_len;
        O2=p2.offset;
        s1=p1.p_pos_list_len/p2.p_pos_list_len;
    }
    else
    {
        q2=p1.p_pos_list;
        I2=p1.p_pos_list_len;
        O2=p1.offset;
        q1=p2.p_pos_list;
        I1=p2.p_pos_list_len;
        O1=p2.offset;
        s1=p2.p_pos_list_len/p1.p_pos_list_len;
    }
	r=(struct pos*)malloc(sizeof(struct pos)*I2);
    uint32_t i1,i2;
    i1=0;
    i2=0;
    lr=0;
    while(i1<I1&&i2<I2)
    {
        if(q1[i1].ref_id>q2[i2].ref_id)
        {
        	i2++;
        }
        else if(q1[i1].ref_id<q2[i2].ref_id)
        {
        	i1++;
        }
        else
        {
        	int64_t a,b;
        	a=q1[i1].pos;
        	b=q2[i2].pos;
        	if(a+O1>b+O2)
        	{
        		i2++;
        	}
        	else if(a+O1<b+O2)
        	{
        		i1++;
        	}
        	else
        	{
                r[lr].pos=(uint32_t)(a+O1);
                r[lr].ref_id=q1[i1].ref_id;
                lr++;
                i1++;
                i2++;
        	}
        }
    }
    if(lr==0)
    {
    	p3->p_pos_list=NULL;
    	p3->p_pos_list_len=0;
    	return;
    }
    else if(lr!=I2)
    {
    	r=(struct pos*)realloc(r,sizeof(struct pos)*lr);
    }
	p3->p_pos_list=r;
	p3->p_pos_list_len=lr;
	p3->offset=0;
}
int32_t mergy_pos_list_multi(vector <struct position_list> L, struct position_list *r)
{
	if(L.size()==1)
	{
	    struct position_list c;
    	c.p_pos_list=(struct pos*)malloc(sizeof(struct pos)*L[0].p_pos_list_len);
    	c.p_pos_list_len=L[0].p_pos_list_len;
    	for(uint32_t i=0;i<L[0].p_pos_list_len;i++)
    	{
    		c.p_pos_list[i]=L[0].p_pos_list[i];
    	}
    	c.offset=L[0].offset;
		r->offset=c.offset;
		r->p_pos_list=c.p_pos_list;
		r->p_pos_list_len=c.p_pos_list_len;
		return -1;
	}
//    sort(L.begin(),L.end(),cmp);
    struct position_list r_tmp;
    mergy_pos_list(L[0],L[1],r);
    if(r->p_pos_list_len==0)
    {
    	r->p_pos_list_len=0;
        r->p_pos_list=NULL;
        return -1;
    }

    for(uint32_t i=2;i<L.size();i++)
    {
//    	cerr << "i:"<<i<<endl;
//    	cerr << "r->p_pos_list_len:" << r->p_pos_list_len << endl;
        if(r->p_pos_list_len==0)
        {
        	r->p_pos_list_len=0;
            r->p_pos_list=NULL;
            return -1;
        }
        else
        {
            r_tmp.offset=r->offset;
            r_tmp.p_pos_list=r->p_pos_list;
            r_tmp.p_pos_list_len=r->p_pos_list_len;
//            if(i==70)
//            {
//            	cout << r->p_pos_list_len << endl;
//            	cout << r->offset << endl;
//            	for(uint32_t kk=0;kk<r->p_pos_list_len;kk++)
//            	{
//            		cout <<(int)r->p_pos_list[kk].ref_id <<":" << r->p_pos_list[kk].pos << endl;
//            	}
//            	cout << L[i].p_pos_list_len << endl;
//            	cout << L[i].offset << endl;
//            	for(uint32_t kk=0;kk<L[i].p_pos_list_len;kk++)
//            	{
//            		cout <<(int)L[i].p_pos_list[kk].ref_id <<":" << L[i].p_pos_list[kk].pos << endl;
//            	}
//
//            }
            mergy_pos_list(L[i],r_tmp,r);
            free(r_tmp.p_pos_list);
        }
    }
    return 0;
}
void combn_pos_list(struct position_list a,struct position_list b,struct position_list* p1)
{
    uint32_t l;
    l=a.p_pos_list_len+b.p_pos_list_len;
    struct pos * p;
    p=(struct pos*)malloc(sizeof(struct pos)*l);

    uint32_t x=0,y=0,z=0;
    while(x<a.p_pos_list_len&&y<b.p_pos_list_len)
    {
    	if(a.p_pos_list[x].ref_id<b.p_pos_list[y].ref_id)
    	{
            p[z].pos=a.p_pos_list[x].pos+a.offset;
            p[z].ref_id=a.p_pos_list[x].ref_id;
            x++;
            z++;
    	}
    	else if(a.p_pos_list[x].ref_id>b.p_pos_list[y].ref_id)
    	{
            p[z].pos=b.p_pos_list[y].pos+b.offset;
            p[z].ref_id=b.p_pos_list[y].ref_id;
            y++;
            z++;
    	}
    	else
    	{
    		if(a.p_pos_list[x].pos+a.offset<b.p_pos_list[y].pos+b.offset)
    		{
                p[z].pos=a.p_pos_list[x].pos+a.offset;
                p[z].ref_id=a.p_pos_list[x].ref_id;
                x++;
                z++;
    		}
    		else if(a.p_pos_list[x].pos+a.offset>b.p_pos_list[y].pos+b.offset)
    		{
    			p[z].pos=b.p_pos_list[y].pos+b.offset;
				p[z].ref_id=b.p_pos_list[y].ref_id;
				y++;
				z++;
    		}
    		else
    		{
                p[z].pos=a.p_pos_list[x].pos+a.offset;
                p[z].ref_id=a.p_pos_list[x].ref_id;
                x++;
                y++;
                z++;
    		}
    	}
    }
    if(x<a.p_pos_list_len)
    {
        for(uint32_t i=x;i<a.p_pos_list_len;i++)
        {
            p[z].pos=a.p_pos_list[x].pos+a.offset;
            p[z].ref_id=a.p_pos_list[x].ref_id;
            x++;
            z++;
        }
    }
    else if(y<b.p_pos_list_len)
    {
        for(uint32_t i=y;i<b.p_pos_list_len;i++)
        {
            p[z].pos=b.p_pos_list[y].pos+b.offset;
            p[z].ref_id=b.p_pos_list[y].ref_id;
            y++;
            z++;
        }
    }

    p=(struct pos*)realloc(p,sizeof(struct pos)*z);

    p1->p_pos_list=p;
    p1->p_pos_list_len=z;
    p1->offset=0;
}
void combn_pos_list_multi(vector <struct position_list> a, struct position_list* p)
{
    struct position_list c;
    struct position_list c_tmp;

    if(a.size()==0)
    {
    	p->offset=0;
    	p->p_pos_list=NULL;
    	p->p_pos_list_len=0;
    	return;
    }

    if(a.size()<2)
    {
    	c.p_pos_list=(struct pos*)malloc(sizeof(struct pos)*a[0].p_pos_list_len);
    	c.p_pos_list_len=a[0].p_pos_list_len;
    	for(uint32_t i=0;i<a[0].p_pos_list_len;i++)
    	{
    		c.p_pos_list[i]=a[0].p_pos_list[i];
    	}
    	c.offset=a[0].offset;
    	p->offset=c.offset;
    	p->p_pos_list=c.p_pos_list;
    	p->p_pos_list_len=c.p_pos_list_len;
        return;
    }
    else
    {
        combn_pos_list(a[0],a[1],&c);
        for(uint32_t i=2;i<a.size();i++)
        {
            c_tmp=c;
            combn_pos_list(a[i],c_tmp,&c);
            free(c_tmp.p_pos_list);
        }
        p->p_pos_list=c.p_pos_list;
        p->p_pos_list_len=c.p_pos_list_len;
        p->offset=c.offset;
    }
}
