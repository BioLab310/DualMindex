/*
 * Locate.h
 *
 *  Created on: Mar 28, 2022
 *      Author: bio
 */

#ifndef LOCATE_H_
#define LOCATE_H_

#include "basic.h"
#include "inputFMindex.h"
#include "posList.h"
#include "bitVec.h"
#include "DMindex_Oper.h"

struct DMindex
{
	sFMindex FMidx;
	char* p_bitVec;
	uint32_t v_bit_gap;
	uint32_t* p_Nb;
	uint32_t* p_Pi;
	struct pos* p_Ps;
	uint32_t kmer_len;

};

struct shift_miniID
{
	uint32_t mini_edge_id;
	int32_t shift_value;
};

struct shift_edge
{
	char * e;
	int32_t shift_value;
};

struct position_list
{
	struct pos * p_pos_list;
    uint32_t p_pos_list_len;
    int32_t offset;
};
void mergy_pos_list(struct position_list p1,struct position_list p2,struct position_list *p3);
int32_t mergy_pos_list_multi(vector <struct position_list> L, struct position_list *r);
void combn_pos_list(struct position_list a,struct position_list b,struct position_list* p1);
void combn_pos_list_multi(vector <struct position_list> a, struct position_list* p);
int32_t Locate_edge(struct DMindex &pDM,char *e,uint32_t label,struct shift_miniID** p_r,uint32_t *p_r_len,int32_t shift_value);
int32_t Locate_path(struct DMindex &pDM,char *h,struct position_list *r);
int32_t Locate_edge_whole(struct DMindex &pDM,char *e,struct position_list *p);
#endif /* LOCATE_H_ */
