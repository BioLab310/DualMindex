/*
 * DMindex_Oper.h
 *
 *  Created on: Feb 28, 2022
 *      Author: bio
 */

#ifndef DMINDEX_OPER_H_
#define DMINDEX_OPER_H_

#include "basic.h"
#include "inputFMindex.h"
#include "exactMatchFMindex.h"

uint32_t cal_edge_id(sFMindex &FMidx,char* e);

uint32_t cal_edge_id_new(sFMindex &FMidx,char* e,uint32_t l);

void cal_SA_oneLeftShift(sFMindex &FMidx,uint32_t l,uint32_t h,char c,uint32_t *x,uint32_t *y);


uint32_t cal_Node_indegree(sFMindex &FMidx,char* n,char* ad);


uint32_t cal_Node_outdegree(sFMindex &FMidx,char* n,char* ad);


uint32_t is_Node_or_Edge_the_member_of_graph(sFMindex &FMidx,char* n,uint32_t* l,uint32_t* h);


#endif /* DMINDEX_OPER_H_ */
