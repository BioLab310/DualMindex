/*
 * bitVec.h
 *
 *  Created on: Feb 27, 2022
 *      Author: bio
 */

#ifndef BITVEC_H_
#define BITVEC_H_

#include "basic.h"
#include "inputFMindex.h"
#include "DMindex_Oper.h"

struct Para_gen_vec
{
	char* p_mnrPath;
	char* v_bit;
	uint32_t kmer_len;
	uint32_t start;
	uint32_t end;
	sFMindex FMidx;
};
uint32_t check_Minimal_Edge(char* e,sFMindex &FMidx);
void genBitVec(char *MNRpath, char **p_bitVec, uint64_t *p_bitVec_len,uint32_t kmer_len,uint32_t thread_num);
void genBitVec_total(char* p_file,uint32_t kmer_len,uint32_t thread_num,uint32_t v_bit_gap);
uint32_t cal_ones(char a);
uint32_t is_minimal_edge(char *p_bitVec,uint32_t id);
int64_t get_minimal_edge_id(char *p_bitVec,uint32_t id,uint32_t v_bit_gap,uint32_t* p_Nb);
#endif /* BITVEC_H_ */
