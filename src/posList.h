/*
 * posList.h
 *
 *  Created on: Mar 23, 2022
 *      Author: bio
 */

#ifndef POSLIST_H_
#define POSLIST_H_

#include "basic.h"
#include "bitVec.h"
#include "DMindex_Oper.h"
#include "inputRef.h"
struct pos
{
	uint8_t ref_id;
	uint32_t pos;
};
struct Para_gen_list
{
	sFMindex FMidx;
	char *seq;
	uint64_t seq_length;
	uint32_t thread_num;
	uint32_t thread_id;
	char *p_bitVec;
	uint32_t v_bit_gap;
	uint32_t* p_Nb;
	uint32_t kmer_len;
	uint8_t ref_i;
//	struct pos **p_pos;
//	uint32_t *p_pos_len;

};
void gen_list(char * ref_name,char *p_bitVec_file,uint32_t v_bit_gap,char* p_Nb_file,uint32_t kmer_len,uint32_t thread_num);
//void gen_list(char * ref_name,char *p_bitVec,uint32_t v_bit_gap,char* p_Nb,uint32_t kmer_len);
void read_Pi(uint32_t** p_Pi,char* p_Pi_file);
void read_Ps(struct pos** p_Ps,char* p_Ps_file);
void read_Nb(uint32_t** p_Nb,char* p_Nb_file);
void read_vb(char** p_bitVec,char* p_bitVec_file);
uint32_t count_minimal_edges(char *p_bitVec_file);
#endif /* POSLIST_H_ */
