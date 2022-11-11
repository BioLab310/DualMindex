/*
 * BitVec_build.h
 *
 *  Created on: Jan 4, 2020
 *      Author: pluto
 */

#ifndef BITVEC_INDEX_H_
#define BITVEC_INDEX_H_


#include "basic.h"
#include <nmmintrin.h>
#include "Bit_sseInstruction.h"
#include "build.h"
#define U64HSET1 0x8000000000000000

struct sFMindexv
{
	__m128i *bvh;
	__m128i *bvl;
	uint32_t bseqlen;
	uint32_t u64num;
	uint32_t dolpos;
	uint32_t *c;
	uint32_t *occa;
	uint32_t *sa;
	uint32_t occ_num;
	uint32_t sa_num;
};

void build_bVec(sFMindexv *pindxv);
void build_occAVec(sFMindexv *pindxv, struct build_para para);
uint32_t calc_OCCV(sFMindexv indxv,struct build_para para,char c,uint32_t pos);
uint32_t calc_OCCVu(sFMindexv indxv,struct build_para para,char c,uint32_t pos, uint32_t methodsel);

#endif /* BITVEC_INDEX_H_ */
