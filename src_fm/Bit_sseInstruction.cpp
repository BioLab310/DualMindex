/*
 * Bit_sseInstruction.cpp
 *
 *  Created on: Dec 31, 2019
 *      Author: pluto
 */

#include "Bit_sseInstruction.h"
#include <nmmintrin.h>
#include <stdio.h>
#include "print.h"

uint64_t CUTVECL[65];	//0的个数和下标相同   将64位的高位截去
uint64_t CUTVECH[65];	//1的个数和下标相同   将64位的低位截去

uint8_t ALLBITSET1_128[16] __aligned__ = {
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255,
/* 255 */255 };

void init_CUTVECL()
{
	CUTVECL[0] = 0l;
	CUTVECL[64] = ~0l;
	for(uint32_t i = 1; i < 64; i++)
	{
		CUTVECL[i] = CUTVECL[64] >> (64-i);
	}
}

void init_CUTVECH()
{
	CUTVECH[0] = 0l;
	CUTVECH[64] = ~0l;
	for(uint32_t i = 1; i < 64; i++)
	{
		CUTVECH[i] = CUTVECH[64] << (64-i);
	}
}

uint32_t popcnt_sse_u64(__m128i reg){
	uint32_t result = 0;
	result += _mm_popcnt_u64(*(uint64_t *)&reg);
	result += _mm_popcnt_u64(*((uint64_t *)&reg + 1));
	return result;
}

__m128i _mm_not_si128(__m128i reg)
{
	return _mm_andnot_si128(reg, *(__m128i *)ALLBITSET1_128);
}

__m128i get_cutm128i(uint32_t cutlen, uint64_t *CUTVEC)
{
	__m128i mcut;
	uint64_t cutmp[2] = {0};
	if(cutlen > 64)
	{
		cutmp[0] = ~0l;
		cutmp[1] = CUTVEC[cutlen - 64];
	}
	else
	{
		cutmp[0] = CUTVEC[cutlen];
		cutmp[1] = 0l;
	}
	mcut = _mm_set_epi64x(cutmp[1],cutmp[0]);
	return mcut;
}

uint32_t calcVec_Anum(__m128i regh, __m128i regl)
{
	__m128i reght, reglt, tmp;
	reght = _mm_not_si128(regh);
	reglt = _mm_not_si128(regl);
	tmp = _mm_and_si128(reght,reglt);
	return popcnt_sse_u64(tmp);
}

uint32_t calcVec_Cnum(__m128i regh, __m128i regl)
{
	__m128i reght, tmp;
	reght = _mm_not_si128(regh);
	tmp = _mm_and_si128(reght,regl);
	return popcnt_sse_u64(tmp);
}

uint32_t calcVec_Gnum(__m128i regh, __m128i regl)
{
	__m128i reglt, tmp;
	reglt = _mm_not_si128(regl);
	tmp = _mm_and_si128(regh,reglt);
	return popcnt_sse_u64(tmp);
}

uint32_t calcVec_Tnum(__m128i regh, __m128i regl)
{
	__m128i tmp;
	tmp = _mm_and_si128(regh,regl);
	return popcnt_sse_u64(tmp);
}

uint32_t calcVec_Acut(__m128i regh, __m128i regl, uint32_t cutlen, uint64_t *CUTVEC)
{
	__m128i reght, reglt, tmp, cutv;
	reght = _mm_not_si128(regh);
	reglt = _mm_not_si128(regl);
	tmp = _mm_and_si128(reght,reglt);
	cutv = get_cutm128i(cutlen, CUTVEC);
	tmp = _mm_and_si128(tmp,cutv);
	return popcnt_sse_u64(tmp);
}

uint32_t calcVec_Ccut(__m128i regh, __m128i regl, uint32_t cutlen, uint64_t *CUTVEC)
{
	__m128i reght, tmp, cutv;
	reght = _mm_not_si128(regh);
	tmp = _mm_and_si128(reght,regl);
	cutv = get_cutm128i(cutlen, CUTVEC);
	tmp = _mm_and_si128(tmp,cutv);
	return popcnt_sse_u64(tmp);
}

uint32_t calcVec_Gcut(__m128i regh, __m128i regl, uint32_t cutlen, uint64_t *CUTVEC)
{
	__m128i reglt, tmp, cutv;
	reglt = _mm_not_si128(regl);
	tmp = _mm_and_si128(regh,reglt);
	cutv = get_cutm128i(cutlen, CUTVEC);
	tmp = _mm_and_si128(tmp,cutv);
	return popcnt_sse_u64(tmp);
}

uint32_t calcVec_Tcut(__m128i regh, __m128i regl, uint32_t cutlen, uint64_t *CUTVEC)
{
	__m128i tmp, cutv;
	tmp = _mm_and_si128(regh,regl);
	cutv = get_cutm128i(cutlen, CUTVEC);
	tmp = _mm_and_si128(tmp,cutv);
	return popcnt_sse_u64(tmp);
}
