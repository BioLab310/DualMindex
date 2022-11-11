/*
 * Bit_sseInstruction.h
 *
 *  Created on: Dec 31, 2019
 *      Author: pluto
 */

#ifndef BIT_SSEINSTRUCTION_H_
#define BIT_SSEINSTRUCTION_H_
#ifndef __aligned__
	#define __aligned__ __attribute__((aligned(16)))
#endif

#include <stdint.h>
#include <nmmintrin.h>

extern uint64_t CUTVECL[];
extern uint64_t CUTVECH[];

void init_CUTVECL();
void init_CUTVECH();
uint32_t popcnt_sse_u64(__m128i reg);
__m128i _mm_not_si128(__m128i reg);
__m128i get_cutm128i(uint32_t cutlen, uint64_t *CUTVEC);
uint32_t calcVec_Anum(__m128i regh, __m128i regl);
uint32_t calcVec_Cnum(__m128i regh, __m128i regl);
uint32_t calcVec_Gnum(__m128i regh, __m128i regl);
uint32_t calcVec_Tnum(__m128i regh, __m128i regl);
uint32_t calcVec_Acut(__m128i regh, __m128i regl, uint32_t cutlen, uint64_t *CUTVEC);
uint32_t calcVec_Ccut(__m128i regh, __m128i regl, uint32_t cutlen, uint64_t *CUTVEC);
uint32_t calcVec_Gcut(__m128i regh, __m128i regl, uint32_t cutlen, uint64_t *CUTVEC);
uint32_t calcVec_Tcut(__m128i regh, __m128i regl, uint32_t cutlen, uint64_t *CUTVEC);


#endif /* BIT_SSEINSTRUCTION_H_ */
