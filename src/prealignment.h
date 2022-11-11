/*
 * prealignment.h
 *
 *  Created on: May 25, 2020
 *      Author: pluto
 */

#ifndef PREALIGNMENT_H_
#define PREALIGNMENT_H_

#include "basic.h"
#include <xmmintrin.h>
#include <tmmintrin.h>
#include <emmintrin.h>
#include <nmmintrin.h>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration.hpp>
#include <boost/preprocessor/arithmetic.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include "exactMatchFMindex.h"
#ifndef __aligned__
	#define __aligned__ __attribute__((aligned(16)))
#endif

#define SSE_BIT_LENGTH		128
#define SSE_BYTE_NUM		BOOST_PP_DIV(SSE_BIT_LENGTH, 8)
#define BASE_SIZE1			1
#define SSE_BASE_NUM1		BOOST_PP_DIV(SSE_BIT_LENGTH, BASE_SIZE1)

//static const uint8_t char_to_uint8_table[256] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

static inline uint8_t char_to_uint8(char c) {
	if(c >= 256)
	{printf("error %s\n", __FUNCTION__);}
  return char_to_uint8_table[(uint8_t)c];
}

int banded_edit_distance(const char *pattern, const char *text, \
		int text_length, const int error_threshold,int * mapping_end_posi);
void sse3_convert2bit1(char *str, uint8_t *bits0, uint8_t *bits1);
int pre_alignment(char *read, uint32_t read_len, char *ref, uint32_t start, uint32_t end, uint32_t max_error);
extern uint8_t read_vec0_t[SSE_BYTE_NUM] __aligned__;
extern uint8_t read_vec1_t[SSE_BYTE_NUM] __aligned__;


#endif /* PREALIGNMENT_H_ */
