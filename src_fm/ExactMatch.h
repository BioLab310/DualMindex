/*
 * ExactMatch.h
 *
 *  Created on: 2019年7月13日
 *      Author: jupiter
 */

#ifndef EXACTMATCH_H_
#define EXACTMATCH_H_

#include "basic.h"
#include "read.h"

static const uint8_t char_to_uint8_table[256] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
												 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
static const char uint8_to_char_table[8] = {'A', 'C', 'G', 'T', 'N', 'A', 'N', 'N'};


uint32_t calc_OCC(const sFMindex &mem,char c,size_t pos);
uint32_t* calc_SArangeSeq(const sFMindex &mem,const char *read);
uint32_t* calc_SArangeChar(const sFMindex &mem,uint32_t *pre, char ch);
uint32_t calc_SA(const sFMindex &mem,size_t pos);

inline char C_next(char c)
{
	return uint8_to_char_table[char_to_uint8_table[c] + 1];
}

inline uint32_t calc_C(const sFMindex &mem,char c)
{
	return mem.c[char_to_uint8_table[c]];
}

inline uint32_t LF_Mapping(const sFMindex &mem,char c,size_t pos)
{
	//last first mapping
	return calc_C(mem,c) + calc_OCC(mem,c,pos) + 1;
}

inline uint32_t LF_Mapping_l(const sFMindex &mem,char c,size_t pos)
{
	//last first mapping
//	return char_to_uint8_table[c] != 4 ? calc_C(mem,c) + calc_OCC(mem,c,pos) + 1 : -1;
	return calc_C(mem,c) + calc_OCC(mem,c,pos) + 1;
}
inline uint32_t LF_Mapping_h(const sFMindex &mem,char c,size_t pos)
{
	//last first mapping
//	return calc_C(mem,c) + calc_OCC(mem,c,pos) + 1 - 1;
//	return char_to_uint8_table[c] != 4 ? calc_C(mem,c) + calc_OCC(mem,c,pos+1) : 0;
	return calc_C(mem,c) + calc_OCC(mem,c,pos+1);
}


#endif /* EXACTMATCH_H_ */
