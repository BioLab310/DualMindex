/*
 * print.h
 *
 *  Created on: Jan 2, 2020
 *      Author: pluto
 */

#ifndef PRINT_H_
#define PRINT_H_


#include <stdint.h>
#include <tmmintrin.h>

void printbytevector(uint8_t *data, int length);
void print128_bit(__m128i var);


#endif /* PRINT_H_ */
