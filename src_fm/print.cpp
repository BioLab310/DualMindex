/*
 * print.cpp
 *
 *  Created on: Jan 2, 2020
 *      Author: pluto
 */


#include "print.h"
#include <stdio.h>

void printbytevector(uint8_t *data, int length) {
	int i;
	for (i = length-1; i >= 0; i--) {
		int m;
		for (m = 7; m >= 0; m--) {	//低位开始打印  低字节-->高字节  第一个字节到第十六个字节
			if (data[i] & (1ULL << m))
				printf("1");
			else
				printf("0");
		}
	}
}

void print128_bit(__m128i var) {
	uint8_t *val = (uint8_t*) &var;
	printbytevector(val, 16);

	printf("\n");

}

