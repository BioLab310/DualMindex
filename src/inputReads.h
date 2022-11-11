/*
 * inputfile.h
 *
 *  Created on: Mar 26, 2021
 *      Author: bio
 */

#ifndef INPUTREADS_H_
#define INPUTREADS_H_
#include "basic.h"

void input_reads(char ** p_reads,uint32_t* num_read,char * p_readfilename,\
		uint32_t start_line,uint32_t task_size,uint32_t len_read);


#endif /* INPUTREADS_H_ */
