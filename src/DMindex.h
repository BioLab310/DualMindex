/*
 * DMindex.h
 *
 *  Created on: Feb 14, 2022
 *      Author: bio
 */

#ifndef DMINDEX_H_
#define DMINDEX_H_

#include "basic.h"
#include "inputRef.h"
#include "Hash.h"
#include "hashT.h"

void Total_Gen_MNRpath(char * pathFile,uint32_t kmer_len,char** nrpath,uint64_t * nrpath_len);
void Total_Gen_MNRpath_from_bigger_one(char * pathFile,uint32_t kmer_len,char** nrpath,uint64_t * nrpath_len);

#endif /* DMINDEX_H_ */
