/*
 * hashT.h
 *
 *  Created on: Feb 16, 2022
 *      Author: bio
 */

#ifndef HASHT_H_
#define HASHT_H_

#include "basic.h"
#define hash_table_size 0x100000


uint32_t hash_function(uint64_t *a,struct bit256KmerPara para);

struct ikmer
{
	uint64_t i;
	uint64_t *kmer;
};

int cmp(uint64_t *a,uint64_t *b,uint32_t l);

//unsigned hash_function(uint64_t *a);


void hash_table_initial(vector< vector<struct ikmer> > &a);


void hash_table_insert(vector< vector<struct ikmer> > aa,uint64_t *node,uint64_t i,bit256KmerPara para);


int hash_table_query(vector< vector<struct ikmer> > aa,uint64_t *node,uint64_t *res_index,vector <string> b,bit256KmerPara para);


void hash_table_delete(vector< vector<struct ikmer> > aa,uint64_t *node,uint64_t res_index,bit256KmerPara para);


void hash_table_destroy(vector< vector<struct ikmer> > a);


#endif /* HASHT_H_ */
