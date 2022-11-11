/*
 * Hash.h
 *
 *  Created on: Jan 24, 2019
 *      Author: bio
 */

#ifndef HASH_H_
#define HASH_H_

#include "basic.h"
#include "BplusTreeBit.h"
#define AddSize 1
#define MinSizeForSearch 100
#define HashSize 0x100000
#define HashFSize 0x1000000
uint32_t bit256hashFFunction(uint64_t *a,\
		struct bit256KmerPara para);
struct NodeBit** bit256initialHashFTable();
struct NodeBit * bit256initialSingleHashFTable();
void bit256insertHashFTable(struct NodeBit** pk,\
		struct nodeBit c_tmp_hashtable,\
		struct bit256KmerPara para);
void bit256freeHashFTable(struct NodeBit **ph,\
		struct bit256KmerPara para);
int64_t getHashFTableValue(struct NodeBit** pk,\
		uint64_t *a,\
		struct bit256KmerPara para);

/***************************hashForBitCharKmer*********************************/
struct bit256KmerPosition{
	uint64_t 	kmer[4];
	uint64_t 	position;
};
struct bit256Hash{
	uint64_t** p_kmer;
	uint64_t** p_pos;
	uint64_t* len;
	uint64_t* maxlen;
};
uint32_t bit256HashFunction(uint64_t* a,struct bit256KmerPara para);
uint32_t bit256HashFunction_left(uint64_t* a,struct bit256KmerPara para);
struct bit256Hash* initial256BitHashTable();
void insert256BitHashTable(struct bit256Hash *ph,uint64_t* kmer,uint64_t pos,bit256KmerPara para);
uint64_t binary256BitSearch(uint64_t * ppo,uint64_t* p,uint64_t len,uint64_t* y,bit256KmerPara para);
uint64_t search256BitHashTable(struct bit256Hash *ph,uint64_t* y,bit256KmerPara para);
void freeHash256BitTable(struct bit256Hash *ph);
/***************************hashForBitCharKmer*********************************/

#endif /* HASH_H_ */
