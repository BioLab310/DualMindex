#ifndef _basic_h_
#define _basic_h_
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include "stdio.h"
#include <fstream>
#include <ctime>
#include "stdlib.h"
#include <math.h>
#include <iostream>
#include <bitset>
#include "assert.h"
#include <pthread.h>
#include <sys/time.h>
#include <unistd.h> 	//access header file
#include <sys/types.h>	//mdkir header file
#include <sys/stat.h> 	//mdkir header file
using namespace std;
using std::bitset;
typedef struct hashValue
{
	int seqID;
	int pos;
}HASHV;
typedef struct range
{
	int start;
	int end;
	int len;

}RANGE;
typedef struct pair
{
	int pos;
	int fre;
}PAIR;
typedef struct M_node
{
	int i;
	int j;
	int k;
	int seg;
	int label;
	int xLen;
}MNode;
typedef struct hashline
{
	int id;
	int dv_pos;
}HASHLINE;
typedef struct group_stratend
{
	int start;
	int end;
}GR_STartEnd;
struct RefFilePath
{
	char **pRefFilePath;
	uint32_t NumberOfPathes;
};
#define _hash_gap_ 1500
#define _hash_gap_short_ 20
#define _MAX_GroupStrLen_ 1024*1024*50
#define _OCC_gap_ 10
#define _SA_gap_ 4

//#define _GENE_
#define _PRESS_OCC_
#define _PRESS_SA_
//#define _URL_
#define _AHP_FILTER_
//#define _WAVE_
#define _RECURED_
#define _testsegmentation_

#ifdef _URL_
	#define _XIGMA_SIZE_ 46
	#define _XIGMA_START_ 78
#else
	#define _XIGMA_SIZE_ 27
	#define _XIGMA_START_ 97
#endif

#define _BufferSize_ 257

void reverseq(char *seq);
void ReadSeq(char **seq1,uint32_t *seq_length,char* p_ref);
void getRefFilePathes(char* pathFile, struct RefFilePath* p);
#endif





