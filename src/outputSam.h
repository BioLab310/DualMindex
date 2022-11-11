/*
 * outputSam.h
 *
 *  Created on: Aug 9, 2021
 *      Author: bio
 */

#ifndef OUTPUTSAM_H_
#define OUTPUTSAM_H_

#include "basic.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"

#define kvec_t(type) struct { size_t n, m; type *a; }  //define A B 用A来表示B

//@abstract Structure for one alignment.
struct sam_t{
	size_t n, m;
	bam1_t* v;
};

struct  kvec_t_uint32_t{
  kvec_t(uint32_t) v;
} ;

struct OutSam{
	samFile *output_sam_file;    //typedef htsFile samFile   samFile作为htsFile的新名字，typedef A B 用B来表示A
	sam_hdr_t *sam_header;
	sam_t *sam_alignment;
};

//处理结构体OutSam的sam_hdr_t的成员变量
void output_sam_header(const char *output_file_path, char *p_ref_path, OutSam *output_sam) ;

void generate_bam1_t(uint8_t edit_distance, kstring_t *MD_tag, uint32_t mapping_start_position, \
		int32_t reference_sequence_index, uint8_t mapping_quality, uint16_t flag, \
		const char *query_name, uint16_t query_name_length, uint32_t *cigar, \
		uint32_t num_cigar_operations, const char *query, const char *query_qual, int32_t query_length, bam1_t *sam_alignment);

void outSam(OutSam *output_sam,const char *output_file_path, char *p_ref_path,uint8_t edit_distance, \
						int mapping_start_position,uint32_t reference_sequence_index, char *read_name, int read_name_length,\
						char* read_sequence,char*read_qual,int read_length);



#endif /* OUTPUTSAM_H_ */
