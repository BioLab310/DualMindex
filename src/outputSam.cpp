/*
 * outputSam.cpp
 *
 *  Created on: Aug 9, 2021
 *      Author: bio
 */

#include "outputSam.h"

//处理结构体OutSam的sam_hdr_t的成员变量
void output_sam_header(const char *output_file_path, char *p_ref_path,
		OutSam *output_sam) {
	//sam_hdr_t sam_header;
	output_sam->sam_header->n_targets = 1; //n_targets: number of reference sequences
	//lengths of the reference sequences
	output_sam->sam_header->target_len = (uint32_t*) malloc(
			output_sam->sam_header->n_targets * sizeof(uint32_t));
	//names of the reference sequences
	output_sam->sam_header->target_name = (char**) malloc(
			output_sam->sam_header->n_targets * sizeof(char*));
	kstring_t header_kstr = { 0, 0, NULL };
	// TODO(Haowen): add PG later
	//ksprintf(&header_kstr, "@PG\tID:FEM\tPN:FEM\tVN:%s\tCL:%s", FEM_PACKAGE, argv[0]);
	//for (int i = 1; i < argc; ++i) {
	//  ksprintf(&header_kstr, " %s", argv[i]);
	//}
	char *p_ref;
	uint32_t ref_length;
	ReadSeq(&p_ref, &ref_length, p_ref_path);
	output_sam->sam_header->target_len[0] = ref_length;
	output_sam->sam_header->target_name[0] = p_ref_path;
	ksprintf(&header_kstr, "@SQ\tSN:%s\tLN:%d\n",
			output_sam->sam_header->target_name[0],
			output_sam->sam_header->target_len[0]);
	//length of the plain text in the header (may be zero if the header has been edited)
	output_sam->sam_header->l_text = ks_len(&header_kstr);
//   plain text (may be NULL if the header has been edited)
	output_sam->sam_header->text = ks_str(&header_kstr);
//   header dictionary
	output_sam->sam_header->sdict = NULL;
//pointer to the extended header struct (internal use only)
	output_sam->sam_header->hrecs = NULL;
//reference count
	output_sam->sam_header->ref_count = 1;
	int htslib_err = sam_hdr_write(output_sam->output_sam_file,
			output_sam->sam_header);
	assert(htslib_err == 0);
}

void generate_bam1_t(uint8_t edit_distance, kstring_t *MD_tag,
		uint32_t mapping_start_position, int32_t reference_sequence_index,
		uint8_t mapping_quality, uint16_t flag, const char *query_name,
		uint16_t query_name_length, uint32_t *cigar,
		uint32_t num_cigar_operations, const char *query,
		const char *query_qual, int32_t query_length, bam1_t *sam_alignment) {
	/*! @typedef
	 *  @abstract Structure for core alignment information.
	 *  @field  pos     0-based leftmost coordinate
	 *  @field  tid     chromosome ID, defined by sam_hdr_t
	 *  @field  bin     bin calculated by bam_reg2bin()
	 *  @field  qual    mapping quality
	 *  @field  l_extranul length of extra NULs between qname & cigar (for alignment)
	 *  @field  flag    bitwise flag
	 *  @field  l_qname length of the query name
	 *  @field  n_cigar number of CIGAR operations
	 *  @field  l_qseq  length of the query sequence (read)
	 *  @field  mtid    chromosome ID of next read in template, defined by sam_hdr_t
	 *  @field  mpos    0-based leftmost coordinate of next read in template
	 *  @field  isize   observed template length ("insert size")
	 */
	sam_alignment->core.pos = mapping_start_position;
	sam_alignment->core.tid = reference_sequence_index;
	sam_alignment->core.qual = mapping_quality;
	sam_alignment->core.l_extranul = 0;
	if ((query_name_length + 1) % 4 != 0) {
		sam_alignment->core.l_extranul = 4 - ((query_name_length + 1) % 4);
	}
	sam_alignment->core.flag = flag;
	sam_alignment->core.l_qname = query_name_length + 1
			+ sam_alignment->core.l_extranul;
	sam_alignment->core.n_cigar = num_cigar_operations;
	sam_alignment->core.l_qseq = query_length;
	sam_alignment->core.mtid = -1;
	sam_alignment->core.mpos = -1;
	sam_alignment->core.isize = 0;

	/*! @typedef
	 @abstract Structure for one alignment.
	 @field  core       core information about the alignment
	 @field  id
	 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
	 @field  l_data     current length of bam1_t::data
	 @field  m_data     maximum length of bam1_t::data
	 @field  mempolicy  memory handling policy, see bam_set_mempolicy()
	 @discussion Notes:
	 1. The data blob should be accessed using bam_get_qname, bam_get_cigar,
	 bam_get_seq, bam_get_qual and bam_get_aux macros.  These returns pointers
	 to the start of each type of data.
	 2. qname is terminated by one to four NULs, so that the following
	 cigar data is 32-bit aligned; core.l_qname includes these trailing NULs,
	 while core.l_extranul counts the excess NULs (so 0 <= l_extranul <= 3).
	 3. Cigar data is encoded 4 bytes per CIGAR operation.
	 See the bam_cigar_* macros for manipulation.
	 4. seq is nibble-encoded according to bam_nt16_table.
	 See the bam_seqi macro for retrieving individual bases.
	 5. Per base qualilties are stored in the Phred scale with no +33 offset.
	 Ie as per the BAM specification and not the SAM ASCII printable method.
	 */
	/*
	 @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
	 8 for T and 15 for N. Two bases are packed in one byte with the base
	 at the higher 4 bits having smaller coordinate on the read. It is
	 recommended to use bam_seqi() macro to get the base.
	 */
	// First calculate the length of data
	sam_alignment->l_data = sam_alignment->core.l_qname
			+ (sam_alignment->core.n_cigar << 2)
			+ ((sam_alignment->core.l_qseq + 1) >> 1)
			+ sam_alignment->core.l_qseq;
	if (sam_alignment->l_data > sam_alignment->m_data) {
		sam_alignment->m_data = sam_alignment->l_data;
		free(sam_alignment->data);
		sam_alignment->data = (uint8_t*) calloc(sam_alignment->m_data,
				sizeof(uint8_t));
	}
	// copy qname
	memcpy(bam_get_qname(sam_alignment), query_name,
			query_name_length * sizeof(char));
	// add NULs after qname and before cigar, let me add one nul at the moment and see if okay.
	bam_get_qname(sam_alignment)[query_name_length] = '\0';
	// copy cigar
	memcpy(bam_get_cigar(sam_alignment), cigar,
			sam_alignment->core.n_cigar * sizeof(uint32_t));
	// set seq
	uint8_t *seq = bam_get_seq(sam_alignment);
	for (size_t i = 0; i < query_length; ++i) {
		bam_set_seqi(seq, i, seq_nt16_table[(uint8_t )query[i]]);
	}
	// copy seq qual
	uint8_t *seq_qual = bam_get_qual(sam_alignment);
	memcpy(seq_qual, query_qual, query_length * sizeof(char));
	// remove +33 offset
	for (int i = 0; i < query_length; ++i) {
		seq_qual[i] -= 33;
	}
	bam_aux_update_int(sam_alignment, "NM", edit_distance);
	bam_aux_update_str(sam_alignment, "MD", ks_len(MD_tag) + 1, ks_str(MD_tag));
}

//one aligment/
/*
 *reference_sequence_index:candidate_position >> 32
 *POSITIVE_DIRECTION 0
 *#define NEGATIVE_DIRECTION 1
 */
void outSam(OutSam *output_sam, const char *output_file_path, char *p_ref_path,
		uint8_t edit_distance, int mapping_start_position,
		uint32_t reference_sequence_index, char *read_name,
		int read_name_length, char* read_sequence, char*read_qual,
		int read_length) {
	//初始化结构体,处理结构体的samFile *output_sam_file这个成员变量;
	output_sam->sam_alignment = (sam_t*) malloc(sizeof(sam_t));
//		kv_init(output_sam->sam_alignment_kvec[i].v);
	output_sam->sam_alignment->m = 0;
	output_sam->sam_alignment->n = 0;
	output_sam->sam_alignment->v = 0;

	output_sam->output_sam_file = sam_open_format(output_file_path, "w", NULL);
	output_sam->sam_header = sam_hdr_init();
	output_sam_header(output_file_path, p_ref_path, output_sam);
	fprintf(stderr, "Initialize struct successfully!\n");

	//把映射结果写入结构体
	sam_t *sam_alignment;
	kstring_t MD_tag = { 0, 0, NULL };
	kvec_t_uint32_t cigar_uint32_t;
//		kv_init(cigar_uint32_t.v);
	cigar_uint32_t.v.a = 0;
	cigar_uint32_t.v.m = 0;
	cigar_uint32_t.v.n = 0;

	uint8_t mapping_quality = 255;
	uint16_t flag = 1 == 0 ? 0 : BAM_FREVERSE;   //mappings[mi].direction
	MD_tag.l = 0;
	//第一个参数参考基因候选位置，第二个参数read，第三个参数是阈值
	//r = banded_edit_distance(p.p_ref+can_pos,p.p_reads+i*p.len_read,p.len_read,p.e,&mes);
	generate_bam1_t(edit_distance, &MD_tag, mapping_start_position,
			reference_sequence_index, mapping_quality, flag, read_name,
			read_name_length, cigar_uint32_t.v.a, cigar_uint32_t.v.n,
			read_sequence, read_qual, read_length,
			output_sam->sam_alignment->v);

	//从结构体写入文件
	// Save the result into the file
	int htslib_err = sam_write1(output_sam->output_sam_file,
			output_sam->sam_header, output_sam->sam_alignment->v);
	assert(htslib_err >= 0);

}



