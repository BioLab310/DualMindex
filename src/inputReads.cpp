/*
 * inputfile.cpp
 *
 *  Created on: Mar 26, 2021
 *      Author: bio
 */

#include "inputReads.h"

void input_reads(char ** p_reads,uint32_t* num_read,char * p_readfilename,\
		uint32_t start_line,uint32_t task_size,uint32_t len_read)
{
	char *seq;
	seq=(char*) malloc (sizeof(char)*task_size*len_read + 1);
	memset(seq, 0, task_size*len_read+1);

	uint32_t buffer_size=256;
	char buffer_line[256] = {0};

	FILE *fp;
	fp = fopen(p_readfilename,"r+");
	if(fp==NULL)
	{
		cout <<"error from input reads: the Read file can not be open!" << endl;
	}

	for(uint32_t i=0;i<start_line;i++)
	{
		fgets(buffer_line,buffer_size,fp);
	}

	uint32_t size_tmp=0;
	uint32_t seq_index=0;
	while (fgets(buffer_line,buffer_size,fp)!=NULL)
	{
		//add code for dealing with reads containing other characters, such as 'N';

		memcpy(seq+seq_index,buffer_line,len_read);
		seq_index+=len_read;
		size_tmp++;
		if(size_tmp==task_size)
		{
			break;
		}
	}

	*p_reads=seq;
	*num_read=size_tmp;
}


