/*
 * inputRef.cpp
 *
 *  Created on: May 8, 2021
 *      Author: bio
 */

#include "inputRef.h"

void ReadSeq_ref(char **seq1,uint64_t *seq_length,char* p_ref)
{
	uint32_t buffer_size=256;
	char buffer_line[256];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
		return;
	}

	uint64_t total_size=0;
	fseek(fp,0,2);
	total_size=ftell(fp);

	char *seq;
	seq=(char*) malloc (sizeof(char)*total_size);

	fseek(fp,0,0);

	uint64_t len=0;
	while (fgets(buffer_line,buffer_size-1,fp)!=NULL)
	{
		if(buffer_line[0]=='>')
			continue;
		else
		{
			for(uint32_t i=0;i<buffer_size;i++)
			{
				if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
				{
					break;
				}
				if(buffer_line[i]>='a')
				{
					buffer_line[i]-=32;
				}
				if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
				{
					buffer_line[i]='A';
				}
				seq[len]=buffer_line[i];
				len++;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	*seq_length=len;
	*seq1=seq;
	cout << "the length of seq is: " << len << endl;
}
