#include "basic.h"

void reverseq(char *seq)
{
	char *front = seq;
	char *rear = seq + strlen(seq) - 1;
	char tmp;
	while(front < rear)
	{
		tmp = *front;
		*front = *rear;
		*rear = tmp;
		front++;
		rear--;
	}
}

void fill_char_with_four_char(uint8_t* current,char*p)
{
	for(uint32_t i=0;i<3;i++)
	{
		switch(p[i])
		{
			case 'A':
				*current=*current<<2;
				break;
			case 'C':
				*current=*current|1;
				*current=*current<<2;
				break;
			case 'G':
				*current=*current|2;
				*current=*current<<2;
				break;
			case 'T':
				*current=*current|3;
				*current=*current<<2;
				break;
			default:
				*current=*current<<2;
				break;
		}
	}
	switch(p[3])
	{
		case 'A':
			break;
		case 'C':
			*current=*current|1;
			break;
		case 'G':
			*current=*current|2;
			break;
		case 'T':
			*current=*current|3;
			break;
		default:
			break;
	}
}
void ReadSeq_bit(uint8_t **seq1,uint64_t *seq_length,char* p_ref)
{
	uint32_t cursize;
	uint32_t maxsize;
	uint32_t addsize;

	maxsize=pow(2,20);
	addsize=pow(2,20);

	uint8_t *seq;//=seq1;
	seq=(uint8_t*) malloc (sizeof(uint8_t)*maxsize);
	cursize=maxsize;

	uint32_t buffer_size=_BufferSize_;
	char buffer_line[_BufferSize_];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint64_t len=0;
	uint32_t cycle=0;
	uint32_t number_left=0;
	uint32_t number_buffer_char;
	char array_left[4];
	while (1)
	{
		if(fgets(buffer_line,buffer_size-number_left,fp)==NULL)
		{
			break;
		}
		if(buffer_line[0]=='>')
		{
			continue;
		}
		else
		{
			char buffer_line_tmp[_BufferSize_];
			if(number_left!=0)
			{
				strcpy(buffer_line_tmp,array_left);
				strcpy(buffer_line_tmp+number_left,buffer_line);
			}
			else
			{
				strcpy(buffer_line_tmp,buffer_line);
			}

			uint32_t buffer_line_char_number=strlen(buffer_line_tmp);
			if(buffer_line_tmp[buffer_line_char_number-1]=='\n')
			{
				buffer_line_char_number--;
			}
			number_left=buffer_line_char_number%4;
			number_buffer_char=buffer_line_char_number/4;
			if(number_left!=0)
			{
				for(uint32_t i=0;i<number_left;i++)
				{
					array_left[i]=buffer_line_tmp[buffer_line_char_number-(number_left-i)];
				}
				array_left[number_left]='\0';
			}

			if(len+number_buffer_char<cursize)
			{
				for(uint32_t i=0;i<number_buffer_char;i++)
				{
					uint8_t tmp=0;

					for(uint32_t j=0;j<4;j++)
					{
						if(buffer_line_tmp[4*i+j]>='a')
						{
							buffer_line_tmp[4*i+j]-=32;
						}
					}
					fill_char_with_four_char(&tmp,buffer_line_tmp+4*i);
					seq[len]=tmp;
					len++;
				}
			}
			else
			{
				seq=(uint8_t*) realloc (seq,sizeof(uint8_t)*(cursize+addsize));
				cursize=cursize+addsize;
				for(uint32_t i=0;i<buffer_line_char_number/4;i++)
				{
					uint8_t tmp=0;

					for(uint32_t j=0;j<4;j++)
					{
						if(buffer_line_tmp[4*i+j]>='a')
						{
							buffer_line_tmp[4*i+j]-=32;
						}
					}
					fill_char_with_four_char(&tmp,buffer_line_tmp+4*i);
					seq[len]=tmp;
					len++;
				}
				cout <<"add 1024*1024 byte for seq: " << cycle++ <<endl;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	*seq_length=len*4+number_left;
	if(number_left!=0)
	{
		for(uint32_t i=number_left;i<4;i++)
		array_left[i]='A';
		uint8_t tmp;
		fill_char_with_four_char(&tmp,array_left);
		seq[len]=tmp;
		len++;
	}

	seq1[0]=seq;

	uint8_t* seq_shift1=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift1[i]=seq[i]<<2;
		seq_shift1[i]=seq_shift1[i]|(seq[i+1]>>6);
	}
	seq_shift1[len-1]=seq[len-1]<<2;
	seq1[1]=seq_shift1;

	uint8_t* seq_shift2=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift2[i]=seq_shift1[i]<<2;
		seq_shift2[i]=seq_shift2[i]|(seq_shift1[i+1]>>6);
	}
	seq_shift2[len-1]=seq_shift1[len-1]<<2;
	seq1[2]=seq_shift2;

	uint8_t* seq_shift3=(uint8_t*) malloc (sizeof(uint8_t)*len);
	for(uint32_t i=0;i<len-1;i++)
	{
		seq_shift3[i]=seq_shift2[i]<<2;
		seq_shift3[i]=seq_shift3[i]|(seq_shift2[i+1]>>6);
	}
	seq_shift3[len-1]=seq_shift2[len-1]<<2;
	seq1[3]=seq_shift3;

	cout << "the length of seq is: " << *seq_length << endl;
}

void ReadSeq(char **seq1,uint32_t *seq_length,char* p_ref)
{
	uint32_t cursize;
	uint32_t maxsize;
	uint32_t addsize;

	maxsize=pow(2,20);
	addsize=pow(2,20);

	char *seq;
	seq=(char*) malloc (sizeof(char)*maxsize);
	cursize=maxsize;

	uint32_t buffer_size=256;
	char buffer_line[256];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint32_t len=0;
	uint32_t cycle=0;
	while (fgets(buffer_line,buffer_size-1,fp)!=NULL)
	{
		if(buffer_line[0]=='>')
			continue;
		else
		{
			if(len+buffer_size<cursize)
			{
				for(uint32_t i=0;i<buffer_size;i++)
				{
					if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
					{
						break;
					}
					if(buffer_line[i]=='?')
					{
						buffer_line[i]='b';
					}
					else if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
					{
						buffer_line[i]='A';
					}
					seq[len]=buffer_line[i];
					len++;
				}
			}
			else
			{
				seq=(char*) realloc (seq,sizeof(char)*(cursize+addsize));
				cursize=cursize+addsize;
				for(uint32_t i=0;i<buffer_size;i++)
				{
					if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
					{
						break;
					}
					if(buffer_line[i]=='?')
					{
						buffer_line[i]='b';
					}
					else if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
					{
						buffer_line[i]='A';
					}
					seq[len]=buffer_line[i];
					len++;
				}
//				cout <<"add 1 MB for seq: " << cycle++ <<endl;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	*seq_length=len;
	*seq1=seq;
	cout << "the length of seq is: " << len << endl;
}
void getRefFilePathes(char* pathFile, struct RefFilePath* p)
{
	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	ifstream int_input;

	cout << pathFile << endl;
	int_input.open(pathFile,ios::in);
	uint64_t total_hash_number=0;
	char hv_tmp[64];
	while(int_input.getline(hv_tmp,63))
	{
		total_hash_number++;
	}
	int_input.close();
	cout << "total ref number is: " <<total_hash_number << endl;

	p->pRefFilePath=(char**)malloc(sizeof(char*)*total_hash_number);
	for(uint32_t i=0;i<total_hash_number;i++)
	{
		p->pRefFilePath[i]=(char*)malloc(sizeof(char)*256);
	}

	uint32_t line_ref_path=0;
	int_input.open(pathFile,ios::in);
	while(int_input.getline(hv_tmp,63))
	{
		strcpy(p->pRefFilePath[line_ref_path],hv_tmp);
		line_ref_path++;
	}
	int_input.close();
	p->NumberOfPathes=line_ref_path;

	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
}
