#include "basic.h"

void kmercpy(uint64_t *dest, const uint64_t * from, uint32_t unitper)
{
	if(dest == NULL || from == NULL)
	{
		printf("kmercpy error\n");
		exit(1);
	}
	for(uint32_t i = 0; i < unitper; i++)
	{
		dest[i] = from[i];
	}
}
void get_para(struct bit256KmerPara *para1,uint32_t kmer_length)
{
	para1->kmer1Len = kmer_length*2;
	para1->remainer1to64 = para1->kmer1Len%64;
	para1->kmer64Len = para1->kmer1Len/64+(para1->remainer1to64?1:0);
	if(para1->remainer1to64==0)
	{
		para1->codefor1 = 0;
	}
	else
	{
		para1->codefor1 = 0;
		for(uint32_t i=0;i<para1->remainer1to64-1;i++)
		{
			para1->codefor1 = para1->codefor1|1;
			para1->codefor1 = para1->codefor1<<1;
		}
		para1->codefor1 = para1->codefor1|1;
	}
}
void ReadSeq(char **seq1,uint32_t *seq_length,const char* p_ref)
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
	char buffer_line[256] = {0};

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
	}

	uint32_t len=0;
	while (fgets(buffer_line,buffer_size,fp)!=NULL)
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
					if(buffer_line[i]>='a')
					{
						buffer_line[i] -= 32;
					}
					if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
					{
						continue;
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
					if(buffer_line[i]>='a')
					{
						buffer_line[i] -= 32;
					}
					if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
					{
						continue;
					}
					seq[len]=buffer_line[i];
					len++;
				}
			}
		}
		memset(buffer_line,0,buffer_size*sizeof(char));
	}
	*seq_length=len;
	*seq1=seq;
//	cout << "the length of seq is: " << len << endl;
}
uint32_t cmp256BitKmer(uint64_t*a,uint64_t*b,uint32_t len)
{
	uint32_t r=2;
	for(uint32_t i=0;i<len;i++)
	{
		if(a[i]<b[i])
		{
			r=0;
			break;
		}
		else if(a[i]>b[i])
		{
			r=1;
			break;
		}
		else
		{
			if(i==len-1)
			{
				r=2;
				break;
			}
		}
	}
	return r;
}
void cal_hash_value_directly_256bit(char *seq,uint64_t * current,\
		struct bit256KmerPara para)
{
	uint64_t tmp;
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<para.kmer64Len;i++)
	{
		char* loop_tmp=k_mer_temp+32*i;
		if(i==para.kmer64Len-1&&para.remainer1to64!=0)
		{
			tmp=0;
			for(uint32_t j=0;j<para.remainer1to64/2-1;j++)
			{
				switch(loop_tmp[j])
				{
					case 'A':
						tmp=tmp<<2;
						break;
					case 'C':
						tmp=tmp|1;
						tmp=tmp<<2;
						break;
					case 'G':
						tmp=tmp|2;
						tmp=tmp<<2;
						break;
					case 'T':
						tmp=tmp|3;
						tmp=tmp<<2;
						break;
					default:
						tmp=tmp<<2;
						break;
				}
			}
			switch(loop_tmp[para.remainer1to64/2-1])
			{
				case 'A':
					break;
				case 'C':
					tmp=tmp|1;
					break;
				case 'G':
					tmp=tmp|2;
					break;
				case 'T':
					tmp=tmp|3;
					break;
				default:
					break;
			}
			current[i]=tmp;
		}
		else
		{
			tmp=0;
			for(uint32_t j=0;j<31;j++)
			{
				switch(loop_tmp[j])
				{
					case 'A':
						tmp=tmp<<2;
						break;
					case 'C':
						tmp=tmp|1;
						tmp=tmp<<2;
						break;
					case 'G':
						tmp=tmp|2;
						tmp=tmp<<2;
						break;
					case 'T':
						tmp=tmp|3;
						tmp=tmp<<2;
						break;
					default:
						tmp=tmp<<2;
						break;
				}
			}
			switch(loop_tmp[31])
			{
				case 'A':
					break;
				case 'C':
					tmp=tmp|1;
					break;
				case 'G':
					tmp=tmp|2;
					break;
				case 'T':
					tmp=tmp|3;
					break;
				default:
					break;
			}
			current[i]=tmp;
		}
	}
}
void cal_hash_value_indirectly_256bit(char *seq,uint64_t* original,uint64_t* current,\
		struct bit256KmerPara para)
{
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<para.kmer64Len-1;i++)
	{
		if(i==para.kmer64Len-2&&para.remainer1to64!=0)
		{
			current[i]=original[i]<<2;
			current[i]=current[i]|(original[i+1]>>(para.remainer1to64-2));
		}
		else
		{
			current[i]=original[i]<<2;
			current[i]=current[i]|(original[i+1]>>62);
		}
	}
	if(para.remainer1to64==0)
	{
		current[para.kmer64Len-1]=original[para.kmer64Len-1]<<2;
		switch(k_mer_temp[para.kmer1Len/2-1])
		{
			case 'A':
				//current=current|0;
				break;
			case 'C':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|1;
				//current=current<<2;
				break;
			case 'G':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|2;
				//current=current<<2;
				break;
			case 'T':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|3 ;
				//current=current<<2;
				break;
			default:
				break;
		}
	}
	else
	{
		current[para.kmer64Len-1]=original[para.kmer64Len-1]<<2;
		current[para.kmer64Len-1]=current[para.kmer64Len-1]&para.codefor1;
		switch(k_mer_temp[para.kmer1Len/2-1])
		{
			case 'A':
				//current=current|0;
				break;
			case 'C':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|1;
				//current=current<<2;
				break;
			case 'G':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|2;
				//current=current<<2;
				break;
			case 'T':
				current[para.kmer64Len-1]=current[para.kmer64Len-1]|3 ;
				//current=current<<2;
				break;
			default:
				break;
		}
	}
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

uint32_t counter(char a)
{
	uint32_t r;
	r=0;
	char x;
	x=1;
	if(a&x==x)
	{
		r++;
	}
	for(uint32_t i=1;i<=7;i++)
	{
		x=x<<1;
		if(a&x==x)
		{
			r++;
		}
	}

	return r;
}
