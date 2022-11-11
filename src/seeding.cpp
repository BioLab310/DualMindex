/*
 * seeding.cpp
 *
 *  Created on: Mar 26, 2021
 *      Author: bio
 */

#include "seeding.h"

uint32_t * read_seed_build(uint32_t num_read,uint32_t e){
	uint32_t * read_sav_seed;
	read_sav_seed = (uint32_t *)malloc(sizeof(uint32_t)*num_read*(e+1));
	return read_sav_seed;
}
uint32_t get_seed_pos(uint32_t len_read,uint32_t e,uint32_t seed_id)
{
	uint32_t p=len_read/(e+1);
	uint32_t q=len_read%(e+1);
	if(seed_id<q)
	{
		return (p+1)*seed_id;
	}
	else
	{
		return (p+1)*q+p*(seed_id-q);
	}

}
void *Parallel_cal_seeds(void *p)
{
	struct bit256KmerPara para_tmp;

	struct para_multi_seed *tmp=(struct para_multi_seed *)p;
	uint64_t * seed_code;
	struct seed_message *seed_message;
	uint8_t *seed_thread_id;
	uint32_t seed_num=1;// the seed_num value is 1 bigger than the actual size of the seeds
	uint32_t seed_64num=tmp->read_len/(tmp->e+1)/32+1+1;
	seed_code=(uint64_t *)malloc(sizeof(uint64_t)*seed_num*seed_64num);
	seed_message=(struct seed_message *)malloc(sizeof(struct seed_message)*seed_num);
	seed_thread_id=(uint8_t *)malloc(sizeof(uint8_t)*seed_num);
	for(uint32_t i=tmp->start;i<=tmp->end;i++)
	{
		struct seed_message message_tmp;
		uint32_t p_SAI[2];
		char seed_tmp[64];
		//i*len_read is the start pos the current read
		for(uint32_t j=0;j<=tmp->e;j++)
		{
			//1)divide
			uint32_t start_pos;
			uint32_t seed_len;
			generate_seed(&start_pos,&seed_len,tmp->read_len,tmp->e,j);
			get_para(&para_tmp,seed_len);
			uint32_t start_tmp;
//			memset(seed_tmp,'\0',strlen(seed_tmp));
			start_tmp=(seed_num-1)*seed_64num;
			seed_code[start_tmp]=seed_len;
			cal_hash_value_directly_256bit(tmp->p_read+i*tmp->read_len+start_pos,\
					seed_code+start_tmp+1,para_tmp);
			seed_code=(uint64_t *)realloc(seed_code,sizeof(uint64_t)*(seed_num+1)*seed_64num);

			message_tmp.p_read_id_length=i;
			message_tmp.seed_len=seed_len;
			strncpy(seed_tmp,tmp->p_read+i*tmp->read_len+start_pos,seed_len);
			seed_tmp[seed_len]='\0';
			calc_SArangeSeq(tmp->FMidx,seed_tmp,p_SAI,p_SAI+1);
			message_tmp.saInterval[0]=p_SAI[0];
			message_tmp.saInterval[1]=p_SAI[1];

			seed_message[seed_num-1]=message_tmp;
			seed_message=(struct seed_message *)realloc(seed_message,sizeof(struct seed_message)*(seed_num+1));

			seed_thread_id[seed_num-1]=(bit256hashFFunction(seed_code+start_tmp+1,para_tmp))%tmp->thread_num;
			seed_thread_id=(uint8_t *)realloc(seed_thread_id,sizeof(uint8_t)*(seed_num+1));

			seed_num=seed_num+1;
		}
	}

	*(tmp->p_seed)=seed_code;
	*(tmp->p_seed_message)=seed_message;
	*(tmp->p_seed_thread_id)=seed_thread_id;
	*(tmp->p_seed_num)=seed_num-1;
}
void *Parallel_sav_seeds(void *p)
{
	struct bit256KmerPara para_tmp;

	struct para_multi_seed_final *tmp=(struct para_multi_seed_final*)p;
	uint32_t num_id=0;
	for(uint32_t i=0;i<tmp->thread_num;i++)
	{
		for(uint32_t j=0;j<tmp->p_seed_num[i];j++)
		{
			if(tmp->seed_thread_id[i][j]==tmp->thread_id)
			{
				get_para(&para_tmp,tmp->seed_message[i][j].seed_len);
				para_tmp.kmer64Len+=1;
				struct nodeBit c_tmp_hashtable;
				c_tmp_hashtable.hashValue=tmp->p_seed[i]+j*para_tmp.kmer64Len;
				int64_t x;
				x=getHashFTableValue(tmp->p_tmp,c_tmp_hashtable.hashValue,para_tmp);
				if(x!=-1)
				{
					(*(tmp->VS_message))[x].p_read_id=(uint32_t*)realloc((*(tmp->VS_message))[x].p_read_id,\
							sizeof(uint32_t)*((*(tmp->VS_message))[x].p_read_id_length+1));
					(*(tmp->VS_message))[x].p_read_id[(*(tmp->VS_message))[x].p_read_id_length]=tmp->seed_message[i][j].p_read_id_length;
					(*(tmp->VS_message))[x].p_read_id_length++;
					*(tmp->read_sav_seed+num_id+j)=x;
				}
				else
				{
					//if the current seed is not on the B+tree
					//save the message of the current seed
					struct seed_message tmp_message;
					tmp_message.e=0;
					tmp_message.seed_len=tmp->seed_message[i][j].seed_len;
					tmp_message.p_read_id=(uint32_t*)malloc(sizeof(uint32_t));
					tmp_message.p_read_id_length=0;
					tmp_message.p_read_id[tmp_message.p_read_id_length]=tmp->seed_message[i][j].p_read_id_length;
					tmp_message.p_read_id_length++;
					tmp_message.saInterval[0] = tmp->seed_message[i][j].saInterval[0];
					tmp_message.saInterval[1] = tmp->seed_message[i][j].saInterval[1];
//					calc_SA
					if(tmp_message.saInterval[0] <=tmp_message.saInterval[1])
					{
						tmp_message.seed_sa=(uint32_t *)malloc(sizeof(uint32_t)*(tmp_message.saInterval[1]-tmp_message.saInterval[0]+1));
						for(uint32_t m = 0; m <= tmp_message.saInterval[1]-tmp_message.saInterval[0]; m++)
						{
							tmp_message.seed_sa[m] = calc_SA(tmp->FMidx,tmp_message.saInterval[0]+m);
						}
					}

					pthread_mutex_lock(tmp->m);
					*(tmp->read_sav_seed+num_id+j)=*(tmp->array_id);
					c_tmp_hashtable.arrayID=*(tmp->array_id);
					*(tmp->array_id)=*(tmp->array_id)+1;
					(*(tmp->VS_message)).push_back(tmp_message);
					pthread_mutex_unlock(tmp->m);
					bit256insertHashFTable(tmp->p_tmp,c_tmp_hashtable,para_tmp);
				}
			}
		}
		num_id+=tmp->p_seed_num[i];
	}
}
void generate_seed(uint32_t*start_pos,uint32_t*seed_len,\
		uint32_t len_read,uint32_t e,uint32_t seed_id)
{
	uint32_t p=len_read/(e+1);
	uint32_t q=len_read%(e+1);
	if(seed_id<q)
	{
		*seed_len=p+1;
		*start_pos=(p+1)*seed_id;
	}
	else
	{
		*seed_len=p;
		*start_pos=(p+1)*q+p*(seed_id-q);
	}
}
struct NodeBit** seed_indexing(char* p_reads,uint32_t num_read,uint32_t len_read,\
		uint32_t e,uint32_t thread_num,vector <struct seed_message> &VS_message,uint32_t *reads_sav_seed)
{
	sFMindex FMidx;
	char *path=(char*)malloc(sizeof(char)*20);
	path=".";
	read_bfile2index(path,FMidx,0);
//	free(path);

	struct bit256KmerPara para_tmp;

	// 1)initialize the Bplus trees
	struct NodeBit** p_tmp;
	p_tmp=bit256initialHashFTable();

	//2)generate seeds and insert them into the BplusTrees
	if(thread_num==1)
	{
		uint64_t arrayId=0;
		for(uint32_t i=0;i<num_read;i++)
		{
			uint64_t * seed_code;
			seed_code=(uint64_t *)malloc(sizeof(uint64_t)*(len_read/(e+1)/32+1+1));//
			//i*len_read is the start pos the current read
			char seed_tmp[64];
			uint32_t p_SAI[2];
			for(uint32_t j=0;j<=e;j++)
			{
				//1)divide
				uint32_t start_pos;
				uint32_t seed_len;
				generate_seed(&start_pos,&seed_len,len_read,e,j);
				get_para(&para_tmp,seed_len);
				seed_code[0]=seed_len;
				cal_hash_value_directly_256bit(p_reads+i*len_read+start_pos,\
						seed_code+1,para_tmp);
				para_tmp.kmer64Len+=1;

				//2)insert
				struct nodeBit c_tmp_hashtable;
				c_tmp_hashtable.hashValue=seed_code;
				int64_t x;
				x=getHashFTableValue(p_tmp,c_tmp_hashtable.hashValue,para_tmp);
				if(x!=-1)
				{
					VS_message[x].p_read_id=(uint32_t*)realloc(VS_message[x].p_read_id,\
							sizeof(uint32_t)*(VS_message[x].p_read_id_length+1));
					VS_message[x].p_read_id[VS_message[x].p_read_id_length]=i;
					VS_message[x].p_read_id_length++;
					*(reads_sav_seed+i*(e+1)+j)=x;
				}
				else
				{
					//if the current seed is not on the B+tree
					c_tmp_hashtable.arrayID=arrayId;
					//save the message of the current seed
//					memset(seed_tmp,'\0',strlen(seed_tmp));
					strncpy(seed_tmp,p_reads+i*len_read+start_pos,seed_len);
					seed_tmp[seed_len]='\0';
					calc_SArangeSeq(FMidx,seed_tmp,p_SAI,p_SAI+1);
					struct seed_message tmp;
					tmp.e=0;
					tmp.seed_len=seed_len;
					tmp.p_read_id=(uint32_t*)malloc(sizeof(uint32_t));
					tmp.p_read_id_length=0;
					tmp.p_read_id[tmp.p_read_id_length]=i;
					tmp.p_read_id_length++;
					tmp.saInterval[0]=p_SAI[0];
					tmp.saInterval[1]=p_SAI[1];
					if(p_SAI[1]>=p_SAI[0])
					{
//						cout << p_SAI[1]-p_SAI[0]+1 << endl;
						tmp.seed_sa=(uint32_t *)malloc(sizeof(uint32_t)*(p_SAI[1]-p_SAI[0]+1));
						for(uint32_t k=0;k<p_SAI[1]-p_SAI[0]+1;k++)
						{
							tmp.seed_sa[k]=calc_SA(FMidx,p_SAI[0]+k);
						}
					}
					VS_message.push_back(tmp);
					bit256insertHashFTable(p_tmp,c_tmp_hashtable,para_tmp);
					*(reads_sav_seed+i*(e+1)+j)=arrayId;
					arrayId++;
				}
			}
		}
	}
	else
	{
		struct NodeBit** p_tmp;
		p_tmp=bit256initialHashFTable();
		pthread_t* t;
		t=(pthread_t*)malloc(sizeof(pthread_t)*thread_num);

		//1)calculate the seeds with multiple threads
		uint32_t p=num_read/thread_num;
		uint32_t q=num_read%thread_num;
		struct para_multi_seed * para;
		para=(struct para_multi_seed*)malloc(sizeof(struct para_multi_seed)*thread_num);
		uint32_t cur=0;
		uint64_t * p_seed[thread_num];
		struct seed_message *seed_message[thread_num];
		uint8_t *seed_thread_id[thread_num];
		uint32_t p_seed_num[thread_num];

		pthread_mutex_t m;
		pthread_mutex_init(&m,NULL);
		uint32_t array_id=0;

		for(uint32_t i=0;i<thread_num;i++)
		{
			para[i].start=cur;
			if(i<q)
			{
				cur=cur+p+1;
			}
			else
			{
				cur=cur+p;
			}
			para[i].FMidx=FMidx;
			para[i].end=cur-1;
			para[i].p_read=p_reads;
			para[i].read_len=len_read;
			para[i].p_seed=&(p_seed[i]);
			para[i].p_seed_message=&(seed_message[i]);
			para[i].p_seed_thread_id=&(seed_thread_id[i]);
			para[i].p_seed_num=p_seed_num+i;
			para[i].thread_id=i;
			para[i].thread_num=thread_num;
			para[i].e=e;
			if(pthread_create(t+i, NULL, Parallel_cal_seeds, (void*)(para+i))!=0)
			{
				cout << "error!" << endl;
			}
		}
		for(uint32_t i=0;i<thread_num;i++)
		{
			pthread_join(t[i], NULL);
		}

		free(para);

//		2)insert the seeds into the B+trees

		struct para_multi_seed_final * para_final;
		para_final=(struct para_multi_seed_final*)malloc(sizeof(struct para_multi_seed_final)*thread_num);
		for(uint32_t i=0;i<thread_num;i++)
		{
			para_final[i].VS_message=&VS_message;
			para_final[i].array_id=&array_id;
			para_final[i].m=&m;
			para_final[i].p_seed=p_seed;
			para_final[i].p_seed_num=p_seed_num;
			para_final[i].p_tmp=p_tmp;
			para_final[i].seed_message = seed_message;
			para_final[i].seed_thread_id = seed_thread_id;
			para_final[i].thread_id = i;
			para_final[i].thread_num = thread_num;
			para_final[i].FMidx = FMidx;
			para_final[i].read_sav_seed=reads_sav_seed;
			if(pthread_create(t+i, NULL, Parallel_sav_seeds, (void*)(para_final+i))!=0)
			{
				cout << "error!" << endl;
			}
		}
		for(uint32_t i=0;i<thread_num;i++)
		{
			pthread_join(t[i], NULL);
		}
		cout << VS_message.size() <<endl;
		free(para_final);
		free(t);
		for(uint32_t i=0;i<thread_num;i++)
		{
			if(p_seed[i]!=NULL)
			{
				free(p_seed[i]);
			}
			if(seed_message[i]!=NULL)
			{
				free(seed_message[i]);
			}
			if(seed_thread_id[i]!=NULL)
			{
				free(seed_thread_id[i]);
			}
		}
	}

	return p_tmp;
}
