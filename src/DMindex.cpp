/*
 * DMindex.cpp
 *
 *  Created on: Feb 14, 2022
 *      Author: bio
 */

#include "DMindex.h"

void cal_kmer(string a, bit256KmerPara para, int label,uint64_t * kmer)
{
	char* p=(char*)malloc(sizeof(char)*a.size());
	for(uint32_t i=0;i<a.size();i++)
	{
		p[i]=a[i];
	}
	if(label==0)
	{
		cal_hash_value_directly_256bit(p,kmer,para);
	}
	else
	{
		cal_hash_value_directly_256bit(p+(a.size()-para.kmer1Len/2),kmer,para);
	}
	free(p);
}
void Gen_NRpaths(char * pathFile,uint32_t kmer_len,char** nrpath,uint64_t * nrpath_len)
{
	bit256KmerPara para;
	get_para(&para,kmer_len+1);

//	bit256KmerPara para_x;
//	get_para(&para_x,kmer_len);
//	string x="GCGAAGGCATGAGGGCACAAA";
//
//	uint64_t kmer_current_x[4];
//
//	for(uint32_t i=0;i<4;i++)
//	{
//		kmer_current_x[i]=0;
//	}
//	cal_kmer(x,para_x,0,kmer_current_x);
//	cout << kmer_current_x[0] << endl;

//1 read sequence names
	struct RefFilePath p_ref_path;
	getRefFilePathes(pathFile, &p_ref_path);
//2 initialize the hash table for saving the (k+1)-mers
	struct bit256Hash* p_root=NULL;
	p_root=initial256BitHashTable();
//3 initialize the temp char* for saving the Current NRpaths
	char * nrp=NULL;
	uint64_t nrp_lenl;
	uint64_t M_nrp_lenl;
	uint32_t add_len=1024;
	nrp=(char*)malloc(sizeof(char)*1024*1024);
	nrp_lenl=0;
	M_nrp_lenl=1024*1024;
//4 loop for the sequences and generate the NRpaths
	char *seq;
	uint64_t seq_length;
	uint64_t seq_k_current[4];
	for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
	{
		cout << ref_i+1 << ":"<< p_ref_path.NumberOfPathes <<endl;
		ReadSeq_ref(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		char *temp_seq;
		temp_seq=seq;
		//loop for the (k+1)-mers
		cal_hash_value_directly_256bit(temp_seq,seq_k_current,para);
		int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
		if(arrayID==0||nrp_lenl==0)
		{
			insert256BitHashTable(p_root,seq_k_current,1,para);
			if(nrp[nrp_lenl-1]=='$'||nrp_lenl==0)
			{
				if(nrp_lenl+kmer_len+1>=M_nrp_lenl)
				{
					nrp=(char*)realloc(nrp,sizeof(char)*(M_nrp_lenl+add_len));
					M_nrp_lenl=M_nrp_lenl+add_len;

				}
				for(uint32_t i=0;i<kmer_len+1;i++)
				{
					nrp[nrp_lenl]=temp_seq[i];
					nrp_lenl++;
				}
			}
			else
			{
				// in this case nrp[nrp_lenl-1]!='$'&& nrp_lenl!=0 && the first (k+1)-mer of a new reference
				if(nrp_lenl+1+kmer_len+1>=M_nrp_lenl)
				{
					nrp=(char*)realloc(nrp,sizeof(char)*(M_nrp_lenl+add_len));
					M_nrp_lenl=M_nrp_lenl+add_len;

				}
				nrp[nrp_lenl]='$';
				nrp_lenl++;
				for(uint32_t i=0;i<kmer_len+1;i++)
				{
					nrp[nrp_lenl]=temp_seq[i];
					nrp_lenl++;
				}
				// the following two line are the original operation
//				nrp[nrp_lenl]=temp_seq[kmer_len];
//				nrp_lenl++;
			}
		}
		else
		{
			if(nrp[nrp_lenl-1]!='$')
			{
				nrp[nrp_lenl]='$';
				nrp_lenl++;
			}
		}

		for(uint32_t ii=1;ii<seq_length-kmer_len;ii++)
		{
			cal_hash_value_indirectly_256bit(temp_seq+ii,seq_k_current,seq_k_current,para);
			int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
			if(arrayID==0||nrp_lenl==0)
			{
				insert256BitHashTable(p_root,seq_k_current,1,para);
				if(nrp[nrp_lenl-1]=='$'||nrp_lenl==0)
				{
					if(nrp_lenl+kmer_len+1>=M_nrp_lenl)
					{
						nrp=(char*)realloc(nrp,sizeof(char)*(M_nrp_lenl+add_len));
						M_nrp_lenl=M_nrp_lenl+add_len;

					}
					for(uint32_t i=0;i<kmer_len+1;i++)
					{
						nrp[nrp_lenl]=temp_seq[ii+i];
						nrp_lenl++;
					}
				}
				else
				{
					if(nrp_lenl+1>=M_nrp_lenl)
					{
						nrp=(char*)realloc(nrp,sizeof(char)*(M_nrp_lenl+add_len));
						M_nrp_lenl=M_nrp_lenl+add_len;

					}
					nrp[nrp_lenl]=temp_seq[ii+kmer_len];
					nrp_lenl++;
				}
			}
			else
			{
				if(nrp[nrp_lenl-1]!='$')
				{
					nrp[nrp_lenl]='$';
					nrp_lenl++;
				}
			}
		}
	}

	*nrpath=nrp;
	*nrpath_len=nrp_lenl;
	cout <<nrp_lenl << endl;
//	cout << nrp << endl;
	freeHash256BitTable(p_root);
}

void Gen_NRpaths_from_biger_one(char * p_file,uint32_t kmer_len,char** nrpath,uint64_t * nrpath_len)
{
	//from aa22 to generate aa18 (example)

	char *MNRpath;
	uint64_t MNRpath_length;

	FILE * file = NULL;
	file = fopen(p_file,"r");
	if(file)
	{
		fseek(file,0,SEEK_END);//将文件内部的指针指向文件末尾
		uint64_t len = ftell(file);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
		MNRpath_length=len;
		rewind(file);//将文件内部的指针重新指向一个流的开头
		MNRpath = (char *)malloc(len+1*sizeof(char));//申请内存空间，len*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
		//用malloc申请的内存是没有初始值的，如果不赋值会导致写入的时候找不到结束标志符而出现内存比实际申请值大，写入数据后面跟随乱码的情况
		memset(MNRpath,0,len+1);//将内存空间都赋值为‘\0’

		fread(MNRpath, len, 1, file);
		for(uint32_t i=0;i<len;i++)
		{
			if(MNRpath[i]!='A'&&MNRpath[i]!='C'&&MNRpath[i]!='G'&&MNRpath[i]!='T'&&MNRpath[i]!='?')
			{
				MNRpath[i]='A';
			}
		}
		fclose(file);
	}

	bit256KmerPara para;
	get_para(&para,kmer_len+1);

//1 read sequence names
//	struct RefFilePath p_ref_path;
//	getRefFilePathes(pathFile, &p_ref_path);
//2 initialize the hash table for saving the (k+1)-mers
	struct bit256Hash* p_root=NULL;
	p_root=initial256BitHashTable();
//3 initialize the temp char* for saving the Current NRpaths
	char * nrp=NULL;
	uint64_t nrp_lenl;
	uint64_t M_nrp_lenl;
	uint32_t add_len=1024;
	nrp=(char*)malloc(sizeof(char)*1024*1024);
	nrp_lenl=0;
	M_nrp_lenl=1024*1024;
//4 loop for the sequences and generate the NRpaths
	char *seq;
	uint32_t seq_length;
	uint64_t seq_k_current[4];
	for(uint64_t ref_i=0;ref_i<MNRpath_length-1;)
	{
//		cout << ref_i+1 << ":"<< p_ref_path.NumberOfPathes <<endl;
//		ReadSeq_ref(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
		uint32_t x=0;
		while(MNRpath[ref_i+x]!='?')
		{
			x++;
		}
		seq=MNRpath+ref_i;
		seq_length=x;

		char *temp_seq;
		temp_seq=seq;
		//loop for the (k+1)-mers
		cal_hash_value_directly_256bit(temp_seq,seq_k_current,para);
		int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
		if(arrayID==0||nrp_lenl==0)
		{
			insert256BitHashTable(p_root,seq_k_current,1,para);
			if(nrp[nrp_lenl-1]=='$'||nrp_lenl==0)
			{
				if(nrp_lenl+kmer_len+1>=M_nrp_lenl)
				{
					nrp=(char*)realloc(nrp,sizeof(char)*(M_nrp_lenl+add_len));
					M_nrp_lenl=M_nrp_lenl+add_len;

				}
				for(uint32_t i=0;i<kmer_len+1;i++)
				{
					nrp[nrp_lenl]=temp_seq[i];
					nrp_lenl++;
				}
			}
			else
			{
				// in this case nrp[nrp_lenl-1]!='$'&& nrp_lenl!=0 && the first (k+1)-mer of a new reference
				if(nrp_lenl+1+kmer_len+1>=M_nrp_lenl)
				{
					nrp=(char*)realloc(nrp,sizeof(char)*(M_nrp_lenl+add_len));
					M_nrp_lenl=M_nrp_lenl+add_len;

				}
				nrp[nrp_lenl]='$';
				nrp_lenl++;
				for(uint32_t i=0;i<kmer_len+1;i++)
				{
					nrp[nrp_lenl]=temp_seq[i];
					nrp_lenl++;
				}
				// the following two line are the original operation
//				nrp[nrp_lenl]=temp_seq[kmer_len];
//				nrp_lenl++;
			}
		}
		else
		{
			if(nrp[nrp_lenl-1]!='$')
			{
				nrp[nrp_lenl]='$';
				nrp_lenl++;
			}
		}

		for(uint32_t ii=1;ii<seq_length-kmer_len;ii++)
		{
			cal_hash_value_indirectly_256bit(temp_seq+ii,seq_k_current,seq_k_current,para);
			int64_t arrayID = search256BitHashTable(p_root,seq_k_current,para);
			if(arrayID==0||nrp_lenl==0)
			{
				insert256BitHashTable(p_root,seq_k_current,1,para);
				if(nrp[nrp_lenl-1]=='$'||nrp_lenl==0)
				{
					if(nrp_lenl+kmer_len+1>=M_nrp_lenl)
					{
						nrp=(char*)realloc(nrp,sizeof(char)*(M_nrp_lenl+add_len));
						M_nrp_lenl=M_nrp_lenl+add_len;

					}
					for(uint32_t i=0;i<kmer_len+1;i++)
					{
						nrp[nrp_lenl]=temp_seq[ii+i];
						nrp_lenl++;
					}
				}
				else
				{
					if(nrp_lenl+1>=M_nrp_lenl)
					{
						nrp=(char*)realloc(nrp,sizeof(char)*(M_nrp_lenl+add_len));
						M_nrp_lenl=M_nrp_lenl+add_len;

					}
					nrp[nrp_lenl]=temp_seq[ii+kmer_len];
					nrp_lenl++;
				}
			}
			else
			{
				if(nrp[nrp_lenl-1]!='$')
				{
					nrp[nrp_lenl]='$';
					nrp_lenl++;
				}
			}
		}
		ref_i=ref_i+x+1;
	}

	*nrpath=nrp;
	*nrpath_len=nrp_lenl;
	cout <<nrp_lenl << endl;
//	cout << nrp << endl;
	freeHash256BitTable(p_root);
}


uint32_t find_lr(string query,vector <string> v_nrp,\
		vector< vector<struct ikmer> > hash_table_l,\
		vector< vector<struct ikmer> > hash_table_r,\
		uint32_t *i, uint32_t *j,bit256KmerPara para)
{
	// 0: not find
	// 1: find left match
	// 2: find right match
	// 3: find double match
	uint64_t kmer_current[4];
	uint64_t kmer_current1[4];
	uint64_t r,r1,r2;
//	cal_kmer(query, para, 0,kmer_current);
//	cal_kmer(query, para, 1,kmer_current1);
//	if(hash_table_query(hash_table_r,kmer_current,&r1,v_nrp,para)!=0)
//	{
//		if(hash_table_query(hash_table_l,kmer_current1,&r2,v_nrp,para)!=0)
//		{
//			*i=r1;
//			*j=r2;
//			r=3;
//		}
//		else
//		{
//			*i=r1;
//			r=1;
//
//		}
//	}
//	else
//	{
//		if(hash_table_query(hash_table_l,kmer_current1,&r2,v_nrp,para)!=0)
//		{
//			*j=r2;
//			r=2;
//
//		}
//		else
//		{
//			r=0;
//		}
//	}
	return r;
}

 void Gen_MNRpath(char* p_nrp,char** p,uint64_t nrpath_len,uint64_t* p_len,uint32_t kmer_len)
{
	bit256KmerPara para;
	get_para(&para,kmer_len);

	vector <string> v_nrp;
	v_nrp.clear();
	string nrp_tmp;
	nrp_tmp.clear();

	uint64_t kmer_current[4];
	uint64_t kmer_current1[4];

	for(uint32_t i=0;i<4;i++)
	{
		kmer_current[i]=0;
		kmer_current1[i]=0;
	}

	//1 splite the nrp into multiple strings
	uint64_t total_len1=0;
	for(uint64_t i=0;i<nrpath_len;i++)
	{
		if(p_nrp[i]=='$')
		{
			if(nrp_tmp.empty())
			{
				cout << "error:the string is an empty string!"<< endl;
			}
			else
			{
				v_nrp.push_back(nrp_tmp);
				total_len1+=nrp_tmp.size();
				nrp_tmp.clear();
			}
		}
		else
		{
			nrp_tmp.push_back(p_nrp[i]);
		}
	}
	if(nrp_tmp.size()!=0)
	{
		v_nrp.push_back(nrp_tmp);
		total_len1+=nrp_tmp.size();
		nrp_tmp.clear();
	}

	cout << "the size of v_nrp is: "<<v_nrp.size() << endl;
	cout << "the total length of strings in v_nrp is: "<< total_len1 << endl;
	//2 free the nrp space
	free(p_nrp);

	//3 merger1
	vector <string> v_nrp1;
	v_nrp1.clear();

	vector< vector<struct ikmer> > hash_table_l;
	vector< vector<struct ikmer> > hash_table_r;

	hash_table_initial(hash_table_l);
	hash_table_initial(hash_table_r);

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);


	for(uint32_t i=0;i<v_nrp.size();i++)
	{
		if(v_nrp1.size()==0)
		{
			v_nrp1.push_back(v_nrp[i]);
			//add v_nrp[i] into the hash table
			cal_kmer(v_nrp[i], para, 0,kmer_current);
//			hash_table_insert(hash_table_l,kmer_current,v_nrp1.size()-1,para);
			struct ikmer a;
			a.i=v_nrp1.size()-1;
			a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
			for(uint32_t x=0;x<para.kmer64Len;x++)
			{
				a.kmer[x]=kmer_current[x];
			}
			hash_table_l[hash_function(kmer_current,para)].push_back(a);
			cal_kmer(v_nrp[i], para, 1,kmer_current);
//			hash_table_insert(hash_table_r,kmer_current,v_nrp1.size()-1,para);
			a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
			for(uint32_t x=0;x<para.kmer64Len;x++)
			{
				a.kmer[x]=kmer_current[x];
			}
			hash_table_r[hash_function(kmer_current,para)].push_back(a);

		}
		else
		{
			cal_kmer(v_nrp[i], para, 0,kmer_current);
			cal_kmer(v_nrp[i], para, 1,kmer_current1);

			//for a new string
			uint64_t x,y;
			uint32_t r,l1,l2;
//			r=find_lr(v_nrp[i],v_nrp1,hash_table_l,hash_table_r,&x, &y,para);
			uint32_t value;
			value=hash_function(kmer_current,para);

			l1=0;
			if(hash_table_r[value].size()!=0)
			{
				unsigned j;
				for(j=0;j<hash_table_r[value].size();j++)
				{
					if(cmp(hash_table_r[value][j].kmer,kmer_current,para.kmer64Len)==1)
					{
						if(!v_nrp1[hash_table_r[value][j].i].empty())
						{
							x=hash_table_r[value][j].i;
							l1=1;
							break;
						}
					}
				}
			}

			value=hash_function(kmer_current1,para);

			l2=0;
			if(hash_table_l[value].size()!=0)
			{
				unsigned j;
				for(j=0;j<hash_table_l[value].size();j++)
				{
					if(cmp(hash_table_l[value][j].kmer,kmer_current1,para.kmer64Len)==1)
					{
						if(!v_nrp1[hash_table_l[value][j].i].empty())
						{
							if(l1!=0)
							{
								if(x!=hash_table_l[value][j].i)
								{
									y=hash_table_l[value][j].i;
									l2=1;
									break;
								}
							}
							else
							{
								y=hash_table_l[value][j].i;
								l2=1;
								break;
							}
						}
					}
				}
			}

			if(l1!=0)
			{
				if(l2!=0)
				{
					r=3;
				}
				else
				{
					r=1;

				}
			}
			else
			{
				if(l2!=0)
				{
					r=2;

				}
				else
				{
					r=0;
				}
			}

//			if(i==53345704)
//			{
//				cout << "r:" << r << endl;
//				cout << "v_nrp[i]:" << v_nrp[i] << endl;
//				cout << "x:" << x << endl;
//				if(r==1||r==3)
//				{
//					cout << "v_nrp1[x]:" << v_nrp1[x] << endl;
//				}
//				cout << "y:" << y << endl;
//				if(r==2||r==3)
//				{
//					cout << "v_nrp1[y]:" << v_nrp1[y] << endl;
//				}
//				cout << "53345704:" << v_nrp1[53345704] << endl;
//
//			}
//			r=0;
			if(r==0)
			{
				v_nrp1.push_back(v_nrp[i]);
				//add v_nrp[i] into the hash table
				cal_kmer(v_nrp[i], para, 0,kmer_current);
//				hash_table_insert(hash_table_l,kmer_current,v_nrp1.size()-1,para);
				struct ikmer a;
				a.i=v_nrp1.size()-1;
				a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t x=0;x<para.kmer64Len;x++)
				{
					a.kmer[x]=kmer_current[x];
				}
				hash_table_l[hash_function(kmer_current,para)].push_back(a);

				cal_kmer(v_nrp[i], para, 1,kmer_current);
//				hash_table_insert(hash_table_r,kmer_current,v_nrp1.size()-1,para);
				a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t x=0;x<para.kmer64Len;x++)
				{
					a.kmer[x]=kmer_current[x];
				}
				hash_table_r[hash_function(kmer_current,para)].push_back(a);

			}
			else if(r==1)
			{
				string tmp;
				tmp.clear();
				tmp=v_nrp1[x];
				tmp=tmp+v_nrp[i].substr(kmer_len);
//				cal_kmer(v_nrp1[x], para, 0,kmer_current);
//				hash_table_delete(hash_table_l,kmer_current,x,para);
//				cal_kmer(v_nrp1[x], para, 1,kmer_current);
//				hash_table_delete(hash_table_r,kmer_current,x,para);
				v_nrp1[x].clear();

				v_nrp1.push_back(tmp);
				//add v_nrp[i] into the hash table
				cal_kmer(tmp, para, 0,kmer_current);
//				hash_table_insert(hash_table_l,kmer_current,v_nrp1.size()-1,para);
				struct ikmer a;
				a.i=v_nrp1.size()-1;
				a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t x=0;x<para.kmer64Len;x++)
				{
					a.kmer[x]=kmer_current[x];
				}
				hash_table_l[hash_function(kmer_current,para)].push_back(a);

				cal_kmer(tmp, para, 1,kmer_current);
//				hash_table_insert(hash_table_r,kmer_current,v_nrp1.size()-1,para);
				a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t x=0;x<para.kmer64Len;x++)
				{
					a.kmer[x]=kmer_current[x];
				}
				hash_table_r[hash_function(kmer_current,para)].push_back(a);
			}
			else if(r==2)
			{
				string tmp;
				tmp.clear();
				tmp=v_nrp[i];
				tmp=tmp+v_nrp1[y].substr(kmer_len);
//				cal_kmer(v_nrp1[y], para, 0,kmer_current);
//				hash_table_delete(hash_table_l,kmer_current,y,para);
//				cal_kmer(v_nrp1[y], para, 1,kmer_current);
//				hash_table_delete(hash_table_r,kmer_current,y,para);
				v_nrp1[y].clear();

				v_nrp1.push_back(tmp);
				//add v_nrp[i] into the hash table
				cal_kmer(tmp, para, 0,kmer_current);
//				hash_table_insert(hash_table_l,kmer_current,v_nrp1.size()-1,para);
				struct ikmer a;
				a.i=v_nrp1.size()-1;
				a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t x=0;x<para.kmer64Len;x++)
				{
					a.kmer[x]=kmer_current[x];
				}
				hash_table_l[hash_function(kmer_current,para)].push_back(a);

				cal_kmer(tmp, para, 1,kmer_current);
//				hash_table_insert(hash_table_r,kmer_current,v_nrp1.size()-1,para);
				a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
				for(uint32_t x=0;x<para.kmer64Len;x++)
				{
					a.kmer[x]=kmer_current[x];
				}
				hash_table_r[hash_function(kmer_current,para)].push_back(a);
			}
			else if(r==3)
			{
				if(x!=y)
				{
					string tmp;
					tmp.clear();
					tmp=v_nrp1[x];
					tmp=tmp+v_nrp[i].substr(kmer_len);
					tmp=tmp+v_nrp1[y].substr(kmer_len);
	//				cal_kmer(v_nrp1[x], para, 0,kmer_current);
	//				hash_table_delete(hash_table_l,kmer_current,x,para);
	//				cal_kmer(v_nrp1[x], para, 1,kmer_current);
	//				hash_table_delete(hash_table_r,kmer_current,x,para);
	//				cal_kmer(v_nrp1[y], para, 0,kmer_current);
	//				hash_table_delete(hash_table_l,kmer_current,y,para);
	//				cal_kmer(v_nrp1[y], para, 1,kmer_current);
	//				hash_table_delete(hash_table_r,kmer_current,y,para);
					v_nrp1[x].clear();
					v_nrp1[y].clear();

					v_nrp1.push_back(tmp);
					//add v_nrp[i] into the hash table
					cal_kmer(tmp, para, 0,kmer_current);
	//				hash_table_insert(hash_table_l,kmer_current,v_nrp1.size()-1,para);
					struct ikmer a;
					a.i=v_nrp1.size()-1;
					a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
					for(uint32_t x=0;x<para.kmer64Len;x++)
					{
						a.kmer[x]=kmer_current[x];
					}
					hash_table_l[hash_function(kmer_current,para)].push_back(a);

					cal_kmer(tmp, para, 1,kmer_current);
	//				hash_table_insert(hash_table_r,kmer_current,v_nrp1.size()-1,para);
					a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
					for(uint32_t x=0;x<para.kmer64Len;x++)
					{
						a.kmer[x]=kmer_current[x];
					}
					hash_table_r[hash_function(kmer_current,para)].push_back(a);
				}
				else
				{
					string tmp;
					tmp.clear();
					tmp=v_nrp[i];
					tmp=tmp+v_nrp1[y].substr(kmer_len);
	//				cal_kmer(v_nrp1[y], para, 0,kmer_current);
	//				hash_table_delete(hash_table_l,kmer_current,y,para);
	//				cal_kmer(v_nrp1[y], para, 1,kmer_current);
	//				hash_table_delete(hash_table_r,kmer_current,y,para);
					v_nrp1[y].clear();

					v_nrp1.push_back(tmp);
					//add v_nrp[i] into the hash table
					cal_kmer(tmp, para, 0,kmer_current);
	//				hash_table_insert(hash_table_l,kmer_current,v_nrp1.size()-1,para);
					struct ikmer a;
					a.i=v_nrp1.size()-1;
					a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
					for(uint32_t x=0;x<para.kmer64Len;x++)
					{
						a.kmer[x]=kmer_current[x];
					}
					hash_table_l[hash_function(kmer_current,para)].push_back(a);

					cal_kmer(tmp, para, 1,kmer_current);
	//				hash_table_insert(hash_table_r,kmer_current,v_nrp1.size()-1,para);
					a.kmer=(uint64_t *)malloc(sizeof(uint64_t)*para.kmer64Len);
					for(uint32_t x=0;x<para.kmer64Len;x++)
					{
						a.kmer[x]=kmer_current[x];
					}
					hash_table_r[hash_function(kmer_current,para)].push_back(a);
				}
			}
		}
	}
	hash_table_destroy(hash_table_l);
	hash_table_destroy(hash_table_r);
	v_nrp.clear();

	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"the second step time cost is: "<<span<<endl;

	gettimeofday(&tvs,NULL);

	vector <string> v_nrp2;
	v_nrp2.clear();

	uint64_t v1=0;
	for(uint32_t i=0;i<v_nrp1.size();i++)
	{
		if(v_nrp1[i].size()!=0)
		{
			v1++;
			cal_kmer(v_nrp1[i], para, 0,kmer_current);
			cal_kmer(v_nrp1[i], para, 1,kmer_current1);
			if(cmp(kmer_current,kmer_current1,para.kmer64Len)==1)
			{
				v_nrp2.push_back(v_nrp1[i]);
				v_nrp1[i].clear();
			}
		}
	}

	cout << "the size of v_nrp1 is: "<<v1 << endl;
	cout << "the size of v_nrp2 is: "<<v_nrp2.size() << endl;
//	cout << "hrhr" << endl;

	uint32_t maxl=0;
	for(uint32_t i=0;i<v_nrp1.size();i++)
	{
		if(v_nrp1[i].size()>maxl)
		{
			maxl=v_nrp1[i].size();
		}
	}
	if(v_nrp2.size()!=0)
	{
		for(uint32_t i=0;i<v_nrp2.size();i++)
		{
			if(v_nrp2[i].size()>maxl)
			{
				maxl=v_nrp2[i].size();
			}
		}
	}

	uint32_t x_tmp;
	x_tmp=0;
	if(v_nrp2.size()!=0)
	{
		struct bit256Hash* p_root=NULL;
		p_root=initial256BitHashTable();
		char* p=(char*)malloc(sizeof(char)*maxl);

		for(uint64_t i=0;i<v_nrp1.size();i++)
		{
			if(v_nrp1[i].size()!=0)
			{
				for(uint32_t j=0;j<v_nrp1[i].size();j++)
				{
					p[j]=v_nrp1[i][j];
				}
				cal_hash_value_directly_256bit(p,kmer_current,para);
				insert256BitHashTable(p_root,kmer_current,i+1,para);
				if(kmer_current[0]==40845413009)
				{
					x_tmp++;
					cout << "x_tmp1:"<< x_tmp << endl;
					uint64_t index;
					index=search256BitHashTable(p_root,kmer_current,para);
					cout << "index:" << index << endl;
				}

				for(uint32_t j=1;j<=v_nrp1[i].size()-kmer_len;j++)
				{
					cal_hash_value_indirectly_256bit(p+j,kmer_current,kmer_current,para);
					insert256BitHashTable(p_root,kmer_current,i+1,para);
					if(kmer_current[0]==40845413009)
					{
						x_tmp++;
						cout << "x_tmp2:"<< x_tmp << endl;
						uint64_t index;
						index=search256BitHashTable(p_root,kmer_current,para);
						cout << "index:" << index << endl;
					}
				}
			}
		}

		uint32_t label_first=0;
		while(1)
		{
			uint32_t ll=0;
			for(uint32_t i=0;i<v_nrp2.size();i++)
			{
				if(v_nrp2[i].size()!=0)
				{
//					cout << v_nrp2[i] << endl;
					cal_kmer(v_nrp2[i], para, 0,kmer_current);
					// find the target string in v_nrp1
					uint64_t index;
					index=search256BitHashTable(p_root,kmer_current,para);

					if(index!=0)
					{
						index--;
						// change the string in v_nrp1
			//			pos=v_nrp1[index].find(v_nrp2[i].substr(0,kmer_len-1));

						if(label_first==0)
						{
							uint32_t pos=0;
							for(uint32_t j=0;j<v_nrp1[index].size();j++)
							{
								p[j]=v_nrp1[index][j];
							}
							cal_hash_value_directly_256bit(p,kmer_current1,para);

							for(uint32_t j=1;j<v_nrp1[index].size()-kmer_len;j++)
							{
								cal_hash_value_indirectly_256bit(p+j,kmer_current1,kmer_current1,para);
								if(cmp(kmer_current,kmer_current1,para.kmer64Len)==1)
								{
									pos=j;
									break;
								}
							}
							if(pos==0)
							{
								cout << "error:no pos is found!" << endl;
								cout << v_nrp1[index] << endl;
								ll=1;
							}
							else
							{
//								cout << "pos is found!" << endl;
								string tmp;
								tmp.clear();
								tmp=v_nrp1[index].substr(0,pos+kmer_len);
								tmp=tmp+v_nrp2[i].substr(kmer_len);
								tmp=tmp+v_nrp1[index].substr(pos+kmer_len);
								v_nrp1[index].clear();
								v_nrp1[index]=tmp;

								if(tmp.size()>maxl)
								{
									maxl=tmp.size();
									p=(char*)realloc(p,sizeof(char)*maxl);
								}

								for(uint32_t j=0;j<v_nrp2[i].size();j++)
								{
									p[j]=v_nrp2[i][j];
								}
								cal_hash_value_directly_256bit(p,kmer_current,para);

								for(uint32_t j=1;j<v_nrp2[i].size()-kmer_len;j++)
								{
									cal_hash_value_indirectly_256bit(p+j,kmer_current,kmer_current,para);
									insert256BitHashTable(p_root,kmer_current,index+1,para);
									if(kmer_current[0]==40845413009)
									{
										x_tmp++;
										cout << "x_tmp3:"<< x_tmp << endl;
										uint64_t index;
										index=search256BitHashTable(p_root,kmer_current,para);
										cout << "index:" << index << endl;
									}
								}
								v_nrp2[i].clear();
							}
						}
						else
						{
							int32_t pos=-1;
							for(uint32_t j=0;j<v_nrp1[index].size();j++)
							{
								p[j]=v_nrp1[index][j];
							}
							cal_hash_value_directly_256bit(p,kmer_current1,para);
							int32_t j=0;
							if(cmp(kmer_current,kmer_current1,para.kmer64Len)==1)
							{
								pos=j;
							}
							else
							{
								for(j=1;j<=v_nrp1[index].size()-kmer_len;j++)
								{
									cal_hash_value_indirectly_256bit(p+j,kmer_current1,kmer_current1,para);
									if(cmp(kmer_current,kmer_current1,para.kmer64Len)==1)
									{
										pos=j;
										break;
									}
								}
							}

							if(pos==-1)
							{
								cout << "error:no pos is found!" << endl;
								cout << v_nrp1[index] << endl;
								ll=1;
							}
							else
							{
								cout << "pos is found!" << endl;
								string tmp;
								tmp.clear();
								tmp=v_nrp1[index].substr(0,pos+kmer_len);
								tmp=tmp+v_nrp2[i].substr(kmer_len);
								tmp=tmp+v_nrp1[index].substr(pos+kmer_len);
								v_nrp1[index].clear();
								v_nrp1[index]=tmp;

								if(tmp.size()>maxl)
								{
									maxl=tmp.size();
									p=(char*)realloc(p,sizeof(char)*maxl);
								}

								for(uint32_t j=0;j<v_nrp2[i].size();j++)
								{
									p[j]=v_nrp2[i][j];
								}
								cal_hash_value_directly_256bit(p,kmer_current,para);

								for(uint32_t j=1;j<v_nrp2[i].size()-kmer_len;j++)
								{
									cal_hash_value_indirectly_256bit(p+j,kmer_current,kmer_current,para);
									insert256BitHashTable(p_root,kmer_current,index+1,para);
									if(kmer_current[0]==40845413009)
									{
										x_tmp++;
										cout << "x_tmp4:"<< x_tmp << endl;
										uint64_t index;
										index=search256BitHashTable(p_root,kmer_current,para);
										cout << "index:" << index << endl;
									}
								}
								v_nrp2[i].clear();
							}
						}
					}
					else
					{
						cout << "error: string is not found!" << endl;
						cout << v_nrp2[i] << endl;
						for(uint32_t j=0;j<v_nrp2[i].size();j++)
						{
							p[j]=v_nrp2[i][j];
						}
						cal_hash_value_directly_256bit(p,kmer_current,para);

						for(uint32_t j=1;j<=v_nrp2[i].size()-kmer_len;j++)
						{
							cal_hash_value_indirectly_256bit(p+j,kmer_current,kmer_current,para);
							index=search256BitHashTable(p_root,kmer_current,para);
							if(index!=0)
							{
								cout << "sloved!" << endl;
								string tmp;
								tmp=v_nrp2[i].substr(j);
								tmp=tmp+v_nrp2[i].substr(kmer_len,j);
								v_nrp2[i]=tmp;
								break;
							}
							else
							{
								cout << "error:not find the inserted node!" << endl;
//								char a;
//								cin >>a;
							}
						}
						ll=1;
					}
				}
			}
			if(ll==0)
			{
				break;
			}
			label_first++;
		}

		free(p);
		v_nrp2.clear();
	}

	uint64_t total_len;
	total_len=0;
	for(uint32_t i=0;i<v_nrp1.size();i++)
	{
		if(v_nrp1[i].size()!=0)
		{
			total_len=total_len+v_nrp1[i].size()+1;
		}
	}
	char * pp;
	pp=(char*)malloc(sizeof(char)*(total_len+1));
	uint64_t jj=0;
	for(uint32_t i=0;i<v_nrp1.size();i++)
	{
		if(v_nrp1[i].size()!=0)
		{
			for(uint32_t j=0;j<v_nrp1[i].size();j++)
			{
				pp[jj]=v_nrp1[i][j];
				jj++;
			}
			pp[jj]='?';
			jj++;
		}
	}
//	pp[jj]='\0';
	*p=pp;
	*p_len=total_len;
	v_nrp1.clear();
	gettimeofday(&tve,NULL);
	span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"the third step time cost is: "<<span<<endl;
}

void Total_Gen_MNRpath(char * pathFile,uint32_t kmer_len,char** nrpath,uint64_t * nrpath_len)
{
	char* nrpath_tmp;
	uint64_t nrpath_len_tmp;

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	cout <<"start..."<<endl;

	Gen_NRpaths(pathFile,kmer_len,&nrpath_tmp,&nrpath_len_tmp);
//	ofstream out_nrp;
//	out_nrp.open("nrp");
//	for(uint32_t i=0;i<nrpath_len_tmp;i++)
//	{
//		out_nrp << nrpath_tmp[i];
//	}

	cout << "end..."<< endl;
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"the first step time cost is: "<<span<<endl;

	gettimeofday(&tvs,NULL);
	cout <<"start..."<<endl;

	Gen_MNRpath(nrpath_tmp,nrpath,nrpath_len_tmp,nrpath_len,kmer_len);

	char *file;
	file=(char*)malloc(sizeof(char)*6);
	file[0]='a';
	file[1]='a';
	char swap[4];
	sprintf(swap,"%d",kmer_len);
	strcpy(file+2,swap);

	FILE *fp;
	fp=fopen(file,"wb+");
	fwrite(*nrpath,sizeof(char),*nrpath_len,fp);
	fclose(fp);
	free(file);

	cout << "end..."<< endl;
	gettimeofday(&tve,NULL);
	span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"the second and third step time cost is: "<<span<<endl;
}

void Total_Gen_MNRpath_from_bigger_one(char * pathFile,uint32_t kmer_len,char** nrpath,uint64_t * nrpath_len)
{
	char* nrpath_tmp;
	uint64_t nrpath_len_tmp;

	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	cout <<"start..."<<endl;

	Gen_NRpaths_from_biger_one(pathFile,kmer_len,&nrpath_tmp,&nrpath_len_tmp);
//	ofstream out_nrp;
//	out_nrp.open("nrp");
//	for(uint32_t i=0;i<nrpath_len_tmp;i++)
//	{
//		out_nrp << nrpath_tmp[i];
//	}

	cout << "end..."<< endl;
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"the first step time cost is: "<<span<<endl;

	gettimeofday(&tvs,NULL);
	cout <<"start..."<<endl;

	Gen_MNRpath(nrpath_tmp,nrpath,nrpath_len_tmp,nrpath_len,kmer_len);

	char *file;
	file=(char*)malloc(sizeof(char)*6);
	file[0]='a';
	file[1]='a';
	char swap[4];
	sprintf(swap,"%d",kmer_len);
	strcpy(file+2,swap);

	FILE *fp;
	fp=fopen(file,"wb+");
	fwrite(*nrpath,sizeof(char),*nrpath_len,fp);
	fclose(fp);
	free(file);

	cout << "end..."<< endl;
	gettimeofday(&tve,NULL);
	span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout <<"the second and third step time cost is: "<<span<<endl;
}
