/*
 * mappers.cpp
 *
 *  Created on: May 8, 2021
 *      Author: bio
 */
#include "mappers.h"
extern uint32_t num_read;
uint64_t candidate_num = 0;

void* Parall_nativeMapper(void *p) {

	struct para_multi_p *tmp = (struct para_multi_p*) p;
	for (uint32_t i = tmp->start; i <= tmp->end; i++) {
		for (uint32_t j = 0; j <= tmp->p->e; j++) {
			uint32_t seed_id_tmp =
					tmp->p->read_sav_seed[i * (tmp->p->e + 1) + j];
			if (tmp->p->VS_message[seed_id_tmp].saInterval[1]
					>= tmp->p->VS_message[seed_id_tmp].saInterval[0]) {
				pthread_mutex_lock(tmp->m);
				//候选位置数量有几个
				uint32_t saInterval_size =
						tmp->p->VS_message[seed_id_tmp].saInterval[1]
								- tmp->p->VS_message[seed_id_tmp].saInterval[0]
								+ 1;

				//全部read的候选位置数量
				candidate_num += saInterval_size;
				pthread_mutex_unlock(tmp->m);
				//得出种子在一段read上的位置
				uint32_t seed_pos = get_seed_pos(tmp->p->len_read, tmp->p->e,
						j);

				//遍历每个候选位置
				for (uint32_t k = 0; k < saInterval_size; k++) {
					uint32_t can_pos = 0;
					if (tmp->p->VS_message[seed_id_tmp].seed_sa[k]
							< seed_pos + tmp->p->e) {
						can_pos = 0;
					} else {
						can_pos = tmp->p->VS_message[seed_id_tmp].seed_sa[k]
								- (seed_pos + tmp->p->e);
					}
					int r;
					int mes;
					//如果任务read数目大于实际的read数目，输出错误的范围
					if (tmp->p->p_reads + i * tmp->p->len_read
							>= tmp->p->p_reads + tmp->p->len_read * num_read) {
						cout << "Wrong Range" << endl;
					}

					//第一个参数参考基因候选位置，第二个参数read，第三个参数是阈值
					r = banded_edit_distance(tmp->p->p_ref + can_pos,
							tmp->p->p_reads + i * tmp->p->len_read,
							tmp->p->len_read, tmp->p->e, &mes);
				}
			}
		}
	}
}
void naiveMapper(struct naiveMapper_para p) {
	cout << "p.task_size:" << p.task_size << endl;
	cout << "p.len_read:" << p.len_read << endl;
	cout << "num_read:" << num_read << endl;
	if (p.thread_num == 1) {
		for (uint32_t i = 0; i < p.task_size; i++) {
			for (uint32_t j = 0; j < p.e + 1; j++) {
				uint32_t seed_id_tmp = p.read_sav_seed[i * (p.e + 1) + j];
				if (p.VS_message[seed_id_tmp].saInterval[1]
						>= p.VS_message[seed_id_tmp].saInterval[0])
				{
					uint32_t saInterval_size =p.VS_message[seed_id_tmp].saInterval[1]
									- p.VS_message[seed_id_tmp].saInterval[0]
									+ 1;

					candidate_num += saInterval_size;
					uint32_t seed_pos = get_seed_pos(p.len_read, p.e, j);
					for (uint32_t k = 0; k < saInterval_size; k++) {
						uint32_t can_pos = 0;
						if (p.VS_message[seed_id_tmp].seed_sa[k]
								< seed_pos + p.e)
						{
							can_pos = 0;
						}
						else
						{
							can_pos = p.VS_message[seed_id_tmp].seed_sa[k]
									- (seed_pos + p.e);
						}

						if (p.VS_message[seed_id_tmp].seed_sa[k]+128
								> p.ref_length)
						{
							can_pos = p.ref_length-128;
						}
						int r;
						int mes;
						if (p.p_reads + i * p.len_read
								>= p.p_reads + p.len_read * num_read) {
							cout << "Wrong Range!" << endl;
						}
						r = banded_edit_distance(p.p_ref + can_pos,
								p.p_reads + i * p.len_read, p.len_read, p.e,
								&mes);

						if(!pre_alignment(p.p_reads, p.len_read, p.p_ref, can_pos, can_pos+2*p.e, p.e))
						{
							cerr << i << endl;
						}
					}
				}
			}
		}
	} else {
		cout << "多线程开始" << endl;

		//多线程
		pthread_t *t;
		t = (pthread_t *) malloc(sizeof(pthread_t) * p.thread_num);
		pthread_mutex_t m;
		pthread_mutex_init(&m,NULL);

		//计算每个线程要处理的read数目
		uint32_t r = p.task_size / p.thread_num;
		uint32_t s = p.task_size % p.thread_num;

		//计算第几个read
		uint32_t cur = 0;

		//构造多线程带的参数
		struct para_multi_p *para;
		para = (struct para_multi_p*) malloc(
				sizeof(struct para_multi_p) * p.thread_num);

		for (uint32_t i = 0; i < p.thread_num; i++)
		{
			para[i].start = cur;
			if (i < s) {
				cur = cur + r + 1;
			} else {
				cur = cur + r;
			}
			para[i].end = cur - 1;
			para[i].p = &p;
			para[i].m = &m;

			if (pthread_create(t + i, NULL, Parall_nativeMapper,
					(void*) (para + i)) != 0) {
				cout << "error!" << endl;
			}
		}
		for (uint32_t i = 0; i < p.thread_num; i++) {
			pthread_join(t[i], NULL);
		}
		free(para);
		free(t);
	}
	cout << "candidate_num:" << candidate_num << endl;
}

