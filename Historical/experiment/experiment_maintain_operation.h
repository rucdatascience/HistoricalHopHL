#pragma once
#include "Historical/graph_with_time_span/graph.h"
#include <boost/heap/fibonacci_heap.hpp>
#include "Historical/graph_with_time_span/two_hop_label.h"
#include <algorithm>
#include <mutex>
#include <map>
#include <CPU/tool_functions/ThreadPool.h>
#include <Historical/experiment/experiment_operation.h>
#include <thread>

namespace experiment
{
	namespace nonhop
	{
		class pair_label
		{ // pair_label2 is stored in NoP
		public:
			int first, second;
			pair_label(int _first, int _second)
			{
				first = _first;
				second = _second;
			}
			bool operator==(const pair_label other) const
			{
				return (first == other.first && second == other.second);
			}
			bool operator<(const pair_label other) const
			{ // used to sort/search pair_label2 in set
				if (first != other.first)
					return first < other.first;
				return second < other.second;
			}
		};

		struct node_for_DIFFUSE
		{
			int index;
			int disx;
			node_for_DIFFUSE() {}
			node_for_DIFFUSE(int _u, int _dis)
			{
				index = _u;
				disx = _dis;
			}
		}; // define the node in the queue
		bool operator<(node_for_DIFFUSE const &x, node_for_DIFFUSE const &y)
		{
			return x.disx > y.disx; // < is the max-heap; > is the min heap
		}
		class affected_label
		{
		public:
			int first, second;
			int dis;
			int t_s, t_e;
			affected_label() {}
			affected_label(int _first, int _second, int _dis)
			{
				first = _first;
				second = _second;
				dis = _dis;
			}
		};

		typedef typename boost::heap::fibonacci_heap<node_for_DIFFUSE>::handle_type handle_t_for_DIFFUSE;
		std::vector<std::vector<std::pair<int, int>>> Dis;
		std::vector<std::vector<int>> Q_value;
		std::vector<std::vector<handle_t_for_DIFFUSE>> Q_handles;

		std::shared_mutex mtx_595_1, mtx_595_2;
		std::vector<std::shared_mutex> mtx_5952(max_N_ID_for_mtx_595);

		void initialize_experiment_global_values_dynamic(int N, int thread_num)
		{
			Dis.resize(thread_num);
			Q_value.resize(thread_num);
			Q_handles.resize(thread_num);
			std::queue<int>().swap(Qid_595);
			for (int i = 0; i < thread_num; i++)
			{
				Dis[i].resize(N, {-1, -1});
				Q_value[i].resize(N, 1e7);
				Q_handles[i].resize(N);
				Qid_595.push(i);
			}
		};
		namespace ruc
		{
			namespace decrease
			{
				void decrease_maintain_step1_batch(std::map<std::pair<int, int>, int> &v_map, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> *CL,
												   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
				{
					for (auto it : v_map)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, CL]
																		  {
								int v1 = it.first.first, v2 = it.first.second;
								int w_new = it.second;
								for (int sl = 0; sl < 2; sl++) {
									if (sl == 1) {
										std::swap(v1, v2);
									}
									for (auto& it : (*L)[v1]) {
										if (it.vertex <= v2 && it.distance + w_new < 2e6 && it.t_e == std::numeric_limits<int>::max()) {
											auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current(*L, it.vertex, v2); // query_result is {distance, common hub}
											if (query_result.first > it.distance + w_new) {
												mtx_595_1.lock();
												CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
												mtx_595_1.unlock();
											}
											else {
												auto search_result = search_sorted_two_hop_label_weight_in_current((*L)[v2], it.vertex);
												// TODO-GPY MAX_VALUE PRUNE 1e7
												if (search_result > it.distance + w_new && search_result < 1e7) {
													mtx_595_1.lock();
													CL->push_back(affected_label{ v2, it.vertex, it.distance + w_new });
													mtx_595_1.unlock();
												}
												if (query_result.second != it.vertex) {
													mtx_5952[v2].lock();
													PPR_TYPE::PPR_insert(*PPR, v2, query_result.second, it.vertex);
													mtx_5952[v2].unlock();
												}
												if (query_result.second != v2) {
													mtx_5952[it.vertex].lock();
													PPR_TYPE::PPR_insert(*PPR, it.vertex, query_result.second, v2);
													mtx_5952[it.vertex].unlock();
												}
											}
										}
									}
								}
								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				};

				void DIFFUSE_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> &CL,
								   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
				{

					// Deduplication
					std::map<std::pair<int, int>, int> CL_edge_map;
					for (auto &it : CL)
					{
						if (CL_edge_map.count({it.first, it.second}) == 0)
						{
							CL_edge_map[{it.first, it.second}] = it.dis;
						}
						else if (CL_edge_map[{it.first, it.second}] > it.dis)
						{
							CL_edge_map[{it.first, it.second}] = it.dis;
						}
					}

					// extract each unique hub v and its (u,dis) list
					std::map<int, std::vector<std::pair<int, int>>> CL_map; // CL_map[v]=(u1,dis1),(u2,dis2)...
					for (auto &it : CL_edge_map)
					{
						int u = it.first.first;
						int v = it.first.second;
						int dis = it.second;
						if (CL_map.count(v) == 0)
						{
							std::vector<std::pair<int, int>> vec_with_hub_v;
							vec_with_hub_v.emplace_back(std::make_pair(u, dis));
							CL_map[v] = vec_with_hub_v;
						}
						else
						{
							std::vector<std::pair<int, int>> vec_with_hub_v = CL_map[v];
							vec_with_hub_v.emplace_back(std::make_pair(u, dis));
							CL_map[v] = vec_with_hub_v;
						}
					}

					std::vector<std::pair<int, std::vector<std::pair<int, int>>>> CL_map_vec(CL_map.begin(), CL_map.end());
					sort(CL_map_vec.begin(), CL_map_vec.end(), [](const std::pair<int, std::vector<std::pair<int, int>>> &a, const std::pair<int, std::vector<std::pair<int, int>>> &b)
						 { return a.first < b.first; });

					// each thread processes one unique hub
					for (auto &it : CL_map_vec)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, &instance_graph, PPR]
																		  {

								mtx_595_1.lock();
								int current_tid = Qid_595.front();
								Qid_595.pop();
								mtx_595_1.unlock();

								int v = it.first;
								std::vector<std::pair<int, int>> vec_with_hub_v = it.second;

								// int u = it.first.first, v = it.first.second;
								// weightTYPE du = it.second;
								mtx_595[v].lock_shared();
								auto Lv = (*L)[v]; // to avoid interlocking
								mtx_595[v].unlock_shared();

								std::vector<int> Dis_changed;
								auto& DIS = Dis[current_tid];
								auto& Q_HANDLES = Q_handles[current_tid];
								auto& Q_VALUE = Q_value[current_tid];
								// int N=instance_graph.size();
								// std::vector<std::pair<weightTYPE, int>> DIS(N, { -1, -1 });
								// std::vector<handle_t_for_DIFFUSE> Q_HANDLES(N);
								// std::vector<weightTYPE> Q_VALUE(N, MAX_VALUE);
								boost::heap::fibonacci_heap<node_for_DIFFUSE> Q;

								for (auto& it : vec_with_hub_v) {
									int u = it.first;
									int du = it.second;
									DIS[u] = { du, v }; // <distance, hub responsible for this distance>
									Dis_changed.push_back(u);
									Q_HANDLES[u] = Q.push(node_for_DIFFUSE(u, du));
									Q_VALUE[u] = du;
								}

								while (!Q.empty()) {

									node_for_DIFFUSE temp2 = Q.top();
									int x = temp2.index;
									int dx = temp2.disx;
									Q.pop();
									Q_VALUE[x] = 1e7;

									mtx_595[x].lock_shared();
									long long int d_old = search_sorted_two_hop_label_weight_and_hub_in_current((*L)[x], v).first;
									mtx_595[x].unlock_shared();
									if (d_old > dx) {
										mtx_595[x].lock();
										insert_sorted_two_hop_label((*L)[x], v, dx, t);
										mtx_595[x].unlock();
									}
									else {
										continue;
										//dx=d_old;
									}

									for (auto& nei : instance_graph[x]) {
										int xnei = nei.first;
										int d_new = dx + nei.second;

										if (v < xnei && d_new < 2e6) {
											if (DIS[xnei].first == -1) {
												mtx_595[xnei].lock_shared();
												DIS[xnei] = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[xnei], Lv);
												mtx_595[xnei].unlock_shared();
												Dis_changed.push_back(xnei);
											}
											if (DIS[xnei].first > d_new) {
												DIS[xnei] = { d_new, v };
												if (Q_VALUE[xnei] >= 1e7) {
													Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
												}
												else {
													Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
												}
												Q_VALUE[xnei] = d_new;
											}
											else {
												mtx_595[xnei].lock_shared();
												auto search_result = search_sorted_two_hop_label_weight_and_hub_in_current((*L)[xnei], v);
												mtx_595[xnei].unlock_shared();
												if (search_result.second != -1 && std::min(search_result.first, Q_VALUE[xnei]) > d_new) {
													if (Q_VALUE[xnei] >= 1e7) {
														Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
													}
													else {
														Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
													}
													Q_VALUE[xnei] = d_new;
												}
												if (DIS[xnei].second != v) {
													mtx_5952[xnei].lock();
													PPR_TYPE::PPR_insert(*PPR, xnei, DIS[xnei].second, v);
													mtx_5952[xnei].unlock();
												}
												if (DIS[xnei].second != xnei) {
													mtx_5952[v].lock();
													PPR_TYPE::PPR_insert(*PPR, v, DIS[xnei].second, xnei);
													mtx_5952[v].unlock();
												}
											}
										}
									}
								}

								for (int i : Dis_changed) {
									DIS[i] = { -1, -1 };
								}

								mtx_595_1.lock();
								Qid_595.push(current_tid);
								mtx_595_1.unlock();

								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				template <typename weight_type>
				void decrease_maintain(graph<weight_type> &instance_graph, two_hop_case_info &mm, std::vector<std::pair<int, int>> &v, std::vector<int> &w_new,
									   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{
					std::map<std::pair<int, int>, int> w_new_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_new_map.count(v[i]) == 0)
						{
							w_new_map[v[i]] = w_new[i];
						}
						else if (w_new_map[v[i]] > w_new[i])
						{
							w_new_map[v[i]] = w_new[i];
						}
					}
					std::vector<affected_label> CL;

					decrease_maintain_step1_batch(w_new_map, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic, time);

					DIFFUSE_batch(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic, time);
				}
			}
			namespace increase
			{
				void SPREAD1_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L,
								   std::vector<affected_label> &al1, std::vector<pair_label> *al2, std::map<std::pair<int, int>, weightTYPE> &w_old_map,
								   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{

					for (auto &it : al1)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, al2, &instance_graph, &w_old_map, time]
																		  {
							std::queue<std::pair<int, weightTYPE> > q; //(u,d)
							int v = it.second;
							q.push(std::pair<int, weightTYPE>(it.first, it.dis));
							while (!q.empty()) {
								int x = q.front().first;
								weightTYPE dx = q.front().second;
								q.pop();
								mtx_595[x].lock();
								insert_sorted_two_hop_label((*L)[x], v, MAX_VALUE, time);
								mtx_595[x].unlock();
								mtx_595_1.lock();
								al2->push_back(pair_label(x, v));
								mtx_595_1.unlock();
								for (auto nei : instance_graph[x]) {
									if (v < nei.first) {
										mtx_595[nei.first].lock();
										weightTYPE search_weight = search_sorted_two_hop_label_weight_in_current((*L)[nei.first], v);
										mtx_595[nei.first].unlock();
										weightTYPE w_old;
										if (w_old_map.count(std::pair<int, int>(x, nei.first)) > 0) {
											w_old = w_old_map[std::pair<int, int>(x, nei.first)];
										}
										else if (w_old_map.count(std::pair<int, int>(nei.first, x)) > 0) {
											w_old = w_old_map[std::pair<int, int>(nei.first, x)];
										}
										else {
											w_old = nei.second;
										}
										if (dx + w_old == search_weight && search_weight < MAX_VALUE) {
											q.push(std::pair<int, weightTYPE>(nei.first, dx + w_old));
										}
									}
								}
							}

							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void SPREAD2_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR,
								   std::vector<pair_label> &al2, std::vector<affected_label> *al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{

					for (auto &it : al2)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al3, &instance_graph]
																		  {

							int v = it.first, u = it.second;
							mtx_5952[v].lock();
							std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
							mtx_5952[v].unlock();
							temp.push_back(u);
							for (auto t : temp) {
								if (v < t) {
									long long d1 = MAX_VALUE;
									for (auto nei : instance_graph[t]) {
										mtx_595[nei.first].lock();
										d1 = std::min(d1, search_sorted_two_hop_label_weight_in_current((*L)[nei.first], v) + (long long)nei.second);
										mtx_595[nei.first].unlock();
									}
									auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current(*L, t, v);
									if (d1 >= 2e6) continue;
									if (query_result.first > d1) { // only add new label when it's absolutely necessary
										mtx_595_1.lock();
										al3->push_back(affected_label(t, v, d1));
										mtx_595_1.unlock();
									}
									else {
										if (query_result.second != v) {
											mtx_5952[t].lock();
											PPR_TYPE::PPR_insert(*PPR, t, query_result.second, v);
											mtx_5952[t].unlock();
										}
										if (query_result.second != t) {
											mtx_5952[v].lock();
											PPR_TYPE::PPR_insert(*PPR, v, query_result.second, t);
											mtx_5952[v].unlock();
										}
									}
								}
								else if (t < v) {
									long long d1 = MAX_VALUE;
									for (auto nei : instance_graph[v]) {
										mtx_595[nei.first].lock();
										d1 = std::min(d1, search_sorted_two_hop_label_weight_in_current((*L)[nei.first], t) + (long long)nei.second);
										mtx_595[nei.first].unlock();
									}
									if (d1 >= 2e6) continue;
									auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current(*L, v, t);
									if (query_result.first > d1) {
										mtx_595_1.lock();
										al3->push_back(affected_label(v, t, d1));
										mtx_595_1.unlock();
									}
									else {
										if (query_result.second != v) {
											mtx_5952[t].lock();
											PPR_TYPE::PPR_insert(*PPR, t, query_result.second, v);
											mtx_5952[t].unlock();
										}
										if (query_result.second != t) {
											mtx_5952[v].lock();
											PPR_TYPE::PPR_insert(*PPR, v, query_result.second, t);
											mtx_5952[v].unlock();
										}
									}
								}
							}

							return 1; }));
					}
					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void SPREAD3_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR, std::vector<affected_label> &al3,
								   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{

					// Deduplication (u,v,dis)
					std::map<std::pair<int, int>, weightTYPE> al3_edge_map;
					for (auto &it : al3)
					{
						if (al3_edge_map.count({it.first, it.second}) == 0)
						{
							al3_edge_map[{it.first, it.second}] = it.dis;
						}
						else if (al3_edge_map[{it.first, it.second}] > it.dis)
						{
							al3_edge_map[{it.first, it.second}] = it.dis;
						}
					}

					// extract each unique hub v and its (u,dis) list
					std::map<int, std::vector<std::pair<int, weightTYPE>>> al3_map; // al3_map[v]=(u1,dis1),(u2,dis2)...
					for (auto &it : al3_edge_map)
					{
						int u = it.first.first;
						int v = it.first.second;
						weightTYPE dis = it.second;
						if (al3_map.count(v) == 0)
						{
							std::vector<std::pair<int, weightTYPE>> vec_with_hub_v;
							vec_with_hub_v.emplace_back(std::make_pair(u, dis));
							al3_map[v] = vec_with_hub_v;
						}
						else
						{
							std::vector<std::pair<int, weightTYPE>> vec_with_hub_v = al3_map[v];
							vec_with_hub_v.emplace_back(std::make_pair(u, dis));
							al3_map[v] = vec_with_hub_v;
						}
					}

					std::vector<std::pair<int, std::vector<std::pair<int, weightTYPE>>>> al3_map_vec(al3_map.begin(), al3_map.end());
					sort(al3_map_vec.begin(), al3_map_vec.end(), [](const std::pair<int, std::vector<std::pair<int, weightTYPE>>> &a, const std::pair<int, std::vector<std::pair<int, weightTYPE>>> &b)
						 { return a.first < b.first; });

					// std::cout<<"SPREAD3_batch"<<std::endl;
					for (auto &it : al3_map_vec)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, &instance_graph, PPR, time]
																		  {

							mtx_595_1.lock();
							int current_tid = Qid_595.front();
							Qid_595.pop();
							mtx_595_1.unlock();
							// al3_map new---
							int v = it.first;
							std::vector<std::pair<int, weightTYPE>> vec_with_hub_v = it.second;

							// al3_map origin---
							// int u = it.first.first, v = it.first.second;
							// weightTYPE du = it.second;
							// al3 ---
							// int u = it.first, v = it.second;
							// weightTYPE du = it.dis;

							mtx_595[v].lock_shared();
							auto Lv = (*L)[v]; // to avoid interlocking
							mtx_595[v].unlock_shared();

							std::vector<int> Dis_changed;
							auto& DIS = Dis[current_tid];
							auto& Q_HANDLES = Q_handles[current_tid];
							auto& Q_VALUE = Q_value[current_tid];
							// int N = instance_graph.size();
							// std::vector<std::pair<weightTYPE, int>> DIS(N, { -1, -1 });
							// std::vector<handle_t_for_DIFFUSE> Q_HANDLES(N);
							// std::vector<weightTYPE> Q_VALUE(N, MAX_VALUE);
							boost::heap::fibonacci_heap<node_for_DIFFUSE> pq;

							for (auto& it : vec_with_hub_v) {
								int u = it.first;
								weightTYPE du = it.second;
								// std::cout<<"u: "<<u<<" v: "<<v<<" du: "<<du<<"\n";
								mtx_595[u].lock_shared();
								auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[u], Lv);
								mtx_595[u].unlock_shared();
								bool flag = false;
								if (query_result.first < du) {
									if (query_result.second != v) {
										mtx_5952[u].lock();
										PPR_TYPE::PPR_insert(*PPR, u, query_result.second, v);
										mtx_5952[u].unlock();
										flag = true;
									}
									if (query_result.second != u) {
										mtx_5952[v].lock();
										PPR_TYPE::PPR_insert(*PPR, v, query_result.second, u);
										mtx_5952[v].unlock();
										flag = true;
									}

									// mtx_595_1.lock();
									// Qid_595.push(current_tid);
									// mtx_595_1.unlock();
									// return 1;
								}

								if (flag == true) {
									continue;
								}

								DIS[u] = { du, v }; // <distance, hub responsible for this distance>
								Dis_changed.push_back(u);
								Q_HANDLES[u] = pq.push(node_for_DIFFUSE(u, du));
								Q_VALUE[u] = du;
							}

							while (!pq.empty()) {
								int x = pq.top().index;
								weightTYPE dx = pq.top().disx;
								pq.pop();
								Q_VALUE[x] = MAX_VALUE;

								mtx_595[x].lock();
								weightTYPE d_old = search_sorted_two_hop_label_weight_in_current((*L)[x], v);
								// std::cout<<"insert: "<<x<<" "<<v<<" "<<dx<<" "<<d_old<<"\n";
								if (dx < d_old) {
									insert_sorted_two_hop_label((*L)[x], v, dx, time);
								}
								else {
									mtx_595[x].unlock();
									continue;
									//dx = d_old;
								}
								mtx_595[x].unlock();

								for (auto nei : instance_graph[x]) {
									int xnei = nei.first;
									long long d_new = (long long)dx + nei.second;
									if (v < xnei && d_new < 2e6) {
										if (DIS[xnei].first == -1) {
											mtx_595[xnei].lock_shared();
											DIS[xnei] = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[xnei], Lv);
											mtx_595[xnei].unlock_shared();
											Dis_changed.push_back(xnei);
										}
										if (DIS[xnei].first > d_new) {
											DIS[xnei] = { d_new, v };
											if (Q_VALUE[xnei] >= MAX_VALUE) {
												Q_HANDLES[xnei] = pq.push(node_for_DIFFUSE(xnei, d_new));
											}
											else {
												pq.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
											}
											Q_VALUE[xnei] = d_new;
										}
										else {
											if (DIS[xnei].second != v) {
												mtx_5952[xnei].lock();
												PPR_TYPE::PPR_insert(*PPR, xnei, DIS[xnei].second, v);
												mtx_5952[xnei].unlock();
											}
											if (DIS[xnei].second != xnei) {
												mtx_5952[v].lock();
												PPR_TYPE::PPR_insert(*PPR, v, DIS[xnei].second, xnei);
												mtx_5952[v].unlock();
											}
										}
									}
								}
							}
							for (int i : Dis_changed) {
								DIS[i] = { -1, -1 };
							}

							mtx_595_1.lock();
							Qid_595.push(current_tid);
							mtx_595_1.unlock();

							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void nonHOP_WeightIncreaseMaintenance_improv_batch(graph<int> &instance_graph, two_hop_case_info &mm, std::vector<std::pair<int, int>> &v, std::vector<int> &w_old_vec,
																   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{
					std::vector<affected_label> al1, al3;
					std::vector<pair_label> al2;
					std::map<std::pair<int, int>, weightTYPE> w_old_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_old_map.count(v[i]) == 0)
						{
							w_old_map[v[i]] = w_old_vec[i];
						}
					}

					for (auto &it : w_old_map)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, &al1, &instance_graph, &mm, &w_old_map]
																		  {
							int v1 = it.first.first;
							int v2 = it.first.second;
							weightTYPE w_old = it.second;
							for (auto it : mm.L[v1]) {
								mtx_595[v2].lock();
								long long search_weight = search_sorted_two_hop_label_weight_in_current(mm.L[v2], it.vertex);
								mtx_595[v2].unlock();
								if (it.vertex <= v2 && search_weight == (long long)it.distance + w_old && search_weight < MAX_VALUE) {
									mtx_595_1.lock();
									al1.push_back(affected_label(v2, it.vertex, it.distance + w_old));
									mtx_595_1.unlock();
								}
							}
							for (auto it : mm.L[v2]) {
								mtx_595[v1].lock();
								long long search_weight = search_sorted_two_hop_label_weight_in_current(mm.L[v1], it.vertex);
								mtx_595[v1].unlock();
								if (it.vertex <= v1 && search_weight == (long long)it.distance + w_old && search_weight < MAX_VALUE) {
									mtx_595_1.lock();
									al1.push_back(affected_label(v1, it.vertex, it.distance + w_old));
									mtx_595_1.unlock();
								}
							}
							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
					SPREAD1_batch(instance_graph, &mm.L, al1, &al2, w_old_map, pool_dynamic, results_dynamic, time);
					SPREAD2_batch(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic, time);
					SPREAD3_batch(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic, time);
				}

			}
		}

		namespace algorithm2021
		{
			namespace decrease
			{

				void ProDecreasep_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR,
										std::vector<affected_label> &CL_curr, std::vector<affected_label> *CL_next, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{
					bool is_debug = false;
					if (CL_next->size() > 100000)
					{
						is_debug = true;
					}
					for (const affected_label &it : CL_curr)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, PPR, CL_next, &instance_graph, is_debug]
																		  {
								int v = it.first, u = it.second;

								mtx_595[u].lock();
								std::vector<two_hop_label> Lu = (*L)[u]; // to avoid interlocking
								//std::cout << "lu address is " << &Lu << std::endl;
								//std::cout << "L[u] address is " << &((*(L))[u]) << std::endl;
								mtx_595[u].unlock();

								for (auto nei : instance_graph[v]) {
									int vnei = nei.first;
									long long dnew = it.dis + nei.second;

									if (u < vnei) {
										mtx_595[vnei].lock();
										auto query_result = graph_weighted_two_hop_extract_2hop_label_by_backup_label((*L)[vnei], Lu); // query_result is {distance, common hub}
										mtx_595[vnei].unlock();
										if (query_result.first.t_s == -1 || (query_result.first.distance + query_result.second.distance) > dnew) {
											mtx_595[vnei].lock();
											//if (is_debug) {
											//	std::cout << "2021 decrease function two hop label1 is from " << query_result.first.t_s << " cost " << query_result.first.distance
											//		<< " hop label2 is from " << query_result.second.t_s << " cost " << query_result.second.distance
											//		<< " and all time is " << query_result.first.distance + query_result.second.distance << " greater than " << dnew << std::endl;
											//}

											insert_sorted_two_hop_label((*L)[vnei], u, dnew, time);
											mtx_595[vnei].unlock();
											mtx_595_1.lock();
											CL_next->push_back(affected_label(vnei, u, dnew));
											mtx_595_1.unlock();
										}
										else {
											mtx_595[vnei].lock();
											two_hop_label search_result = search_sorted_two_hop_label_in_current((*L)[vnei], u);
											mtx_595[vnei].unlock();
											if (search_result.distance < 1e7 && search_result.distance > dnew) {
												mtx_595[vnei].lock();
												// std::cout << "decrease label has better answer : old label is " << vnei << " to " << search_result.vertex << " old value is " << search_result.distance << " to " << dnew << " t_s is " << search_result.t_s << std::endl;
												insert_sorted_two_hop_label((*L)[vnei], search_result.vertex, dnew, time);
												// (*L)[vnei][search_result.second].distance = dnew;
												mtx_595[vnei].unlock();
												mtx_595_1.lock();
												CL_next->push_back(affected_label(vnei, u, dnew));
												mtx_595_1.unlock();
											}
											if (query_result.first.vertex != u) {
												mtx_5952[vnei].lock();
												PPR_TYPE::PPR_insert(*PPR, vnei, query_result.first.vertex, u);
												mtx_5952[vnei].unlock();
											}
											if (query_result.first.vertex != vnei) {
												mtx_5952[u].lock();
												PPR_TYPE::PPR_insert(*PPR, u, query_result.first.vertex, vnei);
												mtx_5952[u].unlock();
											}
										}
									}
								}

								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void decrease_maintain(graph<int> &instance_graph, two_hop_case_info &mm, std::vector<std::pair<int, int>> &v, std::vector<int> &w_new,
									   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{
					std::map<std::pair<int, int>, int> w_new_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_new_map.count(v[i]) == 0)
						{
							w_new_map[v[i]] = w_new[i];
						}
						else if (w_new_map[v[i]] > w_new[i])
						{
							w_new_map[v[i]] = w_new[i];
						}
					}

					std::vector<affected_label> CL_curr, CL_next;

					auto &L = mm.L;
					/*
					the following part does not suit parallel computation:
					the reason is that L is changed below, and as a result, in each following loop, L[v2] or L[v1] is locked at each step,
					which means that following loops cannot be actually parallized
					*/

					for (auto &it : w_new_map)
					{
						int v1 = it.first.first, v2 = it.first.second;
						int w_new = it.second;
						for (int sl = 0; sl < 2; sl++)
						{
							if (sl == 1)
							{
								std::swap(v1, v2);
							}
							for (auto it : L[v1])
							{
								int v = it.vertex;
								long long dis = it.distance + w_new;
								if (v <= v2)
								{
									auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current(L, v, v2); // query_result is {distance, common hub}
									if (query_result.first > dis)
									{
										mtx_595[v2].lock();
										insert_sorted_two_hop_label(L[v2], v, dis, time);
										mtx_595[v2].unlock();
										CL_curr.push_back(affected_label(v2, v, dis));
									}
									else
									{
										auto search_result = search_sorted_two_hop_label_weight_and_hub_in_current(L[v2], v);
										if (search_result.first < 1e7 && search_result.first > dis)
										{
											mtx_595[v2].lock();
											// ����ֱ���滻 ʹ�÷���
											insert_sorted_two_hop_label(L[v2], search_result.second, dis, time);
											mtx_595[v2].unlock();
											// L[v2][search_result.second].distance = dis;
											CL_curr.push_back(affected_label(v2, v, dis));
										}
										if (query_result.second != v)
										{
											PPR_TYPE::PPR_insert(mm.PPR, v2, query_result.second, v);
										}
										if (query_result.second != v2)
										{
											PPR_TYPE::PPR_insert(mm.PPR, v, query_result.second, v2);
										}
									}
								}
							}
						}
					}

					while (CL_curr.size())
					{
						// std::cout << "2021 decrease cl_curr size is " << CL_curr.size() << std::endl;
						ProDecreasep_batch(instance_graph, &mm.L, &mm.PPR, CL_curr, &CL_next, pool_dynamic, results_dynamic, time);
						CL_curr = CL_next;
						std::vector<affected_label>().swap(CL_next);
					}
				}
			}
			namespace increase
			{
				void PI11(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L,
						  std::vector<affected_label> &al1_curr, std::vector<affected_label> *al1_next,
						  std::map<std::pair<int, int>, weightTYPE> &w_old_map,
						  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{
					bool is_debug = false;
					if (al1_curr.size() >= 100000)
					{
						is_debug = true;
					}
					for (auto it : al1_curr)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, al1_next, &instance_graph, &w_old_map, time, is_debug]
																		  {
							for (auto nei : instance_graph[it.first]) {
								mtx_595[nei.first].lock();
								two_hop_label search_weight = search_sorted_two_hop_label_in_current((*L)[nei.first], it.second);
								mtx_595[nei.first].unlock();
								weightTYPE w_old;
								int i = 0;
								if (w_old_map.count(std::pair<int, int>(it.first, nei.first)) > 0) {
									i = 1;
									w_old = w_old_map[std::pair<int, int>(it.first, nei.first)];
								}
								else if (w_old_map.count(std::pair<int, int>(nei.first, it.first)) > 0) {
									i = 2;
									w_old = w_old_map[std::pair<int, int>(nei.first, it.first)];
								}
								else {
									i = 3;
									w_old = nei.second;
								}
								if (it.dis + w_old == search_weight.distance && search_weight.t_s!=time) {
									mtx_595_1.lock();
									/*std::cout << i << " and its t_s is " << search_weight.t_s << std::endl;
									std::cout << nei.first << " to " << it.second << " weight is it.dis " << it.dis << " + w_old " << w_old << "=" << it.dis + w_old << " and search_weight is " << search_weight.distance << std::endl;*/
									// if (is_debug) {
									// 	std::cout << nei.first << " to " << it.second << " weight is it.dis " << it.dis << " + w_old " << w_old << "=" << it.dis + w_old << " and search_weight is " << search_weight.distance << " and time is " <<search_weight.t_s<<std::endl;
									// }
									al1_next->push_back(affected_label(nei.first, it.second, search_weight.distance));
									mtx_595_1.unlock();
								}
							}
							mtx_595[it.first].lock();
							insert_sorted_two_hop_label((*L)[it.first], it.second, MAX_VALUE, time); // this does not change the size of L[it->first] here, so does not need to lock here
							mtx_595[it.first].unlock();
							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void PI12(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR,
						  std::vector<affected_label> &al1_curr, std::vector<pair_label> *al2_next, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{
					for (auto it : al1_curr)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al2_next, &instance_graph, time]
																		  {

							int v = it.first, u = it.second;
							mtx_5952[v].lock();
							std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
							mtx_5952[v].unlock();
							temp.push_back(u);

							mtx_595[v].lock_shared();
							auto Lv = (*L)[v]; // to avoid interlocking
							mtx_595[v].unlock_shared();

							for (auto t : temp) {
								if (v < t) {
									long long d1 = MAX_VALUE;
									for (auto nei : instance_graph[t]) {
										mtx_595[nei.first].lock();
										d1 = std::min(d1, search_sorted_two_hop_label_in_current((*L)[nei.first], v).distance + (long long)nei.second);
										mtx_595[nei.first].unlock();
									}
									mtx_595[t].lock();
									auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[t], Lv);
									mtx_595[t].unlock();
									if (query_result.first > d1) {
										mtx_595[t].lock();
										insert_sorted_two_hop_label((*L)[t], v, d1, time);
										mtx_595[t].unlock();
										mtx_595_1.lock();
										al2_next->push_back(pair_label(t, v));
										mtx_595_1.unlock();
									}
									else {
										if (query_result.second != v) {
											mtx_5952[t].lock();
											PPR_TYPE::PPR_insert(*PPR, t, query_result.second, v);
											mtx_5952[t].unlock();
										}
										if (query_result.second != t) {
											mtx_5952[v].lock();
											PPR_TYPE::PPR_insert(*PPR, v, query_result.second, t);
											mtx_5952[v].unlock();
										}
									}
								}
								if (t < v) {
									long long d1 = MAX_VALUE;
									for (auto nei : instance_graph[v]) {
										mtx_595[nei.first].lock();
										d1 = std::min(d1, search_sorted_two_hop_label_in_current((*L)[nei.first], t).distance + (long long)nei.second);
										mtx_595[nei.first].unlock();
									}
									mtx_595[t].lock();
									auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[t], Lv);
									mtx_595[t].unlock();
									if (query_result.first > d1) {
										mtx_595[v].lock();
										insert_sorted_two_hop_label((*L)[v], t, d1, time);
										mtx_595[v].unlock();
										mtx_595_1.lock();
										al2_next->push_back(pair_label(v, t));
										mtx_595_1.unlock();
									}
									else {
										if (query_result.second != v) {
											mtx_5952[t].lock();
											PPR_TYPE::PPR_insert(*PPR, t, query_result.second, v);
											mtx_5952[t].unlock();
										}
										if (query_result.second != t) {
											mtx_5952[v].lock();
											PPR_TYPE::PPR_insert(*PPR, v, query_result.second, t);
											mtx_5952[v].unlock();
										}
									}
								}
							}

							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void PI22(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR,
						  std::vector<pair_label> &al2_curr, std::vector<pair_label> *al2_next, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{

					for (auto it = al2_curr.begin(); it != al2_curr.end(); it++)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al2_next, &instance_graph, time]
																		  {
							try {
								//if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
								//	throw reach_limit_time_string;
								//}

								mtx_595[it->second].lock_shared();
								auto Lxx = (*L)[it->second]; // to avoid interlocking
								mtx_595[it->second].unlock_shared();

								for (auto nei : instance_graph[it->first]) {
									if (nei.first > it->second) {
										mtx_595[it->first].lock();
										long long search_result = search_sorted_two_hop_label_in_current((*L)[it->first], it->second).distance + (long long)nei.second;
										mtx_595[it->first].unlock();
										mtx_595[nei.first].lock();
										auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[nei.first], Lxx);
										mtx_595[nei.first].unlock();
										if (query_result.first > search_result) {
											mtx_595[nei.first].lock();
											insert_sorted_two_hop_label((*L)[nei.first], it->second, search_result, time);
											mtx_595[nei.first].unlock();
											mtx_595_1.lock();
											al2_next->push_back(pair_label(nei.first, it->second));
											mtx_595_1.unlock();
										}
										else {
											if (query_result.second != it->second) {
												mtx_5952[nei.first].lock();
												PPR_TYPE::PPR_insert(*PPR, nei.first, query_result.second, it->second);
												mtx_5952[nei.first].unlock();
											}
											if (query_result.second != nei.first) {
												mtx_5952[it->second].lock();
												PPR_TYPE::PPR_insert(*PPR, it->second, query_result.second, nei.first);
												mtx_5952[it->second].unlock();
											}
										}
									}
								}
							}
							catch (std::string& e) {
								std::cout << "error is " << e << std::endl;
							}
							return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void increase_maintain(graph<int> &instance_graph, two_hop_case_info &mm, std::vector<std::pair<int, int>> &v, std::vector<weightTYPE> &w_old_vec,
									   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{

					std::map<std::pair<int, int>, weightTYPE> w_old_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_old_map.count(v[i]) == 0)
						{
							w_old_map[v[i]] = w_old_vec[i];
						}
					}

					std::vector<affected_label> al1_curr, al1_next;
					std::vector<pair_label> al2_curr, al2_next;

					for (auto &iter : w_old_map)
					{
						int v1 = iter.first.first;
						int v2 = iter.first.second;
						weightTYPE w_old = iter.second;

						for (auto it : mm.L[v1])
						{
							mtx_595[v2].lock();
							long long search_weight = search_sorted_two_hop_label_weight_in_current(mm.L[v2], it.vertex);
							mtx_595[v2].unlock();
							if (it.vertex <= v2 && search_weight >= (long long)it.distance + w_old && search_weight < MAX_VALUE)
							{
								al1_curr.push_back(affected_label(v2, it.vertex, it.distance + w_old));
							}
						}
						for (auto it : mm.L[v2])
						{
							mtx_595[v1].lock();
							long long search_weight = search_sorted_two_hop_label_weight_in_current(mm.L[v1], it.vertex);
							mtx_595[v1].unlock();
							if (it.vertex <= v1 && search_weight >= (long long)it.distance + w_old && search_weight < MAX_VALUE)
							{
								al1_curr.push_back(affected_label(v1, it.vertex, it.distance + w_old));
							}
						}
					}

					while (al1_curr.size() || al2_curr.size())
					{

						// if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
						//	throw reach_limit_time_string;
						// }
						PI11(instance_graph, &mm.L, al1_curr, &al1_next, w_old_map, pool_dynamic, results_dynamic, time);
						PI12(instance_graph, &mm.L, &mm.PPR, al1_curr, &al2_next, pool_dynamic, results_dynamic, time);
						PI22(instance_graph, &mm.L, &mm.PPR, al2_curr, &al2_next, pool_dynamic, results_dynamic, time);
						// std::cout << "increase 2021 al1_cuur size is " << al1_curr.size() << " al2_curr size is " << al2_curr.size() << " al2_next size is " << al2_next.size() << std::endl;
						al1_curr = al1_next;
						al2_curr = al2_next;
						std::vector<affected_label>().swap(al1_next);
						std::vector<pair_label>().swap(al2_next);
					}
				}

			}
		}

	}
	namespace hop
	{
		int TwoM_value = 2 * 1e6;

		std::shared_mutex mtx_599_1, mtx_599_2;

		std::vector<std::shared_mutex> mtx_ruc_decrease(max_N_ID_for_mtx_599);
		std::vector<std::shared_mutex> mtx_ruc_increase(max_N_ID_for_mtx_599);
		std::vector<std::shared_mutex> mtx_2021_decrease(max_N_ID_for_mtx_599);
		std::vector<std::shared_mutex> mtx_2021_increase(max_N_ID_for_mtx_599);

		std::vector<std::shared_mutex> mtx_5992(max_N_ID_for_mtx_599);
		std::queue<int> Qid_599_v2, Qid_599_v3;
		std::vector<std::vector<std::pair<int, int>>> dist_hop_599_v2, dist_hop_599_v3;
		std::vector<std::vector<std::vector<weightTYPE>>> Q_value;
		class hop_constrained_affected_label
		{
		public:
			int first, second, hop;
			long long int dis;
			hop_constrained_affected_label() {}
			hop_constrained_affected_label(int _first, int _second, int _hop, long long int _dis)
			{
				first = _first;
				second = _second;
				hop = _hop;
				dis = _dis;
			}
		};

		class hop_constrained_pair_label
		{
		public:
			int first, second;
			int hop;
			hop_constrained_pair_label(int _first, int _second, int _hop)
			{
				first = _first;
				second = _second;
				hop = _hop;
			}
			bool operator==(const hop_constrained_pair_label other) const
			{
				return (first == other.first && second == other.second && hop == other.hop);
			}
			bool operator<(const hop_constrained_pair_label other) const
			{ // used to sort/search pair_label2 in set
				if (first != other.first)
					return first < other.first;
				if (second != other.second)
					return second < other.second;
				return hop < other.hop;
			}
		};

		class hop_constrained_label_v2
		{
		public:
			int hub_vertex, hop;
			weightTYPE distance;
			hop_constrained_label_v2(int _vertex, int _hop, weightTYPE _dis)
			{
				hub_vertex = _vertex;
				hop = _hop;
				distance = _dis;
			}
			bool operator==(const hop_constrained_label_v2 other) const
			{
				return (hub_vertex == other.hub_vertex && hop == other.hop && distance == other.distance);
			}
			bool operator<(const hop_constrained_label_v2 other) const
			{ // used to sort/search pair_label2 in set
				if (hub_vertex != other.hub_vertex)
					return hub_vertex < other.hub_vertex;
				if (hop != other.hop)
					return hop < other.hop;
				return distance < other.distance;
			}
		};

		struct hop_constrained_node_for_DIFFUSE
		{
			int index;
			int hop;
			weightTYPE disx;
			hop_constrained_node_for_DIFFUSE() {}
			hop_constrained_node_for_DIFFUSE(int _u, int _hop, weightTYPE _dis)
			{
				index = _u;
				hop = _hop;
				disx = _dis;
			}
		}; // define the node in the queue

		typedef typename boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE>::handle_type hop_constrained_handle_t_for_DIFFUSE;

		bool operator<(hop_constrained_node_for_DIFFUSE const &x, hop_constrained_node_for_DIFFUSE const &y)
		{
			return x.disx > y.disx; // < is the max-heap; > is the min heap
		}

		void initialize_global_values_dynamic_hop_constrained(int N, int thread_num, int upper_k)
		{

			dist_hop_599_v2.resize(thread_num);
			dist_hop_599_v3.resize(thread_num);
			Q_value.resize(thread_num);
			std::queue<int>().swap(Qid_599_v2);
			std::queue<int>().swap(Qid_599_v3);
			for (int i = 0; i < thread_num; i++)
			{

				Qid_599_v2.push(i);
				Qid_599_v3.push(i);
				dist_hop_599_v2[i].resize(N, {-1, 0});
				dist_hop_599_v3[i].resize(N, {-1, 0});
				Q_value[i].resize(N, std::vector<weightTYPE>(upper_k + 1, MAX_VALUE));
			}
		}
		namespace ruc
		{
			namespace decrease
			{
				void decrease_maintain_step1_batch(std::map<std::pair<int, int>, weightTYPE> &v_map, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR, std::vector<hop_constrained_affected_label> *CL,
												   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
				{
					for (auto v_map_item : v_map)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([t, v_map_item, L, PPR, CL]
																		  {

								int v1 = v_map_item.first.first, v2 = v_map_item.first.second;
								weightTYPE w_new = v_map_item.second;
								for (int sl = 0; sl < 2; sl++)
								{
									if (sl == 1)
									{
										std::swap(v1, v2);
									}
									for (auto& it : (*L)[v1])
									{
										if (it.hub_vertex <= v2 && (long long int)it.distance + w_new < TwoM_value && it.t_e == std::numeric_limits<int>::max())
										{
											auto query_result = hop_constrained_extract_distance_and_hub(*L, it.hub_vertex, v2, it.hop + 1); // query_result is {distance, common hub}
											if ((long long int)query_result.first > (long long int)it.distance + w_new)
											{
												mtx_599_1.lock();
												CL->push_back(hop_constrained_affected_label{ v2, it.hub_vertex, it.hop + 1, it.distance + w_new });
												mtx_599_1.unlock();
											}
											else
											{
												auto search_result = search_sorted_hop_constrained_weight_two_hop_label((*L)[v2], it.hub_vertex, it.hop + 1);
												if (search_result < MAX_VALUE && search_result > it.distance + w_new)
												{
													mtx_599_1.lock();
													CL->push_back(hop_constrained_affected_label{ v2, it.hub_vertex, it.hop + 1, it.distance + w_new });
													mtx_599_1.unlock();
												}
												if (query_result.second != -1 && query_result.second != it.hub_vertex)
												{
													mtx_5992[v2].lock();
													PPR_TYPE::PPR_insert(*PPR, v2, query_result.second, it.hub_vertex);
													mtx_5992[v2].unlock();
												}
												if (query_result.second != -1 && query_result.second != v2)
												{
													mtx_5992[it.hub_vertex].lock();
													PPR_TYPE::PPR_insert(*PPR, it.hub_vertex, query_result.second, v2);
													mtx_5992[it.hub_vertex].unlock();
												}
											}
										}
									}
								}
								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void DIFFUSE_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR, std::vector<hop_constrained_affected_label> &CL,
								   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k, int t)
				{
					std::map<hop_constrained_pair_label, weightTYPE> CL_edge_map;
					for (auto &it : CL)
					{
						if (CL_edge_map.count({it.first, it.second, it.hop}) == 0)
						{
							CL_edge_map[{it.first, it.second, it.hop}] = it.dis;
						}
						else if (CL_edge_map[{it.first, it.second, it.hop}] > it.dis)
						{
							CL_edge_map[{it.first, it.second, it.hop}] = it.dis;
						}
					}

					// extract each unique hub v and its (u,hop,dis) list
					std::map<int, std::vector<hop_constrained_label_v2>> CL_map;
					for (auto &it : CL_edge_map)
					{
						int u = it.first.first;
						int v = it.first.second;
						int hop = it.first.hop;
						weightTYPE dis = it.second;
						if (CL_map.count(v) == 0)
						{
							std::vector<hop_constrained_label_v2> vec_with_hub_v;
							hop_constrained_label_v2 tmp(u, hop, dis);
							vec_with_hub_v.emplace_back(tmp);
							CL_map[v] = vec_with_hub_v;
						}
						else
						{
							std::vector<hop_constrained_label_v2> vec_with_hub_v = CL_map[v];
							hop_constrained_label_v2 tmp(u, hop, dis);
							vec_with_hub_v.emplace_back(tmp);
							CL_map[v] = vec_with_hub_v;
						}
					}

					for (auto &it : CL_map)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, &instance_graph, PPR, upper_k]
																		  {
								mtx_599_1.lock();
								int current_tid = Qid_599_v2.front();
								Qid_599_v2.pop();
								mtx_599_1.unlock();

								int v = it.first;
								std::vector<hop_constrained_label_v2> vec_with_hub_v = it.second;

								mtx_ruc_decrease[v].lock_shared();
								auto Lv = (*L)[v]; // to avoid interlocking
								mtx_ruc_decrease[v].unlock_shared();

								std::vector<int> dist_hop_changes;
								auto& dist_hop = dist_hop_599_v2[current_tid];
								boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE> pq;
								std::map<std::pair<int, int>, std::pair<hop_constrained_handle_t_for_DIFFUSE, int>> Q_handle;
								std::vector<int> hubs;
								hubs.resize(instance_graph.size(), -1);
								auto& Q_VALUE = Q_value[current_tid];

								for (auto& it : vec_with_hub_v) {
									int u = it.hub_vertex;
									int h_v = it.hop;
									weightTYPE du = it.distance;

									dist_hop[u] = { du, h_v }; //  {dis, hop}
									dist_hop_changes.push_back(u);
									hop_constrained_node_for_DIFFUSE tmp;
									tmp.index = u;
									tmp.hop = h_v;
									tmp.disx = du;
									Q_handle[{u, h_v}] = { pq.push({tmp}), du }; //{node_for_DIFFUSE_v2,dis}
									// what the meaning of Q_VALUE? mark the data of pq
									if (h_v <= upper_k)
										Q_VALUE[u][h_v] = du;

								}

								while (!pq.empty())
								{
									int x = pq.top().index;
									int xhv = pq.top().hop;
									weightTYPE dx = pq.top().disx;
									pq.pop();
									if (xhv <= upper_k)
										Q_VALUE[x][xhv] = MAX_VALUE;

									mtx_ruc_decrease[x].lock_shared();
									weightTYPE d_old = search_sorted_hop_constrained_weight_two_hop_label((*L)[x], v, xhv);
									mtx_ruc_decrease[x].unlock_shared();
									if (dx >= 0 && dx < d_old)
									{
										mtx_ruc_decrease[x].lock();
										insert_sorted_hop_constrained_two_hop_label((*L)[x], v, xhv, dx, t);
										mtx_ruc_decrease[x].unlock();
									}

									if (xhv + 1 > upper_k)
										continue;

									for (auto nei : instance_graph[x])
									{
										if (dx + nei.second >= TwoM_value)
											continue;
										int xnei = nei.first;
										int hop_nei = xhv + 1;
										long long int d_new = dx + (long long int)nei.second;
										hop_constrained_node_for_DIFFUSE node = { xnei, hop_nei, (weightTYPE)d_new };
										if (v < xnei)
										{
											if (dist_hop[xnei].first == -1)
											{
												// Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
												// Q_VALUE[xnei][hop_nei] = d_new;
												mtx_ruc_decrease[xnei].lock_shared();
												std::pair<int, int> temp_dis = graph_weighted_two_hop_extract_distance_and_hop_by_backup_label((*L)[xnei], Lv, xhv + 1);
												std::pair<int, int> temp_dis_hub = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[xnei], Lv, xhv + 1);
												//std::pair<int, int> temp_dis = hop_constrained_extract_distance_and_hop(*L, xnei, v, xhv + 1);
												mtx_ruc_decrease[xnei].unlock_shared();
												hubs[xnei] = temp_dis_hub.second;

												dist_hop[xnei].first = temp_dis.first;
												dist_hop[xnei].second = temp_dis.second;
												dist_hop_changes.push_back(xnei);
											}

											if (d_new < dist_hop[xnei].first)
											{

												//if (Q_handle.find({xnei, hop_nei}) != Q_handle.end())
												if (Q_VALUE[xnei][hop_nei] < MAX_VALUE)
												{
													if (Q_handle[{xnei, hop_nei}].second > d_new)
													{
														pq.update(Q_handle[{xnei, hop_nei}].first, node);
														Q_handle[{xnei, hop_nei}].second = d_new;
													}
												}
												else
												{
													Q_handle[{xnei, hop_nei}] = { pq.push(node), d_new };
												}
												dist_hop[xnei].first = d_new;
												dist_hop[xnei].second = hop_nei;
												dist_hop_changes.push_back(xnei);
												hubs[xnei] = v;
												Q_VALUE[xnei][hop_nei] = d_new;
											}
											else if (hop_nei < dist_hop[xnei].second)
											{

												//if (Q_handle.find({xnei, hop_nei}) != Q_handle.end())
												if (Q_VALUE[xnei][hop_nei] < MAX_VALUE)
												{
													if (Q_handle[{xnei, hop_nei}].second > d_new)
													{
														pq.update(Q_handle[{xnei, hop_nei}].first, node);
														Q_handle[{xnei, hop_nei}].second = d_new;
													}
												}
												else
												{
													Q_handle[{xnei, hop_nei}] = { pq.push(node), d_new };
												}
												Q_VALUE[xnei][hop_nei] = d_new;
											}


											if (dist_hop[xnei].first < d_new)
											{
												if (hubs[xnei] != -1 && hubs[xnei] != v)
												{
													mtx_5992[xnei].lock();
													PPR_TYPE::PPR_insert(*PPR, xnei, hubs[xnei], v);
													mtx_5992[xnei].unlock();
												}
												if (hubs[xnei] != -1 && hubs[xnei] != xnei)
												{
													mtx_5992[v].lock();
													PPR_TYPE::PPR_insert(*PPR, v, hubs[xnei], xnei);
													mtx_5992[v].unlock();
												}
											}


										}
									}
								}


								for (int i : dist_hop_changes)
								{
									dist_hop[i] = { -1, 0 };
								}

								mtx_599_1.lock();
								Qid_599_v2.push(current_tid);
								mtx_599_1.unlock();
								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void HOP_WeightDecreaseMaintenance_improv_batch(graph<int> &instance_graph, two_hop_case_info &mm,
																std::vector<std::pair<int, int>> &v, std::vector<int> &w_new, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
				{
					std::map<std::pair<int, int>, weightTYPE> w_new_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_new_map.count(v[i]) == 0)
						{
							w_new_map[v[i]] = w_new[i];
						}
						else if (w_new_map[v[i]] > w_new[i])
						{
							w_new_map[v[i]] = w_new[i];
						}
					}
					std::vector<hop_constrained_affected_label> CL;
					decrease_maintain_step1_batch(w_new_map, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic, t);
					std::cout << "ruc decrease CL size is" << CL.size() << std::endl;
					DIFFUSE_batch(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic, mm.upper_k, t);
				}
			}
			namespace increase
			{
				void HOP_maintain_SPREAD1_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L,
												std::vector<hop_constrained_affected_label> &al1, std::vector<hop_constrained_pair_label> *al2, std::map<std::pair<int, int>, int> &w_old_map, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
				{

					for (auto it : al1)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, al2, &instance_graph, &w_old_map]
																		  {
								std::queue<hop_constrained_node_for_DIFFUSE> q; //(u,h_v, d)
								int v = it.second;
								q.push(hop_constrained_node_for_DIFFUSE(it.first, it.hop, it.dis));
								while (!q.empty()) {
									int x = q.front().index;
									int h_x = q.front().hop;
									int dx = q.front().disx;
									q.pop();
									mtx_ruc_increase[x].lock();
									insert_sorted_hop_constrained_two_hop_label((*L)[x], v, h_x, MAX_VALUE, t); // this does not change the size of L[x] here, so does not need to lock here
									mtx_ruc_increase[x].unlock();
									mtx_599_1.lock();
									al2->push_back(hop_constrained_pair_label{x, v, h_x});
									mtx_599_1.unlock();

									for (auto nei : instance_graph[x])
									{
										if (v < nei.first) {
											mtx_ruc_increase[nei.first].lock_shared();
											int search_weight = search_sorted_hop_constrained_weight_two_hop_label((*L)[nei.first], v, h_x + 1);
											mtx_ruc_increase[nei.first].unlock_shared();
											int w_old = nei.second;
											if (w_old_map.count(std::pair<int, int>(x, nei.first)) > 0) {
												w_old = w_old_map[std::pair<int, int>(x, nei.first)];
											}
											else if (w_old_map.count(std::pair<int, int>(nei.first, x)) > 0) {
												w_old = w_old_map[std::pair<int, int>(nei.first, x)];
											}
											else {
												w_old = nei.second;
											}
											if (dx + w_old <= search_weight && search_weight < MAX_VALUE) {
												// mtx_599_2.lock();
												// std::cout << "x is "<< x <<" v is " << v << " x's nei is " << nei.first <<" dx is " << dx << " w_lod is " << w_old <<" search_weight is " << search_weight << std::endl;
												// mtx_599_2.unlock();;
												q.push(hop_constrained_node_for_DIFFUSE(nei.first, h_x + 1, dx + nei.second));
											}

										}
									}
								}

								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void HOP_maintain_SPREAD2_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR,
												std::vector<hop_constrained_pair_label> &al2, std::vector<hop_constrained_affected_label> *al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k)
				{
					for (const auto &it : al2)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([&it, L, PPR, al3, &instance_graph, upper_k]
																		  {
								int v = it.first, u = it.second, h_u = it.hop;
								mtx_5992[v].lock_shared();
								std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
								mtx_5992[v].unlock_shared();
								temp.push_back(u);
								mtx_ruc_increase[v].lock_shared();
								auto Lv = (*L)[v]; // to avoid interlocking
								mtx_ruc_increase[v].unlock_shared();

								for (auto t : temp) {
									if (v < t) {
										long long int d1 = MAX_VALUE;
										int hop_vn = 0;
										for (const auto& nei : instance_graph[t]) {
											//mtx_599[nei.first].lock();
											std::pair<int, int> dis_hop = get_shortest_distance_hop_two_hop_label2((*L)[nei.first], v);
											//mtx_599[nei.first].unlock();
											if (d1 > dis_hop.first + (long long int)nei.second)
											{
												d1 = dis_hop.first + (long long int)nei.second;
												hop_vn = dis_hop.second;
											}

										}
										
										// if(d1 >= TwoM_value)
										// 	continue;
										
										for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++)
										{
											if (hop_i > upper_k)
												break;
											long long int di = MAX_VALUE;
											for (const auto& nei : instance_graph[t]) {
												//mtx_599[nei.first].lock();
												di = std::min(di, search_sorted_hop_constrained_weight_two_hop_label((*L)[nei.first], v, hop_i - 1) + (long long int)nei.second);
												//mtx_599[nei.first].unlock();
											}
											if (di >= TwoM_value)
												continue;
											//mtx_599[t].lock_shared();
											//auto query_result = graph_hash_of_mixed_weightejd_two_hop_v2_extract_distance_no_reduc2(*L, t.first, v, hop_i);
											auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[t], Lv, hop_i);
											//mtx_599[t].unlock_shared();

											if (query_result.first > di) { // only add new label when it's absolutely necessary
												mtx_599_1.lock();
												//cout<<"query_result.first > d1 + 1e-5: "<<t_first<<' '<<v << ' ' << hop_vn+1 << ' ' << d1 << endl;
												al3->push_back(hop_constrained_affected_label{t, v, hop_i, di});
												mtx_599_1.unlock();

											}
											else {
												if (query_result.second != -1 && query_result.second != v) {
													mtx_5992[t].lock();
													PPR_TYPE::PPR_insert(*PPR, t, query_result.second, v);
													mtx_5992[t].unlock();
												}
												if (query_result.second != -1 && query_result.second != t) {
													mtx_5992[v].lock();
													PPR_TYPE::PPR_insert(*PPR, v, query_result.second, t);
													mtx_5992[v].unlock();
												}
											}
										}
									}
									if (t < v) {
										long long int d1 = MAX_VALUE;
										int hop_vn = 0;
										for (const auto& nei : instance_graph[v]) {
											//d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], t_first, t.second) + (int)nei.second);
											//mtx_599[nei.first].lock();
											std::pair<int, int> dis_hop = get_shortest_distance_hop_two_hop_label2((*L)[nei.first], t);
											//mtx_599[nei.first].unlock();
											if (d1 > dis_hop.first + (long long int)nei.second)
											{
												d1 = dis_hop.first + (long long int)nei.second;
												hop_vn = dis_hop.second;
											}

										}

										// if(d1 >= TwoM_value)
										// 	continue;

										for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++) {
											if (hop_i > upper_k)
												break;
											long long int di = MAX_VALUE;
											for (auto nei : instance_graph[v]) {
												//mtx_599[nei.first].lock();
												di = std::min(di, search_sorted_hop_constrained_weight_two_hop_label((*L)[nei.first], t, hop_i - 1) + (long long int)nei.second);
												//mtx_599[nei.first].unlock();
											}

											if (di >= TwoM_value)
												continue;

											//mtx_599[t].lock_shared();
											//auto query_result = graph_hash_of_mixed_weighted_two_hop_v2_extract_distance_no_reduc2(*L, v, t_first, hop_i);
											auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[t], Lv, hop_i);
											//mtx_599[t].unlock_shared();

											if (query_result.first > di) {
												mtx_599_1.lock();
												al3->push_back(hop_constrained_affected_label{v, t, hop_i, di});
												mtx_599_1.unlock();
											}
											else {
												if (query_result.second != -1 && query_result.second != v) {
													mtx_5992[t].lock();
													PPR_TYPE::PPR_insert(*PPR, t, query_result.second, v);
													mtx_5992[t].unlock();
												}
												if (query_result.second != -1 && query_result.second != t) {
													mtx_5992[v].lock();
													PPR_TYPE::PPR_insert(*PPR, v, query_result.second, t);
													mtx_5992[v].unlock();
												}
											}
										}
									}
								}

								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void HOP_maintain_SPREAD3_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR, std::vector<hop_constrained_affected_label> &al3,
												ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k, int t)
				{
					std::map<hop_constrained_pair_label, int> al3_edge_map;
					for (auto &it : al3)
					{
						if (al3_edge_map.count({it.first, it.second, it.hop}) == 0)
						{
							al3_edge_map[{it.first, it.second, it.hop}] = it.dis;
						}
						else if (al3_edge_map[{it.first, it.second, it.hop}] > it.dis)
						{
							al3_edge_map[{it.first, it.second, it.hop}] = it.dis;
						}
					}

					// extract each unique hub v and its (u,dis) list
					std::map<int, std::vector<hop_constrained_label_v2>> al3_map; // al3_map[v]=(u1,hop1,dis1),(u2,hop2,dis2)...
					for (auto &it : al3_edge_map)
					{
						int u = it.first.first;
						int v = it.first.second;
						int hop = it.first.hop;
						int dis = it.second;
						if (al3_map.count(v) == 0)
						{
							std::vector<hop_constrained_label_v2> vec_with_hub_v;
							hop_constrained_label_v2 tmp(u, hop, dis);
							vec_with_hub_v.emplace_back(tmp);
							al3_map[v] = vec_with_hub_v;
						}
						else
						{
							std::vector<hop_constrained_label_v2> vec_with_hub_v = al3_map[v];
							hop_constrained_label_v2 tmp(u, hop, dis);
							vec_with_hub_v.emplace_back(tmp);
							al3_map[v] = vec_with_hub_v;
						}
					}

					for (auto &it : al3_map)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, &instance_graph, PPR, upper_k]
																		  {
								mtx_599_1.lock();
								int current_tid = Qid_599_v2.front();
								Qid_599_v2.pop();
								mtx_599_1.unlock();

								int v = it.first;
								std::vector<hop_constrained_label_v2> vec_with_hub_v = it.second;

								mtx_ruc_increase[v].lock_shared();
								auto Lv = (*L)[v]; // to avoid interlocking
								mtx_ruc_increase[v].unlock_shared();

								std::vector<int> dist_hop_changes;
								auto& dist_hop = dist_hop_599_v2[current_tid];
								boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE> pq;
								std::map<std::pair<int, int>, std::pair<hop_constrained_handle_t_for_DIFFUSE, int>> Q_handle;
								std::vector<int> hubs;
								hubs.resize(instance_graph.size(), -1);
								auto& Q_VALUE = Q_value[current_tid];

								for (auto& it : vec_with_hub_v)
								{
									int u = it.hub_vertex;
									int h_v = it.hop;
									int du = it.distance;

									mtx_ruc_increase[u].lock_shared();
									auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[u], Lv, h_v);
									mtx_ruc_increase[u].unlock_shared();

									bool flag = false;
									if (query_result.first < du)
									{
										if (query_result.second != -1 && query_result.second != v)
										{
											mtx_5992[u].lock();
											PPR_TYPE::PPR_insert(*PPR, u, query_result.second, v);
											mtx_5992[u].unlock();
											flag = true;
										}
										if (query_result.second != -1 && query_result.second != u)
										{
											mtx_5992[v].lock();
											PPR_TYPE::PPR_insert(*PPR, v, query_result.second, u);
											mtx_5992[v].unlock();
											flag = true;
										}
									}

									if (flag == true)
									{
										continue;
									}

									dist_hop[u] = { du, h_v }; //  {dis, hop}
									dist_hop_changes.push_back(u);
									hop_constrained_node_for_DIFFUSE tmp;
									tmp.index = u;
									tmp.hop = h_v;
									tmp.disx = du;
									Q_handle[{u, h_v}] = { pq.push({tmp}), du }; //{hop_constrained_node_for_DIFFUSE,dis}
									Q_VALUE[u][h_v] = du;
								}

								while (!pq.empty())
								{
									int x = pq.top().index;
									int xhv = pq.top().hop;
									int dx = pq.top().disx;
									pq.pop();
									if (xhv <= upper_k)
										Q_VALUE[x][xhv] = MAX_VALUE;

									mtx_ruc_increase[x].lock_shared();
									int d_old = search_sorted_hop_constrained_weight_two_hop_label((*L)[x], v, xhv);
									mtx_ruc_increase[x].unlock_shared();
									if (dx >= 0 && dx < d_old)
									{
										mtx_ruc_increase[x].lock();
										insert_sorted_hop_constrained_two_hop_label((*L)[x], v, xhv, dx, t);
										mtx_ruc_increase[x].unlock();
									}

									if (xhv + 1 > upper_k)
										continue;

									for (auto nei : instance_graph[x])
									{

										if (dx + nei.second >= TwoM_value)
											continue;

										int xnei = nei.first;
										int hop_nei = xhv + 1;
										long long int d_new = dx + (long long int)nei.second;
										hop_constrained_node_for_DIFFUSE node = { xnei, xhv + 1, (weightTYPE)d_new };

										if (v < xnei)
										{

											if (dist_hop[xnei].first == -1)
											{
												//Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
												//Q_VALUE[xnei][hop_nei] = d_new;
												dist_hop[xnei].first = d_new;
												dist_hop[xnei].second = hop_nei;
												dist_hop_changes.push_back(xnei);

												mtx_ruc_increase[xnei].lock_shared();
												std::pair<int, int> tmp = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[xnei], Lv, xhv + 1);
												mtx_ruc_increase[xnei].unlock_shared();
												hubs[xnei] = tmp.second;
											}
											if (d_new < dist_hop[xnei].first)
											{
												//if (Q_handle.find({xnei, node.hop}) != Q_handle.end())
												if (Q_VALUE[xnei][hop_nei] < MAX_VALUE)
												{
													if (Q_handle[{xnei, hop_nei}].second > d_new)
													{
														pq.update(Q_handle[{xnei, hop_nei}].first, node);
														Q_handle[{xnei, hop_nei}].second = d_new;
													}
												}
												else
												{
													Q_handle[{xnei, hop_nei}] = { pq.push(node), d_new };
												}
												dist_hop[xnei].first = d_new;
												dist_hop[xnei].second = hop_nei;
												dist_hop_changes.push_back(xnei);
												hubs[xnei] = v;
												Q_VALUE[xnei][hop_nei] = d_new;
											}
											else if (hop_nei < dist_hop[xnei].second)
											{
												//if (Q_handle.find({xnei, node.hop}) != Q_handle.end())
												if (Q_VALUE[xnei][hop_nei] < MAX_VALUE)
												{
													if (Q_handle[{xnei, hop_nei}].second > d_new)
													{
														pq.update(Q_handle[{xnei, hop_nei}].first, node);
														Q_handle[{xnei, hop_nei}].second = d_new;
													}
												}
												else
												{
													Q_handle[{xnei, hop_nei}] = { pq.push(node), d_new };
												}
												Q_VALUE[xnei][hop_nei] = d_new;
											}

											if (dist_hop[xnei].first < d_new)
											{

												if (hubs[xnei] != -1 && hubs[xnei] != v)
												{
													mtx_5992[xnei].lock();
													PPR_TYPE::PPR_insert(*PPR, xnei, hubs[xnei], v);
													mtx_5992[xnei].unlock();
												}
												if (hubs[xnei] != -1 && hubs[xnei] != xnei)
												{
													mtx_5992[v].lock();
													PPR_TYPE::PPR_insert(*PPR, v, hubs[xnei], xnei);
													mtx_5992[v].unlock();
												}
											}
										}
									}
								}

								for (int i : dist_hop_changes) {
									dist_hop[i] = { -1, 0 };
								}

								mtx_599_1.lock();
								Qid_599_v2.push(current_tid);
								mtx_599_1.unlock();

								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void HOP_WeightIncreaseMaintenance_improv_batch(graph<int> &instance_graph, two_hop_case_info &mm, std::vector<std::pair<int, int>> &v, std::vector<int> &w_old_vec,
																ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
				{

					std::vector<hop_constrained_affected_label> al1, al3;
					std::vector<hop_constrained_pair_label> al2;

					std::map<std::pair<int, int>, int> w_old_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first < v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_old_map.count(v[i]) == 0)
						{
							w_old_map[v[i]] = w_old_vec[i];
						}
					}

					for (auto iter : w_old_map)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([iter, &al1, &instance_graph, &mm, &w_old_map]
																		  {
								int v1 = iter.first.first;
								int v2 = iter.first.second;
								int w_old = iter.second;
								for (const auto& it : mm.L[v1]) {
									int search_weight = search_sorted_hop_constrained_weight_two_hop_label(mm.L[v2], it.hub_vertex, it.hop + 1);
									if (it.hub_vertex <= v2 && search_weight >= (long long int)it.distance + w_old && search_weight < MAX_VALUE && it.t_e == std::numeric_limits<int>::max()) {
										mtx_599_1.lock();
										al1.push_back(hop_constrained_affected_label{v2, it.hub_vertex, it.hop + 1, it.distance + w_old});
										mtx_599_1.unlock();
									}
								}
								for (const auto& it : mm.L[v2]) {
									int search_weight = search_sorted_hop_constrained_weight_two_hop_label(mm.L[v1], it.hub_vertex, it.hop + 1);
									if (it.hub_vertex <= v1 && search_weight >= (long long int)it.distance + w_old && search_weight < MAX_VALUE && it.t_e == std::numeric_limits<int>::max()) {
										mtx_599_1.lock();
										al1.push_back(hop_constrained_affected_label{v1, it.hub_vertex, it.hop + 1, it.distance + w_old});
										mtx_599_1.unlock();
									}
								}
								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
					HOP_maintain_SPREAD1_batch(instance_graph, &mm.L, al1, &al2, w_old_map, pool_dynamic, results_dynamic, t);
					std::cout << "ruc increase al1 size is " << al1.size() << " and al2 size is " << al2.size() << std::endl;
					HOP_maintain_SPREAD2_batch(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic, mm.upper_k);
					std::cout << "ruc increase al2 size is " << al2.size() << " and al3 size is " << al3.size() << std::endl;
					HOP_maintain_SPREAD3_batch(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic, mm.upper_k, t);
				}
			}
		}
		namespace algorithm2021
		{
			namespace decrease
			{
				void ProDecreasep_batch(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR,
										std::vector<hop_constrained_affected_label> &CL_curr, std::vector<hop_constrained_affected_label> *CL_next,
										ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k, int t)
				{

					for (auto it : CL_curr)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, PPR, CL_next, &instance_graph, upper_k]
																		  {

								int v = it.first, u = it.second;

								mtx_2021_decrease[u].lock_shared();
								auto Lu = (*L)[u]; // to avoid interlocking
								mtx_2021_decrease[u].unlock_shared();

								if (it.hop + 1 > upper_k)
									return 1;

								for (auto nei : instance_graph[v]) {
									int vnei = nei.first;
									int hop_u = it.hop;
									long long dnew = it.dis + nei.second;
									if (u < vnei) {
										mtx_2021_decrease[vnei].lock_shared();
										auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[vnei], Lu, hop_u + 1); // query_result is {distance, common hub}
										mtx_2021_decrease[vnei].unlock_shared();
										if ((long long)query_result.first > dnew) {
											mtx_2021_decrease[vnei].lock();
											insert_sorted_hop_constrained_two_hop_label((*L)[vnei], u, hop_u + 1, dnew, t);
											mtx_2021_decrease[vnei].unlock();
											mtx_599_1.lock();
											CL_next->push_back(hop_constrained_affected_label{vnei, u, hop_u + 1, dnew});
											mtx_599_1.unlock();
										}
										else {
											mtx_2021_decrease[vnei].lock_shared();
											auto search_result = search_sorted_hop_constrained_weight_and_index_two_hop_label((*L)[vnei], u, hop_u + 1);
											mtx_2021_decrease[vnei].unlock_shared();
											if (search_result.first < MAX_VALUE && search_result.first > dnew) {
												mtx_2021_decrease[vnei].lock();
												insert_sorted_hop_constrained_two_hop_label((*L)[vnei], u, hop_u + 1, dnew, t);
												// (*L)[vnei][search_result.second].distance = dnew;
												mtx_2021_decrease[vnei].unlock();
												mtx_599_1.lock();
												CL_next->push_back(hop_constrained_affected_label{vnei, u, hop_u + 1, dnew});
												mtx_599_1.unlock();
											}
											if (query_result.second != u) {
												mtx_5992[vnei].lock();
												PPR_TYPE::PPR_insert(*PPR, vnei, query_result.second, u);
												mtx_5992[vnei].unlock();
											}
											if (query_result.second != vnei) {
												mtx_5992[u].lock();
												PPR_TYPE::PPR_insert(*PPR, u, query_result.second, vnei);
												mtx_5992[u].unlock();
											}
										}
									}
								}

								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void HOP_WeightDecrease2021_batch(graph<int> &instance_graph, two_hop_case_info &mm, std::vector<std::pair<int, int>> &v, std::vector<weightTYPE> &w_new,
												  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
				{
					std::map<std::pair<int, int>, weightTYPE> w_new_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_new_map.count(v[i]) == 0)
						{
							w_new_map[v[i]] = w_new[i];
						}
						else if (w_new_map[v[i]] > w_new[i])
						{
							w_new_map[v[i]] = w_new[i];
						}
					}

					std::vector<hop_constrained_affected_label> CL_curr, CL_next;

					auto &L = mm.L;
					/*
					the following part does not suit parallel computation:
					the reason is that L is changed below, and as a result, in each following loop, L[v2] or L[v1] is locked at each step,
					which means that following loops cannot be actually parallized
					*/
					for (auto &it : w_new_map)
					{
						int v1 = it.first.first, v2 = it.first.second;
						weightTYPE w_new = it.second;
						for (int sl = 0; sl < 2; sl++)
						{
							if (sl == 1)
							{
								std::swap(v1, v2);
							}
							for (auto it : L[v1])
							{
								int v = it.hub_vertex;
								int hop_v = it.hop;
								long long dis = it.distance + w_new;
								if (v <= v2 && it.t_e == std::numeric_limits<int>::max())
								{
									auto query_result = hop_constrained_extract_distance_and_hub(L, v, v2, hop_v + 1); // query_result is {distance, common hub}

									if ((long long)query_result.first > dis)
									{
										insert_sorted_hop_constrained_two_hop_label(L[v2], v, hop_v + 1, dis, t);
										CL_curr.push_back(hop_constrained_affected_label(v2, v, hop_v + 1, dis));
									}
									else
									{
										auto search_result = search_sorted_hop_constrained_weight_and_index_two_hop_label(L[v2], v, hop_v + 1);
										if (search_result.first < MAX_VALUE && search_result.first > dis)
										{
											L[v2][search_result.second].distance = dis;
											CL_curr.push_back(hop_constrained_affected_label(v2, v, hop_v + 1, dis));
										}
										if (query_result.second != v)
										{
											PPR_TYPE::PPR_insert(mm.PPR, v2, query_result.second, v);
										}
										if (query_result.second != v2)
										{
											PPR_TYPE::PPR_insert(mm.PPR, v, query_result.second, v2);
										}
									}
								}
							}
						}
					}
					while (CL_curr.size())
					{
						ProDecreasep_batch(instance_graph, &mm.L, &mm.PPR, CL_curr, &CL_next, pool_dynamic, results_dynamic, mm.upper_k, t);
						CL_curr = CL_next;
						std::vector<hop_constrained_affected_label>().swap(CL_next);
					}
				}

			}
			namespace increase
			{
				void PI11(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L,
						  std::vector<hop_constrained_affected_label> &al1_curr, std::vector<hop_constrained_affected_label> *al1_next,
						  std::map<std::pair<int, int>, weightTYPE> &w_old_map,
						  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
				{

					for (auto it : al1_curr)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, al1_next, &instance_graph, &w_old_map]
																		  {
								for (auto nei : instance_graph[it.first]) {
									mtx_2021_increase[nei.first].lock_shared();
									two_hop_label search_weight = search_sorted_hop_constrained_label_two_hop_label((*L)[nei.first], it.second, it.hop + 1);
									mtx_2021_increase[nei.first].unlock_shared();
									weightTYPE w_old = nei.second;
									if (w_old_map.count(std::pair<int, int>(it.first, nei.first)) > 0) {
										w_old = w_old_map[std::pair<int, int>(it.first, nei.first)];
									}
									else if (w_old_map.count(std::pair<int, int>(nei.first, it.first)) > 0) {
										w_old = w_old_map[std::pair<int, int>(nei.first, it.first)];
									}
									else {
										w_old = nei.second;
									}

									if (it.dis + w_old <= search_weight.distance && search_weight.distance < MAX_VALUE && search_weight.t_s != t) {
										mtx_599_1.lock();
										al1_next->push_back(hop_constrained_affected_label{nei.first, it.second, it.hop + 1, it.dis + w_old});
										mtx_599_1.unlock();
									}
								}
								mtx_2021_increase[it.first].lock();
								insert_sorted_hop_constrained_two_hop_label((*L)[it.first], it.second, it.hop, MAX_VALUE, t); // this does not change the size of L[it->first] here, so does not need to lock here
								mtx_2021_increase[it.first].unlock();
								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void PI12(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR,
						  std::vector<hop_constrained_affected_label> &al1_curr, std::vector<hop_constrained_pair_label> *al2_next, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k, int time)
				{

					for (auto it : al1_curr)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, PPR, al2_next, &instance_graph, upper_k]
																		  {

								int v = it.first, u = it.second;
								int hop_u = it.hop;
								mtx_5992[v].lock();
								std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
								mtx_5992[v].unlock();
								temp.push_back(u);

								mtx_2021_increase[v].lock_shared();
								auto Lv = (*L)[v]; // to avoid interlocking
								mtx_2021_increase[v].unlock_shared();

								for (auto t : temp) {

									if (v < t)
									{
										long long d1 = MAX_VALUE;
										int hop_vn = 0;
										for (auto nei : instance_graph[t])
										{
											mtx_2021_increase[nei.first].lock_shared();
											std::pair<weightTYPE, int> dis_hop = get_shortest_distance_hop_two_hop_label2((*L)[nei.first], v);
											mtx_2021_increase[nei.first].unlock_shared();
											if (d1 > dis_hop.first + (long long)nei.second)
											{
												d1 = dis_hop.first + (long long)nei.second;
												hop_vn = dis_hop.second;
											}
										}
										for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++)
										{
											if (hop_i > upper_k)
												break;
											long long di = MAX_VALUE;
											for (auto nei : instance_graph[t]) {
												mtx_2021_increase[nei.first].lock_shared();
												di = std::min(di, search_sorted_hop_constrained_weight_two_hop_label((*L)[nei.first], v, hop_i - 1) + (long long)nei.second);
												mtx_2021_increase[nei.first].unlock_shared();
											}

											mtx_2021_increase[t].lock_shared();
											auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[t], Lv, hop_i);
											mtx_2021_increase[t].unlock_shared();

											if (query_result.first > di) {
												mtx_2021_increase[t].lock();
												insert_sorted_hop_constrained_two_hop_label((*L)[t], v, hop_i, di, time);
												mtx_2021_increase[t].unlock();
												mtx_599_1.lock();
												al2_next->push_back(hop_constrained_pair_label{t, v, hop_i});
												mtx_599_1.unlock();
											}
											else {
												if (query_result.second != -1 && query_result.second != v) {
													mtx_5992[t].lock();
													PPR_TYPE::PPR_insert(*PPR, t, query_result.second, v);
													mtx_5992[t].unlock();
												}
												if (query_result.second != -1 && query_result.second != t) {
													mtx_5992[v].lock();
													PPR_TYPE::PPR_insert(*PPR, v, query_result.second, t);
													mtx_5992[v].unlock();
												}
											}
										}
									}
									if (t < v) {
										long long d1 = MAX_VALUE;
										int hop_vn = 0;
										for (auto nei : instance_graph[v]) {
											mtx_2021_increase[nei.first].lock_shared();
											std::pair<weightTYPE, int> dis_hop = get_shortest_distance_hop_two_hop_label2((*L)[nei.first], t);
											mtx_2021_increase[nei.first].unlock_shared();
											if (d1 > dis_hop.first + (long long)nei.second)
											{
												d1 = dis_hop.first + (long long)nei.second;
												hop_vn = dis_hop.second;
											}
										}

										for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++)
										{
											if (hop_i > upper_k)
												break;
											long long di = MAX_VALUE;
											for (auto nei : instance_graph[v]) {
												mtx_2021_increase[nei.first].lock_shared();
												di = std::min(di, search_sorted_hop_constrained_weight_two_hop_label((*L)[nei.first], t, hop_i - 1) + (long long)nei.second);
												mtx_2021_increase[nei.first].unlock_shared();
											}
											mtx_2021_increase[t].lock_shared();
											auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[t], Lv, hop_i);
											mtx_2021_increase[t].unlock_shared();

											if (query_result.first > di) {
												mtx_2021_increase[v].lock();
												insert_sorted_hop_constrained_two_hop_label((*L)[v], t, hop_i, di, time);
												mtx_2021_increase[v].unlock();
												mtx_599_1.lock();
												al2_next->push_back(hop_constrained_pair_label{v, t, hop_i});
												mtx_599_1.unlock();
											}
											else {
												if (query_result.second != -1 && query_result.second != v) {
													mtx_5992[t].lock();
													PPR_TYPE::PPR_insert(*PPR, t, query_result.second, v);
													mtx_5992[t].unlock();
												}
												if (query_result.second != -1 && query_result.second != t) {
													mtx_5992[v].lock();
													PPR_TYPE::PPR_insert(*PPR, v, query_result.second, t);
													mtx_5992[v].unlock();
												}
											}
										}
									}
								}
								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void PI22(graph<int> &instance_graph, std::vector<std::vector<two_hop_label>> *L, PPR_TYPE::PPR_type *PPR,
						  std::vector<hop_constrained_pair_label> &al2_curr, std::vector<hop_constrained_pair_label> *al2_next, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k, int time)
				{

					for (auto it = al2_curr.begin(); it != al2_curr.end(); it++)
					{
						results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, PPR, al2_next, &instance_graph, upper_k]
																		  {

								mtx_2021_increase[it->second].lock_shared();
								auto Lxx = (*L)[it->second]; // to avoid interlocking
								mtx_2021_increase[it->second].unlock_shared();

								if (it->hop + 1 > upper_k)
									return 1;

								for (auto nei : instance_graph[it->first]) {
									if (nei.first > it->second) {
										mtx_2021_increase[it->first].lock_shared();
										long long search_result = search_sorted_hop_constrained_weight_two_hop_label((*L)[it->first], it->second, it->hop) + (long long)nei.second;
										mtx_2021_increase[it->first].unlock_shared();
										mtx_2021_increase[nei.first].lock_shared();
										auto query_result = graph_weighted_two_hop_extract_distance_and_hub_by_backup_label((*L)[nei.first], Lxx, it->hop + 1);
										mtx_2021_increase[nei.first].unlock_shared();
										if (query_result.first > search_result) {
											mtx_2021_increase[nei.first].lock();
											insert_sorted_hop_constrained_two_hop_label((*L)[nei.first], it->second, it->hop + 1, search_result, time);
											mtx_2021_increase[nei.first].unlock();
											mtx_599_1.lock();
											al2_next->push_back(hop_constrained_pair_label{nei.first, it->second, it->hop + 1});
											mtx_599_1.unlock();
										}
										else {
											if (query_result.second != -1 && query_result.second != it->second) {
												mtx_5992[nei.first].lock();
												PPR_TYPE::PPR_insert(*PPR, nei.first, query_result.second, it->second);
												mtx_5992[nei.first].unlock();
											}
											if (query_result.second != -1 && query_result.second != nei.first) {
												mtx_5992[it->second].lock();
												PPR_TYPE::PPR_insert(*PPR, it->second, query_result.second, nei.first);
												mtx_5992[it->second].unlock();
											}
										}
									}
								}

								return 1; }));
					}

					for (auto &&result : results_dynamic)
					{
						result.get();
					}
					std::vector<std::future<int>>().swap(results_dynamic);
				}

				void HOP_WeightIncrease2021_batch(graph<int> &instance_graph, two_hop_case_info &mm,
												  std::vector<std::pair<int, int>> &v, std::vector<weightTYPE> &w_old_vec,
												  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
				{
					std::map<std::pair<int, int>, weightTYPE> w_old_map;
					int batch_size = v.size();
					for (int i = 0; i < batch_size; i++)
					{
						if (v[i].first > v[i].second)
						{
							std::swap(v[i].first, v[i].second);
						}
						if (w_old_map.count(v[i]) == 0)
						{
							w_old_map[v[i]] = w_old_vec[i];
						}
					}

					std::vector<hop_constrained_affected_label> al1_curr, al1_next;
					std::vector<hop_constrained_pair_label> al2_curr, al2_next;

					for (auto &iter : w_old_map)
					{
						int v1 = iter.first.first;
						int v2 = iter.first.second;
						weightTYPE w_old = iter.second;
						for (auto it : mm.L[v1])
						{
							long long search_weight = search_sorted_hop_constrained_weight_two_hop_label(mm.L[v2], it.hub_vertex, it.hop + 1);
							if (it.hub_vertex <= v2 && search_weight >= (long long)it.distance + w_old && search_weight < MAX_VALUE && it.t_e == std::numeric_limits<int>::max())
							{
								al1_curr.push_back(hop_constrained_affected_label(v2, it.hub_vertex, it.hop + 1, it.distance + w_old));
							}
						}
						for (auto it : mm.L[v2])
						{
							long long search_weight = search_sorted_hop_constrained_weight_two_hop_label(mm.L[v1], it.hub_vertex, it.hop + 1);
							if (it.hub_vertex <= v1 && search_weight >= (long long)it.distance + w_old && search_weight < MAX_VALUE && it.t_e == std::numeric_limits<int>::max())
							{
								al1_curr.push_back(hop_constrained_affected_label(v1, it.hub_vertex, it.hop + 1, it.distance + w_old));
							}
						}
					}
					while (al1_curr.size() || al2_curr.size())
					{
						std::cout << "al1 curr size is " << al1_curr.size() << " al2 cur size is " << al2_curr.size() << std::endl;
						PI11(instance_graph, &mm.L, al1_curr, &al1_next, w_old_map, pool_dynamic, results_dynamic, time);
						std::cout << "al1 next size is " << al1_next.size() << " al2 next size is " << al2_next.size() << std::endl;
						PI12(instance_graph, &mm.L, &mm.PPR, al1_curr, &al2_next, pool_dynamic, results_dynamic, mm.upper_k, time);
						std::cout << "al1 next size is " << al1_next.size() << " al2 next size is " << al2_next.size() << std::endl;
						PI22(instance_graph, &mm.L, &mm.PPR, al2_curr, &al2_next, pool_dynamic, results_dynamic, mm.upper_k, time);
						std::cout << "al1 next size is " << al1_next.size() << " al2 next size is " << al2_next.size() << std::endl;

						al1_curr = al1_next;
						al2_curr = al2_next;
						std::vector<hop_constrained_affected_label>().swap(al1_next);
						std::vector<hop_constrained_pair_label>().swap(al2_next);
					}
				}

			}
		}
	}
}