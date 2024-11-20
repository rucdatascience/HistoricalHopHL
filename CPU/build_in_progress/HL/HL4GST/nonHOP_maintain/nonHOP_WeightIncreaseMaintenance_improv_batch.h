#pragma once

#include <CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_maintain_PLL.h>

#include <queue>
#include <map>

using namespace std;
using namespace nonHop;

namespace nonHop
{
	void SPREAD1_batch(graph_v_of_v<int>& instance_graph, vector<vector<two_hop_label>>* L,
		std::vector<affected_label>& al1, std::vector<pair_label>* al2, std::map<pair<int,int>,weightTYPE >& w_old_map,
		ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic,int time) {

		for (auto &it : al1) {
			results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, al2, &instance_graph, &w_old_map,time] {
				queue<pair<int, weightTYPE> > q; //(u,d)
				int v = it.second;
				q.push(pair<int, weightTYPE>(it.first, it.dis));
				while (!q.empty()) {
					int x = q.front().first;
					weightTYPE dx = q.front().second;
					q.pop();		
					insert_sorted_two_hop_label((*L)[x], v, MAX_VALUE,time); 
					mtx_595_1.lock();
					al2->push_back(pair_label(x, v));
					mtx_595_1.unlock();
					for (auto nei : instance_graph[x]) {
						if (v < nei.first) {
							weightTYPE search_weight = search_sorted_two_hop_label((*L)[nei.first], v);
							weightTYPE w_old;
							if(w_old_map.count(pair<int,int>(x,nei.first))>0){
								w_old=w_old_map[pair<int,int>(x,nei.first)];
							}
							else if(w_old_map.count(pair<int,int>(nei.first,x))>0){
								w_old=w_old_map[pair<int,int>(nei.first,x)];
							}
							else{
								w_old=nei.second;
							}
							if (dx + w_old == search_weight && search_weight < MAX_VALUE) {
								q.push(pair<int, weightTYPE>(nei.first, dx + w_old));
							}
						}
					}
				}

				return 1; }));
		}

		for (auto&& result : results_dynamic) {
			result.get();
		}
		std::vector<std::future<int>>().swap(results_dynamic);
	}

	void SPREAD2_batch(graph_v_of_v<int>& instance_graph, vector<vector<two_hop_label>>* L, PPR_type* PPR,
		std::vector<pair_label>& al2, std::vector<affected_label>* al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic,int time) {

		for (auto &it : al2) {
			results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al3, &instance_graph] {

				int v = it.first, u = it.second;
				mtx_5952[v].lock_shared();
				std::vector<int> temp = PPR_retrieve(*PPR, v, u);
				mtx_5952[v].unlock_shared();
				PPR_binary_operations_insert(temp, u);
				for (auto t : temp) {
					if (v < t) {
						long long d1 = MAX_VALUE;
						for (auto nei : instance_graph[t]) {
							d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], v) + (long long)nei.second);
						}
						auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, t, v);
						if(d1 >= 2e6) continue;
						if (query_result.first > d1) { // only add new label when it's absolutely necessary
							mtx_595_1.lock();
							al3->push_back(affected_label(t, v, d1));
							mtx_595_1.unlock();
						}
						else {
							if (query_result.second != v) {
								mtx_5952[t].lock();
								PPR_insert(*PPR, t, query_result.second, v);
								mtx_5952[t].unlock();
							}
							if (query_result.second != t) {
								mtx_5952[v].lock();
								PPR_insert(*PPR, v, query_result.second, t);
								mtx_5952[v].unlock();
							}
						}
					}
					else if (t < v) {
						long long d1 = MAX_VALUE;
						for (auto nei : instance_graph[v]) {
							d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], t) + (long long)nei.second);
						}
						if(d1 >= 2e6) continue;
						auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, v, t);
						if (query_result.first > d1) {
							mtx_595_1.lock();
							al3->push_back(affected_label(v, t, d1));
							mtx_595_1.unlock();
						}
						else {
							if (query_result.second != v) {
								mtx_5952[t].lock();
								PPR_insert(*PPR, t, query_result.second, v);
								mtx_5952[t].unlock();
							}
							if (query_result.second != t) {
								mtx_5952[v].lock();
								PPR_insert(*PPR, v, query_result.second, t);
								mtx_5952[v].unlock();
							}
						}
					}
				}

				return 1; }));
		}
		for (auto&& result : results_dynamic) {
			result.get();
		}
		std::vector<std::future<int>>().swap(results_dynamic);
	}

	void SPREAD3_batch(graph_v_of_v<int>& instance_graph, vector<vector<two_hop_label>>* L, PPR_type* PPR, std::vector<affected_label>& al3,
	ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic,int time) {

		// Deduplication (u,v,dis)
		std::map<std::pair<int,int>,weightTYPE> al3_edge_map;
		for(auto &it:al3){
			if(al3_edge_map.count({it.first,it.second})==0){
				al3_edge_map[{it.first,it.second}]=it.dis;
			}
			else if(al3_edge_map[{it.first,it.second}]>it.dis){
				al3_edge_map[{it.first,it.second}]=it.dis;
			}
		}

		// extract each unique hub v and its (u,dis) list
		std::map<int,std::vector<std::pair<int,weightTYPE>>> al3_map; // al3_map[v]=(u1,dis1),(u2,dis2)...
		for(auto &it:al3_edge_map){
			int u = it.first.first;
			int v = it.first.second;
			weightTYPE dis = it.second;
			if(al3_map.count(v)==0){
				std::vector<std::pair<int, weightTYPE>> vec_with_hub_v;
				vec_with_hub_v.emplace_back(make_pair(u, dis));
				al3_map[v] = vec_with_hub_v;
			}
			else{
				std::vector<std::pair<int, weightTYPE>> vec_with_hub_v = al3_map[v];
				vec_with_hub_v.emplace_back(make_pair(u, dis));
				al3_map[v] = vec_with_hub_v;
			}
		}

		std::vector<std::pair<int,std::vector<std::pair<int,weightTYPE>>>> al3_map_vec(al3_map.begin(),al3_map.end());
		sort(al3_map_vec.begin(),al3_map_vec.end(),[](const std::pair<int,std::vector<std::pair<int,weightTYPE>>>& a, const std::pair<int,std::vector<std::pair<int,weightTYPE>>>& b){
			return a.first<b.first;
		});

		// std::cout<<"SPREAD3_batch"<<std::endl;
		for (auto &it : al3_map_vec) {
			results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, &instance_graph, PPR,time] {

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

				vector<int> Dis_changed;
				auto& DIS = Dis[current_tid];
				auto& Q_HANDLES = Q_handles[current_tid];
				auto& Q_VALUE = Q_value[current_tid];
				// int N = instance_graph.size();
				// std::vector<std::pair<weightTYPE, int>> DIS(N, { -1, -1 });
				// std::vector<handle_t_for_DIFFUSE> Q_HANDLES(N);
				// std::vector<weightTYPE> Q_VALUE(N, MAX_VALUE);
				boost::heap::fibonacci_heap<node_for_DIFFUSE> pq;

				for(auto &it:vec_with_hub_v){
					int u = it.first;
					weightTYPE du = it.second;
					// std::cout<<"u: "<<u<<" v: "<<v<<" du: "<<du<<"\n";
					mtx_595[u].lock_shared();
					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[u], Lv);
					mtx_595[u].unlock_shared();
					bool flag = false;
					if (query_result.first < du) {
						if (query_result.second != v) {
							mtx_5952[u].lock();
							PPR_insert(*PPR, u, query_result.second, v);
							mtx_5952[u].unlock();
							flag = true;
						}
						if (query_result.second != u) {
							mtx_5952[v].lock();
							PPR_insert(*PPR, v, query_result.second, u);
							mtx_5952[v].unlock();
							flag = true;
						}

						// mtx_595_1.lock();
						// Qid_595.push(current_tid);
						// mtx_595_1.unlock();
						// return 1;
					}
					
					if(flag == true){
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
					weightTYPE d_old = search_sorted_two_hop_label((*L)[x], v);
					// std::cout<<"insert: "<<x<<" "<<v<<" "<<dx<<" "<<d_old<<"\n";
					if (dx < d_old) {
						insert_sorted_two_hop_label((*L)[x], v, dx,time);
					}
					else {
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
								DIS[xnei] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[xnei], Lv);
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
									PPR_insert(*PPR, xnei, DIS[xnei].second, v);
									mtx_5952[xnei].unlock();
								}
								if (DIS[xnei].second != xnei) {
									mtx_5952[v].lock();
									PPR_insert(*PPR, v, DIS[xnei].second, xnei);
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

		for (auto&& result : results_dynamic) {
			result.get();
		}
		std::vector<std::future<int>>().swap(results_dynamic);
	}


	void nonHOP_WeightIncreaseMaintenance_improv_batch(graph_v_of_v<int>& instance_graph, two_hop_case_info& mm, vector<pair<int,int> >& v, vector<int>& w_old_vec,
	ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic,int time) {

		label_operation_times = 0;
		global_query_times = 0;

		std::vector<affected_label> al1, al3;
		std::vector<pair_label> al2;
		std::map<pair<int,int>,weightTYPE > w_old_map;
		int batch_size = v.size();
		for(int i=0;i<batch_size;i++){
			if(v[i].first>v[i].second){
				swap(v[i].first,v[i].second);
			}
			if(w_old_map.count(v[i])==0){
				w_old_map[v[i]]=w_old_vec[i];
			}
		}

		for(auto &it:w_old_map){
			results_dynamic.emplace_back(pool_dynamic.enqueue([it, &al1, &instance_graph, &mm, &w_old_map] {
				int v1=it.first.first;
				int v2=it.first.second;
				weightTYPE w_old=it.second;
				for (auto it : mm.L[v1]) {
					long long search_weight = search_sorted_two_hop_label(mm.L[v2], it.vertex);
					if (it.vertex <= v2 && search_weight == (long long)it.distance + w_old && search_weight < MAX_VALUE) {
						mtx_595_1.lock();
						al1.push_back(affected_label(v2, it.vertex, it.distance + w_old));
						mtx_595_1.unlock();
					}
				}
				for (auto it : mm.L[v2]) {
					long long search_weight = search_sorted_two_hop_label(mm.L[v1], it.vertex);
					if (it.vertex <= v1 && search_weight == (long long)it.distance + w_old && search_weight < MAX_VALUE) {
						mtx_595_1.lock();
						al1.push_back(affected_label(v1, it.vertex, it.distance + w_old));
						mtx_595_1.unlock();
					}
				}
			return 1; }));
		}

		for (auto&& result : results_dynamic) {
			result.get();
		}
		std::vector<std::future<int>>().swap(results_dynamic);

		SPREAD1_batch(instance_graph, &mm.L, al1, &al2, w_old_map, pool_dynamic, results_dynamic, time);
		SPREAD2_batch(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic, time);
		SPREAD3_batch(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic, time);

	}
}