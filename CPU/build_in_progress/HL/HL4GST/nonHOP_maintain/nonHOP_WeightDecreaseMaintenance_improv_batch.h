#pragma once
#include <map>
#include "CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_maintain_PLL.h"
using namespace std;
namespace nonHop
{
	void WeightDecreaseMaintenance_improv_step1_batch(std::map<pair<int, int>, weightTYPE> &v_map, vector<vector<two_hop_label>> *L, PPR_type *PPR, std::vector<affected_label> *CL,
													  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
	{

		for (auto it : v_map)
		{
			results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, CL]
															  {
			int v1 = it.first.first, v2 = it.first.second;
			weightTYPE w_new = it.second;
			for (int sl = 0; sl < 2; sl++) {
				if (sl == 1) {
					swap(v1, v2);
				}
				for (auto &it : (*L)[v1]) {
					if (it.vertex <= v2 && it.distance + w_new < 2e6&& it.t_e == std::numeric_limits<int>::max()) {
						auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, it.vertex, v2); // query_result is {distance, common hub}
						if (query_result.first > it.distance + w_new) {
							mtx_595_1.lock();
							CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
							mtx_595_1.unlock();
						}
						else {
							auto search_result = search_sorted_two_hop_label((*L)[v2], it.vertex);
							if (search_result > it.distance + w_new && search_result < MAX_VALUE) {
								mtx_595_1.lock();
								CL->push_back(affected_label{ v2, it.vertex, it.distance + w_new });
								mtx_595_1.unlock();
							}
							if (query_result.second != it.vertex) {
								mtx_5952[v2].lock();
								PPR_insert(*PPR, v2, query_result.second, it.vertex);
								mtx_5952[v2].unlock();
							}
							if (query_result.second != v2) {
								mtx_5952[it.vertex].lock();
								PPR_insert(*PPR, it.vertex, query_result.second, v2);
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
	}

	void DIFFUSE_batch(graph_v_of_v<int> &instance_graph, vector<vector<two_hop_label>> *L, PPR_type *PPR, std::vector<affected_label> &CL,
					   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
	{

		// Deduplication
		std::map<std::pair<int, int>, weightTYPE> CL_edge_map;
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
		std::map<int, std::vector<std::pair<int, weightTYPE>>> CL_map; // CL_map[v]=(u1,dis1),(u2,dis2)...
		for (auto &it : CL_edge_map)
		{
			int u = it.first.first;
			int v = it.first.second;
			weightTYPE dis = it.second;
			if (CL_map.count(v) == 0)
			{
				std::vector<std::pair<int, weightTYPE>> vec_with_hub_v;
				vec_with_hub_v.emplace_back(make_pair(u, dis));
				CL_map[v] = vec_with_hub_v;
			}
			else
			{
				std::vector<std::pair<int, weightTYPE>> vec_with_hub_v = CL_map[v];
				vec_with_hub_v.emplace_back(make_pair(u, dis));
				CL_map[v] = vec_with_hub_v;
			}
		}

		std::vector<std::pair<int, std::vector<std::pair<int, weightTYPE>>>> CL_map_vec(CL_map.begin(), CL_map.end());
		sort(CL_map_vec.begin(), CL_map_vec.end(), [](const std::pair<int, std::vector<std::pair<int, weightTYPE>>> &a, const std::pair<int, std::vector<std::pair<int, weightTYPE>>> &b)
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
			std::vector<std::pair<int, weightTYPE>> vec_with_hub_v = it.second;

			// int u = it.first.first, v = it.first.second;
			// weightTYPE du = it.second;
			mtx_595[v].lock_shared();
			auto Lv = (*L)[v]; // to avoid interlocking
			mtx_595[v].unlock_shared();

			vector<int> Dis_changed;
			auto& DIS = Dis[current_tid];
			auto& Q_HANDLES = Q_handles[current_tid];
			auto& Q_VALUE = Q_value[current_tid];
			// int N=instance_graph.size();
			// std::vector<std::pair<weightTYPE, int>> DIS(N, { -1, -1 });
			// std::vector<handle_t_for_DIFFUSE> Q_HANDLES(N);
			// std::vector<weightTYPE> Q_VALUE(N, MAX_VALUE);
			boost::heap::fibonacci_heap<node_for_DIFFUSE> Q;

			for(auto &it:vec_with_hub_v){
				int u = it.first;
				weightTYPE du = it.second;
				DIS[u] = { du, v }; // <distance, hub responsible for this distance>
				Dis_changed.push_back(u);
				Q_HANDLES[u] = Q.push(node_for_DIFFUSE(u, du));
				Q_VALUE[u] = du;
			}

			while (!Q.empty()) {

				node_for_DIFFUSE temp2 = Q.top();
				int x = temp2.index;
				weightTYPE dx = temp2.disx;
				Q.pop();
				Q_VALUE[x] = MAX_VALUE;

				mtx_595[x].lock();
				weightTYPE d_old=search_sorted_two_hop_label2((*L)[x], v).first;
				if(d_old>dx){
					insert_sorted_two_hop_label((*L)[x], v, dx,t);
				}
				else{
					continue;
					//dx=d_old;
				}
				mtx_595[x].unlock();

				for (auto& nei : instance_graph[x]) {
					int xnei = nei.first;
					weightTYPE d_new = dx + nei.second;

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
								Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
							}
							else {
								Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
							}
							Q_VALUE[xnei] = d_new;
						}
						else {
							mtx_595[xnei].lock_shared();
							auto search_result = search_sorted_two_hop_label2((*L)[xnei], v);
							mtx_595[xnei].unlock_shared();
							if (search_result.second != -1 && std::min(search_result.first, Q_VALUE[xnei]) > d_new) {
								if (Q_VALUE[xnei] >= MAX_VALUE) {
									Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
								}
								else {
									Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
								}
								Q_VALUE[xnei] = d_new;
							}
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

		for (auto &&result : results_dynamic)
		{
			result.get();
		}
		std::vector<std::future<int>>().swap(results_dynamic);
	}

	void nonHOP_WeightDecreaseMaintenance_improv_batch(graph_v_of_v<int> &instance_graph, two_hop_case_info &mm, std::vector<pair<int, int>> &v, std::vector<int> &w_new,
													   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
	{

		label_operation_times = 0;
		global_query_times = 0;

		std::map<pair<int, int>, weightTYPE> w_new_map;
		int batch_size = v.size();
		for (int i = 0; i < batch_size; i++)
		{
			if (v[i].first > v[i].second)
			{
				swap(v[i].first, v[i].second);
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

		WeightDecreaseMaintenance_improv_step1_batch(w_new_map, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic, time);

		DIFFUSE_batch(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic, time);
	}
}