#pragma once

#include <map>
#include "CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_maintain_PLL.h"
using namespace std;
namespace nonHop
{
	void ProDecreasep_batch(graph_v_of_v<int> &instance_graph, vector<vector<two_hop_label>> *L, PPR_type *PPR,
							std::vector<affected_label> &CL_curr, std::vector<affected_label> *CL_next, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
	{

		for (auto it : CL_curr)
		{
			results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, PPR, CL_next, &instance_graph]
															  {

			if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
				throw reach_limit_time_string;
			}

			int v = it.first, u = it.second;

			mtx_595[u].lock();
			auto Lu = (*L)[u]; // to avoid interlocking
			mtx_595[u].unlock();

			for (auto nei : instance_graph[v]) {
				int vnei = nei.first;
				long long dnew = it.dis + nei.second;

				if (u < vnei) {
					mtx_595[vnei].lock();
					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[vnei], Lu); // query_result is {distance, common hub}
					mtx_595[vnei].unlock();
					if (query_result.first > dnew) {
						mtx_595[vnei].lock();
						insert_sorted_two_hop_label((*L)[vnei], u, dnew, time);
						mtx_595[vnei].unlock();
						mtx_595_1.lock();
						CL_next->push_back(affected_label(vnei, u, dnew));
						mtx_595_1.unlock();
					}
					else {
						mtx_595[vnei].lock();
						auto search_result = search_sorted_two_hop_label2((*L)[vnei], u);
						mtx_595[vnei].unlock();
						if (search_result.first < MAX_VALUE && search_result.first > dnew) {
							mtx_595[vnei].lock();
							(*L)[vnei][search_result.second].distance = dnew;
							mtx_595[vnei].unlock();
							mtx_595_1.lock();
							CL_next->push_back(affected_label(vnei, u, dnew));
							mtx_595_1.unlock();
						}
						if (query_result.second != u) {
							mtx_5952[vnei].lock();
							PPR_insert(*PPR, vnei, query_result.second, u);
							mtx_5952[vnei].unlock();
						}
						if (query_result.second != vnei) {
							mtx_5952[u].lock();
							PPR_insert(*PPR, u, query_result.second, vnei);
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

	void nonHOP_WeightDecrease2021_batch(graph_v_of_v<int> &instance_graph, two_hop_case_info &mm, std::vector<pair<int, int>> &v, std::vector<weightTYPE> &w_new,
										 ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time)
	{

		label_operation_times = 0;
		global_query_times = 0;

		begin_time = std::chrono::high_resolution_clock::now();
		max_run_time_nanosec = _2023algo_max_second * 1e9;

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
			weightTYPE w_new = it.second;
			for (int sl = 0; sl < 2; sl++)
			{
				if (sl == 1)
				{
					swap(v1, v2);
				}
				for (auto it : L[v1])
				{
					int v = it.vertex;
					long long dis = it.distance + w_new;
					if (v <= v2)
					{
						auto query_result = Query2(v, v2); // query_result is {distance, common hub}
						if (query_result.first > dis)
						{
							insert_sorted_two_hop_label(L[v2], v, dis, time);
							CL_curr.push_back(affected_label(v2, v, dis));
						}
						else
						{
							auto search_result = search_sorted_two_hop_label2(L[v2], v);
							if (search_result.first < MAX_VALUE && search_result.first > dis)
							{
								L[v2][search_result.second].distance = dis;
								CL_curr.push_back(affected_label(v2, v, dis));
							}
							if (query_result.second != v)
							{
								PPR_insert(mm.PPR, v2, query_result.second, v);
							}
							if (query_result.second != v2)
							{
								PPR_insert(mm.PPR, v, query_result.second, v2);
							}
						}
					}
				}
			}
		}

		while (CL_curr.size())
		{
			ProDecreasep_batch(instance_graph, &mm.L, &mm.PPR, CL_curr, &CL_next, pool_dynamic, results_dynamic, time);
			CL_curr = CL_next;
			std::vector<affected_label>().swap(CL_next);
		}
	}
}
