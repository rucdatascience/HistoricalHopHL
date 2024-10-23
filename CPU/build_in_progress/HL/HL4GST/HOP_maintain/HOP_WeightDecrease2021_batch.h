#pragma once
using namespace std;
#include "CPU/tool_functions/ThreadPool.h"
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h>
#include <algorithm>
#include <map>

void ProDecreasep_batch(graph_v_of_v<int> &instance_graph, vector<vector<hop_constrained_two_hop_label>> *L, PPR_type *PPR,
                        std::vector<hop_constrained_affected_label> &CL_curr, std::vector<hop_constrained_affected_label> *CL_next,
                        ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k, int t)
{

    for (auto it : CL_curr)
    {
        results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, PPR, CL_next, &instance_graph, upper_k]
                                                          {

            if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
                throw reach_limit_time_string;
            }

			int v = it.first, u = it.second;
    
			mtx_599[u].lock();
			auto Lu = (*L)[u]; // to avoid interlocking
			mtx_599[u].unlock();

            if(it.hop + 1 > upper_k)
                return 1;

            for (auto nei : instance_graph[v]) {
				int vnei = nei.first;
                int hop_u = it.hop;
                long long dnew = it.dis + nei.second;
				if (u < vnei) {
					mtx_599[vnei].lock();
					auto query_result = hop_constrained_extract_distance_and_hub((*L),vnei, u, hop_u + 1); // query_result is {distance, common hub}
                    mtx_599[vnei].unlock();
					if ((long long)query_result.first > dnew) {
						mtx_599[vnei].lock();
						insert_sorted_hop_constrained_two_hop_label((*L)[vnei], u, hop_u + 1, dnew, t);
						mtx_599[vnei].unlock();
						mtx_599_1.lock();
						CL_next->push_back(hop_constrained_affected_label(vnei, u, hop_u + 1, dnew));
						mtx_599_1.unlock();
					}
					else {
						mtx_599[vnei].lock();
						auto search_result = search_sorted_hop_constrained_two_hop_label_and_index((*L)[vnei], u, hop_u + 1);
						mtx_599[vnei].unlock();
						if (search_result.first < MAX_VALUE && search_result.first > dnew) {
							mtx_599[vnei].lock();
							(*L)[vnei][search_result.second].distance = dnew;
							mtx_599[vnei].unlock();
							mtx_599_1.lock();
							CL_next->push_back(hop_constrained_affected_label(vnei, u, hop_u + 1, dnew));
							mtx_599_1.unlock();
						}
						if (query_result.second != u) {
							mtx_5992[vnei].lock();
							PPR_insert(*PPR, vnei, query_result.second, u);
							mtx_5992[vnei].unlock();
						}
						if (query_result.second != vnei) {
							mtx_5992[u].lock();
							PPR_insert(*PPR, u, query_result.second, vnei);
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

void HOP_WeightDecrease2021_batch(graph_v_of_v<int> &instance_graph, hop_constrained_case_info &mm, std::vector<pair<int, int>> &v, std::vector<weightTYPE> &w_new,
                                  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
{

    global_query_times = 0;
    label_operation_times = 0;

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
                swap(v1, v2);
            }
            for (auto it : L[v1])
            {
                int v = it.hub_vertex;
                int hop_v = it.hop;
                long long dis = it.distance + w_new;
                if (v <= v2)
                {
                    auto query_result = hop_constrained_extract_distance_and_hub(L, v, v2, hop_v + 1); // query_result is {distance, common hub}

                    if ((long long)query_result.first > dis)
                    {
                        insert_sorted_hop_constrained_two_hop_label(L[v2], v, hop_v + 1, dis, t);
                        CL_curr.push_back(hop_constrained_affected_label(v2, v, hop_v + 1, dis));
                    }
                    else
                    {
                        auto search_result = search_sorted_hop_constrained_two_hop_label_and_index(L[v2], v, hop_v + 1);
                        if (search_result.first < MAX_VALUE && search_result.first > dis)
                        {
                            L[v2][search_result.second].distance = dis;
                            CL_curr.push_back(hop_constrained_affected_label(v2, v, hop_v + 1, dis));
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
        ProDecreasep_batch(instance_graph, &mm.L, &mm.PPR, CL_curr, &CL_next, pool_dynamic, results_dynamic, mm.upper_k, t);
        CL_curr = CL_next;
        std::vector<hop_constrained_affected_label>().swap(CL_next);
    }
}
