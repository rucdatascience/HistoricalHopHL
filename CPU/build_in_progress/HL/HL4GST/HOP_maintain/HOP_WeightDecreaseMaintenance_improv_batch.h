#pragma once
using namespace std;
#include "CPU/tool_functions/ThreadPool.h"
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h>
#include <algorithm>
#include <map>
void WeightDecreaseMaintenance_improv_step1_batch(std::map<pair<int, int>, weightTYPE> &v_map, vector<vector<hop_constrained_two_hop_label>> *L, PPR_type *PPR, std::vector<hop_constrained_affected_label> *CL,
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
                    swap(v1, v2);
                }
                for (auto it : (*L)[v1])
                {
                    if (it.hub_vertex <= v2 && (long long int)it.distance + w_new < TwoM_value && it.t_e == std::numeric_limits<int>::max())
                    {
                            auto query_result = hop_constrained_extract_distance_and_hub(*L, it.hub_vertex, v2, it.hop + 1); // query_result is {distance, common hub}
                            if ((long long int)query_result.first > (long long int)it.distance + w_new)
                            {
                                mtx_599_1.lock();
                                CL->push_back(hop_constrained_affected_label{v2, it.hub_vertex, it.hop + 1, it.distance + w_new});
                                mtx_599_1.unlock();
                            }
                            else
                            {
                                auto search_result = search_sorted_hop_constrained_two_hop_label((*L)[v2], it.hub_vertex, it.hop + 1);
                                if (search_result < MAX_VALUE && search_result > it.distance + w_new)
                                {
                                    mtx_599_1.lock();
                                    CL->push_back(hop_constrained_affected_label{v2, it.hub_vertex, it.hop + 1, it.distance + w_new});
                                    mtx_599_1.unlock();
                                }
                                if (query_result.second != -1 && query_result.second != it.hub_vertex)
                                {
                                    mtx_5992[v2].lock();
                                    PPR_insert(*PPR, v2, query_result.second, it.hub_vertex);
                                    mtx_5992[v2].unlock();
                                }
                                if (query_result.second != -1 && query_result.second != v2)
                                {
                                    mtx_5992[it.hub_vertex].lock();
                                    PPR_insert(*PPR, it.hub_vertex, query_result.second, v2);
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

void DIFFUSE_batch(graph_v_of_v<int> &instance_graph, vector<vector<hop_constrained_two_hop_label>> *L, PPR_type *PPR, std::vector<hop_constrained_affected_label> &CL,
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

			mtx_599[v].lock_shared();
			auto Lv = (*L)[v]; // to avoid interlocking
			mtx_599[v].unlock_shared();

			vector<int> dist_hop_changes;
			auto & dist_hop = dist_hop_599_v2[current_tid];
			boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE> pq;
			map<pair<int, int>, pair<hop_constrained_handle_t_for_DIFFUSE, int>> Q_handle;
            vector<int> hubs;
			hubs.resize(instance_graph.size(), -1);
            auto& Q_VALUE = Q_value[current_tid];
            
            for(auto &it:vec_with_hub_v){
				int u = it.hub_vertex;
                int h_v = it.hop;
                weightTYPE du = it.distance;

				dist_hop[u] = {du, h_v}; //  {dis, hop}
                dist_hop_changes.push_back(u);
                hop_constrained_node_for_DIFFUSE tmp;
                tmp.index = u;
                tmp.hop = h_v;
                tmp.disx = du;
                Q_handle[{u, h_v}] = {pq.push({tmp}), du}; //{node_for_DIFFUSE_v2,dis}
                // what the meaning of Q_VALUE? mark the data of pq
                if(h_v <= upper_k)
                    Q_VALUE[u][h_v] = du;
                
            }
            
			while (!pq.empty())
			{
				int x = pq.top().index;
				int xhv = pq.top().hop;
				weightTYPE dx = pq.top().disx;
				pq.pop();
                if(xhv <= upper_k)
                    Q_VALUE[x][xhv] = MAX_VALUE;

				mtx_599[x].lock();
				weightTYPE d_old = search_sorted_hop_constrained_two_hop_label((*L)[x], v, xhv);
				if (dx >= 0 && dx < d_old)
				{
                    insert_sorted_hop_constrained_two_hop_label((*L)[x], v, xhv, dx,t);
				}
				mtx_599[x].unlock();
                
                if(xhv + 1 > upper_k)
					continue;

				for (auto nei : instance_graph[x])
				{
					
                    if (dx + nei.second >= TwoM_value)
					    continue;
				
					int xnei = nei.first;
                    int hop_nei = xhv + 1;
                    long long int d_new = dx + (long long int)nei.second;
					hop_constrained_node_for_DIFFUSE node = {xnei, hop_nei, (weightTYPE)d_new};

					if (v < xnei)
					{

						if (dist_hop[xnei].first == -1)
						{
							
							// Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                            // Q_VALUE[xnei][hop_nei] = d_new;
                            mtx_599[xnei].lock_shared();
							// std::pair<int, int> tmp = hop_constrained_extract_distance_and_hub_2((*L)[xnei], Lv, xhv + 1); 
							std::pair<int, int> temp_dis = hop_constrained_extract_distance_and_hop(*L,xnei, v, xhv + 1); 
							mtx_599[xnei].unlock_shared();
							// hubs[xnei] = tmp.second;

                            dist_hop[xnei].first = temp_dis.first;
                            dist_hop[xnei].second = temp_dis.second;							
                            dist_hop_changes.push_back(xnei);
						}
                        
                        if (d_new < dist_hop[xnei].first)
						{
                            
                            //if (Q_handle.find({xnei, hop_nei}) != Q_handle.end())
							if(Q_VALUE[xnei][hop_nei] < MAX_VALUE)
                            {
								if (Q_handle[{xnei, hop_nei}].second > d_new)
								{
									pq.update(Q_handle[{xnei, hop_nei}].first, node);
									Q_handle[{xnei, hop_nei}].second = d_new;
								}
							}
							else
							{
								Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
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
							if(Q_VALUE[xnei][hop_nei] < MAX_VALUE)
                            {
								if (Q_handle[{xnei, hop_nei}].second > d_new)
								{
									pq.update(Q_handle[{xnei, hop_nei}].first, node);
									Q_handle[{xnei, hop_nei}].second = d_new;
								}
							}
							else
							{
								Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
							}
                            Q_VALUE[xnei][hop_nei] = d_new;
						}

						
                        if(dist_hop[xnei].first < d_new)
                        {
                            if (hubs[xnei] != -1 && hubs[xnei] != v)
							{
								mtx_5992[xnei].lock();
								PPR_insert(*PPR, xnei, hubs[xnei], v);
								mtx_5992[xnei].unlock();
							}
							if (hubs[xnei] != -1 && hubs[xnei] != xnei)
							{
								mtx_5992[v].lock();
								PPR_insert(*PPR, v, hubs[xnei], xnei);
								mtx_5992[v].unlock();
							}
                        }
							

					}
				}
			}


			for(int i : dist_hop_changes)
			{
				dist_hop[i] = {-1, 0};
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

void HOP_WeightDecreaseMaintenance_improv_batch(graph_v_of_v<int> &instance_graph, hop_constrained_case_info &mm,
                                                std::vector<pair<int, int>> &v, std::vector<int> &w_new, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
{

    global_query_times = 0;
    label_operation_times = 0;

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
    std::vector<hop_constrained_affected_label> CL;
    mm.mark_time("start get the affected label");
    WeightDecreaseMaintenance_improv_step1_batch(w_new_map, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic, t);
    mm.mark_time("end get the affected label and current CL size is" + std::to_string(CL.size()));
    mm.mark_time("start diffuse");
    DIFFUSE_batch(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic, mm.upper_k, t);
    mm.mark_time("end diffuse");
}
