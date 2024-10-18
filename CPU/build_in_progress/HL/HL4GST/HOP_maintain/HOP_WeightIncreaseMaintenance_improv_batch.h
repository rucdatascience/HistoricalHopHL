#pragma once
using namespace std;
#include "CPU/tool_functions/ThreadPool.h"
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h>
#include <map>
#include <algorithm>

// #define MAX_VALUE std::numeric_limits<int>::max()

void HOP_maintain_SPREAD1_batch(graph_v_of_v<int> &instance_graph, vector<vector<hop_constrained_two_hop_label>> *L,
								std::vector<hop_constrained_affected_label> &al1, std::vector<hop_constrained_pair_label> *al2, std::map<pair<int, int>, int> &w_old_map, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
{

	for (auto it : al1)
	{
		results_dynamic.emplace_back(pool_dynamic.enqueue([t, it, L, al2, &instance_graph, &w_old_map]
														  {
				queue<hop_constrained_node_for_DIFFUSE> q; //(u,h_v, d)
				int v = it.second;
				q.push(hop_constrained_node_for_DIFFUSE(it.first, it.hop, it.dis));
				while (!q.empty()) {
					int x = q.front().index;
					int h_x = q.front().hop;
					int dx = q.front().disx;
					q.pop();
					insert_sorted_hop_constrained_two_hop_label((*L)[x], v, h_x, MAX_VALUE,t); // this does not change the size of L[x] here, so does not need to lock here
					mtx_599_1.lock();
					al2->push_back(hop_constrained_pair_label(x, v, h_x));
					mtx_599_1.unlock();
	
					for (auto nei : instance_graph[x])
					{
						if (v < nei.first) {
							int search_weight = search_sorted_hop_constrained_two_hop_label((*L)[nei.first], v, h_x + 1);
							int w_old = nei.second;
							if (w_old_map.count(pair<int, int>(x, nei.first)) > 0) {
								w_old = w_old_map[pair<int, int>(x, nei.first)];
							}
							else if (w_old_map.count(pair<int, int>(nei.first, x)) > 0) {
								w_old = w_old_map[pair<int, int>(nei.first, x)];
							}
							else {
								w_old = nei.second;
							}
							if (dx + w_old <= search_weight && search_weight < MAX_VALUE) {
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

void HOP_maintain_SPREAD2_batch(graph_v_of_v<int> &instance_graph, vector<vector<hop_constrained_two_hop_label>> *L, PPR_type *PPR,
								std::vector<hop_constrained_pair_label> &al2, std::vector<hop_constrained_affected_label> *al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k)
{

	for (auto it : al2)
	{
		results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al3, &instance_graph, upper_k]
														  {
				int v = it.first, u = it.second, h_u = it.hop;
				mtx_5992[v].lock_shared();
				std::vector<int> temp = PPR_retrieve(*PPR, v, u);
				mtx_5992[v].unlock_shared();
				PPR_binary_operations_insert(temp, u);

				mtx_599[v].lock();
				auto Lv = (*L)[v]; // to avoid interlocking
				mtx_599[v].unlock();

				for (auto t : temp) {
					if (v < t) {
						long long int d1 = MAX_VALUE;
						int hop_vn = 0;
						for (auto nei : instance_graph[t]) {
                            //mtx_599[nei.first].lock();
							pair<int, int> dis_hop = get_shortest_distance_hop_two_hop_label2((*L)[nei.first], v);
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
                            if(hop_i > upper_k)
								break;
							long long int di = MAX_VALUE;
							for (auto nei : instance_graph[t]) {
								//mtx_599[nei.first].lock();
								di = std::min(di, search_sorted_hop_constrained_two_hop_label((*L)[nei.first], v, hop_i - 1) + (long long int)nei.second);
								//mtx_599[nei.first].unlock();
							}
							if(di >= TwoM_value)
								continue;
							//mtx_599[t].lock_shared();
							//auto query_result = graph_hash_of_mixed_weighted_two_hop_v2_extract_distance_no_reduc2(*L, t.first, v, hop_i);
							auto query_result = hop_constrained_extract_distance_and_hub(*L, t, v, hop_i);
							//mtx_599[t].unlock_shared();

							if (query_result.first > di) { // only add new label when it's absolutely necessary
								mtx_599_1.lock();
								//cout<<"query_result.first > d1 + 1e-5: "<<t_first<<' '<<v << ' ' << hop_vn+1 << ' ' << d1 << endl;
								al3->push_back(hop_constrained_affected_label(t, v, hop_i, di));
								mtx_599_1.unlock();
			
							}
							else {
								if (query_result.second != -1 && query_result.second != v) {
									mtx_5992[t].lock();
									PPR_insert(*PPR, t, query_result.second, v);
									mtx_5992[t].unlock();
								}
								if (query_result.second != -1 && query_result.second != t) {
									mtx_5992[v].lock();
									PPR_insert(*PPR, v, query_result.second, t);
									mtx_5992[v].unlock();
								}
							}
						}
					}
					if (t < v) {
						long long int d1 = MAX_VALUE;
						int hop_vn = 0;
						for (auto nei : instance_graph[v]) {
							//d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], t_first, t.second) + (int)nei.second);
							//mtx_599[nei.first].lock();
							pair<int, int> dis_hop = get_shortest_distance_hop_two_hop_label2((*L)[nei.first], t);
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
                            if(hop_i > upper_k)
								break;
							long long int di = MAX_VALUE;
							for (auto nei : instance_graph[v]) {
								//mtx_599[nei.first].lock();
								di = std::min(di, search_sorted_hop_constrained_two_hop_label((*L)[nei.first], t, hop_i - 1) + (long long int)nei.second);
								//mtx_599[nei.first].unlock();
							}

							if(di >= TwoM_value)
								continue;

							//mtx_599[t].lock_shared();
							//auto query_result = graph_hash_of_mixed_weighted_two_hop_v2_extract_distance_no_reduc2(*L, v, t_first, hop_i);
							auto query_result = hop_constrained_extract_distance_and_hub(*L,v, t, hop_i);
							//mtx_599[t].unlock_shared();

							if (query_result.first > di) {
								mtx_599_1.lock();
								al3->push_back(hop_constrained_affected_label(v, t, hop_i, di));
								mtx_599_1.unlock();
							}
							else {
								if (query_result.second != -1 && query_result.second != v) {
									mtx_5992[t].lock();
									PPR_insert(*PPR, t, query_result.second, v);
									mtx_5992[t].unlock();
								}
								if (query_result.second != -1 && query_result.second != t) {
									mtx_5992[v].lock();
									PPR_insert(*PPR, v, query_result.second, t);
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

void HOP_maintain_SPREAD3_batch(graph_v_of_v<int> &instance_graph, vector<vector<hop_constrained_two_hop_label>> *L, PPR_type *PPR, std::vector<hop_constrained_affected_label> &al3,
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

			mtx_599[v].lock_shared();
			auto Lv = (*L)[v]; // to avoid interlocking
			mtx_599[v].unlock_shared();

			vector<int> dist_hop_changes;
			auto &dist_hop = dist_hop_599_v2[current_tid];
			boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE> pq;
			map<pair<int, int>, pair<hop_constrained_handle_t_for_DIFFUSE, int>> Q_handle;
			vector<int> hubs;
			hubs.resize(instance_graph.size(), -1);
			auto& Q_VALUE = Q_value[current_tid];

			for (auto &it : vec_with_hub_v)
			{
				int u = it.hub_vertex;
				int h_v = it.hop;
				int du = it.distance;

				mtx_599[u].lock_shared();
				auto query_result = hop_constrained_extract_distance_and_hub(*L,u,v, h_v);
				mtx_599[u].unlock_shared();

				bool flag = false;
				if (query_result.first < du)
				{
					if (query_result.second != -1 && query_result.second != v)
					{
						mtx_5992[u].lock();
						PPR_insert(*PPR, u, query_result.second, v);
						mtx_5992[u].unlock();
						flag = true;
					}
					if (query_result.second != -1 && query_result.second != u)
					{
						mtx_5992[v].lock();
						PPR_insert(*PPR, v, query_result.second, u);
						mtx_5992[v].unlock();
						flag = true;
					}
				}

				if (flag == true)
				{
					continue;
				}

				dist_hop[u] = {du, h_v}; //  {dis, hop}
				dist_hop_changes.push_back(u);
				hop_constrained_node_for_DIFFUSE tmp;
				tmp.index = u;
				tmp.hop = h_v;
				tmp.disx = du;
				Q_handle[{u, h_v}] = {pq.push({tmp}), du}; //{hop_constrained_node_for_DIFFUSE,dis}
				Q_VALUE[u][h_v] = du;
			}

			while (!pq.empty())
			{
				int x = pq.top().index;
				int xhv = pq.top().hop;
				int dx = pq.top().disx;
				pq.pop();
				if(xhv <= upper_k)
                    Q_VALUE[x][xhv] = MAX_VALUE;

				mtx_599[x].lock();
				int d_old = search_sorted_hop_constrained_two_hop_label((*L)[x], v, xhv);
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
					hop_constrained_node_for_DIFFUSE node = {xnei, xhv + 1, (weightTYPE)d_new};

					if (v < xnei)
					{

						if (dist_hop[xnei].first == -1)
						{
							Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                            Q_VALUE[xnei][hop_nei] = d_new;
                            dist_hop[xnei].first = d_new;
							dist_hop[xnei].second = hop_nei;
							dist_hop_changes.push_back(xnei);

                            mtx_599[xnei].lock_shared();
							std::pair<int, int> tmp = hop_constrained_extract_distance_and_hub(*L,xnei,v, xhv + 1); 
							mtx_599[xnei].unlock_shared();
							hubs[xnei] = tmp.second;
						}
						if (d_new < dist_hop[xnei].first)
						{
							//if (Q_handle.find({xnei, node.hop}) != Q_handle.end())
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
							//if (Q_handle.find({xnei, node.hop}) != Q_handle.end())
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

						if (dist_hop[xnei].first < d_new)
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

			for (int i : dist_hop_changes) {
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

void HOP_WeightIncreaseMaintenance_improv_batch(graph_v_of_v<int> &instance_graph, hop_constrained_case_info &mm, vector<pair<int, int>> &v, vector<int> &w_old_vec,
												ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int t)
{

	global_query_times = 0;
	label_operation_times = 0;

	std::vector<hop_constrained_affected_label> al1, al3;
	std::vector<hop_constrained_pair_label> al2;

	std::map<pair<int, int>, int> w_old_map;
	int batch_size = v.size();
	for (int i = 0; i < batch_size; i++)
	{
		if (v[i].first < v[i].second)
		{
			swap(v[i].first, v[i].second);
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
				for (auto it : mm.L[v1]) {
					int search_weight = search_sorted_hop_constrained_two_hop_label(mm.L[v2], it.hub_vertex, it.hop + 1);
					if (it.hub_vertex <= v2 && search_weight >= (long long int)it.distance + w_old && search_weight < MAX_VALUE) {
						mtx_599_1.lock();
						al1.push_back(hop_constrained_affected_label(v2, it.hub_vertex, it.hop + 1, it.distance + w_old));
						mtx_599_1.unlock();
					}
				}
				for (auto it : mm.L[v2]) {
				    int search_weight = search_sorted_hop_constrained_two_hop_label(mm.L[v1], it.hub_vertex, it.hop + 1);
					if (it.hub_vertex <= v1 && search_weight >= (long long int)it.distance + w_old && search_weight < MAX_VALUE) {
						mtx_599_1.lock();
						al1.push_back(hop_constrained_affected_label(v1, it.hub_vertex, it.hop + 1, it.distance + w_old));
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
    mm.mark_time("start increase");
	HOP_maintain_SPREAD1_batch(instance_graph, &mm.L, al1, &al2, w_old_map, pool_dynamic, results_dynamic, t);
	HOP_maintain_SPREAD2_batch(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic, mm.upper_k);
	HOP_maintain_SPREAD3_batch(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic, mm.upper_k, t);
	mm.mark_time("end increase");
}
