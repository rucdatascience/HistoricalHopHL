#pragma once

#include <CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_maintain_PLL.h>
using namespace nonHop;
namespace nonHop
{
	void PI11(graph_v_of_v<int>& instance_graph, vector<vector<two_hop_label>>* L,
		std::vector<affected_label>& al1_curr, std::vector<affected_label>* al1_next,
		std::map<pair<int, int>, weightTYPE >& w_old_map,
		ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic, int time) {

		for (auto it : al1_curr) {
			results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, al1_next, &instance_graph, &w_old_map, time] {
				for (auto nei : instance_graph[it.first]) {
					mtx_595[nei.first].lock();
					two_hop_label search_weight = search_sorted_two_hop_label_entity_not_realTime((*L)[nei.first], it.second, time);
					mtx_595[nei.first].unlock();
					weightTYPE w_old;
					int i = 0;
					if (w_old_map.count(pair<int, int>(it.first, nei.first)) > 0) {
						i = 1;
						w_old = w_old_map[pair<int, int>(it.first, nei.first)];
					}
					else if (w_old_map.count(pair<int, int>(nei.first, it.first)) > 0) {
						i = 2;
						w_old = w_old_map[pair<int, int>(nei.first, it.first)];
					}
					else {
						i = 3;
						w_old = nei.second;
					}
					if (it.dis + w_old == search_weight.distance) {
						mtx_595_1.lock();
						/*std::cout << i << " and its t_s is " << search_weight.t_s << std::endl;
						std::cout << nei.first << " to " << it.second << " weight is it.dis " << it.dis << " + w_old " << w_old << "=" << it.dis + w_old << " and search_weight is " << search_weight.distance << std::endl;*/
						al1_next->push_back(affected_label(nei.first, it.second, search_weight.distance));
						mtx_595_1.unlock();
					}
				}
				mtx_595[it.first].lock();
				insert_sorted_two_hop_label((*L)[it.first], it.second, MAX_VALUE, time); // this does not change the size of L[it->first] here, so does not need to lock here
				mtx_595[it.first].unlock();
				return 1; }));
		}

		for (auto&& result : results_dynamic) {
			result.get();
		}
		std::vector<std::future<int>>().swap(results_dynamic);
	}

	void PI12(graph_v_of_v<int>& instance_graph, vector<vector<two_hop_label>>* L, PPR_type* PPR,
		std::vector<affected_label>& al1_curr, std::vector<pair_label>* al2_next, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic, int time) {
		std::cout << "al1_curr length is " << al1_curr.size() << std::endl;
		for (auto it : al1_curr) {
			results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al2_next, &instance_graph, time] {

				int v = it.first, u = it.second;
				mtx_5952[v].lock();
				std::vector<int> temp = PPR_retrieve(*PPR, v, u);
				mtx_5952[v].unlock();
				PPR_binary_operations_insert(temp, u);

				mtx_595[v].lock();
				auto Lv = (*L)[v]; // to avoid interlocking
				mtx_595[v].unlock();

				for (auto t : temp) {
					if (v < t) {
						long long d1 = MAX_VALUE;
						for (auto nei : instance_graph[t]) {
							mtx_595[nei.first].lock();
							d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], v) + (long long)nei.second);
							mtx_595[nei.first].unlock();
						}
						mtx_595[t].lock();
						auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[t], Lv);
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
					if (t < v) {
						long long d1 = MAX_VALUE;
						for (auto nei : instance_graph[v]) {
							mtx_595[nei.first].lock();
							d1 = min(d1, search_sorted_two_hop_label((*L)[nei.first], t) + (long long)nei.second);
							mtx_595[nei.first].unlock();
						}
						mtx_595[t].lock();
						auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[t], Lv);
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

	void PI22(graph_v_of_v<int>& instance_graph, vector<vector<two_hop_label>>* L, PPR_type* PPR,
		std::vector<pair_label>& al2_curr, std::vector<pair_label>* al2_next, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic, int time) {

		for (auto it = al2_curr.begin(); it != al2_curr.end(); it++) {
			results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, al2_next, &instance_graph, time] {
				try {
					//if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
					//	throw reach_limit_time_string;
					//}

					mtx_595[it->second].lock();
					auto Lxx = (*L)[it->second]; // to avoid interlocking
					mtx_595[it->second].unlock();

					for (auto nei : instance_graph[it->first]) {
						if (nei.first > it->second) {
							mtx_595[it->first].lock();
							long long search_result = search_sorted_two_hop_label((*L)[it->first], it->second) + (long long)nei.second;
							mtx_595[it->first].unlock();
							mtx_595[nei.first].lock();
							auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4((*L)[nei.first], Lxx);
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
									PPR_insert(*PPR, nei.first, query_result.second, it->second);
									mtx_5952[nei.first].unlock();
								}
								if (query_result.second != nei.first) {
									mtx_5952[it->second].lock();
									PPR_insert(*PPR, it->second, query_result.second, nei.first);
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

		for (auto&& result : results_dynamic) {
			result.get();
		}
		std::vector<std::future<int>>().swap(results_dynamic);
	}

	void nonHOP_WeightIncrease2021_batch(graph_v_of_v<int>& instance_graph, two_hop_case_info& mm, vector<pair<int, int> >& v, vector<weightTYPE>& w_old_vec,
		ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic, int time) {

		label_operation_times = 0;
		global_query_times = 0;

		begin_time = std::chrono::high_resolution_clock::now();
		max_run_time_nanosec = _2023algo_max_second * 1e9;

		std::map<pair<int, int>, weightTYPE> w_old_map;
		int batch_size = v.size();
		for (int i = 0; i < batch_size; i++) {
			if (v[i].first > v[i].second) {
				swap(v[i].first, v[i].second);
			}
			if (w_old_map.count(v[i]) == 0) {
				w_old_map[v[i]] = w_old_vec[i];
			}
		}

		std::vector<affected_label> al1_curr, al1_next;
		std::vector<pair_label> al2_curr, al2_next;

		for (auto& iter : w_old_map) {
			int v1 = iter.first.first;
			int v2 = iter.first.second;
			weightTYPE w_old = iter.second;

			for (auto it : mm.L[v1]) {
				mtx_595[v2].lock();
				long long search_weight = search_sorted_two_hop_label(mm.L[v2], it.vertex);
				mtx_595[v2].unlock();
				if (it.vertex <= v2 && search_weight >= (long long)it.distance + w_old && search_weight < MAX_VALUE) {
					al1_curr.push_back(affected_label(v2, it.vertex, it.distance + w_old));
				}
			}
			for (auto it : mm.L[v2]) {
				mtx_595[v1].lock();
				long long search_weight = search_sorted_two_hop_label(mm.L[v1], it.vertex);
				mtx_595[v1].unlock();
				if (it.vertex <= v1 && search_weight >= (long long)it.distance + w_old && search_weight < MAX_VALUE) {
					al1_curr.push_back(affected_label(v1, it.vertex, it.distance + w_old));
				}
			}
		}



		while (al1_curr.size() || al2_curr.size()) {

			//if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
			//	throw reach_limit_time_string;
			//}
			std::cout << "2021 increase pl11" << std::endl;
			PI11(instance_graph, &mm.L, al1_curr, &al1_next, w_old_map, pool_dynamic, results_dynamic, time);
			std::cout << "2021 increase pl12" << std::endl;
			PI12(instance_graph, &mm.L, &mm.PPR, al1_curr, &al2_next, pool_dynamic, results_dynamic, time);
			std::cout << "2021 increase pl22" << std::endl;
			PI22(instance_graph, &mm.L, &mm.PPR, al2_curr, &al2_next, pool_dynamic, results_dynamic, time);
			al1_curr = al1_next;
			al2_curr = al2_next;
			std::vector<affected_label>().swap(al1_next);
			std::vector<pair_label>().swap(al2_next);
		}
	}
}

