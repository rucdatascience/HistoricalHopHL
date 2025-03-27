#pragma once
#include "Historical/graph_with_time_span/graph.h"
#include <filesystem>
#include "Historical/experiment/experiment_config.h"
#include "Historical/graph_with_time_span/two_hop_label.h"
#include <boost/heap/fibonacci_heap.hpp>
#include <shared_mutex>
#include <CPU/tool_functions/ThreadPool.h>
#include "Historical/utils/BinaryPersistence.h"
#include <fstream>

namespace experiment
{
	template <typename weight_type>
	void read_graph(graph<weight_type> &graph, ExperimentConfig &config)
	{
		std::string readPath = config.data_source.string();
		std::string line_content;
		boost::random::uniform_int_distribution<> random_weight = boost::random::uniform_int_distribution<>(1, 100);
		int v_num = 0;
		std::ifstream myfile(readPath);
		if (myfile.is_open())
		{
			while (getline(myfile, line_content))
			{
				std::vector<std::string> Parsed_content = experiment::parse_string(line_content, "\t");

				if (!Parsed_content[0].compare("#"))
				{
					if (!Parsed_content[1].compare("Nodes"))
					{
						v_num = std::stoi(Parsed_content[2]);
						graph.resize(v_num);
					}
				}
				else
				{
					int v1 = std::stoi(Parsed_content[0]);
					int v2 = std::stoi(Parsed_content[1]);
					int w = random_weight(boost_random_time_seed);
					graph.add_edge(v1, v2, w);
				}
			}
		}
	};

	template <typename weight_type>
	class iteration_info
	{
	private:
		int _v_num;
		int _iteration;
		int _change_num;
		int _upper;
		int _lower;
		boost::random::uniform_int_distribution<> _random_v;
		boost::random::uniform_int_distribution<> _random_weight;
		graph<weight_type> instance_graph;

	public:
		std::vector<std::queue<change_edge_info>> q_list;
		iteration_info(){

		};
		iteration_info(int v_num, int iteration, int change_num, int upper, int lower, graph<weight_type> graph) : _v_num(v_num), _iteration(iteration), _change_num(change_num), _upper(upper), _lower(lower), instance_graph(graph)
		{
			this->_random_v = boost::random::uniform_int_distribution<>(0, this->_v_num);
			this->_random_weight = boost::random::uniform_int_distribution<>(this->_lower, this->_upper);
			q_list = std::vector<std::queue<change_edge_info>>(this->_iteration + 1, std::queue<change_edge_info>());
		};

		void build_random_change()
		{
			std::map<std::pair<int, int>, int> pair2dis;
			for (int i = 1; i <= this->_iteration; i++)
			{
				int j = 0;
				while (j < this->_change_num)
				{
					int index_i = this->_random_v(boost_random_time_seed);
					if (instance_graph[index_i].size() == 0)
					{
						continue;
					}
					boost::random::uniform_int_distribution<> dis_inner(0, instance_graph[index_i].size() - 1);
					int index_j_relatively = dis_inner(boost_random_time_seed);
					int i_j_weight = this->_random_weight(boost_random_time_seed);
					int index_j = instance_graph[index_i][index_j_relatively].first;
					if (index_i > index_j)
					{
						std::swap(index_i, index_j);
					}
					change_edge_info info = {index_i, index_j, i_j_weight, i};
					q_list[i].push(info);
					++j;
				}
			}
		}

		void serialize(std::ofstream &out) const
		{
			saveBinary(out, _v_num);
			saveBinary(out, _iteration);
			saveBinary(out, _change_num);
			saveBinary(out, _upper);
			saveBinary(out, _lower);
			saveBinary(out, instance_graph);

			for (const auto &queue : q_list)
			{
				size_t size = queue.size();
				saveBinary(out, size);
				std::queue<change_edge_info> temp = queue;
				while (!temp.empty())
				{
					saveBinary(out, temp.front().v1);
					saveBinary(out, temp.front().v2);
					saveBinary(out, temp.front().time);
					saveBinary(out, temp.front().weight);
					temp.pop();
				}
			}
		}
		void deserialize(std::ifstream &in)
		{
			loadBinary(in, _v_num);
			loadBinary(in, _iteration);
			loadBinary(in, _change_num);
			loadBinary(in, _upper);
			loadBinary(in, _lower);
			loadBinary(in, instance_graph);

			q_list.resize(this->_iteration + 1, std::queue<change_edge_info>());
			for (auto &queue : q_list)
			{
				size_t size;
				loadBinary(in, size);
				for (size_t i = 0; i < size; i++)
				{
					change_edge_info info;
					loadBinary(in, info.v1);
					loadBinary(in, info.v2);
					loadBinary(in, info.time);
					loadBinary(in, info.weight);
					queue.push(info);
				}
			}
		}
	};

	template <typename weight_type>
	class BinarySerializer<iteration_info<weight_type>>
	{
	public:
		static void saveBinary(std::ofstream &out, const iteration_info<weight_type> &info)
		{
			info.serialize(out);
		}

		static void loadBinary(std::ifstream &in, iteration_info<weight_type> &info)
		{
			info.deserialize(in);
		}
	};

	namespace nonhop
	{
		int globalLabelSize = 0;
		int globalLabelCleanSize = 0;
		int globalPprSize = 0;
		int max_N_ID_for_mtx_595 = 1e7;
		std::vector<std::shared_mutex> mtx_595(max_N_ID_for_mtx_595);
		std::vector<std::shared_mutex> ppr_595(max_N_ID_for_mtx_595);
		std::vector<std::vector<two_hop_label>> L_temp_595;
		PPR_TYPE::PPR_type PPR_595;
		std::vector<std::vector<int>> P_dij_595;
		std::vector<std::vector<int>> T_dij_595;
		typedef typename boost::heap::fibonacci_heap<two_hop_label>::handle_type PLL_handle_t_for_sp;
		std::vector<std::vector<PLL_handle_t_for_sp>> Q_handles_595;
		std::queue<int> Qid_595;
		std::vector<std::vector<two_hop_label>> Lv_final;
		void PLL_dij_function(int v_k, graph<int> &input_graph)
		{
			mtx_595[max_N_ID_for_mtx_595 - 1].lock();
			auto startTime = std::chrono::steady_clock::now();
			int used_id = Qid_595.front();
			Qid_595.pop();
			mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

			std::vector<int> P_changed_vertices, T_changed_vertices;
			std::vector<int> &T_dij = T_dij_595[used_id], P_dij = P_dij_595[used_id];

			std::vector<PLL_handle_t_for_sp> &Q_handles = Q_handles_595[used_id];

			boost::heap::fibonacci_heap<two_hop_label> Q;
			two_hop_label node(0);
			node.vertex = v_k;
			node.distance = 0;
			Q_handles[v_k] = Q.push(node);
			P_dij[v_k] = 0;
			P_changed_vertices.push_back(v_k);
			mtx_595[v_k].lock_shared();
			for (const auto &xx : L_temp_595[v_k])
			{
				int L_v_k_i_vertex = xx.vertex;
				T_dij[L_v_k_i_vertex] = xx.distance;
				T_changed_vertices.push_back(L_v_k_i_vertex);
			}
			mtx_595[v_k].unlock_shared();
			int new_label_num = 0;
			int count = 0;
			while (Q.size())
			{
				count++;
				node = Q.top();
				Q.pop();
				int u = node.vertex;

				int P_u = node.distance;
				
				if(v_k > u){
					continue;
				}

				int query_v_k_u = std::numeric_limits<int>::max();
				int common_hub_for_query_v_k_u = 0;
				mtx_595[u].lock_shared(); // put lock in for loop is very slow
				for (const auto &xx : L_temp_595[u])
				{
					long long int dis = xx.distance + (long long int)T_dij[xx.vertex]; // long long int is to avoid overflow
					if (query_v_k_u > dis)
					{
						query_v_k_u = dis;
						common_hub_for_query_v_k_u = xx.vertex;
					}
				}
				mtx_595[u].unlock_shared();

				if (P_u < query_v_k_u)
				{
					node.vertex = v_k;
					node.distance = P_u;

					mtx_595[u].lock();
					L_temp_595[u].push_back(node);
					++globalLabelSize;
					mtx_595[u].unlock();
					new_label_num++;

					for (const auto &xx : input_graph.ADJs[u])
					{
						int adj_v = xx.first, ec = xx.second;
						if (P_dij[adj_v] == std::numeric_limits<int>::max())
						{
							node.vertex = adj_v;
							node.distance = P_u + ec;

							Q_handles[adj_v] = Q.push(node);
							P_dij[adj_v] = node.distance;
							P_changed_vertices.push_back(adj_v);
						}
						else if (P_dij[adj_v] > P_u + ec)
						{
							node.vertex = adj_v;
							node.distance = P_u + ec;
							Q.update(Q_handles[adj_v], node);
							P_dij[adj_v] = node.distance;
						}
					}
				}
				else
				{
					if (common_hub_for_query_v_k_u != v_k)
					{
						ppr_595[u].lock();
						PPR_TYPE::PPR_insert(PPR_595, u, common_hub_for_query_v_k_u, v_k);
						ppr_595[u].unlock();
					}
					if (common_hub_for_query_v_k_u != u)
					{
						ppr_595[v_k].lock();
						PPR_TYPE::PPR_insert(PPR_595, v_k, common_hub_for_query_v_k_u, u);
						ppr_595[v_k].unlock();
					}
				}
			}

			for (auto i : P_changed_vertices)
			{
				P_dij[i] = std::numeric_limits<int>::max(); // reverse-allocate P values
			}
			for (auto i : T_changed_vertices)
			{
				T_dij[i] = std::numeric_limits<int>::max(); // reverse-allocate T values
			}

			mtx_595[max_N_ID_for_mtx_595 - 1].lock();
			auto endTime = std::chrono::steady_clock::now();
			std::cout << "print pll v_k: " << v_k << " time cost is " << std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count() << std::endl;
			Qid_595.push(used_id);
			mtx_595[max_N_ID_for_mtx_595 - 1].unlock();
		}

		std::vector<std::vector<two_hop_label>> sortL(int num_of_threads)
		{

			/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/

			int N = L_temp_595.size();
			std::vector<std::vector<two_hop_label>> output_L(N);

			/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/
			ThreadPool pool(num_of_threads);
			std::vector<std::future<int>> results; // return typename: xxx
			for (int v_k = 0; v_k < N; v_k++)
			{
				results.emplace_back(
					pool.enqueue([&output_L, v_k] { // pass const type value j to thread; [] can be empty
						sort(L_temp_595[v_k].begin(), L_temp_595[v_k].end(), compare_two_hop_label_small_to_large);
						std::vector<two_hop_label>(L_temp_595[v_k]).swap(L_temp_595[v_k]); // swap�ͷ�vector�ж���ռ�
						output_L[v_k] = L_temp_595[v_k];
						std::vector<two_hop_label>().swap(L_temp_595[v_k]); // clear new labels for RAM efficiency

						return 1; // return to results; the return type must be the same with results
					}));
			}
			for (auto &&result : results)
				result.get(); // all threads finish here
			results.clear();

			return output_L;
		}

		void clean_L(two_hop_case_info &case_info, int thread_num)
		{

			auto &L = case_info.L;
			int N = L.size();

			ThreadPool pool(thread_num);
			std::vector<std::future<int>> results;

			for (int v = 0; v < N; v++)
			{
				results.emplace_back(
					pool.enqueue([v, &L] { // pass const type value j to thread; [] can be empty
						mtx_595[max_N_ID_for_mtx_595 - 1].lock();
						int used_id = Qid_595.front();
						Qid_595.pop();
						mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

						std::vector<two_hop_label> &Lv_final_inner = Lv_final[v];

						std::vector<two_hop_label> &Lv = L[v];

						auto &T = T_dij_595[used_id];

						for (const auto &Lvi : Lv)
						{
							int u = Lvi.vertex;
							if (v == u)
							{
								Lv_final_inner.push_back(two_hop_label(Lvi));
								++globalLabelCleanSize;
								T[v] = Lvi.distance;
								continue;
							}
							const auto &Lu = L[u];

							int min_dis = std::numeric_limits<int>::max();
							for (const auto &label : Lu)
							{
								long long int query_dis = label.distance + (long long int)T[label.vertex];
								if (query_dis < min_dis)
								{
									min_dis = query_dis;
								}
							}

							if (min_dis > Lvi.distance)
							{
								Lv_final_inner.push_back(two_hop_label(Lvi));
								++globalLabelCleanSize;
								T[u] = Lvi.distance;
							}
						}

						for (const auto &label : Lv_final_inner)
						{
							T[label.vertex] = std::numeric_limits<int>::max();
						}

						mtx_595[max_N_ID_for_mtx_595 - 1].lock();
						Qid_595.push(used_id);
						std::cout << "print pll v: " << v  << std::endl;
						mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

						return 1; // return to results; the return type must be the same with results
					}));
			}

			for (auto &&result : results)
				result.get(); // all threads finish here
			case_info.L = std::move(Lv_final);
			results.clear();
		}

		void PLL_clear_global_values()
		{
			std::vector<std::vector<two_hop_label>>().swap(L_temp_595);
			std::vector<std::vector<two_hop_label>>().swap(Lv_final);
			PPR_TYPE::PPR_type().swap(PPR_595);
			std::queue<int>().swap(Qid_595);
			std::vector<std::vector<int>>().swap(P_dij_595);
			std::vector<std::vector<int>>().swap(T_dij_595);
			std::vector<std::vector<PLL_handle_t_for_sp>>().swap(Q_handles_595);
		}

		template <typename weight_type>
		void pll(graph<weight_type> &input_graph, nonhop::two_hop_case_info &case_info)
		{
			//----------------------------------- step 1: initialization ------------------------------------------------------------------
			timer.startSubtask("step 1: initialization");
			int num_of_threads = case_info.thread_num;
			int N = input_graph.ADJs.size();

			L_temp_595.resize(N);
			PPR_595.resize(N);
			Lv_final.resize(N);
			timer.endSubtask();
			//---------------------------------------------------------------------------------------------------------------------------------------

			//----------------------------------------------- step 2: generate labels ---------------------------------------------------------------
			timer.startSubtask("step 2: generate labels");
			{
				// to save RAM of ThreadPool
				/*seaching shortest paths*/
				ThreadPool pool(num_of_threads);
				std::vector<std::future<int>> results; // return typename: xxx
				P_dij_595.resize(num_of_threads);
				T_dij_595.resize(num_of_threads);
				Q_handles_595.resize(num_of_threads);
				std::queue<int>().swap(Qid_595);
				for (int i = 0; i < num_of_threads; i++)
				{
					P_dij_595[i].resize(N, std::numeric_limits<int>::max());
					T_dij_595[i].resize(N, std::numeric_limits<int>::max());
					Q_handles_595[i].resize(N);
					Qid_595.push(i);
				}

				int last_check_vID = N - 1;

				for (int v_k = 0; v_k <= last_check_vID; v_k++)
				{
					results.emplace_back(
						pool.enqueue([v_k, &input_graph, last_check_vID] { // pass const type value j to thread; [] can be empty
							PLL_dij_function(v_k, input_graph);
							return 1; // return to results; the return type must be the same with results
						}));
				}
				for (auto &&result : results)
					result.get(); // all threads finish here
				results.clear();
			}
			timer.endSubtask();
			//---------------------------------------------------------------------------------------------------------------------------------------

			//----------------------------------------------- step 3: sortL ---------------------------------------------------------------
			timer.startSubtask("step 3: sortL");
			case_info.L = sortL(num_of_threads);
			case_info.PPR = PPR_595;
			timer.endSubtask();
			//---------------------------------------------------------------------------------------------------------------------------------------

			//----------------------------------------------- step 3: canonical_repair ---------------------------------------------------------------
			timer.startSubtask("step 4: canonical_repair");
			clean_L(case_info, num_of_threads);
			timer.endSubtask();
			//---------------------------------------------------------------------------------------------------------------------------------------
			PLL_clear_global_values();
		}

	}

	namespace hop
	{
		int globalLabelSize = 0;
		int globalLabelCleanSize = 0;
		int globalPprSize = 0;
		int max_N_ID_for_mtx_599 = 1e7;
		std::queue<int> Qid_599;
		std::vector<std::shared_mutex> mtx_599(max_N_ID_for_mtx_599);
		std::vector<std::shared_mutex> ppr_599(max_N_ID_for_mtx_599);

		int global_upper_k = 0;

		template <typename weight_type>
		graph<weight_type> ideal_graph_599;

		typedef typename boost::heap::fibonacci_heap<two_hop_label>::handle_type hop_constrained_node_handle;
		std::vector<std::vector<two_hop_label>> L_temp_599;
		std::vector<std::vector<two_hop_label>> Lv_final_599;
		PPR_TYPE::PPR_type PPR_599;
		std::vector<std::vector<std::vector<std::pair<int, int>>>> Temp_L_vk_599;
		std::vector<std::vector<std::pair<int, int>>> dist_hop_599;
		std::vector<std::vector<std::vector<std::pair<hop_constrained_node_handle, int>>>> Q_handle_priorities_599;
		std::vector<std::vector<std::vector<int>>> Vh_599;

		/* use asynchronous tasks to generate L labels and the parameter v_k should be the index of the endpoint */
		template <typename weight_type>
		void HSDL_thread_function(int v_k)
		{
			/* get unique thread id */
			/* critical section obtain array index  */
			mtx_599[max_N_ID_for_mtx_599 - 1].lock();
			auto startTime = std::chrono::steady_clock::now();
			int used_id = Qid_599.front();
			Qid_599.pop();
			mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

			/* store the Temp L and dist_hop in current thread*/
			std::vector<int> Temp_L_vk_changes, dist_hop_changes;
			/* Temp_L_vk stores the dest_vertex_id and distance and hop */
			auto &Temp_L_vk = Temp_L_vk_599[used_id];
			auto &dist_hop = dist_hop_599[used_id]; // record the minimum distance (and the corresponding hop) of a searched vertex in Q
			std::vector<std::pair<int, int>> Q_handle_priorities_changes;
			/* get the label list in current thread*/
			auto &Q_handle_priorities = Q_handle_priorities_599[used_id];

			/* a class contains information about destination vertex, hop count, and cost in the priority queue */
			boost::heap::fibonacci_heap<two_hop_label> Q;

			/* generate a label for the starting vertex itself */
			two_hop_label node;
			node.hub_vertex = v_k;
			node.hop = 0;
			node.distance = 0;
			Q_handle_priorities[v_k][0] = {Q.push({node}), node.distance};
			Q_handle_priorities_changes.push_back({v_k, 0});
			double costQuery = 0;
			double costUpdate = 0;
			double costPPR = 0;
			// size_t size = 0;
			/* Temp_L_vk_599 stores the label (dist and hop) of vertex v_k */
			mtx_599[v_k].lock_shared();
			/* root is vk-> vk->obj info -> vector<obj> -> index-> vertexId obj-><distance,hop> */
			for (auto &xx : L_temp_599[v_k])
			{
				int L_vk_vertex = xx.hub_vertex;
				Temp_L_vk[L_vk_vertex].push_back({xx.distance, xx.hop});
				Temp_L_vk_changes.push_back(L_vk_vertex);
			}
			mtx_599[v_k].unlock_shared();
			/*  dist_hop_599 stores the shortest distance from vk to any other vertices with its hop_cst,
				note that the hop_cst is determined by the shortest distance */
			dist_hop[v_k] = {0, 0};
			dist_hop_changes.push_back(v_k);
			while (Q.size() > 0)
			{
				// size = std::max(size, Q.size());
				/* poll the vertex from heap.In other words, poll the vertex with the minimal cost */
				node = Q.top();
				Q.pop();

				/* current node, u, which is the node generating the labels*/
				int u = node.hub_vertex;

				if (v_k > u)
				{
					continue;
				}

				int u_hop = node.hop;
				int P_u = node.distance;
				int common_hub_for_query_v_k_u = -1;
				int query_v_k_u = std::numeric_limits<int>::max();
				auto startTime1 = std::chrono::steady_clock::now();
				mtx_599[u].lock_shared();
				for (auto &xx : L_temp_599[u])
				{
					int common_v = xx.hub_vertex;
					for (auto &yy : Temp_L_vk[common_v])
					{
						long long int dis_opt = (long long int)xx.distance + yy.first;
						if (xx.hop + yy.second <= u_hop)
						{
							long long int dis = (long long int)xx.distance + yy.first;
							if (query_v_k_u > dis)
							{
								query_v_k_u = dis;
								common_hub_for_query_v_k_u = xx.hub_vertex;
							}
						}
					}
				}
				mtx_599[u].unlock_shared();
				auto endTime1 = std::chrono::steady_clock::now();
				costQuery += std::chrono::duration_cast<std::chrono::duration<double>>(endTime1 - startTime1).count();
				if (P_u < query_v_k_u)
				{
					node.hub_vertex = v_k;
					node.hop = u_hop;
					node.distance = P_u;
					mtx_599[u].lock();
					L_temp_599[u].push_back(node);
					++globalLabelSize;
					mtx_599[u].unlock();

					if (u_hop + 1 > global_upper_k)
					{
						continue;
					}

					/* update adj */
					/* Traverse neighboring nodes */
					for (auto &xx : ideal_graph_599<weight_type>[u])
					{
						/* adh_v is the neighborhood and the ec is the distance from u to ajd_v*/
						int adj_v = xx.first, ec = xx.second;

						/* update node info */
						node.hub_vertex = adj_v;
						node.distance = P_u + ec;
						node.hop = u_hop + 1;

						auto &yy = Q_handle_priorities[adj_v][node.hop];

						if (yy.second <= node.distance)
						{ // adj_v has been reached with a smaller distance and the same hop
							continue;
						}
						/* the vertex has not been visited*/
						if (dist_hop[adj_v].first == std::numeric_limits<int>::max())
						{ // adj_v has not been reached
							auto updateStartTime = std::chrono::steady_clock::now();
							yy = {Q.push({node}), node.distance};
							Q_handle_priorities_changes.push_back({adj_v, node.hop});
							dist_hop[adj_v].first = node.distance;
							dist_hop[adj_v].second = node.hop;
							dist_hop_changes.push_back(adj_v);
							auto updateEndTime = std::chrono::steady_clock::now();
							costUpdate += std::chrono::duration_cast<std::chrono::duration<double>>(updateEndTime - updateStartTime).count();
						}
						else
						{
							if (node.distance < dist_hop[adj_v].first)
							{ // adj_v has been reached with a less distance
								auto updateStartTime = std::chrono::steady_clock::now();
								if (yy.second != std::numeric_limits<int>::max())
								{
									Q.update(yy.first, node);
									yy.second = node.distance;
								}
								else
								{
									yy = {Q.push(node), node.distance};
									Q_handle_priorities_changes.push_back({adj_v, node.hop});
								}
								dist_hop[adj_v].first = node.distance;
								dist_hop[adj_v].second = node.hop;
								auto updateEndTime = std::chrono::steady_clock::now();
								costUpdate += std::chrono::duration_cast<std::chrono::duration<double>>(updateEndTime - updateStartTime).count();
							}
							else if (node.hop < dist_hop[adj_v].second)
							{ // adj_v has been reached with a less hop
								auto updateStartTime = std::chrono::steady_clock::now();
								if (yy.second != std::numeric_limits<int>::max())
								{
									Q.update(yy.first, node);
									yy.second = node.distance;
								}
								else
								{
									yy = {Q.push(node), node.distance};
									Q_handle_priorities_changes.push_back({adj_v, node.hop});
								}
								auto updateEndTime = std::chrono::steady_clock::now();
								costUpdate += std::chrono::duration_cast<std::chrono::duration<double>>(updateEndTime - updateStartTime).count();
							}
						}
					}
				}
				else
				{
					auto pprStartTime = std::chrono::steady_clock::now();
					/* add v_k into PPR(u,common_hub_for_query_v_k_u), and add u into PPR(v_k,common_hub_for_query_v_k_u)*/
					if (common_hub_for_query_v_k_u != v_k)
					{
						ppr_599[u].lock();
						PPR_TYPE::PPR_insert(PPR_599, u, common_hub_for_query_v_k_u, v_k);
						++globalPprSize;
						ppr_599[u].unlock();
					}
					if (common_hub_for_query_v_k_u != u)
					{
						ppr_599[v_k].lock();
						PPR_TYPE::PPR_insert(PPR_599, v_k, common_hub_for_query_v_k_u, u);
						++globalPprSize;
						ppr_599[v_k].unlock();
					}
					auto pprEndTime = std::chrono::steady_clock::now();
					costPPR += std::chrono::duration_cast<std::chrono::duration<double>>(pprEndTime - pprStartTime).count();
				}
			}
			auto restoreStartTime = std::chrono::steady_clock::now();
			for (auto &xx : Temp_L_vk_changes)
			{
				std::vector<std::pair<int, int>>().swap(Temp_L_vk[xx]);
			}
			for (auto &xx : dist_hop_changes)
			{
				dist_hop[xx] = {std::numeric_limits<int>::max(), 0};
			}
			hop_constrained_node_handle handle_x;
			for (auto &xx : Q_handle_priorities_changes)
			{
				Q_handle_priorities[xx.first][xx.second] = {handle_x, std::numeric_limits<int>::max()};
			}
			auto restoreEndTime = std::chrono::steady_clock::now();

			// mtx_599[v_k].lock();
			// std::vector<two_hop_label>(L_temp_599[v_k]).swap(L_temp_599[v_k]);
			// mtx_599[v_k].unlock();

			mtx_599[max_N_ID_for_mtx_599 - 1].lock();
			auto endTime = std::chrono::steady_clock::now();
			double cost = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count();
			std::cout << "print pll v_k: " << v_k << " time cost is " << cost << std::endl;
			if (cost > 50)
			{
				std::cout << " query opt by L cost is " << costQuery << std::endl;
				std::cout << " update cost is " << costUpdate << std::endl;
				std::cout << " ppr update cost is " << costPPR << std::endl;
				std::cout << " restore update cost is " << std::chrono::duration_cast<std::chrono::duration<double>>(restoreEndTime - restoreStartTime).count() << std::endl;
				// std::cout << " Q size is " << size << std::endl;
			}
			Qid_599.push(used_id);
			mtx_599[max_N_ID_for_mtx_599 - 1].unlock();
		}

		std::vector<std::vector<two_hop_label>> hop_constrained_sortL(int num_of_threads)
		{

			/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/

			int N = L_temp_599.size();
			std::vector<std::vector<two_hop_label>> output_L(N);

			/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/
			ThreadPool pool(num_of_threads);
			std::vector<std::future<int>> results; // return typename: xxx
			for (int v_k = 0; v_k < N; v_k++)
			{
				results.emplace_back(
					pool.enqueue([&output_L, v_k] { // pass const type value j to thread; [] can be empty
						sort(L_temp_599[v_k].begin(), L_temp_599[v_k].end(), compare_hop_constrained_two_hop_label);
						std::vector<two_hop_label>(L_temp_599[v_k]).swap(L_temp_599[v_k]); // ʹ��vector��swap�Ż��ڴ�ռ�ã��ͷŶ���Ŀռ�
						output_L[v_k] = L_temp_599[v_k];
						std::vector<two_hop_label>().swap(L_temp_599[v_k]); // clear new labels for RAM efficiency

						return 1; // return to results; the return type must be the same with results
					}));
			}
			for (auto &&result : results)
				result.get(); // all threads finish here

			return output_L;
		}

		/*canonical_repair*/
		void hop_constrained_clean_L(two_hop_case_info &case_info, int thread_num)
		{

			auto &L = case_info.L;
			int N = L.size();

			ThreadPool pool(thread_num);
			std::vector<std::future<int>> results;
			/* test the correctness of async */
			// vector<int> list;
			// list.push_back(4);
			// list.push_back(5);
			// list.push_back(0);
			// list.push_back(1);
			// list.push_back(2);
			// list.push_back(3);
			// for (int v = 0; v < N; v++)
			// for (int index = 0; index < N; index++)
			for (int v = 0; v < N; v++)
			{
				results.emplace_back(
					pool.enqueue([v, &L] { // pass const type value j to thread; [] can be empty
						mtx_599[max_N_ID_for_mtx_599 - 1].lock();
						auto startTime = std::chrono::steady_clock::now();
						int used_id = Qid_599.front();
						Qid_599.pop();
						mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

						std::vector<two_hop_label> &Lv_final = Lv_final_599[v];

						/**
						 * get the L result of the current vertex
						 */

						std::vector<two_hop_label> &Lv = L[v];

						/**
						 * the temp_L in this thread
						 */
						auto &T = Temp_L_vk_599[used_id];

						/**
						 * Traverse the L-list of the current vertex
						 */
						for (const auto &Lvi : Lv)
						{
							int u = Lvi.hub_vertex;
							int u_hop = Lvi.hop;

							/**
							 * Traverse the L on the opposite vertex of the current label.
							 */
							const auto &Lu = L[u];

							/**
							 * traverse downward from the perfectly correct first vertex
							 */
							int min_dis = std::numeric_limits<int>::max();
							for (auto &label1 : Lu)
							{
								for (auto &label2 : T[label1.hub_vertex])
								{
									if (label1.hop + label2.second <= u_hop)
									{
										long long int query_dis = label1.distance + (long long int)label2.first;
										if (query_dis < min_dis)
										{
											min_dis = query_dis;
										}
									}
								}
							}

							if (min_dis > Lvi.distance)
							{
								Lv_final.push_back(two_hop_label(Lvi));
								++globalLabelCleanSize;
								T[u].push_back({Lvi.distance, Lvi.hop});
							}
						}

						for (const auto &label : Lv_final)
						{
							std::vector<std::pair<int, int>>().swap(T[label.hub_vertex]);
						}

						mtx_599[max_N_ID_for_mtx_599 - 1].lock();
						Qid_599.push(used_id);
						auto endTime = std::chrono::steady_clock::now();
						std::cout << "print pll v: " << v << " time cost is " << std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count() << std::endl;
						mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

						return 1; // return to results; the return type must be the same with results
					}));
			}

			for (auto &&result : results)
				result.get(); // all threads finish here
			std::cout << "start move" << std::endl;
			case_info.L = std::move(Lv_final_599);
			std::cout << "end move" << std::endl;
			results.clear();
		}

		template <typename weight_type>
		void hop_constrained_clear_global_values()
		{
			std::vector<std::vector<two_hop_label>>().swap(L_temp_599);
			std::vector<std::vector<two_hop_label>>().swap(Lv_final_599);
			ideal_graph_599<weight_type>.clear();
			std::vector<std::vector<std::vector<std::pair<int, int>>>>().swap(Temp_L_vk_599);
			std::vector<std::vector<std::pair<int, int>>>().swap(dist_hop_599);
			std::vector<std::vector<std::vector<std::pair<hop_constrained_node_handle, int>>>>().swap(Q_handle_priorities_599);
			std::vector<std::vector<std::vector<int>>>().swap(Vh_599);
			std::queue<int>().swap(Qid_599);
			PPR_TYPE::PPR_type().swap(PPR_599);
		}

		template <typename weight_type>
		void pll(graph<weight_type> &graph, hop::two_hop_case_info &case_info)
		{
			//----------------------------------- step 1: initialization -----------------------------------
			timer.startSubtask("step 1: initialization");
			int N = graph.size();
			/* store the L Label and PPR*/
			L_temp_599.resize(N);
			Lv_final_599.resize(N);
			PPR_599.resize(N);

			int num_of_threads = case_info.thread_num;
			ThreadPool pool(num_of_threads);
			std::vector<std::future<int>> results;
			timer.endSubtask();
			//----------------------------------------------- step 2: generate labels ---------------------------------------------------------------
			/**
			 * Use Temp_L_vk_599 to mark the positional relationship between the iterated nodes
			 * and the nodes with already generated labels. This is done to reduce the process of
			 * traversing the L labels of the iterated nodes. Additionally, register multithreaded
			 * tasks and retrieve the results
			 */
			timer.startSubtask("step 2: generate labels");
			global_upper_k = case_info.upper_k;
			ideal_graph_599<weight_type> = graph;
			Temp_L_vk_599.resize(num_of_threads);
			dist_hop_599.resize(num_of_threads);
			Q_handle_priorities_599.resize(num_of_threads);
			Vh_599.resize(num_of_threads);
			hop_constrained_node_handle handle_x;
			for (int i = 0; i < num_of_threads; i++)
			{
				Temp_L_vk_599[i].resize(N);
				dist_hop_599[i].resize(N, {std::numeric_limits<int>::max(), 0});
				Q_handle_priorities_599[i].resize(N);
				for (int j = 0; j < N; j++)
				{
					Q_handle_priorities_599[i][j].resize(global_upper_k + 1, {handle_x, std::numeric_limits<int>::max()});
				}
				Vh_599[i].resize(global_upper_k + 2);
				Qid_599.push(i);
			}

			int last_check_vID = N - 1;

			for (int v_k = 0; v_k <= last_check_vID; v_k++)
			{
				results.emplace_back(
					pool.enqueue([v_k]
								 {
							HSDL_thread_function<weight_type>(v_k);
							return 1; }));
			}

			for (auto &&result : results)
				result.get();
			timer.endSubtask();
			//----------------------------------------------- step 3: sortL ---------------------------------------------------------------
			timer.startSubtask("step 3: sortL");
			case_info.L = L_temp_599;
			case_info.L = hop_constrained_sortL(num_of_threads);
			case_info.PPR = PPR_599;
			timer.endSubtask();
			//----------------------------------------------- step 4: canonical_repair---------------------------------------------------------------
			timer.startSubtask("step 4: canonical_repair");
			hop_constrained_clean_L(case_info, num_of_threads);
			timer.endSubtask();
			//---------------------------------------------------------------------------------------------------------------------------------------
			hop_constrained_clear_global_values<weight_type>();
			std::cout << "end" << std::endl;
		}

	}
}
