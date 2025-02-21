#include "Historical/graph_with_time_span/graph.h"
#include <filesystem>
#include "Historical/experiment/experiment_config.h"
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "Historical/graph_with_time_span/two_hop_label.h"
#include <boost/heap/fibonacci_heap.hpp>
#include <shared_mutex>
#include <CPU/tool_functions/ThreadPool.h>

namespace experiment {
	boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };
	int max_N_ID_for_mtx_595 = 1e7;
	std::vector<std::shared_mutex> mtx_595(max_N_ID_for_mtx_595);
	template <typename weight_type>
	void read_graph(graph<weight_type>& graph, ExperimentConfig& config) {
		std::string readPath = config.data_source.string();
		std::string line_content;
		boost::random::uniform_int_distribution<> random_v;
		boost::random::uniform_int_distribution<> random_weight = boost::random::uniform_int_distribution<>(0, 200);
		int v_num = 0;
		// 读取文件
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
						random_v = boost::random::uniform_int_distribution<>(0, v_num);
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
	namespace nonhop {
		std::vector<std::vector<two_hop_label>> L_temp_595;
		PPR_TYPE::PPR_type PPR_595;
		std::vector<std::vector<int>> P_dij_595;
		std::vector<std::vector<int>> T_dij_595;
		typedef typename boost::heap::fibonacci_heap<two_hop_label>::handle_type PLL_handle_t_for_sp;
		std::vector<std::vector<PLL_handle_t_for_sp>> Q_handles_595;
		std::queue<int> Qid_595;

		void PLL_dij_function(int v_k, graph<int>& input_graph)
		{
			mtx_595[max_N_ID_for_mtx_595 - 1].lock();
			int used_id = Qid_595.front();
			Qid_595.pop();
			mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

			std::vector<int> P_changed_vertices, T_changed_vertices;
			std::vector<int>& T_dij = T_dij_595[used_id], P_dij = P_dij_595[used_id];
			std::vector<PLL_handle_t_for_sp>& Q_handles = Q_handles_595[used_id];

			boost::heap::fibonacci_heap<two_hop_label> Q;
			two_hop_label node(0);
			node.vertex = v_k;
			node.distance = 0;
			Q_handles[v_k] = Q.push(node);
			P_dij[v_k] = 0;
			P_changed_vertices.push_back(v_k);

			mtx_595[v_k].lock_shared();
			for (auto xx : L_temp_595[v_k])
			{ // 因为v-k的标签在从自己出发的过程中不会发生改变，并且在求query的过程中每次都会用到，所以可以提前取出来放在T数组，节省后面查找的时间
				int L_v_k_i_vertex = xx.vertex;
				T_dij[L_v_k_i_vertex] = xx.distance; // allocate T values for L_temp_595[v_k]
				T_changed_vertices.push_back(L_v_k_i_vertex);
			}
			mtx_595[v_k].unlock_shared();

			int new_label_num = 0;

			while (Q.size())
			{

				node = Q.top();
				Q.pop();
				int u = node.vertex;

				int P_u = node.distance;

				int query_v_k_u = std::numeric_limits<int>::max();
				int common_hub_for_query_v_k_u = 0;
				mtx_595[u].lock_shared(); // put lock in for loop is very slow
				for (auto xx : L_temp_595[u])
				{
					long long int dis = xx.distance + (long long int)T_dij[xx.vertex]; // long long int is to avoid overflow
					if (query_v_k_u > dis)
					{
						query_v_k_u = dis;
						common_hub_for_query_v_k_u = xx.vertex;
					}
				} // 求query的值
				mtx_595[u].unlock_shared();

				if (P_u < query_v_k_u)
				{
					node.vertex = v_k;
					node.distance = P_u;

					mtx_595[u].lock();
					L_temp_595[u].push_back(node); // 并行时L_temp_595[u]里面的标签不一定是按照vertex ID排好序的，但是因为什么query时用了T_dij_595的trick，没必要让L_temp_595[u]里面的标签排好序
					mtx_595[u].unlock();
					new_label_num++;

					/*下面是dij更新邻接点的过程，同时更新优先队列和距离*/
					for (auto xx : input_graph.ADJs[u])
					{
						int adj_v = xx.first, ec = xx.second;
						mtx_595[adj_v].lock();
						if (P_dij[adj_v] == std::numeric_limits<int>::max())
						{ // 尚未到达的点
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
						mtx_595[adj_v].unlock();
					}
				}
				//else if (PLL_dynamic_generate_PPR)
				else
				{
					if (common_hub_for_query_v_k_u != v_k)
					{
						mtx_595[u].lock();
						PPR_TYPE::PPR_insert(PPR_595, u, common_hub_for_query_v_k_u, v_k);
						mtx_595[u].unlock();
					}
					if (common_hub_for_query_v_k_u != u)
					{
						mtx_595[v_k].lock();
						PPR_TYPE::PPR_insert(PPR_595, v_k, common_hub_for_query_v_k_u, u);
						mtx_595[v_k].unlock();
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
						std::vector<two_hop_label>(L_temp_595[v_k]).swap(L_temp_595[v_k]); // swap释放vector中多余空间
						output_L[v_k] = L_temp_595[v_k];
						std::vector<two_hop_label>().swap(L_temp_595[v_k]); // clear new labels for RAM efficiency

						return 1; // return to results; the return type must be the same with results
						}));
			}
			for (auto&& result : results)
				result.get(); // all threads finish here
			results.clear();

			return output_L;
		}

		void clean_L(two_hop_case_info& case_info, int thread_num)
		{

			auto& L = case_info.L;
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

						std::vector<two_hop_label> Lv_final;

						mtx_595[v].lock_shared();
						std::vector<two_hop_label> Lv = L[v];
						mtx_595[v].unlock_shared();

						auto& T = T_dij_595[used_id];

						for (auto Lvi : Lv)
						{
							int u = Lvi.vertex;
							if (v == u)
							{
								Lv_final.push_back(Lvi);
								T[v] = Lvi.distance;
								continue;
							}
							mtx_595[u].lock_shared();
							auto Lu = L[u];
							mtx_595[u].unlock_shared();

							int min_dis = std::numeric_limits<int>::max();
							for (auto label : Lu)
							{
								long long int query_dis = label.distance + (long long int)T[label.vertex];
								if (query_dis < min_dis)
								{
									min_dis = query_dis;
								}
							}

							if (min_dis > Lvi.distance)
							{
								Lv_final.push_back(Lvi);
								T[u] = Lvi.distance;
							}
						}

						for (auto label : Lv_final)
						{
							T[label.vertex] = std::numeric_limits<int>::max();
						}

						mtx_595[v].lock();
						L[v] = Lv_final;
						L[v].shrink_to_fit();
						mtx_595[v].unlock();

						mtx_595[max_N_ID_for_mtx_595 - 1].lock();
						Qid_595.push(used_id);
						mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

						return 1; // return to results; the return type must be the same with results
						}));
			}

			for (auto&& result : results)
				result.get(); // all threads finish here
			results.clear();
		}

		void PLL_clear_global_values()
		{
			std::vector<std::vector<two_hop_label>>().swap(L_temp_595);
			PPR_TYPE::PPR_type().swap(PPR_595);
			std::queue<int>().swap(Qid_595);
			std::vector<std::vector<int>>().swap(P_dij_595);
			std::vector<std::vector<int>>().swap(T_dij_595);
			std::vector<std::vector<PLL_handle_t_for_sp>>().swap(Q_handles_595);
		}

		template <typename weight_type>
		void pll(graph<weight_type>& input_graph, nonhop::two_hop_case_info& case_info) {
			//----------------------------------- step 1: initialization ------------------------------------------------------------------
			timer.startSubtask("step 1: initialization");
			int num_of_threads = case_info.thread_num;
			int N = input_graph.ADJs.size();

			L_temp_595.resize(N);
			PPR_595.resize(N);
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
				for (auto&& result : results)
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

	namespace hop {
		template <typename weight_type>
		void pll(graph<weight_type>& graph, hop::two_hop_case_info& case_info) {

		}
	}

}