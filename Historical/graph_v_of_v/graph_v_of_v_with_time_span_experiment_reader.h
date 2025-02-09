#pragma once
#include <vector>
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span_non_hop_constrained.h"
#include "Historical/experiment/mark_timer.h"

#include <filesystem>
#include <queue>
#include <fstream>
#include <numeric>
#include <functional>

void PLL_experiment_clear_global_values()
{
	this_parallel_PLL_is_running_595 = false;
	vector<vector<two_hop_label>>().swap(L_temp_595);
	PPR_type().swap(PPR_595);
	queue<int>().swap(Qid_595);
	vector<vector<int>>().swap(P_dij_595);
	vector<vector<int>>().swap(T_dij_595);
	vector<vector<PLL_handle_t_for_sp>>().swap(Q_handles_595);
};
void initialize_experiment_global_values_dynamic(int N, int thread_num)
{
	Dis.resize(thread_num);
	Q_value.resize(thread_num);
	Q_handles.resize(thread_num);
	queue<int>().swap(Qid_595);
	for (int i = 0; i < thread_num; i++)
	{
		Dis[i].resize(N, { -1, -1 });
		Q_value[i].resize(N, MAX_VALUE);
		Q_handles[i].resize(N);
		Qid_595.push(i);
	}
};

namespace fs = std::filesystem;
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };
struct change_edge_info
{
	int v1;
	int v2;
	int weight;
	int time;
};

class experiment_config
{

private:
	const fs::path experiment_path;
	std::ofstream outputFile;
	//std::ofstream resultOutputFile;
	const fs::path save_path;
	int iteration;
	const int change_num;
	const int upper;
	const int lower;
	const bool is_debug;
	int v_num = 0;
	int e_num = 0;
	boost::random::uniform_int_distribution<> random_v;
	boost::random::uniform_int_distribution<> random_weight;
	// 保存每一个time slot的变化队列
	vector<std::queue<change_edge_info>> q_list;

	two_hop_case_info mm;
	mark_timer mm_mark_timer;

	two_hop_case_info mm2021;
	mark_timer mm2021_mark_timer;
	//void txt_save_label(vector<vector<two_hop_label>>& labels, fs::path _save_path) {
	//	std::ofstream label_save_path;
	//	label_save_path.precision(10);
	//	label_save_path.setf(std::ios::fixed);
	//	label_save_path.setf(std::ios::showpoint);
	//	label_save_path.open(_save_path.generic_u8string());
	//	for (int i = 0; i < labels.size(); i++) {
	//		for (const two_hop_label& label : labels[i]) {
	//			this->outputFile << i << " " << label.vertex << " " << label.distance << " " << label.t_s << " " << label.t_s << "\n";
	//		}
	//	}
	//}

	//void txt_read_label(fs::path _read_path, vector<vector<two_hop_label>> labels) {
	//	std::ifstream myfile(_read_path);
	//}

	void txt_read_base()
	{
		std::string readPath = this->is_debug ? this->save_path.generic_u8string() : this->experiment_path.generic_u8string();
		std::string line_content;
		graph_v_of_v<int> instance_graph;
		// 读取文件
		std::ifstream myfile(readPath);
		if (myfile.is_open())
		{
			while (getline(myfile, line_content))
			{
				if (this->is_debug)
				{
					std::vector<std::string> Parsed_content = parse_string(line_content, " ");
					if (!Parsed_content[0].compare("time"))
					{
						this->iteration = std::stoi(Parsed_content[1]);
						q_list = vector<std::queue<change_edge_info>>(this->iteration + 1, std::queue<change_edge_info>());
						continue;
					}
					if (!Parsed_content[0].compare("vertex"))
					{
						this->v_num = std::stoi(Parsed_content[1]);
						instance_graph.resize(this->v_num);
						continue;
					}
					if (!Parsed_content[0].compare("EOF"))
					{
						break;
					}
					int v1 = std::stoi(Parsed_content[0]);
					int v2 = std::stoi(Parsed_content[1]);
					int w = std::stoi(Parsed_content[2]);
					int time = std::stoi(Parsed_content[3]);
					if (time > 0)
					{
						this->q_list[time].push({ v1, v2, w, time });
					}
					else
					{
						instance_graph.add_edge(v1, v2, w);
					}
				}
				else
				{
					std::vector<std::string> Parsed_content = parse_string(line_content, "\t");

					if (!Parsed_content[0].compare("#"))
					{
						if (!Parsed_content[1].compare("Nodes"))
						{
							instance_graph.resize(std::stoi(Parsed_content[2]));
							this->v_num = std::stoi(Parsed_content[2]);
							this->outputFile << "vertex" << " " << this->v_num << std::endl;
							this->random_v = boost::random::uniform_int_distribution<>(0, this->v_num - 1);
						}
						else if (!Parsed_content[1].compare("Edges"))
						{
							this->e_num = std::stoi(Parsed_content[2]);
							// weight_type weight_upper_limit, weight_type weight_lower_limit
						}
					}
					else
					{
						int v1 = std::stoi(Parsed_content[0]);
						int v2 = std::stoi(Parsed_content[1]);
						int w = this->random_weight(boost_random_time_seed);
						instance_graph.add_edge(v1, v2, w);
						this->txt_save(v1, v2, w, 0);
					}
				}
			}

			this->graph_with_time_span = graph_v_of_v_with_time_span(this->v_num, this->e_num, this->upper, this->lower, this->iteration);
			this->graph_with_time_span.add_graph_time(instance_graph, 0);
			this->graphs.push_back(instance_graph);
			if (!this->is_debug)
			{
				// 迭代指定次数 生成随机改变的边的数组 并保存到硬盘
				for (int i = 1; i <= this->iteration; i++)
				{
					std::map<pair<int, int>, int> diff;
					std::cout << i << std::endl;
					for (int j = 0; j < this->change_num; j++)
					{
						int index_i = this->random_v(boost_random_time_seed);
						boost::random::uniform_int_distribution<> dis_inner(0, this->graphs[0].ADJs[index_i].size() - 1);
						int index_j = dis_inner(boost_random_time_seed);
						int i_j_weight = this->random_weight(boost_random_time_seed);
						diff[{index_i, index_j}] = i_j_weight;
						q_list[i].push({ index_i, index_j, i_j_weight, i });
						// 持久化
						txt_save(index_i, index_j, i_j_weight, i);
					}
					for (auto& iter : diff) {
						int v1 = iter.first.first;
						int v2 = iter.first.second;
						int w = iter.second;

						q_list[i].push({ v1, v2, w, i });
						// 持久化
						txt_save(v1, v2, w, i);
					}
				}
			}
			myfile.close(); // close the file
			if (!this->is_debug)
			{
				this->txt_close();
			}
		}
		else
		{
			std::cout << "Unable to open file " << readPath << std::endl
				<< "Please check the file location or file name." << std::endl; // throw an error message
			getchar();                                                                // keep the console window
			exit(1);                                                                  // end the program
		}
	}
	void txt_save(int v1, int v2, int w, int time)
	{
		this->outputFile << v1 << " " << v2 << " " << w << " " << time << "\n";
	};
	void txt_close()
	{
		this->outputFile << "EOF" << std::endl;
		this->outputFile.close();
	}

public:
	graph_v_of_v_with_time_span graph_with_time_span;
	vector<graph_v_of_v<int>> graphs;
	experiment_config(fs::path _experiment_path, fs::path _save_path, int _iteration, int _change_num, bool _is_debug, int _upper, int _lower) : experiment_path(_experiment_path), save_path(_save_path), iteration(_iteration), change_num(_change_num), is_debug(_is_debug), upper(_upper), lower(_lower)
	{
		random_weight = boost::random::uniform_int_distribution<>(lower, upper);
		if (!this->is_debug)
		{
			q_list = vector<std::queue<change_edge_info>>(this->iteration + 1, std::queue<change_edge_info>());
			this->outputFile.precision(10);
			this->outputFile.setf(std::ios::fixed);
			this->outputFile.setf(std::ios::showpoint);
			this->outputFile.open(_save_path.generic_u8string());
			this->outputFile << "time" << " " << this->iteration << std::endl;
		}
	}

	int init()
	{
		mm.max_labal_byte_size = 6e9;
		mm.max_run_time_seconds = 1e4;
		mm.use_2M_prune = 1;
		mm.use_rank_prune = 1;
		mm.use_canonical_repair = 1;
		mm.thread_num = 10;
		mm2021.max_labal_byte_size = 6e9;
		mm2021.max_run_time_seconds = 1e4;
		mm2021.use_2M_prune = 1;
		mm2021.use_rank_prune = 1;
		mm2021.use_canonical_repair = 1;
		mm2021.thread_num = 10;
		// 读取图的数据到指定的对象中 graphs为原始图的列表 graph_with_time_span为时序图的对象
		this->txt_read_base();
		return 0;
	}
	int process()
	{
		// 初始化nonhop
		initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
		mm_mark_timer.mark();
		PLL(this->graphs[0], this->mm);
		mm_mark_timer.add();
		mm_mark_timer.push();
		std::cout << "PLL1 finished" << std::endl;
		PLL_experiment_clear_global_values();
		//initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
		//mm2021_mark_timer.mark();
		//PLL(this->graphs[0], this->mm2021);
		//mm2021_mark_timer.add();
		//mm2021_mark_timer.push();
		//std::cout << "PLL2 finished" << std::endl;
		PLL_experiment_clear_global_values();
		// 动态维护
		int time = 0;
		vector<pair<int, int>> path_decrease;
		vector<int> weight_decrease;

		vector<pair<int, int>> path_increase;
		vector<int> weight_increase;
		vector<int> weight_old_increase;

		ThreadPool pool_dynamic(mm.thread_num);
		std::vector<std::future<int>> results_dynamic;

		for (int i = 1; i <= this->iteration; i++)
		{
			std::cout << "iteration " << i << std::endl;
			initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
			std::queue<change_edge_info> q = this->q_list[i];
			graph_v_of_v<int> instance_graph_temp(this->graphs[i - 1]);
			while (!q.empty())
			{
				change_edge_info info = q.front();
				q.pop();
				// 1. 读取数据
				int v1 = info.v1;
				int v2 = info.v2;
				int weight = info.weight;
				if (this->graphs[i - 1].ADJs[v1][v2].second < weight)
				{
					// increase
					int old_weight = this->graphs[i - 1].ADJs[v1][v2].second;
					// instance_graph_temp.add_edge(v1, this->graphs[i - 1].ADJs[v1][v2].first, weight);
					path_increase.push_back({ v1, this->graphs[i - 1].ADJs[v1][v2].first });
					weight_old_increase.push_back(old_weight);
					weight_increase.push_back(weight);
				}
				else if (this->graphs[i - 1].ADJs[v1][v2].second > weight)
				{
					// instance_graph_temp.add_edge(v1, this->graphs[i - 1].ADJs[v1][v2].first, weight);
					path_decrease.push_back({ v1, this->graphs[i - 1].ADJs[v1][v2].first });
					weight_decrease.push_back(weight);
				}
				if (path_decrease.size() > mm.thread_num)
				{
					for (int i = 0; i < path_decrease.size(); i++)
					{
						int v1 = path_decrease[i].first;
						int v2 = path_decrease[i].second;
						int w = weight_decrease[i];
						instance_graph_temp.add_edge(v1, v2, w);
					}
					initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
					std::cout << "decrease ruc maintain" << std::endl;
					mm_mark_timer.mark();
					nonHOP_WeightDecreaseMaintenance_improv_batch(instance_graph_temp, mm, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
					mm_mark_timer.add();
					std::cout << "decrease ruc maintain finished" << std::endl;
					//initialize_experiment_global_values_dynamic(this->v_num, this->mm2021.thread_num);
					//mm2021_mark_timer.mark();
					//std::cout << "decrease 2021 maintain" << std::endl;
					//nonHOP_WeightDecrease2021_batch(instance_graph_temp, mm2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
					//mm2021_mark_timer.add();
					//std::cout << "decrease 2021 maintain finished" << std::endl;
					vector<pair<int, int>>().swap(path_decrease);
					vector<int>().swap(weight_decrease);
				}
				if (path_increase.size() > mm.thread_num)
				{
					for (int i = 0; i < path_increase.size(); i++)
					{
						int v1 = path_increase[i].first;
						int v2 = path_increase[i].second;
						int w = weight_increase[i];
						int w_old = weight_old_increase[i];
						instance_graph_temp.add_edge(v1, v2, w);
					}
					initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
					std::cout << "increase ruc maintain" << std::endl;
					mm_mark_timer.mark();
					nonHOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, mm, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
					mm_mark_timer.add();
					std::cout << "increase ruc maintain finished" << std::endl;
					//initialize_experiment_global_values_dynamic(this->v_num, this->mm2021.thread_num);
					//mm2021_mark_timer.mark();
					//std::cout << "increase 2021 maintain" << std::endl;
					//nonHOP_WeightIncrease2021_batch(instance_graph_temp, mm2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
					//mm2021_mark_timer.add();
					//std::cout << "increase 2021 maintain finished" << std::endl;
					vector<pair<int, int>>().swap(path_increase);
					vector<int>().swap(weight_increase);
					vector<int>().swap(weight_old_increase);
				}
			}
			if (path_decrease.size() > 0)
			{
				for (int i = 0; i < path_decrease.size(); i++)
				{
					int v1 = path_decrease[i].first;
					int v2 = path_decrease[i].second;
					int w = weight_decrease[i];
					instance_graph_temp.add_edge(v1, v2, w);
				}
				initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
				std::cout << "decrease ruc maintain" << std::endl;
				mm_mark_timer.mark();
				nonHOP_WeightDecreaseMaintenance_improv_batch(instance_graph_temp, mm, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
				mm_mark_timer.add();
				std::cout << "decrease ruc maintain finished" << std::endl;
				//initialize_experiment_global_values_dynamic(this->v_num, this->mm2021.thread_num);
				//mm2021_mark_timer.mark();
				//std::cout << "decrease 2021 maintain" << std::endl;
				//nonHOP_WeightDecrease2021_batch(instance_graph_temp, mm2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
				//mm2021_mark_timer.add();
				//std::cout << "decrease 2021 maintain finished" << std::endl;
				vector<pair<int, int>>().swap(path_decrease);
				vector<int>().swap(weight_decrease);
			}
			if (path_increase.size() > 0)
			{
				for (int i = 0; i < path_increase.size(); i++)
				{
					int v1 = path_increase[i].first;
					int v2 = path_increase[i].second;
					int w = weight_increase[i];
					int w_old = weight_old_increase[i];
					instance_graph_temp.add_edge(v1, v2, w);
				}
				initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
				std::cout << "increase ruc maintain finished" << std::endl;
				mm_mark_timer.mark();
				nonHOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, mm, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
				mm_mark_timer.add();
				//initialize_experiment_global_values_dynamic(this->v_num, this->mm2021.thread_num);
				//mm2021_mark_timer.mark();
				//std::cout << "increase 2021 maintain finished" << std::endl;
				//nonHOP_WeightIncrease2021_batch(instance_graph_temp, mm2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
				//mm2021_mark_timer.add();
				vector<pair<int, int>>().swap(path_increase);
				vector<int>().swap(weight_increase);
				vector<int>().swap(weight_old_increase);
			}
			mm_mark_timer.push();
			this->graphs.push_back(instance_graph_temp);
			PLL_experiment_clear_global_values();
		}
		return 0;
	}

	int process_2021() {
		// 初始化nonhop
		/*initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
		mm_mark_timer.mark();
		PLL(this->graphs[0], this->mm);
		mm_mark_timer.add();
		mm_mark_timer.push();
		std::cout << "PLL1 finished" << std::endl;
		PLL_experiment_clear_global_values();*/
		initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
		mm2021_mark_timer.mark();
		PLL(this->graphs[0], this->mm2021);
		mm2021_mark_timer.add();
		mm2021_mark_timer.push();
		std::cout << "PLL2 finished" << std::endl;
		PLL_experiment_clear_global_values();
		// 动态维护
		int time = 0;
		vector<pair<int, int>> path_decrease;
		vector<int> weight_decrease;

		vector<pair<int, int>> path_increase;
		vector<int> weight_increase;
		vector<int> weight_old_increase;

		ThreadPool pool_dynamic(mm.thread_num);
		std::vector<std::future<int>> results_dynamic;

		for (int i = 1; i <= this->iteration; i++)
		{
			std::cout << "============iteration " << i << "==============" << std::endl;
			initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
			std::queue<change_edge_info> q = this->q_list[i];
			graph_v_of_v<int> instance_graph_temp(this->graphs[i - 1]);
			while (!q.empty())
			{
				change_edge_info info = q.front();
				q.pop();
				// 1. 读取数据
				int v1 = info.v1;
				int v2 = info.v2;
				int weight = info.weight;
				if (this->graphs[i - 1].ADJs[v1][v2].second < weight)
				{
					// increase
					int old_weight = this->graphs[i - 1].ADJs[v1][v2].second;
					// instance_graph_temp.add_edge(v1, this->graphs[i - 1].ADJs[v1][v2].first, weight);
					path_increase.push_back({ v1, this->graphs[i - 1].ADJs[v1][v2].first });
					weight_old_increase.push_back(old_weight);
					weight_increase.push_back(weight);
				}
				else if (this->graphs[i - 1].ADJs[v1][v2].second > weight)
				{
					// instance_graph_temp.add_edge(v1, this->graphs[i - 1].ADJs[v1][v2].first, weight);
					path_decrease.push_back({ v1, this->graphs[i - 1].ADJs[v1][v2].first });
					weight_decrease.push_back(weight);
				}
				if (path_decrease.size() > mm.thread_num)
				{
					for (int i = 0; i < path_decrease.size(); i++)
					{
						int v1 = path_decrease[i].first;
						int v2 = path_decrease[i].second;
						int w = weight_decrease[i];
						instance_graph_temp.add_edge(v1, v2, w);
					}
					//initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
					//mm_mark_timer.mark();
					//std::cout << "decrease ruc maintain" << std::endl;
					//nonHOP_WeightDecreaseMaintenance_improv_batch(instance_graph_temp, mm, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
					//mm_mark_timer.add();
					//std::cout << "decrease ruc maintain finished" << std::endl;
					initialize_experiment_global_values_dynamic(this->v_num, this->mm2021.thread_num);
					std::cout << "decrease 2021 maintain" << std::endl;
					mm2021_mark_timer.mark();
					nonHOP_WeightDecrease2021_batch(instance_graph_temp, mm2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
					mm2021_mark_timer.add();
					std::cout << "decrease 2021 maintain finished" << std::endl;
					vector<pair<int, int>>().swap(path_decrease);
					vector<int>().swap(weight_decrease);
				}
				if (path_increase.size() > mm.thread_num)
				{
					for (int i = 0; i < path_increase.size(); i++)
					{
						int v1 = path_increase[i].first;
						int v2 = path_increase[i].second;
						int w = weight_increase[i];
						int w_old = weight_old_increase[i];
						instance_graph_temp.add_edge(v1, v2, w);
					}
					//initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
					//mm_mark_timer.mark();
					//std::cout << "increase ruc maintain" << std::endl;
					//nonHOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, mm, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
					//mm_mark_timer.add();
					//std::cout << "increase ruc maintain finished" << std::endl;
					initialize_experiment_global_values_dynamic(this->v_num, this->mm2021.thread_num);
					std::cout << "increase 2021 maintain" << std::endl;
					mm2021_mark_timer.mark();
					nonHOP_WeightIncrease2021_batch(instance_graph_temp, mm2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
					mm2021_mark_timer.add();
					std::cout << "increase 2021 maintain finished" << std::endl;
					vector<pair<int, int>>().swap(path_increase);
					vector<int>().swap(weight_increase);
					vector<int>().swap(weight_old_increase);
				}
			}
			if (path_decrease.size() > 0)
			{
				for (int i = 0; i < path_decrease.size(); i++)
				{
					int v1 = path_decrease[i].first;
					int v2 = path_decrease[i].second;
					int w = weight_decrease[i];
					instance_graph_temp.add_edge(v1, v2, w);
				}
				//initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
				//mm_mark_timer.mark();
				//std::cout << "decrease ruc maintain" << std::endl;
				//nonHOP_WeightDecreaseMaintenance_improv_batch(instance_graph_temp, mm, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
				//mm_mark_timer.add();
				//std::cout << "decrease ruc maintain finished" << std::endl;
				initialize_experiment_global_values_dynamic(this->v_num, this->mm2021.thread_num);
				std::cout << "decrease 2021 maintain" << std::endl;
				mm2021_mark_timer.mark();
				nonHOP_WeightDecrease2021_batch(instance_graph_temp, mm2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
				mm2021_mark_timer.add();
				std::cout << "decrease 2021 maintain finished" << std::endl;
				vector<pair<int, int>>().swap(path_decrease);
				vector<int>().swap(weight_decrease);
			}
			if (path_increase.size() > 0)
			{
				for (int i = 0; i < path_increase.size(); i++)
				{
					int v1 = path_increase[i].first;
					int v2 = path_increase[i].second;
					int w = weight_increase[i];
					int w_old = weight_old_increase[i];
					instance_graph_temp.add_edge(v1, v2, w);
				}
				//initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
				//mm_mark_timer.mark();
				//std::cout << "increase ruc maintain finished" << std::endl;
				//nonHOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, mm, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
				//mm_mark_timer.add();
				initialize_experiment_global_values_dynamic(this->v_num, this->mm2021.thread_num);
				std::cout << "increase 2021 maintain finished" << std::endl;
				mm2021_mark_timer.mark();
				nonHOP_WeightIncrease2021_batch(instance_graph_temp, mm2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
				mm2021_mark_timer.add();
				vector<pair<int, int>>().swap(path_increase);
				vector<int>().swap(weight_increase);
				vector<int>().swap(weight_old_increase);
			}
			mm2021_mark_timer.push();
			this->graphs.push_back(instance_graph_temp);
			PLL_experiment_clear_global_values();
		}
		return 0;
	}

	int print_experiment_result()
	{

		std::cout << "the 2024 maintain algorithm L size is " << mm.compute_L_byte_size() + mm.compute_PPR_byte_size() << std::endl;
		std::cout << "the 2021 maintain algorithm L size is " << mm2021.compute_L_byte_size() + mm2021.compute_PPR_byte_size() << std::endl;
		std::cout << "2021 result is" << mm2021.query(26100, 28900, 40, 60) << std::endl;
		std::cout << "ruc result is" << mm.query(26100, 28900, 40, 60) << std::endl;
		double slot0_2021 = 0;
		double slot1_2021 = 0;
		double slot0_2024 = 0;
		double slot1_2024 = 0;
		int pre = 0;
		int after = 0;
		for (int i = 1; i <= this->iteration; i++) {
			if (i < ceil(this->iteration / 2)) {
				pre++;
				slot0_2021 += this->mm2021_mark_timer.get_experiment_time()[i];
				slot0_2024 += this->mm_mark_timer.get_experiment_time()[i];
			}
			else {
				after++;
				slot1_2021 += this->mm2021_mark_timer.get_experiment_time()[i];
				slot1_2024 += this->mm_mark_timer.get_experiment_time()[i];
			}
		}
		std::cout << "In the 2024 algorithm, the index construction time is " << this->mm_mark_timer.get_experiment_time()[0]
			<< ", the maintenance time for slot0 is "
			<< slot0_2024 / pre
			<< ", and the maintenance time for slot1 is "
			<< slot1_2024 / after
			<< "." << std::endl;
		std::cout << "In the 2021 algorithm, the index construction time is" << this->mm2021_mark_timer.get_experiment_time()[0]
			<< ", the maintenance time for slot0 is "
			<< slot0_2021 / pre
			<< ", and the maintenance time for slot1 is "
			<< slot1_2021 / after
			<< std::endl;
		return 0;
	}

	int close()
	{
		return 0;
	}
};
