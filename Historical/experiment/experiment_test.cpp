#include <iostream>
#include <string>
#include <filesystem>
#include <Historical/utils/ExecutionTimer.h>
experiment::ExecutionTimer timer;
#include "Historical/utils/BinaryPersistence.h"
#include "Historical/graph_with_time_span/two_hop_label.h"
#include "Historical/experiment/experiment_operation.h"
#include "Historical/experiment/experiment_maintain_operation.h"
#include "Historical/experiment/experiment_config.h"
#include "Historical/graph_with_time_span/graph_search_baseline.h"
int main(int argc, char *argv[])
{
	try
	{
		experiment::ExperimentConfig config = experiment::parse_arguments(argc, argv);

		std::cout << "Mode: " << (config.mode == experiment::GENERATE_LABEL ? "Generate Label" : (config.mode == experiment::MAINTAIN_LABEL ? "Maintain Label" : "QueryResult")) << "\n"
				  << "Threads: " << config.threads << "\n"
				  << "Data Source: " << config.data_source << "\n"
				  << "Save Path: " << config.save_path << "\n"
				  << "Hop Limit (k): " << config.hop_limit << "\n";

		if (config.mode == experiment::MAINTAIN_LABEL)
		{
			std::cout << "Iterations: " << config.iterations << "\n"
					  << "Change Count: " << config.change_count << "\n"
					  << "Max Value: " << config.max_value << "\n"
					  << "Min Value: " << config.min_value << "\n";
		}
		if (config.mode == experiment::GENERATE_LABEL)
		{
			experiment::graph<int> instance_graph;
			timer.startTask("generate graph and 2hop label " + std::to_string(config.hop_limit));
			timer.startSubtask("generate graph " + std::to_string(config.hop_limit));
			experiment::read_graph(instance_graph, config);
			timer.endSubtask();
			experiment::graph_with_time_span<int> graph_time;
			graph_time.add_graph_time(instance_graph, 0);
			std::filesystem::path saveDir = std::filesystem::path(config.save_path);
			std::filesystem::create_directory(saveDir);
			if (config.hop_limit != 0)
			{
				// k-constrained pll
				std::string graph_res_filename = "binary_hop_constrained_" + std::to_string(config.hop_limit) + "_graph";
				std::string experiment_res_filename = "GENERATE_LABEL__hop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::filesystem::path graphPath = saveDir.string() + "//" + graph_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_res_filename;
				experiment::hop::two_hop_case_info hop_info;
				hop_info.thread_num = config.threads;
				hop_info.upper_k = config.hop_limit;
				timer.startSubtask("generate 2hop label " + std::to_string(config.hop_limit) + " hop constrained");
				experiment::hop::pll(instance_graph, hop_info);
				timer.endSubtask();
				// hop_info.print_L();
				std::ofstream FILE_GRAPH(graphPath.string(), std::ios::out | std::ofstream::binary);
				experiment::saveBinary(FILE_GRAPH, instance_graph);
				experiment::saveBinary(FILE_GRAPH, graph_time);
				experiment::saveBinary(FILE_GRAPH, hop_info);
				FILE_GRAPH.close();
				std::ofstream outFile;
				outFile.precision(6);
				outFile.setf(std::ios::fixed);
				outFile.setf(std::ios::showpoint);
				outFile.open(resultPath.string());
				timer.writeStatsToFile(outFile);
				hop_info.record_all_details_stream(outFile);
				outFile.close();
			}
			else
			{
				std::string graph_res_filename = "binary_nonhop_constrained_" + std::to_string(config.hop_limit) + "_graph";
				std::string experiment_res_filename = "GENERATE_LABEL_nonhop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::filesystem::path graphPath = saveDir.string() + "//" + graph_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_res_filename;
				experiment::nonhop::two_hop_case_info hop_info;
				hop_info.thread_num = config.threads;
				timer.startSubtask("generate graph and 2hop label " + std::to_string(config.hop_limit) + " nonhop constrained");
				experiment::nonhop::pll(instance_graph, hop_info);
				timer.endSubtask();
				// hop_info.print_L();
				std::ofstream FILE_GRAPH(graphPath.string(), std::ios::out | std::ofstream::binary);
				experiment::saveBinary(FILE_GRAPH, instance_graph);
				experiment::saveBinary(FILE_GRAPH, graph_time);
				experiment::saveBinary(FILE_GRAPH, hop_info);
				FILE_GRAPH.close();
				std::ofstream outFile;
				outFile.precision(6);
				outFile.setf(std::ios::fixed);
				outFile.setf(std::ios::showpoint);
				outFile.open(resultPath.string());
				timer.writeStatsToFile(outFile);
				hop_info.record_all_details_stream(outFile);
				outFile.close();
			}
		}
		else if (config.mode == experiment::MAINTAIN_LABEL)
		{
			experiment::ExecutionTimer timer_ruc;
			experiment::ExecutionTimer timer_2021;
			experiment::ExecutionTimer timer_baseline1;
			experiment::ExecutionTimer timer_baseline2;
			std::vector<experiment::graph<int>> graph_list;
			experiment::graph<int> init_graph;
			experiment::graph_with_time_span<int> graph_time;
			std::filesystem::path saveDir = std::filesystem::path(config.save_path);
			std::filesystem::create_directory(saveDir);
			timer_ruc.startTask("maintain graph and 2hop label " + std::to_string(config.hop_limit) + (config.hop_limit == 0 ? "nonhop_constrained" : "hop constrained"));
			timer_2021.startTask("maintain graph and 2hop label " + std::to_string(config.hop_limit) + (config.hop_limit == 0 ? "nonhop_constrained" : "hop constrained"));
			timer_baseline1.startTask("maintain graph by saving every graph");
			timer_baseline2.startTask("maintain graph by saving edge info with time label");
			timer_ruc.startSubtask("step-1 read graph and original 2hop label");
			timer_2021.startSubtask("step-1 read graph and original 2hop label");
			if (config.hop_limit != 0)
			{
				experiment::hop::two_hop_case_info hop_info;
				std::string graph_res_filename = "binary_hop_constrained_" + std::to_string(config.hop_limit) + "_graph";
				std::string dataSource = config.data_source.string() + "//" + graph_res_filename;

				std::string hop_label_res_filename = "binary_hop_constrained_" + std::to_string(config.hop_limit) + "_2_hop_label_info";
				std::string experiment_res_filename = "MAINTAIN_LABEL_hop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::string change_info_res_filename = "change_info_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::filesystem::path hopLabelPath = saveDir.string() + "//" + hop_label_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_res_filename;
				std::filesystem::path changePath = saveDir.string() + "//" + change_info_res_filename;

				std::ifstream FILE_GRAPH(dataSource, std::ios::in | std::ifstream::binary);
				experiment::loadBinary(FILE_GRAPH, init_graph);
				experiment::loadBinary(FILE_GRAPH, graph_time);
				experiment::loadBinary(FILE_GRAPH, hop_info);
				hop_info.thread_num = config.threads;
				hop_info.upper_k = config.hop_limit;
				experiment::hop::two_hop_case_info hop_info_2021;
				hop_info_2021 = hop_info;
				graph_list.push_back(init_graph);
				timer_ruc.endSubtask();
				timer_2021.endSubtask();

				experiment::iteration_info<int> change_info(
					graph_time.v_num, config.iterations, config.change_count, config.max_value, config.min_value, init_graph);
				change_info.build_random_change();
				std::ofstream CHANGE_PATH_STREAM(changePath.string(), std::ios::out | std::ofstream::binary);
				experiment::saveBinary(CHANGE_PATH_STREAM,change_info);

				// std::ifstream CHANGE_PATH_STREAM(changePath.string(), std::ios::in | std::ifstream::binary);
				// experiment::loadBinary(CHANGE_PATH_STREAM, change_info);
				std::vector<std::pair<int, int>> path_decrease;
				std::map<std::pair<int, int>, int> path2Index4Decrease;
				std::vector<int> weight_decrease;

				std::vector<std::pair<int, int>> path_increase;
				std::vector<int> weight_increase;
				std::vector<int> weight_old_increase;
				std::map<std::pair<int, int>, int> path2Index4Increase;

				ThreadPool pool_dynamic(hop_info.thread_num);
				std::vector<std::future<int>> results_dynamic;
				timer_ruc.startSubtask("step-2 maintain 2 hop label 0 - " + std::to_string(config.iterations / 2));
				timer_2021.startSubtask("step-2 maintain 2 hop label 0 - " + std::to_string(config.iterations / 2));
				timer_baseline1.startSubtask("step-1 save graph from 0 - " + std::to_string(config.iterations / 2));
				timer_baseline2.startSubtask("step-1 save graph with time span edge from 0 - " + std::to_string(config.iterations / 2));
				for (int i = 1; i <= config.iterations; i++)
				{
					if (i == config.iterations / 2 + 1)
					{
						timer_ruc.startSubtask("step-3 maintain 2 hop label " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_2021.startSubtask("step-3 maintain 2 hop label " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_baseline1.startSubtask("step-2 save graph from " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_baseline2.startSubtask("step-2 save graph with time span edge from " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
					}
					std::cout << "iteration " << i << std::endl;
					std::queue<experiment::change_edge_info> q = change_info.q_list[i];
					experiment::graph<int> instance_graph_temp = graph_list[i - 1];
					while (!q.empty())
					{
						experiment::change_edge_info info = q.front();
						q.pop();
						// 1. ��ȡ����
						int v1 = info.v1;
						int v2 = info.v2;
						int weight = info.weight;
						if (instance_graph_temp.ADJs[v1][v2].second < weight)
						{
							auto pairPathV = std::make_pair(v1, instance_graph_temp.ADJs[v1][v2].first);
							auto it = path2Index4Increase.find(pairPathV);
							// increase
							if (it == path2Index4Increase.end())
							{
								int old_weight = instance_graph_temp.ADJs[v1][v2].second;
								path_increase.push_back(pairPathV);
								weight_old_increase.push_back(old_weight);
								weight_increase.push_back(weight);
								path2Index4Increase[pairPathV] = weight_increase.size() - 1;
							}
							else
							{
								weight_increase[path2Index4Increase[pairPathV]] = weight;
							}
						}
						else if (instance_graph_temp.ADJs[v1][v2].second > weight)
						{
							auto pairPathV = std::make_pair(v1, instance_graph_temp.ADJs[v1][v2].first);
							auto it = path2Index4Decrease.find(pairPathV);
							if (it == path2Index4Decrease.end())
							{
								path_decrease.push_back({v1, instance_graph_temp.ADJs[v1][v2].first});
								weight_decrease.push_back(weight);
								path2Index4Decrease[pairPathV] = weight_decrease.size() - 1;
							}
							else
							{
								weight_decrease[path2Index4Decrease[pairPathV]] = weight;
							}
						}
						if (path_decrease.size() > hop_info.thread_num)
						{
							for (int i = 0; i < path_decrease.size(); i++)
							{
								int v1 = path_decrease[i].first;
								int v2 = path_decrease[i].second;
								int w = weight_decrease[i];
								instance_graph_temp.add_edge(v1, v2, w);
							}
							std::cout << "decrease ruc maintain" << std::endl;
							experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
							timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc decrease maintain");
							experiment::hop::ruc::decrease::HOP_WeightDecreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
							timer_ruc.endSubtask();
							std::cout << "decrease 2021 maintain" << std::endl;
							experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
							timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 decrease maintain");
							experiment::hop::algorithm2021::decrease::HOP_WeightDecrease2021_batch(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
							timer_2021.endSubtask();
							std::vector<std::pair<int, int>>().swap(path_decrease);
							std::vector<int>().swap(weight_decrease);
							std::map<std::pair<int, int>, int>().swap(path2Index4Decrease);
						}
						if (path_increase.size() > hop_info.thread_num)
						{
							for (int i = 0; i < path_increase.size(); i++)
							{
								int v1 = path_increase[i].first;
								int v2 = path_increase[i].second;
								int w = weight_increase[i];
								int w_old = weight_old_increase[i];
								instance_graph_temp.add_edge(v1, v2, w);
							}
							std::cout << "increase ruc maintain" << std::endl;
							experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
							timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc increase maintain");
							experiment::hop::ruc::increase::HOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
							timer_ruc.endSubtask();
							std::cout << "increase 2021 maintain" << std::endl;
							experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
							timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 increase maintain");
							experiment::hop::algorithm2021::increase::HOP_WeightIncrease2021_batch(instance_graph_temp, hop_info_2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
							timer_2021.endSubtask();
							std::vector<std::pair<int, int>>().swap(path_increase);
							std::vector<int>().swap(weight_increase);
							std::vector<int>().swap(weight_old_increase);
							std::map<std::pair<int, int>, int>().swap(path2Index4Increase);
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
						std::cout << "decrease ruc maintain" << std::endl;
						experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
						timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc decrease maintain");
						experiment::hop::ruc::decrease::HOP_WeightDecreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
						timer_ruc.endSubtask();
						std::cout << "decrease 2021 maintain" << std::endl;
						experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
						timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 decrease maintain");
						experiment::hop::algorithm2021::decrease::HOP_WeightDecrease2021_batch(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
						timer_2021.endSubtask();
						std::vector<std::pair<int, int>>().swap(path_decrease);
						std::vector<int>().swap(weight_decrease);
						std::map<std::pair<int, int>, int>().swap(path2Index4Decrease);
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
						std::cout << "increase ruc maintain" << std::endl;
						experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
						timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc increase maintain");
						experiment::hop::ruc::increase::HOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
						timer_ruc.endSubtask();
						std::cout << "increase 2021 maintain" << std::endl;
						experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
						timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 increase maintain");
						experiment::hop::algorithm2021::increase::HOP_WeightIncrease2021_batch(instance_graph_temp, hop_info_2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
						timer_2021.endSubtask();
						std::vector<std::pair<int, int>>().swap(path_increase);
						std::vector<int>().swap(weight_increase);
						std::vector<int>().swap(weight_old_increase);
						std::map<std::pair<int, int>, int>().swap(path2Index4Increase);
					}
					timer_baseline1.startSubtask("save graph " + std::to_string(i));
					graph_list.push_back(instance_graph_temp);
					timer_baseline1.endSubtask();
					timer_baseline2.startSubtask("save graph with time span label " + std::to_string(i));
					graph_time.add_graph_time(instance_graph_temp, i);
					timer_baseline2.endSubtask();
					if (i == config.iterations / 2 || i == config.iterations)
					{
						timer_ruc.endSubtask();
						timer_2021.endSubtask();
						timer_baseline1.endSubtask();
						timer_baseline2.endSubtask();
					}
					std::cout <<"current L size is " << hop_info.compute_label_bit_size() << std::endl;
				}
				std::cout << "finish maintain label" << std::endl;
				std::ofstream FILE_HOP_LABEL(hopLabelPath.string(), std::ios::out | std::ofstream::binary);
				experiment::saveBinary(FILE_HOP_LABEL, graph_list);
				experiment::saveBinary(FILE_HOP_LABEL, graph_time);
				experiment::saveBinary(FILE_HOP_LABEL, hop_info);
				experiment::saveBinary(FILE_HOP_LABEL, hop_info_2021);
				FILE_HOP_LABEL.close();
				std::ofstream outFile;
				outFile.precision(6);
				outFile.setf(std::ios::fixed);
				outFile.setf(std::ios::showpoint);
				outFile.open(resultPath.string());
				outFile << "========================ruc maintain======================" << std::endl;
				timer_ruc.writeStatsToFile(outFile);
				hop_info.record_all_details_stream(outFile);
				outFile << "========================2021 maintain=====================" << std::endl;
				timer_2021.writeStatsToFile(outFile);
				hop_info_2021.record_all_details_stream(outFile);
				outFile << "========================baseline1=========================" << std::endl;
				timer_baseline1.writeStatsToFile(outFile);
				long long int graph_list_size = 0;
				for (const auto &graph_instance : graph_list)
				{
					graph_list_size += graph_instance.computeSize();
				}
				outFile << "graph list size is " << graph_list_size << std::endl;
				outFile << "========================baseline2=========================" << std::endl;
				timer_baseline2.writeStatsToFile(outFile);
				graph_time.record_all_details_stream(outFile);
				outFile.close();
				std::cout <<"all finished" <<std::endl;
			}
			else if (config.hop_limit == 0)
			{
				experiment::nonhop::two_hop_case_info hop_info;
				std::string graph_res_filename = "binary_nonhop_constrained_" + std::to_string(config.hop_limit) + "_graph";
				std::string dataSource = config.data_source.string() + "//" + graph_res_filename;

				std::string hop_label_res_filename = "binary_nonhop_constrained_" + std::to_string(config.hop_limit) + "_2_hop_label_info";
				std::string experiment_MAINTAIN_LABEL_res_filename = "MAINTAIN_LABEL_nonhop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::filesystem::path hopLabelPath = saveDir.string() + "//" + hop_label_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_MAINTAIN_LABEL_res_filename;

				std::ifstream FILE_GRAPH(dataSource, std::ios::in | std::ifstream::binary);
				experiment::loadBinary(FILE_GRAPH, init_graph);
				experiment::loadBinary(FILE_GRAPH, graph_time);
				experiment::loadBinary(FILE_GRAPH, hop_info);
				hop_info.thread_num = config.threads;
				experiment::nonhop::two_hop_case_info hop_info_2021;
				hop_info_2021 = hop_info;
				std::cout << "L size is " << hop_info.compute_L_size() << std::endl;
				graph_list.push_back(init_graph);
				timer_ruc.endSubtask();
				timer_2021.endSubtask();
				experiment::iteration_info<int> change_info(
					graph_time.v_num, config.iterations, config.change_count, config.max_value, config.min_value, init_graph);
				change_info.build_random_change();

				std::vector<std::pair<int, int>> path_decrease;
				std::map<std::pair<int, int>, int> path2Index4Decrease;
				std::vector<int> weight_decrease;

				std::vector<std::pair<int, int>> path_increase;
				std::vector<int> weight_increase;
				std::vector<int> weight_old_increase;
				std::map<std::pair<int, int>, int> path2Index4Increase;

				ThreadPool pool_dynamic(hop_info.thread_num);
				std::vector<std::future<int>> results_dynamic;
				timer_ruc.startSubtask("step-2 maintain 2 hop label 0 - " + std::to_string(config.iterations / 2));
				timer_2021.startSubtask("step-2 maintain 2 hop label 0 - " + std::to_string(config.iterations / 2));
				timer_baseline1.startSubtask("step-1 save graph from 0 - " + std::to_string(config.iterations / 2));
				timer_baseline2.startSubtask("step-1 save graph with time span edge from 0 - " + std::to_string(config.iterations / 2));
				for (int i = 1; i <= config.iterations; i++)
				{
					if (i == config.iterations / 2 + 1)
					{
						timer_ruc.startSubtask("step-3 maintain 2 hop label " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_2021.startSubtask("step-3 maintain 2 hop label " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_baseline1.startSubtask("step-2 save graph from " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_baseline2.startSubtask("step-2 save graph with time span edge from " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
					}
					std::cout << "iteration " << i << std::endl;
					std::queue<experiment::change_edge_info> q = change_info.q_list[i];
					experiment::graph<int> instance_graph_temp = graph_list[i - 1];
					while (!q.empty())
					{
						experiment::change_edge_info info = q.front();
						q.pop();
						// 1. ��ȡ����
						int v1 = info.v1;
						int v2 = info.v2;
						int weight = info.weight;
						if (instance_graph_temp.ADJs[v1][v2].second < weight)
						{
							auto pairPathV = std::make_pair(v1, v2);
							auto it = path2Index4Increase.find(pairPathV);
							// increase
							if (it == path2Index4Increase.end())
							{
								int old_weight = instance_graph_temp.ADJs[v1][v2].second;
								path_increase.push_back(pairPathV);
								weight_old_increase.push_back(old_weight);
								weight_increase.push_back(weight);
								path2Index4Increase[pairPathV] = weight_increase.size() - 1;
							}
							else
							{
								weight_increase[path2Index4Increase[pairPathV]] = weight;
							}
						}
						else if (instance_graph_temp.ADJs[v1][v2].second > weight)
						{
							auto pairPathV = std::make_pair(v1, v2);
							auto it = path2Index4Decrease.find(pairPathV);
							if (it == path2Index4Decrease.end())
							{
								path_decrease.push_back(pairPathV);
								weight_decrease.push_back(weight);
								path2Index4Decrease[pairPathV] = weight_decrease.size() - 1;
							}
							else
							{
								weight_decrease[path2Index4Decrease[pairPathV]] = weight;
							}
						}
						if (path_decrease.size() > hop_info.thread_num)
						{
							for (int i = 0; i < path_decrease.size(); i++)
							{
								int v1 = path_decrease[i].first;
								int v2 = path_decrease[i].second;
								int w = weight_decrease[i];
								// std::cout << "decrease old weight is " << instance_graph_temp[v1][v2].second << " new weight is " << w << std::endl;
								instance_graph_temp[v1][v2].second = w;
							}
							std::cout << "decrease ruc maintain" << std::endl;
							experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
							timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc decrease maintain");
							experiment::nonhop::ruc::decrease::decrease_maintain(instance_graph_temp, hop_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
							timer_ruc.endSubtask();
							std::cout << "decrease 2021 maintain" << std::endl;
							experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
							timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 decrease maintain");
							experiment::nonhop::algorithm2021::decrease::decrease_maintain(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
							timer_2021.endSubtask();
							std::vector<std::pair<int, int>>().swap(path_decrease);
							std::vector<int>().swap(weight_decrease);
							std::map<std::pair<int, int>, int>().swap(path2Index4Decrease);
						}
						if (path_increase.size() > hop_info.thread_num)
						{
							for (int i = 0; i < path_increase.size(); i++)
							{
								int v1 = path_increase[i].first;
								int v2 = path_increase[i].second;
								int w = weight_increase[i];
								int w_old = weight_old_increase[i];
								// std::cout << "increase old weight is " << instance_graph_temp[v1][v2].second << " new weight is " << w << std::endl;
								instance_graph_temp[v1][v2].second = w;
							}
							std::cout << "increase ruc maintain" << std::endl;
							experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
							timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc increase maintain");
							experiment::nonhop::ruc::increase::nonHOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
							timer_ruc.endSubtask();
							std::cout << "increase 2021 maintain" << std::endl;
							experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
							timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 increase maintain");
							experiment::nonhop::algorithm2021::increase::increase_maintain(instance_graph_temp, hop_info_2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
							timer_2021.endSubtask();
							std::vector<std::pair<int, int>>().swap(path_increase);
							std::vector<int>().swap(weight_increase);
							std::vector<int>().swap(weight_old_increase);
							std::map<std::pair<int, int>, int>().swap(path2Index4Increase);
						}
					}
					if (path_decrease.size() > 0)
					{
						for (int i = 0; i < path_decrease.size(); i++)
						{
							int v1 = path_decrease[i].first;
							int v2 = path_decrease[i].second;
							int w = weight_decrease[i];
							// std::cout << "decrease old weight is " << instance_graph_temp[v1][v2].second << " new weight is " << w << std::endl;
							instance_graph_temp[v1][v2].second = w;
						}
						std::cout << "decrease ruc maintain" << std::endl;
						experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
						timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc decrease maintain");
						experiment::nonhop::ruc::decrease::decrease_maintain(instance_graph_temp, hop_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
						timer_ruc.endSubtask();
						std::cout << "decrease 2021 maintain" << std::endl;
						experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
						timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 decrease maintain");
						experiment::nonhop::algorithm2021::decrease::decrease_maintain(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
						timer_2021.endSubtask();
						std::vector<std::pair<int, int>>().swap(path_decrease);
						std::vector<int>().swap(weight_decrease);
						std::map<std::pair<int, int>, int>().swap(path2Index4Decrease);
					}
					if (path_increase.size() > 0)
					{
						for (int i = 0; i < path_increase.size(); i++)
						{
							int v1 = path_increase[i].first;
							int v2 = path_increase[i].second;
							int w = weight_increase[i];
							int w_old = weight_old_increase[i];
							// std::cout << "increase old weight is " << instance_graph_temp[v1][v2].second << " new weight is " << w << std::endl;
							instance_graph_temp[v1][v2].second = w;
						}
						std::cout << "increase ruc maintain" << std::endl;
						experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
						timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc increase maintain");
						experiment::nonhop::ruc::increase::nonHOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
						timer_ruc.endSubtask();
						std::cout << "increase 2021 maintain" << std::endl;
						experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
						timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 increase maintain");
						experiment::nonhop::algorithm2021::increase::increase_maintain(instance_graph_temp, hop_info_2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
						timer_2021.endSubtask();
						std::vector<std::pair<int, int>>().swap(path_increase);
						std::vector<int>().swap(weight_increase);
						std::vector<int>().swap(weight_old_increase);
					}
					timer_baseline1.startSubtask("save graph " + std::to_string(i));
					graph_list.push_back(instance_graph_temp);
					timer_baseline1.endSubtask();
					timer_baseline2.startSubtask("save graph with time span label " + std::to_string(i));
					graph_time.add_graph_time(instance_graph_temp, i);
					timer_baseline2.endSubtask();
					if (i == config.iterations / 2 || i == config.iterations)
					{
						timer_ruc.endSubtask();
						timer_2021.endSubtask();
						timer_baseline1.endSubtask();
						timer_baseline2.endSubtask();
					}
				}
				std::ofstream FILE_HOP_LABEL(hopLabelPath.string(), std::ios::out | std::ofstream::binary);
				experiment::saveBinary(FILE_HOP_LABEL, graph_list);
				experiment::saveBinary(FILE_HOP_LABEL, graph_time);
				experiment::saveBinary(FILE_HOP_LABEL, hop_info);
				experiment::saveBinary(FILE_HOP_LABEL, hop_info_2021);
				std::ofstream outFile;
				outFile.precision(6);
				outFile.setf(std::ios::fixed);
				outFile.setf(std::ios::showpoint);
				outFile.open(resultPath.string());
				outFile << "========================ruc maintain======================" << std::endl;
				timer_ruc.writeStatsToFile(outFile);
				hop_info.record_all_details_stream(outFile);
				outFile << "========================2021 maintain=====================" << std::endl;
				timer_2021.writeStatsToFile(outFile);
				hop_info_2021.record_all_details_stream(outFile);
				outFile << "========================baseline1=========================" << std::endl;
				timer_baseline1.writeStatsToFile(outFile);
				long long int graph_list_size = 0;
				for (const auto &graph_instance : graph_list)
				{
					graph_list_size += graph_instance.computeSize();
				}
				outFile << "graph list size is " << graph_list_size << std::endl;
				outFile << "========================baseline2=========================" << std::endl;
				timer_baseline2.writeStatsToFile(outFile);
				graph_time.record_all_details_stream(outFile);
				FILE_HOP_LABEL.close();
				outFile.close();
			}
		}
		else if (config.mode == experiment::QUERY_RESULT)
		{
			// random src and dest
			timer.startTask("query shorest path distance");
			std::vector<experiment::graph<int>> graph_list;
			experiment::graph_with_time_span<int> graph_time;
			if (config.hop_limit != 0)
			{
				std::string experiment_QUERY_RESULT_res_filename = "QUERY_RESULT_nonhop_constrained_" + std::to_string(config.hop_limit) + "_result.txt";
				std::string data_from_filename = "binary_hop_constrained_" + std::to_string(config.hop_limit) + "_2_hop_label_info";
				std::string dataSource = config.data_source.string() + "//" + data_from_filename;
				std::string savePath = config.data_source.string() + "//" + experiment_QUERY_RESULT_res_filename;
				std::ifstream FILE_GRAPH(dataSource, std::ios::in | std::ifstream::binary);

				experiment::hop::two_hop_case_info hop_info;
				experiment::hop::two_hop_case_info hop_info_2021;
				experiment::loadBinary(FILE_GRAPH, graph_list);
				experiment::loadBinary(FILE_GRAPH, graph_time);
				experiment::loadBinary(FILE_GRAPH, hop_info);
				experiment::loadBinary(FILE_GRAPH, hop_info_2021);
				int v_num = graph_time.v_num;
				int time = graph_time.time_max;
				int hop = config.hop_limit;
				boost::random::uniform_int_distribution<> _random_v = boost::random::uniform_int_distribution<>(0, v_num);
				boost::random::uniform_int_distribution<> _random_time = boost::random::uniform_int_distribution<>(0, time);
				boost::random::uniform_int_distribution<> _random_hop = boost::random::uniform_int_distribution<>(0, hop);
				std::ofstream outFile;
				outFile.precision(6);
				outFile.setf(std::ios::fixed);
				outFile.setf(std::ios::showpoint);
				outFile.open(savePath);
				for (int i = 0; i < config.change_count; i++)
				{
					int index_i = _random_v(boost_random_time_seed);
					int index_j = _random_v(boost_random_time_seed);
					int t_1 = _random_time(boost_random_time_seed);
					int t_2 = _random_time(boost_random_time_seed);
					int hop = _random_hop(boost_random_time_seed);
					if (t_1 > t_2)
					{
						std::swap(t_1, t_2);
					}
					timer.startSubtask("====iteration " + std::to_string(i) + " query result info====");
					timer.startSubtask("baseline 1: traverse each time graph");
					int resb1 = experiment::hop::dijkstra_iterator(graph_list, index_i, index_j, t_1, t_2, hop);
					timer.endSubtask();
					timer.startSubtask("baseline 2: traverse graph with time span");
					int resb2 = experiment::hop::search_shortest_path_in_period_time_naive(graph_time, index_i, index_j, t_1, t_2, hop);
					timer.endSubtask();
					timer.startSubtask("search result by ruc maintain algorithm");
					int ruc_res = hop_info.query(index_i, index_j, t_1, t_2, hop);
					timer.endSubtask();
					timer.startSubtask("search result by 2021 maintain algorithm");
					int res_2021 = hop_info_2021.query(index_i, index_j, t_1, t_2, hop);
					timer.endSubtask();
					outFile << resb1 << ":" << resb2 << ":" << ruc_res << ":" << res_2021 << std::endl;
					timer.endSubtask();
				}
				timer.writeStatsToFile(outFile);
				outFile.close();
			}
			else
			{
				std::string experiment_QUERY_RESULT_res_filename = "QUERY_RESULT_nonhop_constrained_" + std::to_string(config.hop_limit) + "_result.txt";
				std::string data_from_filename = "binary_nonhop_constrained_" + std::to_string(config.hop_limit) + "_2_hop_label_info";
				std::string dataSource = config.data_source.string() + "//" + data_from_filename;
				std::string savePath = config.data_source.string() + "//" + experiment_QUERY_RESULT_res_filename;
				std::cout << "save path is " << savePath << std::endl;
				std::ifstream FILE_GRAPH(dataSource, std::ios::in | std::ifstream::binary);
				experiment::nonhop::two_hop_case_info hop_info;
				experiment::nonhop::two_hop_case_info hop_info_2021;
				experiment::loadBinary(FILE_GRAPH, graph_list);
				experiment::loadBinary(FILE_GRAPH, graph_time);
				experiment::loadBinary(FILE_GRAPH, hop_info);
				experiment::loadBinary(FILE_GRAPH, hop_info_2021);
				int v_num = graph_time.v_num;
				int time = graph_time.time_max;
				boost::random::uniform_int_distribution<> _random_v = boost::random::uniform_int_distribution<>(0, v_num);
				boost::random::uniform_int_distribution<> _random_time = boost::random::uniform_int_distribution<>(0, time);
				std::ofstream outFile;
				outFile.precision(6);
				outFile.setf(std::ios::fixed);
				outFile.setf(std::ios::showpoint);
				outFile.open(savePath);
				for (int i = 0; i < config.change_count; i++)
				{
					int index_i = _random_v(boost_random_time_seed);
					int index_j = _random_v(boost_random_time_seed);
					int t_1 = _random_time(boost_random_time_seed);
					int t_2 = _random_time(boost_random_time_seed);
					if (t_1 > t_2)
					{
						std::swap(t_1, t_2);
					}
					timer.startSubtask("====iteration " + std::to_string(i) + " query result info====");
					timer.startSubtask("baseline 1: traverse each time graph");
					int resb1 = experiment::nonhop::dijkstra_iterator(graph_list, index_i, index_j, t_1, t_2);
					timer.endSubtask();
					timer.startSubtask("baseline 2: traverse graph with time span");
					int resb2 = experiment::nonhop::search_shortest_path_in_period_time_naive(graph_time, index_i, index_j, t_1, t_2);
					timer.endSubtask();
					timer.startSubtask("search result by ruc maintain algorithm");
					int ruc_res = hop_info.query(index_i, index_j, t_1, t_2);
					timer.endSubtask();
					timer.startSubtask("search result by 2021 maintain algorithm");
					int res_2021 = hop_info_2021.query(index_i, index_j, t_1, t_2);
					timer.endSubtask();
					outFile << resb1 << ":" << resb2 << ":" << ruc_res << ":" << res_2021 << std::endl;
					timer.endSubtask();
				}
				timer.writeStatsToFile(outFile);
				outFile.close();
			}
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Error: " << ex.what() << "\n";
		exit(EXIT_FAILURE);
	}
	std::cout << "finished" << std::endl;
	exit(EXIT_SUCCESS);
}