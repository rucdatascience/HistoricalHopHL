#include <iostream>
#include <string>
#include <filesystem>
#include <Historical/utils/ExecutionTimer.h>
experiment::ExecutionTimer timer;
experiment::ExecutionTimer timer_ruc;
experiment::ExecutionTimer timer_2021;
experiment::ExecutionTimer timer_baseline1;
experiment::ExecutionTimer timer_baseline2;
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
			instance_graph.graph_v_of_v_update_vertexIDs_by_degrees_large_to_small();
			std::cout << instance_graph.ADJs[0].size() << std::endl;
			std::cout << instance_graph.ADJs[1].size() << std::endl;
			timer.endSubtask();
			experiment::graph_with_time_span graph_time;
			graph_time.add_graph_time(instance_graph, 0);
			std::filesystem::path saveDir = std::filesystem::path(config.save_path);
			std::filesystem::create_directory(saveDir);
			if (config.hop_limit != 0)
			{
				// k-constrained pll
				std::string graph_res_filename = "binary_hop_constrained_" + std::to_string(config.hop_limit) + "_graph";
				std::string experiment_res_filename = "GENERATE_LABEL__hop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result_new.txt";
				std::filesystem::path graphPath = saveDir.string() + "//" + graph_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_res_filename;
				experiment::hop::two_hop_case_info hop_info;
				hop_info.thread_num = config.threads;
				hop_info.upper_k = config.hop_limit;
				timer.startSubtask("generate 2hop label " + std::to_string(config.hop_limit) + " hop constrained");
				experiment::hop::pll(instance_graph, hop_info);
				std::cout << "1" << std::endl;
				timer.endSubtask();
				std::cout << "finish pll " << std::endl;
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
				outFile << "pre L size is " << experiment::hop::globalLabelSize << " clean L size is " << experiment::hop::globalLabelCleanSize << " ppr size is " << experiment::hop::globalPprSize;
				outFile.close();
			}
			else
			{
				experiment::nonhop::two_hop_case_info hop_info;
				experiment::nonhop::two_hop_case_info hop_info_new;
				std::string graph_res_filename = "binary_nonhop_constrained_" + std::to_string(config.hop_limit) + "_graph";
				std::string experiment_res_filename = "GENERATE_LABEL_nonhop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result_old.txt";
				std::filesystem::path graphPath = saveDir.string() + "//" + graph_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_res_filename;
				std::filesystem::path resultGraphPath = saveDir.string() + "//" + "test_save_graph";

				hop_info.thread_num = config.threads;
				timer.startSubtask("generate graph and 2hop label " + std::to_string(config.hop_limit) + " nonhop constrained");
				experiment::nonhop::pll(instance_graph, hop_info);
				timer.endSubtask();
				// hop_info.print_L();
				std::ofstream FILE_GRAPH(graphPath.string(), std::ios::out | std::ofstream::binary | std::ios::app);
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
				outFile << "pre L size is " << experiment::nonhop::globalLabelSize << " clean L size is " << experiment::nonhop::globalLabelCleanSize << " ppr size is " << experiment::nonhop::globalPprSize;
				outFile.close();
			}
		}
		else if (config.mode == experiment::MAINTAIN_LABEL)
		{
			std::vector<experiment::graph<int>> graph_list;
			experiment::graph<int> init_graph;
			experiment::graph_with_time_span graph_time;
			std::filesystem::path saveDir = std::filesystem::path(config.save_path);
			std::filesystem::create_directory(saveDir);
			timer_ruc.startTask("maintain graph and 2hop label " + std::to_string(config.hop_limit) + (config.hop_limit == 0 ? "nonhop_constrained" : "hop constrained"));
			timer_2021.startTask("maintain graph and 2hop label " + std::to_string(config.hop_limit) + (config.hop_limit == 0 ? "nonhop_constrained" : "hop constrained"));
			timer_baseline1.startTask("maintain graph by saving every graph");
			timer_baseline2.startTask("maintain graph by saving edge info with time label");
			timer_ruc.startSubtask("step-1 read graph and original 2hop label");
			timer_2021.startSubtask("step-1 read graph and original 2hop label");
			long long int half_ruc_size = 0;
			long long int half_2021_size = 0;
			long long int half_baseline1_size = 0;
			long long int half_baseline2_size = 0;
			if (config.hop_limit != 0)
			{
				double rucTimeCostAll = 0;
				double a2021TimeCostAll = 0;
				double base1TimeCostAll = 0;
				double base2TimeCostAll = 0;
				double rucTimeCostHalf = 0;
				double a2021TimeCostHalf = 0;
				double base1TimeCostHalf = 0;
				double base2TimeCostHalf = 0;
				experiment::hop::two_hop_case_info hop_info;
				std::string graph_res_filename = "binary_hop_constrained_" + std::to_string(config.hop_limit) + "_graph";
				std::string dataSource = config.data_source.string() + "//" + graph_res_filename;

				std::string hop_label_res_filename = "binary_hop_constrained_" + std::to_string(config.hop_limit) + "_2_hop_label_info";
				std::string experiment_res_filename = "MAINTAIN_LABEL_hop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::string change_info_res_filename = "change_info_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::string change_info_res_filename_detail = "change_info_detail" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::filesystem::path hopLabelPath = saveDir.string() + "//" + hop_label_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_res_filename;
				std::filesystem::path changePath = saveDir.string() + "//" + change_info_res_filename;
				std::filesystem::path changePathDetail = saveDir.string() + "//" + change_info_res_filename_detail;

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
				std::ofstream CHANGE_DETAIL_PATH_STREAM(changePathDetail.string(), std::ios::out | std::ofstream::binary);
				experiment::saveBinary(CHANGE_PATH_STREAM, change_info);
				CHANGE_PATH_STREAM.close();
				change_info.toString(CHANGE_DETAIL_PATH_STREAM);

				// experiment::iteration_info<int> change_info;
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
						half_ruc_size = hop_info.compute_label_bit_size();
						half_2021_size = hop_info_2021.compute_label_bit_size();
						for (const auto &graph_instance : graph_list)
						{
							half_baseline1_size += graph_instance.computeSize();
						};
						half_baseline2_size = graph_time.computeSize();
						timer_ruc.startSubtask("step-3 maintain 2 hop label " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_2021.startSubtask("step-3 maintain 2 hop label " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_baseline1.startSubtask("step-2 save graph from " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_baseline2.startSubtask("step-2 save graph with time span edge from " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
					}
					timer_ruc.startSubtask("start " + std::to_string(i) + " iteration");
					timer_2021.startSubtask("start " + std::to_string(i) + " iteration");
					timer_baseline1.startSubtask("start " + std::to_string(i) + " iteration");
					timer_baseline2.startSubtask("start " + std::to_string(i) + " iteration");
					std::cout << "iteration " << i << std::endl;
					std::queue<experiment::change_edge_info> q = change_info.q_list[i];
					experiment::graph<int> instance_graph_temp = graph_list[i - 1];
					while (!q.empty())
					{
						experiment::change_edge_info info = q.front();
						q.pop();
						int v1 = info.v1;
						int v2 = info.v2;
						int weight = info.weight;
						int old_weight = sorted_vector_binary_operations_search_weight(instance_graph_temp.ADJs[v1], v2);
						if (old_weight < weight)
						{
							auto pairPathV = std::make_pair(v1, v2);
							auto it = path2Index4Increase.find(pairPathV);
							// increase
							if (it == path2Index4Increase.end())
							{
								path_increase.push_back(pairPathV);
								weight_increase.push_back(weight);
								path2Index4Increase[pairPathV] = weight_increase.size() - 1;
							}
							else
							{
								weight_increase[path2Index4Increase[pairPathV]] = weight;
							}
						}
						else if (old_weight > weight)
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
							for (int index = 0; index < path_decrease.size(); index++)
							{
								int v1 = path_decrease[index].first;
								int v2 = path_decrease[index].second;
								int w = weight_decrease[index];
								int old_w = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
								std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << old_w << std::endl;
								timer_baseline1.startSubtask("modify baseline1 edge weight");
								instance_graph_temp.add_edge(v1, v2, w);
								base1TimeCostAll+= timer_baseline1.endSubtask();
								timer_baseline2.startSubtask("modify baseline2 edge weight");
								graph_time.add_edge(v1, v2, w, i);
								base2TimeCostAll+= timer_baseline2.endSubtask();
							}
							std::cout << "decrease ruc maintain" << std::endl;
							experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
							timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc decrease maintain");
							experiment::hop::ruc::decrease::HOP_WeightDecreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
							rucTimeCostAll += timer_ruc.endSubtask();
							std::cout << "decrease 2021 maintain" << std::endl;
							experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
							timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 decrease maintain");
							experiment::hop::algorithm2021::decrease::HOP_WeightDecrease2021_batch(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
							a2021TimeCostAll += timer_2021.endSubtask();
							std::vector<std::pair<int, int>>().swap(path_decrease);
							std::vector<int>().swap(weight_decrease);
							std::map<std::pair<int, int>, int>().swap(path2Index4Decrease);
						}
						if (path_increase.size() > hop_info.thread_num)
						{
							for (int index = 0; index < path_increase.size(); index++)
							{
								int v1 = path_increase[index].first;
								int v2 = path_increase[index].second;
								int w = weight_increase[index];
								int old_w = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
								weight_old_increase.push_back(old_w);
								std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << old_w << std::endl;
								timer_baseline1.startSubtask("modify baseline1 edge weight");
								instance_graph_temp.add_edge(v1, v2, w);
								base1TimeCostAll += timer_baseline1.endSubtask();
								timer_baseline2.startSubtask("modify baseline2 edge weight");
								graph_time.add_edge(v1, v2, w, i);
								base2TimeCostAll+= timer_baseline2.endSubtask();
							}
							std::cout << "increase 2021 maintain" << std::endl;
							experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
							timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 increase maintain");
							experiment::hop::algorithm2021::increase::HOP_WeightIncrease2021_batch(instance_graph_temp, hop_info_2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
							a2021TimeCostAll += timer_2021.endSubtask();
							std::cout << "increase ruc maintain" << std::endl;
							experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
							timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc increase maintain");
							experiment::hop::ruc::increase::HOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
							rucTimeCostAll += timer_ruc.endSubtask();
							std::vector<std::pair<int, int>>().swap(path_increase);
							std::vector<int>().swap(weight_increase);
							std::vector<int>().swap(weight_old_increase);
							std::map<std::pair<int, int>, int>().swap(path2Index4Increase);
						}
					}
					if (path_decrease.size() > 0)
					{
						for (int index = 0; index < path_decrease.size(); index++)
						{
							int v1 = path_decrease[index].first;
							int v2 = path_decrease[index].second;
							int w = weight_decrease[index];
							int old_w = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
							std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << old_w << std::endl;
							timer_baseline1.startSubtask("modify baseline1 edge weight");
							instance_graph_temp.add_edge(v1, v2, w);
							base1TimeCostAll += timer_baseline1.endSubtask();
							timer_baseline2.startSubtask("modify baseline2 edge weight");
							graph_time.add_edge(v1, v2, w, i);
							base2TimeCostAll += timer_baseline2.endSubtask();
						}
						std::cout << "decrease ruc maintain" << std::endl;
						experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
						timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc decrease maintain");
						experiment::hop::ruc::decrease::HOP_WeightDecreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
						rucTimeCostAll += timer_ruc.endSubtask();
						std::cout << "decrease 2021 maintain" << std::endl;
						experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
						timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 decrease maintain");
						experiment::hop::algorithm2021::decrease::HOP_WeightDecrease2021_batch(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
						a2021TimeCostAll += timer_2021.endSubtask();
						std::vector<std::pair<int, int>>().swap(path_decrease);
						std::vector<int>().swap(weight_decrease);
						std::map<std::pair<int, int>, int>().swap(path2Index4Decrease);
					}
					if (path_increase.size() > 0)
					{
						for (int index = 0; index < path_increase.size(); index++)
						{
							int v1 = path_increase[index].first;
							int v2 = path_increase[index].second;
							int w = weight_increase[index];
							int old_w = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
							weight_old_increase.push_back(old_w);
							std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << old_w << std::endl;
							timer_baseline1.startSubtask("modify baseline1 edge weight");
							instance_graph_temp.add_edge(v1, v2, w);
							base1TimeCostAll += timer_baseline1.endSubtask();
							timer_baseline2.startSubtask("modify baseline2 edge weight");
							graph_time.add_edge(v1, v2, w, i);
							base2TimeCostAll += timer_baseline2.endSubtask();
						}
						std::cout << "increase 2021 maintain" << std::endl;
						experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
						timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 increase maintain");
						experiment::hop::algorithm2021::increase::HOP_WeightIncrease2021_batch(instance_graph_temp, hop_info_2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
						a2021TimeCostAll += timer_2021.endSubtask();
						std::cout << "increase ruc maintain" << std::endl;
						experiment::hop::initialize_global_values_dynamic_hop_constrained(instance_graph_temp.size(), hop_info.thread_num, hop_info.upper_k);
						timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc increase maintain");
						experiment::hop::ruc::increase::HOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
						rucTimeCostAll += timer_ruc.endSubtask();
						std::vector<std::pair<int, int>>().swap(path_increase);
						std::vector<int>().swap(weight_increase);
						std::vector<int>().swap(weight_old_increase);
						std::map<std::pair<int, int>, int>().swap(path2Index4Increase);
					}
					graph_list.push_back(instance_graph_temp);
					timer_ruc.endSubtask();
					timer_2021.endSubtask();
					timer_baseline1.endSubtask();
					timer_baseline2.endSubtask();
					if (i == config.iterations / 2 || i == config.iterations)
					{
						if(i==(config.iterations / 2)){
							rucTimeCostHalf = rucTimeCostAll;
							a2021TimeCostHalf = a2021TimeCostAll;
							base1TimeCostHalf = base1TimeCostAll;
							base2TimeCostHalf = base2TimeCostAll;
						}
						timer_ruc.endSubtask();
						timer_2021.endSubtask();
						timer_baseline1.endSubtask();
						timer_baseline2.endSubtask();
					}
					// std::cout <<" time graph size is "<< graph_time.computeSize() <<" when time is "<< i << std::endl;
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
				outFile << "========================2021 maintain=====================" << std::endl;
				timer_2021.writeStatsToFile(outFile);
				outFile << "========================baseline1=========================" << std::endl;
				timer_baseline1.writeStatsToFile(outFile);
				long long int graph_list_size = 0;
				for (const auto &graph_instance : graph_list)
				{
					graph_list_size += graph_instance.computeSize();
				}
				outFile << "========================baseline2=========================" << std::endl;
				timer_baseline2.writeStatsToFile(outFile);
				outFile << "========================solt1=========================" << std::endl;
				outFile << "baseline1 mem cost " << half_baseline1_size << std::endl;
				outFile << "baseline2 mem cost " << half_baseline2_size << std::endl;
				outFile << "2021 mem cost " << half_2021_size << std::endl;
				outFile << "ruc1 mem cost " << half_ruc_size << std::endl;
				outFile << "========================solt2=========================" << std::endl;
				outFile << "baseline1 mem cost " << graph_list_size << std::endl;
				outFile << "baseline2 mem cost " << graph_time.computeSize() << std::endl;
				outFile << "2021 mem cost " << hop_info_2021.compute_label_bit_size() << std::endl;
				outFile << "ruc1 mem cost " << hop_info.compute_label_bit_size() << std::endl;
				outFile << "========================solt1===============================" << std::endl;
				outFile << "baseline1 cost " << base1TimeCostHalf << std::endl;
				outFile << "baseline2 cost " << base2TimeCostHalf << std::endl;
				outFile << "2021 cost " << a2021TimeCostHalf << std::endl;
				outFile << "ruc1 cost " << rucTimeCostHalf << std::endl;
				outFile << "========================solt2===============================" << std::endl;
				outFile << "baseline1 cost " << base1TimeCostAll-base1TimeCostHalf << std::endl;
				outFile << "baseline2 cost " << base2TimeCostAll-base2TimeCostHalf << std::endl;
				outFile << "2021 cost " << a2021TimeCostAll-a2021TimeCostHalf << std::endl;
				outFile << "ruc1 cost " << rucTimeCostAll-rucTimeCostHalf << std::endl;
				outFile.close();
				std::cout << "all finished" << std::endl;
			}
			else if (config.hop_limit == 0)
			{
				double rucTimeCostAll = 0;
				double a2021TimeCostAll = 0;
				double base1TimeCostAll = 0;
				double base2TimeCostAll = 0;
				double rucTimeCostHalf = 0;
				double a2021TimeCostHalf = 0;
				double base1TimeCostHalf = 0;
				double base2TimeCostHalf = 0;
				experiment::nonhop::two_hop_case_info hop_info;
				std::string graph_res_filename = "binary_nonhop_constrained_" + std::to_string(config.hop_limit) + "_graph";
				std::string dataSource = config.data_source.string() + "//" + graph_res_filename;

				std::string hop_label_res_filename = "binary_nonhop_constrained_" + std::to_string(config.hop_limit) + "_2_hop_label_info";
				std::string experiment_MAINTAIN_LABEL_res_filename = "MAINTAIN_LABEL_nonhop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::string change_info_res_filename = "change_info_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::string change_info_res_filename_detail = "change_info_detail" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				
				std::filesystem::path hopLabelPath = saveDir.string() + "//" + hop_label_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_MAINTAIN_LABEL_res_filename;
				std::filesystem::path resultGraphPrePath = saveDir.string() + "//" + "test_read_graph";
				std::filesystem::path changePath = saveDir.string() + "//" + change_info_res_filename;
				std::filesystem::path changePathDetail = saveDir.string() + "//" + change_info_res_filename_detail;

				long long int half_ruc_size = 0;
				long long int half_2021_size = 0;
				long long int half_baseline1_size = 0;
				long long int half_baseline2_size = 0;

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
				std::ofstream CHANGE_PATH_STREAM(changePath.string(), std::ios::out | std::ofstream::binary);
				std::ofstream CHANGE_DETAIL_PATH_STREAM(changePathDetail.string(), std::ios::out | std::ofstream::binary);
				experiment::saveBinary(CHANGE_PATH_STREAM, change_info);
				CHANGE_PATH_STREAM.close();
				change_info.toString(CHANGE_DETAIL_PATH_STREAM);

				// experiment::iteration_info<int> change_info;
				// std::ifstream CHANGE_PATH_STREAM(changePath.string(), std::ios::in | std::ofstream::binary);
				// experiment::loadBinary(CHANGE_PATH_STREAM, change_info);
				// CHANGE_PATH_STREAM.close();

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
						half_ruc_size = hop_info.compute_L_byte_size();
						half_2021_size = hop_info_2021.compute_L_byte_size();
						for (const auto &graph_instance : graph_list)
						{
							half_baseline1_size += graph_instance.computeSize();
						};
						half_baseline2_size = graph_time.computeSize();
						timer_ruc.startSubtask("step-3 maintain 2 hop label " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_2021.startSubtask("step-3 maintain 2 hop label " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_baseline1.startSubtask("step-2 save graph from " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
						timer_baseline2.startSubtask("step-2 save graph with time span edge from " + std::to_string(config.iterations / 2 + 1) + " - " + std::to_string(config.iterations));
					}
					timer_ruc.startSubtask("start " + std::to_string(i) + " iteration");
					timer_2021.startSubtask("start " + std::to_string(i) + " iteration");
					timer_baseline1.startSubtask("start " + std::to_string(i) + " iteration");
					timer_baseline2.startSubtask("start " + std::to_string(i) + " iteration");

					std::cout << "iteration " << i << std::endl;
					std::queue<experiment::change_edge_info> q = change_info.q_list[i];
					experiment::graph<int> instance_graph_temp = graph_list[i - 1];
					while (!q.empty())
					{
						experiment::change_edge_info info = q.front();
						q.pop();
						int v1 = info.v1;
						int v2 = info.v2;
						int weight = info.weight;
						int old_weight = sorted_vector_binary_operations_search_weight(instance_graph_temp.ADJs[v1], v2);
						if (old_weight < weight)
						{
							auto pairPathV = std::make_pair(v1, v2);
							auto it = path2Index4Increase.find(pairPathV);
							// increase
							if (it == path2Index4Increase.end())
							{
								path_increase.push_back(pairPathV);
								weight_increase.push_back(weight);
								path2Index4Increase[pairPathV] = weight_increase.size() - 1;
							}
							else
							{
								weight_increase[path2Index4Increase[pairPathV]] = weight;
							}
						}
						else if (old_weight > weight)
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
						if (path_decrease.size() > hop_info.thread_num * 10)
						{
							for (int index = 0; index < path_decrease.size(); index++)
							{
								int v1 = path_decrease[index].first;
								int v2 = path_decrease[index].second;
								int w = weight_decrease[index];
								int w_old = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
								std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << w_old << std::endl;
								timer_baseline1.startSubtask("modify baseline1 edge weight");
								instance_graph_temp.add_edge(v1, v2, w);
								base1TimeCostAll += timer_baseline1.endSubtask();
								timer_baseline2.startSubtask("modify baseline2 edge weight");
								graph_time.add_edge(v1, v2, w, i);
								base2TimeCostAll+= timer_baseline2.endSubtask();
							}
							std::cout << "decrease ruc maintain" << std::endl;
							experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
							timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc decrease maintain");
							experiment::nonhop::ruc::decrease::decrease_maintain(instance_graph_temp, hop_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
							rucTimeCostAll+= timer_ruc.endSubtask();
							std::cout << "decrease 2021 maintain" << std::endl;
							experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
							timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 decrease maintain");
							experiment::nonhop::algorithm2021::decrease::decrease_maintain(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
							a2021TimeCostAll += timer_2021.endSubtask();
							std::vector<std::pair<int, int>>().swap(path_decrease);
							std::vector<int>().swap(weight_decrease);
							std::map<std::pair<int, int>, int>().swap(path2Index4Decrease);
						}
						if (path_increase.size() > hop_info.thread_num * 10)
						{
							for (int index = 0; index < path_increase.size(); index++)
							{
								int v1 = path_increase[index].first;
								int v2 = path_increase[index].second;
								int w = weight_increase[index];
								int w_old = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
								weight_old_increase.push_back(w_old);
								std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << w_old << std::endl;
								timer_baseline1.startSubtask("modify baseline1 edge weight");
								instance_graph_temp.add_edge(v1, v2, w);
								base1TimeCostAll += timer_baseline1.endSubtask();
								timer_baseline2.startSubtask("modify baseline2 edge weight");
								graph_time.add_edge(v1, v2, w, i);
								base2TimeCostAll += timer_baseline2.endSubtask();
							}
							std::cout << "increase ruc maintain" << std::endl;
							experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
							timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc increase maintain");
							experiment::nonhop::ruc::increase::nonHOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
							rucTimeCostAll += timer_ruc.endSubtask();
							std::cout << "increase 2021 maintain" << std::endl;
							experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
							timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 increase maintain");
							experiment::nonhop::algorithm2021::increase::increase_maintain(instance_graph_temp, hop_info_2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
							a2021TimeCostAll += timer_2021.endSubtask();
							std::vector<std::pair<int, int>>().swap(path_increase);
							std::vector<int>().swap(weight_increase);
							std::vector<int>().swap(weight_old_increase);
							std::map<std::pair<int, int>, int>().swap(path2Index4Increase);
						}
					}
					if (path_decrease.size() > 0)
					{
						for (int index = 0; index < path_decrease.size(); index++)
						{
							int v1 = path_decrease[index].first;
							int v2 = path_decrease[index].second;
							int w = weight_decrease[index];
							int w_old = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
							std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << w_old << std::endl;
							timer_baseline1.startSubtask("modify baseline1 edge weight");
							instance_graph_temp.add_edge(v1, v2, w);
							base1TimeCostAll += timer_baseline1.endSubtask();
							timer_baseline2.startSubtask("modify baseline2 edge weight");
							graph_time.add_edge(v1, v2, w, i);
							base2TimeCostAll += timer_baseline2.endSubtask();
						}
						std::cout << "decrease ruc maintain" << std::endl;
						experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
						timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc decrease maintain");
						auto start = std::chrono::steady_clock::now();
						experiment::nonhop::ruc::decrease::decrease_maintain(instance_graph_temp, hop_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
						auto end = std::chrono::steady_clock::now();
						auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
						std::cout <<"iteration " + std::to_string(i) + " algorithm ruc decrease maintain cost " << duration << std::endl;
						rucTimeCostAll += timer_ruc.endSubtask();
						std::cout << "decrease 2021 maintain" << std::endl;
						experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
						timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 decrease maintain");
						experiment::nonhop::algorithm2021::decrease::decrease_maintain(instance_graph_temp, hop_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, i);
						a2021TimeCostAll += timer_2021.endSubtask();
						std::vector<std::pair<int, int>>().swap(path_decrease);
						std::vector<int>().swap(weight_decrease);
						std::map<std::pair<int, int>, int>().swap(path2Index4Decrease);
					}
					if (path_increase.size() > 0)
					{
						for (int index = 0; index < path_increase.size(); index++)
						{
							int v1 = path_increase[index].first;
							int v2 = path_increase[index].second;
							int w = weight_increase[index];
							int w_old = sorted_vector_binary_operations_search_weight(instance_graph_temp[v1], v2);
							weight_old_increase.push_back(w_old);
							std::cout << "from " << v1 << " to " << v2 << " w " << w << " old_w is " << w_old << std::endl;
							timer_baseline1.startSubtask("modify baseline1 edge weight");
							instance_graph_temp.add_edge(v1, v2, w);
							base1TimeCostAll += timer_baseline1.endSubtask();
							timer_baseline2.startSubtask("modify baseline2 edge weight");
							graph_time.add_edge(v1, v2, w, i);
							base2TimeCostAll += timer_baseline2.endSubtask();
						}
						std::cout << "increase ruc maintain" << std::endl;
						experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
						timer_ruc.startSubtask("iteration " + std::to_string(i) + " algorithm ruc increase maintain");
						experiment::nonhop::ruc::increase::nonHOP_WeightIncreaseMaintenance_improv_batch(instance_graph_temp, hop_info, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
						rucTimeCostAll += timer_ruc.endSubtask();
						std::cout << "increase 2021 maintain" << std::endl;
						experiment::nonhop::initialize_experiment_global_values_dynamic(instance_graph_temp.size(), hop_info.thread_num);
						timer_2021.startSubtask("iteration " + std::to_string(i) + " algorithm 2021 increase maintain");
						experiment::nonhop::algorithm2021::increase::increase_maintain(instance_graph_temp, hop_info_2021, path_increase, weight_old_increase, pool_dynamic, results_dynamic, i);
						a2021TimeCostAll += timer_2021.endSubtask();
						std::vector<std::pair<int, int>>().swap(path_increase);
						std::vector<int>().swap(weight_increase);
						std::vector<int>().swap(weight_old_increase);
						std::map<std::pair<int, int>, int>().swap(path2Index4Increase);
					}
					graph_list.push_back(instance_graph_temp);
					timer_ruc.endSubtask();
					timer_2021.endSubtask();
					timer_baseline1.endSubtask();
					timer_baseline2.endSubtask();
					if (i == config.iterations / 2 || i == config.iterations)
					{
						if(i==(config.iterations / 2)){
							rucTimeCostHalf = rucTimeCostAll;
							a2021TimeCostHalf = a2021TimeCostAll;
							base1TimeCostHalf = base1TimeCostAll;
							base2TimeCostHalf = base2TimeCostAll;
						}
						timer_ruc.endSubtask();
						timer_2021.endSubtask();
						timer_baseline1.endSubtask();
						timer_baseline2.endSubtask();
					}
					// std::cout <<" time graph size is "<< graph_time.computeSize() <<" when time is "<< i << std::endl;
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
				outFile << "========================2021 maintain=====================" << std::endl;
				timer_2021.writeStatsToFile(outFile);
				outFile << "========================baseline1=========================" << std::endl;
				timer_baseline1.writeStatsToFile(outFile);
				long long int graph_list_size = 0;
				for (const auto &graph_instance : graph_list)
				{
					graph_list_size += graph_instance.computeSize();
				}
				outFile << "========================baseline2=========================" << std::endl;
				timer_baseline2.writeStatsToFile(outFile);
				outFile << "========================solt1=========================" << std::endl;
				outFile << "baseline1 mem cost " << half_baseline1_size << std::endl;
				outFile << "baseline2 mem cost " << half_baseline2_size << std::endl;
				outFile << "2021 mem cost " << half_2021_size << std::endl;
				outFile << "ruc1 mem cost " << half_ruc_size << std::endl;
				outFile << "========================solt2=========================" << std::endl;
				outFile << "baseline1 mem cost " << graph_list_size << std::endl;
				outFile << "baseline2 mem cost " << graph_time.computeSize() << std::endl;
				outFile << "2021 mem cost " << hop_info_2021.compute_L_byte_size() << std::endl;
				outFile << "ruc1 mem cost " << hop_info.compute_L_byte_size() << std::endl;
				outFile << "========================solt1===============================" << std::endl;
				outFile << "baseline1 cost " << base1TimeCostHalf << std::endl;
				outFile << "baseline2 cost " << base2TimeCostHalf << std::endl;
				outFile << "2021 cost " << a2021TimeCostHalf << std::endl;
				outFile << "ruc1 cost " << rucTimeCostHalf << std::endl;
				outFile << "========================solt2===============================" << std::endl;
				outFile << "baseline1 cost " << base1TimeCostAll-base1TimeCostHalf << std::endl;
				outFile << "baseline2 cost " << base2TimeCostAll-base2TimeCostHalf << std::endl;
				outFile << "2021 cost " << a2021TimeCostAll-a2021TimeCostHalf << std::endl;
				outFile << "ruc1 cost " << rucTimeCostAll-rucTimeCostHalf << std::endl;
				outFile.close();
			}
		}
		else if (config.mode == experiment::QUERY_RESULT)
		{
			// random src and dest
			std::vector<experiment::graph<int>> graph_list;
			experiment::graph_with_time_span graph_time;
			std::mutex outMutex;
			ThreadPool pool(config.threads);
			std::vector<std::future<int>> results; // return typename: xxx
			if (config.hop_limit != 0)
			{
				double rucTimeCostAll = 0;
				double a2021TimeCostAll = 0;
				double base1TimeCostAll = 0;
				double base2TimeCostAll = 0;
				std::string experiment_QUERY_RESULT_res_filename = "QUERY_RESULT_hop_constrained_" + std::to_string(config.hop_limit) + "_result.txt";
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
				boost::random::uniform_int_distribution<> _random_hop = boost::random::uniform_int_distribution<>(1, hop);
				std::ofstream outFile;
				outFile.precision(6);
				outFile.setf(std::ios::fixed);
				outFile.setf(std::ios::showpoint);
				outFile.open(savePath);
				for (int i = 0; i < config.change_count; i++)
				{
					results.emplace_back(
						pool.enqueue([&graph_list, &graph_time, &hop_info, &hop_info_2021, i, &outFile, &outMutex, &_random_v, &_random_time, &_random_hop, &rucTimeCostAll, &a2021TimeCostAll, &base1TimeCostAll, &base2TimeCostAll]
									 {
							experiment::ExecutionTimer timerQuery;
							timerQuery.startTask("query shorest path distance");
							int index_i = _random_v(boost_random_time_seed);
							int index_j = _random_v(boost_random_time_seed);
							int t_1 = _random_time(boost_random_time_seed);
							int t_2 = _random_time(boost_random_time_seed);
							int hop = _random_hop(boost_random_time_seed);
							if (t_1 > t_2)
							{
								std::swap(t_1, t_2);
							}
							timerQuery.startSubtask("====iteration " + std::to_string(i) + " query result info====");
							timerQuery.startSubtask("baseline 1: traverse each time graph");
							int resb1 = experiment::hop::dijkstra_iterator(graph_list, index_i, index_j, t_1, t_2, hop);
							auto base1Time = timerQuery.endSubtask();
							outMutex.lock();
							base1TimeCostAll+=base1Time;
							outMutex.unlock();
							timerQuery.startSubtask("baseline 2: traverse graph with time span");
							int resb2 = experiment::hop::search_shortest_path_in_period_time_naive(graph_time, index_i, index_j, hop, t_1, t_2);
							auto base2Time =timerQuery.endSubtask();
							outMutex.lock();
							base2TimeCostAll+=base2Time;
							outMutex.unlock();
							timerQuery.startSubtask("search result by ruc maintain algorithm");
							int ruc_res = hop_info.query(index_i, index_j, t_1, t_2, hop);
							auto time_ruc = timerQuery.endSubtask();
							outMutex.lock();
							rucTimeCostAll+=time_ruc;
							outMutex.unlock();
							timerQuery.startSubtask("search result by 2021 maintain algorithm");
							int res_2021 = hop_info_2021.query(index_i, index_j, t_1, t_2, hop);
							auto time_2021 = timerQuery.endSubtask();
							outMutex.lock();
							a2021TimeCostAll+=time_2021;
							outMutex.unlock();
							timerQuery.endSubtask();
							outMutex.lock();
							outFile << "from " << index_i << " to " << index_j << " between " << t_1 << " and " << t_2 << " by " << hop << std::endl;
							outFile << resb1 << ":" << resb2 << ":" << ruc_res << ":" << res_2021 << std::endl;
							timerQuery.writeStatsToFile(outFile);
							outMutex.unlock();
							return 1; }));
				}
				for (auto &&result : results)
					result.get(); // all threads finish here
				results.clear();
				outFile << "ruc1 cost " << rucTimeCostAll << std::endl;
				outFile << "2021 cost " << a2021TimeCostAll << std::endl;
				outFile << "baseline1 cost " << base1TimeCostAll << std::endl;
				outFile << "baseline2 cost " << base2TimeCostAll << std::endl;
				outFile.close();
			}
			else
			{
				double rucTimeCostAll = 0;
				double a2021TimeCostAll = 0;
				double base1TimeCostAll = 0;
				double base2TimeCostAll = 0;
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
					results.emplace_back(
						pool.enqueue([&graph_list, &graph_time, &hop_info, &hop_info_2021, i, &outFile, &outMutex, &_random_v, &_random_time,
									  &base1TimeCostAll, &base2TimeCostAll, &rucTimeCostAll, &a2021TimeCostAll]
									 {
							experiment::ExecutionTimer timerQuery;
							timerQuery.startTask("query shorest path distance");
							int index_i = _random_v(boost_random_time_seed);
							int index_j = _random_v(boost_random_time_seed);
							int t_1 = _random_time(boost_random_time_seed);
							int t_2 = _random_time(boost_random_time_seed);
							if (t_1 > t_2)
							{
								std::swap(t_1, t_2);
							}
							timerQuery.startSubtask("====iteration " + std::to_string(i) + " query result info====");
							timerQuery.startSubtask("baseline 1: traverse each time graph");
							int resb1 = experiment::nonhop::dijkstra_iterator(graph_list, index_i, index_j, t_1, t_2);
							auto base1Time = timerQuery.endSubtask();
							outMutex.lock();
							base1TimeCostAll+=base1Time;
							outMutex.unlock();
							timerQuery.startSubtask("baseline 2: traverse graph with time span");
							int resb2 = experiment::nonhop::search_shortest_path_in_period_time_naive(graph_time, index_i, index_j, t_1, t_2);
							auto base2Time =timerQuery.endSubtask();
							outMutex.lock();
							base2TimeCostAll+=base2Time;
							outMutex.unlock();
							timerQuery.startSubtask("search result by ruc maintain algorithm");
							int ruc_res = hop_info.query(index_i, index_j, t_1, t_2);
							auto time_ruc = timerQuery.endSubtask();
							outMutex.lock();
							rucTimeCostAll+=time_ruc;
							outMutex.unlock();
							timerQuery.startSubtask("search result by 2021 maintain algorithm");
							int res_2021 = hop_info_2021.query(index_i, index_j, t_1, t_2);
							auto time_2021 = timerQuery.endSubtask();
							outMutex.lock();
							a2021TimeCostAll+=time_2021;
							outMutex.unlock();
							timerQuery.endSubtask();
							outMutex.lock();
							outFile << "from " << index_i << " to " << index_j << " between " << t_1 << " and " << t_2 << std::endl;
							outFile << resb1 << ":" << resb2 << ":" << ruc_res << ":" << res_2021 << std::endl;
							timerQuery.writeStatsToFile(outFile);
							outMutex.unlock();
							return 1; }));
				}
				for (auto &&result : results)
					result.get(); // all threads finish here
				outFile << "ruc1 cost " << rucTimeCostAll << std::endl;
				outFile << "2021 cost " << a2021TimeCostAll << std::endl;
				outFile << "baseline1 cost " << base1TimeCostAll << std::endl;
				outFile << "baseline2 cost " << base2TimeCostAll << std::endl;
				outFile.close();
				results.clear();
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