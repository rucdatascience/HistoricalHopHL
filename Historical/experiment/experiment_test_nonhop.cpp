#include <iostream>
#include <string>
#include <filesystem>
#include "Historical/utils/BinaryPersistence.h"
#include "Historical/experiment/experiment_operation.h"
#include "Historical/experiment/experiment_config.h"
#include "Historical/graph_with_time_span/two_hop_label.h"
#include <Historical/utils/ExecutionTimer .h>
experiment::ExecutionTimer timer;
int main(int argc, char* argv[]) {
	try {
		experiment::ExperimentConfig config = experiment::parse_arguments(argc, argv);

		std::cout << "Mode: " << (config.mode == experiment::ExperimentConfig::GENERATE_LABEL ? "Generate Label" : "Maintain Label") << "\n"
			<< "Threads: " << config.threads << "\n"
			<< "Data Source: " << config.data_source << "\n"
			<< "Save Path: " << config.save_path << "\n"
			<< "Hop Limit (k): " << config.hop_limit << "\n";

		if (config.mode == experiment::ExperimentConfig::MAINTAIN_LABEL) {
			std::cout << "Iterations: " << config.iterations << "\n"
				<< "Change Count: " << config.change_count << "\n"
				<< "Max Value: " << config.max_value << "\n"
				<< "Min Value: " << config.min_value << "\n";
		}
		if (config.mode == experiment::ExperimentConfig::GENERATE_LABEL) {
			experiment::graph<int> instance_graph;
			experiment::read_graph(instance_graph, config);
			experiment::graph_with_time_span<int> graph_time;
			graph_time.add_graph_time(instance_graph, 0);
			std::vector<experiment::graph<int>> graph_list(1);
			graph_list.push_back(instance_graph);
			std::filesystem::path saveDir = std::filesystem::path(config.save_path);
			std::filesystem::create_directory(saveDir);
			if (config.hop_limit != 0) {
				//k-constrained pll
				std::string graph_res_filename = "binary_hop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_graph";
				std::string experiment_res_filename = "hop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::filesystem::path graphPath = saveDir.string() + "//" + graph_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_res_filename;
				experiment::hop::two_hop_case_info hop_info;
				hop_info.thread_num = config.threads;
				hop_info.upper_k = config.hop_limit;
				timer.startTask("generate graph and 2hop label " + std::to_string(config.hop_limit) + " hop constrained");
				experiment::hop::pll(instance_graph, hop_info);
				//hop_info.print_L();
				std::ofstream FILE(resultPath.string(), std::ios::out | std::ifstream::binary);
				experiment::saveBinary(FILE, graph_list);
				experiment::saveBinary(FILE, graph_time);
				experiment::saveBinary(FILE, hop_info);
				FILE.close();
				timer.writeStatsToFile(resultPath.string());
			}
			else {
				std::string graph_res_filename = "binary_nonhop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_graph";
				std::string experiment_res_filename = "nonhop_constrained_" + std::to_string(config.hop_limit) + "_" + std::to_string(config.threads) + "_threads_result.txt";
				std::filesystem::path graphPath = saveDir.string() + "//" + graph_res_filename;
				std::filesystem::path resultPath = saveDir.string() + "//" + experiment_res_filename;
				experiment::nonhop::two_hop_case_info hop_info;
				hop_info.thread_num = config.threads;
				timer.startTask("generate graph and 2hop label " + std::to_string(config.hop_limit) + " hop constrained");
				experiment::nonhop::pll(instance_graph, hop_info);
				//hop_info.print_L();
				//experiment::experiment_1_generate_pll_nonhop_result<int> res(graph_list, graph_time, hop_info);
				std::ofstream FILE_GRAPH(graphPath.string(), std::ios::out | std::ifstream::binary);
				experiment::saveBinary(FILE_GRAPH, graph_list);
				experiment::saveBinary(FILE_GRAPH, graph_time);
				experiment::saveBinary(FILE_GRAPH, hop_info);
				FILE_GRAPH.close();
				timer.writeStatsToFile(resultPath.string());
			}
		}
	}
	catch (const std::exception& ex) {
		std::cerr << "Error: " << ex.what() << "\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}