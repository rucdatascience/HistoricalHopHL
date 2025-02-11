#include <cstdlib>
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span_experiment_reader.h"

void print_usage() {
	std::cout << "Usage: experiment -m <iterations> -c <change_count> -p <save_path> -t <threads> -f <data_source> -max <max_value> -min <min_value> [-d] [-L]" << std::endl;
	std::exit(EXIT_FAILURE);
}

ExperimentConfig parse_arguments(int argc, char* argv[]) {
	if (argc < 15) {
		print_usage();
	}

	ExperimentConfig config;
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "-m" && i + 1 < argc) {
			config.iterations = std::stoi(argv[++i]);
			if (config.iterations <= 0) {
				std::cerr << "Error: Iterations (-m) must be greater than 0." << std::endl;
				print_usage();
			}
		}
		else if (arg == "-c" && i + 1 < argc) {
			config.change_count = std::stoi(argv[++i]);
			if (config.change_count <= 0) {
				std::cerr << "Error: Change count (-c) must be greater than 0." << std::endl;
				print_usage();
			}
		}
		else if (arg == "-p" && i + 1 < argc) {
			config.save_path = std::filesystem::path(argv[++i]);
		}
		else if (arg == "-t" && i + 1 < argc) {
			config.thread_count = std::stoi(argv[++i]);
		}
		else if (arg == "-f" && i + 1 < argc) {
			config.data_source = argv[++i];
		}
		else if (arg == "-max" && i + 1 < argc) {
			config.max_value = std::stoi(argv[++i]);
			if (config.max_value <= 0) {
				std::cerr << "Error: Max value (-max) must be greater than 0." << std::endl;
				print_usage();
			}
		}
		else if (arg == "-min" && i + 1 < argc) {
			config.min_value = std::stoi(argv[++i]);
			if (config.min_value <= 0) {
				std::cerr << "Error: Min value (-min) must be greater than 0." << std::endl;
				print_usage();
			}
		}
		else if (arg == "-d") {
			config.debug = true;
		}
		else if (arg == "-l") {
			config.markL = true;
		}
		else {
			std::cerr << "Unknown or malformed argument: " << arg << std::endl;
			print_usage();
		}
	}
	return config;
}

int main(int argc, char* argv[]) {
	ExperimentConfig config = parse_arguments(argc, argv);

	std::cout << "Iterations: " << config.iterations << "\n"
		<< "Change Count: " << config.change_count << "\n"
		<< "Save Path: " << config.save_path << "\n"
		<< "Threads: " << config.thread_count << "\n"
		<< "Data Source: " << config.data_source << "\n"
		<< "Debug Mode: " << (config.debug ? "Enabled" : "Disabled") << "\n"
		<< "Max edge weight:" << config.max_value << "\n"
		<< "Min edge weight:" << config.min_value << "\n";

	experiment_config experiment(config);
	experiment.init();
	experiment.process_2021();
	experiment.process();
	experiment.print_experiment_result();
	experiment.close();

	return 0;
}
