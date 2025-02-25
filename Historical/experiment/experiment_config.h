#pragma once
#include <filesystem>
#include <vector>
#include <queue>
#include "argparse/argparse.hpp"
#include "Historical/graph_with_time_span/graph_with_time_span.h"
#include "Historical/graph_with_time_span/graph.h"
#include "Historical/utils/BinaryPersistence.h"
#include "Historical/graph_with_time_span/two_hop_label.h"
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
namespace experiment {
	enum Mode { GENERATE_LABEL, MAINTAIN_LABEL, QUERY_RESULT };

	struct ExperimentConfig {
		enum Mode mode;
		int threads = 0;
		std::filesystem::path data_source;
		std::filesystem::path save_path;
		int hop_limit = 0;

		// maintain-label param
		int iterations = 0;
		int change_count = 0;
		int max_value = 0;
		int min_value = 0;
	};

	ExperimentConfig parse_arguments(int argc, char* argv[]) {
		argparse::ArgumentParser program("experiment");

		// generate-label
		argparse::ArgumentParser generate_label("generate-label");
		generate_label.add_argument("-t", "--threads").required().scan<'i', int>();
		generate_label.add_argument("-f", "--data_source").required();
		generate_label.add_argument("-p", "--save_path").required();
		generate_label.add_argument("-k", "--hop_limit").required().default_value(0).scan<'i', int>();

		// maintain-label
		argparse::ArgumentParser maintain_label("maintain-label");
		maintain_label.add_argument("-t", "--threads").required().scan<'i', int>();
		maintain_label.add_argument("-f", "--data_source").required();
		maintain_label.add_argument("-p", "--save_path").required();
		maintain_label.add_argument("-k", "--hop_limit").required().scan<'i', int>();
		maintain_label.add_argument("-m", "--iterations").required().scan<'i', int>();
		maintain_label.add_argument("-c", "--change_count").required().scan<'i', int>();
		maintain_label.add_argument("-max", "--max_value").required().scan<'i', int>();
		maintain_label.add_argument("-min", "--min_value").required().scan<'i', int>();

		//query-result
		argparse::ArgumentParser query_label("query-result");
		query_label.add_argument("-f", "--data_source").required();
		query_label.add_argument("-c", "--search_count").required().scan<'i', int>();
		query_label.add_argument("-k", "--hop_limit").required().scan<'i', int>();

		program.add_subparser(generate_label);
		program.add_subparser(maintain_label);
		program.add_subparser(query_label);

		try {
			program.parse_args(argc, argv);
		}
		catch (const std::runtime_error& err) {
			std::cerr << "Error: " << err.what() << "\n";
			std::cerr << program;
			exit(EXIT_FAILURE);
		}

		ExperimentConfig config;
		if (program.is_subcommand_used("generate-label")) {
			config.mode = GENERATE_LABEL;
			config.threads = generate_label.get<int>("-t");
			config.data_source = generate_label.get<std::string>("-f");
			config.save_path = generate_label.get<std::string>("-p");
			config.hop_limit = generate_label.get<int>("-k");
			if (config.hop_limit < 0) {
				throw std::invalid_argument("Error: hop_constrained (-k) must be >= 0.");
			}
		}
		else if (program.is_subcommand_used("maintain-label")) {
			config.mode = MAINTAIN_LABEL;
			config.threads = maintain_label.get<int>("-t");
			config.data_source = maintain_label.get<std::string>("-f");
			config.save_path = maintain_label.get<std::string>("-p");
			config.hop_limit = maintain_label.get<int>("-k");
			config.iterations = maintain_label.get<int>("-m");
			config.change_count = maintain_label.get<int>("-c");
			config.max_value = maintain_label.get<int>("-max");
			config.min_value = maintain_label.get<int>("-min");

			if (config.max_value <= 0 || config.min_value <= 0) {
				throw std::invalid_argument("Error: max_value and min_value must be greater than 0.");
			}
		}
		else if (program.is_subcommand_used("query-result")) {
			config.mode = QUERY_RESULT;
			config.data_source = query_label.get<std::string>("-f");
			config.change_count = query_label.get<int>("-c");
			config.hop_limit = query_label.get<int>("-k");
		}
		else {
			std::cerr << "Error: Unknown subcommand.\n";
			std::cerr << program;
			exit(EXIT_FAILURE);
		}

		return config;
	}

	std::vector<std::string> parse_string(std::string parse_target, std::string delimiter)
	{
		std::vector<std::string> Parsed_content;
		size_t pos = 0;
		std::string token;
		while ((pos = parse_target.find(delimiter)) != std::string::npos) {
			// find(const string& str, size_t pos = 0) function returns the position of the first occurrence of str in the string, or npos if the string is not found.
			token = parse_target.substr(0, pos);
			// The substr(size_t pos = 0, size_t n = npos) function returns a substring of the object, starting at position pos and of length npos
			Parsed_content.push_back(token); // store the subtr to the list
			parse_target.erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
		}
		Parsed_content.push_back(parse_target); // store the subtr to the list

		return Parsed_content;
	}

	struct change_edge_info
	{
		int v1;
		int v2;
		int weight;
		int time;
	};
	template <typename weight_type>
	class iteration_info {
	private:
		const int _v_num;
		const int _iteration;
		const int _change_num;
		const int _upper;
		const int _lower;
		boost::random::uniform_int_distribution<> _random_v;
		boost::random::uniform_int_distribution<> _random_weight;
		graph <weight_type> instance_graph;
	public:
		// 保存每一个time slot的变化队列
		std::vector<std::queue<change_edge_info>> q_list;
		iteration_info(int v_num, int iteration, int change_num, int upper, int lower, graph<weight_type> graph) :_v_num(v_num), _iteration(iteration), _change_num(change_num), _upper(upper), _lower(lower), instance_graph(graph) {
			this->_random_v = boost::random::uniform_int_distribution<>(0, this->_v_num);
			this->_random_weight = boost::random::uniform_int_distribution<>(this->_lower, this->_upper);
			q_list = std::vector<std::queue<change_edge_info>>(this->_iteration + 1, std::queue<change_edge_info>());
		}

		void build_random_change() {
			std::map<std::pair<int, int>, int> pair2dis;
			for (int i = 1; i <= this->_iteration; i++)
			{
				for (int j = 0; j < this->_change_num; j++)
				{
					int index_i = this->_random_v(boost_random_time_seed);
					boost::random::uniform_int_distribution<> dis_inner(0, instance_graph[index_i].size() - 1);
					int index_j = dis_inner(boost_random_time_seed);
					int i_j_weight = this->_random_weight(boost_random_time_seed);
					change_edge_info info = { index_i, index_j, i_j_weight, i };
					std::pair index = std::make_pair(index_j, index_j);
					q_list[i].push(info);
				}
			}
		}
	};
}