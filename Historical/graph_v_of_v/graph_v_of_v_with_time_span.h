#pragma once

#include <iostream>
#include <chrono>
#include "CPU/tool_functions/sorted_vector_binary_operations.h"
#include "CPU/graph_v_of_v/graph_v_of_v.h"
#include "CPU/graph_v_of_v/graph_v_of_v_update_vertexIDs_by_degrees_large_to_small.h"
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightDecreaseMaintenance_improv_batch.h"
#include "CPU/tool_functions/ThreadPool.h"
#include <boost/heap/fibonacci_heap.hpp>
template <typename weight_type> // weight_type may be int, long long int, float, double...
class EdgeInfo
{
public:
	int vertex;
	weight_type weight;
	int startTimeLabel;
	int endTimeLabel;
	EdgeInfo() : vertex(-1) {}
	EdgeInfo(int vertex, weight_type weight, int startTimeLabel)
		: vertex(vertex), weight(weight), startTimeLabel(startTimeLabel)
	{
		endTimeLabel = INT_MAX;
	}
};
template <typename weight_type>
struct compare_tuple
{
	bool operator()(const tuple<int, weight_type, int> &lhs, const tuple<int, weight_type, int> &rhs) const
	{
		return get<1>(lhs) > get<1>(rhs);
	}
};

/**
 * define a graph where each edge has an associated time span
 */
template <typename weight_type> // weight_type may be int, long long int, float, double...
class graph_v_of_v_with_time_span
{
public:
	/*
	this class only suits ideal vertex IDs: from 0 to V-1;

	this class is for undirected and edge-weighted-time-span-label graph
	*/
	// i - 0,j - 0->pair->first=1,second = vector:{<1,23,0,0> <1,30,1,1> <1,11,2,3> <1,8,4,6> <1,17,7,2147483647>}
	vector<vector<pair<int, vector<EdgeInfo<weight_type>>>>> ADJs;

	/*constructors*/
	graph_v_of_v_with_time_span() {}
	graph_v_of_v_with_time_span(int n, int e, weight_type weight_upper_limit, weight_type weight_lower_limit) : v_num(n), e_num(e), weight_dis(weight_lower_limit, weight_upper_limit), ADJs(n)
	{
		if (weight_lower_limit > weight_upper_limit || weight_lower_limit < 0 || weight_upper_limit < 0)
		{
			cout << "weight_lower_limit should be greater than 0 and weight_upper_limit shoule be greater 0 and the weight_lower_limit should be less than weight_upper_limit"
				 << "weight_upper_limit is :" << weight_upper_limit
				 << "weight_lower_limit is :" << weight_lower_limit << endl;
			exit(1);
		}
	}
	int size()
	{
		return ADJs.size();
	}
	std::vector<EdgeInfo<weight_type>> &operator[](int i)
	{
		return ADJs[i];
	}
	/*class member functions*/
	inline weight_type search_shortest_path_in_period_time_naive(int u, int v, int k, int startTime, int endTime);

	inline void print();
	inline vector<graph_v_of_v<weight_type>> graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(int change_num, int maintain_percent, hop_constrained_case_info &info);

	inline void txt_save(std::string save_name);
	inline vector<graph_v_of_v<weight_type>> txt_read(std::string save_name, hop_constrained_case_info &info);

private:
	/* the maximum of time*/
	int time_max;
	/* the number of vertices */
	int v_num;
	/* the number of edges */
	int e_num;
	/* weighted random number generator */
	uniform_int_distribution<weight_type> weight_dis;

	inline void add_edge(int, int, weight_type, int);

	inline void add_graph_time(graph_v_of_v<weight_type>, int);

	inline void clear();
};

/*class member functions*/
/**
 * param
 * 	@u is the source vertex,
 *  @v is the target vertex,
 * 	@startTime,
 * 	@endTime
 */
template <typename weight_type>
weight_type graph_v_of_v_with_time_span<weight_type>::search_shortest_path_in_period_time_naive(int u, int v, int k, int startTime, int endTime)
{
	double res = __DBL_MAX__;
	int N = this->v_num;
	std::vector<double> dist(N);
	std::vector<bool> visited(N);
	boost::heap::fibonacci_heap<tuple<int, double, int>, boost::heap::compare<compare_tuple<weight_type>>> queue;
	for (int queryTime = startTime; queryTime <= endTime; queryTime++)
	{
		dist.assign(N, __DBL_MAX__);
		visited.assign(N, false);
		queue.clear();
		dist[u] = 0;
		queue.push({u, 0, 0});
		while (queue.size() > 0)
		{
			int hop = get<2>(queue.top());
			int vertexBase = get<0>(queue.top());
			queue.pop();
			if (vertexBase == v)
			{
				res = min(res, dist[vertexBase]);
				break;
			}
			if (visited[vertexBase])
				continue;
			visited[vertexBase] = true;

			if (hop == k)
			{
				continue;
			}

			for (const auto &vertices : this->ADJs[vertexBase])
			{
				int next = vertices.first;
				for (const auto &edge_info_time_span : vertices.second)
				{
					if (edge_info_time_span.startTimeLabel <= queryTime && edge_info_time_span.endTimeLabel >= queryTime)
					{
						if (dist[vertexBase] + edge_info_time_span.weight < dist[next])
						{
							dist[next] = dist[vertexBase] + edge_info_time_span.weight;
							queue.push({next, dist[next], hop + 1});
						}
						break;
					}
				}
			}
		}
		cout << "naive" << res << endl;
	}
	return res == __DBL_MAX__ ? -1 : res;
}

template <typename weight_type>
void graph_v_of_v_with_time_span<weight_type>::print()
{
	std::cout << "graph_v_of_v_with_time_span_print:" << std::endl;
	int size = this->ADJs.size();
	for (int i = 0; i < size; i++)
	{
		std::cout << "Vertex " << i << " Adj List: " << endl;
		for (const auto &edges : ADJs[i])
		{
			int v_id = edges.first;
			cout << "\t";
			for (const auto &info : edges.second)
			{
				std::cout << "<" << v_id << "," << info.weight << "," << info.startTimeLabel << "," << info.endTimeLabel << "> ";
			}
			cout << endl;
		}
	}
	std::cout << "graph_v_of_v_with_time_span_print END" << std::endl;
}

template <typename weight_type>
vector<graph_v_of_v<weight_type>> graph_v_of_v_with_time_span<weight_type>::graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(int change_num, int maintain_percent, hop_constrained_case_info &case_info)
{
	if (maintain_percent > 10 || maintain_percent < 0)
	{
		cout << "maintain_percent should be between 0 and 10: " << maintain_percent << endl;
		exit(1);
	}
	if (change_num < 0)
	{
		cout << "the change_num should be greater than or equal to 0" << endl;
	}
	this->time_max = change_num;
	boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};
	graph_v_of_v<weight_type> instance_graph;
	instance_graph = graph_v_of_v_generate_random_graph<weight_type>(this->v_num, this->e_num, this->weight_dis.min(), this->weight_dis.max(), 1, boost_random_time_seed);
	vector<int> is_mock(instance_graph.size());
	for (int i = 0; i < instance_graph.size(); i++)
	{
		is_mock[i] = false;
	}
	instance_graph = graph_v_of_v_update_vertexIDs_by_degrees_large_to_small_mock(instance_graph, is_mock);

	// initialize the label
	hop_constrained_two_hop_labels_generation(instance_graph, case_info);
	ThreadPool pool_dynamic(case_info.thread_num);
	std::vector<std::future<int>> results_dynamic;

	add_graph_time(instance_graph, 0);
	vector<graph_v_of_v<weight_type>> res;
	res.push_back(instance_graph);
	uniform_int_distribution<> dis(1, 10);
	int index = 1;
	int N = instance_graph.ADJs.size();
	while (index <= change_num)
	{
		vector<pair<int, int>> path;
		vector<int> weight;
		for (int i = 0; i < N; i++)
		{
			auto &vertices = instance_graph.ADJs[i];
			for (auto &edges : vertices)
			{
				if (edges.first > i && dis(boost_random_time_seed) > maintain_percent)
				{
					auto next_value = weight_dis(boost_random_time_seed);
					// only edge weight decrease
					if (edges.second > next_value)
					{
						edges.second = next_value;
						int to_position = sorted_vector_binary_operations_search_position(instance_graph.ADJs[edges.first], i);
						instance_graph.ADJs[edges.first][to_position].second = edges.second;
						add_edge(i, edges.first, edges.second, index);
						add_edge(edges.first, i, edges.second, index);
						// maintain the label
						path.push_back({edges.first, i});
						weight.push_back(next_value);
						// 5 is batch_size,it should be definable
						if (path.size() >= 10)
						{
							HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, path, weight, pool_dynamic, results_dynamic, index);
							vector<pair<int, int>>().swap(path);
							vector<int>().swap(weight);
						}
					}
				}
			}
		}
		if (path.size() > 0)
		{
			HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, path, weight, pool_dynamic, results_dynamic, index);
			vector<pair<int, int>>().swap(path);
			vector<int>().swap(weight);
		}
		res.push_back(instance_graph);
		++index;
	}
	return res;
}

template <typename weight_type>
inline void graph_v_of_v_with_time_span<weight_type>::txt_save(std::string save_name)
{
	std::ofstream outputFile;
	outputFile.precision(10);
	outputFile.setf(std::ios::fixed);
	outputFile.setf(std::ios::showpoint);
	outputFile.open(save_name);

	outputFile << "|V|= " << this->v_num << std::endl;
	outputFile << "|E|= " << this->e_num << std::endl;
	outputFile << "|time|= " << this->time_max << std::endl;
	outputFile << std::endl;

	for (int index = 0; index <= this->time_max; index++)
	{
		outputFile << "time " << index << std::endl;
		for (int i = 0; i < this->v_num; i++)
		{
			for (int j = 0; j < this->ADJs[i].size(); j++)
			{
				if (i < this->ADJs[i][j].first)
				{
					for (int k = 0; this->ADJs[i][j].second[k].startTimeLabel <= index; k++)
					{
						if (this->ADJs[i][j].second[k].startTimeLabel == index)
						{
							outputFile << "Edge " << i << " " << ADJs[i][j].first << " " << this->ADJs[i][j].second[k].weight << "\n";
							break;
						}
					}
				}
			}
		}
		outputFile << std::endl;
	}
	outputFile << "EOF" << std::endl;
	outputFile.close();
}

template <typename weight_type>
inline vector<graph_v_of_v<weight_type>> graph_v_of_v_with_time_span<weight_type>::txt_read(std::string save_name, hop_constrained_case_info &case_info)
{
	this->clear();
	std::string line_content;
	int current_time = -1;
	vector<graph_v_of_v<weight_type>> res;
	graph_v_of_v<weight_type> instance_graph;
	ThreadPool pool_dynamic(case_info.thread_num);
	std::vector<std::future<int>> results_dynamic;
	vector<pair<int, int>> path;
	vector<int> weight;

	std::ifstream myfile(save_name); // open the file
	if (myfile.is_open())			 // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			std::vector<std::string> Parsed_content = parse_string(line_content, " ");

			if (!Parsed_content[0].compare("|V|=")) // when it's equal, compare returns 0
			{
				this->v_num = std::stoi(Parsed_content[1]);
				instance_graph.ADJs.resize(this->v_num);
				ADJs.resize(std::stoi(Parsed_content[1]));
				initialize_global_values_dynamic_hop_constrained(this->v_num, case_info.thread_num, case_info.upper_k);
			}
			else if (!Parsed_content[0].compare("|E|="))
			{
				this->e_num = std::stoi(Parsed_content[1]);
			}
			else if (!Parsed_content[0].compare("|time|="))
			{
				this->time_max = std::stoi(Parsed_content[1]);
			}
			else if (!Parsed_content[0].compare("time"))
			{
				current_time = std::stoi(Parsed_content[1]);
			}
			else if (!Parsed_content[0].compare("Edge"))
			{
				int v1 = std::stoi(Parsed_content[1]);
				int v2 = std::stoi(Parsed_content[2]);
				weight_type ec = std::stod(Parsed_content[3]);
				if (current_time == 0)
				{
					// initiate the graph
					instance_graph.add_edge(v1, v2, ec);
				}
				else
				{
					instance_graph.add_edge(v1, v2, ec);
					add_edge(v1, v2, ec, current_time);
					add_edge(v2, v1, ec, current_time);
					// TODO current modification operations are only decrease
					// maintain the label
					path.push_back({v1, v2});
					weight.push_back(ec);

					// 5 is batch_size,it should be definable
					if (path.size() >= case_info.thread_num * 3)
					{
						HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, path, weight, pool_dynamic, results_dynamic, current_time);
						vector<pair<int, int>>().swap(path);
						vector<int>().swap(weight);
					}
				}
			}
			else if (Parsed_content.size() == 1 && Parsed_content[0] == "")
			{
				if (path.size() > 0)
				{
					HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, path, weight, pool_dynamic, results_dynamic, current_time);
					vector<pair<int, int>>().swap(path);
					vector<int>().swap(weight);
				}
				if (current_time >= 0)
				{
					if (current_time == 0)
					{
						hop_constrained_two_hop_labels_generation(instance_graph, case_info);
						add_graph_time(instance_graph, 0);
					}
					res.push_back(instance_graph);
				}
			}
		}
		myfile.close(); // close the file
		return res;
	}
	else
	{
		std::cout << "Unable to open file " << save_name << std::endl
				  << "Please check the file location or file name." << std::endl; // throw an error message
		getchar();																  // keep the console window
		exit(1);																  // end the program
	}
}

template <typename weight_type>
void graph_v_of_v_with_time_span<weight_type>::add_edge(int e1, int e2, weight_type ec, int time)
{
	/* initialize a graph with a time span */
	if (time == 0)
	{
		this->ADJs[e1].push_back({e2, {EdgeInfo(e2, ec, time)}});
		this->ADJs[e2].push_back({e1, {EdgeInfo(e1, ec, time)}});
	}
	else
	{
		int index_e1 = sorted_vector_binary_operations_search_position<vector<EdgeInfo<weight_type>>>(this->ADJs[e1], e2);
		int index_e2 = sorted_vector_binary_operations_search_position<vector<EdgeInfo<weight_type>>>(this->ADJs[e2], e1);
		if (this->ADJs[e1][index_e1].second.back().weight != ec)
		{
			this->ADJs[e1][index_e1].second.back().endTimeLabel = time - 1;
			this->ADJs[e2][index_e2].second.back().endTimeLabel = time - 1;
			this->ADJs[e1][index_e1].second.push_back(EdgeInfo(e2, ec, time));
			this->ADJs[e2][index_e2].second.push_back(EdgeInfo(e1, ec, time));
		}
	}
}

template <typename weight_type>
void graph_v_of_v_with_time_span<weight_type>::add_graph_time(graph_v_of_v<weight_type> graph, int time)
{
	int N = graph.size();
	for (int i = 0; i < N; i++)
	{
		std::vector<std::pair<int, weight_type>> list = graph.ADJs[i];
		for (const auto &edges : list)
		{
			if (edges.first < i)
			{
				continue;
			}
			this->add_edge(i, edges.first, edges.second, time);
		}
	}
}

template <typename weight_type>
inline void graph_v_of_v_with_time_span<weight_type>::clear()
{
	vector<vector<pair<int, vector<EdgeInfo<weight_type>>>>>().swap(this->ADJs);
	this->e_num = 0;
	this->v_num = 0;
	this->time_max = 0;
}

template <typename weight_type>
weight_type dijkstra(graph_v_of_v<weight_type> graph, int u, int v, int k)
{
	std::vector<double> dist(graph.size(), __DBL_MAX__);
	std::vector<bool> visited(graph.size(), false);
	boost::heap::fibonacci_heap<tuple<int, weight_type, int>, boost::heap::compare<compare_tuple<weight_type>>> queue;
	dist[u] = 0;
	queue.push({u, 0, 0});
	while (queue.size() > 0)
	{
		tuple<int, weight_type, int> top = queue.top();
		int hop = get<2>(top);
		int vertexBase = get<0>(top);
		queue.pop();
		if (vertexBase == v)
		{
			return dist[vertexBase];
		}
		if (visited[vertexBase])
			continue;
		visited[vertexBase] = true;
		if (hop == k)
		{
			continue;
		}
		for (const auto &edge : graph[vertexBase])
		{
			int next = edge.first;
			int weight = edge.second;

			if (dist[vertexBase] + weight < dist[next])
			{
				dist[next] = dist[vertexBase] + weight;
				queue.push({next, dist[next], hop + 1});
			}
		}
	}
	return -1;
}

template <typename weight_type>
int dijkstra_iterator(vector<graph_v_of_v<weight_type>> list, int u, int v, int k)
{
	int res = INT_MAX;
	auto start_time = std::chrono::high_resolution_clock::now();
	for (graph_v_of_v<int> graph : list)
	{
		res = min(res, dijkstra(graph, u, v, k));
		cout << "dijkstra" << res << endl;
	}
	auto endTime = std::chrono::high_resolution_clock::now();
	double runtime_n_iterate_dijkstra = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - start_time).count() / 1e9;
	std::cout << runtime_n_iterate_dijkstra << endl;
	return res;
}
