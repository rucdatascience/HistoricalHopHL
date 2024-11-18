#pragma once

#include <iostream>
#include <chrono>
#include "CPU/tool_functions/sorted_vector_binary_operations.h"
#include "CPU/graph_v_of_v/graph_v_of_v.h"
#include "CPU/graph_v_of_v/graph_v_of_v_update_vertexIDs_by_degrees_large_to_small.h"
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightDecreaseMaintenance_improv_batch.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightDecrease2021_batch.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightIncrease2021_batch.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightIncreaseMaintenance_improv_batch.h"

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
		if (get<1>(lhs) == get<1>(rhs))
		{
			return get<2>(lhs) > get<2>(rhs);
		}
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
	virtual weight_type search_shortest_path_in_period_time_naive(int u, int v, int k, int startTime, int endTime) = 0;

	virtual vector<graph_v_of_v<weight_type>> graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(int change_num, int decreate_time, int increase_time, float change_ratio);

	virtual vector<graph_v_of_v<weight_type>> txt_read(std::string save_name);

	inline void print();

	inline void txt_save(std::string save_name);

protected:
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

	inline void process(graph_v_of_v<weight_type> &instance_graph, vector<pair<int, int>> &path, vector<int> &weight, int t);
};

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
					vector<EdgeInfo<weight_type>> list = ADJs[i][j].second;
					for (int k = 0; k < list.size() && list[k].startTimeLabel <= index; k++)
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
inline void graph_v_of_v_with_time_span<weight_type>::process(graph_v_of_v<weight_type> &instance_graph, vector<pair<int, int>> &path, vector<int> &weight, int t)
{
	for (int index = 0; index < path.size(); index++)
	{
		add_edge(path[index].first, path[index].second, weight[index], t);
		instance_graph.add_edge(path[index].first, path[index].second, weight[index]);
	}
}

template <typename weight_type>
weight_type dijkstra(graph_v_of_v<weight_type> &graph, int u, int v, int k)
{
	std::vector<weight_type> dist(graph.size(), std::numeric_limits<weight_type>::max());
	std::vector<int> hop_list(graph.size(), std::numeric_limits<int>::max());
	boost::heap::fibonacci_heap<std::tuple<int, weight_type, int>, boost::heap::compare<compare_tuple<weight_type>>> queue;

	dist[u] = 0;
	hop_list[u] = 0;
	queue.push({u, 0, 0});
	int res = __INT_MAX__;
	while (!queue.empty())
	{
		auto top = queue.top();
		int vertexBase = std::get<0>(top);
		weight_type currentDist = std::get<1>(top);
		int hop = std::get<2>(top);
		queue.pop();

		if (vertexBase == v)
		{
			res = min(res, currentDist);
		}
		if (hop >= k)
		{
			continue;
		}

		for (const auto &edge : graph[vertexBase])
		{
			int next = edge.first;
			weight_type weight = edge.second;

			weight_type newDist = currentDist + weight;

			if (newDist < dist[next] || (hop + 1 < hop_list[next]))
			{
				dist[next] = newDist;
				hop_list[next] = hop + 1;
				queue.push({next, newDist, hop + 1});
			}
		}
	}

	return res;
}

template <typename weight_type>
int dijkstra_iterator(vector<graph_v_of_v<weight_type>> list, int u, int v, int k)
{
	int res = INT_MAX;
	auto start_time = std::chrono::high_resolution_clock::now();
	for (graph_v_of_v<int> graph : list)
	{
		res = min(res, dijkstra(graph, u, v, k));
		// cout << "dijkstra" << res << endl;
	}
	auto endTime = std::chrono::high_resolution_clock::now();
	double runtime_n_iterate_dijkstra = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - start_time).count() / 1e9;
	std::cout << "dijkstra query time" << runtime_n_iterate_dijkstra << endl;
	return res;
};
