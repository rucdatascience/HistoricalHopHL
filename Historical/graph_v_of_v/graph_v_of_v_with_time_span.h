#pragma once

#include <vector>
#include <iostream>
#include <vector>
#include "CPU/tool_functions/sorted_vector_binary_operations.h"
#include "CPU/graph_v_of_v/graph_v_of_v.h"
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
struct compare_pair
{
	bool operator()(const pair<int, weight_type> &lhs, const pair<int, weight_type> &rhs) const
	{
		return lhs.second > rhs.second;
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
	inline weight_type search_shortest_path_in_period_time_naive(int, int, int, int);

	inline void print();
	inline vector<graph_v_of_v<weight_type>> graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(int change_num, int maintain_percent);

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
weight_type graph_v_of_v_with_time_span<weight_type>::search_shortest_path_in_period_time_naive(int u, int v, int startTime, int endTime)
{
	double res = __DBL_MAX__;
	int N = this->v_num;
	std::vector<double> dist(N);
	std::vector<bool> visited(N);
	boost::heap::fibonacci_heap<pair<int, double>, boost::heap::compare<compare_pair<weight_type>>> queue;
	for (int queryTime = startTime; queryTime <= endTime; queryTime++)
	{
		dist.assign(N, __DBL_MAX__);
		visited.assign(N, false);
		queue.clear();
		dist[u] = 0;
		queue.push({u, 0});
		while (queue.size() > 0)
		{
			int vertexBase = queue.top().first;
			queue.pop();
			if (vertexBase == v)
			{
				res = min(res, dist[vertexBase]);
			}
			if (visited[vertexBase])
				continue;
			visited[vertexBase] = true;

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
							queue.push({next, dist[next]});
						}
						break;
					}
				}
			}
		}
		// cout << "naive" << res << endl;
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
		std::cout << "Vertex " << i << " Adj List: ";
		for (const auto &edges : ADJs[i])
		{
			int v_id = edges.first;
			for (const auto &info : edges.second)
			{
				std::cout << "<" << v_id << "," << info.weight << "," << info.startTimeLabel << "," << info.endTimeLabel << "> ";
			}
		}
		std::cout << std::endl;
	}
	std::cout << "graph_v_of_v_with_time_span_print END" << std::endl;
}

template <typename weight_type>
vector<graph_v_of_v<weight_type>> graph_v_of_v_with_time_span<weight_type>::graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(int change_num, int maintain_percent)
{
	if (maintain_percent > 10 || maintain_percent < 0)
	{
		cout << "maintain_percent should be between 0 and 10: " << maintain_percent << endl;
		exit(1);
	}
	if (change_num < 0)
	{
		cout << "the time_max should be greater than or equal to 0" << endl;
	}
	this->time_max = change_num;
	boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};
	graph_v_of_v<weight_type> instance_graph;
	instance_graph = graph_v_of_v_generate_random_graph<weight_type>(this->v_num, this->e_num, this->weight_dis.min(), this->weight_dis.max(), 1, boost_random_time_seed);
	add_graph_time(instance_graph, 0);
	vector<graph_v_of_v<weight_type>> res;
	res.push_back(instance_graph);
	uniform_int_distribution<> dis(1, 10);
	int index = 1;
	int N = instance_graph.ADJs.size();
	while (index <= change_num)
	{
		for (int i = 0; i < N; i++)
		{
			auto &vertices = instance_graph.ADJs[i];
			for (auto &edges : vertices)
			{
				if (edges.first > i && dis(boost_random_time_seed) > maintain_percent)
				{
					edges.second = weight_dis(boost_random_time_seed);
					int to_position = sorted_vector_binary_operations_search_position(instance_graph.ADJs[edges.first], i);
					instance_graph.ADJs[edges.first][to_position].second = edges.second;
					add_edge(i, edges.first, edges.second, index);
					add_edge(edges.first, i, edges.second, index);
				}
			}
		}
		res.push_back(instance_graph);
		++index;
	}
	return res;
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
weight_type dijkstra(graph_v_of_v<weight_type> graph, int u, int v)
{
	std::vector<double> dist(graph.size(), __DBL_MAX__);
	std::vector<bool> visited(graph.size(), false);
	boost::heap::fibonacci_heap<pair<int, weight_type>, boost::heap::compare<compare_pair<weight_type>>> queue;
	dist[u] = 0;
	queue.push({u, 0});
	while (queue.size() > 0)
	{
		int vertexBase = queue.top().first;
		queue.pop();
		if (vertexBase == v)
		{
			return dist[vertexBase];
		}
		if (visited[vertexBase])
			continue;
		visited[vertexBase] = true;

		for (const auto &edge : graph[vertexBase])
		{
			int next = edge.first;
			int weight = edge.second;

			if (dist[vertexBase] + weight < dist[next])
			{
				dist[next] = dist[vertexBase] + weight;
				queue.push({next, dist[next]});
			}
		}
	}
	return -1;
}

template <typename weight_type>
int dijkstra_iterator(vector<graph_v_of_v<weight_type>> list, int u, int v)
{
	int res = INT_MAX;
	for (graph_v_of_v<int> graph : list)
	{
		res = min(res, dijkstra(graph, u, v));
		// cout << "dijkstra" << res << endl;
	}
	return res;
}