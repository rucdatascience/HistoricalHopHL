#pragma once

#include <vector>
#include <iostream>
#include <vector>
#include "Historical/tool_functions/sorted_vector_binary_operations.h"
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
	EdgeInfo(int vertex, weight_type weight, int startTimeLabel, int endTimeLabel)
		: vertex(vertex), weight(weight), startTimeLabel(startTimeLabel), endTimeLabel(endTimeLabel) {}
	bool operator>(int v_id) const
	{
		return vertex > v_id;
	}
	bool operator==(int v_id) const
	{
		return vertex == v_id;
	}
};

template <typename weight_type> // weight_type may be int, long long int, float, double...
class NodeWithTimeSpan
{
public:
	int vertex;
	weight_type weight;
	int startTimeLabel;
	int endTimeLabel;
	NodeWithTimeSpan(int vertex, weight_type weight, int startTimeLabel, int endTimeLabel)
		: vertex(vertex), weight(weight), startTimeLabel(startTimeLabel), endTimeLabel(endTimeLabel) {}
	NodeWithTimeSpan(int vertex, weight_type weight) : vertex(vertex), weight(weight) {}
	bool operator>(weight_type w) const
	{
		return weight > w;
	}
	bool operator==(weight_type w) const
	{
		return weight == w;
	}
	bool operator<(weight_type w) const
	{
		return weight < w;
	}
};

struct compare_node
{
	template <typename weight_type>
	bool operator()(const NodeWithTimeSpan<weight_type> &lhs, const NodeWithTimeSpan<weight_type> &rhs) const
	{
		return lhs.weight > rhs.weight;
	}
};

struct compare_pair
{
	bool operator()(const pair<int, int> &lhs, const pair<int, int> &rhs) const
	{
		return lhs.second > rhs.second;
	}
};

template <typename weight_type> // weight_type may be int, long long int, float, double...
class graph_v_of_v_with_time_span
{
public:
	/*
	this class only suits ideal vertex IDs: from 0 to V-1;

	this class is for undirected and edge-weighted-time-span-label graph
	*/
	std::vector<std::vector<EdgeInfo<weight_type>>> ADJs;

	/*constructors*/
	graph_v_of_v_with_time_span() {}
	graph_v_of_v_with_time_span(int n)
	{
		ADJs.resize(n); // initialize n vertices
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
	inline void add_edge(int, int, weight_type, int);
	inline void add_graph_time(graph_v_of_v<weight_type>, int);
	inline int search_shortest_path_in_period_time(int, int, int, int);
	inline void print();
	// inline int dijkstra_iterator(vector<graph_v_of_v<weight_type>>, int, int);
};

/*class member functions*/

template <typename weight_type>
void graph_v_of_v_with_time_span<weight_type>::add_edge(int e1, int e2, weight_type ec, int time)
{

	/*we assume that the size of g is larger than e1 or e2;
	 this function can update edge weight; there will be no redundent edge*/

	/*
	Add the edges (e1,e2) and (e2,e1) with the weight ec
	When the edge exists, it will update its weight.
	Time complexity:
		O(log n) When edge already exists in graph
		O(n) When edge doesn't exist in graph
	*/
	int index = sorted_vector_binary_operations(this->ADJs[e1], e2);
	/* e2 is the target vertex, and its weight has not changed,and its endTime being the previous time , which means we only need to extend the range of the timespan */
	if (this->ADJs[e1].size() != index && this->ADJs[e1][index].vertex == e2 && this->ADJs[e1][index].weight == ec && this->ADJs[e1][index].endTimeLabel == time - 1)
	{
		this->ADJs[e1][index].endTimeLabel = time;
	}
	else if (this->ADJs[e1].size() == index || this->ADJs[e1][index].vertex != e2)
	{
		EdgeInfo<weight_type> edgeInfo(e2, ec, time, time);
		this->ADJs[e1].insert(ADJs[e1].begin() + index, edgeInfo);
	}
	else
	{
		EdgeInfo<weight_type> edgeInfo(e2, ec, time, time);
		this->ADJs[e1].insert(this->ADJs[e1].begin() + index + 1, edgeInfo);
	}
	index = sorted_vector_binary_operations(this->ADJs[e2], e1);
	if (this->ADJs[e2].size() != index && this->ADJs[e2][index].vertex == e1 && this->ADJs[e2][index].weight == ec && this->ADJs[e2][index].endTimeLabel == time - 1)
	{
		this->ADJs[e2][index].endTimeLabel = time;
	}
	else if (this->ADJs[e2].size() == index || this->ADJs[e2][index].vertex != e1)
	{
		EdgeInfo<weight_type> edgeInfo(e1, ec, time, time);
		this->ADJs[e2].insert(ADJs[e2].begin() + index, edgeInfo);
	}
	else
	{
		EdgeInfo<weight_type> edgeInfo(e1, ec, time, time);
		this->ADJs[e2].insert(ADJs[e2].begin() + index + 1, edgeInfo);
	}
}

template <typename weight_type>
void graph_v_of_v_with_time_span<weight_type>::add_graph_time(graph_v_of_v<weight_type> graph, int time)
{
	int N = graph.size();
	if (N > this->ADJs.size())
	{
		this->ADJs.resize(N);
	}
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

/**
 *
 *
 */
/**
 * param
 * 	@u is the source vertex,
 *  @v is the target vertex,
 * 	@startTime,
 * 	@endTime
 */
template <typename weight_type>
int graph_v_of_v_with_time_span<weight_type>::search_shortest_path_in_period_time(int u, int v, int startTime, int endTime)
{
	if (v > ADJs.size() || u > ADJs.size() || u < 0 || v < 0)
	{
		std::cout << "u and v shoule be between 0 and N" << endl;
		return -1;
	}
	/* a class contains information about destination vertex, hop count, and cost in the priority queue */
	boost::heap::fibonacci_heap<NodeWithTimeSpan<weight_type>, boost::heap::compare<compare_node>> queue;
	NodeWithTimeSpan node(u, 0, startTime, endTime);
	queue.push(node);
	int res = -1;

	vector<vector<int>> visited;
	int N = this->ADJs.size();
	visited.resize(N);
	for (int i = 0; i < N; i++)
	{
		visited[i].resize(endTime - startTime + 1);
	}
	while (queue.size() > 0)
	{
		node = queue.top();
		queue.pop();
		if (node.vertex == v)
		{
			return node.weight;
		}
		int vertexBase = node.vertex;
		int curWeight = node.weight;
		int startTimeLabel = node.startTimeLabel;
		int endTimeLabel = node.endTimeLabel;
		if (visited[vertexBase][startTimeLabel] != 0 && visited[vertexBase][endTimeLabel] != 0)
		{
			continue;
		}
		for (int index = startTimeLabel - startTime; index <= endTimeLabel - startTime; index++)
		{
			visited[vertexBase][index] = 1;
		}
		for (const auto &edges : this->ADJs[vertexBase])
		{
			if (!(edges.startTimeLabel > endTimeLabel || edges.endTimeLabel < startTimeLabel))
			{
				node.vertex = edges.vertex;
				node.weight = curWeight + edges.weight;
				node.startTimeLabel = max(edges.startTimeLabel, startTime);
				node.endTimeLabel = min(edges.endTimeLabel, endTimeLabel);
				queue.push(node);
			}
		}
	}
	return res;
}

template <typename weight_type>
void graph_v_of_v_with_time_span<weight_type>::print()
{
	std::cout << "graph_v_of_v_with_time_span_print:" << std::endl;
	int size = this->ADJs.size();
	for (int i = 0; i < size; i++)
	{
		std::cout << "Vertex " << i << " Adj List: ";
		int v_size = ADJs[i].size();
		for (int j = 0; j < v_size; j++)
		{
			std::cout << "<" << ADJs[i][j].vertex << "," << ADJs[i][j].weight << "," << ADJs[i][j].startTimeLabel << "," << ADJs[i][j].endTimeLabel << "> ";
		}
		std::cout << std::endl;
	}
	std::cout << "graph_v_of_v_with_time_span_print END" << std::endl;
}

int dijkstra(graph_v_of_v<int> graph, int u, int v)
{
	std::vector<int> dist(graph.size(), INT_MAX);
	std::vector<bool> visited(graph.size(), false);
	boost::heap::fibonacci_heap<pair<int, int>, boost::heap::compare<compare_pair>> queue;
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

int dijkstra_iterator(vector<graph_v_of_v<int>> list, int u, int v)
{
	int res = INT_MAX;
	for (graph_v_of_v<int> graph : list)
	{
		res = min(res, dijkstra(graph, u, v));
	}
	return res;
}