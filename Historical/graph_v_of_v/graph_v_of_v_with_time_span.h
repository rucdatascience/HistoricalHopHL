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
	EdgeInfo() : vertex(-1) {}
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
	std::vector<std::vector<EdgeInfo<weight_type> *>> ADJs;
	std::vector<std::vector<EdgeInfo<weight_type> *>> preLabels;
	/*constructors*/
	graph_v_of_v_with_time_span() {}
	graph_v_of_v_with_time_span(int n) : ADJs(n), preLabels(n, std::vector<EdgeInfo<weight_type> *>(n, nullptr))
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
	~graph_v_of_v_with_time_span()
	{
		for (auto &adjList : ADJs)
		{
			for (auto edge : adjList)
			{
				delete edge;
			}
		}
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
	if (preLabels[e1][e2] != nullptr && preLabels[e1][e2]->endTimeLabel == time - 1 && preLabels[e1][e2]->weight == ec)
	{
		preLabels[e1][e2]->endTimeLabel = time;
		preLabels[e2][e1]->endTimeLabel = time;
	}
	else
	{
		EdgeInfo<weight_type> *edgeInfoE1 = new EdgeInfo<weight_type>(e2, ec, time, time);
		this->ADJs[e1].push_back(edgeInfoE1);
		EdgeInfo<weight_type> *edgeInfoE2 = new EdgeInfo<weight_type>(e1, ec, time, time);
		this->ADJs[e2].push_back(edgeInfoE2);
		preLabels[e1][e2] = edgeInfoE1;
		preLabels[e2][e1] = edgeInfoE2;
	}
	// if (this->ADJs[e1].size() != index && this->ADJs[e1][index].vertex == e2 && this->ADJs[e1][index].weight == ec && this->ADJs[e1][index].endTimeLabel == time - 1)
	// {
	// 	this->ADJs[e1][index].endTimeLabel = time;
	// }
	// else if (this->ADJs[e1].size() == index || this->ADJs[e1][index].vertex != e2)
	// {
	// 	EdgeInfo<weight_type> edgeInfo(e2, ec, time, time);
	// 	this->ADJs[e1].insert(ADJs[e1].begin() + index, edgeInfo);
	// }
	// else
	// {
	// 	EdgeInfo<weight_type> edgeInfo(e2, ec, time, time);
	// 	this->ADJs[e1].insert(this->ADJs[e1].begin() + index + 1, edgeInfo);
	// }
	// index = sorted_vector_binary_operations(this->ADJs[e2], e1);
	// if (this->ADJs[e2].size() != index && this->ADJs[e2][index].vertex == e1 && this->ADJs[e2][index].weight == ec && this->ADJs[e2][index].endTimeLabel == time - 1)
	// {
	// 	this->ADJs[e2][index].endTimeLabel = time;
	// }
	// else if (this->ADJs[e2].size() == index || this->ADJs[e2][index].vertex != e1)
	// {
	// 	EdgeInfo<weight_type> edgeInfo(e1, ec, time, time);
	// 	this->ADJs[e2].insert(ADJs[e2].begin() + index, edgeInfo);
	// }
	// else
	// {
	// 	EdgeInfo<weight_type> edgeInfo(e1, ec, time, time);
	// 	this->ADJs[e2].insert(ADJs[e2].begin() + index + 1, edgeInfo);
	// }
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
	int N = this->ADJs.size();
	/**
	 * i -> vertex
	 * j -> index of time_span
	 * pair->time_span
	 */
	vector<vector<pair<int, int>>> visited(N);
	/* a class contains information about destination vertex, hop count, and cost in the priority queue */
	boost::heap::fibonacci_heap<NodeWithTimeSpan<weight_type>, boost::heap::compare<compare_node>> queue;
	NodeWithTimeSpan node(u, 0, startTime, endTime);
	queue.push(node);
	int res = -1;

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
		/* the edges with uncalculated results */
		vector<pair<int, int>> temp;
		vector<pair<int, int>> iteraotr_items;
		int pre = startTimeLabel, after = endTimeLabel;
		int index = startTimeLabel;
		int i = 0;
		int len = visited[vertexBase].size();

		while (i < len && visited[vertexBase][i].second < startTimeLabel)
		{
			temp.push_back(visited[vertexBase][i]);
			++i;
		}
		if (i < len && startTimeLabel >= visited[vertexBase][i].first && endTimeLabel <= visited[vertexBase][i].second)
		{
			continue;
		}
		while (i < len && visited[vertexBase][i].first <= endTimeLabel)
		{
			if (index < visited[vertexBase][i].first)
			{
				iteraotr_items.push_back({index, visited[vertexBase][i].first - 1});
			}
			index = visited[vertexBase][i].second + 1;
			pre = min(pre, visited[vertexBase][i].first);
			after = max(after, visited[vertexBase][i].second);
			++i;
		}
		if (i == len && index <= endTimeLabel)
		{
			iteraotr_items.push_back({index, endTimeLabel});
		}
		else if (i < len)
		{
			iteraotr_items.push_back({index, visited[vertexBase][i].first});
		}
		temp.push_back({pre, after});
		while (i < len)
		{
			temp.push_back(visited[vertexBase][i]);
			i++;
		}
		visited[vertexBase] = temp;
		for (const auto &edges : this->ADJs[vertexBase])
		{
			for (const auto &dfs_time_span : iteraotr_items)
			{
				if (!(edges->startTimeLabel > dfs_time_span.second || edges->endTimeLabel < dfs_time_span.first))
				{
					node.vertex = edges->vertex;
					node.weight = curWeight + edges->weight;
					node.startTimeLabel = max(edges->startTimeLabel, dfs_time_span.first);
					node.endTimeLabel = min(edges->endTimeLabel, dfs_time_span.second);
					queue.push(node);
				}
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
			std::cout << "<" << ADJs[i][j]->vertex << "," << ADJs[i][j]->weight << "," << ADJs[i][j]->startTimeLabel << "," << ADJs[i][j]->endTimeLabel << "> ";
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