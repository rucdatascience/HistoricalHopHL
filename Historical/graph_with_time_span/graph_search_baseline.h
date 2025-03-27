#pragma once
#include <limits>
#include "Historical/graph_with_time_span/graph_with_time_span.h"
#include "Historical/graph_with_time_span/graph.h"
#include <boost/heap/fibonacci_heap.hpp>
#include <algorithm>

namespace experiment {
	namespace nonhop {
		struct compare_pair
		{
			bool operator()(const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) const
			{
				if (lhs.first == rhs.first)
				{
					return lhs.second > rhs.second;
				}
				return lhs.first > rhs.first;
			}
		};

		double search_shortest_path_in_period_time_naive(graph_with_time_span& graph, int u, int v, int startTime, int endTime) {
			double res = std::numeric_limits<double>::max();
			int N = graph.v_num;
			boost::heap::fibonacci_heap<std::pair<int, int>, boost::heap::compare<compare_pair>> queue;
			for (int queryTime = startTime; queryTime <= endTime; queryTime++)
			{
				std::vector<int> dist(N, std::numeric_limits<int>::max());
				queue.clear();
				dist[u] = 0;
				queue.push({ u, 0 });
				while (queue.size() > 0)
				{
					int vertexBase = queue.top().first;
					double currentDist = queue.top().second;
					queue.pop();
					if (vertexBase == v)
					{
						res = std::min(res, currentDist);
					}
					for (const auto& vertices : graph.ADJs[vertexBase])
					{
						int next = vertices.first;
						for (const auto& edge_info_time_span : vertices.second)
						{
							if (edge_info_time_span.startTimeLabel <= queryTime && edge_info_time_span.endTimeLabel >= queryTime)
							{
								int newDist = currentDist + edge_info_time_span.weight;
								if (newDist < dist[next])
								{
									dist[next] = newDist;
									queue.push({ next, newDist });
								}
								break;
							}
						}
					}
				}
			}
			return res;
		}

		int dijkstra(graph<int>& graph, int u, int v)
		{
			std::vector<int> dist(graph.size(), std::numeric_limits<int>::max());
			boost::heap::fibonacci_heap<std::pair<int, int>, boost::heap::compare<compare_pair>> queue;

			dist[u] = 0;
			queue.push({ u, 0 });
			int res = std::numeric_limits<int>::max();
			while (!queue.empty())
			{
				auto top = queue.top();
				int vertexBase = std::get<0>(top);
				int currentDist = std::get<1>(top);
				queue.pop();

				if (vertexBase == v)
				{
					res = std::min(res, currentDist);
				}
				for (const auto& edge : graph[vertexBase])
				{
					int next = edge.first;
					int weight = edge.second;
					int newDist = currentDist + weight;

					if (newDist < dist[next])
					{
						dist[next] = newDist;
						queue.push({ next, newDist });
					}
				}
			}
			return res;
		}
		
		int dijkstra_iterator(std::vector<graph<int>>& list, int u, int v, int ts, int te)
		{
			int res = INT_MAX;
			for (int i = 0; i < list.size(); i++) {
				if (i >= ts && i <= te) {
					res = std::min(res, dijkstra(list[i], u, v));
				}
			}
			return res;
		};
	}

	namespace hop {
		// k-constrianed
		struct compare_tuple
		{
			bool operator()(const std::tuple<int, int, int>& lhs, const std::tuple<int, int, int>& rhs) const
			{
				if (std::get<1>(lhs) == std::get<1>(rhs))
				{
					return std::get<2>(lhs) > std::get<2>(rhs);
				}
				return std::get<1>(lhs) > std::get<1>(rhs);
			}
		};

		double search_shortest_path_in_period_time_naive(graph_with_time_span& graph, int u, int v, int k, int startTime, int endTime)
		{
			double res = std::numeric_limits<double>::max();
			int N = graph.v_num;
			boost::heap::fibonacci_heap<std::tuple<int, int, int>, boost::heap::compare<compare_tuple>> queue;
			for (int queryTime = startTime; queryTime <= endTime; queryTime++)
			{
				std::vector<int> dist(N, std::numeric_limits<int>::max());
				std::vector<int> hop_list(N, std::numeric_limits<int>::max());
				queue.clear();
				dist[u] = 0;
				hop_list[u] = 0;
				queue.push({ u, 0, 0 });
				while (queue.size() > 0)
				{
					int hop = std::get<2>(queue.top());
					int vertexBase = std::get<0>(queue.top());
					int currentDist = std::get<1>(queue.top());
					queue.pop();
					if (vertexBase == v)
					{
						res = (currentDist < res) ? currentDist : res;
					}
					if (hop == k)
					{
						continue;
					}

					for (const auto& vertices : graph.ADJs[vertexBase])
					{
						int next = vertices.first;
						for (const auto& edge_info_time_span : vertices.second)
						{
							if (edge_info_time_span.startTimeLabel <= queryTime && edge_info_time_span.endTimeLabel >= queryTime)
							{
								int newDist = currentDist + edge_info_time_span.weight;
								if (newDist < dist[next] || (hop + 1) < hop_list[next])
								{
									dist[next] = newDist;
									hop_list[next] = hop + 1;
									queue.push({ next, newDist, hop + 1 });
								}
								break;
							}
						}
					}
				}
				// cout << "naive" << res << endl;
			}
			return res;
		};

		int dijkstra(graph<int>& graph, int u, int v, int k)
		{
			std::vector<int> dist(graph.size(), std::numeric_limits<int>::max());
			std::vector<int> hop_list(graph.size(), std::numeric_limits<int>::max());
			boost::heap::fibonacci_heap<std::tuple<int, int, int>, boost::heap::compare<compare_tuple>> queue;

			dist[u] = 0;
			hop_list[u] = 0;
			queue.push({ u, 0, 0 });
			int res = std::numeric_limits<int>::max();
			while (!queue.empty())
			{
				auto top = queue.top();
				int vertexBase = std::get<0>(top);
				int currentDist = std::get<1>(top);
				int hop = std::get<2>(top);
				queue.pop();

				if (vertexBase == v)
				{
					res = std::min(res, currentDist);
				}
				if (hop >= k)
				{
					continue;
				}

				for (const auto& edge : graph[vertexBase])
				{
					int next = edge.first;
					int weight = edge.second;

					int newDist = currentDist + weight;

					if (newDist < dist[next] || (hop + 1 < hop_list[next]))
					{
						dist[next] = newDist;
						hop_list[next] = hop + 1;
						queue.push({ next, newDist, hop + 1 });
					}
				}
			}

			return res;
		};

		int dijkstra_iterator(std::vector<graph<int>>& list, int u, int v, int ts, int te, int k)
		{
			int res = INT_MAX;
			for (int i = 0; i < list.size(); i++) {
				if (i >= ts && i <= te) {
					res = std::min(res, dijkstra(list[i], u, v, k));
				}
			}
			return res;
		};
	}

}