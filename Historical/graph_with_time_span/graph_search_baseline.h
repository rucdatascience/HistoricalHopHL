#pragma once
#include <limits>
#include "Historical/graph_with_time_span/graph_with_time_span.h"
#include "Historical/graph_with_time_span/graph.h"
#include <boost/heap/fibonacci_heap.hpp>
#include <algorithm>

namespace experiment {
	namespace nonhop {
		template <typename weight_type>
		struct compare_pair
		{
			bool operator()(const std::pair<int, weight_type>& lhs, const std::pair<int, weight_type>& rhs) const
			{
				if (lhs.first == rhs.first)
				{
					return lhs.second > rhs.second;
				}
				return lhs.first > rhs.first;
			}
		};
		template <typename weight_type>
		double search_shortest_path_in_period_time_naive(graph_with_time_span<weight_type>& graph, int u, int v, int startTime, int endTime) {
			double res = std::numeric_limits<double>::max();
			int N = graph.v_num;
			boost::heap::fibonacci_heap<std::pair<int, weight_type>, boost::heap::compare<compare_pair<weight_type>>> queue;
			for (int queryTime = startTime; queryTime <= endTime; queryTime++)
			{
				std::vector<weight_type> dist(N, std::numeric_limits<weight_type>::max());
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
								weight_type newDist = currentDist + edge_info_time_span.weight;
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

		template <typename weight_type>
		weight_type dijkstra(graph<weight_type>& graph, int u, int v)
		{
			std::vector<weight_type> dist(graph.size(), std::numeric_limits<weight_type>::max());
			boost::heap::fibonacci_heap<std::pair<int, weight_type>, boost::heap::compare<compare_pair<weight_type>>> queue;

			dist[u] = 0;
			queue.push({ u, 0 });
			int res = std::numeric_limits<int>::max();
			while (!queue.empty())
			{
				auto top = queue.top();
				int vertexBase = std::get<0>(top);
				weight_type currentDist = std::get<1>(top);
				queue.pop();

				if (vertexBase == v)
				{
					res = std::min(res, currentDist);
				}
				for (const auto& edge : graph[vertexBase])
				{
					int next = edge.first;
					weight_type weight = edge.second;
					weight_type newDist = currentDist + weight;

					if (newDist < dist[next])
					{
						dist[next] = newDist;
						queue.push({ next, newDist });
					}
				}
			}
			return res;
		}
		template <typename weight_type>
		int dijkstra_iterator(std::vector<graph<weight_type>>& list, int u, int v, int ts, int te)
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
		template <typename weight_type>
		struct compare_tuple
		{
			bool operator()(const std::tuple<int, weight_type, int>& lhs, const std::tuple<int, weight_type, int>& rhs) const
			{
				if (std::get<1>(lhs) == std::get<1>(rhs))
				{
					return std::get<2>(lhs) > std::get<2>(rhs);
				}
				return std::get<1>(lhs) > std::get<1>(rhs);
			}
		};
		template <typename weight_type>
		double search_shortest_path_in_period_time_naive(graph_with_time_span<weight_type>& graph, int u, int v, int k, int startTime, int endTime)
		{
			double res = std::numeric_limits<double>::max();
			int N = graph.v_num;
			boost::heap::fibonacci_heap<std::tuple<int, weight_type, int>, boost::heap::compare<compare_tuple<weight_type>>> queue;
			for (int queryTime = startTime; queryTime <= endTime; queryTime++)
			{
				std::vector<weight_type> dist(N, std::numeric_limits<weight_type>::max());
				std::vector<int> hop_list(N, std::numeric_limits<int>::max());
				queue.clear();
				dist[u] = 0;
				hop_list[u] = 0;
				queue.push({ u, 0, 0 });
				while (queue.size() > 0)
				{
					int hop = std::get<2>(queue.top());
					int vertexBase = std::get<0>(queue.top());
					weight_type currentDist = std::get<1>(queue.top());
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
								weight_type newDist = currentDist + edge_info_time_span.weight;
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
		template <typename weight_type>
		int dijkstra(graph<weight_type>& graph, int u, int v, int k)
		{
			std::vector<weight_type> dist(graph.size(), std::numeric_limits<weight_type>::max());
			std::vector<int> hop_list(graph.size(), std::numeric_limits<int>::max());
			boost::heap::fibonacci_heap<std::tuple<int, weight_type, int>, boost::heap::compare<compare_tuple<weight_type>>> queue;

			dist[u] = 0;
			hop_list[u] = 0;
			queue.push({ u, 0, 0 });
			int res = std::numeric_limits<int>::max();
			while (!queue.empty())
			{
				auto top = queue.top();
				int vertexBase = std::get<0>(top);
				weight_type currentDist = std::get<1>(top);
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
					weight_type weight = edge.second;

					weight_type newDist = currentDist + weight;

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

		template <typename weight_type>
		int dijkstra_iterator(std::vector<graph<weight_type>>& list, int u, int v, int ts, int te, int k)
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