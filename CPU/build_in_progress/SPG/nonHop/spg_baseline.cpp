//#include <iostream>
//#include <vector>
//#include <queue>
//#include <climits>
//#include <unordered_map>
//#include <algorithm>
//#include <set>
//
//using namespace std;
//
//// 使用邻接表表示图
//typedef pair<int, int> Edge;  // (目标节点, 权重)
//typedef unordered_map<int, vector<Edge>> Graph;
//
//// 定义全局变量
//unordered_map<int, int> dist;  // 存储最短距离
//unordered_map<int, vector<int>> predecessors;  // 存储每个节点的前驱节点列表
//unordered_map<int, bool> visited;  // 存储节点是否被访问过
//
//// 读取图的函数
//Graph read_graph() {
//	Graph graph;
//	// 可以根据需要修改图的输入方式，这里是硬编码示例
//	graph[1] = { {2, 3}, {4, 2}, {5, 1}, {6, 1} };
//	graph[2] = { {1, 3}, {3, 4}, {8, 1}, {9, 2} };
//	graph[3] = { {2, 4}, {4, 1}, {12, 1}, {13, 2} };
//	graph[4] = { {1, 2}, {3, 1}, {14, 2} };
//	graph[5] = { {1, 1}, {6, 1}, {14, 2} };
//	graph[6] = { {1, 1}, {5, 1}, {7, 1} };
//	graph[7] = { {6, 1}, {8, 2} };
//	graph[8] = { {2, 1}, {7, 2}, {9, 1} };
//	graph[9] = { {2, 2}, {8, 1}, {10, 3} };
//	graph[10] = { {9, 3}, {11, 1} };
//	graph[11] = { {10, 1}, {12, 1} };
//	graph[12] = { {3, 1}, {11, 1} };
//	graph[13] = { {3, 2}, {14, 2} };
//	graph[14] = { {4, 2}, {5, 2}, {13, 2} };
//
//
//	return graph;
//}
//
//// Dijkstra算法：计算最短路径并记录多个前驱节点
//void dijkstra_all_paths(const Graph& graph, int start, int end) {
//	priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
//	pq.push({ 0, start });
//	dist[start] = 0;
//
//	// 初始化visited数组，所有节点初始状态为未访问
//	for (const auto& node : graph) {
//		visited[node.first] = false;  // 将每个节点的访问状态设为false
//	}
//
//
//	while (!pq.empty()) {
//		int current_dist = pq.top().first;
//		int current_node = pq.top().second;
//		pq.pop();
//		//标记为已访问
//		visited[current_node] = true;
//
//		//目标点pop后就可以退出
//		if (current_node == end) {
//			return;
//		}
//
//
//		// 遍历所有邻接节点
//		for (const auto& edge : graph.at(current_node)) {
//			int neighbor = edge.first;
//			int weight = edge.second;
//			int new_dist = current_dist + weight;
//
//			if (visited[neighbor] == false) {//对于没有访问的邻居
//				if (new_dist < dist[neighbor]) {
//					dist[neighbor] = new_dist;
//					predecessors[neighbor].clear();  // 如果发现更短路径，清除之前的前驱
//					predecessors[neighbor].push_back(current_node);
//					pq.push({ new_dist, neighbor });
//				}
//				else if (new_dist == dist[neighbor]) {
//					predecessors[neighbor].push_back(current_node);  // 记录多个前驱节点
//				}
//			}
//
//		}
//	}
//}
//
//// 打印predecessors的内容
//void print_predecessors() {
//	cout << "predecessors:\n";
//	for (const auto& entry : predecessors) {
//		cout << "节点 " << entry.first << ": ";
//		for (int pre : entry.second) {
//			cout << pre << " ";
//		}
//		cout << endl;
//	}
//}
//
//// 恢复从源点到目标点的所有最短路径
//void find_all_paths(int node, vector<int>& path, vector<vector<int>>& all_paths) {
//	if (predecessors[node].empty()) {
//		path.push_back(node);
//		reverse(path.begin(), path.end());
//		all_paths.push_back(path);  // 存储一条路径
//		reverse(path.begin(), path.end());
//		path.pop_back();  // 回溯，移除最后一个节点
//		return;
//	}
//
//	// 对于每个前驱节点递归地找到路径
//	path.push_back(node);  // 当前节点加入路径
//	for (int pre : predecessors[node]) {
//		find_all_paths(pre, path, all_paths);
//	}
//	path.pop_back();  // 回溯，移除当前节点
//}
//
//// 构建从源点到目标点的最短路径子图（邻接表）
//// 确保每条边同时出现在两个节点的邻接列表中
//Graph build_subgraph(const Graph& graph, const vector<vector<int>>& all_paths) {
//	Graph subgraph;
//
//	// 遍历所有路径，将路径中的边加入子图
//	for (const auto& path : all_paths) {
//		for (size_t i = 0; i < path.size() - 1; ++i) {
//			int u = path[i];
//			int v = path[i + 1];
//			int weight = -1;
//
//			// 查找这条边的权重
//			for (const auto& edge : graph.at(u)) {
//				if (edge.first == v) {
//					weight = edge.second;
//					break;
//				}
//			}
//			if (weight != -1) {
//				// 将边加入起点的邻接表
//				auto it_u = find_if(subgraph[u].begin(), subgraph[u].end(),
//					[v](const Edge& e) { return e.first == v; });
//				if (it_u == subgraph[u].end()) {
//					subgraph[u].push_back({ v, weight });
//				}
//
//				// 将边加入终点的邻接表（无向图特性）
//				auto it_v = find_if(subgraph[v].begin(), subgraph[v].end(),
//					[u](const Edge& e) { return e.first == u; });
//				if (it_v == subgraph[v].end()) {
//					subgraph[v].push_back({ u, weight });
//				}
//			}
//		}
//	}
//	return subgraph;
//}
//
//// 打印子图的邻接表
//void print_subgraph(const Graph& subgraph) {
//	cout << "子图的邻接表：\n";
//	for (const auto& node : subgraph) {
//		cout << node.first << ": ";
//		for (const auto& edge : node.second) {
//			cout << "(" << edge.first << ", " << edge.second << ") ";
//		}
//		cout << endl;
//	}
//}
//
////int main() {
////    // 读取图
////    Graph graph = read_graph();
////
////    // 用户输入源点和目标点
////    int s = 9; // 示例源点
////    int t = 4; // 示例目标点
////
////    // 初始化dist
////    for (const auto& node : graph) {
////        dist[node.first] = INT_MAX;
////    }
////
////    // 运行Dijkstra算法
////    dijkstra_all_paths(graph, s, t);
////
////    //打印前驱信息
////    print_predecessors();
////
////    // 输出最短路径长度
////    cout << "从 " << s << " 到 " << t << " 的最短路径长度为: " << dist[t] << endl;
////
////    // 输出从s到t的所有最短路径
////    vector<vector<int>> all_paths;
////    vector<int> path;
////    find_all_paths(t, path, all_paths);
////
////    cout << "从 " << s << " 到 " << t << " 的所有最短路径：\n";
////    for (const auto& p : all_paths) {
////        for (int node : p) {
////            cout << node << " ";
////        }
////        cout << endl;
////    }
////
////    // 构建并输出从s到t的最短路径子图
////    Graph subgraph = build_subgraph(graph, all_paths);
////    print_subgraph(subgraph);
////
////    return 0;
////}