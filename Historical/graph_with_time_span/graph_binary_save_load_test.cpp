//#include "Historical/graph_with_time_span/graph_with_time_span.h"
//#include "Historical/graph_with_time_span/graph.h"
//#include "Historical/graph_with_time_span/graph_search_baseline.h"
//
//int main() {
//	int M = 100;
//	experiment::graph<int> instance_graph(10);
//	experiment::graph_with_time_span<int> instance_graph_list;
//	instance_graph.add_edge(0, 1, 5);
//	instance_graph.add_edge(0, 2, M);
//	instance_graph.add_edge(0, 4, 1);
//	instance_graph.add_edge(0, 5, 7);
//	instance_graph.add_edge(0, 6, 3);
//	instance_graph.add_edge(0, 8, M);
//	// 1-x
//	instance_graph.add_edge(1, 2, M);
//	instance_graph.add_edge(1, 3, 15);
//	instance_graph.add_edge(1, 5, 2);
//	// 2-x
//	instance_graph.add_edge(2, 3, M);
//	instance_graph.add_edge(2, 4, M);
//	// 3-x
//	instance_graph.add_edge(3, 4, 3);
//	instance_graph.add_edge(3, 7, M);
//	// 4-x
//	instance_graph.add_edge(4, 9, M);
//	// 5-x
//	instance_graph.add_edge(5, 6, 3);
//	instance_graph.add_edge(5, 7, M);
//	// 6-x
//	instance_graph.add_edge(6, 8, M);
//
//	instance_graph_list.add_graph_time(instance_graph, 0);
//
//	instance_graph.print();
//	instance_graph_list.print();
//
//
//	std::ofstream FILE("test_binary_save_graph", std::ios::out | std::ifstream::binary);
//	std::ofstream FILE_LIST("test_binary_save_graph_list", std::ios::out | std::ifstream::binary);
//	std::ifstream FILE1("test_binary_save_graph", std::ios::in | std::ifstream::binary);
//	std::ifstream FILE_LIST1("test_binary_save_graph_list", std::ios::in | std::ifstream::binary);
//
//	instance_graph.serialize(FILE);
//	instance_graph_list.serialize(FILE_LIST);
//	FILE.close();
//	FILE_LIST.close();
//	experiment::graph<int> instance_graph_new;
//	experiment::graph_with_time_span<int> instance_graph_list_new;
//	instance_graph_new.deserialize(FILE1);
//	instance_graph_list_new.deserialize(FILE_LIST1);
//
//	std::cout << "================" << std::endl;
//	instance_graph_new.print();
//	instance_graph_list_new.print();
//
//	double res = experiment::hop::search_shortest_path_in_period_time_naive<int>(instance_graph_list_new, 0, 4, 3, 0, 0);
//	printf("from %d to %d weight is %lf\n", 0, 4, res);
//	std::cout << "finished" << std::endl;
//}