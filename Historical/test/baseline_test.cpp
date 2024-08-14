using namespace std;
#include <chrono>
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"

int testBaseLineAndBaseline2()
{
    // generate a ramdom graph
    // graph_v_of_v_with_time_span<int> graphWithTimeSpan(5, 10, 30, 6);
    // vector<graph_v_of_v<int>> graphs = graphWithTimeSpan.graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(10, 6);
    // generate a larger ramdom graph
    graph_v_of_v_with_time_span<int> graphWithTimeSpan(5000, 100000, 30, 6);
    vector<graph_v_of_v<int>> graphs = graphWithTimeSpan.graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(100, 3);
    // print the graphs
    // for (graph_v_of_v<int> graph : graphs)
    // {
    //     graph.print();
    // }
    // graphWithTimeSpan.print();

    // // dijkstra_iterator baseline 1
    vector<graph_v_of_v<int>> subsequence(graphs.begin() + 3, graphs.begin() + 7 + 1);
    auto start_time = std::chrono::high_resolution_clock::now();
    int res_n_iterate_dijkstra = dijkstra_iterator(subsequence, 1, 4);
    auto endTime = std::chrono::high_resolution_clock::now();
    double runtime_n_iterate_dijkstra = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - start_time).count() / 1e9;
    // std::cout << res_n_iterate_dijkstra << std::endl;

    // dfs to calculate the shortest path baseline 2
    auto start_time_base_line_2 = std::chrono::high_resolution_clock::now();
    int res_base_line_with_span = graphWithTimeSpan.search_shortest_path_in_period_time_naive(1, 4, 3, 7);
    auto end_time_base_line_2 = std::chrono::high_resolution_clock::now();
    double runtime_base_line_with_span = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_base_line_2 - start_time_base_line_2).count() / 1e9;

    // // print the result
    std::cout << res_n_iterate_dijkstra  << ":" << res_base_line_with_span << std::endl;
    std::cout << runtime_n_iterate_dijkstra << ":" << runtime_base_line_with_span << std::endl;
    // debug
    // graphWithTimeSpan.print();
    // for (graph_v_of_v<int> graph : graphs)
    // {
    //     graph.print();
    // }

    // // // dijkstra_iterator baseline 1
    // int res_base_line_with_span_debug = graphWithTimeSpan.search_shortest_path_in_period_time_naive(1, 4, 3, 7);
    return 0;
}
int main()
{
    testBaseLineAndBaseline2();
    exit(0);
}