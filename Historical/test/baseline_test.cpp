using namespace std;
#include <chrono>
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"

int testBaseLineAndBaseline2()
{
    // generate a ramdom graph
    vector<graph_v_of_v<int>> graphs = generate_in_same_edge_different_weight(5, 10, 10, 30, 6, 9);
    graph_v_of_v_with_time_span<int> graphWithTimeSpan;
    // generate the time span
    for (int i = 0; i < graphs.size(); i++)
    {
        graphWithTimeSpan.add_graph_time(graphs[i], i);
    }
    graphWithTimeSpan.print();

    // dijkstra_iterator baseline 1
    vector<graph_v_of_v<int>> subsequence(graphs.begin() + 3, graphs.begin() + 7);
    auto start_time = std::chrono::high_resolution_clock::now();
    int resIterator = dijkstra_iterator(subsequence, 1, 4);
    auto endTime = std::chrono::high_resolution_clock::now();
    double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - start_time).count();

    // dfs to calculate the shortest path baseline 2
    start_time = std::chrono::high_resolution_clock::now();
    int res = graphWithTimeSpan.search_shortest_path_in_period_time(1, 4, 3, 7);
    endTime = std::chrono::high_resolution_clock::now();
    double runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - start_time).count();

    // print the result
    std::cout << res << ":" << resIterator << std::endl;
    std::cout << runningtime1 << ":" << runningtime2 << std::endl;
    return 0;
}