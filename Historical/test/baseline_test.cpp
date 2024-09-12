using namespace std;
#include <chrono>
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"

int testBaseLineAndBaseline2()
{
    // query param
    int source = 1, target = 4;
    int queryStartTime = 2, queryEndTime = 4;
    int k = 10;
    // generate a random graph
    // int v_num = 6, e_num = 6;
    // int upper = 30, lower = 6;
    // int change_num = 5, maintain_percent = 7;

    // generate a larger random graph
    int v_num = 20, e_num = 50;
    int upper = 100, lower = 1;
    int change_num = 5, maintain_percent = 3;

    // initialize the 2-hop label with time span
    hop_constrained_case_info mm;
    mm.upper_k = k;
    mm.max_bit_size = 6e9;
    mm.use_2M_prune = 1;
    mm.use_rank_prune = 1; // set true
    mm.use_2023WWW_generation = 0;
    mm.use_canonical_repair = 1;
    mm.max_run_time_seconds = 1e2;
    mm.thread_num = 1;

    graph_v_of_v_with_time_span<int> graphWithTimeSpan(v_num, e_num, upper, lower);
    graph_hop_constrained_two_hop_label_time_span two_hop_label_with_time_span(v_num);
    vector<graph_v_of_v<int>> graphs = graphWithTimeSpan.graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(change_num, maintain_percent);

    two_hop_label_with_time_span.process(graphs, mm);
    // print the graphs
    // for (graph_v_of_v<int> graph : graphs)
    // {
    //     graph.print();
    // }
    // graphWithTimeSpan.print();

    // dijkstra_iterator baseline 1
    vector<graph_v_of_v<int>> subsequence(graphs.begin() + queryStartTime, graphs.begin() + queryEndTime + 1);
    int res_n_iterate_dijkstra = dijkstra_iterator(subsequence, source, target, k);

    // dfs to calculate the shortest path baseline 2
    auto start_time_base_line_2 = std::chrono::high_resolution_clock::now();
    int res_base_line_with_span = graphWithTimeSpan.search_shortest_path_in_period_time_naive(source, target, k, queryStartTime, queryEndTime);
    auto end_time_base_line_2 = std::chrono::high_resolution_clock::now();
    double runtime_base_line_with_span = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_base_line_2 - start_time_base_line_2).count() / 1e9;
    std::cout << runtime_base_line_with_span << std::endl;
    std::cout << res_n_iterate_dijkstra << ":" << res_base_line_with_span << std::endl;
    // (test) the dsp using 2-hop labeling with time_span
    auto start_time_base_line_two_hop_label_test = std::chrono::high_resolution_clock::now();
    int res_two_hop_label = two_hop_label_with_time_span.test_query_method(source, target, queryStartTime, queryEndTime, k);
    auto end_time_base_line_two_hop_label_test = std::chrono::high_resolution_clock::now();
    double runtime_base_line_with_two_hop_label_test = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_base_line_two_hop_label_test - start_time_base_line_two_hop_label_test).count() / 1e9;
    std::cout << runtime_base_line_with_two_hop_label_test << std::endl;
    // // print the result
    std::cout << res_n_iterate_dijkstra << ":" << res_base_line_with_span << ":" << res_two_hop_label << std::endl;

    // debug
    two_hop_label_with_time_span.print_L_by_index(source);
    two_hop_label_with_time_span.print_L_by_index(target);
    for (auto &graph : subsequence)
    {
        hop_constrained_two_hop_labels_generation(graph, mm);
        mm.print_L_vk(source);
        mm.print_L_vk(target);
        mm.clear_labels();
    }
    graph_hop_constrained_two_hop_label_time_span two_hop_label_with_time_span_debug(v_num);
    two_hop_label_with_time_span_debug.process(graphs, mm);
    return 0;
}
int main()
{
    testBaseLineAndBaseline2();
    exit(0);
}