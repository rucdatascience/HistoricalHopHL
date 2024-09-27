using namespace std;
#include <chrono>
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels.h"

int testBaseLineAndBaseline2()
{

    try
    {
        // query param
        int source = 900, target = 999;
        int queryStartTime = 1, queryEndTime = 1;
        int k = 10;
        // generate a random graph
        // int v_num = 6, e_num = 6;
        // int upper = 30, lower = 6;
        // int change_num = 5, maintain_percent = 7;

        // generate a larger random graph
        int v_num = 10000, e_num = 20000;
        int upper = 20, lower = 1;
        int change_num = 1, decrease_time = 30;
        float change_ratio = 0.3;

        // initialize the 2-hop label with time span
        hop_constrained_case_info mm;
        mm.upper_k = k;
        mm.max_bit_size = 6e9;
        mm.use_2M_prune = 1;
        mm.use_rank_prune = 1; // set true
        mm.use_2023WWW_generation = 0;
        mm.use_canonical_repair = 0;
        mm.max_run_time_seconds = 1e2;
        mm.thread_num = 10;

        bool use_save_read = false;
        bool use_2_hop_label = true;
        graph_v_of_v_with_time_span<int> graph_with_time_span;
        vector<graph_v_of_v<int>> graphs;
        if (use_save_read)
        {
            graph_with_time_span = graph_v_of_v_with_time_span<int>();
            graphs = graph_with_time_span.txt_read("time-graph.txt", mm);
            if (graph_with_time_span.size() < source || graph_with_time_span.size() < target)
            {
                cout << "vertex is out of range" << endl;
                return 0;
            }
        }
        else
        {
            initialize_global_values_dynamic_hop_constrained(v_num, mm.thread_num, mm.upper_k);
            graph_with_time_span = graph_v_of_v_with_time_span<int>(v_num, e_num, upper, lower);
            graphs = graph_with_time_span.graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(change_num, decrease_time, change_ratio, mm);
            // graph_with_time_span.txt_save("time-graph.txt");
        }
        // print the graphs
        // for (graph_v_of_v<int> graph : graphs)
        // {
        //     graph.print();
        // }
        // graphWithTimeSpan.print();

        // dijkstra_iterator baseline 1
        if (queryStartTime < 0 || queryEndTime < queryStartTime || queryEndTime > change_num)
        {
            std::cerr << "error query time" << std::endl;
            return 1;
        }
        vector<graph_v_of_v<int>> subsequence(graphs.begin() + queryStartTime, graphs.begin() + queryEndTime + 1);
        int res_n_iterate_dijkstra = dijkstra_iterator(subsequence, source, target, k);

        // dfs to calculate the shortest path baseline 2
        auto start_time_base_line_2 = std::chrono::high_resolution_clock::now();
        int res_base_line_with_span = graph_with_time_span.search_shortest_path_in_period_time_naive(source, target, k, queryStartTime, queryEndTime);
        auto end_time_base_line_2 = std::chrono::high_resolution_clock::now();
        double runtime_base_line_with_span = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_base_line_2 - start_time_base_line_2).count() / 1e9;
        std::cout << runtime_base_line_with_span << std::endl;
        std::cout << res_n_iterate_dijkstra << ":" << res_base_line_with_span << std::endl;
        if (use_2_hop_label)
        {
            auto start_time_base_line_2 = std::chrono::high_resolution_clock::now();
            // mm.print_L();
            int res = mm.query(source, target, queryStartTime, queryEndTime, k);
            std::cout << res_n_iterate_dijkstra << ":" << res_base_line_with_span << ":" << res << std::endl;
        }
        return 0;
    }
    catch (...)
    {
        std::cerr << "Error: "<< std::endl;
        return 1;
    }
}
int main()
{
    testBaseLineAndBaseline2();
    exit(0);
}