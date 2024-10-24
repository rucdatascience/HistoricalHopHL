using namespace std;
#include <chrono>
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels.h"
int testBaseLineAndBaseline2()
{
    int iterator = 10;
    int index = 0;
    while (index < iterator)
    {
        ++index;
        try
        {
            // query param
            int source = 90, target = 99;
            int queryStartTime = 4, queryEndTime = 8;
            int k = 10;
            // generate a random graph
            // int v_num = 10, e_num = 20;
            // int upper = 20, lower = 1;
            // int change_num = 2, decrease_time = 0, increase_time = 5;
            // float change_ratio = 0.3;
            // generate a larger random graph
            int v_num = 1000, e_num = 10000;
            int upper = 100, lower = 80;
            int change_num = 10, decrease_time = 10, increase_time = 10;
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
            mm.source = source;
            mm.target = target;
            mm.t_s = queryStartTime;
            mm.t_e = queryEndTime;

            hop_constrained_case_info mm2021;
            mm2021.upper_k = k;
            mm2021.max_bit_size = 6e9;
            mm2021.use_2M_prune = 1;
            mm2021.use_rank_prune = 1; // set true
            mm2021.use_2023WWW_generation = 0;
            mm2021.use_canonical_repair = 0;
            mm2021.max_run_time_seconds = 1e2;
            mm2021.thread_num = 10;
            mm2021.source = source;
            mm2021.target = target;
            mm2021.t_s = queryStartTime;
            mm2021.t_e = queryEndTime;

            bool use_save_read = false;
            bool use_2_hop_label = true;
            graph_v_of_v_with_time_span<int> graph_with_time_span;
            vector<graph_v_of_v<int>> graphs;
            mm.begin_timing();
            mm2021.begin_timing();
            if (use_save_read)
            {
                graph_with_time_span = graph_v_of_v_with_time_span<int>();
                graphs = graph_with_time_span.txt_read("time-graph.txt", mm, mm2021);
                // graphs = graph_with_time_span.txt_read("time-graph-2024-10-09-1729.txt", mm);
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
                mm.mark_time("initialize_global_values_dynamic_hop_constrained");
                graphs = graph_with_time_span.graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(change_num, decrease_time, increase_time, change_ratio, mm, mm2021);
                graph_with_time_span.txt_save("time-graph.txt");
            }
            std::cout << "maintain runtime is " << mm.get_maintain_time() << std::endl;
            std::cout << "2021 maintain runtime is " << mm2021.get_maintain_time() << std::endl;

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
            std::cout << "query time of graph with time_span :" << runtime_base_line_with_span << std::endl;
            std::cout << res_n_iterate_dijkstra << ":" << res_base_line_with_span << std::endl;
            if (use_2_hop_label)
            {
                auto start_time_2_hop_label = std::chrono::high_resolution_clock::now();
                int res = mm.query(source, target, queryStartTime, queryEndTime, k);
                auto end_time_2_hop_label = std::chrono::high_resolution_clock::now();
                double runtime_2_hop_label_with_span = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_2_hop_label - start_time_2_hop_label).count() / 1e9;

                auto start_time_2_hop_label_2021 = std::chrono::high_resolution_clock::now();
                int res_2021 = mm2021.query(source, target, queryStartTime, queryEndTime, k);
                auto end_time_2_hop_label_2021 = std::chrono::high_resolution_clock::now();
                double runtime_2_hop_label_with_span_2021 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_2_hop_label_2021 - start_time_2_hop_label_2021).count() / 1e9;

                std::cout << "query time is " << runtime_2_hop_label_with_span << std::endl;
                std::cout << "2021 query time is " << runtime_2_hop_label_with_span_2021 << std::endl;
                std::cout << res_n_iterate_dijkstra << ":" << res_base_line_with_span << ":" << res << ":" << res_2021 << std::endl;
            }
        }
        catch (const char *c)
        {
            std::cerr << "Error: " << c << std::endl;
        }
    }
    return 0;
}
int main()
{
    testBaseLineAndBaseline2();
    exit(0);
}