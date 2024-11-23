using namespace std;
#include <chrono>
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span_non_hop_constrained.h"
using namespace nonHop;

int testBaseLineAndBaseline2()
{
    int iterator = 10;
    int index = 0;
    timer.is_debug = false;
    while (index < iterator)
    {
        timer.begin_timing();
        ++index;
        try
        {
            // query param
            int source = 7, target = 9;
            int queryStartTime = 1, queryEndTime = 9;
            // generate a random graph
            // int v_num = 10, e_num = 20;
            // int upper = 20, lower = 1;
            // int change_num = 2, decrease_time = 0, increase_time = 5;
            // float change_ratio = 0.3;
            // generate a larger random graph
            int v_num = 5, e_num = 10;
            int upper = 100, lower = 50;
            int change_num = 10, decrease_time = 3, increase_time = 3;
            float change_ratio = 0.3;

            graph_v_of_v_with_time_span_non_hop_constrained<int> graph_with_time_span_non_hop_constrained;
            vector<graph_v_of_v<int>> graphs;

            // initialize the 2-hop label with time span
            two_hop_case_info mm;
            mm.max_labal_byte_size = 6e9;
            mm.max_run_time_seconds = 1e4;
            mm.use_2M_prune = 1;
            mm.use_rank_prune = 1;
            mm.use_canonical_repair = 1;
            mm.thread_num = 1;
            mm.source = source;
            mm.target = target;
            mm.t_s = queryStartTime;
            mm.t_e = queryEndTime;

            two_hop_case_info mm2021;
            mm.max_labal_byte_size = 6e9;
            mm.max_run_time_seconds = 1e4;
            mm.use_2M_prune = 1;
            mm.use_rank_prune = 1;
            mm.use_canonical_repair = 1;
            mm.thread_num = 5;
            mm.source = source;
            mm.target = target;
            mm.t_s = queryStartTime;
            mm.t_e = queryEndTime;

            bool use_save_read = true;
            bool use_2_hop_label = true;

            if (use_save_read)
            {
                graph_with_time_span_non_hop_constrained = graph_v_of_v_with_time_span_non_hop_constrained<int>();
                graphs = graph_with_time_span_non_hop_constrained.txt_read("time-graph-nonhop.txt", mm, mm2021);
                // graphs = graph_with_time_span.txt_read("time-graph-2024-10-09-1729.txt", mm);
                if (graph_with_time_span_non_hop_constrained.size() < source || graph_with_time_span_non_hop_constrained.size() < target)
                {
                    cout << "vertex is out of range" << endl;
                    return 0;
                }
            }
            else
            {
                graph_with_time_span_non_hop_constrained = graph_v_of_v_with_time_span_non_hop_constrained<int>(v_num, e_num, upper, lower);
                graphs = graph_with_time_span_non_hop_constrained.graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(change_num, decrease_time, increase_time, change_ratio, mm, mm2021);
                graph_with_time_span_non_hop_constrained.txt_save("time-graph-nonhop.txt");
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
            int res_n_iterate_dijkstra = dijkstra_iterator(subsequence, source, target);
            // int temp = graph_with_time_span_non_hop_constrained.search_shortest_path_in_period_time_naive(2, 3, 6, 6);
            // dfs to calculate the shortest path baseline 2
            auto start_time_base_line_2 = std::chrono::high_resolution_clock::now();
            int res_base_line_with_span = graph_with_time_span_non_hop_constrained.search_shortest_path_in_period_time_naive(source, target, queryStartTime, queryEndTime);
            auto end_time_base_line_2 = std::chrono::high_resolution_clock::now();
            double runtime_base_line_with_span = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_base_line_2 - start_time_base_line_2).count() / 1e9;
            if (use_2_hop_label)
            {
                // mm.print_L();
                // mm2021.print_L();
                auto start_time_2_hop_label = std::chrono::high_resolution_clock::now();
                int res = mm.query(source, target, queryStartTime, queryEndTime);
                auto end_time_2_hop_label = std::chrono::high_resolution_clock::now();
                double runtime_2_hop_label_with_span = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_2_hop_label - start_time_2_hop_label).count() / 1e9;

                auto start_time_2_hop_label_2021 = std::chrono::high_resolution_clock::now();
                int res_2021 = mm2021.query(source, target, queryStartTime, queryEndTime);
                auto end_time_2_hop_label_2021 = std::chrono::high_resolution_clock::now();
                double runtime_2_hop_label_with_span_2021 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time_2_hop_label_2021 - start_time_2_hop_label_2021).count() / 1e9;

                std::cout << "query time is " << runtime_2_hop_label_with_span << std::endl;
                std::cout << "2021 query time is " << runtime_2_hop_label_with_span_2021 << std::endl;
                std::cout << res_n_iterate_dijkstra << ":" << res_base_line_with_span << ":" << res << ":" << res_2021 << std::endl;

                if (!(res_n_iterate_dijkstra == res_base_line_with_span && res == res_2021 && res_n_iterate_dijkstra == res))
                {
                    throw "error result .please check the algorithm";
                }
            }
            else
            {
                std::cout << "query time of graph with time_span :" << runtime_base_line_with_span << std::endl;
            }
        }
        catch (const char *c)
        {
            std::cerr << "Error: " << c << std::endl;
            return 1;
        }
    }
    return 0;
}
int main()
{
    testBaseLineAndBaseline2();
    exit(0);
}