#include "Historical/hop_label/k_hop_constrained_two_hop_label_time_span.h"
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"

#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h"
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"

#include <chrono>
#include <limits.h>
int main(int argc, char const *argv[])
{
    boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};
    graph_v_of_v<int> instance_graph;
    int v_num = 500, e_num = 10000;
    int upper = 30, lower = 6;
    // int change_num = 100;
    instance_graph = graph_v_of_v_generate_random_graph<int>(v_num, e_num, lower, upper, 1, boost_random_time_seed);

    int res = dijkstra(instance_graph, 1, 3, 100);
    auto two_hop_label_with_time_span = graph_hop_constrained_two_hop_label_time_span(v_num);
    hop_constrained_case_info mm;
    mm.upper_k = 5;
    mm.max_bit_size = 6e9;
    mm.use_2M_prune = 1;
    mm.use_rank_prune = 1; // set true
    mm.use_2023WWW_generation = 0;
    mm.use_canonical_repair = 0;
    mm.max_run_time_seconds = 1e2;
    mm.thread_num = 1;
    two_hop_label_with_time_span.initiate_2_hop_label(instance_graph, mm);
    int res_2 = two_hop_label_with_time_span.query(1, 3, 100);
    std::cout << "res : " << res << "|res_2 : " << res_2 << std::endl;
    system("pause");
}
