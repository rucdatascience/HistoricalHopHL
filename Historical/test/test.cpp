#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"

#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h"
#include "CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h"

#include <chrono>
#include <limits.h>
int main(int argc, char const *argv[])
{
    int M = 100;
    graph_v_of_v<int> instance_graph(10);
    // 0-x
    instance_graph.add_edge(0, 1, 5);
    instance_graph.add_edge(0, 2, M);
    instance_graph.add_edge(0, 4, 1);
    instance_graph.add_edge(0, 5, 7);
    instance_graph.add_edge(0, 6, 3);
    instance_graph.add_edge(0, 8, M);
    // 1-x
    instance_graph.add_edge(1, 2, M);
    instance_graph.add_edge(1, 3, 15);
    instance_graph.add_edge(1, 5, 2);
    // 2-x
    instance_graph.add_edge(2, 3, M);
    instance_graph.add_edge(2, 4, M);
    // 3-x
    instance_graph.add_edge(3, 4, 3);
    instance_graph.add_edge(3, 7, M);
    // 4-x
    instance_graph.add_edge(4, 9, M);
    // 5-x
    instance_graph.add_edge(5, 6, 3);
    instance_graph.add_edge(5, 7, M);
    // 6-x
    instance_graph.add_edge(6, 8, M);
    auto two_hop_label_with_time_span = graph_hop_constrained_two_hop_label_time_span(10);
    hop_constrained_case_info mm;
    mm.upper_k = 5;
    mm.max_bit_size = 6e9;
    mm.use_2M_prune = 1;
    mm.use_rank_prune = 1; // set true
    mm.use_2023WWW_generation = 0;
    mm.use_canonical_repair = 0;
    mm.max_run_time_seconds = 100;
    mm.thread_num = 1;
    two_hop_label_with_time_span.initiate_2_hop_label(instance_graph, mm);
    two_hop_label_with_time_span.print_L();
    two_hop_label_with_time_span.add_new_edge_or_weight_decrease(3, 5, 6);
    two_hop_label_with_time_span.add_new_edge_or_weight_decrease(4, 6, 3);
    two_hop_label_with_time_span.add_time();
    std::cout << "after" << std::endl;
    two_hop_label_with_time_span.print_L();
}
