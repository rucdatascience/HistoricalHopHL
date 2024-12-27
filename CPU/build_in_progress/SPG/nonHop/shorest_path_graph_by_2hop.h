#include "CPU\build_in_progress\HL\HL4GST\nonHOP_maintain\nonHOP_maintain_two_hop_labels.h"
#include "CPU\graph_v_of_v\graph_v_of_v.h"
#include "Historical\graph_v_of_v\graph_v_of_v_with_time_span_non_hop_constrained.h"
using namespace nonHop;

template <typename WEIGHT_TYPE>
vector<int> process(int u, int v, two_hop_case_info &info, graph_v_of_v_with_time_span_non_hop_constrained<WEIGHT_TYPE> &graph_info)
{
    // SPG的结果集 res
    // 查询u-v的最短距离 shorest_cost
    // (TODO: 并行) 1. 处理u和v的每一个label
    // 1.入队u,v 入队的信息是info(vertex, cost_else)
    // 2.队列不为空
}

void diffuse(int u, int target, int cost, graph_v_of_v<weightTYPE> &graph, two_hop_case_info &info, ThreadPool &pool, std::vector<std::future<vector<int>>> &results_dynamic)
{

}

void search_shorest_path_graph_partial(int u, int v, int cost)
{
}
