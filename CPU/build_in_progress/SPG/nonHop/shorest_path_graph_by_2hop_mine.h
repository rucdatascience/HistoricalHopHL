#pragma once
#include "CPU\build_in_progress\HL\HL4GST\nonHOP_maintain\nonHOP_maintain_two_hop_labels.h"
#include "CPU\graph_v_of_v\graph_v_of_v.h"
#include "Historical\graph_v_of_v\graph_v_of_v_with_time_span_non_hop_constrained.h"
using namespace nonHop;

// define the node in the queue
struct node_for_SPG_diffuse
{
    // 当前的节点
    int index;
    // 已经走了的距离
    weightTYPE disx;
    // 已经走了的距离
    vector<weightTYPE> hubElse;
    // hub
    vector<int> hub;
    // 最短路径的另一侧
    int target;
    // mode-0 直接查找 必须包含 不包含可以试做不合法
    // mode-1 也是直接查找，可以不包含 不包含要他的真正目的
    int mode;
    node_for_SPG_diffuse() {}
    node_for_SPG_diffuse(int _u, weightTYPE _dis, vector<weightTYPE> _hubElse, vector<int> _hub, int _target, int _mode) : index(_u), disx(_dis), hubElse(_hubElse), hub(_hub), target(_target), mode(_mode) {}
};

std::vector<std::pair<int, int>> nullPair = {{-1, -1}};
std::vector<std::pair<int, int>> samePair = {{0, 0}};
bool operator<(node_for_SPG_diffuse const &x, node_for_SPG_diffuse const &y)
{
    if (x.mode != y.mode)
    {
        return x.mode > y.mode;
    }
    return x.index > y.index; // < is the max-heap; > is the min heap
}

typedef typename boost::heap::fibonacci_heap<node_for_SPG_diffuse>::handle_type handle_t_for_SPG_diffuse; // pairing heap has a similar speed with fibonacci_heap here

template <typename WEIGHT_TYPE>
vector<int> process(int u, int v, int t_s, int t_e, two_hop_case_info &info, graph_v_of_v_with_time_span_non_hop_constrained<WEIGHT_TYPE> &graph_info)
{
    std::pair<std::vector<std::pair<int, int>>, std::vector<int>> dis2hub = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2_find_all_hub(info.L, u, v);
    // 查询u-v的最短距离 shortest_cost
    if (dis2hub.first == nullPair)
    {
        return dis2hub.second;
    }
    // 如果查询的是一个点 直接返回
    if (dis2hub.first == samePair)
    {
        return dis2hub.second;
    }
    int shortest_path_dis = dis2hub.first[0].first + dis2hub.first[0].second;
    boost::heap::fibonacci_heap<node_for_SPG_diffuse> Q;
    std::vector<handle_t_for_SPG_diffuse> Q_handles(graph_info.size());
    std::vector<int> status(graph_info.size(), 0);
    // SPG的结果集 res
    std::vector<int> res;
    int mark = u + v;
    // hub一定是结果
    for (const int &hub : dis2hub.second)
    {
        status[hub] = mark;
        res.push_back(hub);
    }
    if (dis2hub.second.size() == 1 && dis2hub.second[0] == u)
    {

        Q_handles[v] = Q.push(node_for_SPG_diffuse(v, 0, {shortest_path_dis}, dis2hub.second, u, 0));
    }
    else if (dis2hub.second.size() == 1 && dis2hub.second[0] == v)
    {
        Q_handles[u] = Q.push(node_for_SPG_diffuse(u, 0, {shortest_path_dis}, dis2hub.second, v, 0));
    }
    else
    {
        // commonhub is a high rank vertex
        vector<weightTYPE> costU;
        vector<weightTYPE> costV;
        for (int i = 0; i < dis2hub.first.size(); i++)
        {
            costU.push_back(dis2hub.first[i].first);
            costV.push_back(dis2hub.first[i].second);
        }
        Q_handles[u] = Q.push(node_for_SPG_diffuse(u, 0, costU, dis2hub.second, v, 1));
        Q_handles[v] = Q.push(node_for_SPG_diffuse(v, 0, costV, dis2hub.second, u, 1));
    }
    // 1.2 队列不为空 出队X并遍历邻居(TODO: 此处也可以并行 可以使用乐观锁)
    while (!Q.empty())
    {
        // (TODO: 并行)
        node_for_SPG_diffuse node = Q.top();
        Q.pop();
        res.push_back(node.index);
        int x = node.index;
        for (const pair<int, vector<EdgeInfo<weightTYPE>>> &pair : graph_info[x])
        {
            int v_id = pair.first;
            // 遍历邻居
            // 1. mode-0 直接二分找hub 如果符合则入堆 答案加入 否则直接剪枝
            // 2. mode-1 先二分查找hub
            // 2.1 如果二分hub找到了 并且符合 则入堆 答案加入
            // 2.2 如果没找到或者不符合, 则判断target 是否优先级比自己高。
            //     如果高则在当前的中二分target, 如果找到了target 并且符合
            //       则入堆，并且hub是target mode是0
            //     否则则在直接剪枝
            for (const EdgeInfo<weightTYPE> &edge : pair.second)
            {
                if (status[edge.vertex] == mark - node.target || status[edge.vertex] == mark)
                {
                    continue;
                }
                if (max(edge.startTimeLabel, t_s) <= min(edge.endTimeLabel, t_e))
                {
                    if (node.mode == 0)
                    {
                        // mode-0: 二分查找 hub
                        bool isValid = false;
                        for (int i = 0; i < node.hub.size() && edge.vertex > node.hub[i] && !isValid; i++)
                        {
                            isValid |= search_sorted_two_hop_label_specify_time_span_cost(info.L[edge.vertex], node.hub[i], shortest_path_dis - (node.disx + edge.weight), t_s, t_e);
                        }
                        if (isValid)
                        {
                            Q_handles[edge.vertex] = Q.push(node_for_SPG_diffuse(edge.vertex, node.disx + edge.weight, node.hubElse, dis2hub.second, node.target, 0));
                        }
                        else
                        {
                            status[edge.vertex] = mark;
                        }
                    }
                    else if (node.mode == 1)
                    {
                        // mode-1: 二分查找 hub
                        bool isValid = false;
                        for (int i = 0; i < node.hub.size() && edge.vertex > node.hub[i] && !isValid; i++)
                        {
                            isValid |= search_sorted_two_hop_label_specify_time_span_cost(info.L[edge.vertex], node.hub[i], node.hubElse[i] - (node.disx + edge.weight), t_s, t_e);
                        }
                        if (isValid)
                        {
                            Q_handles[edge.vertex] = Q.push(node_for_SPG_diffuse(edge.vertex, node.disx + edge.weight, node.hubElse, dis2hub.second, node.target, 1));
                            status[edge.vertex] = mark;
                        }
                        else
                        {
                            isValid = search_sorted_two_hop_label_specify_time_span_cost(info.L[max(edge.vertex, node.target)], min(edge.vertex, node.target), shortest_path_dis - (node.disx + edge.weight), t_s, t_e);
                            if (isValid)
                            {
                                Q_handles[edge.vertex] = Q.push(node_for_SPG_diffuse(edge.vertex, node.disx + edge.weight, node.hubElse, {node.target}, node.target, 0));
                                status[edge.vertex] = mark;
                            }
                            else
                            {
                                status[edge.vertex] = mark - node.target;
                            }
                        }
                    }
                }
            }
        }
    }

    return res;
}
