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
    // 到这个hub剩余的距离
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
boost::heap::fibonacci_heap<node_for_SPG_diffuse> Q_SPG;
std::vector<handle_t_for_SPG_diffuse> Q_SPG_handles;
std::vector<int> status;

template <typename WEIGHT_TYPE>
void preInit(graph_v_of_v_with_time_span_non_hop_constrained<WEIGHT_TYPE> &graph_info)
{
    boost::heap::fibonacci_heap<node_for_SPG_diffuse>().swap(Q_SPG);
    std::vector<handle_t_for_SPG_diffuse>(graph_info.size()).swap(Q_SPG_handles);
    std::vector<int>(graph_info.size(), 0).swap(status);
}

void addNodeToQ(int u, int v, two_hop_case_info &info, int cost, vector<int> &res)
{
    std::pair<std::vector<std::pair<int, int>>, std::vector<int>> dis2hub = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2_find_all_hub(info.L, u, v);

    // // 查询u-v的最短距离 shortest_cost
    // if (dis2hub.first == nullPair)
    // {
    //     return dis2hub.second;
    // }
    if (cost != (dis2hub.first[0].first + dis2hub.first[0].second))
    {
        return;
    }
    int shortest_path_dis = dis2hub.first[0].first + dis2hub.first[0].second;
    // hub一定是结果
    for (const int &hub : dis2hub.second)
    {
        status[hub] = 1;
        res.push_back(hub);
    }
    if (dis2hub.second.size() == 1 && dis2hub.second[0] == u)
    {

        Q_SPG_handles[v] = Q_SPG.push(node_for_SPG_diffuse(v, 0, {shortest_path_dis}, dis2hub.second, u, 0));
    }
    else if (dis2hub.second.size() == 1 && dis2hub.second[0] == v)
    {
        Q_SPG_handles[u] = Q_SPG.push(node_for_SPG_diffuse(u, 0, {shortest_path_dis}, dis2hub.second, v, 0));
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
        Q_SPG_handles[u] = Q_SPG.push(node_for_SPG_diffuse(u, 0, costU, dis2hub.second, v, 1));
        Q_SPG_handles[v] = Q_SPG.push(node_for_SPG_diffuse(v, 0, costV, dis2hub.second, u, 1));
    }
}

template <typename WEIGHT_TYPE>
vector<int> process(int u, int v, int t_s, int t_e, two_hop_case_info &info, graph_v_of_v_with_time_span_non_hop_constrained<WEIGHT_TYPE> &graph_info)
{
    preInit(graph_info);
    std::vector<int> res;
    int costAll = info.query(u, v, t_s, t_e);
    addNodeToQ(u, v, info, costAll, res);
    // SPG的结果集 res
    // 1.2 队列不为空 出队X并遍历邻居(TODO: 此处也可以并行 可以使用乐观锁)
    while (!Q_SPG.empty())
    {
        // (TODO: 并行)
        node_for_SPG_diffuse node = Q_SPG.top();
        Q_SPG.pop();
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
                if (status[edge.vertex] == 1)
                {
                    continue;
                }
                if (max(edge.startTimeLabel, t_s) <= min(edge.endTimeLabel, t_e))
                {
                    if (node.mode == 0)
                    {
                        // mode-0: 二分查找 hub
                        bool isValid = false;
                        vector<int> hubElse;
                        vector<int> hub;
                        node_for_SPG_diffuse tmp(edge.vertex, node.disx + edge.weight, hubElse, hub, node.target, 1);
                        for (int i = 0; i < node.hub.size() && edge.vertex > node.hub[i]; i++)
                        {
                            bool tmp_valid = search_sorted_two_hop_label_specify_time_span_cost(info.L[edge.vertex], node.hub[i], node.hubElse[i] - (node.disx + edge.weight), t_s, t_e);
                            if (tmp_valid)
                            {
                                isValid |= tmp_valid;
                                hubElse.push_back(node.hubElse[i] - (node.disx + edge.weight));
                                hub.push_back(node.hub[i]);
                            }
                        }
                        if (isValid)
                        {
                            Q_SPG_handles[edge.vertex] = Q_SPG.push(tmp);
                        }
                        status[edge.vertex] = 1;
                    }
                    else if (node.mode == 1)
                    {
                        // mode-1: 二分查找 hub
                        bool isValid = false;
                        vector<int> hubElse;
                        vector<int> hub;
                        node_for_SPG_diffuse tmp(edge.vertex, node.disx + edge.weight, hubElse, hub, node.target, 1);
                        for (int i = 0; i < node.hub.size() && edge.vertex > node.hub[i]; i++)
                        {
                            bool tmp_valid = search_sorted_two_hop_label_specify_time_span_cost(info.L[edge.vertex], node.hub[i], node.hubElse[i] - (node.disx + edge.weight), t_s, t_e);
                            if (tmp_valid)
                            {
                                isValid |= tmp_valid;
                                hubElse.push_back(node.hubElse[i] - (node.disx + edge.weight));
                                hub.push_back(node.hub[i]);
                            }
                            else
                            {
                                addNodeToQ(edge.vertex, node.target, info, costAll - (node.disx + edge.weight), res);
                            }
                        }
                        if (isValid)
                        {
                            Q_SPG_handles[edge.vertex] = Q_SPG.push(tmp);
                        }
                        status[edge.vertex] = 1;
                    }
                }
            }
        }
    }
    return res;
}
