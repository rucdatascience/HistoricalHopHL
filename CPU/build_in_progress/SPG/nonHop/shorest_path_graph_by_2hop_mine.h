#include "CPU\build_in_progress\HL\HL4GST\nonHOP_maintain\nonHOP_maintain_two_hop_labels.h"
#include "CPU\graph_v_of_v\graph_v_of_v.h"
#include "Historical\graph_v_of_v\graph_v_of_v_with_time_span_non_hop_constrained.h"
using namespace nonHop;

// define the node in the queue
struct node_for_SPG_diffuse
{
    // 当前的节点
    int index;
    // 剩余距离
    weightTYPE disx;
    // hub
    int hub;
    // 最短路径的另一侧
    int target;
    // mode-0 直接查找 必须包含 不包含可以试做不合法
    // mode-1 也是直接查找，可以不包含 不包含要他的真正目的
    int mode;
    node_for_SPG_diffuse() {}
    node_for_SPG_diffuse(int _u, weightTYPE _dis, int _hub, int _target, int _mode) : index(_u), disx(_dis), hub(_hub), target(_target), mode(_mode) {}
};

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
    std::pair<int, vector<int>> dis2hub = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2_find_all_hub(info.L, u, v);
    // 查询u-v的最短距离 shorest_cost
    if (dis2hub.first == std::numeric_limits<weightTYPE>::max())
    {
        return dis2hub.second;
    }
    // 如果查询的是一个点 直接返回
    if (dis2hub.first == 0)
    {
        return dis2hub.second;
    }
    // SPG的结果集 res
    std::vector<int> res;
    // hub一定是结果
    for (const int &hub : dis2hub.second)
    {
        res.push_back(hub);
    }
    boost::heap::fibonacci_heap<node_for_SPG_diffuse> Q;
    vector<handle_t_for_SPG_diffuse> Q_handles(graph_info.size());
    // 1. 如果hub是不是或者v mode是false 反之mode是true
    if (dis2hub.second != u && dis2hub.second != v)
    {
        // commonhub is a high rank vertex
        Q_handles[u] = Q.push(node_for_SPG_diffuse(u, dis2hub.first, dis2hub.second, v, 1));
        Q_handles[v] = Q.push(node_for_SPG_diffuse(v, dis2hub.first, dis2hub.second, u, 1));
    }
    else if (dis2hub.second == u)
    {
         
        res.push_back(u);
        res.push_back(v);
    }
    else
    {
        Q_handles[u] = Q.push(node_for_SPG_diffuse(u, dis2hub.first, dis2hub.second, v, 0));
    }
    // 1.2 队列不为空 出队X并遍历邻居(TODO: 此处也可以并行 可以使用乐观锁)
    while (!Q.empty())
    {
        // (TODO: 并行)
        node_for_SPG_diffuse node = Q.top();
        Q.pop();
        int x = node.index;
        for (const auto &edges : graph_info[x])
        {
            int v_id = edges.first;
            // 遍历邻居
            // 1. mode-0 直接二分找hub 如果符合则入堆 答案加入 否则直接剪枝
            // 2. mode-1 先二分查找hub
            // 2.1 如果二分hub找到了 并且符合 则入堆 答案加入
            // 2.2 如果没找到或者不符合, 则判断target 是否优先级比自己高。
            //     如果高则在当前的中二分target, 如果找到了target 并且符合
            //       则入堆，并且hub是target mode是0
            //     否则则在直接剪枝
            for (const EdgeInfo<WEIGHT_TYPE> &edge : edges.second)
            {
                if (max(edge.startTimeLabel, t_s) < min(edge.endTimeLabel, t_e))
                {
                    if (node.mode == 0)
                    {
                        // mode-0: 二分查找 hub
                        int cost = search_sorted_two_hop_label(info[edge.vertex], node.hub);
                        if (cost == node.disx - edge.weight)
                        {
                            addToHeap(hub);   // 入堆
                            addToAnswer(hub); // 答案加入
                        }
                        else
                        {
                            continue; // 剪枝
                        }
                    }
                    else if (mode == 1)
                    {
                        // mode-1: 二分查找 hub
                        auto hub = binarySearchHub(info);
                        if (hub && isValidHub(hub, t_s, t_e))
                        {
                            addToHeap(hub);   // 入堆
                            addToAnswer(hub); // 答案加入
                        }
                        else
                        {
                            // 没找到 hub 或不符合，判断 target 优先级
                            if (isHigherPriority(info.target, info.currentNode))
                            {
                                auto target = binarySearchTarget(info);
                                if (target && isValidTarget(target, t_s, t_e))
                                {
                                    addToHeap(target);            // 入堆
                                    setHubToTarget(info, target); // hub 是 target
                                    setMode(info, 0);             // mode 设为 0
                                    addToAnswer(target);          // 答案加入
                                }
                                else
                                {
                                    continue; // 剪枝
                                }
                            }
                            else
                            {
                                continue; // 剪枝
                            }
                        }
                    }
                }
            }
        }
    }
}

void diffuse(int u, int target, int cost, graph_v_of_v<weightTYPE> &graph, two_hop_case_info &info, ThreadPool &pool, std::vector<std::future<vector<int>>> &results_dynamic)
{
}

void search_shorest_path_graph_partial(int u, int v, int cost)
{
}
