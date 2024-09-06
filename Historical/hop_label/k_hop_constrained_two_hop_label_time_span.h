using namespace std;
#include <vector>
#include <list>
#include <iostream>

#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h"
#include "CPU/graph_v_of_v/graph_v_of_v.h"
#include "CPU/graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_binary_operations.h"
class compair_node
{
public:
    int vertex;
    int hop;
    int hub;
    int distance;
    compair_node(int vertex, int hop, int distance) : vertex(vertex), hop(hop), distance(distance)
    {
    }
    compair_node(int vertex, int hub, int hop, int distance) : vertex(vertex), hub(hub), hop(hop), distance(distance)
    {
    }
};

struct compare_diffuse
{
    // vertex,hop,distance
    bool operator()(const compair_node *lhs, const compair_node *rhs) const
    {
        return lhs->hop > rhs->hop;
    }
};

class hop_constrained_two_hop_label_time_span
{
public:
    int hub_vertex, hop;
    int distance;
    int start_time_label, end_time_label;
    hop_constrained_two_hop_label_time_span(int hub_vertex, int hop, int distance, int start_time_label)
        : hub_vertex(hub_vertex), hop(hop), distance(distance), start_time_label(start_time_label)
    {
        end_time_label = INT_MAX;
    }
    hop_constrained_two_hop_label_time_span(int hub_vertex, int hop, int distance, int start_time_label, int end_time_label)
        : hub_vertex(hub_vertex), hop(hop), distance(distance), start_time_label(start_time_label), end_time_label(end_time_label)
    {
    }
    bool operator<(hop_constrained_two_hop_label_time_span const &other)
    {
        if (this->hub_vertex != other.hub_vertex)
        {
            return this->hub_vertex < other.hub_vertex;
        }
        else if (this->end_time_label != other.end_time_label)
        {
            return this->end_time_label < other.end_time_label;
        }
        else if (this->distance != other.distance)
        {
            return this->distance < other.distance;
        }
        else
        {
            return this->hop < other.hop;
        }
    }
};

class graph_hop_constrained_two_hop_label_time_span
{

private:
    int current_time = 0;
    // 邻接表 一次性赋值 后续考虑修改的问题
    vector<vector<pair<int, int>>> AJDS;
    void get_affect_label(int u, int v, int weight, vector<compair_node> &CL);

    void diffuse(vector<compair_node> &CL);

public:
    /**
     * i -> v
     * j -> v_to_label_list
     * k -> label_time_vector
     */
    vector<vector<vector<hop_constrained_two_hop_label_time_span>>> L;

    long int query(int u, int v, int hop);

    graph_hop_constrained_two_hop_label_time_span(int v) : L(v)
    {
    }
    void add_time()
    {
        ++current_time;
    }
    inline void initiate_2_hop_label(graph_v_of_v<int> &graph, hop_constrained_case_info &case_info);
    inline void add_new_edge_or_weight_decrease(int change_sourece, int change_target, int weight);
    inline int test_query_method(int source, int target, int start_time, int end_time, int k);
};

inline long int graph_hop_constrained_two_hop_label_time_span::query(int u, int v, int hop)
{
    long int res = std::numeric_limits<long int>().max();
    int i = 0, j = 0;
    int i_size = this->L[u].size();
    int j_size = this->L[v].size();
    while (i < i_size && j < j_size)
    {
        if (L[u][i].back().hub_vertex < L[v][j].back().hub_vertex)
        {
            ++i;
        }
        else if (L[u][i].back().hub_vertex > L[v][j].back().hub_vertex)
        {
            ++j;
        }
        else
        {
            for (const auto &label_u : L[u][i])
            {
                if (label_u.end_time_label != INT_MAX)
                {
                    continue;
                }
                for (const auto &label_v : L[v][j])
                {
                    if (label_v.end_time_label != INT_MAX)
                    {
                        continue;
                    }
                    if (label_u.hop + label_v.hop <= hop)
                    {
                        res = min(res, (long int)label_u.distance + label_v.distance);
                    }
                }
            }
            ++i;
        }
    }
    return res;
}

inline void graph_hop_constrained_two_hop_label_time_span::get_affect_label(int change_sourece, int change_target, int weight, vector<compair_node> &CL)
{
    for (const auto &label_list : L[change_sourece])
    {
        if (label_list.back().hub_vertex < change_target)
        {
            continue;
        }
        for (const auto &label : label_list)
        {
            if (label.end_time_label == INT_MAX)
            {
                if (this->query(label.hub_vertex, change_target, label.hop + 1) > label.distance + weight)
                {
                    CL.push_back({change_target, label.hub_vertex, label.hop + 1, label.distance + weight});
                }
                else
                {
                    for (const auto &label_list_target : L[change_target])
                    {
                        if (label_list_target.back().hub_vertex != label.hub_vertex)
                        {
                            continue;
                        }
                        else
                        {
                            for (const auto &label_target : label_list_target)
                            {
                                if (label_target.end_time_label == INT_MAX && label_target.hop <= label.hop + 1 && label_target.distance > label.distance + weight)
                                {
                                    CL.push_back({change_target, label.hub_vertex, label.hop + 1, label.distance + weight});
                                }
                            }
                        }
                    }
                    // TODO 修改PPR
                }
            }
        }
    }
}

inline void graph_hop_constrained_two_hop_label_time_span::diffuse(vector<compair_node> &CL)
{
    vector<pair<int, int>> H(this->L.size(), {-1, -1});
    boost::heap::fibonacci_heap<compair_node *, boost::heap::compare<compare_diffuse>> queue;
    auto start = CL.begin();
    auto end = CL.end();
    while (start < end)
    {
        int v = (*start).hub;
        auto start_temp = start;
        while ((*start_temp).hub == v)
        {
            H[(*start_temp).vertex] = {(*start_temp).distance, (*start_temp).hop};
            ++start_temp;
            compair_node temp = {(*start_temp).vertex, (*start_temp).hop, (*start_temp).distance};
            queue.push(&temp);
        }
        start = start_temp;
        while (!queue.empty())
        {
            auto item = queue.top();
            queue.pop();
            int x = item->vertex;
            int hopxv = item->hop;
            int disxv = item->distance;
            for (auto &label_list : L[x])
            {
                if (label_list.back().hub_vertex == v)
                {
                    bool is_found = false;
                    for (auto &label : label_list)
                    {
                        // TODO这里实际还是遍历了
                        if (label.end_time_label == INT_MAX && label.hop >= hopxv && label.distance > disxv)
                        {
                            label.end_time_label = current_time - 1;
                            is_found = true;
                        }
                    }
                    if (is_found)
                    {
                        label_list.push_back(hop_constrained_two_hop_label_time_span(x, v, disxv, current_time));
                    }
                }
            }
            // 这里没写k
            // 遍历邻居
            for (std::pair<int, int> edge : this->AJDS[x])
            {
                int xn = edge.first;
                if (v > xn)
                {
                    if (H[xn].first != -1)
                    {
                        H[xn].first = query(xn, v, hopxv + 1);
                    }
                    int dis_x_xn = sorted_vector_binary_operations_search_weight(this->AJDS[xn], x);
                    if (H[xn].first > disxv + this->AJDS[xn][x].second + dis_x_xn)
                    {
                        H[xn].first = disxv + this->AJDS[xn][x].second + dis_x_xn;
                        H[xn].second = hopxv + 1;
                        for (const auto element : queue)
                        {
                            if (element->vertex == xn)
                            {
                                element->hop = hopxv + 1;
                                element->distance = disxv + this->AJDS[xn][x].second + dis_x_xn;
                            }
                        }
                    }
                    else
                    {
                        if (hopxv + 1 < H[xn].second)
                        {
                            for (const auto element : queue)
                            {
                                if (element->vertex == xn)
                                {
                                    element->hop = hopxv + 1;
                                    element->distance = disxv + this->AJDS[xn][x].second + dis_x_xn;
                                }
                            }
                        }
                        // TODO PPR
                    }
                }
            }
        }
    }
}

void modify_value_in_q(boost::heap::fibonacci_heap<boost::heap::fibonacci_heap<tuple<int, int, int>, boost::heap::compare<compare_diffuse>>> &q, int value)
{
}

inline void graph_hop_constrained_two_hop_label_time_span::initiate_2_hop_label(graph_v_of_v<int> &graph, hop_constrained_case_info &case_info)
{
    int index = 0;
    hop_constrained_two_hop_labels_generation(graph, case_info);
    vector<vector<hop_constrained_two_hop_label>> result = case_info.L;
    for (const auto &label_list : result)
    {
        for (const auto &label : label_list)
        {
            if (!this->L[index].empty() && this->L[index].back().back().hub_vertex == label.hub_vertex)
            {
                this->L[index].back().push_back({hop_constrained_two_hop_label_time_span(label.hub_vertex, label.hop, label.distance, current_time)});
            }
            else
            {
                this->L[index].push_back({hop_constrained_two_hop_label_time_span(label.hub_vertex, label.hop, label.distance, current_time)});
            }
        }
        ++index;
    }
    this->AJDS = graph.ADJs;
    add_time();
}

inline void graph_hop_constrained_two_hop_label_time_span::add_new_edge_or_weight_decrease(int change_sourece, int change_target, int weight)
{
    vector<compair_node> CL;
    get_affect_label(change_sourece, change_target, weight, CL);
    get_affect_label(change_target, change_sourece, weight, CL);
    sort(CL.begin(), CL.end(), [](compair_node x, compair_node y)
         { return x.hub < y.hub; });
    diffuse(CL);
    add_time();
}

int graph_hop_constrained_two_hop_label_time_span::test_query_method(int source, int target, int start_time, int end_time, int k)
{
    int res = INT_MAX;
    auto &M_source = this->L[source];
    auto &M_target = this->L[target];
    vector<vector<hop_constrained_two_hop_label_time_span>>::iterator it_source = M_source.begin();
    vector<vector<hop_constrained_two_hop_label_time_span>>::iterator it_target = M_target.begin();
    while (it_source != M_source.end() && it_target != M_target.end())
    {
        if ((*it_source).back().hub_vertex > (*it_target).back().hub_vertex)
        {
            ++it_target;
        }
        else if ((*it_source).back().hub_vertex < (*it_target).back().hub_vertex)
        {
            ++it_source;
        }
        else
        {
            int i = 0, j = 0;

            vector<hop_constrained_two_hop_label_time_span> &source_L = (*it_source);
            vector<hop_constrained_two_hop_label_time_span> &target_L = (*it_target);
            int max_i = source_L.size(), max_j = target_L.size();
            while (i < max_i && j < max_j)
            {
                int i_num = 0, j_num = 0;
                while (!(source_L[i].start_time_label <= end_time && source_L[i].end_time_label >= start_time))
                {
                    ++i;
                }
                while (i + i_num < max_i && source_L[i + i_num].start_time_label <= end_time && source_L[i + i_num].end_time_label >= start_time)
                {
                    ++i_num;
                }
                while (!(target_L[j].start_time_label <= end_time && target_L[j].end_time_label >= start_time))
                {
                    ++j;
                }
                while (j + j_num < max_j && target_L[j + j_num].start_time_label <= end_time && target_L[j + j_num].end_time_label >= start_time)
                {
                    ++j_num;
                }
                for (int index_i = 0; index_i < i_num; index_i++)
                {
                    for (int index_j = 0; index_j < j_num; index_j++)
                    {
                        if (max(source_L[i + index_i].start_time_label, target_L[j + index_j].start_time_label) <= min(source_L[i + index_i].end_time_label, target_L[j + index_j].end_time_label))
                        {
                            if (source_L[i + index_i].hop + target_L[j + index_j].hop <= k)
                            {
                                res = min(res, source_L[i + index_i].distance + target_L[j + index_j].distance);
                            }
                        }
                    }
                }
                i += i_num;
                j += j_num;
            }
            ++it_target;
        }
    }

    return res;
}