using namespace std;
#include <vector>
#include <list>
#include <iostream>

#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h"
#include "CPU/graph_v_of_v/graph_v_of_v.h"
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

public:
    /**
     * i -> v
     * j -> v_to_label_list
     * k -> label_time_vector
     */
    vector<list<vector<hop_constrained_two_hop_label_time_span>>> L;
    graph_hop_constrained_two_hop_label_time_span(int v) : L(v)
    {
    }
    inline void add_new_hop_label(vector<vector<hop_constrained_two_hop_label>> labels, int time);
    inline void sort_L();
    inline int test_query_method(int source, int target, int start_time, int end_time, int k);
};

void graph_hop_constrained_two_hop_label_time_span::add_new_hop_label(vector<vector<hop_constrained_two_hop_label>> labels, int time)
{
    for (int i = 0; i < labels.size(); i++)
    {
        list<vector<hop_constrained_two_hop_label_time_span>>::iterator it = this->L[i].begin();
        for (int j = 0; j < labels[i].size();)
        {
            while (it != L[i].end() && (*it).back().hub_vertex < labels[i][j].hub_vertex)
            {
                ++it;
            }
            if (it == L[i].end() || (*it).back().hub_vertex > labels[i][j].hub_vertex)
            {
                it = L[i].insert(it, {hop_constrained_two_hop_label_time_span(labels[i][j].hub_vertex, labels[i][j].hop, labels[i][j].distance, time, time)});
                ++j;
            }
            else if ((*it).back().hub_vertex == labels[i][j].hub_vertex)
            {
                int found = false;
                for (auto &label : (*it))
                {
                    if (label.distance == labels[i][j].distance && label.hop == labels[i][j].hop && label.end_time_label == time - 1)
                    {
                        label.end_time_label = time;
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    (*it).push_back(hop_constrained_two_hop_label_time_span(labels[i][j].hub_vertex, labels[i][j].hop, labels[i][j].distance, time, time));
                }
                ++j;
            }
            else
            {
                ++it;
            }
        }
    }
    sort_L();
}

void graph_hop_constrained_two_hop_label_time_span::sort_L()
{
    for (auto &next_to_list : L)
    {
        for (auto &time_labels : next_to_list)
        {
            sort(time_labels.begin(), time_labels.end());
        }
    }
}

int graph_hop_constrained_two_hop_label_time_span::test_query_method(int source, int target, int start_time, int end_time, int k)
{
    int res = INT_MAX;
    auto &M_source = this->L[source];
    auto &M_target = this->L[target];
    list<vector<hop_constrained_two_hop_label_time_span>>::iterator it_source = M_source.begin();
    list<vector<hop_constrained_two_hop_label_time_span>>::iterator it_target = M_target.begin();
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