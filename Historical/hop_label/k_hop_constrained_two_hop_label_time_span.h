using namespace std;
#include <vector>
#include <list>
#include <iostream>

#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels.h"
#include "CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h"
#include "CPU/graph_v_of_v/graph_v_of_v.h"
#include "CPU/graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_binary_operations.h"
#include <map>

vector<vector<pair<int, int>>> dist_hop_599_v2_test, dist_hop_599_v3_test;
queue<int> Qid_599_v2_test, Qid_599_v3_test;
vector<vector<vector<int>>> Q_value_test;
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

void WeightDecreaseMaintenance_improv_step1_batch(std::map<pair<int, int>, weightTYPE> &v_map, vector<vector<hop_constrained_two_hop_label>> *L, PPR_type *PPR, std::vector<hop_constrained_affected_label> *CL,
                                                  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{
    for (auto it : v_map)
    {
        results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, PPR, CL]
                                                          {
            int v1 = it.first.first, v2 = it.first.second;
            weightTYPE w_new = it.second;
            for (int sl = 0; sl < 2; sl++)
            {
                if (sl == 1)
                {
                    swap(v1, v2);
                }
                for (auto it : (*L)[v1])
                {
                    if (it.hub_vertex <= v2 && (long long int)it.distance + w_new < TwoM_value)
                    {
                            auto query_result = hop_constrained_extract_distance_and_hub(*L, it.hub_vertex, v2, it.hop + 1); // query_result is {distance, common hub}
                            if ((long long int)query_result.first > (long long int)it.distance + w_new)
                            {
                                mtx_599_1.lock();
                                CL->push_back(hop_constrained_affected_label{v2, it.hub_vertex, it.hop + 1, it.distance + w_new});
                                mtx_599_1.unlock();
                            }
                            else
                            {
                                auto search_result = search_sorted_hop_constrained_two_hop_label((*L)[v2], it.hub_vertex, it.hop + 1);
                                if (search_result < MAX_VALUE && search_result > it.distance + w_new)
                                {
                                    mtx_599_1.lock();
                                    CL->push_back(hop_constrained_affected_label{v2, it.hub_vertex, it.hop + 1, it.distance + w_new});
                                    mtx_599_1.unlock();
                                }
                                if (query_result.second != -1 && query_result.second != it.hub_vertex)
                                {
                                    mtx_5992[v2].lock();
                                    PPR_insert(*PPR, v2, query_result.second, it.hub_vertex);
                                    mtx_5992[v2].unlock();
                                }
                                if (query_result.second != -1 && query_result.second != v2)
                                {
                                    mtx_5992[it.hub_vertex].lock();
                                    PPR_insert(*PPR, it.hub_vertex, query_result.second, v2);
                                    mtx_5992[it.hub_vertex].unlock();
                                }
                            }
					}
                    }
                }
                return 1; }));
    }

    for (auto &&result : results_dynamic)
    {
        result.get();
    }
    std::vector<std::future<int>>().swap(results_dynamic);
}

void DIFFUSE_batch(graph_v_of_v<int> &instance_graph, vector<vector<hop_constrained_two_hop_label>> *L, PPR_type *PPR, std::vector<hop_constrained_affected_label> &CL,
                   ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int upper_k)
{
    std::map<hop_constrained_pair_label, weightTYPE> CL_edge_map;
    for (auto &it : CL)
    {
        if (CL_edge_map.count({it.first, it.second, it.hop}) == 0)
        {
            CL_edge_map[{it.first, it.second, it.hop}] = it.dis;
        }
        else if (CL_edge_map[{it.first, it.second, it.hop}] > it.dis)
        {
            CL_edge_map[{it.first, it.second, it.hop}] = it.dis;
        }
    }

    // extract each unique hub v and its (u,hop,dis) list
    std::map<int, std::vector<hop_constrained_label_v2>> CL_map;
    for (auto &it : CL_edge_map)
    {
        int u = it.first.first;
        int v = it.first.second;
        int hop = it.first.hop;
        weightTYPE dis = it.second;
        if (CL_map.count(v) == 0)
        {
            std::vector<hop_constrained_label_v2> vec_with_hub_v;
            hop_constrained_label_v2 tmp(u, hop, dis);
            vec_with_hub_v.emplace_back(tmp);
            CL_map[v] = vec_with_hub_v;
        }
        else
        {
            std::vector<hop_constrained_label_v2> vec_with_hub_v = CL_map[v];
            hop_constrained_label_v2 tmp(u, hop, dis);
            vec_with_hub_v.emplace_back(tmp);
            CL_map[v] = vec_with_hub_v;
        }
    }

    for (auto &it : CL_map)
    {
        results_dynamic.emplace_back(pool_dynamic.enqueue([it, L, &instance_graph, PPR, upper_k]
                                                          {
			mtx_599_1.lock();
			int current_tid = Qid_599_v2_test.front();
			Qid_599_v2.pop();
			mtx_599_1.unlock();

            int v = it.first;
            std::vector<hop_constrained_label_v2> vec_with_hub_v = it.second;

			mtx_599[v].lock_shared();
			auto Lv = (*L)[v]; // to avoid interlocking
			mtx_599[v].unlock_shared();

			vector<int> dist_hop_changes;
            // rede
			auto & dist_hop = dist_hop_599_v2_test[current_tid];
			boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE> pq;
			map<pair<int, int>, pair<hop_constrained_handle_t_for_DIFFUSE, int>> Q_handle;
            vector<int> hubs;
			hubs.resize(instance_graph.size(), -1);
            auto& Q_VALUE = Q_value_test[current_tid];
            
            for(auto &it:vec_with_hub_v){
				int u = it.hub_vertex;
                int h_v = it.hop;
                weightTYPE du = it.distance;

				dist_hop[u] = {du, h_v}; //  {dis, hop}
                dist_hop_changes.push_back(u);
                hop_constrained_node_for_DIFFUSE tmp;
                tmp.index = u;
                tmp.hop = h_v;
                tmp.disx = du;
                Q_handle[{u, h_v}] = {pq.push({tmp}), du}; //{node_for_DIFFUSE_v2,dis}
                if(h_v <= upper_k)
                    Q_VALUE[u][h_v] = du;
                
            }
            

			while (!pq.empty())
			{

				int x = pq.top().index;
				int xhv = pq.top().hop;
				weightTYPE dx = pq.top().disx;
				pq.pop();
                if(xhv <= upper_k)
                    Q_VALUE[x][xhv] = MAX_VALUE;

				mtx_599[x].lock();
				weightTYPE d_old = search_sorted_hop_constrained_two_hop_label((*L)[x], v, xhv);
				if (dx > 0 && dx < d_old)
				{
                    insert_sorted_hop_constrained_two_hop_label((*L)[x], v, xhv, dx);
				}
				mtx_599[x].unlock();
                
                if(xhv + 1 > upper_k)
					continue;

				for (auto nei : instance_graph[x])
				{
					
                    if (dx + nei.second >= TwoM_value)
					    continue;
				
					int xnei = nei.first;
                    int hop_nei = xhv + 1;
                    long long int d_new = dx + (long long int)nei.second;
					hop_constrained_node_for_DIFFUSE node = {xnei, hop_nei, (weightTYPE)d_new};

					if (v < xnei)
					{

						if (dist_hop[xnei].first == -1)
						{
							
							Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                            Q_VALUE[xnei][hop_nei] = d_new;
                            dist_hop[xnei].first = d_new;
							dist_hop[xnei].second = hop_nei;
							dist_hop_changes.push_back(xnei);

                            mtx_599[xnei].lock_shared();
							std::pair<int, int> tmp = hop_constrained_extract_distance_and_hub_2((*L)[xnei], Lv, xhv + 1); 
							mtx_599[xnei].unlock_shared();
							hubs[xnei] = tmp.second;
						}
                        
                        if (d_new < dist_hop[xnei].first)
						{
                            
                            //if (Q_handle.find({xnei, hop_nei}) != Q_handle.end())
							if(Q_VALUE[xnei][hop_nei] < MAX_VALUE)
                            {
								if (Q_handle[{xnei, hop_nei}].second > d_new)
								{
									pq.update(Q_handle[{xnei, hop_nei}].first, node);
									Q_handle[{xnei, hop_nei}].second = d_new;
								}
							}
							else
							{
								Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
							}
							dist_hop[xnei].first = d_new;
							dist_hop[xnei].second = hop_nei;
							dist_hop_changes.push_back(xnei);
                            hubs[xnei] = v;
                            Q_VALUE[xnei][hop_nei] = d_new;
                        }
						else if (hop_nei < dist_hop[xnei].second)
						{
                            
							//if (Q_handle.find({xnei, hop_nei}) != Q_handle.end())
							if(Q_VALUE[xnei][hop_nei] < MAX_VALUE)
                            {
								if (Q_handle[{xnei, hop_nei}].second > d_new)
								{
									pq.update(Q_handle[{xnei, hop_nei}].first, node);
									Q_handle[{xnei, hop_nei}].second = d_new;
								}
							}
							else
							{
								Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
							}
                            Q_VALUE[xnei][hop_nei] = d_new;
						}

						
                        if(dist_hop[xnei].first < d_new)
                        {
                            if (hubs[xnei] != -1 && hubs[xnei] != v)
							{
								mtx_5992[xnei].lock();
								PPR_insert(*PPR, xnei, hubs[xnei], v);
								mtx_5992[xnei].unlock();
							}
							if (hubs[xnei] != -1 && hubs[xnei] != xnei)
							{
								mtx_5992[v].lock();
								PPR_insert(*PPR, v, hubs[xnei], xnei);
								mtx_5992[v].unlock();
							}
                        }
							

					}
				}
			}


			for(int i : dist_hop_changes)
			{
				dist_hop[i] = {-1, 0};
			}

			mtx_599_1.lock();
			Qid_599_v2.push(current_tid);
			mtx_599_1.unlock();
			return 1; }));
    }

    for (auto &&result : results_dynamic)
    {
        result.get();
    }
    std::vector<std::future<int>>().swap(results_dynamic);
}

void HOP_WeightDecreaseMaintenance_improv_batch(graph_v_of_v<int> &instance_graph, hop_constrained_case_info &mm,
                                                std::vector<pair<int, int>> &v, std::vector<int> &w_new, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

    global_query_times = 0;
    label_operation_times = 0;

    std::map<pair<int, int>, int> w_new_map;
    int batch_size = v.size();
    for (int i = 0; i < batch_size; i++)
    {
        if (v[i].first > v[i].second)
        {
            swap(v[i].first, v[i].second);
        }
        if (w_new_map.count(v[i]) == 0)
        {
            w_new_map[v[i]] = w_new[i];
        }
        else if (w_new_map[v[i]] > w_new[i])
        {
            w_new_map[v[i]] = w_new[i];
        }
    }
    std::vector<hop_constrained_affected_label> CL;

    WeightDecreaseMaintenance_improv_step1_batch(w_new_map, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic);
    DIFFUSE_batch(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic, mm.upper_k);
}

class graph_hop_constrained_two_hop_label_time_span
{

private:
    int current_time = 0;
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

    long int query(int u, int v, int &hop);

    graph_hop_constrained_two_hop_label_time_span(int v) : L(v)
    {
    }
    void add_time()
    {
        ++current_time;
    }
    inline void initiate_2_hop_label(graph_v_of_v<int> &graph, hop_constrained_case_info &case_info);
    inline void process(vector<graph_v_of_v<int>> graphs, hop_constrained_case_info &case_info);
    inline void add_new_edge_or_weight_decrease(int change_sourece, int change_target, int weight);
    inline int test_query_method(int source, int target, int start_time, int end_time, int k);
    inline void print_L();
    inline void print_L_by_index(int);
};

inline long int graph_hop_constrained_two_hop_label_time_span::query(int u, int v, int &hop)
{
    if (u == v)
    {
        return 0;
    }
    long int res = std::numeric_limits<long int>().max();
    int i = 0, j = 0;
    int i_size = this->L[u].size();
    int j_size = this->L[v].size();
    int hop_upper = hop;
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
                    if (label_u.hop + label_v.hop <= hop_upper)
                    {
                        res = min(res, (long int)(label_u.distance + label_v.distance));
                        hop = label_u.hop + label_v.hop;
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
        if (label_list.back().hub_vertex > change_target)
        {
            continue;
        }
        for (const auto &label : label_list)
        {
            if (label.end_time_label == INT_MAX)
            {
                int hop_upper = label.hop + 1;
                if (this->query(label.hub_vertex, change_target, hop_upper) > label.distance + weight)
                {
                    CL.push_back({change_target, label.hub_vertex, label.hop + 1, label.distance + weight});
                    break;
                }
                else
                {
                    bool is_found = false;
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
                                    is_found = true;
                                    break;
                                }
                            }
                            if (is_found)
                            {
                                break;
                            }
                        }
                    }
                    if (is_found)
                    {
                        break;
                    }
                    // TODO 修改PPR
                }
            }
        }
    }
}

inline void graph_hop_constrained_two_hop_label_time_span::diffuse(vector<compair_node> &CL)
{
    boost::heap::fibonacci_heap<compair_node *, boost::heap::compare<compare_diffuse>> queue;
    auto start = CL.begin();
    auto end = CL.end();
    while (start < end)
    {
        vector<pair<int, int>> H(this->L.size(), {-1, -1});
        queue.clear();
        int v = (*start).hub;
        auto start_temp = start;
        while ((*start_temp).hub == v)
        {
            H[(*start_temp).vertex] = {(*start_temp).distance, (*start_temp).hop};
            compair_node temp = {(*start_temp).vertex, (*start_temp).hop, (*start_temp).distance};
            queue.push(&temp);
            ++start_temp;
        }
        start = start_temp;
        while (!queue.empty())
        {
            auto item = queue.top();
            queue.pop();
            int x = item->vertex;
            int hopxv = item->hop;
            int disxv = item->distance;
            int index = 0;
            for (auto &label_list : L[x])
            {
                if (label_list.back().hub_vertex < v)
                {
                    ++index;
                    continue;
                }
                else if (label_list.back().hub_vertex == v)
                {
                    /**
                     * 感觉有问题 这里的假设
                     */
                    /**
                     * There are four situations regarding the old value compared to the new value:
                     * h large, d large: The old value is redundant; it can be retained but will never be considered as a result.
                     * h large, d small: The old value should be retained as it contains the correct answer.
                     * h small, d large: The old value should be retained as it contains the correct answer.
                     * h small, d small: The new value is redundant. It can be retained but will not be the part of result and might even be deleted.
                     */
                    // bool is_found = false;
                    // for (auto &label : label_list)
                    // {
                    //     // TODO这里实际还是遍历了
                    //     if (label.end_time_label == INT_MAX && label.hop >= hopxv && label.distance > disxv)
                    //     {
                    //         label.end_time_label = current_time - 1;
                    //         is_found = true;
                    //     }
                    // }
                    // if (is_found)
                    // {
                    //     label_list.push_back(hop_constrained_two_hop_label_time_span(x, v, disxv, current_time));
                    // }
                    label_list.push_back(hop_constrained_two_hop_label_time_span(v, hopxv, disxv, current_time));
                    break;
                }
                else
                {
                    vector<hop_constrained_two_hop_label_time_span> temp;
                    temp.push_back(hop_constrained_two_hop_label_time_span(v, hopxv, disxv, current_time));
                    L[x].insert(L[x].begin() + index, temp);
                    break;
                }
            }
            // 这里没写k
            // 遍历邻居
            for (std::pair<int, int> edge : this->AJDS[x])
            {
                int xn = edge.first;
                if (v < xn)
                {
                    if (H[xn].first == -1)
                    {
                        int hop_upper = hopxv + 1;
                        H[xn].first = query(xn, v, hop_upper);
                        H[xn].second = hop_upper;
                    }
                    int dis_x_xn = sorted_vector_binary_operations_search_weight(this->AJDS[xn], x);
                    if (H[xn].first > disxv + dis_x_xn)
                    {
                        H[xn].first = disxv + dis_x_xn;
                        H[xn].second = hopxv + 1;
                        bool is_found = false;
                        for (auto it = queue.begin(); it != queue.end(); ++it)
                        {
                            auto item = (*it);
                            if (item->vertex == xn)
                            {
                                item->hop = hopxv + 1;
                                item->distance = disxv + dis_x_xn;
                                is_found = true;
                                break;
                            }
                        }
                        if (!is_found)
                        {
                            compair_node temp = {xn, H[xn].second, H[xn].first};
                            queue.push(&temp);
                        }
                    }
                    else
                    {
                        bool is_found = false;
                        if (hopxv + 1 < H[xn].second)
                        {
                            for (const auto element : queue)
                            {
                                if (element->vertex == xn)
                                {
                                    is_found = true;
                                    element->hop = hopxv + 1;
                                    element->distance = disxv + dis_x_xn;
                                }
                            }
                            if (!is_found)
                            {
                                compair_node temp = {xn, hopxv + 1, disxv + dis_x_xn};
                                queue.push(&temp);
                            }
                        }
                        // TODO PPR
                    }
                }
            }
        }
    }
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
                this->L[index].back().push_back(hop_constrained_two_hop_label_time_span(label.hub_vertex, label.hop, label.distance, current_time));
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

inline void graph_hop_constrained_two_hop_label_time_span::process(vector<graph_v_of_v<int>> graphs, hop_constrained_case_info &info)
{
    int batch_size = 10;
    this->initiate_2_hop_label(graphs[0], info);
    dist_hop_599_v2_test.resize(info.thread_num);
    dist_hop_599_v3_test.resize(info.thread_num);
    queue<int>().swap(Qid_599_v2_test);
    queue<int>().swap(Qid_599_v3_test);
    Q_value_test.resize(info.thread_num);
    for (int index = 0; index < info.thread_num; index++)
    {
        dist_hop_599_v2_test[index].resize(this->AJDS.size(), {-1, 0});
        dist_hop_599_v3_test[index].resize(this->AJDS.size(), {-1, 0});
        Qid_599_v2_test.push(index);
        Qid_599_v3_test.push(index);
        Q_value_test[index].resize(this->AJDS.size(), vector<int>(10 + 1, MAX_VALUE));
    }
    ThreadPool pool_dynamic(info.thread_num);
    std::vector<std::future<int>> results_dynamic;
    vector<pair<int, int>> vertex_pair_decrease;
    vector<int> weight_decrease;
    vector<pair<int, int>> vertex_pair_increase;
    vector<int> weight_increase;
    for (int index = 1; index < graphs.size(); index++)
    {
        vector<pair<int, int>>().swap(vertex_pair_decrease);
        vector<int>().swap(weight_decrease);
        graph_v_of_v<int> &cur_graph = graphs[index];
        for (int i = 0; i < cur_graph.ADJs.size(); i++)
        {
            for (int j = 0; j < cur_graph.ADJs[i].size(); j++)
            {
                if (i < cur_graph.ADJs[i][j].first && cur_graph.ADJs[i][j].second < graphs[index - 1].ADJs[i][j].second)
                {
                    std::cout << "DePLL:vertex " << i << "-> vertex " << cur_graph.ADJs[i][j].first << " from " << graphs[index - 1].ADJs[i][j].second << " to " << cur_graph.ADJs[i][j].second << std::endl;
                    // store the changing data
                    vertex_pair_decrease.push_back({i, cur_graph.ADJs[i][j].first});
                    weight_decrease.push_back(cur_graph.ADJs[i][j].second);
                    if (weight_decrease.size() >= batch_size)
                    {
                        HOP_WeightDecreaseMaintenance_improv_batch(cur_graph, info, vertex_pair_decrease, weight_decrease, pool_dynamic, results_dynamic);
                    }
                }
                else if (i < cur_graph.ADJs[i][j].first && cur_graph.ADJs[i][j].second > graphs[index - 1].ADJs[i][j].second)
                {
                    std::cout << "InPLL:vertex " << i << "-> vertex " << cur_graph.ADJs[i][j].first << " from " << graphs[index - 1].ADJs[i][j].second << " to " << cur_graph.ADJs[i][j].second << std::endl;
                    // store the changing data
                    vertex_pair_increase.push_back({i, cur_graph.ADJs[i][j].first});
                    weight_increase.push_back(cur_graph.ADJs[i][j].second);
                    if (weight_decrease.size() >= batch_size)
                    {
                        // async inpll
                    }
                }
            }
        }
        if (weight_decrease.size() > 0)
        {
            HOP_WeightDecreaseMaintenance_improv_batch(cur_graph, info, vertex_pair_decrease, weight_decrease, pool_dynamic, results_dynamic);
        }
        if (weight_increase.size() > 0)
        {
            // async inpll
        }
    }
}

inline void graph_hop_constrained_two_hop_label_time_span::add_new_edge_or_weight_decrease(int change_sourece, int change_target, int weight)
{
    vector<compair_node> CL;
    sorted_vector_binary_operations_insert(this->AJDS[change_sourece], change_target, weight);
    sorted_vector_binary_operations_insert(this->AJDS[change_target], change_sourece, weight);
    auto time5 = std::chrono::high_resolution_clock::now();
    get_affect_label(change_sourece, change_target, weight, CL);
    auto time6 = std::chrono::high_resolution_clock::now();
    auto time7 = std::chrono::high_resolution_clock::now();
    get_affect_label(change_target, change_sourece, weight, CL);
    auto time8 = std::chrono::high_resolution_clock::now();
    double section3 = std::chrono::duration_cast<std::chrono::nanoseconds>(time6 - time5).count() / 1e9;
    double section4 = std::chrono::duration_cast<std::chrono::nanoseconds>(time8 - time7).count() / 1e9;
    std::cout << "secton3 and section 4 time cost is " << section3 << ":" << section4 << std::endl;
    sort(CL.begin(), CL.end(), [](auto x, auto y)
         { return x.hub < y.hub; });
    auto time9 = std::chrono::high_resolution_clock::now();
    diffuse(CL);
    auto time10 = std::chrono::high_resolution_clock::now();
    double section5 = std::chrono::duration_cast<std::chrono::nanoseconds>(time10 - time9).count() / 1e9;
    std::cout << "secton5 time cost is " << section5 << std::endl;
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

inline void graph_hop_constrained_two_hop_label_time_span::print_L()
{
    /**
     * i -> v
     * j -> v_to_label_list
     * k -> label_time_vector
     */
    vector<vector<vector<hop_constrained_two_hop_label_time_span>>> L;
    int i = 0;
    for (const auto &i_next_edge : this->L)
    {
        std::cout << "vertex" << i << " :" << std::endl;
        for (const auto &i_j_edge : i_next_edge)
        {
            int j = i_j_edge.back().hub_vertex;
            for (const auto &edge_by_time : i_j_edge)
            {
                std::cout << "\t(" << j << "," << edge_by_time.hop << "," << edge_by_time.distance << "," << edge_by_time.start_time_label << "," << edge_by_time.end_time_label << ")";
            }
            std::cout << std::endl;
        }
        ++i;
    }
}

inline void graph_hop_constrained_two_hop_label_time_span::print_L_by_index(int index)
{
    /**
     * i -> v
     * j -> v_to_label_list
     * k -> label_time_vector
     */
    auto i_next_edge = this->L[index];
    std::cout << "vertex" << index << " :" << std::endl;
    for (const auto &i_j_edge : i_next_edge)
    {
        int j = i_j_edge.back().hub_vertex;
        for (const auto &edge_by_time : i_j_edge)
        {
            std::cout << "\t(" << j << "," << edge_by_time.hop << "," << edge_by_time.distance << "," << edge_by_time.start_time_label << "," << edge_by_time.end_time_label << ")";
        }
        std::cout << std::endl;
    }
}
