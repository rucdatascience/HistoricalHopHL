#pragma once
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_PPR.h>

#define weightTYPE int

/* label format */
class hop_constrained_two_hop_label
{
public:
    int hub_vertex, hop;
    weightTYPE distance;
    int t_s, t_e;
    // hop_constrained_two_hop_label() {}
    // hop_constrained_two_hop_label(int _vertex, int _hop, int _dis)
    // {
    //     hub_vertex = _vertex;
    //     hop = _hop;
    //     distance = _dis;
    // }
    hop_constrained_two_hop_label()
        : t_s(0), t_e(std::numeric_limits<int>::max()) {}
};

class hop_constrained_case_info
{
public:
    /*hop bounded*/
    int thread_num = 1;
    int upper_k = 0;
    bool use_2M_prune = false;
    bool use_2023WWW_generation = false;
    bool use_canonical_repair = 1;
    bool use_rank_prune = 1;

    /*running time records*/
    double time_initialization = 0;
    double time_generate_labels = 0;
    double time_sortL = 0;
    double time_canonical_repair = 0;
    double time_total = 0;

    /*running limits*/
    long long int max_bit_size = 1e12;
    double max_run_time_seconds = 1e12;

    /*labels*/
    vector<vector<hop_constrained_two_hop_label>> L;
    PPR_type PPR;

    double label_size_before_canonical_repair, label_size_after_canonical_repair, canonical_repair_remove_label_ratio;

    long long int compute_label_bit_size()
    {
        long long int size = 0;
        for (auto &xx : L)
        {
            size = size + xx.size() * sizeof(hop_constrained_two_hop_label);
        }
        return size;
    }

    /*clear labels*/
    void clear_labels()
    {
        vector<vector<hop_constrained_two_hop_label>>().swap(L);
        PPR_type().swap(PPR);
    }

    void print_L()
    {
        int index = 0;
        cout << "print_L: (hub_vertex, hop, distance)" << endl;
        for (auto &xx : L)
        {
            cout << "vertex " << index++ << ": ";
            for (auto &yy : xx)
            {
                cout << "(" << yy.hub_vertex << "," << yy.hop << "," << yy.distance << "," << yy.t_s << "," << yy.t_e << ")";
            }
            cout << endl;
        }
    }

    void print_PPR()
    {
        cout << "print_PPR:" << endl;
        for (int i = 0; i < PPR.size(); i++)
        {
            for (int j = 0; j < PPR[i].size(); j++)
            {
                cout << "PPR(" << i << "," << PPR[i][j].first << "): ";
                for (int k = 0; k < PPR[i][j].second.size(); k++)
                {
                    cout << PPR[i][j].second[k] << " ";
                }
                cout << endl;
            }
        }
    }

    void print_L_vk(int v_k)
    {
        for (auto it = L[v_k].begin(); it != L[v_k].end(); it++)
        {
            cout << "<" << it->hub_vertex << "," << it->distance << "," << it->hop << ">";
        }
        cout << endl;
    }

    /*record_all_details*/
    void record_all_details(string save_name)
    {
        ofstream outputFile;
        outputFile.precision(6);
        outputFile.setf(ios::fixed);
        outputFile.setf(ios::showpoint);
        outputFile.open(save_name + ".txt");

        outputFile << "hop_constrained_case_info:" << endl;
        outputFile << "thread_num=" << thread_num << endl;
        outputFile << "upper_k=" << upper_k << endl;
        outputFile << "use_2M_prune=" << use_2M_prune << endl;
        outputFile << "use_2023WWW_generation=" << use_2023WWW_generation << endl;
        outputFile << "use_canonical_repair=" << use_canonical_repair << endl;

        outputFile << "time_initialization=" << time_initialization << endl;
        outputFile << "time_generate_labels=" << time_generate_labels << endl;
        outputFile << "time_sortL=" << time_sortL << endl;
        outputFile << "time_canonical_repair=" << time_canonical_repair << endl;
        outputFile << "time_total=" << time_total << endl;

        outputFile << "max_bit_size=" << max_bit_size << endl;
        outputFile << "max_run_time_seconds=" << max_run_time_seconds << endl;

        outputFile << "label_size_before_canonical_repair=" << label_size_before_canonical_repair << endl;
        outputFile << "label_size_after_canonical_repair=" << label_size_after_canonical_repair << endl;
        outputFile << "canonical_repair_remove_label_ratio=" << canonical_repair_remove_label_ratio << endl;

        outputFile << "compute_label_bit_size()=" << compute_label_bit_size() << endl;

        outputFile.close();
    }
};

#include <chrono>
int _2023algo_max_second = 60;
auto begin_time = std::chrono::high_resolution_clock::now();
double max_run_time_nanosec;
string reach_limit_time_string = "reach limit time";

long long int global_query_times = 0;
long long int label_operation_times = 0;

/**
 * ort the vertices with the following priorities: first by vertex ID,
 * then by the number of hops, and finally by distance, all in ascending order.
 */
bool compare_hop_constrained_two_hop_label(hop_constrained_two_hop_label &i, hop_constrained_two_hop_label &j)
{
    if (i.hub_vertex != j.hub_vertex)
    {
        return i.hub_vertex < j.hub_vertex;
    }
    else if (i.t_e != j.t_e)
    {
        return i.t_e > j.t_e;
    }
    else if (i.hop != j.hop)
    {
        return i.hop < j.hop;
    }
    else if (i.t_s != j.t_s)
    {
        return i.t_s < j.t_s;
    }
    else
    {
        return i.distance < j.distance;
    }
}

void insert_sorted_hop_constrained_two_hop_label(std::vector<hop_constrained_two_hop_label> &input_vector, int key, int hop, weightTYPE new_distance, int t)
{
    label_operation_times++;
    int left = 0, right = input_vector.size() - 1;

    while (left <= right)
    {
        int mid = left + ((right - left) / 2);

        // Step 1: 先匹配 hub_vertex
        if (input_vector[mid].hub_vertex == key)
        {
            // Step 2: 匹配 t_e 是否为无穷大
            if (input_vector[mid].t_e == std::numeric_limits<int>::max())
            {
                // Step 3: 匹配 hop
                if (input_vector[mid].hop == hop)
                {
                    // 找到匹配的标签，先生成旧标签，再更新
                    hop_constrained_two_hop_label old_label = input_vector[mid];
                    old_label.t_e = t; // 将旧标签的 t_e 设置为当前时间 t
                    input_vector[mid].distance = new_distance;
                    input_vector[mid].t_s = t; // 更新当前标签的 t_s

                    // Step 4: 旧标签需要按排序规则插入
                    int insert_pos = mid + 1;
                    while (insert_pos < input_vector.size() &&
                           (input_vector[insert_pos].hub_vertex == old_label.hub_vertex &&
                            (input_vector[insert_pos].t_e > old_label.t_e ||
                             (input_vector[insert_pos].t_e == old_label.t_e && input_vector[insert_pos].hop < old_label.hop))))
                    {
                        insert_pos++;
                    }

                    // 插入旧标签到合适的位置
                    input_vector.insert(input_vector.begin() + insert_pos, old_label);

                    return; // 更新完成，退出函数
                }
                else if (input_vector[mid].hop < hop)
                {
                    left = mid + 1;
                }
                else
                {
                    right = mid - 1;
                }
            }
            else
            {
                right = mid - 1; // t_e 不是无穷大，继续向左找
            }
        }
        else if (input_vector[mid].hub_vertex > key)
        {
            right = mid - 1;
        }
        else
        {
            left = mid + 1;
        }
    }

    // Step 5: 如果没找到，则插入新标签
    hop_constrained_two_hop_label new_label;
    new_label.hub_vertex = key;
    new_label.hop = hop;
    new_label.distance = new_distance;
    new_label.t_s = t;
    new_label.t_e = std::numeric_limits<int>::max(); // 新标签的 t_e 是无穷大

    // 找到插入点
    input_vector.insert(input_vector.begin() + left, new_label); // 在 left 位置插入新标签
}

weightTYPE search_sorted_hop_constrained_two_hop_label(std::vector<hop_constrained_two_hop_label> &input_vector, int key, int hop)
{
    label_operation_times++;
    int left = 0, right = input_vector.size() - 1;

    while (left <= right)
    {
        int mid = left + ((right - left) / 2); // mid is between left and right (may be equal)

        // Step 1: 先匹配 hub_vertex
        if (input_vector[mid].hub_vertex == key)
        {
            // Step 2: 在匹配 hub_vertex 后，检查 t_e 是否为无穷大
            if (input_vector[mid].t_e == std::numeric_limits<int>::max())
            {
                // Step 3: 再检查 hop
                if (input_vector[mid].hop == hop)
                {
                    return input_vector[mid].distance; // 找到符合条件的标签，返回 distance
                }
                else if (input_vector[mid].hop < hop)
                {
                    left = mid + 1; // hop 太小，向右找
                }
                else
                {
                    right = mid - 1; // hop 太大，向左找
                }
            }
            else
            {
                right = mid - 1; // t_e 不是无穷大，向左继续查找
            }
        }
        else if (input_vector[mid].hub_vertex > key)
        {
            right = mid - 1; // hub_vertex 太大，向左找
        }
        else
        {
            left = mid + 1; // hub_vertex 太小，向右找
        }
    }

    return std::numeric_limits<int>::max(); // 没找到符合条件的标签，返回无穷大
}

pair<weightTYPE, int> search_sorted_hop_constrained_two_hop_label_2(std::vector<hop_constrained_two_hop_label> &input_vector, int key, int hop)
{
    label_operation_times++;
    int left = 0, right = input_vector.size() - 1;

    while (left <= right)
    {
        int mid = left + ((right - left) / 2); // mid is between left and right (may be equal);
        if (input_vector[mid].hub_vertex == key && input_vector[mid].hop == hop)
        {
            return {input_vector[mid].distance, mid};
        }
        else if (input_vector[mid].hub_vertex > key || (input_vector[mid].hub_vertex == key && input_vector[mid].hop > hop))
        {
            right = mid - 1;
        }
        else
        {
            left = mid + 1;
        }
    }

    return {std::numeric_limits<int>::max(), -1};
}

pair<weightTYPE, int> get_shortest_distance_hop_two_hop_label(std::vector<hop_constrained_two_hop_label> &input_vector, int key)
{

    label_operation_times++;
    int idx = 0, right = input_vector.size() - 1;
    weightTYPE mindis = std::numeric_limits<int>::max();
    int hop_val = 0;
    while (idx <= right)
    {
        if (input_vector[idx].hub_vertex > key)
            break;
        if (input_vector[idx].hub_vertex == key)
        {
            if (input_vector[idx].distance < mindis)
            {
                mindis = input_vector[idx].distance;
                hop_val = input_vector[idx].hop;
            }
        }
        idx++;
    }

    return {mindis, hop_val};
}

class hop_constrained_affected_label
{
public:
    int first, second, hop;
    weightTYPE dis;
    hop_constrained_affected_label() {}
    hop_constrained_affected_label(int _first, int _second, int _hop, weightTYPE _dis)
    {
        first = _first;
        second = _second;
        hop = _hop;
        dis = _dis;
    }
};

class hop_constrained_pair_label
{
public:
    int first, second;
    int hop;
    hop_constrained_pair_label(int _first, int _second, int _hop)
    {
        first = _first;
        second = _second;
        hop = _hop;
    }
    bool operator==(const hop_constrained_pair_label other) const
    {
        return (first == other.first && second == other.second && hop == other.hop);
    }
    bool operator<(const hop_constrained_pair_label other) const
    { // used to sort/search pair_label2 in set
        if (first != other.first)
            return first < other.first;
        if (second != other.second)
            return second < other.second;
        return hop < other.hop;
    }
};

class hop_constrained_label_v2
{
public:
    int hub_vertex, hop;
    weightTYPE distance;
    hop_constrained_label_v2(int _vertex, int _hop, weightTYPE _dis)
    {
        hub_vertex = _vertex;
        hop = _hop;
        distance = _dis;
    }
    bool operator==(const hop_constrained_label_v2 other) const
    {
        return (hub_vertex == other.hub_vertex && hop == other.hop && distance == other.distance);
    }
    bool operator<(const hop_constrained_label_v2 other) const
    { // used to sort/search pair_label2 in set
        if (hub_vertex != other.hub_vertex)
            return hub_vertex < other.hub_vertex;
        if (hop != other.hop)
            return hop < other.hop;
        return distance < other.distance;
    }
};

struct hop_constrained_node_for_DIFFUSE
{
    int index;
    int hop;
    weightTYPE disx;
    hop_constrained_node_for_DIFFUSE() {}
    hop_constrained_node_for_DIFFUSE(int _u, int _hop, weightTYPE _dis)
    {
        index = _u;
        hop = _hop;
        disx = _dis;
    }
}; // define the node in the queue

bool operator<(hop_constrained_node_for_DIFFUSE const &x, hop_constrained_node_for_DIFFUSE const &y)
{
    return x.disx > y.disx; // < is the max-heap; > is the min heap
}

weightTYPE hop_constrained_extract_distance(vector<vector<hop_constrained_two_hop_label>> &L, int source, int terminal, int hop_cst)
{

    /*return std::numeric_limits<int>::max() is not connected*/

    if (hop_cst < 0)
    {
        return std::numeric_limits<int>::max();
    }
    if (source == terminal)
    {
        return 0;
    }
    else if (hop_cst == 0)
    {
        return std::numeric_limits<int>::max();
    }

    int distance = std::numeric_limits<int>::max();
    auto vector1_check_pointer = L[source].begin();
    auto vector2_check_pointer = L[terminal].begin();
    auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();

    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
        {

            auto vector1_end = vector1_check_pointer;
            while (vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end != pointer_L_s_end)
            {
                vector1_end++;
            }
            auto vector2_end = vector2_check_pointer;
            while (vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end != pointer_L_t_end)
            {
                vector2_end++;
            }

            for (auto vector1_begin = vector1_check_pointer; vector1_begin != vector1_end; vector1_begin++)
            {
                // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
                for (auto vector2_begin = vector2_check_pointer; vector2_begin != vector2_end; vector2_begin++)
                {
                    // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
                    if (vector1_begin->hop + vector2_begin->hop <= hop_cst)
                    {
                        long long int dis = (long long int)vector1_begin->distance + vector2_begin->distance;
                        if (distance > dis)
                        {
                            distance = dis;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }

            vector1_check_pointer = vector1_end;
            vector2_check_pointer = vector2_end;
        }
        else if (vector1_check_pointer->hub_vertex > vector2_check_pointer->hub_vertex)
        {
            vector2_check_pointer++;
        }
        else
        {
            vector1_check_pointer++;
        }
    }

    return distance;
}

weightTYPE hop_constrained_extract_distance(vector<vector<hop_constrained_two_hop_label>> &L, int source, int terminal, int hop_cst, int t_s, int t_e)
{

    /*return std::numeric_limits<int>::max() is not connected*/

    if (hop_cst < 0)
    {
        return std::numeric_limits<int>::max();
    }
    if (source == terminal)
    {
        return 0;
    }
    else if (hop_cst == 0)
    {
        return std::numeric_limits<int>::max();
    }

    int distance = std::numeric_limits<int>::max();
    auto vector1_check_pointer = L[source].begin();
    auto vector2_check_pointer = L[terminal].begin();
    auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();

    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        while (vector1_check_pointer != pointer_L_s_end && vector1_check_pointer->t_e != std::numeric_limits<int>::max())
        {
            vector1_check_pointer++;
        }
        if (vector1_check_pointer == pointer_L_s_end)
        {
            break;
        }
        while (vector2_check_pointer != pointer_L_t_end && vector2_check_pointer->t_e != std::numeric_limits<int>::max())
        {
            vector2_check_pointer++;
        }
        if (vector2_check_pointer == pointer_L_t_end)
        {
            break;
        }
        if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
        {

            auto vector1_end = vector1_check_pointer;
            while (vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end != pointer_L_s_end)
            {
                vector1_end++;
            }
            auto vector2_end = vector2_check_pointer;
            while (vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end != pointer_L_t_end)
            {
                vector2_end++;
            }

            for (auto vector1_begin = vector1_check_pointer; vector1_begin != vector1_end; vector1_begin++)
            {
                // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
                for (auto vector2_begin = vector2_check_pointer; vector2_begin != vector2_end; vector2_begin++)
                {
                    // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
                    if (vector1_begin->hop + vector2_begin->hop <= hop_cst && max(t_s, max(vector1_begin->t_s, vector2_begin->t_s)) <= min(t_e, min(vector1_begin->t_e, vector2_begin->t_e)))
                    {
                        long long int dis = (long long int)vector1_begin->distance + vector2_begin->distance;
                        if (distance > dis)
                        {
                            distance = dis;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }

            vector1_check_pointer = vector1_end;
            vector2_check_pointer = vector2_end;
        }
        else if (vector1_check_pointer->hub_vertex > vector2_check_pointer->hub_vertex)
        {
            vector2_check_pointer++;
        }
        else
        {
            vector1_check_pointer++;
        }
    }

    return distance;
}

pair<weightTYPE, int> hop_constrained_extract_distance_and_hub(vector<vector<hop_constrained_two_hop_label>> &L, int source, int terminal, int hop_cst)
{
    global_query_times++;

    /*return std::numeric_limits<int>::max() is not connected*/

    if (hop_cst < 0)
    {
        return {std::numeric_limits<int>::max(), -1};
    }
    if (source == terminal)
    {
        return {0, -1};
    }
    else if (hop_cst == 0)
    {
        return {std::numeric_limits<int>::max(), -1};
    }

    int distance = std::numeric_limits<int>::max();
    int common_hub = -1;
    auto vector1_check_pointer = L[source].begin();
    auto vector2_check_pointer = L[terminal].begin();
    auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();

    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
        {

            auto vector1_end = vector1_check_pointer;
            while (vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end->t_e == std::numeric_limits<int>::max() && vector1_end != pointer_L_s_end)
            {
                vector1_end++;
            }
            auto vector2_end = vector2_check_pointer;
            while (vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end->t_e == std::numeric_limits<int>::max() && vector2_end != pointer_L_t_end)
            {
                vector2_end++;
            }

            for (auto vector1_begin = vector1_check_pointer; vector1_begin != vector1_end; vector1_begin++)
            {
                // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
                for (auto vector2_begin = vector2_check_pointer; vector2_begin != vector2_end; vector2_begin++)
                {
                    // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
                    if (vector1_begin->hop + vector2_begin->hop <= hop_cst)
                    {
                        long long int dis = (long long int)vector1_begin->distance + vector2_begin->distance;
                        if (distance > dis)
                        {
                            distance = dis;
                            common_hub = vector1_check_pointer->hub_vertex;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }

            vector1_check_pointer = vector1_end;
            vector2_check_pointer = vector2_end;
        }
        else if (vector1_check_pointer->hub_vertex > vector2_check_pointer->hub_vertex)
        {
            vector2_check_pointer++;
        }
        else
        {
            vector1_check_pointer++;
        }
    }

    return {distance, common_hub};
}

pair<weightTYPE, int> hop_constrained_extract_distance_and_hub(vector<vector<hop_constrained_two_hop_label>> &L, int source, int terminal, int hop_cst, int t_s, int t_e)
{
    global_query_times++;

    /*return std::numeric_limits<int>::max() is not connected*/

    if (hop_cst < 0)
    {
        return {std::numeric_limits<int>::max(), -1};
    }
    if (source == terminal)
    {
        return {0, -1};
    }
    else if (hop_cst == 0)
    {
        return {std::numeric_limits<int>::max(), -1};
    }

    int distance = std::numeric_limits<int>::max();
    int common_hub = -1;
    auto vector1_check_pointer = L[source].begin();
    auto vector2_check_pointer = L[terminal].begin();
    auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();

    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        while (vector1_check_pointer != pointer_L_s_end && vector1_check_pointer->t_e != std::numeric_limits<int>::max())
        {
            vector1_check_pointer++;
        }
        if (vector1_check_pointer == pointer_L_s_end)
        {
            break;
        }
        while (vector2_check_pointer != pointer_L_t_end && vector2_check_pointer->t_e != std::numeric_limits<int>::max())
        {
            vector2_check_pointer++;
        }
        if (vector2_check_pointer == pointer_L_t_end)
        {
            break;
        }
        if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
        {

            auto vector1_end = vector1_check_pointer;
            while (vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end->t_e == std::numeric_limits<int>::max() && vector1_end != pointer_L_s_end)
            {
                vector1_end++;
            }
            auto vector2_end = vector2_check_pointer;
            while (vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end->t_e == std::numeric_limits<int>::max() && vector2_end != pointer_L_t_end)
            {
                vector2_end++;
            }

            for (auto vector1_begin = vector1_check_pointer; vector1_begin != vector1_end; vector1_begin++)
            {
                // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
                for (auto vector2_begin = vector2_check_pointer; vector2_begin != vector2_end; vector2_begin++)
                {
                    // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
                    if (vector1_begin->hop + vector2_begin->hop <= hop_cst && max(t_s, max(vector1_begin->t_s, vector2_begin->t_s)) <= min(t_e, min(vector1_begin->t_e, vector2_begin->t_e)))
                    {
                        long long int dis = (long long int)vector1_begin->distance + vector2_begin->distance;
                        if (distance > dis)
                        {
                            distance = dis;
                            common_hub = vector1_check_pointer->hub_vertex;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }

            vector1_check_pointer = vector1_end;
            vector2_check_pointer = vector2_end;
        }
        else if (vector1_check_pointer->hub_vertex > vector2_check_pointer->hub_vertex)
        {
            vector2_check_pointer++;
        }
        else
        {
            vector1_check_pointer++;
        }
    }

    return {distance, common_hub};
}

weightTYPE hop_constrained_extract_distance_2(vector<hop_constrained_two_hop_label> &L_s, vector<hop_constrained_two_hop_label> &L_t, int hop_cst)
{

    /*return std::numeric_limits<double>::max() is not connected*/
    if (hop_cst < 0)
    {
        return std::numeric_limits<int>::max();
    }
    else if (hop_cst == 0)
    {
        return std::numeric_limits<int>::max();
    }

    int distance = std::numeric_limits<int>::max(); // if disconnected, return this large value

    auto vector1_check_pointer = L_s.begin();
    auto vector2_check_pointer = L_t.begin();
    auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
        {
            if (vector1_check_pointer->hop + vector2_check_pointer->hop <= hop_cst)
            {
                long long int dis = (long long int)vector1_check_pointer->distance + vector2_check_pointer->distance;
                if (distance > dis)
                {
                    distance = dis;
                }
            }

            vector1_check_pointer++;
        }
        else if (vector1_check_pointer->hub_vertex > vector2_check_pointer->hub_vertex)
        {
            vector2_check_pointer++;
        }
        else
        {
            vector1_check_pointer++;
        }
    }

    return distance;
}

pair<weightTYPE, int> hop_constrained_extract_distance_and_hub_2(vector<hop_constrained_two_hop_label> &L_s, vector<hop_constrained_two_hop_label> &L_t, int hop_cst)
{

    global_query_times++;

    /*return std::numeric_limits<double>::max() is not connected*/
    if (hop_cst < 0)
    {
        return {std::numeric_limits<int>::max(), -1};
    }
    else if (hop_cst == 0)
    {
        return {std::numeric_limits<int>::max(), -1};
    }

    int distance = std::numeric_limits<int>::max(); // if disconnected, return this large value
    int common_hub = -1;

    auto vector1_check_pointer = L_s.begin();
    auto vector2_check_pointer = L_t.begin();
    auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
    while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
    {
        if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
        {
            if (vector1_check_pointer->hop + vector2_check_pointer->hop <= hop_cst)
            {
                long long int dis = (long long int)vector1_check_pointer->distance + vector2_check_pointer->distance;
                if (distance > dis)
                {
                    distance = dis;
                    common_hub = vector1_check_pointer->hub_vertex;
                }
            }

            vector1_check_pointer++;
        }
        else if (vector1_check_pointer->hub_vertex > vector2_check_pointer->hub_vertex)
        {
            vector2_check_pointer++;
        }
        else
        {
            vector1_check_pointer++;
        }
    }

    return {distance, common_hub};
}
