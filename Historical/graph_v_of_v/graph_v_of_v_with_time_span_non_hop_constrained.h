#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"
#include "CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_maintain_two_hop_labels.h"
#include "CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_maintain_PLL.h"
#include "CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_WeightDecreaseMaintenance_improv_batch.h"
#include "CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_WeightDecrease2021_batch.h"

using namespace std;
using namespace nonHop;
template <typename weight_type>
struct compare_pair
{
    bool operator()(const pair<int, weight_type> &lhs, const pair<int, weight_type> &rhs) const
    {
        if (lhs.first == rhs.first)
        {
            return lhs.second > rhs.second;
        }
        return lhs.first > rhs.first;
    }
};
template <typename weight_type> // weight_type may be int, long long int, float, double...
class graph_v_of_v_with_time_span_non_hop_constrained : public graph_v_of_v_with_time_span<weight_type>
{
private:
public:
    graph_v_of_v_with_time_span_non_hop_constrained() : graph_v_of_v_with_time_span<weight_type>() {};
    graph_v_of_v_with_time_span_non_hop_constrained(int n, int e, weight_type weight_upper_limit, weight_type weight_lower_limit) : graph_v_of_v_with_time_span<weight_type>(n, e, weight_upper_limit, weight_lower_limit) {};

    /**
     * param
     * 	@u is the source vertex,
     *  @v is the target vertex,
     * 	@startTime,
     * 	@endTime
     */
    weight_type search_shortest_path_in_period_time_naive(int u, int v, int startTime, int endTime)
    {
        weight_type res = std::numeric_limits<weight_type>::max();
        int N = this->v_num;

        boost::heap::fibonacci_heap<pair<int, weight_type>, boost::heap::compare<compare_pair<weight_type>>> queue;
        for (int queryTime = startTime; queryTime <= endTime; queryTime++)
        {
            std::vector<weight_type> dist(N, std::numeric_limits<weight_type>::max());
            queue.clear();
            dist[u] = 0;
            queue.push({u, 0});
            while (queue.size() > 0)
            {
                int vertexBase = queue.top().first;
                weight_type currentDist = queue.top().second;
                queue.pop();
                if (vertexBase == v)
                {
                    res = min(res, currentDist);
                }
                for (const auto &vertices : this->ADJs[vertexBase])
                {
                    int next = vertices.first;
                    for (const auto &edge_info_time_span : vertices.second)
                    {
                        if (edge_info_time_span.startTimeLabel <= queryTime && edge_info_time_span.endTimeLabel >= queryTime)
                        {
                            weight_type newDist = currentDist + edge_info_time_span.weight;
                            if (newDist < dist[next])
                            {
                                dist[next] = newDist;
                                queue.push({next, newDist});
                            }
                            break;
                        }
                    }
                }
            }
            // cout << "naive" << res << endl;
        }
        return res;
    };

    vector<graph_v_of_v<weight_type>> graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(int change_num, int decrease_time, int increase_time, float change_ratio, two_hop_case_info &case_info, two_hop_case_info &case_info_2021)
    {
        if (change_num < 0)
        {
            cout << "the change_num should be greater than or equal to 0" << endl;
        }
        this->time_max = change_num;
        boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};
        graph_v_of_v<weight_type> instance_graph;
        instance_graph = graph_v_of_v_generate_random_graph<weight_type>(this->v_num, this->e_num, this->weight_dis.min(), this->weight_dis.max(), 1, boost_random_time_seed);
        vector<int> is_mock(instance_graph.size());
        for (int i = 0; i < instance_graph.size(); i++)
        {
            is_mock[i] = false;
        }
        instance_graph = graph_v_of_v_update_vertexIDs_by_degrees_large_to_small_mock(instance_graph, is_mock);
        timer.mark_time("====time 0====");
        PLL(instance_graph, case_info);
        PLL(instance_graph, case_info_2021);
        timer.mark_time("initialize_global_values_dynamic_non_hop_constrained");
        timer.mark_time("====maintain process====");
        ThreadPool pool_dynamic(case_info.thread_num);
        std::vector<std::future<int>> results_dynamic;

        this->add_graph_time(instance_graph, 0);
        vector<graph_v_of_v<weight_type>> res;
        res.push_back(instance_graph);
        uniform_int_distribution<> dis(0, this->v_num);
        int index = 1;
        int N = instance_graph.ADJs.size();
        while (index <= change_num)
        {
            int current_decrease_time = decrease_time;
            int current_increase_time = increase_time;
            timer.mark_time("====time " + std::to_string(index) + "====");
            vector<pair<int, int>> path;
            vector<int> weight;
            int i, j;
            while (current_decrease_time > 0)
            {
                current_decrease_time--;
                i = dis(boost_random_time_seed);
                if (instance_graph.ADJs[i].size() == 0)
                {
                    continue;
                }
                uniform_int_distribution<> dis_inner(0, instance_graph.ADJs[i].size() - 1);
                j = dis_inner(boost_random_time_seed);
                int next_value = (instance_graph.ADJs[i][j].second) * (1 - change_ratio);
                if (next_value > instance_graph.ADJs[i][j].second)
                {
                    cout << "error in decrease" << endl;
                }
                if (next_value == instance_graph.ADJs[i][j].second || next_value < 0)
                {
                    continue;
                }
                instance_graph.add_edge(i, instance_graph.ADJs[i][j].first, next_value);
                this->add_edge(i, instance_graph.ADJs[i][j].first, next_value, index);
                path.push_back({i, instance_graph.ADJs[i][j].first});
                weight.push_back(next_value);

                if (path.size() > case_info.thread_num)
                {
                    auto time1 = std::chrono::high_resolution_clock::now();
                    nonHOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, path, weight, pool_dynamic, results_dynamic, index);
                    auto time2 = std::chrono::high_resolution_clock::now();
                    case_info.time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                    auto time3 = std::chrono::high_resolution_clock::now();
                    nonHOP_WeightDecrease2021_batch(instance_graph, case_info_2021, path, weight, pool_dynamic, results_dynamic, index);
                    auto time4 = std::chrono::high_resolution_clock::now();
                    case_info_2021.time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);
                    vector<pair<int, int>>().swap(path);
                    vector<int>().swap(weight);
                }
            }
            if (path.size() > 0)
            {
                auto time1 = std::chrono::high_resolution_clock::now();
                nonHOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, path, weight, pool_dynamic, results_dynamic, index);
                auto time2 = std::chrono::high_resolution_clock::now();
                case_info.time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                auto time3 = std::chrono::high_resolution_clock::now();
                nonHOP_WeightDecrease2021_batch(instance_graph, case_info_2021, path, weight, pool_dynamic, results_dynamic, index);
                auto time4 = std::chrono::high_resolution_clock::now();
                case_info_2021.time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);
                vector<pair<int, int>>().swap(path);
                vector<int>().swap(weight);
            }

            while (current_increase_time > 0)
            {
                current_increase_time--;
                i = dis(boost_random_time_seed);
                if (instance_graph.ADJs[i].size() == 0)
                {
                    continue;
                }
                uniform_int_distribution<> dis_inner(0, instance_graph.ADJs[i].size() - 1);
                j = dis_inner(boost_random_time_seed);
                int old_value = instance_graph.ADJs[i][j].second;
                int next_value = (instance_graph.ADJs[i][j].second) * (1 + change_ratio);
                if (next_value < instance_graph.ADJs[i][j].second)
                {
                    cout << "error in increase" << endl;
                }
                if (next_value == instance_graph.ADJs[i][j].second)
                {
                    continue;
                }
                instance_graph.add_edge(i, instance_graph.ADJs[i][j].first, next_value);
                this->add_edge(i, instance_graph.ADJs[i][j].first, next_value, index);
                path.push_back({i, instance_graph.ADJs[i][j].first});
                weight.push_back(old_value);

                if (path.size() > case_info.thread_num)
                {
                    auto time1 = std::chrono::high_resolution_clock::now();
                    // HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, case_info, path, weight, pool_dynamic, results_dynamic, index);
                    auto time2 = std::chrono::high_resolution_clock::now();
                    case_info.time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                    auto time3 = std::chrono::high_resolution_clock::now();
                    // HOP_WeightIncrease2021_batch(instance_graph, case_info_2021, path, weight, pool_dynamic, results_dynamic, index);
                    auto time4 = std::chrono::high_resolution_clock::now();
                    case_info_2021.time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                    vector<pair<int, int>>().swap(path);
                    vector<int>().swap(weight);
                }
            }
            if (path.size() > 0)
            {
                auto time1 = std::chrono::high_resolution_clock::now();
                // HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, case_info, path, weight, pool_dynamic, results_dynamic, index);
                auto time2 = std::chrono::high_resolution_clock::now();
                case_info.time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                auto time3 = std::chrono::high_resolution_clock::now();
                // HOP_WeightIncrease2021_batch(instance_graph, case_info_2021, path, weight, pool_dynamic, results_dynamic, index);
                auto time4 = std::chrono::high_resolution_clock::now();
                case_info_2021.time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                vector<pair<int, int>>().swap(path);
                vector<int>().swap(weight);
            }
            // case_info.print_L();
            res.push_back(instance_graph);
            ++index;
        }
        return res;
    };

    vector<graph_v_of_v<weight_type>> txt_read(std::string save_name, two_hop_case_info &case_info, two_hop_case_info &case_info_2021)
    {
        this->clear();
        std::string line_content;
        int current_time = -1;
        vector<graph_v_of_v<weight_type>> res;
        graph_v_of_v<weight_type> instance_graph;
        ThreadPool pool_dynamic(case_info.thread_num);
        std::vector<std::future<int>> results_dynamic;
        vector<pair<int, int>> path_decrease;
        vector<weight_type> weight_decrease;
        vector<pair<int, int>> path_increase;
        vector<weight_type> weight_increase;
        vector<weight_type> old_weight_increase;

        std::ifstream myfile(save_name); // open the file
        if (myfile.is_open())            // if the file is opened successfully
        {
            while (getline(myfile, line_content)) // read file line by line
            {
                std::vector<std::string> Parsed_content = parse_string(line_content, " ");

                if (!Parsed_content[0].compare("|V|=")) // when it's equal, compare returns 0
                {
                    this->v_num = std::stoi(Parsed_content[1]);
                    instance_graph.ADJs.resize(this->v_num);
                    this->ADJs.resize(std::stoi(Parsed_content[1]));
                    initialize_global_values_dynamic(this->v_num, case_info.thread_num);
                }
                else if (!Parsed_content[0].compare("|E|="))
                {
                    this->e_num = std::stoi(Parsed_content[1]);
                }
                else if (!Parsed_content[0].compare("|time|="))
                {
                    this->time_max = std::stoi(Parsed_content[1]);
                }
                else if (!Parsed_content[0].compare("time"))
                {
                    current_time = std::stoi(Parsed_content[1]);
                    timer.mark_time("====time " + Parsed_content[1] + "====");
                    // std::cout << "====time " << Parsed_content[1] << "====" << endl;
                }
                else if (!Parsed_content[0].compare("Edge"))
                {
                    int v1 = std::stoi(Parsed_content[1]);
                    int v2 = std::stoi(Parsed_content[2]);
                    weight_type ec = std::stod(Parsed_content[3]);
                    if (current_time == 0)
                    {
                        // initiate the graph
                        instance_graph.add_edge(v1, v2, ec);
                    }
                    else
                    {
                        // maintain the label
                        int old_ec = instance_graph.edge_weight(v1, v2);
                        if (old_ec > ec)
                        {
                            path_decrease.push_back({v1, v2});
                            weight_decrease.push_back(ec);
                        }
                        else if (old_ec < ec)
                        {
                            path_increase.push_back({v1, v2});
                            weight_increase.push_back(ec);
                            old_weight_increase.push_back(old_ec);
                        }
                        if (path_decrease.size() >= case_info.thread_num)
                        {
                            this->process(instance_graph, path_decrease, weight_decrease, current_time);

                            auto time1 = std::chrono::high_resolution_clock::now();
                            nonHOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, current_time);
                            auto time2 = std::chrono::high_resolution_clock::now();
                            case_info.time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                            auto time3 = std::chrono::high_resolution_clock::now();
                            nonHOP_WeightDecrease2021_batch(instance_graph, case_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, current_time);
                            auto time4 = std::chrono::high_resolution_clock::now();
                            case_info_2021.time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                            vector<pair<int, int>>().swap(path_decrease);
                            vector<int>().swap(weight_decrease);
                        }
                        if (path_increase.size() >= case_info.thread_num)
                        {
                            this->process(instance_graph, path_increase, weight_increase, current_time);

                            auto time1 = std::chrono::high_resolution_clock::now();
                            // HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, case_info, path_increase, old_weight_increase, pool_dynamic, results_dynamic, current_time);
                            auto time2 = std::chrono::high_resolution_clock::now();
                            case_info.time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                            auto time3 = std::chrono::high_resolution_clock::now();
                            // HOP_WeightIncrease2021_batch(instance_graph, case_info_2021, path_increase, old_weight_increase, pool_dynamic, results_dynamic, current_time);
                            auto time4 = std::chrono::high_resolution_clock::now();
                            case_info_2021.time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                            vector<pair<int, int>>().swap(path_increase);
                            vector<int>().swap(weight_increase);
                            vector<int>().swap(old_weight_increase);
                        }
                    }
                }
                else if (Parsed_content.size() == 1 && Parsed_content[0] == "")
                {
                    if (path_decrease.size() > 0)
                    {
                        this->process(instance_graph, path_decrease, weight_decrease, current_time);

                        auto time1 = std::chrono::high_resolution_clock::now();
                        nonHOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, current_time);
                        auto time2 = std::chrono::high_resolution_clock::now();
                        case_info.time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                        auto time3 = std::chrono::high_resolution_clock::now();
                        nonHOP_WeightDecrease2021_batch(instance_graph, case_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, current_time);
                        auto time4 = std::chrono::high_resolution_clock::now();
                        case_info_2021.time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                        vector<pair<int, int>>().swap(path_decrease);
                        vector<int>().swap(weight_decrease);
                    }
                    if (path_increase.size() > 0)
                    {
                        this->process(instance_graph, path_increase, weight_increase, current_time);

                        auto time1 = std::chrono::high_resolution_clock::now();
                        // HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, case_info, path_increase, old_weight_increase, pool_dynamic, results_dynamic, current_time);
                        auto time2 = std::chrono::high_resolution_clock::now();
                        case_info.time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                        auto time3 = std::chrono::high_resolution_clock::now();
                        // HOP_WeightIncrease2021_batch(instance_graph, case_info_2021, path_increase, old_weight_increase, pool_dynamic, results_dynamic, current_time);
                        auto time4 = std::chrono::high_resolution_clock::now();
                        case_info_2021.time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                        vector<pair<int, int>>().swap(path_increase);
                        vector<int>().swap(weight_increase);
                        vector<int>().swap(old_weight_increase);
                    }
                    if (current_time >= 0)
                    {
                        if (current_time == 0)
                        {
                            PLL(instance_graph, case_info);
                            PLL(instance_graph, case_info_2021);
                            timer.mark_time("initialize the 2-hop label");
                            this->add_graph_time(instance_graph, 0);
                        }
                        res.push_back(instance_graph);
                    }
                }
            }
            myfile.close(); // close the file
            return res;
        }
        else
        {
            std::cout << "Unable to open file " << save_name << std::endl
                      << "Please check the file location or file name." << std::endl; // throw an error message
            getchar();                                                                // keep the console window
            exit(1);                                                                  // end the program
        }
    };
};

template <typename weight_type>
weight_type dijkstra(graph_v_of_v<weight_type> &graph, int u, int v)
{
    std::vector<weight_type> dist(graph.size(), std::numeric_limits<weight_type>::max());
    boost::heap::fibonacci_heap<std::pair<int, weight_type>, boost::heap::compare<compare_pair<weight_type>>> queue;

    dist[u] = 0;
    queue.push({u, 0});
    int res = __INT_MAX__;
    while (!queue.empty())
    {
        auto top = queue.top();
        int vertexBase = std::get<0>(top);
        weight_type currentDist = std::get<1>(top);
        queue.pop();

        if (vertexBase == v)
        {
            res = min(res, currentDist);
        }
        for (const auto &edge : graph[vertexBase])
        {
            int next = edge.first;
            weight_type weight = edge.second;
            weight_type newDist = currentDist + weight;

            if (newDist < dist[next])
            {
                dist[next] = newDist;
                queue.push({next, newDist});
            }
        }
    }
    return res;
}

template <typename weight_type>
int dijkstra_iterator(vector<graph_v_of_v<weight_type>> list, int u, int v)
{
    int res = INT_MAX;
    auto start_time = std::chrono::high_resolution_clock::now();
    for (graph_v_of_v<int> graph : list)
    {
        res = min(res, dijkstra(graph, u, v));
        // cout << "dijkstra" << res << endl;
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    double runtime_n_iterate_dijkstra = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - start_time).count() / 1e9;
    std::cout << "dijkstra query time" << runtime_n_iterate_dijkstra << endl;
    return res;
};