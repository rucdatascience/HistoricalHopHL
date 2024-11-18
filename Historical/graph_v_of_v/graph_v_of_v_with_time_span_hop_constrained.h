#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"
using namespace std;

template <typename weight_type> // weight_type may be int, long long int, float, double...
class graph_v_of_v_with_time_span_hop_constrained : protected graph_v_of_v_with_time_span
{
private:
    hop_constrained_case_info *case_info;
    hop_constrained_case_info *case_info_2021;

public:
    graph_v_of_v_with_time_span_hop_constrained(hop_constrained_case_info &info, hop_constrained_case_info &info2021) : graph_v_of_v_with_time_span(), case_info(&info), case_info_2021(&info2021) {};
    graph_v_of_v_with_time_span_hop_constrained(int n, int e, weight_type weight_upper_limit, weight_type weight_lower_limit, hop_constrained_case_info &info, hop_constrained_case_info &info2021) : graph_v_of_v_with_time_span(n, e, weight_upper_limit, weight_lower_limit), case_info(&info), case_info_2021(&info2021) {};

    /*class member functions*/
    /**
     * param
     * 	@u is the source vertex,
     *  @v is the target vertex,
     * 	@startTime,
     * 	@endTime
     */
    weight_type search_shortest_path_in_period_time_naive(int u, int v, int k, int startTime, int endTime) override
    {
        weight_type res = std::numeric_limits<weight_type>::max();
        int N = this->v_num;

        boost::heap::fibonacci_heap<tuple<int, weight_type, int>, boost::heap::compare<compare_tuple<weight_type>>> queue;
        for (int queryTime = startTime; queryTime <= endTime; queryTime++)
        {
            std::vector<weight_type> dist(N, std::numeric_limits<weight_type>::max());
            std::vector<int> hop_list(N, std::numeric_limits<int>::max());
            queue.clear();
            dist[u] = 0;
            hop_list[u] = 0;
            queue.push({u, 0, 0});
            while (queue.size() > 0)
            {
                int hop = get<2>(queue.top());
                int vertexBase = get<0>(queue.top());
                weight_type currentDist = std::get<1>(queue.top());
                queue.pop();
                if (vertexBase == v)
                {
                    res = min(res, currentDist);
                }
                if (hop == k)
                {
                    continue;
                }

                for (const auto &vertices : this->ADJs[vertexBase])
                {
                    int next = vertices.first;
                    for (const auto &edge_info_time_span : vertices.second)
                    {
                        if (edge_info_time_span.startTimeLabel <= queryTime && edge_info_time_span.endTimeLabel >= queryTime)
                        {
                            weight_type newDist = currentDist + edge_info_time_span.weight;
                            if (newDist < dist[next] || (hop + 1) < hop_list[next])
                            {
                                dist[next] = newDist;
                                hop_list[next] = hop + 1;
                                queue.push({next, newDist, hop + 1});
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

    vector<graph_v_of_v<weight_type>> graph_v_of_v_with_time_span<weight_type>::graph_v_of_v_generate_random_graph_with_same_edges_of_different_weight(int change_num, int decrease_time, int increase_time, float change_ratio) override
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
        hop_constrained_two_hop_labels_generation(instance_graph, *case_info);
        hop_constrained_two_hop_labels_generation(instance_graph, *case_info_2021);
        timer.mark_time("initialize_global_values_dynamic_hop_constrained");
        timer.mark_time("====maintain process====");
        ThreadPool pool_dynamic(case_info->thread_num);
        std::vector<std::future<int>> results_dynamic;

        add_graph_time(instance_graph, 0);
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

                if (path.size() > case_info->thread_num)
                {
                    auto time1 = std::chrono::high_resolution_clock::now();
                    HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, *case_info, path, weight, pool_dynamic, results_dynamic, index);
                    auto time2 = std::chrono::high_resolution_clock::now();
                    case_info->time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                    auto time3 = std::chrono::high_resolution_clock::now();
                    HOP_WeightDecrease2021_batch(instance_graph, *case_info_2021, path, weight, pool_dynamic, results_dynamic, index);
                    auto time4 = std::chrono::high_resolution_clock::now();
                    case_info_2021->time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);
                    vector<pair<int, int>>().swap(path);
                    vector<int>().swap(weight);
                }
            }
            if (path.size() > 0)
            {
                auto time1 = std::chrono::high_resolution_clock::now();
                HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, *case_info, path, weight, pool_dynamic, results_dynamic, index);
                auto time2 = std::chrono::high_resolution_clock::now();
                case_info->time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                auto time3 = std::chrono::high_resolution_clock::now();
                HOP_WeightDecrease2021_batch(instance_graph, *case_info_2021, path, weight, pool_dynamic, results_dynamic, index);
                auto time4 = std::chrono::high_resolution_clock::now();
                case_info_2021->time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);
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
                    HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, *case_info, path, weight, pool_dynamic, results_dynamic, index);
                    auto time2 = std::chrono::high_resolution_clock::now();
                    case_info->time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                    auto time3 = std::chrono::high_resolution_clock::now();
                    HOP_WeightIncrease2021_batch(instance_graph, *case_info_2021, path, weight, pool_dynamic, results_dynamic, index);
                    auto time4 = std::chrono::high_resolution_clock::now();
                    case_info_2021->time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                    vector<pair<int, int>>().swap(path);
                    vector<int>().swap(weight);
                }
            }
            if (path.size() > 0)
            {
                auto time1 = std::chrono::high_resolution_clock::now();
                HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, *case_info, path, weight, pool_dynamic, results_dynamic, index);
                auto time2 = std::chrono::high_resolution_clock::now();
                case_info->time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                auto time3 = std::chrono::high_resolution_clock::now();
                HOP_WeightIncrease2021_batch(instance_graph, *case_info_2021, path, weight, pool_dynamic, results_dynamic, index);
                auto time4 = std::chrono::high_resolution_clock::now();
                case_info_2021->time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                vector<pair<int, int>>().swap(path);
                vector<int>().swap(weight);
            }
            // case_info.print_L();
            res.push_back(instance_graph);
            ++index;
        }
        return res;
    };

    vector<graph_v_of_v<weight_type>> graph_v_of_v_with_time_span<weight_type>::txt_read(std::string save_name)
    {
        this->clear();
        std::string line_content;
        int current_time = -1;
        vector<graph_v_of_v<weight_type>> res;
        graph_v_of_v<weight_type> instance_graph;
        ThreadPool pool_dynamic(case_info->thread_num);
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
                    ADJs.resize(std::stoi(Parsed_content[1]));
                    initialize_global_values_dynamic_hop_constrained(this->v_num, case_info->thread_num, case_info->upper_k);
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
                        if (path_decrease.size() >= case_info->thread_num)
                        {
                            this->process(instance_graph, path_decrease, weight_decrease, current_time);

                            auto time1 = std::chrono::high_resolution_clock::now();
                            HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, *case_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, current_time);
                            auto time2 = std::chrono::high_resolution_clock::now();
                            case_info->time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                            auto time3 = std::chrono::high_resolution_clock::now();
                            HOP_WeightDecrease2021_batch(instance_graph, *case_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, current_time);
                            auto time4 = std::chrono::high_resolution_clock::now();
                            case_info_2021->time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                            vector<pair<int, int>>().swap(path_decrease);
                            vector<int>().swap(weight_decrease);
                        }
                        if (path_increase.size() >= case_info.thread_num)
                        {
                            this->process(instance_graph, path_increase, weight_increase, current_time);

                            auto time1 = std::chrono::high_resolution_clock::now();
                            HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, *case_info, path_increase, old_weight_increase, pool_dynamic, results_dynamic, current_time);
                            auto time2 = std::chrono::high_resolution_clock::now();
                            case_info->time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                            auto time3 = std::chrono::high_resolution_clock::now();
                            HOP_WeightIncrease2021_batch(instance_graph, *case_info_2021, path_increase, old_weight_increase, pool_dynamic, results_dynamic, current_time);
                            auto time4 = std::chrono::high_resolution_clock::now();
                            case_info_2021->time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

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
                        HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, *case_info, path_decrease, weight_decrease, pool_dynamic, results_dynamic, current_time);
                        auto time2 = std::chrono::high_resolution_clock::now();
                        case_info->time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                        auto time3 = std::chrono::high_resolution_clock::now();
                        HOP_WeightDecrease2021_batch(instance_graph, *case_info_2021, path_decrease, weight_decrease, pool_dynamic, results_dynamic, current_time);
                        auto time4 = std::chrono::high_resolution_clock::now();
                        case_info_2021->time_decrease.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                        vector<pair<int, int>>().swap(path_decrease);
                        vector<int>().swap(weight_decrease);
                    }
                    if (path_increase.size() > 0)
                    {
                        this->process(instance_graph, path_increase, weight_increase, current_time);

                        auto time1 = std::chrono::high_resolution_clock::now();
                        HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, *case_info, path_increase, old_weight_increase, pool_dynamic, results_dynamic, current_time);
                        auto time2 = std::chrono::high_resolution_clock::now();
                        case_info->time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
                        auto time3 = std::chrono::high_resolution_clock::now();
                        HOP_WeightIncrease2021_batch(instance_graph, *case_info_2021, path_increase, old_weight_increase, pool_dynamic, results_dynamic, current_time);
                        auto time4 = std::chrono::high_resolution_clock::now();
                        case_info_2021->time_increase.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(time4 - time3).count() / 1e9);

                        vector<pair<int, int>>().swap(path_increase);
                        vector<int>().swap(weight_increase);
                        vector<int>().swap(old_weight_increase);
                    }
                    if (current_time >= 0)
                    {
                        if (current_time == 0)
                        {
                            hop_constrained_two_hop_labels_generation(instance_graph, *case_info);
                            hop_constrained_two_hop_labels_generation(instance_graph, *case_info_2021);
                            timer.mark_time("initialize the 2-hop label");
                            add_graph_time(instance_graph, 0);
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