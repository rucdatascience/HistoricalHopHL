#pragma once

/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------

#include <iostream>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_test.h>


int main()
{
    HOP_maintain_test();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/rucgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)


*/
#include <CPU/graph_v_of_v/graph_v_of_v.h>
#include <CPU/graph_v_of_v/graph_v_of_v_update_vertexIDs_by_degrees_large_to_small.h>
#include <CPU/graph_v_of_v/graph_v_of_v_generate_random_graph.h>
#include <CPU/graph_v_of_v/graph_v_of_v_shortest_paths.h>
#include <CPU/graph_v_of_v/graph_v_of_v_hop_constrained_shortest_distance.h>
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels_generation.h>
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightDecreaseMaintenance_improv_batch.h>
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightIncreaseMaintenance_improv_batch.h>
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightDecrease2021_batch.h>
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_WeightIncrease2021_batch.h>
#include <CPU/text_mining/print_items.h>
#include <CPU/text_mining/binary_save_read_vector.h>

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};

void add_vertex_groups(graph_v_of_v<int> &instance_graph, int group_num)
{

    double dummy_edge_probability = 0.2;
    boost::random::uniform_int_distribution<> dist{static_cast<int>(1), static_cast<int>(100)};

    /* get the number of vertex in the graph */
    int N = instance_graph.size();

    /* add the number of vertices */
    instance_graph.ADJs.resize(N + group_num);
    for (int i = N; i < N + group_num; i++)
    {
        for (int j = 0; j < N; j++)
        {
            /* in 20% of cases, connect the group vertex to the existing vertex with a dummy edge */
            if ((double)dist(boost_random_time_seed) / 100 < dummy_edge_probability)
            {
                instance_graph.add_edge(i, j, 1e6); // add a dummy edge
            }
        }
    }
}

void hop_constrained_check_correctness(hop_constrained_case_info &case_info, graph_v_of_v<int> &instance_graph,
                                       int iteration_source_times, int iteration_terminal_times, int upper_k)
{
    // boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};

    boost::random::uniform_int_distribution<> vertex_range{static_cast<int>(0), static_cast<int>(instance_graph.ADJs.size() - 1)};
    boost::random::uniform_int_distribution<> hop_range{static_cast<int>(1), static_cast<int>(upper_k)};

    for (int yy = 0; yy < iteration_source_times; yy++)
    {
        int source = vertex_range(boost_random_time_seed);

        while (is_mock[source])
        {
            source = vertex_range(boost_random_time_seed);
        }

        std::vector<int> distances(instance_graph.size());

        int hop_cst = hop_range(boost_random_time_seed);

        graph_v_of_v_hop_constrained_shortest_distance(instance_graph, source, hop_cst, distances);

        for (int xx = 0; xx < iteration_terminal_times; xx++)
        {
            int terminal = vertex_range(boost_random_time_seed);

            while (is_mock[terminal])
            {
                terminal = vertex_range(boost_random_time_seed);
            }

            int query_dis = hop_constrained_extract_distance(case_info.L, source, terminal, hop_cst);

            if (abs(query_dis - distances[terminal]) > 1e-4 && (query_dis < TwoM_value || distances[terminal] < TwoM_value))
            {
                hop_constrained_case_info case_info2;

                case_info2.upper_k = case_info.upper_k;
                case_info2.use_2M_prune = case_info.use_2M_prune;
                case_info2.use_rank_prune = case_info.use_rank_prune;
                case_info2.use_2023WWW_generation = case_info.use_2023WWW_generation;
                case_info2.use_canonical_repair = case_info.use_canonical_repair;
                case_info2.max_run_time_seconds = case_info.max_run_time_seconds;
                case_info2.thread_num = case_info.thread_num;

                hop_constrained_two_hop_labels_generation(instance_graph, case_info2);
                // instance_graph.print();
                // case_info.print_L();
                cout << "source = " << source << endl;
                cout << "terminal = " << terminal << endl;
                cout << "hop_cst = " << hop_cst << endl;
                cout << "query_dis = " << query_dis << endl;
                cout << "distances[terminal] = " << distances[terminal] << endl;
                cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
                cout << "source vector:" << endl;
                case_info.print_L_vk(source);
                cout << "terminal vector:" << endl;
                case_info.print_L_vk(terminal);
                cout << "source correct label:" << endl;
                case_info2.print_L_vk(source);
                cout << "terminal correct label:" << endl;
                case_info2.print_L_vk(terminal);
                getchar();
            }
        }
    }
}

void graph_change_and_label_maintenance(graph_v_of_v<int> &instance_graph, hop_constrained_case_info &case_info,
                                        int V, int weightIncrease_time, int weightDecrease_time, double weightChange_ratio, double &avg_maintain_time, int batch_size, int insert_and_delete_edge, int ec_min, int ec_max)
{
    ThreadPool pool_dynamic(case_info.thread_num);
    std::vector<std::future<int>> results_dynamic;

    while (weightIncrease_time + weightDecrease_time)
    {

        /*randomly select an edge*/
        vector<pair<int, int>> selected_edge_vec;
        vector<int> selected_edge_weight_vec;
        vector<int> new_edge_weight_vec;
        /* generate random vertex and weight which be stored to vector */
        if (!insert_and_delete_edge)
        {
            int ct = 0;
            while (ct < batch_size)
            {
                boost::random::uniform_int_distribution<> dist_v1{static_cast<int>(0), static_cast<int>(V - 1)};
                int v1 = dist_v1(boost_random_time_seed);

                if (instance_graph[v1].size() > 0)
                {
                    boost::random::uniform_int_distribution<> dist_v2{static_cast<int>(0), static_cast<int>(instance_graph[v1].size() - 1)};
                    int rand = dist_v2(boost_random_time_seed);
                    int v2 = instance_graph[v1][rand].first;
                    int wt = instance_graph[v1][rand].second;
                    if (is_mock[v2] || wt >= 1e6)
                        continue;

                    selected_edge_vec.push_back({v1, v2});
                    selected_edge_weight_vec.push_back(wt);
                    ct++;
                }
            }
        }

        /*change weight*/
        if (weightIncrease_time >= weightDecrease_time)
        {
            weightIncrease_time--;
            bool next_flag = false;
            if (!insert_and_delete_edge)
            {
                for (int i = 0; i < batch_size; i++)
                {
                    int new_ec = selected_edge_weight_vec[i] * (1 + weightChange_ratio);
                    if (new_ec > 1e6 || new_ec == selected_edge_weight_vec[i])
                    {
                        weightIncrease_time++;
                        next_flag = true;
                        break;
                    }

                    new_edge_weight_vec.push_back(new_ec);
                }
                if (next_flag)
                    continue;
            }
            else
            {
                int ct = 0;
                while (ct < batch_size)
                {
                    boost::random::uniform_int_distribution<> dist_v1{static_cast<int>(0), static_cast<int>(V - 1)};
                    int v1 = dist_v1(boost_random_time_seed);

                    if (instance_graph[v1].size() > 0)
                    {
                        boost::random::uniform_int_distribution<> dist_v2{static_cast<int>(0), static_cast<int>(instance_graph[v1].size() - 1)};
                        int rand = dist_v2(boost_random_time_seed);
                        int v2 = instance_graph[v1][rand].first;
                        int wt = instance_graph[v1][rand].second;
                        if (is_mock[v2] || wt >= 1e6)
                            continue;
                        if (v1 > v2)
                            swap(v1, v2);
                        if (find(selected_edge_vec.begin(), selected_edge_vec.end(), make_pair(v1, v2)) != selected_edge_vec.end())
                            continue;
                        selected_edge_vec.push_back({v1, v2});
                        selected_edge_weight_vec.push_back(wt);
                        new_edge_weight_vec.push_back(MAX_VALUE);
                        // std::cout<<v1<<" "<<v2<<" : "<<wt<<" -> "<<MAX_VALUE<<std::endl;
                        ct++;
                    }
                }
            }

            for (int i = 0; i < batch_size; i++)
            {
                instance_graph.add_edge(selected_edge_vec[i].first, selected_edge_vec[i].second, new_edge_weight_vec[i]);
            }

            auto begin = std::chrono::high_resolution_clock::now();

            /*maintain labels*/
            // HOP_WeightIncrease2021_batch(instance_graph, case_info, selected_edge_vec, selected_edge_weight_vec, pool_dynamic, results_dynamic);
            HOP_WeightIncreaseMaintenance_improv_batch(instance_graph, case_info, selected_edge_vec, selected_edge_weight_vec, pool_dynamic, results_dynamic);
            //    cout << "1ec change " << selected_edge.first << " " << selected_edge.second << " " << new_ec << endl;
            //    mm.print_L();
            //    mm.print_PPR();

            auto end = std::chrono::high_resolution_clock::now();
            avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        }
        else
        {
            weightDecrease_time--;
            bool next_flag = false;
            if (!insert_and_delete_edge)
            {
                for (int i = 0; i < batch_size; i++)
                {
                    // int new_ec = 0; // delete an edge
                    weightTYPE new_ec = selected_edge_weight_vec[i] * (1 - weightChange_ratio);
                    if (new_ec < 1)
                    {
                        weightDecrease_time++;
                        next_flag = true;
                        break;
                    }

                    new_edge_weight_vec.push_back(new_ec);
                }
                if (next_flag)
                    continue;
            }
            else
            {
                int ct = 0;
                while (ct < batch_size)
                {
                    boost::random::uniform_int_distribution<> dist_v1{static_cast<int>(0), static_cast<int>(V - 1)};
                    int v1 = dist_v1(boost_random_time_seed);

                    if (instance_graph[v1].size() < instance_graph.size())
                    {
                        boost::random::uniform_int_distribution<> dist_v2{static_cast<int>(0), static_cast<int>(V - 1)};
                        int v2 = dist_v2(boost_random_time_seed);
                        int wt = instance_graph.edge_weight(v1, v2);
                        if (is_mock[v2] || wt <= 1e6)
                            continue;
                        if (v1 > v2)
                            swap(v1, v2);
                        if (find(selected_edge_vec.begin(), selected_edge_vec.end(), make_pair(v1, v2)) != selected_edge_vec.end())
                            continue;
                        selected_edge_vec.push_back({v1, v2});
                        selected_edge_weight_vec.push_back(wt);
                        boost::random::uniform_int_distribution<> dist{static_cast<int>(ec_min), static_cast<int>(ec_max)};
                        int new_ec = dist(boost_random_time_seed);
                        new_edge_weight_vec.push_back(new_ec);
                        // std::cout<<v1<<" "<<v2<<" : "<<wt<<" -> "<<new_ec<<std::endl;
                        ct++;
                    }
                }
            }

            for (int i = 0; i < batch_size; i++)
            {
                instance_graph.add_edge(selected_edge_vec[i].first, selected_edge_vec[i].second, new_edge_weight_vec[i]);
            }
            auto begin = std::chrono::high_resolution_clock::now();

            int t = 1;
            /*maintain labels*/
            // HOP_WeightDecrease2021_batch(instance_graph, case_info, selected_edge_vec, new_edge_weight_vec, pool_dynamic, results_dynamic);
            HOP_WeightDecreaseMaintenance_improv_batch(instance_graph, case_info, selected_edge_vec, new_edge_weight_vec, pool_dynamic, results_dynamic, t);
            //       cout << "2ec change " << selected_edge.first << " " << selected_edge.second << " " << new_ec << endl;
            //       mm.print_L();
            //       mm.print_PPR();

            auto end = std::chrono::high_resolution_clock::now();
            avg_maintain_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        }
    }
}

void HOP_maintain_test()
{

    /* problem parameters */
    int iteration_graph_times = 1e4, iteration_source_times = 10, iteration_terminal_times = 10;
    int V = 100, E = 500, group_num = 0;
    int ec_min = 1, ec_max = 20;

    int weightIncrease_time = 20, weightDecrease_time = 0;
    double weightChange_ratio = 0.2;
    int batch_size = 10;

    /* result info */
    double avg_index_time = 0, avg_index_size_per_v = 0, avg_maintain_time = 0;
    double avg_canonical_repair_remove_label_ratio = 0;

    bool insert_and_delete_edge = 1;

    bool generate_new_graph = 1;

    /* hop bounded info */
    hop_constrained_case_info mm;
    mm.upper_k = 5;
    mm.max_bit_size = 6e9;
    mm.use_2M_prune = 1;
    mm.use_rank_prune = 1; // set true
    mm.use_2023WWW_generation = 0;
    mm.use_canonical_repair = 0;
    mm.max_run_time_seconds = 1e2;
    mm.thread_num = 10;

    /* iteration */
    for (int i = 0; i < iteration_graph_times; i++)
    {
        cout << "iteration " << i << endl;

        graph_v_of_v<int> instance_graph;
        if (generate_new_graph == 1)
        {
            /* generate the graph*/
            instance_graph = graph_v_of_v_generate_random_graph<int>(V, E, ec_min, ec_max, 1, boost_random_time_seed);
            /*add vertex groups. In fact ,these steps will never be executed*/
            if (group_num > 0)
            {
                add_vertex_groups(instance_graph, group_num);
            }
            is_mock.resize(V + group_num);
            for (int j = 0; j < V; j++)
            {
                is_mock[j] = false;
            }
            for (int j = 0; j < group_num; j++)
            {
                is_mock[V + j] = true;
            }
            /* sort vertices in descending order by degree and re-generate the graph*/
            instance_graph = graph_v_of_v_update_vertexIDs_by_degrees_large_to_small_mock(instance_graph, is_mock); // sort vertices
            instance_graph.txt_save("simple_iterative_tests.txt");
            binary_save_vector("simple_iterative_tests_is_mock.txt", is_mock);
        }
        else
        {
            instance_graph.txt_read("simple_iterative_tests.txt");
            binary_read_vector("simple_iterative_tests_is_mock.txt", is_mock);
        }

        // instance_graph.print();

        auto begin = std::chrono::high_resolution_clock::now();
        try
        {
            /* generate the two hop labels with k hop constraints.the k is 5 in this program*/
            hop_constrained_two_hop_labels_generation(instance_graph, mm);
            if (0)
            {
                // graph_hash_of_mixed_weighted_print(instance_graph);
                instance_graph.print();
                mm.print_L();
                mm.print_PPR();
            }
        }
        catch (string s)
        {
            cout << s << endl;
            hop_constrained_clear_global_values();
            continue;
        }
        auto end = std::chrono::high_resolution_clock::now();
        double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
        avg_index_time = avg_index_time + runningtime / iteration_graph_times;
        avg_canonical_repair_remove_label_ratio += mm.canonical_repair_remove_label_ratio / iteration_graph_times;

        /*dynamic maintenance*/
        initialize_global_values_dynamic_hop_constrained(V + group_num, mm.thread_num, mm.upper_k);
        graph_change_and_label_maintenance(instance_graph, mm, V, weightIncrease_time, weightDecrease_time, weightChange_ratio, avg_maintain_time, batch_size, insert_and_delete_edge, ec_min, ec_max);

        /*debug*/
        if (0)
        {
            instance_graph.print();
            mm.print_L();
        }

        hop_constrained_check_correctness(mm, instance_graph, iteration_source_times, iteration_terminal_times, mm.upper_k);

        long long int index_size = 0;
        for (auto it = mm.L.begin(); it != mm.L.end(); it++)
        {
            index_size = index_size + (*it).size();
        }
        avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;

        mm.clear_labels();
    }

    cout << "avg_index_time: " << avg_index_time << "s" << endl;
    cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
    cout << "avg_maintain_time: " << (double)avg_maintain_time / iteration_graph_times << "s" << endl;
    cout << "avg_canonical_repair_remove_label_ratio: " << avg_canonical_repair_remove_label_ratio << endl;
}
