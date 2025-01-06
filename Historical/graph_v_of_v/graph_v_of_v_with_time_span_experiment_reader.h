#pragma once
#include <vector>
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span_non_hop_constrained.h"
#include "Historical/experiment/mark_timer.h"

#include <filesystem>
#include <queue>
#include <fstream>

void PLL_clear_global_values()
{
    this_parallel_PLL_is_running_595 = false;
    vector<vector<two_hop_label>>().swap(L_temp_595);
    PPR_type().swap(PPR_595);
    queue<int>().swap(Qid_595);
    vector<vector<int>>().swap(P_dij_595);
    vector<vector<int>>().swap(T_dij_595);
    vector<vector<PLL_handle_t_for_sp>>().swap(Q_handles_595);
}
void initialize_experiment_global_values_dynamic(int N, int thread_num)
{
    Dis.resize(thread_num);
    Q_value.resize(thread_num);
    Q_handles.resize(thread_num);
    queue<int>().swap(Qid_595);
    for (int i = 0; i < thread_num; i++)
    {
        Dis[i].resize(N, {-1, -1});
        Q_value[i].resize(N, MAX_VALUE);
        Q_handles[i].resize(N);
        Qid_595.push(i);
    }
}

namespace fs = std::filesystem;
boost::random::mt19937 boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};
struct change_edge_info
{
    int v1;
    int v2;
    int weight;
    int time;
};

class experiment_config
{

private:
    const fs::path experiment_path;
    std::ofstream outputFile;
    const fs::path save_path;
    const int iteration;
    const int change_num;
    const int upper;
    const int lower;
    const bool is_debug;
    int v_num = 0;
    int e_num = 0;
    uniform_int_distribution<> random_v;
    uniform_int_distribution<> random_weight;
    std::queue<change_edge_info> q;

    two_hop_case_info mm;
    mark_timer mm_mark_timer;

    two_hop_case_info mm2021;
    mark_timer mm2021_mark_timer;

    void txt_read_base()
    {
        std::string readPath = this->is_debug ? this->save_path.generic_u8string() : this->experiment_path.generic_u8string();
        std::string line_content;
        graph_v_of_v<int> instance_graph;
        // 读取文件
        std::ifstream myfile(readPath);
        if (myfile.is_open())
        {
            while (getline(myfile, line_content))
            {
                if (this->is_debug)
                {
                    std::vector<std::string> Parsed_content = parse_string(line_content, " ");
                    int v1 = std::stoi(Parsed_content[0]);
                    int v2 = std::stoi(Parsed_content[1]);
                    int w = std::stoi(Parsed_content[2]);
                    int time = std::stoi(Parsed_content[3]);
                    if (time > 0)
                    {
                        this->q.push({v1, v2, w, time});
                    }
                    else
                    {
                        instance_graph.add_edge(v1, v2, w);
                    }
                }
                else
                {
                    std::vector<std::string> Parsed_content = parse_string(line_content, "\t");

                    if (!Parsed_content[0].compare("#"))
                    {
                        if (!Parsed_content[1].compare("Nodes"))
                        {
                            instance_graph.resize(std::stoi(Parsed_content[2]));
                            this->v_num = std::stoi(Parsed_content[2]);
                            this->random_v = uniform_int_distribution<>(0, this->v_num);
                        }
                        else if (!Parsed_content[1].compare("Edges"))
                        {
                            this->e_num = std::stoi(Parsed_content[2]);
                            // weight_type weight_upper_limit, weight_type weight_lower_limit
                        }
                    }
                    else
                    {
                        int v1 = std::stoi(Parsed_content[0]);
                        int v2 = std::stoi(Parsed_content[1]);
                        int w = this->random_weight(boost_random_time_seed);
                        instance_graph.add_edge(v1, v2, w);
                        this->txt_save(v1, v2, w, 0);
                    }
                }
            }

            this->graph_with_time_span = graph_v_of_v_with_time_span<int>(this->v_num, this->e_num, this->upper, this->lower);
            this->graph_with_time_span.add_graph_time(instance_graph, 0);
            this->graphs.push_back(instance_graph);
            if (!this->is_debug)
            {
                // 迭代指定次数 生成随机改变的边的数组 并保存到硬盘
                for (int i = 1; i < this->iteration; i++)
                {
                    std::cout << i << std::endl;
                    for (int j = 0; j < this->change_num; j++)
                    {
                        int index_i = this->random_v(boost_random_time_seed);
                        uniform_int_distribution<> dis_inner(0, this->graphs[0].ADJs[index_i].size() - 1);
                        int index_j = dis_inner(boost_random_time_seed);
                        int i_j_weight = this->random_weight(boost_random_time_seed);
                        q.push({index_i, index_j, i_j_weight, i});
                        // 持久化
                        txt_save(index_i, index_j, i_j_weight, i);
                    }
                }
            }
            myfile.close(); // close the file
        }
        else
        {
            std::cout << "Unable to open file " << readPath << std::endl
                      << "Please check the file location or file name." << std::endl; // throw an error message
            getchar();                                                                // keep the console window
            exit(1);                                                                  // end the program
        }
    }
    void txt_save(int v1, int v2, int w, int time)
    {
        this->outputFile << v1 << " " << v2 << " " << w << " " << time << "\n";
    };
    void txt_close()
    {
        this->outputFile << "EOF" << std::endl;
        this->outputFile.close();
    }

public:
    graph_v_of_v_with_time_span<int> graph_with_time_span;
    vector<graph_v_of_v<int>> graphs;
    experiment_config(fs::path _experiment_path, fs::path _save_path, int _iteration, int _change_num, bool _is_debug, int _upper, int _lower) : experiment_path(_experiment_path), save_path(_save_path), iteration(_iteration), change_num(_change_num), is_debug(_is_debug), upper(_upper), lower(_lower)
    {
        random_weight = uniform_int_distribution<>(lower, upper);
        this->outputFile.precision(10);
        this->outputFile.setf(std::ios::fixed);
        this->outputFile.setf(std::ios::showpoint);
        this->outputFile.open(_save_path.generic_u8string());
        this->outputFile << "time" << "\t" << this->iteration << std::endl;
    }
    int init()
    {
        mm.max_labal_byte_size = 6e9;
        mm.max_run_time_seconds = 1e4;
        mm.use_2M_prune = 1;
        mm.use_rank_prune = 1;
        mm.use_canonical_repair = 1;
        mm.thread_num = 1;
        mm2021.max_labal_byte_size = 6e9;
        mm2021.max_run_time_seconds = 1e4;
        mm2021.use_2M_prune = 1;
        mm2021.use_rank_prune = 1;
        mm2021.use_canonical_repair = 1;
        mm2021.thread_num = 5;
        return 0;
    }
    int process()
    {
        // 读取图的数据到指定的对象中 graphs为原始图的列表 graph_with_time_span为时序图的对象
        this->txt_read_base();
        // 初始化nonhop
        initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
        mm_mark_timer.mark();
        PLL(this->graphs[0], this->mm);
        mm_mark_timer.add();
        mm_mark_timer.push();
        initialize_experiment_global_values_dynamic(this->v_num, this->mm.thread_num);
        mm2021_mark_timer.mark();
        PLL(this->graphs[0], this->mm2021);
        mm2021_mark_timer.add();
        mm2021_mark_timer.push();
        // 动态维护
        int time = 0;
        while (!this->q.empty())
        {
            change_edge_info edge_info = this->q.front();
            q.pop();
            //TODO 把维护代码集成进来
        }

        return 0;
    }
    int close()
    {
        this->txt_close();
        return 0;
    }
};
