#pragma once
#include <vector>
#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span.h"

#include <filesystem>
namespace fs = std::filesystem;

class experiment_config
{

private:
    const fs::path experiment_path;
    const int iteration;
    const int change_num;
    const int upper;
    const int lower;
    const bool is_debug;
    int v_num = 0;
    int e_num = 0;

    void txt_read_base(std::string save_name)
    {
        std::string line_content;
        graph_v_of_v<int> instance_graph;
        // 读取文件
        std::ifstream myfile(save_name);
        if (myfile.is_open())
        {
            while (getline(myfile, line_content))
            {
                std::vector<std::string> Parsed_content = parse_string(line_content, "\t");

                if (!Parsed_content[0].compare("#"))
                {
                    if (!Parsed_content[1].compare("Nodes"))
                    {
                        instance_graph.resize(std::stoi(Parsed_content[2]));
                        this->v_num = std::stoi(Parsed_content[2]);
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
                    instance_graph.add_edge(v1, v2, 1);
                }
            }
            this->graph_with_time_span = graph_v_of_v_with_time_span<int>(this->v_num, this->e_num, this->upper, this->lower);

            this->graph_with_time_span.add_graph_time(instance_graph, 0);

            this->graphs.push_back(instance_graph);
            myfile.close(); // close the file
        }
        else
        {
            std::cout << "Unable to open file " << save_name << std::endl
                      << "Please check the file location or file name." << std::endl; // throw an error message
            getchar();                                                                // keep the console window
            exit(1);                                                                  // end the program
        }
    }

public:
    graph_v_of_v_with_time_span<int> graph_with_time_span;
    vector<graph_v_of_v<int>> graphs;
    experiment_config(fs::path _experiment_path, int _iteration, int _change_num, bool _is_debug, int _upper, int _lower) : experiment_path(_experiment_path), iteration(_iteration), change_num(_change_num), is_debug(_is_debug), upper(_upper), lower(_lower)
    {
    }

    int process()
    {
        // 迭代指定次数 生成随机改变的边的数组 并保存到硬盘
        // 如果debug模式的话 读取上一次的结果
        std::string name = this->experiment_path.u8string();
        this->txt_read_base(name);

        // 初始化索引

        // 动态维护
    }
};
