#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include "Historical/experiment/experiment_config.h"
#include <CPU/tool_functions/sorted_vector_binary_operations.h>
#include <CPU/text_mining/binary_save_read_vector_of_vectors.h>
#include "Historical/utils/BinaryPersistence.h"

namespace experiment
{
	bool sortEdgeById(const std::pair<int, int> &i, std::pair<int, int> &j)
	{
		/*< is nearly 10 times slower than >*/
		return i.first < j.first; // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
	}

	bool compare_graph_v_of_v_update_vertexIDs_by_degrees_large_to_small(const std::pair<int, int> &i, std::pair<int, int> &j)
	{
		/*< is nearly 10 times slower than >*/
		return i.second > j.second; // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
	}

	template <typename weight_type> // weight_type may be int, long long int, float, double...
	class graph
	{
	public:
		/*
		this class only suits ideal vertex IDs: from 0 to V-1;

		this class is for undirected and edge-weighted graph
		*/
		std::vector<std::vector<std::pair<int, weight_type>>> ADJs;

		/*constructors*/
		graph() {}
		graph(int n)
		{
			ADJs.resize(n); // initialize n vertices
		}
		int size()
		{
			return ADJs.size();
		}

		void resize(int n)
		{
			ADJs.resize(n); // initialize n vertices
		}

		std::vector<std::pair<int, weight_type>> &operator[](int i)
		{
			return ADJs[i];
		}

		long long int computeSize() const
		{
			long long int res = 0;
			for (const auto &item_first : this->ADJs)
			{
				for (const auto &item_second : item_first)
				{
					res += sizeof(int);
					res += sizeof(weight_type);
				}
			}
			return res;
		}

		/*class member functions*/
		void add_edge(int e1, int e2, weight_type ec)
		{
			/*we assume that the size of g is larger than e1 or e2;
			 this function can update edge weight; there will be no redundent edge*/

			/*
			Add the edges (e1,e2) and (e2,e1) with the weight ec
			When the edge exists, it will update its weight.
			Time complexity:
				O(log n) When edge already exists in graph
				O(n) When edge doesn't exist in graph
			*/
			sorted_vector_binary_operations_insert(ADJs[e1], e2, ec);
			sorted_vector_binary_operations_insert(ADJs[e2], e1, ec);
		}
		void remove_edge(int e1, int e2)
		{

			/*we assume that the size of g is larger than e1 or e2*/
			/*
			 Remove the edges (e1,e2) and (e2,e1)
			 If the edge does not exist, it will do nothing.
			 Time complexity: O(n)
			*/

			sorted_vector_binary_operations_erase(ADJs[e1], e2);
			sorted_vector_binary_operations_erase(ADJs[e2], e1);
		}

		void remove_all_adjacent_edges(int v)
		{

			for (auto it = ADJs[v].begin(); it != ADJs[v].end(); it++)
			{
				sorted_vector_binary_operations_erase(ADJs[it->first], v);
			}

			std::vector<std::pair<int, weight_type>>().swap(ADJs[v]);
		}

		bool contain_edge(int e1, int e2) const
		{

			/*
			Return true if graph contain edge (e1,e2)
			Time complexity: O(logn)
			*/

			return sorted_vector_binary_operations_search(ADJs[e1], e2);
		}

		weight_type edge_weight(int e1, int e2) const
		{

			/*
			Return the weight of edge (e1,e2)
			If the edge does not exist, return std::numeric_limits<double>::max()
			Time complexity: O(logn)
			*/

			return sorted_vector_binary_operations_search_weight(ADJs[e1], e2);
		}

		long long int edge_number() const
		{

			/*
			Returns the number of edges in the figure
			(e1,e2) and (e2,e1) will be counted only once
			Time complexity: O(n)
			*/

			int num = 0;
			for (const auto &it : ADJs)
			{
				num = num + it.size();
			}

			return num / 2;
		}
		void print() const
		{
			std::cout << "graph_print:" << std::endl;
			int size = ADJs.size();
			for (int i = 0; i < size; i++)
			{
				std::cout << "Vertex " << i << " Adj List: ";
				int v_size = ADJs[i].size();
				for (int j = 0; j < v_size; j++)
				{
					std::cout << "<" << ADJs[i][j].first << "," << ADJs[i][j].second << "> ";
				}
				std::cout << std::endl;
			}
			std::cout << "graph_v_of_v_print END" << std::endl;
		}

		void clear()
		{

			return std::vector<std::vector<std::pair<int, weight_type>>>().swap(ADJs);
		}

		int degree(int v) const
		{
			return ADJs[v].size();
		}

		int search_adjv_by_weight(int e1, weight_type ec) const
		{

			for (auto &xx : ADJs[e1])
			{
				if (xx.second == ec)
				{
					return xx.first;
				}
			}

			return -1;
		}

		void txt_save(std::string save_name) const
		{

			std::ofstream outputFile;
			outputFile.precision(10);
			outputFile.setf(std::ios::fixed);
			outputFile.setf(std::ios::showpoint);
			outputFile.open(save_name);

			outputFile << "|V|= " << ADJs.size() << std::endl;
			outputFile << "|E|= " << graph<weight_type>::edge_number() << std::endl;
			outputFile << std::endl;

			int size = ADJs.size();
			for (int i = 0; i < size; i++)
			{
				int v_size = ADJs[i].size();
				for (int j = 0; j < v_size; j++)
				{
					if (i < ADJs[i][j].first)
					{
						outputFile << "Edge " << i << " " << ADJs[i][j].first << " " << ADJs[i][j].second << '\n';
					}
				}
			}
			outputFile << std::endl;

			outputFile << "EOF" << std::endl;
		}

		void txt_read(std::string save_name)
		{

			graph<weight_type>::clear();

			std::string line_content;
			std::ifstream myfile(save_name); // open the file
			if (myfile.is_open())			 // if the file is opened successfully
			{
				while (getline(myfile, line_content)) // read file line by line
				{
					std::vector<std::string> Parsed_content = experiment::parse_string(line_content, " ");

					if (!Parsed_content[0].compare("|V|=")) // when it's equal, compare returns 0
					{
						ADJs.resize(std::stoi(Parsed_content[1]));
					}
					else if (!Parsed_content[0].compare("Edge"))
					{
						int v1 = std::stoi(Parsed_content[1]);
						int v2 = std::stoi(Parsed_content[2]);
						weight_type ec = std::stod(Parsed_content[3]);
						graph<weight_type>::add_edge(v1, v2, ec);
					}
				}

				myfile.close(); // close the file
			}
			else
			{
				std::cout << "Unable to open file " << save_name << std::endl
						  << "Please check the file location or file name." << std::endl; // throw an error message
				getchar();																  // keep the console window
				exit(1);																  // end the program
			}
		}

		void graph_v_of_v_update_vertexIDs_by_degrees_large_to_small()
		{
			int N = this->ADJs.size();

			std::vector<std::pair<int, int>> sorted_vertices;
			for (int i = 0; i < N; i++)
			{
				sorted_vertices.push_back({i, this->ADJs[i].size()});
			}
			std::sort(sorted_vertices.begin(), sorted_vertices.end(), compare_graph_v_of_v_update_vertexIDs_by_degrees_large_to_small);
			std::vector<int> vertexID_old_to_new(N);
			for (int i = 0; i < N; i++)
			{
				vertexID_old_to_new[sorted_vertices[i].first] = i;
			}
			for (int i = 0; i < N; i++)
			{
				std::vector<std::pair<int, weight_type>> &edge_info = this->ADJs.at(i);
				for (std::pair<int, weight_type> &edge : edge_info)
				{
					edge.first = vertexID_old_to_new[edge.first];
				}
				std::sort(edge_info.begin(), edge_info.end(), sortEdgeById);
			}
			for (int i = 0; i < N; i++)
			{
				if (vertexID_old_to_new[i] < i)
				{
					std::swap(this->ADJs[i], this->ADJs[vertexID_old_to_new[i]]);
				}
			}
		}

		void serialize(std::ofstream &out) const
		{
			saveBinary(out, this->ADJs);
		}

		void deserialize(std::ifstream &in)
		{
			loadBinary(in, this->ADJs);
		}
	};

	template <typename weight_type>
	class BinarySerializer<graph<weight_type>>
	{
	public:
		static void saveBinary(std::ofstream &out, const graph<weight_type> &vec)
		{
			vec.serialize(out);
		}

		static void loadBinary(std::ifstream &in, graph<weight_type> &vec)
		{
			vec.deserialize(in);
		}
	};
}