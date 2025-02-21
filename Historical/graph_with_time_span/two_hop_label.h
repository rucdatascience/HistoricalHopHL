#pragma once
#include <climits>
#include <limits>
#include <vector>
#include <iostream>
#include "Historical/utils/BinaryPersistence.h"
#include <algorithm>
#include "Historical/utils/vector_operations.h"
using WEIGHT_TYPE = long;
namespace experiment {
	namespace PPR_TYPE {
		using PPR_type = std::vector<std::vector<std::pair<int, std::vector<int>>>>;
		template <typename T>
		void saveBinary(std::ofstream& out, const PPR_type& vec) {
			size_t size = vec.size();
			saveBinary(out, size);
			for (const auto& item : vec) {
				saveBinary(out, item);
			}
		}

		template <typename T>
		void loadBinary(std::ifstream& in, PPR_type& vec) {
			size_t size;
			loadBinary(in, size);
			std::vector<T>().swap(vec);
			vec.resize(size);
			for (auto& item : vec) {
				loadBinary(in, item);
			}
		}

		int PPR_binary_operations_insert(std::vector<int>& input_vector, int key)
		{

			int left = 0, right = input_vector.size() - 1;

			while (left <= right) // it will be skept when input_vector.size() == 0
			{
				int mid = left + ((right - left) / 2); // mid is between left and right (may be equal);
				if (input_vector[mid] == key)
				{
					return mid;
				}
				else if (input_vector[mid] > key)
				{
					right = mid - 1; // the elements after right are always either empty, or have larger keys than input key
				}
				else
				{
					left = mid + 1; // the elements before left are always either empty, or have smaller keys than input key
				}
			}

			/*the following code is used when key is not in vector, i.e., left > right, specifically, left = right + 1;
			the elements before left are always either empty, or have smaller keys than input key;
			the elements after right are always either empty, or have larger keys than input key;
			so, the input key should be insert between right and left at this moment*/
			input_vector.insert(input_vector.begin() + left, key);
			return left;
		}

		void PPR_insert(PPR_type& PPR, int v1, int v2, int v3)
		{

			/*add v3 into PPR(v1, v2)*/

			int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
			if (pos == -1)
			{
				std::vector<int> x = { v3 };
				graph_hash_of_mixed_weighted_binary_operations_insert(PPR[v1], v2, x);
			}
			else
			{
				PPR_binary_operations_insert(PPR[v1][pos].second, v3);
			}
		}

		std::vector<int> PPR_retrieve(PPR_type& PPR, int v1, int v2)
		{

			/*retrieve PPR(v1, v2)*/

			int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
			if (pos == -1)
			{
				std::vector<int> x;
				return x;
			}
			else
			{
				return PPR[v1][pos].second;
			}
		}

		void PPR_replace(PPR_type& PPR, int v1, int v2, std::vector<int>& loads)
		{

			/*replace PPR(v1, v2) = loads*/

			int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
			if (pos == -1)
			{
				graph_hash_of_mixed_weighted_binary_operations_insert(PPR[v1], v2, loads);
			}
			else
			{
				PPR[v1][pos].second = loads;
			}
		}

		void PPR_erase(PPR_type& PPR, int v1, int v2, int v3)
		{
			int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
			for (auto it = PPR[v1][pos].second.begin(); it != PPR[v1][pos].second.end(); it++)
			{
				if (*it == v3)
				{
					PPR[v1][pos].second.erase(it);
					break;
				}
			}
		}

	}

	namespace nonhop {
		class two_hop_label
		{
		public:
			int vertex;
			WEIGHT_TYPE distance;
			int t_s, t_e;
			two_hop_label() {
				t_s = 0;
				t_e = INT_MAX;
				vertex = std::numeric_limits<int>::max();
				distance = std::numeric_limits<WEIGHT_TYPE>::max();
			}
			two_hop_label(int start_time)
			{
				t_s = start_time;
				t_e = INT_MAX;
				vertex = std::numeric_limits<int>::max();
				distance = std::numeric_limits<WEIGHT_TYPE>::max();
			}
			void serialize(std::ofstream& out) const {
				experiment::saveBinary(out, vertex);
				experiment::saveBinary(out, distance);
				experiment::saveBinary(out, t_s);
				experiment::saveBinary(out, t_e);
			}

			void deserialize(std::ifstream& in) {
				experiment::loadBinary(in, vertex);
				experiment::loadBinary(in, distance);
				experiment::loadBinary(in, t_s);
				experiment::loadBinary(in, t_e);
			}
		};
		bool operator<(two_hop_label const& x, two_hop_label const& y)
		{
			return x.distance > y.distance; // < is the max-heap; > is the min heap
		}
		bool compare_two_hop_label_small_to_large(two_hop_label& i, two_hop_label& j)
		{
			if (i.t_e != j.t_e)
				return i.t_e > j.t_e;	// t_e降序
			return i.vertex < j.vertex; // < is from small to big; > is from big to small
		};
		class two_hop_case_info {
		public:
			int thread_num = 1;

			/*labels*/
			std::vector<std::vector<two_hop_label>> L;
			PPR_TYPE::PPR_type PPR;

			void serialize(std::ofstream& out) const {
				experiment::saveBinary(out, thread_num);
				experiment::saveBinary(out, L);
				experiment::saveBinary(out, PPR);
			}

			void deserialize(std::ifstream& in) {
				experiment::loadBinary(in, thread_num);
				experiment::loadBinary(in, L);
				experiment::loadBinary(in, PPR);
			}


			/*clear labels*/
			void clear_labels()
			{
				std::vector<std::vector<two_hop_label>>().swap(L);
				PPR_TYPE::PPR_type().swap(PPR);
			}

			/*compute label size; this should equal label_size_after_canonical_repair when use_canonical_repair==true*/
			long long int compute_L_byte_size()
			{
				long long int size = 0;
				for (auto it = L.begin(); it != L.end(); it++)
				{
					size = size + (*it).size() * sizeof(two_hop_label); // 12 byte per two_hop_label
				}
				return size;
			}

			long long int compute_PPR_byte_size()
			{
				long long int size = 0;
				for (int i = 0; i < PPR.size(); i++)
				{
					for (int j = 0; j < PPR[i].size(); j++)
					{
						size = size + (PPR[i][j].second.size() + 1) * sizeof(int); // + 1 ��Ӧ PPR[i][j].first
					}
				}
				return size;
			}

			/*printing*/
			void print_L()
			{
				std::cout << "print_L:" << std::endl;
				for (int i = 0; i < L.size(); i++)
				{
					std::cout << "L[" << i << "]=";
					for (int j = 0; j < L[i].size(); j++)
					{
						std::cout << "{" << L[i][j].vertex << "," << L[i][j].distance << "," << L[i][j].t_s << "," << L[i][j].t_e << "}";
					}
					std::cout << std::endl;
				}
			}
			void print_PPR()
			{
				std::cout << "print_PPR:" << std::endl;
				for (int i = 0; i < PPR.size(); i++)
				{
					for (int j = 0; j < PPR[i].size(); j++)
					{
						std::cout << "PPR(" << i << "," << PPR[i][j].first << "): ";
						for (int k = 0; k < PPR[i][j].second.size(); k++)
						{
							std::cout << PPR[i][j].second[k] << " ";
						}
						std::cout << std::endl;
					}
				}
			}

			/*record_all_details*/
			void record_all_details(std::string save_name)
			{
				std::ofstream outputFile;
				outputFile.precision(6);
				outputFile.setf(std::ios::fixed);
				outputFile.setf(std::ios::showpoint);
				outputFile.open(save_name + ".txt");

				outputFile << "PLL info:" << std::endl;
				outputFile << "thread_num=" << thread_num << std::endl;
				outputFile << "compute_label_byte_size()=" << compute_L_byte_size() << std::endl;

				outputFile.close();
			}

			long long int query(int source, int terminal, int t_s, int t_e)
			{
				if (source == terminal)
				{
					return 0;
				}

				int distance = std::numeric_limits<int>::max();
				auto vector1_check_pointer = L[source].begin();
				auto vector2_check_pointer = L[terminal].begin();
				auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();

				for (auto vector1_begin = vector1_check_pointer; vector1_begin != pointer_L_s_end; vector1_begin++)
				{
					// cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
					for (auto vector2_begin = vector2_check_pointer; vector2_begin != pointer_L_t_end; vector2_begin++)
					{
						// cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
						if (vector1_begin->vertex == vector2_begin->vertex && std::max(vector1_begin->t_s, std::max(vector2_begin->t_s, t_s)) <= std::min(vector1_begin->t_e, std::min(vector2_begin->t_e, t_e)))
						{
							long long int dis = (long long int)vector1_begin->distance + vector2_begin->distance;
							if (distance > dis)
							{
								distance = dis;
								// cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance <<  ") " << endl;
								// cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance <<  ") " << endl;
							}
						}
					}
				}
				return distance;
			}
		};
	};
	namespace hop {
		class two_hop_label
		{
		public:
			int hub_vertex, hop;
			WEIGHT_TYPE distance;
			int t_s, t_e;
			// hop_constrained_two_hop_label() {}
			// hop_constrained_two_hop_label(int _vertex, int _hop, int _dis)
			// {
			//     hub_vertex = _vertex;
			//     hop = _hop;
			//     distance = _dis;
			// }
			two_hop_label()
				: t_s(0), t_e(std::numeric_limits<int>::max()) {
				hub_vertex = 0;
				hop = 0;
				distance = std::numeric_limits<WEIGHT_TYPE>::max();
			}
			void serialize(std::ofstream& out) const {
				experiment::saveBinary(out, hub_vertex);
				experiment::saveBinary(out, hop);
				experiment::saveBinary(out, distance);
				experiment::saveBinary(out, t_s);
				experiment::saveBinary(out, t_e);
			}

			void deserialize(std::ifstream& in) {
				experiment::loadBinary(in, hub_vertex);
				experiment::loadBinary(in, hop);
				experiment::loadBinary(in, distance);
				experiment::loadBinary(in, t_s);
				experiment::loadBinary(in, t_e);
			}
		};


		class two_hop_case_info {
		public:
			/*hop bounded*/
			int thread_num = 1;
			int upper_k = 0;

			/*labels*/
			std::vector<std::vector<two_hop_label>> L;
			PPR_TYPE::PPR_type PPR;

			void serialize(std::ofstream& out) const {
				experiment::saveBinary(out, thread_num);
				experiment::saveBinary(out, upper_k);
				experiment::saveBinary(out, L);
				experiment::saveBinary(out, PPR);
			}

			void deserialize(std::ifstream& in) {
				experiment::loadBinary(in, thread_num);
				experiment::loadBinary(in, upper_k);
				experiment::loadBinary(in, L);
				experiment::loadBinary(in, PPR);
			}

			long long int compute_label_bit_size()
			{
				long long int size = 0;
				for (auto& xx : L)
				{
					size = size + xx.size() * sizeof(two_hop_label);
				}
				return size;
			}

			/*clear labels*/
			void clear_labels()
			{
				std::vector<std::vector<two_hop_label>>().swap(L);
				PPR_TYPE::PPR_type().swap(PPR);
			}

			void print_L()
			{
				int index = 0;
				std::cout << "print_L: (hub_vertex, hop, distance)" << std::endl;
				for (auto& xx : L)
				{
					std::cout << "vertex " << index++ << ": ";
					for (auto& yy : xx)
					{
						std::cout << "(" << yy.hub_vertex << "," << yy.hop << "," << yy.distance << "," << yy.t_s << "," << yy.t_e << ")";
					}
					std::cout << std::endl;
				}
			}

			void print_PPR()
			{
				std::cout << "print_PPR:" << std::endl;
				for (int i = 0; i < PPR.size(); i++)
				{
					for (int j = 0; j < PPR[i].size(); j++)
					{
						std::cout << "PPR(" << i << "," << PPR[i][j].first << "): ";
						for (int k = 0; k < PPR[i][j].second.size(); k++)
						{
							std::cout << PPR[i][j].second[k] << " ";
						}
						std::cout << std::endl;
					}
				}
			}

			void print_L_vk(int v_k)
			{
				for (auto it = L[v_k].begin(); it != L[v_k].end(); it++)
				{
					std::cout << "(" << it->hub_vertex << "," << it->hop << "," << it->distance << "," << it->t_s << "," << it->t_e << ")";
				}
				std::cout << std::endl;
			}

			/*record_all_details*/
			void record_all_details(std::string save_name)
			{
				std::ofstream outputFile;
				outputFile.precision(6);
				outputFile.setf(std::ios::fixed);
				outputFile.setf(std::ios::showpoint);
				outputFile.open(save_name + ".txt");

				outputFile << "hop_constrained_case_info:" << std::endl;
				outputFile << "thread_num=" << thread_num << std::endl;
				outputFile << "upper_k=" << upper_k << std::endl;

				outputFile << "compute_label_bit_size()=" << compute_label_bit_size() << std::endl;

				outputFile.close();
			}

			long long int query(int source, int terminal, int t_s, int t_e, int hop_cst)
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

				for (auto vector1_begin = vector1_check_pointer; vector1_begin != pointer_L_s_end; vector1_begin++)
				{
					// cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
					for (auto vector2_begin = vector2_check_pointer; vector2_begin != pointer_L_t_end; vector2_begin++)
					{
						// cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
						if (vector1_begin->hub_vertex == vector2_begin->hub_vertex && vector1_begin->hop + vector2_begin->hop <= hop_cst && std::max(vector1_begin->t_s, std::max(vector2_begin->t_s, t_s)) <= std::min(vector1_begin->t_e, std::min(vector2_begin->t_e, t_e)))
						{
							long long int dis = (long long int)vector1_begin->distance + vector2_begin->distance;
							if (distance > dis)
							{
								distance = dis;
								// cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance <<  ") " << endl;
								// cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance <<  ") " << endl;
							}
						}
					}
				}
				return distance;
			}
		};
	}
	void saveBinary(std::ofstream& out, const nonhop::two_hop_label& data) {
		data.serialize(out);
	}

	void loadBinary(std::ifstream& in, nonhop::two_hop_label& data) {
		data.deserialize(in);
	}

	void saveBinary(std::ofstream& out, const hop::two_hop_label& data) {
		data.serialize(out);
	}

	void loadBinary(std::ifstream& in, hop::two_hop_label& data) {
		data.deserialize(in);
	}

	void saveBinary(std::ofstream& out, const hop::two_hop_case_info& data) {
		data.serialize(out);
	}

	void loadBinary(std::ifstream& in, hop::two_hop_case_info& data) {
		data.deserialize(in);
	}

	void saveBinary(std::ofstream& out, const nonhop::two_hop_case_info& data) {
		data.serialize(out);
	}

	void loadBinary(std::ifstream& in, nonhop::two_hop_case_info& data) {
		data.deserialize(in);
	}

}