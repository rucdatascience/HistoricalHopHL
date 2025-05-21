#pragma once
#include <climits>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include "Historical/utils/BinaryPersistence.h"
#include <algorithm>
#include "Historical/utils/vector_operations.h"
using WEIGHT_TYPE = long;
#define weightTYPE int
#define MAX_VALUE 1e7
namespace experiment
{
	static long long PPR_INSERT_RUC = 0;
	static long long PPR_INSERT_2021 = 0;
	namespace PPR_TYPE
	{
		using PPR_type = std::vector<std::vector<std::pair<int, std::vector<int>>>>;

		long long int getSize(PPR_type &PPR)
		{
			long long int res = 0;
			for (const auto &array : PPR)
			{
				for (const auto &arrayInner : array)
				{
					res += arrayInner.second.size();
				}
			}
			return res;
		};

		int PPR_binary_operations_insert(std::vector<int> &input_vector, int key)
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

		void PPR_insert(PPR_type &PPR, int v1, int v2, int v3)
		{

			/*add v3 into PPR(v1, v2)*/

			int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
			if (pos == -1)
			{
				std::vector<int> x = {v3};
				graph_hash_of_mixed_weighted_binary_operations_insert(PPR[v1], v2, x);
			}
			else
			{
				PPR_binary_operations_insert(PPR[v1][pos].second, v3);
			}
		}
		void PPR_insert_mark(PPR_type &PPR, int v1, int v2, int v3,bool isRuc)
		{
			if(isRuc){
				PPR_INSERT_RUC++;
			}else{
				PPR_INSERT_2021++;
			}
			/*add v3 into PPR(v1, v2)*/

			int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(PPR[v1], v2);
			if (pos == -1)
			{
				std::vector<int> x = {v3};
				graph_hash_of_mixed_weighted_binary_operations_insert(PPR[v1], v2, x);
			}
			else
			{
				PPR_binary_operations_insert(PPR[v1][pos].second, v3);
			}
		}

		std::vector<int> PPR_retrieve(PPR_type &PPR, int v1, int v2)
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

		void PPR_replace(PPR_type &PPR, int v1, int v2, std::vector<int> &loads)
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

		void PPR_erase(PPR_type &PPR, int v1, int v2, int v3)
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

	namespace nonhop
	{
		static long long ruc_query_count = 0;
		static long long a2021_query_count = 0;
		static long long ruc_label_query_count = 0;
		static long long a2021_label_query_count = 0;
		static long long ruc_label_insert_count = 0;
		static long long a2021_label_insert_count = 0;
		class two_hop_label
		{
		public:
			int vertex;
			WEIGHT_TYPE distance;
			int t_s, t_e;
			bool operator==(const two_hop_label &other) const
			{
				return vertex == other.vertex && distance == other.distance && t_s == other.t_s && t_e == other.t_e;
			}
			two_hop_label()
			{
				t_s = 0;
				t_e = INT_MAX;
				vertex = std::numeric_limits<int>::max();
				distance = std::numeric_limits<WEIGHT_TYPE>::max();
			}

			two_hop_label(const two_hop_label &other)
			{
				t_s = other.t_s;
				t_e = other.t_e;
				vertex = other.vertex;
				distance = other.distance;
			}

			two_hop_label(int start_time)
			{
				t_s = start_time;
				t_e = INT_MAX;
				vertex = std::numeric_limits<int>::max();
				distance = std::numeric_limits<WEIGHT_TYPE>::max();
			}
			void serialize(std::ofstream &out) const
			{
				experiment::saveBinary(out, vertex);
				experiment::saveBinary(out, distance);
				experiment::saveBinary(out, t_s);
				experiment::saveBinary(out, t_e);
			}

			void deserialize(std::ifstream &in)
			{
				experiment::loadBinary(in, vertex);
				experiment::loadBinary(in, distance);
				experiment::loadBinary(in, t_s);
				experiment::loadBinary(in, t_e);
			}
		};
		bool operator<(two_hop_label const &x, two_hop_label const &y)
		{
			return x.distance > y.distance; // < is the max-heap; > is the min heap
		}
		bool compare_two_hop_label_small_to_large(two_hop_label &i, two_hop_label &j)
		{
			if (i.t_e != j.t_e)
				return i.t_e > j.t_e;	// t_e降序
			return i.vertex < j.vertex; // < is from small to big; > is from big to small
		};

		// method of 2hop label
		std::pair<int, int> graph_weighted_two_hop_extract_distance_and_hub_in_current(std::vector<std::vector<two_hop_label>> &L, int source, int terminal)
		{
			/*return std::numeric_limits<double>::max() is not connected*/

			if (source == terminal)
			{
				return {0, source};
			}

			long long int distance = std::numeric_limits<long long int>::max(); // if disconnected, return this large value
			int common_hub;

			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					long long int dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						common_hub = vector1_check_pointer->vertex;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
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

		std::pair<int, int> graph_weighted_two_hop_extract_distance_and_hub_in_current_mark(std::vector<std::vector<two_hop_label>> &L, int source, int terminal,bool isRuc)
		{
			/*return std::numeric_limits<double>::max() is not connected*/
			if(isRuc){
				ruc_query_count++;
			}else{
				a2021_query_count++;
			}
			if (source == terminal)
			{
				return {0, source};
			}

			long long int distance = std::numeric_limits<long long int>::max(); // if disconnected, return this large value
			int common_hub;

			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					long long int dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						common_hub = vector1_check_pointer->vertex;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
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

		std::pair<int, int> graph_weighted_two_hop_extract_distance_and_hub_by_backup_label(std::vector<two_hop_label> &L_s, std::vector<two_hop_label> &L_t)
		{

			/*return std::numeric_limits<double>::max() is not connected*/
			int distance = std::numeric_limits<int>::max(); // if disconnected, return this large value
			int common_hub = -1;

			auto vector1_check_pointer = L_s.begin();
			auto vector2_check_pointer = L_t.begin();
			auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					int dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						common_hub = vector1_check_pointer->vertex;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
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
		std::pair<int, int> graph_weighted_two_hop_extract_distance_and_hub_by_backup_label_mark(std::vector<two_hop_label> &L_s, std::vector<two_hop_label> &L_t,bool isRuc)
		{
			if(isRuc){
				ruc_query_count++;
			}else{
				a2021_query_count++;
			}
			/*return std::numeric_limits<double>::max() is not connected*/
			int distance = std::numeric_limits<int>::max(); // if disconnected, return this large value
			int common_hub = -1;

			auto vector1_check_pointer = L_s.begin();
			auto vector2_check_pointer = L_t.begin();
			auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					int dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						common_hub = vector1_check_pointer->vertex;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
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

		std::pair<two_hop_label, two_hop_label> graph_weighted_two_hop_extract_2hop_label_by_backup_label(std::vector<two_hop_label> &L_s, std::vector<two_hop_label> &L_t)
		{

			/*return std::numeric_limits<double>::max() is not connected*/

			int distance = std::numeric_limits<int>::max(); // if disconnected, return this large value
			int common_hub;
			two_hop_label res1 = two_hop_label{-1};
			two_hop_label res2 = two_hop_label{-1};
			auto vector1_check_pointer = L_s.begin();
			auto vector2_check_pointer = L_t.begin();
			auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					int dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						common_hub = vector1_check_pointer->vertex;
						res1 = *vector1_check_pointer;
						res2 = *vector2_check_pointer;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
				{
					vector2_check_pointer++;
				}
				else
				{
					vector1_check_pointer++;
				}
			}

			return {res1, res2};
		}

		std::pair<two_hop_label, two_hop_label> graph_weighted_two_hop_extract_2hop_label_by_backup_label_mark(std::vector<two_hop_label> &L_s, std::vector<two_hop_label> &L_t,bool isRuc)
		{
			if(isRuc){
				ruc_query_count++;
			}else{
				a2021_query_count++;
			}
			/*return std::numeric_limits<double>::max() is not connected*/

			int distance = std::numeric_limits<int>::max(); // if disconnected, return this large value
			int common_hub;
			two_hop_label res1 = two_hop_label{-1};
			two_hop_label res2 = two_hop_label{-1};
			auto vector1_check_pointer = L_s.begin();
			auto vector2_check_pointer = L_t.begin();
			auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					int dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						common_hub = vector1_check_pointer->vertex;
						res1 = *vector1_check_pointer;
						res2 = *vector2_check_pointer;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
				{
					vector2_check_pointer++;
				}
				else
				{
					vector1_check_pointer++;
				}
			}

			return {res1, res2};
		}


		std::pair<two_hop_label, two_hop_label> graph_weighted_two_hop_extract_2hop_label_by_backup_label_not_real_time(std::vector<two_hop_label> &L_s, std::vector<two_hop_label> &L_t, int time)
		{

			/*return std::numeric_limits<double>::max() is not connected*/

			int distance = std::numeric_limits<int>::max(); // if disconnected, return this large value
			int common_hub;
			two_hop_label res1 = two_hop_label{-1};
			two_hop_label res2 = two_hop_label{-1};
			auto vector1_check_pointer = L_s.begin();
			auto vector2_check_pointer = L_t.begin();
			auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					if (vector1_check_pointer->t_s == time)
					{
						vector1_check_pointer++;
					}
					else if (vector2_check_pointer->t_s == time)
					{
						vector2_check_pointer++;
					}
					int dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						common_hub = vector1_check_pointer->vertex;
						res1 = *vector1_check_pointer;
						res2 = *vector2_check_pointer;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
				{
					vector2_check_pointer++;
				}
				else
				{
					vector1_check_pointer++;
				}
			}

			return {res1, res2};
		}

		int search_sorted_two_hop_label_weight_in_current(std::vector<two_hop_label> &input_vector, int key)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key)
					{
						return input_vector[mid].distance;
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}

			return std::numeric_limits<int>::max();
		};
		int search_sorted_two_hop_label_weight_in_current_mark(std::vector<two_hop_label> &input_vector, int key,bool isRuc)
		{
			if(isRuc){
				ruc_label_query_count++;
			}else{
				a2021_label_query_count++;
			}
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key)
					{
						return input_vector[mid].distance;
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}

			return std::numeric_limits<int>::max();
		};
		
		std::pair<int, int> search_sorted_two_hop_label_weight_and_hub_in_current(std::vector<two_hop_label> &input_vector, int key)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key)
					{
						return {input_vector[mid].distance, mid};
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}

			return {std::numeric_limits<int>::max(), -1};
		}

		std::pair<int, int> search_sorted_two_hop_label_weight_and_hub_in_current_mark(std::vector<two_hop_label> &input_vector, int key,bool isRuc)
		{
			if(isRuc){
				ruc_label_query_count++;
			}else{
				a2021_label_query_count++;
			}
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key)
					{
						return {input_vector[mid].distance, mid};
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}

			return {std::numeric_limits<int>::max(), -1};
		}


		two_hop_label search_sorted_two_hop_label_in_current(std::vector<two_hop_label> &input_vector, int key)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key)
					{
						return input_vector[mid];
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}
			auto res = two_hop_label(-1);
			res.distance = std::numeric_limits<int>::max();
			res.vertex = -1;
			return res;
		}
		two_hop_label search_sorted_two_hop_label_in_current_mark(std::vector<two_hop_label> &input_vector, int key,bool isRuc)
		{
			if(isRuc){
				ruc_label_query_count ++;
			}else{
				a2021_label_query_count++;
			}
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key)
					{
						return input_vector[mid];
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}
			auto res = two_hop_label(-1);
			res.distance = std::numeric_limits<int>::max();
			res.vertex = -1;
			return res;
		}


		two_hop_label search_sorted_two_hop_label_entity_not_realTime(std::vector<two_hop_label> &input_vector, int key, int time)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key && input_vector[mid].t_s != time)
					{
						return input_vector[mid];
					}
					else if (input_vector[mid].vertex == key && input_vector[mid].t_s == time)
					{
						break;
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}
			auto res = two_hop_label(-1);
			res.distance = std::numeric_limits<int>::max();
			res.vertex = -1;
			return res;
		}

		void insert_sorted_two_hop_label(std::vector<two_hop_label> &input_vector, int key, int value, int time)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key)
					{
						two_hop_label old_label = input_vector[mid];
						old_label.t_e = time - 1;

						input_vector[mid].distance = value;
						input_vector[mid].t_s = time;
						if (old_label.distance != MAX_VALUE)
						{
							int insert_left = mid + 1, insert_right = input_vector.size() - 1;

							while (insert_left <= insert_right)
							{
								int insert_mid = insert_left + ((insert_right - insert_left) / 2);

								if (input_vector[insert_mid].t_e < time)
								{
									insert_right = insert_mid - 1;
								}
								else if (input_vector[insert_mid].t_e == time)
								{
									if (input_vector[insert_mid].vertex > key)
									{
										insert_right = insert_mid - 1;
									}
									else if (input_vector[insert_mid].vertex < key)
									{
										insert_left = insert_mid + 1;
									}
									else
									{
										insert_left = insert_mid;
										break;
									}
								}
								else
								{
									insert_left = insert_mid + 1;
								}
							}
							input_vector.insert(input_vector.begin() + insert_left, old_label);
						}
						return;
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}
			if (value != MAX_VALUE)
			{
				two_hop_label new_label(time);
				new_label.vertex = key;
				new_label.distance = value;

				input_vector.insert(input_vector.begin() + left, new_label);
			}
		}
		void insert_sorted_two_hop_label_mark(std::vector<two_hop_label> &input_vector, int key, int value, int time,bool isRuc)
		{
			if(isRuc){
				ruc_label_insert_count++;
			}else{
				a2021_label_insert_count++;
			}
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].vertex == key)
					{
						two_hop_label old_label = input_vector[mid];
						old_label.t_e = time - 1;

						input_vector[mid].distance = value;
						input_vector[mid].t_s = time;
						if (old_label.distance != MAX_VALUE)
						{
							int insert_left = mid + 1, insert_right = input_vector.size() - 1;

							while (insert_left <= insert_right)
							{
								int insert_mid = insert_left + ((insert_right - insert_left) / 2);

								if (input_vector[insert_mid].t_e < time)
								{
									insert_right = insert_mid - 1;
								}
								else if (input_vector[insert_mid].t_e == time)
								{
									if (input_vector[insert_mid].vertex > key)
									{
										insert_right = insert_mid - 1;
									}
									else if (input_vector[insert_mid].vertex < key)
									{
										insert_left = insert_mid + 1;
									}
									else
									{
										insert_left = insert_mid;
										break;
									}
								}
								else
								{
									insert_left = insert_mid + 1;
								}
							}
							input_vector.insert(input_vector.begin() + insert_left, old_label);
						}
						return;
					}
					else if (input_vector[mid].vertex < key)
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
					right = mid - 1;
				}
			}
			if (value != MAX_VALUE)
			{
				two_hop_label new_label(time);
				new_label.vertex = key;
				new_label.distance = value;

				input_vector.insert(input_vector.begin() + left, new_label);
			}
		}


		class two_hop_case_info
		{
		public:
			int thread_num = 1;

			/*labels*/
			std::vector<std::vector<two_hop_label>> L;
			PPR_TYPE::PPR_type PPR;
			bool operator==(const two_hop_case_info &other) const
			{
				if (thread_num != other.thread_num || L != other.L || PPR != other.PPR)
				{
					return false;
				}
				return true;
			}
			void serialize(std::ofstream &out) const
			{
				experiment::saveBinary(out, thread_num);
				size_t size = L.size();
				BinarySerializer<size_t>::saveBinary(out, size);
				for (auto &item : L)
				{
					BinarySerializer<std::vector<experiment::nonhop::two_hop_label>>::saveBinary(out, item);
					out.flush();
				}
				experiment::saveBinary(out, PPR);
			}

			void deserialize(std::ifstream &in)
			{
				experiment::loadBinary(in, thread_num);
				size_t size;
				BinarySerializer<size_t>::loadBinary(in, size);
				L.resize(size);
				for (auto &item : L)
				{
					BinarySerializer<std::vector<experiment::nonhop::two_hop_label>>::loadBinary(in, item);
				}
				experiment::loadBinary(in, PPR);
			}

			long long int compute_L_size()
			{
				long long int res = 0;
				for (const auto &L_info : L)
				{
					for (const auto &inner : L_info)
					{
						++res;
					}
				}
				return res;
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
						// TODO-GPY 如果正确要修改这里
						// size = size + (PPR[i][j].second.size() + 1) * sizeof(int);
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
			// void print_PPR()
			// {
			// 	std::cout << "print_PPR:" << std::endl;
			// 	for (int i = 0; i < PPR.size(); i++)
			// 	{
			// 		for (int j = 0; j < PPR[i].size(); j++)
			// 		{
			// 			std::cout << "PPR(" << i << "," << PPR[i][j].first << "): ";
			// 			for (int k = 0; k < PPR[i][j].second.size(); k++)
			// 			{
			// 				std::cout << PPR[i][j].second[k] << " ";
			// 			}
			// 			std::cout << std::endl;
			// 		}
			// 	}
			// }

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

			void record_all_details_stream(std::ofstream &outputFile)
			{
				outputFile << "PLL info:" << std::endl;
				outputFile << "thread_num=" << thread_num << std::endl;
				outputFile << "compute_label_byte_size()=" << compute_L_byte_size() << std::endl;
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
	namespace hop
	{
		class two_hop_label
		{
		public:
			int hub_vertex, hop;
			WEIGHT_TYPE distance;
			int t_s, t_e;

			two_hop_label(const two_hop_label &other)
			{
				t_s = other.t_s;
				t_e = other.t_e;
				hub_vertex = other.hub_vertex;
				hop = other.hop;
				distance = other.distance;
			}
			// hop_constrained_two_hop_label() {}
			// hop_constrained_two_hop_label(int _vertex, int _hop, int _dis)
			// {
			//     hub_vertex = _vertex;
			//     hop = _hop;
			//     distance = _dis;
			// }
			two_hop_label()
				: t_s(0), t_e(std::numeric_limits<int>::max())
			{
				hub_vertex = 0;
				hop = 0;
				distance = std::numeric_limits<WEIGHT_TYPE>::max();
			}
			void serialize(std::ofstream &out) const
			{
				experiment::saveBinary(out, hub_vertex);
				experiment::saveBinary(out, hop);
				experiment::saveBinary(out, distance);
				experiment::saveBinary(out, t_s);
				experiment::saveBinary(out, t_e);
			}

			void deserialize(std::ifstream &in)
			{
				experiment::loadBinary(in, hub_vertex);
				experiment::loadBinary(in, hop);
				experiment::loadBinary(in, distance);
				experiment::loadBinary(in, t_s);
				experiment::loadBinary(in, t_e);
			}
		};

		bool compare_hop_constrained_two_hop_label(two_hop_label &i, two_hop_label &j)
		{
			if (i.t_e != j.t_e)
			{
				return i.t_e > j.t_e;
			}
			else if (i.hub_vertex != j.hub_vertex)
			{
				return i.hub_vertex < j.hub_vertex;
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

		bool operator<(two_hop_label const &x, two_hop_label const &y)
		{
			if (x.distance != y.distance)
			{
				return x.distance > y.distance; // < is the max-heap; > is the min heap
			}
			else
			{
				return x.hop > y.hop; // < is the max-heap; > is the min heap
			}
		}

		std::pair<weightTYPE, int> hop_constrained_extract_distance_and_hub(std::vector<std::vector<two_hop_label>> &L, int source, int terminal, int hop_cst)
		{

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
			auto pointer_L_s_end = vector1_check_pointer, pointer_L_t_end = vector2_check_pointer;
			while (pointer_L_s_end != L[source].end() && pointer_L_s_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_s_end++;
			}
			while (pointer_L_t_end != L[terminal].end() && pointer_L_t_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_t_end++;
			}

			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
				{
					auto vector1_end = vector1_check_pointer;
					while (vector1_end != pointer_L_s_end && vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end->t_e == std::numeric_limits<int>::max())
					{
						vector1_end++;
					}
					auto vector2_end = vector2_check_pointer;
					while (vector2_end != pointer_L_t_end && vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end->t_e == std::numeric_limits<int>::max())
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

		std::pair<weightTYPE, int> hop_constrained_extract_distance_and_hop(std::vector<std::vector<two_hop_label>> &L, int source, int terminal, int hop_cst)
		{
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
			int hop = std::numeric_limits<int>::max();
			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = vector1_check_pointer, pointer_L_t_end = vector2_check_pointer;
			while (pointer_L_s_end != L[source].end() && pointer_L_s_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_s_end++;
			}
			while (pointer_L_t_end != L[terminal].end() && pointer_L_t_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_t_end++;
			}

			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
				{
					auto vector1_end = vector1_check_pointer;
					while (vector1_end != pointer_L_s_end && vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end->t_e == std::numeric_limits<int>::max())
					{
						vector1_end++;
					}
					auto vector2_end = vector2_check_pointer;
					while (vector2_end != pointer_L_t_end && vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end->t_e == std::numeric_limits<int>::max())
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
									hop = vector1_begin->hop + vector2_begin->hop;
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

			return {distance, hop};
		}

		std::tuple<weightTYPE, int, int> graph_weighted_two_hop_extract_distance_and_hop_and_hub_by_backup_label(std::vector<two_hop_label> &L_s, std::vector<two_hop_label> &L_t, int hop_cst)
		{
			/*return std::numeric_limits<int>::max() is not connected*/

			if (hop_cst < 0)
			{
				return {std::numeric_limits<int>::max(), -1, -1};
			}
			else if (hop_cst == 0)
			{
				return {std::numeric_limits<int>::max(), -1, -1};
			}

			int distance = std::numeric_limits<int>::max();
			int hop = std::numeric_limits<int>::max();
			int hub = -1;
			auto vector1_check_pointer = L_s.begin();
			auto vector2_check_pointer = L_t.begin();
			auto pointer_L_s_end = vector1_check_pointer, pointer_L_t_end = vector2_check_pointer;
			while (pointer_L_s_end != L_s.end() && pointer_L_s_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_s_end++;
			}
			while (pointer_L_t_end != L_t.end() && pointer_L_t_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_t_end++;
			}

			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
				{
					auto vector1_end = vector1_check_pointer;
					while (vector1_end != pointer_L_s_end && vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end->t_e == std::numeric_limits<int>::max())
					{
						vector1_end++;
					}
					auto vector2_end = vector2_check_pointer;
					while (vector2_end != pointer_L_t_end && vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end->t_e == std::numeric_limits<int>::max())
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
									hop = vector1_begin->hop + vector2_begin->hop;
									hub = vector1_begin->hub_vertex;
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

			return {distance, hop, hub};
		}

		std::pair<weightTYPE, int> graph_weighted_two_hop_extract_distance_and_hop_by_backup_label(std::vector<two_hop_label> &L_s, std::vector<two_hop_label> &L_t, int hop_cst)
		{
			/*return std::numeric_limits<int>::max() is not connected*/

			if (hop_cst < 0)
			{
				return {std::numeric_limits<int>::max(), -1};
			}
			else if (hop_cst == 0)
			{
				return {std::numeric_limits<int>::max(), -1};
			}

			int distance = std::numeric_limits<int>::max();
			int hop = std::numeric_limits<int>::max();
			auto vector1_check_pointer = L_s.begin();
			auto vector2_check_pointer = L_t.begin();
			auto pointer_L_s_end = vector1_check_pointer, pointer_L_t_end = vector2_check_pointer;
			while (pointer_L_s_end != L_s.end() && pointer_L_s_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_s_end++;
			}
			while (pointer_L_t_end != L_t.end() && pointer_L_t_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_t_end++;
			}

			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
				{
					auto vector1_end = vector1_check_pointer;
					while (vector1_end != pointer_L_s_end && vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end->t_e == std::numeric_limits<int>::max())
					{
						vector1_end++;
					}
					auto vector2_end = vector2_check_pointer;
					while (vector2_end != pointer_L_t_end && vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end->t_e == std::numeric_limits<int>::max())
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
									hop = vector1_begin->hop + vector2_begin->hop;
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

			return {distance, hop};
		}

		std::pair<weightTYPE, int> graph_weighted_two_hop_extract_distance_and_hub_by_backup_label(std::vector<two_hop_label> &L_s, std::vector<two_hop_label> &L_t, int hop_cst)
		{
			/*return std::numeric_limits<int>::max() is not connected*/

			if (hop_cst < 0)
			{
				return {std::numeric_limits<int>::max(), -1};
			}
			else if (hop_cst == 0)
			{
				return {std::numeric_limits<int>::max(), -1};
			}

			int distance = std::numeric_limits<int>::max();
			int hub = -1;
			auto vector1_check_pointer = L_s.begin();
			auto vector2_check_pointer = L_t.begin();
			auto pointer_L_s_end = vector1_check_pointer, pointer_L_t_end = vector2_check_pointer;
			while (pointer_L_s_end != L_s.end() && pointer_L_s_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_s_end++;
			}
			while (pointer_L_t_end != L_t.end() && pointer_L_t_end->t_e == std::numeric_limits<int>::max())
			{
				pointer_L_t_end++;
			}
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex)
				{
					auto vector1_end = vector1_check_pointer;
					while (vector1_end != pointer_L_s_end && vector1_check_pointer->hub_vertex == vector1_end->hub_vertex && vector1_end->t_e == std::numeric_limits<int>::max())
					{
						vector1_end++;
					}
					auto vector2_end = vector2_check_pointer;
					while (vector2_end != pointer_L_t_end && vector2_check_pointer->hub_vertex == vector2_end->hub_vertex && vector2_end->t_e == std::numeric_limits<int>::max())
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
									hub = vector1_begin->hub_vertex;
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

			return {distance, hub};
		}

		weightTYPE search_sorted_hop_constrained_weight_two_hop_label(std::vector<two_hop_label> &input_vector, int key, int hop)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].hub_vertex == key)
					{
						if (input_vector[mid].hop == hop)
						{
							return input_vector[mid].distance;
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
					else if (input_vector[mid].hub_vertex < key)
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
					right = mid - 1;
				}
			}

			return std::numeric_limits<int>::max();
		}

		two_hop_label search_sorted_hop_constrained_label_two_hop_label(std::vector<two_hop_label> &input_vector, int key, int hop)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].hub_vertex == key)
					{
						if (input_vector[mid].hop == hop)
						{
							return input_vector[mid];
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
					else if (input_vector[mid].hub_vertex < key)
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
					right = mid - 1;
				}
			}
			two_hop_label res;
			return res;
		}

		std::pair<weightTYPE, int> get_shortest_distance_hop_two_hop_label2(std::vector<two_hop_label> &input_vector, int key)
		{
			int left = 0, right = input_vector.size() - 1;
			weightTYPE mindis = std::numeric_limits<int>::max();
			int hop_val = 0;

			while (left <= right)
			{
				int mid = (right - left) / 2 + left;
				if (input_vector[mid].t_e != std::numeric_limits<int>::max())
				{
					right = mid - 1;
				}
				else
				{
					if (input_vector[mid].hub_vertex < key)
					{
						left = mid + 1;
					}
					else if (input_vector[mid].hub_vertex > key)
					{
						right = mid - 1;
					}
					else
					{
						mindis = input_vector[mid].distance;
						hop_val = input_vector[mid].hop;
						left = mid + 1;
					}
				}
			}

			return {mindis, hop_val};
		}

		std::pair<weightTYPE, int> get_shortest_distance_hop_two_hop_label2(std::vector<two_hop_label> &input_vector, int key, int hop_k)
		{
			int left = 0, right = input_vector.size() - 1;
			weightTYPE mindis = std::numeric_limits<int>::max();
			int hop_val = 0;

			while (left <= right)
			{
				int mid = (right - left) / 2 + left;
				if (input_vector[mid].t_e != std::numeric_limits<int>::max())
				{
					right = mid - 1;
				}
				else
				{
					if (input_vector[mid].hub_vertex < key)
					{
						left = mid + 1;
					}
					else if (input_vector[mid].hub_vertex > key)
					{
						right = mid - 1;
					}
					else
					{
						mindis = input_vector[mid].distance;
						hop_val = input_vector[mid].hop;
						left = mid + 1;
						if (hop_val == hop_k)
						{
							return {mindis, hop_val};
						}
					}
				}
			}

			return {mindis, hop_val};
		}

		void insert_sorted_hop_constrained_two_hop_label(std::vector<two_hop_label> &input_vector, int key, int hop, weightTYPE new_distance, int t)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].hub_vertex == key)
					{
						if (input_vector[mid].hop == hop)
						{
							two_hop_label old_label = input_vector[mid];
							old_label.t_e = t - 1;
							if (input_vector[mid].distance < new_distance && new_distance != MAX_VALUE)
							{
								std::cout << "error input " << std::endl;
								return;
							}
							input_vector[mid].distance = new_distance;
							input_vector[mid].t_s = t;
							if (old_label.t_s != t && old_label.distance != MAX_VALUE)
							{
								int insert_left = mid + 1, insert_right = input_vector.size() - 1;

								while (insert_left <= insert_right)
								{
									int insert_mid = insert_left + ((insert_right - insert_left) / 2);

									if (input_vector[insert_mid].t_e < t - 1)
									{
										insert_right = insert_mid - 1;
									}
									else if (input_vector[insert_mid].t_e == t - 1)
									{
										if (input_vector[insert_mid].hub_vertex > key ||
											(input_vector[insert_mid].hub_vertex == key && input_vector[insert_mid].hop > hop))
										{
											insert_right = insert_mid - 1;
										}
										else
										{
											insert_left = insert_mid + 1;
										}
									}
									else
									{
										insert_left = insert_mid + 1;
									}
								}

								input_vector.insert(input_vector.begin() + insert_left, old_label);
							}
							return;
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
					else if (input_vector[mid].hub_vertex < key)
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
					right = mid - 1;
				}
			}
			if (new_distance != MAX_VALUE)
			{
				two_hop_label new_label;
				new_label.hub_vertex = key;
				new_label.hop = hop;
				new_label.distance = new_distance;
				new_label.t_s = t;
				new_label.t_e = std::numeric_limits<int>::max();

				input_vector.insert(input_vector.begin() + left, new_label);
			}
		}

		std::pair<weightTYPE, int> search_sorted_hop_constrained_weight_and_index_two_hop_label(std::vector<two_hop_label> &input_vector, int key, int hop)
		{
			int left = 0, right = input_vector.size() - 1;

			while (left <= right)
			{
				int mid = left + ((right - left) / 2);

				if (input_vector[mid].t_e == std::numeric_limits<int>::max())
				{
					if (input_vector[mid].hub_vertex == key)
					{
						if (input_vector[mid].hop == hop)
						{
							return {input_vector[mid].distance, mid};
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
					else if (input_vector[mid].hub_vertex < key)
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
					right = mid - 1;
				}
			}

			return {std::numeric_limits<int>::max(), -1};
		}

		class two_hop_case_info
		{
		public:
			/*hop bounded*/
			int thread_num = 1;
			int upper_k = 0;

			/*labels*/
			std::vector<std::vector<two_hop_label>> L;
			PPR_TYPE::PPR_type PPR;

			void serialize(std::ofstream &out) const
			{
				experiment::saveBinary(out, thread_num);
				experiment::saveBinary(out, upper_k);
				size_t size = L.size();
				BinarySerializer<size_t>::saveBinary(out, size);
				for (auto &item : L)
				{
					BinarySerializer<std::vector<experiment::hop::two_hop_label>>::saveBinary(out, item);
					out.flush();
				}
				experiment::saveBinary(out, PPR);
			}

			void deserialize(std::ifstream &in)
			{
				experiment::loadBinary(in, thread_num);
				experiment::loadBinary(in, upper_k);
				size_t size;
				BinarySerializer<size_t>::loadBinary(in, size);
				L.resize(size);
				for (auto &item : L)
				{
					BinarySerializer<std::vector<experiment::hop::two_hop_label>>::loadBinary(in, item);
				}
				experiment::loadBinary(in, PPR);
			}

			long long int compute_label_bit_size()
			{
				long long int size = 0;
				for (auto &xx : L)
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
				for (auto &xx : L)
				{
					std::cout << "vertex " << index++ << ": ";
					for (auto &yy : xx)
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

			void record_all_details_stream(std::ofstream &outputFile)
			{
				outputFile << "hop_constrained_case_info:" << std::endl;
				outputFile << "thread_num=" << thread_num << std::endl;
				outputFile << "upper_k=" << upper_k << std::endl;

				outputFile << "compute_label_bit_size()=" << compute_label_bit_size() << std::endl;
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

	template <>
	class BinarySerializer<nonhop::two_hop_label>
	{
	public:
		static void saveBinary(std::ofstream &out, const nonhop::two_hop_label &label)
		{
			label.serialize(out);
		}

		static void loadBinary(std::ifstream &in, nonhop::two_hop_label &label)
		{
			label.deserialize(in);
		}
	};
	template <>
	class BinarySerializer<hop::two_hop_label>
	{
	public:
		static void saveBinary(std::ofstream &out, const hop::two_hop_label &label)
		{
			label.serialize(out);
		}

		static void loadBinary(std::ifstream &in, hop::two_hop_label &label)
		{
			label.deserialize(in);
		}
	};

	template <>
	class BinarySerializer<hop::two_hop_case_info>
	{
	public:
		static void saveBinary(std::ofstream &out, const hop::two_hop_case_info &info)
		{
			info.serialize(out);
		}

		static void loadBinary(std::ifstream &in, hop::two_hop_case_info &info)
		{
			info.deserialize(in);
		}
	};

	template <>
	class BinarySerializer<nonhop::two_hop_case_info>
	{
	public:
		static void saveBinary(std::ofstream &out, const nonhop::two_hop_case_info &vec)
		{
			vec.serialize(out);
		}

		static void loadBinary(std::ifstream &in, nonhop::two_hop_case_info &vec)
		{
			vec.deserialize(in);
		}
	};
}