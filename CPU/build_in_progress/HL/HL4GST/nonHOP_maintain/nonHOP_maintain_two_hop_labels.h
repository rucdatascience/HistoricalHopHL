#pragma once
#include <vector>
#include <string>
#include <shared_mutex>
#include <boost/heap/fibonacci_heap.hpp>
#include <chrono>
#include <algorithm>
#include "CPU/tool_functions/ThreadPool.h"
#include "CPU/build_in_progress/HL/HL4GST/nonHOP_maintain/nonHOP_PPR.h"
#include "Historical/graph_v_of_v/Timer.h"
using namespace std;
#define weightTYPE int
namespace nonHop
{
	class two_hop_label
	{
	public:
		int vertex, distance;
		int t_s, t_e;
		two_hop_label(int start_time)
		{
			t_s = start_time;
			t_e = INT_MAX;
			vertex = std::numeric_limits<int>::max();
			distance = std::numeric_limits<int>::max();
		}
	};

	/*global values that should be cleared after usig PLL*/
	string reach_limit_error_string_MB = "reach limit error MB";
	string reach_limit_error_string_time = "reach limit error time";
	long long int max_labal_num_595;
	long long int labal_num_595;
	auto begin_time_595 = std::chrono::high_resolution_clock::now();
	double max_run_time_nanoseconds_595;
	bool this_parallel_PLL_is_running_595 = false;
	vector<vector<two_hop_label>> L_temp_595;
	PPR_type PPR_595;
	std::shared_mutex mtx_595_1, mtx_595_2;
	int max_N_ID_for_mtx_595 = 1e7;							 // this is the max N to run
	vector<std::shared_mutex> mtx_595(max_N_ID_for_mtx_595); // std::mutex has no copy or move constructor, while std::vector::resize() requires that; you cannot resize mtx; moreover, do not change mtx to a pointer and then points to local values, it is very slow!!
	vector<int> is_mock(max_N_ID_for_mtx_595);
	queue<int> Qid_595; // IDs of available elements of P T
	vector<vector<int>> P_dij_595;
	vector<vector<int>> T_dij_595;
	long long int label_size_before_canonical_repair_595 = 0;
	long long int label_size_after_canonical_repair_595 = 0;
	bool global_use_2M_prune = false, global_use_rank_prune = true;
	int TwoM_value = 2 * 1e6; // suppose that dummy edge has a weight of 1e6

	bool operator<(two_hop_label const& x, two_hop_label const& y)
	{
		return x.distance > y.distance; // < is the max-heap; > is the min heap
	}
	typedef typename boost::heap::fibonacci_heap<two_hop_label>::handle_type PLL_handle_t_for_sp;
	vector<vector<PLL_handle_t_for_sp>> Q_handles_595;

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

	class two_hop_case_info
	{
	public:
		/*use info*/
		bool use_2M_prune = false;
		bool use_rank_prune = true;
		bool use_canonical_repair = false;
		int thread_num = 1;

		/*canonical_repair info*/
		long long int label_size_before_canonical_repair = 0;
		long long int label_size_after_canonical_repair = 0;
		double canonical_repair_remove_label_ratio = 0;

		/*running time records*/
		double time_initialization = 0;
		double time_generate_labels = 0;
		double time_sortL = 0;
		double time_canonical_repair = 0;
		double time_total = 0;

		/*running limits*/
		long long int max_labal_byte_size = 1e12;
		double max_run_time_seconds = 1e12;

		/**
		 * query param
		 */
		int source;
		int target;
		int t_s;
		int t_e;

		/* array of maintenance algorithm running time each time*/
		vector<double> time_decrease;
		vector<double> time_increase;

		double get_maintain_time()
		{
			double res = 0;
			for (double time : this->time_decrease)
			{
				res += time;
			}
			for (double time : this->time_increase)
			{
				res += time;
			}
			return res;
		}

		/*labels*/
		vector<vector<two_hop_label>> L;
		PPR_type PPR;

		/*clear labels*/
		void clear_labels()
		{
			vector<vector<two_hop_label>>().swap(L);
			PPR_type().swap(PPR);
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
			cout << "print_L:" << endl;
			for (int i = 0; i < L.size(); i++)
			{
				cout << "L[" << i << "]=";
				for (int j = 0; j < L[i].size(); j++)
				{
					cout << "{" << L[i][j].vertex << "," << L[i][j].distance << "," << L[i][j].t_s << "," << L[i][j].t_e << "}";
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

		void print_times()
		{
			cout << "print_times:" << endl;
			cout << "time_initialization: " << time_initialization << "s" << endl;
			cout << "time_generate_labels: " << time_generate_labels << "s" << endl;
			cout << "time_sortL: " << time_sortL << "s" << endl;
			cout << "time_canonical_repair: " << time_canonical_repair << "s" << endl;
		}

		/*record_all_details*/
		void record_all_details(string save_name)
		{
			ofstream outputFile;
			outputFile.precision(6);
			outputFile.setf(ios::fixed);
			outputFile.setf(ios::showpoint);
			outputFile.open(save_name + ".txt");

			outputFile << "PLL info:" << endl;
			outputFile << "thread_num=" << thread_num << endl;
			outputFile << "use_2M_prune=" << use_2M_prune << endl;
			outputFile << "use_rank_prune=" << use_rank_prune << endl;
			outputFile << "use_canonical_repair=" << use_canonical_repair << endl;

			outputFile << "label_size_before_canonical_repair=" << label_size_before_canonical_repair << endl;
			outputFile << "label_size_after_canonical_repair=" << label_size_after_canonical_repair << endl;
			outputFile << "canonical_repair_remove_label_ratio=" << canonical_repair_remove_label_ratio << endl;

			outputFile << "time_initialization=" << time_initialization << endl;
			outputFile << "time_generate_labels=" << time_generate_labels << endl;
			outputFile << "time_sortL=" << time_sortL << endl;
			outputFile << "time_canonical_repair=" << time_canonical_repair << endl;
			outputFile << "time_total=" << time_total << endl;

			outputFile << "max_labal_byte_size=" << max_labal_byte_size << endl;
			outputFile << "max_run_time_seconds=" << max_run_time_seconds << endl;

			outputFile << "compute_label_byte_size()=" << compute_L_byte_size() << endl;

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
					if (vector1_begin->vertex == vector2_begin->vertex && max(vector1_begin->t_s, max(vector2_begin->t_s, t_s)) <= min(vector1_begin->t_e, min(vector2_begin->t_e, t_e)))
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

	int _2023algo_max_second = 60;

	auto begin_time = std::chrono::high_resolution_clock::now();
	double max_run_time_nanosec;
	string reach_limit_time_string = "reach limit time in WeightIncrease";

	/*common functions*/

	long long int label_operation_times = 0;

	bool compare_two_hop_label_small_to_large(two_hop_label& i, two_hop_label& j)
	{
		if (i.t_e != j.t_e)
			return i.t_e > j.t_e;	// t_e降序
		return i.vertex < j.vertex; // < is from small to big; > is from big to small
	}

	void insert_sorted_two_hop_label(std::vector<two_hop_label>& input_vector, int key, weightTYPE value, int time)
	{
		label_operation_times++;
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
					return;
				}
				else if (input_vector[mid].vertex < key)
				{
					left = mid + 1; // hub_vertex 太小，继续向右找
				}
				else
				{
					right = mid - 1; // hub_vertex 太大，向左找
				}
			}
			else
			{
				right = mid - 1; // t_e 不是无穷大，向左继续查找
			}
		}

		// Step 6: 如果没找到符合条件的标签，则插入新标签
		two_hop_label new_label(time);
		new_label.vertex = key;
		new_label.distance = value;

		// 找到插入点，确保排序规则
		input_vector.insert(input_vector.begin() + left, new_label); // 在 left 位置插入新标签
	}

	bool search_sorted_two_hop_label_specify_time_span_cost(std::vector<two_hop_label>& input_vector, int key, int cost, int t_s, int t_e)
	{
		// TODO: 这里使用的是傻瓜方法 之后修复 可以降到对数级别
		for (int i = 0; i < input_vector.size(); i++)
		{
			if (input_vector[i].vertex == key && input_vector[i].distance == cost && (input_vector[i].t_s, t_s) <= min(input_vector[i].t_e, t_e))
			{
				return true;
			}
		}
		return false;
	}

	two_hop_label search_sorted_two_hop_label_entity(std::vector<two_hop_label>& input_vector, int key)
	{
		label_operation_times++;
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

	two_hop_label search_sorted_two_hop_label_entity_not_realTime(std::vector<two_hop_label>& input_vector, int key, int time)
	{
		label_operation_times++;
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
				else if (input_vector[mid].t_s == time) {
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

	weightTYPE search_sorted_two_hop_label(std::vector<two_hop_label>& input_vector, int key)
	{

		label_operation_times++;
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
	}

	pair<weightTYPE, int> search_sorted_two_hop_label2(std::vector<two_hop_label>& input_vector, int key)
	{

		label_operation_times++;
		int left = 0, right = input_vector.size() - 1;

		while (left <= right)
		{
			int mid = left + ((right - left) / 2);

			if (input_vector[mid].t_e == std::numeric_limits<int>::max())
			{
				if (input_vector[mid].vertex == key)
				{
					return { input_vector[mid].distance, mid };
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

		return { std::numeric_limits<int>::max(), -1 };
	}

#define Query(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, x, y)	 // reduction is not used here
#define Query2(x, y) graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(L, x, y) // reduction is not used here
	// #define MAX_VALUE std::numeric_limits<int>::max()
#define MAX_VALUE 1e7

	class affected_label
	{
	public:
		int first, second;
		weightTYPE dis;
		int t_s, t_e;
		affected_label() {}
		affected_label(int _first, int _second, weightTYPE _dis)
		{
			first = _first;
			second = _second;
			dis = _dis;
		}
	};

	class pair_label
	{ // pair_label2 is stored in NoP
	public:
		int first, second;
		pair_label(int _first, int _second)
		{
			first = _first;
			second = _second;
		}
		bool operator==(const pair_label other) const
		{
			return (first == other.first && second == other.second);
		}
		bool operator<(const pair_label other) const
		{ // used to sort/search pair_label2 in set
			if (first != other.first)
				return first < other.first;
			return second < other.second;
		}
	};

	struct node_for_DIFFUSE
	{
		int index;
		weightTYPE disx;
		node_for_DIFFUSE() {}
		node_for_DIFFUSE(int _u, weightTYPE _dis)
		{
			index = _u;
			disx = _dis;
		}
	}; // define the node in the queue

	bool operator<(node_for_DIFFUSE const& x, node_for_DIFFUSE const& y)
	{
		return x.disx > y.disx; // < is the max-heap; > is the min heap
	}
	typedef typename boost::heap::fibonacci_heap<node_for_DIFFUSE>::handle_type handle_t_for_DIFFUSE; // pairing heap has a similar speed with fibonacci_heap here

	vector<std::shared_mutex> mtx_5952(max_N_ID_for_mtx_595);
	vector<vector<pair<weightTYPE, int>>> Dis;
	vector<vector<weightTYPE>> Q_value;
	vector<vector<handle_t_for_DIFFUSE>> Q_handles;

	void initialize_global_values_dynamic(int N, int thread_num)
	{
		Dis.resize(thread_num);
		Q_value.resize(thread_num);
		Q_handles.resize(thread_num);
		queue<int>().swap(Qid_595);
		for (int i = 0; i < thread_num; i++)
		{
			Dis[i].resize(N, { -1, -1 });
			Q_value[i].resize(N, MAX_VALUE);
			Q_handles[i].resize(N);
			Qid_595.push(i);
		}
	}

	/*codes for querying distances*/

	long long int global_query_times = 0;

	weightTYPE graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(vector<vector<two_hop_label>>& L, int source, int terminal)
	{

		global_query_times++;

		/*return std::numeric_limits<double>::max() is not connected*/

		if (source == terminal)
		{
			return 0;
		}

		weightTYPE distance = std::numeric_limits<weightTYPE>::max(); // if disconnected, return this large value

		auto vector1_check_pointer = L[source].begin();
		auto vector2_check_pointer = L[terminal].begin();
		auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
		{
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
			{
				weightTYPE dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
				if (distance > dis)
				{
					distance = dis;
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

		return distance;
	}

	pair<weightTYPE, int> graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(vector<vector<two_hop_label>>& L, int source, int terminal)
	{

		global_query_times++;

		/*return std::numeric_limits<double>::max() is not connected*/

		if (source == terminal)
		{
			return { 0, source };
		}

		weightTYPE distance = std::numeric_limits<weightTYPE>::max(); // if disconnected, return this large value
		int common_hub;

		auto vector1_check_pointer = L[source].begin();
		auto vector2_check_pointer = L[terminal].begin();
		auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
		{
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
			{
				weightTYPE dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
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

		return { distance, common_hub };
	}

	pair<vector<pair<weightTYPE, weightTYPE>>, vector<int>> graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2_find_all_hub(vector<vector<two_hop_label>>& L, int source, int terminal)
	{

		global_query_times++;

		/*return std::numeric_limits<double>::max() is not connected*/

		vector<int> common_hub;
		if (source == terminal)
		{
			common_hub.push_back(source);
			return { {{0, 0}}, common_hub };
		}

		vector<pair<weightTYPE, weightTYPE>> distanceList = { {-1, -1} }; // if disconnected, return this large value
		weightTYPE distance = std::numeric_limits<weightTYPE>::max();
		auto vector1_check_pointer = L[source].begin();
		auto vector2_check_pointer = L[terminal].begin();
		auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
		{
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
			{
				weightTYPE dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
				if (distance > dis)
				{
					distance = dis;
					vector<pair<weightTYPE, weightTYPE>>().swap(distanceList);
					vector<pair<weightTYPE, weightTYPE>>{{vector1_check_pointer->distance, vector2_check_pointer->distance}}.swap(distanceList);
					vector<int>().swap(common_hub);
					common_hub.push_back(vector1_check_pointer->vertex);
				}
				else if (distance == dis)
				{
					distanceList.push_back({ vector1_check_pointer->distance, vector2_check_pointer->distance });
					common_hub.push_back(vector1_check_pointer->vertex);
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

		return { distanceList, common_hub };
	}

	weightTYPE graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc3(vector<two_hop_label>& L_s, vector<two_hop_label>& L_t)
	{

		global_query_times++;

		/*return std::numeric_limits<double>::max() is not connected*/

		weightTYPE distance = std::numeric_limits<weightTYPE>::max(); // if disconnected, return this large value

		auto vector1_check_pointer = L_s.begin();
		auto vector2_check_pointer = L_t.begin();
		auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
		{
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
			{
				weightTYPE dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
				if (distance > dis)
				{
					distance = dis;
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

		return distance;
	}

	pair<two_hop_label, two_hop_label> graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4_two_hop_label(vector<two_hop_label>& L_s, vector<two_hop_label>& L_t)
	{

		global_query_times++;

		/*return std::numeric_limits<double>::max() is not connected*/

		weightTYPE distance = std::numeric_limits<weightTYPE>::max(); // if disconnected, return this large value
		int common_hub;
		two_hop_label res1 = two_hop_label{ -1 };
		two_hop_label res2 = two_hop_label{ -1 };
		auto vector1_check_pointer = L_s.begin();
		auto vector2_check_pointer = L_t.begin();
		auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
		{
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
			{
				weightTYPE dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
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

		return { res1, res2 };
	}

	pair<weightTYPE, int> graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc4(vector<two_hop_label>& L_s, vector<two_hop_label>& L_t)
	{

		global_query_times++;

		/*return std::numeric_limits<double>::max() is not connected*/

		weightTYPE distance = std::numeric_limits<weightTYPE>::max(); // if disconnected, return this large value
		int common_hub;

		auto vector1_check_pointer = L_s.begin();
		auto vector2_check_pointer = L_t.begin();
		auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end && vector1_check_pointer->t_e == std::numeric_limits<int>::max() && vector2_check_pointer->t_e == std::numeric_limits<int>::max())
		{
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
			{
				weightTYPE dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
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

		return { distance, common_hub };
	}

	/*querying distances or paths*/

	int extract_distance(two_hop_case_info& case_info, int source, int terminal)
	{

		auto& L = case_info.L;

		/*return std::numeric_limits<double>::max() is not connected*/

		if (source == terminal)
		{
			return 0;
		}

		int distance = std::numeric_limits<int>::max(); // if disconnected, return this large value
		auto vector1_check_pointer = L[source].begin();
		auto vector2_check_pointer = L[terminal].begin();
		auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
		{
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
			{
				int dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
				if (distance > dis)
				{
					distance = dis;
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

		return distance;
	}
}