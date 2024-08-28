#pragma once
#include <boost/heap/fibonacci_heap.hpp>
#include <CPU/build_in_progress/HL/HL4GST/HOP_maintain/HOP_maintain_hop_constrained_two_hop_labels.h>
#include <shared_mutex>
#include <CPU/tool_functions/ThreadPool.h>
#include <CPU/graph_v_of_v/graph_v_of_v.h>

#define MAX_VALUE 1e7
int TwoM_value = 2 * 1e6; // suppose that dummy edge has a weight of 1e6
bool global_use_2M_prune = false, global_use_rank_prune = true;
string reach_limit_error_string_MB = "reach limit error MB";
string reach_limit_error_string_time = "reach limit error time";
bool HSDL_dynamic_generate_PPR = true;

/*unique code for this file: 599*/
long long int max_labal_size_599;
/* the number of labels that during the generation process must not exceed the max_labal_size_599 */
long long int labal_size_599;
int max_N_ID_for_mtx_599 = 1e7;
vector<int> is_mock(max_N_ID_for_mtx_599);
/* the maximum runtime in nanoseconds is controlled by program parameters */
double max_run_time_nanoseconds_599;
int global_upper_k;
long long int label_size_before_canonical_repair_599, label_size_after_canonical_repair_599;
/* the time at program startup */
auto begin_time_599 = std::chrono::high_resolution_clock::now();
vector<std::shared_timed_mutex> mtx_599(max_N_ID_for_mtx_599);
graph_v_of_v<int> ideal_graph_599;
vector<vector<hop_constrained_two_hop_label>> L_temp_599;
/**
 * index:[i][j][k]
 * i: threadId j:dest_vertexId k:labelId
 */
vector<vector<vector<pair<int, int>>>> Temp_L_vk_599;
vector<vector<pair<int, int>>> dist_hop_599, dist_hop_599_v2, dist_hop_599_v3;
vector<vector<vector<weightTYPE>>> Q_value;
vector<vector<vector<int>>> Vh_599;
/* mark the process id in the queue and the id ranges from 0 to threadNum*/
queue<int> Qid_599, Qid_599_v2, Qid_599_v3;
PPR_type PPR_599;

typedef typename boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE>::handle_type hop_constrained_handle_t_for_DIFFUSE; // pairing heap has a similar speed with fibonacci_heap here

std::shared_mutex mtx_599_1, mtx_599_2;
vector<std::shared_mutex> mtx_5992(max_N_ID_for_mtx_599);

void initialize_global_values_dynamic_hop_constrained(int N, int thread_num, int upper_k)
{

	dist_hop_599_v2.resize(thread_num);
	dist_hop_599_v3.resize(thread_num);
	Q_value.resize(thread_num);
	queue<int>().swap(Qid_599_v2);
	queue<int>().swap(Qid_599_v3);
	for (int i = 0; i < thread_num; i++)
	{

		Qid_599_v2.push(i);
		Qid_599_v3.push(i);
		dist_hop_599_v2[i].resize(N, {-1, 0});
		dist_hop_599_v3[i].resize(N, {-1, 0});
		Q_value[i].resize(N, vector<weightTYPE>(upper_k + 1, MAX_VALUE));
	}
}

/* override the operator for hop_constrained_two_hop_label to use the Fibonacci minimum heap */
bool operator<(hop_constrained_two_hop_label const &x, hop_constrained_two_hop_label const &y)
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
/**
 * a minimum fibonacci heap for the hop-constrained 2-hop label
 */
typedef typename boost::heap::fibonacci_heap<hop_constrained_two_hop_label>::handle_type hop_constrained_node_handle;
/**
 * index[i][j][k];
 * i: threadNum;
 * j: vertexId
 * k: hop
 * value:label_node,dist
 * in the BFS iteration process for each vertex,the vector will be reset by changed vector
 */
vector<vector<vector<pair<hop_constrained_node_handle, int>>>> Q_handle_priorities_599;

void hop_constrained_clear_global_values()
{
	vector<vector<hop_constrained_two_hop_label>>().swap(L_temp_599);
	ideal_graph_599.clear();
	vector<vector<vector<pair<int, int>>>>().swap(Temp_L_vk_599);
	vector<vector<pair<int, int>>>().swap(dist_hop_599);
	vector<vector<vector<pair<hop_constrained_node_handle, int>>>>().swap(Q_handle_priorities_599);
	vector<vector<vector<int>>>().swap(Vh_599);
	queue<int>().swap(Qid_599);
	PPR_type().swap(PPR_599);
}

/* use asynchronous tasks to generate L labels and the parameter v_k should be the index of the endpoint */
void HSDL_thread_function(int v_k)
{
	// cout << "HSDL_thread_function" << endl;
	if (labal_size_599 > max_labal_size_599)
	{
		throw reach_limit_error_string_MB;
	}
	/* exceeding the maximum time limit*/
	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_599).count() > max_run_time_nanoseconds_599)
	{
		throw reach_limit_error_string_time;
	}

	/* get unique thread id */
	/* critical section obtain array index  */
	mtx_599[max_N_ID_for_mtx_599 - 1].lock();
	int used_id = Qid_599.front();
	Qid_599.pop();
	mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

	/* store the Temp L and dist_hop in current thread*/
	vector<int> Temp_L_vk_changes, dist_hop_changes;
	/* Temp_L_vk stores the dest_vertex_id and distance and hop */
	auto &Temp_L_vk = Temp_L_vk_599[used_id];
	auto &dist_hop = dist_hop_599[used_id]; // record the minimum distance (and the corresponding hop) of a searched vertex in Q
	vector<pair<int, int>> Q_handle_priorities_changes;
	/* get the label list in current thread*/
	auto &Q_handle_priorities = Q_handle_priorities_599[used_id];

	long long int new_label_num = 0;

	/* a class contains information about destination vertex, hop count, and cost in the priority queue */
	boost::heap::fibonacci_heap<hop_constrained_two_hop_label> Q;

	/* generate a label for the starting vertex itself */
	hop_constrained_two_hop_label node;
	node.hub_vertex = v_k;
	node.hop = 0;
	node.distance = 0;
	Q_handle_priorities[v_k][0] = {Q.push({node}), node.distance};
	Q_handle_priorities_changes.push_back({v_k, 0});

	/* Temp_L_vk_599 stores the label (dist and hop) of vertex v_k */
	mtx_599[v_k].lock();
	L_temp_599[v_k].push_back(node);
	new_label_num++;
	/* root is vk-> vk->obj info -> vector<obj> -> index-> vertexId obj-><distance,hop> */
	for (auto &xx : L_temp_599[v_k])
	{
		int L_vk_vertex = xx.hub_vertex;
		Temp_L_vk[L_vk_vertex].push_back({xx.distance, xx.hop});
		Temp_L_vk_changes.push_back(L_vk_vertex);
	}
	mtx_599[v_k].unlock();

	/*  dist_hop_599 stores the shortest distance from vk to any other vertices with its hop_cst,
		note that the hop_cst is determined by the shortest distance */
	dist_hop[v_k] = {0, 0};
	dist_hop_changes.push_back(v_k);

	while (Q.size() > 0)
	{
		/* poll the vertex from heap.In other words, poll the vertex with the minimal cost */
		node = Q.top();
		Q.pop();

		/* current node, u, which is the node generating the labels*/
		int u = node.hub_vertex;

		/* it doesn't judge nodes with degrees larger than the current node for pruning */
		if (global_use_rank_prune && v_k > u)
		{
			continue;
		}

		// if (!global_use_2M_prune) {
		//	if (global_use_rank_prune && v_k > u) {
		//		continue;
		//	}
		// }
		// else {
		//	if (v_k > u && !is_mock[u]) {
		//		continue;
		//	}
		// }

		int u_hop = node.hop;
		int P_u = node.distance;
		int common_hub_for_query_v_k_u = -1;
		int common_hub_for_query_v_k_u_opt = -1;
		int query_v_k_u = std::numeric_limits<int>::max();
		int query_v_k_u_opt = std::numeric_limits<int>::max();

		mtx_599[u].lock();
		for (auto &xx : L_temp_599[u])
		{
			int common_v = xx.hub_vertex;
			for (auto &yy : Temp_L_vk[common_v])
			{
				long long int dis_opt = (long long int)xx.distance + yy.first;
				if (query_v_k_u_opt > dis_opt)
				{
					query_v_k_u_opt = dis_opt;
					common_hub_for_query_v_k_u_opt = xx.hub_vertex;
				}
				if (xx.hop + yy.second <= u_hop)
				{
					long long int dis = (long long int)xx.distance + yy.first;
					if (query_v_k_u > dis)
					{
						query_v_k_u = dis;
						common_hub_for_query_v_k_u = xx.hub_vertex;
					}
					break;
				}
			}
		}
		mtx_599[u].unlock();

		// if (u == 4 && v_k == 0 && u_hop == 3 && P_u == 8) {
		//	cout << "xx: " << P_u << " " << query_v_k_u << endl;
		// }

		// Either a shorter path to the current label exists, or the node is v_k
		if (P_u < query_v_k_u || query_v_k_u == 0)
		{ // query_v_k_u == 0 is to start the while loop by searching neighbors of v_k

			if (P_u < query_v_k_u)
			{
				node.hub_vertex = v_k;
				node.hop = u_hop;
				node.distance = P_u;
				mtx_599[u].lock();
				L_temp_599[u].push_back(node);
				mtx_599[u].unlock();
				new_label_num++;
			}

			// if (HSDL_dynamic_generate_PPR && query_v_k_u_opt < P_u)
			// {
			// 	if (common_hub_for_query_v_k_u_opt != v_k)
			// 	{
			// 		mtx_599[u].lock();
			// 		PPR_insert(PPR_599, u, common_hub_for_query_v_k_u_opt, v_k);
			// 		mtx_599[u].unlock();
			// 	}
			// 	if (common_hub_for_query_v_k_u_opt != u)
			// 	{
			// 		mtx_599[v_k].lock();
			// 		PPR_insert(PPR_599, v_k, common_hub_for_query_v_k_u_opt, u);
			// 		mtx_599[v_k].unlock();
			// 	}
			// }

			if (u_hop + 1 > global_upper_k)
			{
				continue;
			}

			/* update adj */
			/* Traverse neighboring nodes */
			for (auto &xx : ideal_graph_599[u])
			{
				/* adh_v is the neighborhood and the ec is the distance from u to ajd_v*/
				int adj_v = xx.first, ec = xx.second;

				if (global_use_2M_prune && P_u + ec >= TwoM_value)
				{
					continue;
				}

				/* update node info */
				node.hub_vertex = adj_v;
				node.distance = P_u + ec;
				node.hop = u_hop + 1;

				auto &yy = Q_handle_priorities[adj_v][node.hop];

				/*directly using the following codes without dist_hop is OK, but is slower; dist_hop is a pruning technique without increasing the time complexity*/
				// if (yy.second != std::numeric_limits<int>::max()) {
				//	if (yy.second > node.distance) {
				//		Q.update(yy.first, node);
				//		yy.second = node.distance;
				//	}
				// }
				// else {
				//	yy = { Q.push(node), node.distance };
				//	Q_handle_priorities_changes.push_back({ adj_v, node.hop });
				// }

				if (yy.second <= node.distance)
				{ // adj_v has been reached with a smaller distance and the same hop
					continue;
				}
				/* the vertex has not been visited*/
				if (dist_hop[adj_v].first == std::numeric_limits<int>::max())
				{ // adj_v has not been reached
					yy = {Q.push({node}), node.distance};
					Q_handle_priorities_changes.push_back({adj_v, node.hop});
					dist_hop[adj_v].first = node.distance;
					dist_hop[adj_v].second = node.hop;
					dist_hop_changes.push_back(adj_v);
				}
				else
				{
					if (node.distance < dist_hop[adj_v].first)
					{ // adj_v has been reached with a less distance
						if (yy.second != std::numeric_limits<int>::max())
						{
							Q.update(yy.first, node);
							yy.second = node.distance;
						}
						else
						{
							yy = {Q.push(node), node.distance};
							Q_handle_priorities_changes.push_back({adj_v, node.hop});
						}
						dist_hop[adj_v].first = node.distance;
						dist_hop[adj_v].second = node.hop;
					}
					else if (node.hop < dist_hop[adj_v].second)
					{ // adj_v has been reached with a less hop
						if (yy.second != std::numeric_limits<int>::max())
						{
							Q.update(yy.first, node);
							yy.second = node.distance;
						}
						else
						{
							yy = {Q.push(node), node.distance};
							Q_handle_priorities_changes.push_back({adj_v, node.hop});
						}
					}
				}
			}
		}
		else if (HSDL_dynamic_generate_PPR)
		{
			/* add v_k into PPR(u,common_hub_for_query_v_k_u), and add u into PPR(v_k,common_hub_for_query_v_k_u)*/
			if (common_hub_for_query_v_k_u != v_k)
			{
				mtx_599[u].lock();
				PPR_insert(PPR_599, u, common_hub_for_query_v_k_u, v_k);
				mtx_599[u].unlock();
			}
			if (common_hub_for_query_v_k_u != u)
			{
				mtx_599[v_k].lock();
				PPR_insert(PPR_599, v_k, common_hub_for_query_v_k_u, u);
				mtx_599[v_k].unlock();
			}
		}
	}

	for (auto &xx : Temp_L_vk_changes)
	{
		vector<pair<int, int>>().swap(Temp_L_vk[xx]);
	}
	for (auto &xx : dist_hop_changes)
	{
		dist_hop[xx] = {std::numeric_limits<int>::max(), 0};
	}
	hop_constrained_node_handle handle_x;
	for (auto &xx : Q_handle_priorities_changes)
	{
		Q_handle_priorities[xx.first][xx.second] = {handle_x, std::numeric_limits<int>::max()};
	}

	mtx_599[v_k].lock();
	vector<hop_constrained_two_hop_label>(L_temp_599[v_k]).swap(L_temp_599[v_k]);
	mtx_599[v_k].unlock();

	mtx_599[max_N_ID_for_mtx_599 - 1].lock();
	Qid_599.push(used_id);
	labal_size_599 = labal_size_599 + new_label_num;
	mtx_599[max_N_ID_for_mtx_599 - 1].unlock();
}

void _2023WWW_thread_function(int v_k)
{

	if (labal_size_599 > max_labal_size_599)
	{
		throw reach_limit_error_string_MB;
	}
	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_599).count() > max_run_time_nanoseconds_599)
	{
		throw reach_limit_error_string_time;
	}

	/* get unique thread id */
	mtx_599[max_N_ID_for_mtx_599 - 1].lock();
	int used_id = Qid_599.front();
	Qid_599.pop();
	mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

	vector<int> Temp_L_vk_changes, dist_hop_changes;
	auto &Temp_L_vk = Temp_L_vk_599[used_id];
	auto &dist_hop = dist_hop_599[used_id]; // record {dis, predecessor}
	auto &Vh = Vh_599[used_id];

	long long int new_label_num = 0;

	hop_constrained_two_hop_label node;
	node.hub_vertex = v_k;
	node.hop = 0;
	node.distance = 0;

	/* Temp_L_vk_599 stores the label (dist and hop) of vertex v_k */
	mtx_599[v_k].lock();
	// L_temp_599[v_k].push_back(node); new_label_num++;
	for (auto &xx : L_temp_599[v_k])
	{
		int L_vk_vertex = xx.hub_vertex;
		Temp_L_vk[L_vk_vertex].push_back({xx.distance, xx.hop});
		Temp_L_vk_changes.push_back(L_vk_vertex);
	}
	mtx_599[v_k].unlock();

	Vh[0].push_back(v_k);

	dist_hop[v_k] = {0, v_k};
	dist_hop_changes.push_back(v_k);

	vector<tuple<int, int, int>> dh_updates;

	for (int h = 0; h <= global_upper_k; h++)
	{

		for (auto &xx : dh_updates)
		{
			if (dist_hop[get<0>(xx)].first > get<1>(xx))
			{
				dist_hop[get<0>(xx)] = {get<1>(xx), get<2>(xx)};
				dist_hop_changes.push_back(get<0>(xx));
			}
		}
		vector<tuple<int, int, int>>().swap(dh_updates);

		for (auto u : Vh[h])
		{
			int P_u = dist_hop[u].first;

			if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_599).count() > max_run_time_nanoseconds_599)
			{
				throw reach_limit_error_string_time;
			}

			// if (u < v_k) { // rank pruning ò�Ʋ��������Լ���_2023WWW_thread_function��_2023WWW_thread_function֮����������ΪBFS�����˹�����������Ҫ����ϴ��label������ıȷ�����Ķ༸��һ����������
			//	continue;
			// }

			int query_v_k_u = std::numeric_limits<int>::max();
			mtx_599[u].lock();
			for (auto &xx : L_temp_599[u])
			{
				int common_v = xx.hub_vertex;
				for (auto &yy : Temp_L_vk[common_v])
				{
					if (xx.hop + yy.second <= h)
					{
						long long int dis = (long long int)xx.distance + yy.first;
						if (query_v_k_u > dis)
						{
							query_v_k_u = dis;
						}
					}
				}
			}
			mtx_599[u].unlock();

			if (P_u < query_v_k_u)
			{

				node.hub_vertex = v_k;
				node.hop = h;
				node.distance = P_u;
				mtx_599[u].lock();
				L_temp_599[u].push_back(node);
				mtx_599[u].unlock();
				new_label_num++;

				/* update adj */
				for (auto &xx : ideal_graph_599[u])
				{
					int adj_v = xx.first, ec = xx.second;
					if (P_u + ec < dist_hop[adj_v].first)
					{
						Vh[h + 1].push_back(adj_v);
						dh_updates.push_back({adj_v, P_u + ec, u});
					}
				}
			}
		}
	}

	for (auto &xx : Temp_L_vk_changes)
	{
		vector<pair<int, int>>().swap(Temp_L_vk[xx]);
	}
	for (int i = 0; i <= global_upper_k; i++)
	{
		vector<int>().swap(Vh[i]);
	}
	for (auto &xx : dist_hop_changes)
	{
		dist_hop[xx] = {std::numeric_limits<int>::max(), 0};
	}

	mtx_599[v_k].lock();
	vector<hop_constrained_two_hop_label>(L_temp_599[v_k]).swap(L_temp_599[v_k]);
	mtx_599[v_k].unlock();

	mtx_599[max_N_ID_for_mtx_599 - 1].lock();
	Qid_599.push(used_id);
	labal_size_599 = labal_size_599 + new_label_num;
	mtx_599[max_N_ID_for_mtx_599 - 1].unlock();
}

/*sortL*/
// bool compare_hop_constrained_two_hop_label(hop_constrained_two_hop_label &i, hop_constrained_two_hop_label &j)
// {
// 	if (i.hub_vertex != j.hub_vertex)
// 	{
// 		return i.hub_vertex < j.hub_vertex;
// 	}
// 	else if (i.hop != j.hop)
// 	{
// 		return i.hop < j.hop;
// 	}
// 	else
// 	{
// 		return i.distance < j.distance;
// 	}
// }
vector<vector<hop_constrained_two_hop_label>> hop_constrained_sortL(int num_of_threads)
{

	/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/

	int N = L_temp_599.size();
	vector<vector<hop_constrained_two_hop_label>> output_L(N);

	/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/
	ThreadPool pool(num_of_threads);
	std::vector<std::future<int>> results; // return typename: xxx
	for (int v_k = 0; v_k < N; v_k++)
	{
		results.emplace_back(
			pool.enqueue([&output_L, v_k] 
			{ // pass const type value j to thread; [] can be empty
				sort(L_temp_599[v_k].begin(), L_temp_599[v_k].end(), compare_hop_constrained_two_hop_label);
				vector<hop_constrained_two_hop_label>(L_temp_599[v_k]).swap(L_temp_599[v_k]); // 使用vector的swap优化内存占用，释放多余的空间
				output_L[v_k] = L_temp_599[v_k];
				vector<hop_constrained_two_hop_label>().swap(L_temp_599[v_k]); // clear new labels for RAM efficiency

				return 1; // return to results; the return type must be the same with results
			}));
	}
	for (auto &&result : results)
		result.get(); // all threads finish here

	return output_L;
}

/*canonical_repair*/
void hop_constrained_clean_L(hop_constrained_case_info &case_info, int thread_num)
{

	auto &L = case_info.L;
	int N = L.size();
	label_size_before_canonical_repair_599 = 0;
	label_size_after_canonical_repair_599 = 0;

	ThreadPool pool(thread_num);
	std::vector<std::future<int>> results;
	/* test the correctness of async */
	// vector<int> list;
	// list.push_back(4);
	// list.push_back(5);
	// list.push_back(0);
	// list.push_back(1);
	// list.push_back(2);
	// list.push_back(3); 
	// for (int v = 0; v < N; v++)
	// for (int index = 0; index < N; index++)
	for(int v=0;v<N;v++)
	{
		results.emplace_back(
			pool.enqueue([v, &L] { // pass const type value j to thread; [] can be empty
				mtx_599[max_N_ID_for_mtx_599 - 1].lock();
				int used_id = Qid_599.front();
				Qid_599.pop();
				mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

				vector<hop_constrained_two_hop_label> Lv_final;

				/**
				 * get the L result of the current vertex
				*/
				mtx_599[v].lock_shared();
				vector<hop_constrained_two_hop_label> Lv = L[v];
				mtx_599[v].unlock_shared();
				label_size_before_canonical_repair_599 += Lv.size();

				/**
				 * the temp_L in this thread
				*/
				auto &T = Temp_L_vk_599[used_id];

				/**
				 * Traverse the L-list of the current vertex
				*/
				for (auto Lvi : Lv)
				{
					int u = Lvi.hub_vertex;
					int u_hop = Lvi.hop;

					/**
					 * Traverse the L on the opposite vertex of the current label.
					*/
					mtx_599[u].lock_shared();
					auto Lu = L[u];
					mtx_599[u].unlock_shared();

					/**
					 * traverse downward from the perfectly correct first vertex
					*/
					int min_dis = std::numeric_limits<int>::max();
					for (auto &label1 : Lu)
					{
						for (auto &label2 : T[label1.hub_vertex])
						{
							if (label1.hop + label2.second <= u_hop)
							{
								long long int query_dis = label1.distance + (long long int)label2.first;
								if (query_dis < min_dis)
								{
									min_dis = query_dis;
								}
							}
						}
					}

					if (min_dis > Lvi.distance)
					{
						Lv_final.push_back(Lvi);
						T[u].push_back({Lvi.distance, Lvi.hop});
					}
				}

				for (auto label : Lv_final)
				{
					vector<pair<int, int>>().swap(T[label.hub_vertex]);
				}

				mtx_599[v].lock();
				L[v] = Lv_final;
				L[v].shrink_to_fit();
				mtx_599[v].unlock();
				label_size_after_canonical_repair_599 += Lv_final.size();

				mtx_599[max_N_ID_for_mtx_599 - 1].lock();
				Qid_599.push(used_id);
				mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

				return 1; // return to results; the return type must be the same with results
			}));
	}

	for (auto &&result : results)
		result.get(); // all threads finish here
	results.clear();

	case_info.label_size_before_canonical_repair = label_size_before_canonical_repair_599;
	case_info.label_size_after_canonical_repair = label_size_after_canonical_repair_599;
	case_info.canonical_repair_remove_label_ratio = (double)(label_size_before_canonical_repair_599 - label_size_after_canonical_repair_599) / label_size_before_canonical_repair_599;
}

void hop_constrained_two_hop_labels_generation(graph_v_of_v<int> &input_graph, hop_constrained_case_info &case_info)
{

	//----------------------------------- step 1: initialization -----------------------------------
	/* generate containers to store L labels and PPR, and create a thread pool */
	auto begin = std::chrono::high_resolution_clock::now();
	/*  */
	labal_size_599 = 0;
	begin_time_599 = std::chrono::high_resolution_clock::now();
	/* max run time in nanoseconds. 100*1e9 in this program*/
	max_run_time_nanoseconds_599 = case_info.max_run_time_seconds * 1e9;
	// TODO why the max_bit_size is 6e9
	max_labal_size_599 = case_info.max_bit_size / sizeof(hop_constrained_two_hop_label);

	int N = input_graph.size();
	/* store the L Label and PPR*/
	L_temp_599.resize(N);
	PPR_599.resize(N);
	if (N > max_N_ID_for_mtx_599)
	{
		cout << "N > max_N_ID_for_mtx_599!" << endl;
		exit(1);
	}

	int num_of_threads = case_info.thread_num;
	ThreadPool pool(num_of_threads);
	std::vector<std::future<int>> results;

	ideal_graph_599 = input_graph;
	global_use_2M_prune = case_info.use_2M_prune;
	global_use_rank_prune = case_info.use_rank_prune;

	auto end = std::chrono::high_resolution_clock::now();
	case_info.time_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;

	//----------------------------------------------- step 2: generate labels ---------------------------------------------------------------
	/** 
	 * Use Temp_L_vk_599 to mark the positional relationship between the iterated nodes 
	 * and the nodes with already generated labels. This is done to reduce the process of
	 * traversing the L labels of the iterated nodes. Additionally, register multithreaded
	 * tasks and retrieve the results 
	 */ 
	begin = std::chrono::high_resolution_clock::now();

	global_upper_k = case_info.upper_k == 0 ? std::numeric_limits<int>::max() : case_info.upper_k;

	Temp_L_vk_599.resize(num_of_threads);
	dist_hop_599.resize(num_of_threads);
	Q_handle_priorities_599.resize(num_of_threads);
	Vh_599.resize(num_of_threads);
	hop_constrained_node_handle handle_x;
	for (int i = 0; i < num_of_threads; i++)
	{
		Temp_L_vk_599[i].resize(N);
		dist_hop_599[i].resize(N, {std::numeric_limits<int>::max(), 0});
		Q_handle_priorities_599[i].resize(N);
		for (int j = 0; j < N; j++)
		{
			Q_handle_priorities_599[i][j].resize(global_upper_k + 1, {handle_x, std::numeric_limits<int>::max()});
		}
		Vh_599[i].resize(global_upper_k + 2);
		Qid_599.push(i);
	}
	if (case_info.use_2023WWW_generation)
	{
		for (int v_k = 0; v_k < N; v_k++)
		{
			results.emplace_back(
				pool.enqueue([v_k]
							 {
					_2023WWW_thread_function(v_k);
					return 1; }));
		}
	}
	else
	{
		int last_check_vID = N - 1;
		if (global_use_2M_prune)
		{
			for (int v_k = N - 1; v_k >= 0; v_k--)
			{
				if (is_mock[v_k])
				{
					last_check_vID = v_k;
					break;
				}
			}
		}

		for (int v_k = 0; v_k <= last_check_vID; v_k++)
		{
			// if (global_use_2M_prune && is_mock[v_k]) {
			//	continue;
			// }
			results.emplace_back(
				pool.enqueue([v_k]
							 {
					HSDL_thread_function(v_k);
					return 1; }));
		}
	}
	for (auto &&result : results)
		result.get();

	end = std::chrono::high_resolution_clock::now();
	case_info.time_generate_labels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	//----------------------------------------------- step 3: sortL---------------------------------------------------------------
	begin = std::chrono::high_resolution_clock::now();
	case_info.L = L_temp_599;
	case_info.L = hop_constrained_sortL(num_of_threads);
	case_info.PPR = PPR_599;

	end = std::chrono::high_resolution_clock::now();
	case_info.time_sortL = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	//----------------------------------------------- step 4: canonical_repair---------------------------------------------------------------
	begin = std::chrono::high_resolution_clock::now();
	if (case_info.use_canonical_repair)
	{
		hop_constrained_clean_L(case_info, num_of_threads);
	}

	end = std::chrono::high_resolution_clock::now();
	case_info.time_canonical_repair = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	//---------------------------------------------------------------------------------------------------------------------------------------

	case_info.time_total = case_info.time_initialization + case_info.time_generate_labels + case_info.time_sortL + case_info.time_canonical_repair;

	hop_constrained_clear_global_values();
}