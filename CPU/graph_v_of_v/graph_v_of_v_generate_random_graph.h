#pragma once
#include <boost/random.hpp>
/*this function generates a random graph with vertex and edge weights, and this graph
may not be connected.*/
#include "CPU/graph_v_of_v/graph_v_of_v.h"

template <typename weight_type>
graph_v_of_v<weight_type> graph_v_of_v_generate_random_graph(long long int V, long long int E, double ec_min, double ec_max,
															 int input_precision, boost::random::mt19937 &boost_random_time_seed)
{ // must use long long int for large V E

	/*time complexity: O(|V||E|)*/

	/* randomly generate the weight of edge */
	double precision = std::pow(10, input_precision);
	boost::random::uniform_int_distribution<> dist_ec{static_cast<int>(ec_min * precision), static_cast<int>(ec_max * precision)};

	/*time complexity: O(|V|)*/
	graph_v_of_v<weight_type> random_graph(V); // generate vertices

	/*add edges to random_graph*/
	long long int max_E = V * (V - 1) / 2; // must use long long int for large V
	/* If this graph is a complete graph  */
	if (E == max_E)
	{ // complete graphs
		/*time complexity: O(|V|+|E|)*/
		for (int i = 0; i < V; i++)
		{
			for (int j = 0; j < i; j++)
			{
				weight_type new_cost = (weight_type)dist_ec(boost_random_time_seed) / precision; // generate ec
				random_graph.add_edge(i, j, new_cost);
			}
		}
	}
	/* the number of e is a wrong value*/
	else if (E > max_E)
	{
		std::cout << "E: " << E << std::endl;
		std::cout << "V * (V - 1) / 2: " << max_E << std::endl;
		std::cout << "E > V * (V - 1) / 2 in graph_v_of_v_generate_random_graph!" << '\n';
		exit(1);
	}
	else
	{ // incomplete graphs

		/*time complexity: O(|V|)*/
		/* init a vector to store the set of vertices that may potentially have edges added */
		std::vector<int> not_full_vertices; // vertices without a full degree
		for (int i = 0; i < V; i++)
		{
			not_full_vertices.push_back(i);
		}

		/*time complexity: O(|V||E|)*/
		int edge_num = 0;
		while (edge_num < E)
		{
			boost::random::uniform_int_distribution<> dist_id{static_cast<int>(0), static_cast<int>(not_full_vertices.size() - 1)};
			int RAND = dist_id(boost_random_time_seed); // generate int random number  0, not_full_vertices.size()-1
			if (random_graph.degree(not_full_vertices[RAND]) < V - 1)
			{ // randomly select a vertex and this is a vertex without a full degree

				/*time complexity: O(|V|)*/
				/* generate an increasing vector<int> from zero */
				std::vector<int> unchecked(V);
				std::iota(std::begin(unchecked), std::end(unchecked), 0);
				bool added = false;
				while (added == false)
				{
					/* set the range of random numbers for edges */
					boost::random::uniform_int_distribution<> dist_id2{static_cast<int>(0), static_cast<int>(unchecked.size() - 1)};
					/* randomly select a edge */
					int x = dist_id2(boost_random_time_seed);
					int j = unchecked[x];
					/* if this edge does not point to itself (RAND) and does not exist in the graph*/
					if (not_full_vertices[RAND] != j &&
						random_graph.contain_edge(not_full_vertices[RAND], j) == 0)
					{
						// This edge does not exist
						/* calculate the weight */
						weight_type new_cost = (weight_type)dist_ec(boost_random_time_seed) / precision; // generate ec
						/* add the edge with weight to graph*/
						random_graph.add_edge(not_full_vertices[RAND], j, new_cost); // add a new edge
						edge_num++;
						added = true;
						break; // break after adding one edge
					}
					else
					{
						unchecked.erase(unchecked.begin() + x);
					}
				}
			}
			else
			{ // this is a vertex with a full degree
				not_full_vertices.erase(not_full_vertices.begin() + RAND);
			}
		}
	}

	return random_graph;
}