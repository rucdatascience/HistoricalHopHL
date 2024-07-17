#include "graph_v_of_v/graph_v_of_v.h"
#include "graph_v_of_v/graph_v_of_v_generate_random_graph.h"
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

void test_generate(){
    int V = 100, E = 500;
    int ec_min = 1, ec_max = 20;
    graph_v_of_v<int> instance_graph;
    instance_graph = graph_v_of_v_generate_random_graph<int>(V, E, ec_min, ec_max, 1, boost_random_time_seed);
    instance_graph.print();
}