# HistoricalHopHL
## How to use HistoricalHopHL
<b>Normally, there is an example in each header file to explain how to use the codes in this file. If there is not, you can add a mark, e.g., "/\*an example is needed\*/", in the end of this file, and then pull a request. I will add an example ASAP.</b>
You can download the whole repository into a rucgraph folder, and add the path of this folder when compiling cpp codes.
### An example of using rucgraph on a Windows
> 1. Download the whole HistoricalHopHL folder
> 2. Given a cpp file named as "try_generate_random_graph.cpp", the contents in which are
> ```
> using namespace std;
> 
> #include "GPU/graph_v_of_v/graph_v_of_v_generate_random_graph_test.h"
> 
> int main()
> {
> 	/**
> 		Test the algorithm for generating undirected weighted graphs
>	*/
> 	test_generate();
> }
> ```
> ,where "graph_v_of_v_generate_random_graph_test.h" is a header file in HistoricalHopHL. In the terminal, compile and run the above codes using the following commands:
> ```
> $ g++ -std=c++17 -g ${workspaceFolder}\try**.cpp -o ${workspaceFolder}\try***.exe -I {boost_lnclude_path} -I {std_lnclude_path} -I {workspaceFolder}
> ```
> , where "-I {workspaceFolder}" 、 "-I{boost_lnclude_path}" and  "-I{std_lnclude_path}" are to add the path of the workfolder and stdlib and boost folder when compiling. Then, you successfully run an example application that generates a weighted undirected graph.