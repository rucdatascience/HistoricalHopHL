#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span_experiment_reader.h"

int main()
{
    try
    {
        experiment_config experiment(
            "E:\\project\\postgraduate\\HistoricalHopHL\\Historical\\experiment\\Email-Enron.txt",
            "E:\\project\\postgraduate\\HistoricalHopHL\\Historical\\experiment\\Email-Enron-saving.txt",
            10,
            30,
            true,
            200, 1);
        experiment.init();
        experiment.process();
        experiment.print_experiment_result();
        experiment.close();
    }
    catch (const char *error)
    {
        std::cout << std::string(error)<<"test" << std::endl;
    }
}
