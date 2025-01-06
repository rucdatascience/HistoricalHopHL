#include "Historical/graph_v_of_v/graph_v_of_v_with_time_span_experiment_reader.h"

int main()
{
    experiment_config experiment(
        "D:\\project\\postgraduate\\graph\\HistoricalHopHL\\Historical\\experiment\\Email-Enron.txt",
        "D:\\project\\postgraduate\\graph\\HistoricalHopHL\\Historical\\experiment\\Email-Enron-saving.txt",
        10,
        30,
        false,
        200, 1);
        experiment.init();
        experiment.process();
        experiment.close();
}
