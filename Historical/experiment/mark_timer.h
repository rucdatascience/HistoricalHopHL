#pragma once
#include <vector>
class mark_timer
{
private:
    // 0 is init
    // 1 is the slot time1
    // .....
    std::vector<double> experiment_time;
    std::chrono::_V2::system_clock::time_point time1;
    double time_cost = 0;
    int addTimeMark(double time)
    {
        time_cost = time_cost + (time);
        return 0;
    }
    int pushTimeMarkToVector(double time)
    {
        this->experiment_time.push_back(time);
        return 0;
    }

public:
    int mark()
    {
        this->time1 = std::chrono::high_resolution_clock::now();
        return 0;
    }

    int add()
    {
        std::chrono::_V2::system_clock::time_point time2 = std::chrono::high_resolution_clock::now();
        addTimeMark(std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() / 1e9);
        return 0;
    }

    int push()
    {
        pushTimeMarkToVector(this->time_cost);
        this->time_cost = 0;
        return 0;
    }

    std::vector<double> get_experiment_time()
    {
        return this->experiment_time;
    }

    mark_timer(){
    };
    ~mark_timer() {

    };
};
