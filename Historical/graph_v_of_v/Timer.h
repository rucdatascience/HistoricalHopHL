#pragma once
#include <chrono>
#include <string>
#include <iostream>

class Timer
{
public:
    bool is_debug = false;
    std::chrono::_V2::system_clock::time_point start_time;

    void begin_timing()
    {
        start_time = std::chrono::high_resolution_clock::now();
    }
    void mark_time(std::string current_step)
    {
        if (is_debug)
        {
            auto now = std::chrono::high_resolution_clock::now();
            double runtime_base_line_with_span = std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_time).count() / 1e9;
            std::cout << current_step << ": " << runtime_base_line_with_span << std::endl;
        }
    }
};
static Timer timer;
