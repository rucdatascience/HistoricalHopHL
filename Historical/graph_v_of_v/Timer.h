#pragma once
#include <chrono>
#include <string>
#include <iostream>

class Timer
{
public:
	bool is_debug = false;
	std::chrono::steady_clock::time_point start_time;

	void begin_timing()
	{
		start_time = std::chrono::steady_clock::now();
	}

	void mark_time(const std::string& current_step)
	{
		if (is_debug)
		{
			auto now = std::chrono::steady_clock::now();
			double elapsed_time = std::chrono::duration<double>(now - start_time).count();
			std::cout << current_step << ": " << elapsed_time << "s" << std::endl;
		}
	}
};

static Timer timer;
