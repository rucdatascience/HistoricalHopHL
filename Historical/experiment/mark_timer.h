#pragma once
#include <vector>
#include <chrono>

class mark_timer {
private:
	std::vector<double> experiment_time;
	std::chrono::steady_clock::time_point time1;
	double time_cost = 0;

	void addTimeMark(double time) {
		time_cost += time;
	}

	void pushTimeMarkToVector(double time) {
		experiment_time.push_back(time);
	}

public:
	void mark() {
		time1 = std::chrono::steady_clock::now();
	}

	void add() {
		auto time2 = std::chrono::steady_clock::now();
		addTimeMark(std::chrono::duration<double>(time2 - time1).count());
	}

	void push() {
		pushTimeMarkToVector(time_cost);
		time_cost = 0;
	}

	std::vector<double> get_experiment_time() const {
		return experiment_time;
	}

	mark_timer() = default;
	~mark_timer() = default;
};
