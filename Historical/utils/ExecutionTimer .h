#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <chrono>
#include <map>
#include <iomanip>
#include <sstream>
namespace experiment {

	struct Subtask {
		std::string name;
		std::chrono::steady_clock::time_point startTime;
		std::chrono::steady_clock::time_point endTime;
		Subtask() : startTime(std::chrono::steady_clock::now()) {
			name = "task";
		}
		Subtask(const std::string& taskName)
			: name(taskName), startTime(std::chrono::steady_clock::now()) {
		}

		void end() {
			endTime = std::chrono::steady_clock::now();
		}

		double getDuration() const {
			return std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count();
		}
	};

	class ExecutionTimer {
	public:
		void startTask(const std::string& taskName) {
			// Start a new main task, do not start timer
			currentMainTaskName = taskName;
			currentSubtask.clear();
			std::map<std::string, Subtask>().swap(allTasks);
		}

		void startSubtask(const std::string& subtaskName) {
			if (currentSubtask.empty()) {
				// If no current subtask is running, start a new subtask under the current task
				currentSubtask = subtaskName;
				if (allTasks.find(subtaskName) == allTasks.end()) {
					allTasks[subtaskName] = Subtask(subtaskName);
				}
				else {
					std::cerr << "Error: any subtask has the same name for " << subtaskName << std::endl;
				}
			}
			else {
				std::cerr << "Error: subtask is already running under the current task." << std::endl;
			}
		}

		void endSubtask() {
			if (!currentSubtask.empty()) {
				// End the current subtask
				allTasks[currentSubtask].end();
				currentSubtask.clear();
			}
			else {
				std::cerr << "Error: No subtask is currently running.\n";
			}
		}

		double getTaskDuration() {
			double totalDuration = 0;
			for (const auto& subtask : allTasks) {
				totalDuration += subtask.second.getDuration();
			}
			return totalDuration;
		}

		void printStats() {
			for (const auto& task : allTasks) {
				printTaskStats(task.first);
			}
		}
		void writeStatsToFile(std::string& filename) {
			std::ofstream outFile;
			outFile.precision(6);
			outFile.setf(std::ios::fixed);
			outFile.setf(std::ios::showpoint);
			outFile.open(filename);
			if (!outFile.is_open()) {
				std::cerr << "Error opening file for writing: " << filename << std::endl;
				return;
			}
			outFile << "Task Timing Statistics:" << std::endl;
			outFile << "Task: " << currentMainTaskName << ", Total Time: " << this->getTaskDuration() << " seconds" << std::endl;
			for (const auto& task : allTasks) {
				outFile << "    Subtask: " << task.first << ", Time: " << std::to_string(task.second.getDuration()) << " seconds" << std::endl;
			}

			outFile.close();
		}

	private:
		void printTaskStats(const std::string& taskName) {
			std::cout << "Main Task: " << taskName << "\n";
			for (const auto& subtask : allTasks) {
				const auto& sub = subtask.second;
				std::cout << "  Subtask: " << subtask.first << " - "
					<< std::fixed << std::setprecision(3) << sub.getDuration() << "s\n";
			}
		}

	private:
		std::string currentMainTaskName;
		std::string currentSubtask;
		std::map<std::string, Subtask> allTasks;
	};
}