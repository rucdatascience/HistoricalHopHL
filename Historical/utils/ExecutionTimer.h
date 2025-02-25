#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <memory>

namespace experiment {

	struct Subtask {
		std::string name;
		std::chrono::steady_clock::time_point startTime;
		std::chrono::steady_clock::time_point endTime;
		std::map<std::string, std::shared_ptr<Subtask>> subtasks;
		std::weak_ptr<Subtask> parent;

		Subtask() : startTime(std::chrono::steady_clock::now()) {
			name = "task";
		}

		Subtask(const std::string& taskName, std::weak_ptr<Subtask> parentTask = std::weak_ptr<Subtask>())
			: name(taskName), startTime(std::chrono::steady_clock::now()), parent(parentTask) {
		}

		void end() {
			endTime = std::chrono::steady_clock::now();
		}

		double getDuration() const {
			return std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count();
		}

		double getTotalDuration() const {
			double totalDuration = getDuration();
			for (const auto& subtask : subtasks) {
				totalDuration += subtask.second->getTotalDuration();
			}
			return totalDuration;
		}
	};

	class ExecutionTimer {
	public:
		void startTask(const std::string& taskName) {
			// Start a new main task, do not start timer
			currentMainTaskName = taskName;
			currentSubtask = nullptr; // Reset current subtask
			allTasks.clear();
		}

		void startSubtask(const std::string& subtaskName) {
			if (currentSubtask) {
				// Adding a subtask to the current task
				auto& currentTask = *currentSubtask;
				currentTask.subtasks[subtaskName] = std::make_shared<Subtask>(subtaskName, currentSubtask);
				currentSubtask = currentTask.subtasks[subtaskName]; // Set the new subtask as the current one
			}
			else {
				// If no current subtask is running, start a new subtask under the current task
				allTasks[subtaskName] = std::make_shared<Subtask>(subtaskName);
				currentSubtask = allTasks[subtaskName]; // Set it as the current task
			}
		}

		void endSubtask() {
			if (currentSubtask) {
				// End the current subtask
				currentSubtask->end();

				// Move to parent task (if any)
				auto parentTask = currentSubtask->parent.lock();
				if (parentTask) {
					currentSubtask = parentTask;
				}
				else {
					currentSubtask = nullptr;
				}
			}
			else {
				std::cerr << "Error: No subtask is currently running.\n";
			}
		}

		double getTaskDuration() {
			double totalDuration = 0;
			for (const auto& task : allTasks) {
				totalDuration += task.second->getTotalDuration();
			}
			return totalDuration;
		}

		void printStats() {
			std::cout << "Task Timing Statistics:" << std::endl;
			printTaskStats(currentMainTaskName, allTasks, 0);
		}

		void writeStatsToFile(std::ofstream& outFile) {
			if (!outFile.is_open()) {
				std::cerr << "Error opening file for writing: " << std::endl;
				return;
			}

			outFile.precision(6);
			outFile.setf(std::ios::fixed);
			outFile.setf(std::ios::showpoint);
			outFile << "Task Timing Statistics:" << std::endl;
			writeStatsToFileRecursively(outFile, currentMainTaskName, allTasks, 0);
		}

		void writeStatsToFileRecursively(std::ofstream& outFile, const std::string& taskName, const std::map<std::string, std::shared_ptr<Subtask>>& tasks, int level) {
			for (const auto& task : tasks) {
				outFile << std::string(level * 4, ' ') << "-> " << task.first << ", Time: " << task.second->getDuration() << " seconds" << std::endl;
				writeStatsToFileRecursively(outFile, task.first, task.second->subtasks, level + 1);
			}
		}

	private:
		void printTaskStats(const std::string& taskName, const std::map<std::string, std::shared_ptr<Subtask>>& tasks, int level) {
			for (const auto& task : tasks) {

				std::cout << std::string(level * 4, ' ') << "-> " << task.first << ": "
					<< std::fixed << std::setprecision(3) << task.second->getDuration() << "s" << std::endl;
				printTaskStats(task.first, task.second->subtasks, level + 1);
			}
		}

	private:
		std::string currentMainTaskName;
		std::shared_ptr<Subtask> currentSubtask;
		std::map<std::string, std::shared_ptr<Subtask>> allTasks;
	};
}
