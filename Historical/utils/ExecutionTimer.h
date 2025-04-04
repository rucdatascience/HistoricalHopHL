#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <memory>

namespace experiment
{

	struct Subtask
	{
		std::string name;
		std::chrono::steady_clock::time_point startTime;
		std::chrono::steady_clock::time_point endTime;
		std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> subtasks;
		std::weak_ptr<Subtask> parent;
		bool isEnded = false;

		Subtask() : startTime(std::chrono::steady_clock::now())
		{
			name = "task";
		}

		Subtask(const std::string &taskName, std::weak_ptr<Subtask> parentTask = std::weak_ptr<Subtask>())
			: name(taskName), startTime(std::chrono::steady_clock::now()), parent(parentTask)
		{
		}

		void end()
		{
			endTime = std::chrono::steady_clock::now();
			isEnded = true;
		}

		double getDuration() const
		{
			if (!isEnded)
				return 0.0;
			return std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count();
		}

		double getTotalDuration() const
		{
			if (!subtasks.empty())
			{
				double totalDuration = 0.0;
				for (const auto &subtask : subtasks)
				{
					totalDuration += subtask.second->getTotalDuration();
				}
				return totalDuration;
			}
			return getDuration();
		}
	};

	class ExecutionTimer
	{
	public:
		void startTask(const std::string &taskName)
		{
			currentMainTaskName = taskName;
			currentSubtask = nullptr;
			allTasks.clear();
		}

		void startSubtask(const std::string &subtaskName)
		{
			auto newSubtask = std::make_shared<Subtask>(subtaskName, currentSubtask);

			if (currentSubtask)
			{
				currentSubtask->subtasks.emplace_back(subtaskName, newSubtask);
			}
			else
			{
				allTasks.emplace_back(subtaskName, newSubtask);
			}
			currentSubtask = newSubtask;
		}

		long long int endSubtask()
		{
			if (currentSubtask)
			{
				currentSubtask->end();
				long long res = currentSubtask->getTotalDuration();
				auto parentTask = currentSubtask->parent.lock();
				currentSubtask = parentTask ? parentTask : nullptr;
				return res;
			}
			else
			{
				std::cerr << "Error: No subtask is currently running.\n";
				return 0;
			}
		}

		double getTaskDuration()
		{
			double totalDuration = 0;
			for (const auto &task : allTasks)
			{
				totalDuration += task.second->getTotalDuration();
			}
			return totalDuration;
		}

		void printStats()
		{
			std::cout << "Task Timing Statistics: " << getTaskDuration() << " seconds\n";
			printTaskStats(currentMainTaskName, allTasks, 0);
		}

		void writeStatsToFile(std::ofstream &outFile)
		{
			if (!outFile.is_open())
			{
				std::cerr << "Error opening file for writing." << std::endl;
				return;
			}

			outFile.precision(6);
			outFile.setf(std::ios::fixed);
			outFile.setf(std::ios::showpoint);
			outFile << "Task Timing Statistics:" << std::endl;
			writeStatsToFileRecursively(outFile, currentMainTaskName, allTasks, 0);
		}

	private:
		void writeStatsToFileRecursively(std::ofstream &outFile, const std::string &taskName, const std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> &tasks, int level)
		{
			for (const auto &task : tasks)
			{
				outFile << std::string(level * 4, ' ') << "-> " << task.first << ", Time: "
						<< task.second->getTotalDuration() << " seconds" << std::endl;
				writeStatsToFileRecursively(outFile, task.first, task.second->subtasks, level + 1);
			}
		}

		void printTaskStats(const std::string &taskName, const std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> &tasks, int level)
		{
			for (const auto &task : tasks)
			{
				std::cout << std::string(level * 4, ' ') << "-> " << task.first << ": "
						  << std::fixed << std::setprecision(3) << task.second->getTotalDuration() << "s" << std::endl;
				printTaskStats(task.first, task.second->subtasks, level + 1);
			}
		}

	private:
		std::string currentMainTaskName;
		std::shared_ptr<Subtask> currentSubtask;
		std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> allTasks;
	};

}
