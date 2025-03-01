#include "ExecutionTimer.h"
#include <thread>

int main()
{
    experiment::ExecutionTimer timer;
    timer.startTask("Main Task");

    timer.startSubtask("Subtask 1");
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    timer.startSubtask("Subtask 1.1");
    std::this_thread::sleep_for(std::chrono::milliseconds(300));
    timer.endSubtask(); // 结束 Subtask 1.1
    timer.startSubtask("Subtask 1.2");
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    timer.endSubtask(); // 结束 Subtask 1.2
    timer.endSubtask(); // 结束 Subtask 1

    timer.startSubtask("Subtask 2");
    std::this_thread::sleep_for(std::chrono::milliseconds(400));
    timer.startSubtask("Subtask 2.1");
    std::this_thread::sleep_for(std::chrono::milliseconds(300));
    timer.endSubtask(); // 结束 Subtask 1.1
    timer.startSubtask("Subtask .2");
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    timer.endSubtask(); // 结束 Subtask 1.2
    timer.endSubtask(); // 结束 Subtask 2

    timer.printStats();

    return 0;
}