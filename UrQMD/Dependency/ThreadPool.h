#ifndef THREADPOOL_H
#define THREADPOOL_H
#define MULTITHREAD 1
#include <iostream>
#include <thread>
#include <vector>
#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>

class ThreadPool
{
public:
    static ThreadPool &getInstance()
    {
        static ThreadPool instance(std::thread::hardware_concurrency());
        return instance;
    }
    ~ThreadPool()
    {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }

        condition.notify_all();

        for (auto &worker : workers)
        {
            worker.join(); // Wait for all threads to complete
        }
    }
    void enqueue(std::function<void()> task)
    {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            tasks.push(std::move(task));
        }
        condition.notify_one(); // Notify one thread to start the task
    }

    void waitForCompletion()
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        condition.wait(lock, [this]()
                       { return tasks.empty() && activeTasks == 0; });
    }
    std::mutex queueMutex;
private:
    ThreadPool(size_t numThreads) : stop(false), activeTasks(0)
    {
        for (size_t i = 0; i < numThreads; ++i)
        {
            workers.push_back(std::thread([this]()
                                          {
            while (true) {
                std::function<void()> task;

                {
                    std::unique_lock<std::mutex> lock(queueMutex);
                    condition.wait(lock, [this]() { return stop || !tasks.empty(); });

                    if (stop && tasks.empty())
                        return;

                    task = std::move(tasks.front());
                    tasks.pop();
                    ++activeTasks;  // Increment active task counter
                }

                task();  // Execute the task

                {
                    std::unique_lock<std::mutex> lock(queueMutex);
                    --activeTasks;  // Decrement active task counter
                    if (tasks.empty() && activeTasks == 0) {
                        // Notify the main thread if all tasks are completed
                        condition.notify_all();
                    }
                }
            } }));
        }
    }
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::condition_variable condition;
    bool stop;
    size_t activeTasks;
};
#endif // THREADPOOL_H