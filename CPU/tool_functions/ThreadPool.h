/*this is from https://github.com/progschj/ThreadPool 

an explanation: https://www.cnblogs.com/chenleideblog/p/12915534.html
*/


#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool {
public:
    ThreadPool(size_t); // �̳߳صĹ��캯��
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> std::future<typename std::result_of<F(Args...)>::type>; // ���������ӵ��̳߳ص����������
    ~ThreadPool(); // �̳߳ص���������
private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers; // ���ڴ���̵߳����飬��vector��������
    // the task queue
    std::queue< std::function<void()> > tasks; // ���ڴ������Ķ��У���queue���н��б��档��������Ϊstd::function<void()>����Ϊstd::function��ͨ�ö�̬������װ������������������д�ŵ���һ��������
    
    // synchronization
    std::mutex queue_mutex; // һ������������еĻ��������ڲ�����������߳�ȡ��������Ҫ�������������а�ȫ����
    std::condition_variable condition; // һ������֪ͨ�߳��������״̬����������������������֪ͨ�߳̿���ִ�У��������wait״̬
    bool stop; // ��ʶ�̳߳ص�״̬�����ڹ����������ж��̳߳�״̬���˽�
};
 
// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads) // ���캯������Ϊinline�����ղ���threads��ʾ�̳߳���Ҫ�������ٸ��̡߳�
    :   stop(false) // stop��ʼΪfalse������ʾ�̳߳������š�
{
    for(size_t i = 0;i<threads;++i) // ����forѭ�������δ���threads���̣߳��������߳�����workers�С�
        workers.emplace_back( 
            /*
            ��vector�У�emplace_back()��Ա������������������β������һ����������Ч����push_back()һ����������������΢���죬
            ��emplace_back(args)�з���Ķ���Ĳ�������push_back(OBJ(args))�з�����Ƕ��󡣼�emplace_back()ֱ�����������Դ���Ĳ���ֱ�ӵ��ö���Ĺ��캯�������µĶ���
            ��push_back()���ȵ��ö���Ĺ��캯������һ����ʱ�����ٽ���ʱ���󿽱��������ڴ��С�
            
            lambda����ʽ�ĸ�ʽΪ��
            [ ���� ] ( �β� ) ˵����(��ѡ) �쳣˵�� attr -> �������� { ������ }
            ��������lambda����ʽΪ [ ���� ] { ������ } ���͡������lambda����ʽ���Ǵ����̵߳�ִ�к���.       
            */
            [this] // ��lambda����ʽ�����̳߳�ָ��this�����ں�������ʹ�ã������̳߳س�Ա����stop��tasks�ȣ�    
            {
                for(;;) // for(;;)Ϊһ����ѭ������ʾÿ���̶߳��ᷴ������ִ�У�����ʵÿ���̳߳��е��̶߳���������
                {
                    std::function<void()> task; // ��ѭ���У��ȴ���һ����װvoid()������std::function����task�����ڽ��պ�������������е�������ʵ����

                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex); // ��{}�ڣ�queue_mutex������״̬

                        /*�����ʾ���̳߳���ֹͣ������������в�Ϊ�գ��򲻻���뵽wait״̬��
                          ���ڸտ�ʼ�����̳߳أ��̳߳ر�ʾδֹͣ�����������Ϊ�գ�����ÿ���̶߳�����뵽wait״̬��
                          ������������������֪ͨ���߳̾ͻ�������½���
                        */
                        this->condition.wait(lock,
                            [this]{ return this->stop || !this->tasks.empty(); });

                        /*���̳߳��Ѿ�ֹͣ���������Ϊ�գ���return���������߳�������ѭ�������׳��̳߳أ��ڴ�֮ǰ��ÿ���߳�Ҫô��wait״̬��Ҫô��ִ�������task*/
                        if(this->stop && this->tasks.empty())
                            return;

                        /*
                        ����������еĵ�һ��������task��ǣ�Ȼ����������и����񵯳������˴��߳�ʵ�ڻ������������еĻ�����������½��еģ����������������̺߳�
                        �߳���wait�����ڵõ���������еĻ������Ż��������ִ�С���������ֻ����һ���߳��õ����񣬲��ᷢ����ȺЧӦ��
                        ���˳���{ }�����Ƕ�������е����ӵ���Ҳ�ͷ��ˣ�Ȼ�����ǵ��߳̾Ϳ���ִ�������õ�������task�ˣ�ִ�����֮���߳��ֽ�������ѭ����
                        */
                        task = std::move(this->tasks.front());
                        this->tasks.pop(); // task����ѻ���queue��
                    }

                    task();
                }
            }
        );
}

// add new work item to the pool
/*
equeue��һ��ģ�庯�����������β�ΪF��Args������class... Args��ʾ��������βΡ�
auto�����Զ��Ƶ���equeue�ķ������ͣ��������β�Ϊ(F&& f, Args&&... args)������&&��ʾ��ֵ���á���ʾ����һ��F���͵�f�������ɸ�Args���͵�args��

typename std::result_of<F(Args...)>::type   //�����ArgsΪ������F�ĺ������͵ķ�������
std::future<typename std::result_of<F(Args...)>::type> //std::future���������첽�����Ľ��
���շ��ص��Ƿ���std::future�е�F(Args��)�������͵��첽ִ�н����
*/
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) 
    -> std::future<typename std::result_of<F(Args...)>::type> // ��ʾ�������ͣ���lambda����ʽ�еı�ʾ����һ����
{
    using return_type = typename std::result_of<F(Args...)>::type; // �����ArgsΪ������F�ĺ������͵ķ�������

    auto task = std::make_shared< std::packaged_task<return_type()> >(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        
    std::future<return_type> res = task->get_future(); // res�б���������Ϊreturn_type�ı�������task�첽ִ����ϲſ��Խ�ֵ�����ȥ
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace([task](){ (*task)(); });
    }
    condition.notify_one(); // //�������������к���Ҫȥ����һ���߳�
    return res; // //���߳�ִ����ϣ���ִ�еĽ������
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool()
{
    /*
    �����������У��ȶ���������м�������ֹͣ�������Ϊtrue������������ʹ���µĲ����������Ҳ��ִ��ʧ��
    */
    {
        std::unique_lock<std::mutex> lock(queue_mutex); // ��{}�ڣ�queue_mutex������״̬
        stop = true;
    }
    /*ʹ�������������������̣߳������̶߳�������ִ��*/
    condition.notify_all();

    /*��stop����Ϊtrue�����������Ϊ��ʱ����Ӧ���߳̽�������ѭ������
    ��ÿ���߳�����Ϊjoin���ȵ�ÿ���߳̽�����Ϻ����߳����˳���*/
    for(std::thread &worker: workers)
        worker.join();
}

#endif











/*Example:

#include <iostream>
#include <GPU/tool_functions/ThreadPool.h>
using namespace std;

class example_class
{
public:
    int a;
    double b;
};
example_class example_function(example_class x)
{
    return x;
}

void ThreadPool_example()
{
    ThreadPool pool(4);	// ����һ���̳߳أ������߳�Ϊ4
    std::vector<std::future<example_class>> results; // return typename: example_class; ������߳�ִ�н��
    for (int i = 0; i < 10; ++i)
    { // ����10������
        int j = i + 10;
        results.emplace_back(  // ����ÿ���첽���
            pool.enqueue([j] { // ��ÿ��������뵽��������У�lambda����ʽ�� pass const type value j to thread; [] can be empty
                example_class a;
                a.a = j;
                return example_function(a); // return to results; the return type must be the same with results
            }));
    }
    for (auto &&result : results)			// һ��ȡ��������results�еĽ��
        std::cout << result.get().a << ' '; // result.get() makes sure this thread has been finished here;
    results.clear(); // future get value can only be called once. get��results�����future�Ͳ������ˣ�clear֮��results���ܱ��ٴ�ʹ��
    std::cout << std::endl;
    // if result.get() is pair<int, string>, then you cannot use result.get().first = result.get().second
}

int main()
{
    ThreadPool_example();
}

*/



/*Example: multi_thread write a vector (use std lock mechanism)  


#include <iostream>
#include <GPU/tool_functions/ThreadPool.h>
#include <mutex>
using namespace std;

int vector_size = 3;
vector<vector<int>> vectors(vector_size);
vector<std::mutex> mtx(vector_size); // ����vectors
void thread_function(int ID, int value)
{
    mtx[ID].lock(); // only one thread can lock mtx[ID] here, until mtx[ID] is unlocked
    vectors[ID].push_back(value);
    mtx[ID].unlock();
}
void ThreadPool_example()
{
    ThreadPool pool(5);	// use 5 threads
    std::vector<std::future<int>> results; // return typename: xxx
    for (int i = 0; i < vector_size; ++i)
    {
        for (int j = 0; j < 1e1; j++)
        {
            results.emplace_back(
                pool.enqueue([i, j] { // pass const type value j to thread; [] can be empty
                    thread_function(i, j);
                    return 1; // return to results; the return type must be the same with results
                }));
        }
    }
    for (auto &&result : results)
        result.get(); // all threads finish here
    for (int i = 0; i < vector_size; ++i)
    {
        cout << "vectors[" << i << "].size(): " << vectors[i].size() << endl;
        for (int x = 0; x < vectors[i].size(); x++)
        {
            std::cout << vectors[i][x] << ' ';
        }
        std::cout << std::endl;
    }
}
int main(){ThreadPool_example();}


*/
