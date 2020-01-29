#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>

class Timer
{
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    double duration;
    
public:
    Timer() 
    {
        start = std::chrono::system_clock::now();
    }
    ~Timer() 
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        duration = elapsed.count();
        std::cout << "Elapsed Time: " << duration << "s\n";
    }
};

#endif