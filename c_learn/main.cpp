#include <iostream>
#include <array>
#include <vector>
#include <thread>
#include <mutex>

#include "cross_vaidation.h"
#include "timer.h"


int main()
{
    auto stat = std::make_shared<Linear_regression>();
  
    std::array<std::thread, 8> workers;

    {
        Timer t;
        for (auto& t : workers)
            t = std::thread(cross_validation_bench, stat);
        
        for (auto& t : workers)
            t.join();
    }
    return 0;
}