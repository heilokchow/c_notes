#include <iostream>
#include <array>
#include <vector>
#include <thread>
#include <mutex>
#include "omp.h"

#include "cross_vaidation.h"
#include "timer.h"

static long num_steps = 100000000;
double step;
#define NUM_THREADS 8

int main()
{
    auto stat = std::make_shared<Linear_regression>(10000, 500);
    
    omp_set_num_threads(NUM_THREADS);


    {
        Timer t;
#pragma omp parallel for schedule (dynamic)
    for (int i = 0; i < 100; i++) {
        MatrixXd X = MatrixXd::Random(2000, 2000);
        MatrixXd Y = MatrixXd::Random(2000, 2000);
        MatrixXd C;
        C.noalias() = X * Y;
        printf("pos: %d\n", i);
        }
    }

    /*OMP guide code
    step = 1.0 / static_cast<double>(num_steps);
    double pi = 0.0, sum = 0.0;
    omp_set_num_threads(NUM_THREADS);

    double start = omp_get_wtime();
#pragma omp parallel
    {
        int i = omp_get_thread_num(), j;
        int k = omp_get_num_threads();
        printf("%d, %d\n", i, k);
        double x, sum_par = 0.0;
        for ( j = i; j < num_steps; j = j + k) {
            x = (j + 0.5) * step;
            sum_par += 4.0 / (1.0 + x * x);
        }
#pragma omp critical
        sum += sum_par;
    }
    double end = omp_get_wtime();

    pi = sum * step;

    printf("%f, Pi: %.6f", end - start, pi);*/

    /*Use std::thread for multuithreading
    auto stat = std::make_shared<Linear_regression>(100,10);
  
    std::array<std::thread, 8> workers;

    {
        Timer t;
        for (auto& t : workers)
            t = std::thread(cross_validation_bench, stat);
        
        for (auto& t : workers)
            t.join();
    }*/

    return 0;
}