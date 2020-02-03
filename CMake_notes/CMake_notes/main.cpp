#ifdef MAIN

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
    /*OMP PI Example
    double pi = 0.0, sum = 0.0;
    double step = 1.0 / static_cast<double>(num_steps);

    double start = omp_get_wtime();

    #pragma omp parallel
    {
        double x;
        #pragma omp for reduction (+:sum)
            for (int i = 0; i < num_steps; i++) {
                x = (i + 0.5) * step;
                sum += 4.0 / (1.0 + x * x);
            }
    }

    double end = omp_get_wtime();

    pi = sum * step;
    printf("%f, Pi: %.6f\n", end - start, pi);*/

    auto stat = std::make_shared<Linear_regression>(10000, 1000);
    
    omp_set_num_threads(NUM_THREADS);
    int flag = 1;

    {
        Timer t;
        #pragma omp parallel for schedule (dynamic)
        for (int i = 0; i < 100; i++) {
            int id = omp_get_thread_num();
            cross_validation_bench_omp(stat, i, id);
        }
    }

    /*OMP guide code Module 1
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

#else

#include <omp.h>
#ifdef APPLE
#include <stdlib.h>
#else
#include <memory>
#endif
#include <stdio.h>

#define N        10000000

/* Some random number constants from numerical recipies */
#define SEED       2531
#define RAND_MULT  1366
#define RAND_ADD   150889
#define RAND_MOD   714025
int randy = SEED;

/* function to fill an array with random numbers */
void fill_rand(int length, double* a)
{
    int i;
    for (i = 0; i < length; i++) {
        randy = (RAND_MULT * randy + RAND_ADD) % RAND_MOD;
        *(a + i) = ((double)randy) / ((double)RAND_MOD);
    }
}

/* function to sum the elements of an array */
double Sum_array(int length, double* a)
{
    int i;  double sum = 0.0;
    for (i = 0; i < length; i++)  sum += *(a + i);
    return sum;
}

int main()
{
    double* A, sum, runtime;
    int flag = 0;

    A = (double*)malloc(N * sizeof(double));

    runtime = omp_get_wtime();

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            fill_rand(N, A);        // Producer: fill an array of data
            #pragma omp flush
            flag = 1;
            #pragma omp flush (flag)
        }

        #pragma omp section
        {
            #pragma omp flush (flag)
            while (flag == 0) {
                #pragma omp flush (flag)
            }
            #pragma omp flush
            sum = Sum_array(N, A);  // Consumer: sum the array
        }
    }

    runtime = omp_get_wtime() - runtime;

    printf(" In %f seconds, The sum is %f \n", runtime, sum);
}

#endif