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
    printf("%f, Pi: %.6f", end - start, pi);

    /*OMP guide code Module 5
    auto stat = std::make_shared<Linear_regression>(4000, 500);
    
    omp_set_num_threads(NUM_THREADS);


    {
        Timer t;
#pragma omp parallel for schedule (dynamic)
        for (int i = 0; i < 100; i++) {
            cross_validation_bench_omp(stat, i);
        }
    }*/

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

# define NPOINTS 1000
# define MAXITER 1000


struct d_complex {
    double r;
    double i;
};

void testpoint(d_complex&);

struct d_complex c;
int numoutside = 0;
omp_lock_t _lock;

int main() {
    int i, j;
    double area, error, eps = 1.0e-5;
    omp_init_lock(&_lock);

    //   Loop over grid of points in the complex plane which contains the Mandelbrot set,
    //   testing each point to see whether it is inside or outside the set.

#pragma omp parallel for default(shared) private(c) firstprivate(eps)
    for (i = 0; i < NPOINTS; i++) {
        for (j = 0; j < NPOINTS; j++) {
            c.r = -2.0 + 2.5 * (double)(i) / (double)(NPOINTS)+eps;
            c.i = 1.125 * (double)(j) / (double)(NPOINTS)+eps;
            testpoint(c);
        }
    }

    omp_destroy_lock(&_lock);

    // Calculate area of set and error estimate and output the results
    
    area = 2.0 * 2.5 * 1.125 * (double)(NPOINTS * NPOINTS - numoutside) / (double)(NPOINTS * NPOINTS);
    error = area / (double)NPOINTS;

    printf("Area of Mandlebrot set = %12.8f +/- %12.8f, %d\n", area, error, numoutside);
    printf("Correct answer should be around 1.510659\n");

}

void testpoint(d_complex& c) {

    // Does the iteration z=z*z+c, until |z| > 2 when point is known to be outside set
    // If loop count reaches MAXITER, point is considered to be inside the set

    struct d_complex z;
    int iter;
    double temp;

    z = c;
    printf("%d, %d\n", z.i, z.r);
    for (iter = 0; iter < MAXITER; iter++) {
        temp = (z.r * z.r) - (z.i * z.i) + c.r;
        z.i = z.r * z.i * 2 + c.i;
        z.r = temp;
        if ((z.r * z.r + z.i * z.i) > 4.0) {
            omp_set_lock(&_lock);
            numoutside++;
            omp_unset_lock(&_lock);
            break;
        }
    }

}

#endif