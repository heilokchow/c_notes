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

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#ifndef N
#define N 20
#endif
#ifndef FS
#define FS 25
#endif

struct node {
    int data;
    int fibdata;
    struct node* next;
};

int fib(int n) {
    int x, y;
    if (n < 2) {
        return (n);
    }
    else {
        x = fib(n - 1);
        y = fib(n - 2);
        return (x + y);
    }
}

void processwork(struct node* p)
{
    int n;
    n = p->data;
    p->fibdata = fib(n);
}

struct node* init_list(struct node* p) {
    int i;
    struct node* head = NULL;
    struct node* temp = NULL;

    head = (node*)malloc(sizeof(struct node));
    p = head;
    p->data = FS;
    p->fibdata = 0;
    for (i = 0; i < N; i++) {
        temp = (node*)malloc(sizeof(struct node));
        p->next = temp;
        p = temp;
        p->data = FS + i + 1;
        p->fibdata = i + 1;
    }
    p->next = NULL;
    return head;
}

int main(int argc, char* argv[]) {
    double start, end;
    struct node* p = NULL;
    struct node* temp = NULL;
    struct node* head = NULL;

    printf("Process linked list\n");
    printf("  Each linked list node will be processed by function 'processwork()'\n");
    printf("  Each ll node will compute %d fibonacci numbers beginning with %d\n", N, FS);

    p = init_list(p);
    head = p;
    
    start = omp_get_wtime();
#pragma omp parallel
#pragma omp single
    {
        for (temp = head; temp; temp = temp->next)
#pragma omp task firstprivate(temp)
            processwork(temp);
    }
    end = omp_get_wtime();
    p = head;
    while (p != NULL) {
        printf("%d : %d\n", p->data, p->fibdata);
        temp = p->next;
        free(p);
        p = temp;
    }
    free(p);

    printf("Compute Time: %f seconds\n", end - start);

    printf("Hello World\n");

    return 0;
}


#endif