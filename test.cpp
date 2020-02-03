#include <iostream>
#include <omp.h>
#include "dynamics.hpp"


int main()
{
    pcg32 RNG(23, 42);
    struct trial trial;
    initTrial(trial, 10, 2, 2.0, RNG, NULL, NULL);
    std::cout << "MAXN=" << MAXN << std::endl;
    std::cout << "MAXK=" << MAXK << std::endl;
    std::cout << "sizeof(trial)[MB]=" << sizeof(trial)/1e3 << std::endl;

    int na = 20;
    int trials = 100;
    omp_set_num_threads(3);
#pragma omp parallel default(none) shared(na,trials)
    {
#pragma omp single
        {
            printf("Number of threads: %d\n", omp_get_num_threads());
        }
        printf("Hi from thread %d\n", omp_get_thread_num());
#pragma omp for
        for (int i=0; i<na; i++) {
            printf("  For: thread %d\n", omp_get_thread_num());
        }
    }

    return 0;
}
