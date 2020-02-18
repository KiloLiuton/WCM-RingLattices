#include <iostream>
#include <omp.h>
#include <functional>
#include "dynamics.hpp"
#define maxEle(x) ((x[0]>x[1])&&(x[0]>x[2])) ? 0 : ((x[1]>x[0])&&(x[1]>x[2])) ? 1 : ((x[2]>x[0])&&(x[2]>x[1])) ? 2 : -1

void customInitStates(struct trial &t) {
    initWave(t, 6, false);
    t.states[14] = 1;
    t.states[15] = 2;
    t.states[16] = 0;
}

int main()
{
    pcg32 RNG(23, 42);
    struct trial trial;
    std::cout << "MAXN=" << MAXN << std::endl;
    std::cout << "MAXK=" << MAXK << std::endl;
    std::cout << "sizeof(trial)[MB]=" << sizeof(trial)/1e3 << std::endl;

    std::function<void(struct trial &)> initStates = NULL;
    initStates = customInitStates;
    int N = 91;
    int K = 2;
    initTrial(trial, N, K, 2.8, RNG, NULL, initStates);

    for (int i=0; i<N; i++) {
        printf("%d", trial.states[i]);
    }
    printf("\n");
    double psi = psiOP(trial);
    int cycl = cycles(trial, 7);
    printf("psi=%f  cycles=%d\n", psi, cycl);

    return 0;
}
