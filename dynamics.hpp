#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP
#include <bits/stdint-uintn.h>
#include <sys/types.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>
#include "pcg_random.hpp"
#define MAXN 25001
#define MAXK 12500
#define MAX_TABLE_SIZE 4*MAXK+1
#define SIN_PHI1 0.8660254037844387

struct trial {
    int N, N0, N1, N2;
    int K;
    pcg32 rng;
    std::uniform_real_distribution<double> uniform{0.0, 1.0};
    uint8_t states[MAXN] = {0};
    int32_t deltas[MAXN] = {0};
    double rates[MAXN] = {0};
    double rates_sum;
    double t;
    double table[MAX_TABLE_SIZE] = {0};
};

void initTable(struct trial &trial, double a)
{
    for (int i=0; i<4*trial.K+1; i++) {
        trial.table[i] = exp(a*(i-2*trial.K)/(2*trial.K));
    }
}

double getRate(struct trial &trial, int delta)
{
    return trial.table[delta + 2*trial.K];
}

int getDelta(int i, struct trial &trial)
{
    uint8_t statei = trial.states[i];
    uint8_t statej;
    int delta = 0;
    for (int k=1; k<=trial.K; k++) {
        statej = trial.states[((i+k) >= trial.N) ? i+k-trial.N : i+k];
        if (statej == statei)            delta--;
        else if (statej == (statei+1)%3) delta++;

        statej = trial.states[((i-k) < 0) ? i-k+trial.N : i-k];
        if (statej == statei)            delta--;
        else if (statej == (statei+1)%3) delta++;
    }
    return delta;
}

void initUniform(
        struct trial &trial,
        int size,
        int range,
        double a,
        uint32_t seed,
        uint32_t stream
        )
{
    /*Initialize the system with all phase values to 0*/
    trial.N = size;
    trial.N0 = trial.N, trial.N1 = 0, trial.N2 = 0;
    trial.K = range;
    trial.rates_sum = 0;
    trial.rng.seed(seed, stream);
    initTable(trial, a);
    for (int i=0; i<trial.N; i++) {
        trial.states[i] = 0;
        trial.deltas[i] = -2*trial.K;
        double rate = getRate(trial, -2*trial.K);
        trial.rates[i] = rate;
        trial.rates_sum += rate;
    }
    trial.t = 0;
}

void initWave(
        struct trial &trial,
        int size,
        int range,
        double a,
        uint32_t seed,
        uint32_t stream,
        int num_waves
        )
{
    /*Initialize the system with `num_waves` domains of constant phase*/
    trial.N = size;
    trial.N0 = trial.N1 = trial.N2 = 0;
    trial.K = range;
    trial.rates_sum = 0;
    trial.rng.seed(seed, stream);
    initTable(trial, a);
    int tmp = size / num_waves;
    for (int i=0; i<size; i++) {
        uint8_t s = (i / tmp) % 3;
        trial.states[i] = s;
        switch (s) {
            case 0:
                trial.N0++;
                break;
            case 1:
                trial.N1++;
                break;
            case 2:
                trial.N2++;
                break;
        }
    }
    for (int i=0; i<size; i++) {
        trial.deltas[i] = getDelta(i, trial);
        trial.rates[i] = getRate(trial, trial.deltas[i]);
    }
    trial.t = 0;
}

inline int sampleRate(struct trial &trial)
{
    bool searching = true;
    int index = trial.N-1;
    double tmp = trial.uniform(trial.rng);
    double thold = tmp*trial.rates_sum;
    double partial_sum = 0;
    for (int i=0; i<trial.N; i++) {
        partial_sum += trial.rates[i];
        if (searching && (partial_sum > thold)) {
            index = i;
            searching = false;
        }
    }
    trial.rates_sum = partial_sum;
    return index;
}

inline void updateSite(struct trial &trial, int index)
{
    trial.t += 1./trial.rates_sum;
    uint8_t prev_state = trial.states[index];
    switch (prev_state) {
        case 0:
            trial.N0--;
            trial.N1++;
            break;
        case 1:
            trial.N1--;
            trial.N2++;
            break;
        case 2:
            trial.N2--;
            trial.N0++;
            break;
        default:
            printf("Invalid state found in updateSite\n");
            exit(EXIT_FAILURE);
    }
    trial.states[index] = (trial.states[index]+1) % 3;
    trial.deltas[index] = 0;
    for (int i=1; i<=trial.K; i++) {
        // update CW and CCW neighbors
        int nbrs[2] = {
            (index+i)>=trial.N ? index+i-trial.N : index+i,
            (index-i)<0 ? index-i+trial.N : index-i
        };
        for (int j=0; j<2; j++) {
            if (trial.states[nbrs[j]] == trial.states[index]) {
                trial.deltas[index] -= 1;
                trial.deltas[nbrs[j]] -= 1;
            } else if (trial.states[nbrs[j]] == prev_state) {
                trial.deltas[nbrs[j]] += 2;
            } else {
                trial.deltas[index] += 1;
                trial.deltas[nbrs[j]] -= 1;
            }
            trial.rates_sum -= trial.rates[nbrs[j]];
            trial.rates[nbrs[j]] = getRate(trial, trial.deltas[nbrs[j]]);
            trial.rates_sum += trial.rates[nbrs[j]];
        }
    }
    trial.rates_sum -= trial.rates[index];
    trial.rates[index] = getRate(trial, trial.deltas[index]);
    trial.rates_sum += trial.rates[index];
}

inline int update(struct trial &trial)
{
    int index = sampleRate(trial);
    updateSite(trial, index);
    return index;
}

double kuramotoOP(struct trial const &trial) {
    return (double) sqrt(
            trial.N0*trial.N0
            + trial.N1*trial.N1
            + trial.N2*trial.N2
            - trial.N1*trial.N2
            - trial.N0*trial.N1
            - trial.N0*trial.N2) / trial.N;
}

double psiOP(struct trial const &trial) {
    double s1 = 0, s2 = 0;
    for (int i = 0; i < trial.N; i++) {
        uint8_t s = trial.states[i];
        switch (s) {
       case 0:
            s1 += trial.rates[i];
            break;
        case 1:
            s1 += -0.5 * trial.rates[i];
            s2 += SIN_PHI1 * trial.rates[i];
            break;
        case 2:
            s1 += -0.5 * trial.rates[i];
            s2 += -SIN_PHI1 * trial.rates[i];
            break;
        default:
            std::cout << (int) s << " <- Unknown state in psiOP!\n";
            break;
        }
    }
    return sqrt(std::pow(s1, 2.0) + std::pow(s2, 2.0)) / trial.N;
}

#endif
