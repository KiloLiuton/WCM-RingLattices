#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP
#include <bits/stdint-uintn.h>
#include <cstdlib>
#include <sys/types.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>
#include <functional>
#include "pcg_random.hpp"

#ifndef MAXN
#define MAXN 25001
#endif
#define MAXK MAXN/2 - (1 - MAXN % 2)
#define MAX_TABLE_SIZE 4*MAXK+1
#define SIN_PHI1 0.8660254037844387
#define maxEle(x) ((x[0]>x[1])&&(x[0]>x[2])) ? 0 : ((x[1]>x[0])&&(x[1]>x[2])) ? 1 : ((x[2]>x[0])&&(x[2]>x[1])) ? 2 : -1

struct trial {
    int N, N0, N1, N2;
    int K;
    double a;
    double t, dt;
    pcg32 rng;
    double rates_sum;
    std::uniform_real_distribution<double> uniform{0.0, 1.0};
    double  nat_freqs[MAXN]       = {0};
    uint8_t states[MAXN]          = {0};
    int32_t deltas[MAXN]          = {0};
    double  rates[MAXN]           = {0};
    double  table[MAX_TABLE_SIZE] = {0};
};

inline double getRate(int delta, struct trial &trial)
{
    return trial.table[delta + 2*trial.K];
}


void initTable(struct trial &trial)
{
    /*Requires trial.K and trial.a to be initialized*/
    for (int i=0; i<4*trial.K+1; i++) {
        trial.table[i] = exp(trial.a*(i-2*trial.K)/(2*trial.K));
    }
}

int getDelta(int i, struct trial &trial)
{
    uint8_t statei = trial.states[i];
    uint8_t statej;
    int delta = 0;
    for (int k=1; k<=trial.K; k++) {
        // CCW neighbor
        int nbr_index = ((i+k) >= trial.N) ? i+k-trial.N : i+k;
        statej = trial.states[nbr_index];
        if (statej == statei)            delta--;
        else if (statej == (statei+1)%3) delta++;

        // CW neighbor
        nbr_index = ((i-k) < 0) ? i-k+trial.N : i-k;
        statej = trial.states[nbr_index];
        if (statej == statei)            delta--;
        else if (statej == (statei+1)%3) delta++;
    }
    return delta;
}

void initDeltas(struct trial &trial)
{
    /*Requires trial.states and trial.N to be initialized*/
    for (int i=0; i<trial.N; i++) trial.deltas[i] = getDelta(i, trial);
}

void initRates(struct trial &trial)
{
    /*Requires trial.deltas and trial.N to be initialized*/
    trial.rates_sum = 0;
    for (int i=0; i<trial.N; i++) {
        double rate = trial.nat_freqs[i] * getRate(trial.deltas[i], trial);
        trial.rates[i] = rate;
        trial.rates_sum += rate;
    }
}

void initGaussianNaturalFreqs(
        struct trial &trial,
        double mean,
        double stddev
        )
{
    /*Requires trial.rng and trial.N to be initialized*/
    std::uniform_real_distribution<double> Normal(mean, stddev);
    for (int i=0; i<trial.N; i++) {
        trial.nat_freqs[i] = Normal(trial.rng);
    }
}

void initUniform(struct trial &trial, uint8_t phase)
{
    /*Initialize the system with all phase values to phase.
     * Requires trial.N to be initialized*/
    for (int i=0; i<trial.N; i++) trial.states[i] = phase;
}

void initWave(struct trial &trial, int num_waves, bool reversed)
{
    /*Initialize the system with `num_waves` domains of constant phase*/
    int tmp = trial.N / num_waves;
    for (int i=0; i<trial.N; i++) {
        uint8_t s = (i / tmp) % 3;
        if (reversed) {
            trial.states[trial.N - 1 - i] = s;
        } else {
            trial.states[i] = s;
        }
    }
}

void initPops(struct trial &trial)
{
    /*Requires trial.states and trial.N to be initialized*/
    trial.N0 = trial.N1 = trial.N2 = 0;
    for (int i=0; i<trial.N; i++) {
        uint8_t j = trial.states[i];
        switch (j) {
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
}

void initTrial(
        struct trial &trial,
        int size,
        int range,
        double coupling,
        pcg32 &RNG,
        std::function<void(struct trial &)> initNaturalFreqs=NULL,
        std::function<void(struct trial &)> initStates=NULL
        )
{
    if (initNaturalFreqs == NULL) {
        initNaturalFreqs = [](struct trial &trial){
            initGaussianNaturalFreqs(trial, 1.0, 0.1);
        };
    }
    if (initStates == NULL) {
        initStates = [](struct trial &trial){ initUniform(trial, 0); };
    }
    trial.N = size;
    trial.K = range;
    trial.a = coupling,
    trial.t = 0;
    trial.dt = 0;
    trial.rng = RNG;
    initTable(trial);
    initStates(trial);
    initPops(trial);
    initNaturalFreqs(trial);
    initDeltas(trial);
    initRates(trial);
}

inline int sampleRate(struct trial &trial)
{
    bool searching = true;
    int index = trial.N-1;
    double thold = trial.uniform(trial.rng) * trial.rates_sum;
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

inline void updateSite(struct trial &trial, int i)
{
    uint8_t prev_state = trial.states[i];
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
    trial.states[i] = (trial.states[i]+1) % 3;
    trial.deltas[i] = 0;
    for (int k=1; k<=trial.K; k++) {
        // get indexes for CCW and CW neighbors and update
        int nbrs[2] = {
            (i+k)>=trial.N ? i+k-trial.N : i+k,
            (i-k)<0 ? i-k+trial.N : i-k
        };
        for (auto j : nbrs) {
            if (trial.states[j] == trial.states[i]) {
                trial.deltas[i] -= 1;
                trial.deltas[j] -= 1;
            } else if (trial.states[j] == prev_state) {
                trial.deltas[j] += 2;
            } else {
                trial.deltas[i] += 1;
                trial.deltas[j] -= 1;
            }
            double newrate = trial.nat_freqs[j] * getRate(trial.deltas[j], trial);
            trial.rates_sum += newrate - trial.rates[j];
            trial.rates[j] = newrate;
        }
    }
    double newrate = trial.nat_freqs[i] * getRate(trial.deltas[i], trial);
    trial.rates_sum += newrate - trial.rates[i];
    trial.rates[i] = newrate;
    trial.dt = 1./trial.rates_sum;
    trial.t += trial.dt;
}

inline int update(struct trial &trial)
{
    int index = sampleRate(trial);
    updateSite(trial, index);
    return index;
}

inline double kuramotoOP(struct trial const &trial)
{
    return (double) sqrt(
            trial.N0*trial.N0
            + trial.N1*trial.N1
            + trial.N2*trial.N2
            - trial.N1*trial.N2
            - trial.N0*trial.N1
            - trial.N0*trial.N2) / trial.N;
}

inline void psiOPandCycles(struct trial const &trial, double result[2], int const window_size)
{
    double sum1 = 0, sum2 = 0;
    uint8_t curr_state;
    int pops[3] = {0, 0, 0};
    int curr_window_state, prev_window_state = -1;
    int count = 0;
    int cycles = 0;
    for (int i = 0; i < trial.N; i++) {
        curr_state = trial.states[i];
        switch (curr_state) {
       case 0:
            sum1 += trial.rates[i];
            break;
        case 1:
            sum1 += -0.5 * trial.rates[i];
            sum2 += SIN_PHI1 * trial.rates[i];
            break;
        case 2:
            sum1 += -0.5 * trial.rates[i];
            sum2 += -SIN_PHI1 * trial.rates[i];
            break;
        default:
            std::cout << (int) curr_state << " <- Unknown state in psiOP!\n";
            break;
        }

        if ( window_size > 0 ) {
        // to get the cycles we take the median of the previous `window_size` phase values in order to avoid small fluctuations.
        // `window_size` should be set to one sixth of the size of the minimum wave size (that is, half the size of third of a wave)
            pops[curr_state]++;
            count++;
            if ( (count == window_size) || (i == (trial.N-1)) ) {
                curr_window_state = maxEle(pops);
                // if a window is undefined (but not the first window) repeat the state of the previous window
                if ( curr_window_state == -1 ) {
                    if ( prev_window_state == -1 ) curr_window_state = prev_window_state = 0;
                    curr_window_state = prev_window_state;
                }
                if ( (prev_window_state != curr_window_state) ) {
                    if ( prev_window_state==0 ) {
                        if      ( curr_window_state==1 ) cycles++;
                        else if ( curr_window_state==2 ) cycles--;
                    } else if ( prev_window_state==1 ) {
                        if      ( curr_window_state==2 ) cycles++;
                        else if ( curr_window_state==0 ) cycles--;
                    } else if ( prev_window_state==2 ) {
                        if      ( curr_window_state==0 ) cycles++;
                        else if ( curr_window_state==1 ) cycles--;
                    }
                }
                prev_window_state = curr_window_state;
                count   = 0;
                pops[0] = 0;
                pops[1] = 0;
                pops[2] = 0;
            }
        }
    }
    result[0] = sqrt(pow(sum1, 2.0) + pow(sum2, 2.0)) / trial.N;
    result[1] = cycles;
}

inline int cycles(struct trial const &trial, int window_size)
{
    int pops[3] = {0, 0, 0};
    int curr_state;
    int curr_window_state;
    int prev_window_state = -1;
    int count = 0;
    int cycles = 0;
    for (int i=0; i<trial.N; i++) {
        curr_state = trial.states[i];
        pops[curr_state]++;
        count++;
        if ( (count == window_size) || (i == (trial.N-1)) ) {
            curr_window_state = maxEle(pops);
            // if a window is undefined (but not the first window) repeat the state of the previous window
            if ( curr_window_state == -1 ) {
                if ( prev_window_state == -1 ) curr_window_state = prev_window_state = 0;
                curr_window_state = prev_window_state;
            }
            if ( (prev_window_state != curr_window_state) ) {
                if ( prev_window_state==0 ) {
                    if      ( curr_window_state==1 ) cycles++;
                    else if ( curr_window_state==2 ) cycles--;
                } else if ( prev_window_state==1 ) {
                    if      ( curr_window_state==2 ) cycles++;
                    else if ( curr_window_state==0 ) cycles--;
                } else if ( prev_window_state==2 ) {
                    if      ( curr_window_state==0 ) cycles++;
                    else if ( curr_window_state==1 ) cycles--;
                }
            }
            prev_window_state = curr_window_state;
            count   = 0;
            pops[0] = 0;
            pops[1] = 0;
            pops[2] = 0;
        }
    }
    return cycles;

    //int s1, s2;
    //int pops[3];
    //int count = 0;
    //int cycles = 0;
    //for (int i=0; i<trial.N; i++) {
    //    s1 = trial.states[i];
    //    int j = (i < (trial.N-1)) ? i+1 : 0;
    //    s2 = trial.states[j];
    //    pops[s2]++;
    //    count++;
    //    if (count == window_size) {
    //        s2 = maxEle(pops);
    //        if ( s1!=s2 ) {
    //            if        ( s1==0 ) {
    //                if      ( s2==1 ) cycles++;
    //                else if ( s2==2 ) cycles--;
    //            } else if ( s1==1 ) {
    //                if      ( s2==2 ) cycles++;
    //                else if ( s2==0 ) cycles--;
    //            } else if ( s1==2 ) {
    //                if      ( s2==0 ) cycles++;
    //                else if ( s2==1 ) cycles--;
    //            }
    //        }

    //        count = 0;
    //        pops[0]  = 0;
    //        pops[1]  = 0;
    //        pops[2]  = 0;
    //    }
    //}
    //return cycles;
}

#endif
