#include <bits/types/struct_timespec.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <getopt.h>
#include <ctime>
#include <map>
#include <memory>
#include <string>
#include <fstream>
#include <functional>
#include "dynamics.hpp"
#include "logging.hpp"

pcg32 medianRNG;
int median(int n0, int n1, int n2)
{
    // return the phase value of the majority, selecting randomly on ties.
    if (n0 > n2 && n0 > n1)       return 0;
    if (n1 > n2 && n1 > n0)       return 1;
    if (n2 > n0 && n2 > n1)       return 2;
    if ((n0 == n1) && (n1 == n2)) return pcg_extras::bounded_rand(medianRNG, 3);
    if (n0 == n1)                 return pcg_extras::bounded_rand(medianRNG, 2);
    if (n1 == n2)                 return pcg_extras::bounded_rand(medianRNG, 2) + 1;
    if (n0 == n2)                 return pcg_extras::bounded_rand(medianRNG, 2) ? 0 : 2;
    
    std::cout << "oops, this shouldn't have happened\n";
    return 0;
}

void compressAndLogPhases(FILE *file, struct trial &trial, int window_size)
{
    if (window_size > 1) {
        int compress_counter, n0, n1, n2;
        compress_counter = 0;
        n0 = n1 = n2 = 0;
        for (int j=0; j<trial.N; j++) {
            if (trial.states[j] == 0)      n0++;
            else if (trial.states[j] == 1) n1++;
            else if (trial.states[j] == 2) n2++;
            compress_counter++;
            if (compress_counter == window_size) {
                fprintf(file, "%d,", median(n0, n1, n2));
                n0 = n1 = n2 = 0;
                compress_counter = 0;
            }
        }
        if (compress_counter > 0) {
            fprintf(file, "%d,%f\n", median(n0, n1, n2), trial.t);
        } else {
            fprintf(file, "%f\n", trial.t);
        }
    } else if (window_size == 1) {
        for (int j=0; j<trial.N; j++) {
            fprintf(file, "%d,", trial.states[j]);
        }
        fprintf(file, "%f\n", trial.t);
    }
}

int main(int argc, char *argv[])
{
    int optN = 200;
    int optK = 20;
    double opta = 1.8;
    double gmean = 1.0;
    double gstddev = 0.1;
    int iters = 5000;
    int seed = 23;
    int stream = 42;
    int log_phases = 0;
    int log_both = 0;
    int log_interval = 1;
    std::string filename = "";
    std::string ic = "wave2";
    int num_waves;

    char opt;
    int opt_idx = 0;
    std::string sopt = "n:k:a:i:s:r:o:c:";
    static struct option lopt[] = {
        {"size",          required_argument, 0, 'n'},
        {"range",         required_argument, 0, 'k'},
        {"coupling",      required_argument, 0, 'a'},
        {"gmean",         required_argument, 0, '4'},
        {"gstddev",       required_argument, 0, '5'},
        {"iters",         required_argument, 0, 'i'},
        {"log-phases",    required_argument, 0,  1 },
        {"log-both",      required_argument, 0,  2 },
        {"log-interval",  required_argument, 0,  3 },
        {"random-seed",   required_argument, 0, 's'},
        {"random-stream", required_argument, 0, 'r'},
        {"outfile",       required_argument, 0, 'o'},
        {"ic",            required_argument, 0, 'c'},
        {0,               0,                 0,  0 }
    };
    while ((opt=getopt_long(argc, argv, sopt.c_str(), lopt, &opt_idx)) != -1) {
        switch (opt) {
        case 'n':
            if ((size_t) atoi(optarg) > MAXN) {
                printf("  -n option cannot exceed MAXN=%d\n", MAXN);
                exit(EXIT_FAILURE);
            }
            optN = atoi(optarg);
            break;
        case 'k':
            if ((size_t) atoi(optarg) > MAXK) {
                printf("  -k option cannot exceed MAXN=%d\n", MAXK);
                exit(EXIT_FAILURE);
            }
            optK = atoi(optarg);
            break;
        case 'a':
            opta = atof(optarg);
            break;
        case 1:
            log_phases = abs(atoi(optarg));
            break;
        case 2:
            log_both = abs(atoi(optarg));
            break;
        case 3:
            log_interval = abs(atoi(optarg));
            if (log_interval < 1) log_interval = 1;
            break;
        case 4:
            gmean = atof(optarg);
            break;
        case 5:
            gstddev = atof(optarg);
            break;
        case 'i':
            iters = atoi(optarg);
            break;
        case 's':
            seed = atoi(optarg);
            break;
        case 'r':
            stream = atoi(optarg);
            break;
        case 'o':
            filename = std::string(optarg);
            break;
        case 'c':
            ic = std::string(optarg);
            break;
        }
    }
    if (optK > MAXK) {
        printf("-k --range option cannot exceed %d\n", MAXK);
        exit(EXIT_FAILURE);
    } else if (optK > optN/2) {
        printf("-k --range option cannot exceed N/2\n");
        exit(EXIT_FAILURE);
    }
    if (optN > MAXN) {
        printf("-n --size option cannot exceed N/2\n");
        exit(EXIT_FAILURE);
    }
    if (filename.empty()) {
        const char* prefix = log_phases ? "phasetrial" : "trial";
        char suffix[16];
        sprintf(suffix, "a-%.6f", opta);
        filename = defaultFilePath(optN, optK, prefix, suffix);
    }
    std::ostringstream metadata;
    metadata << "N=" << optN << " K=" << optK << " a=" << opta
             << " gmean=" << gmean << " gstddev=" << gstddev
             << " ic=" << ic.c_str() << " log-interval=" << log_interval
             << " log-phases=" << log_phases
             << " iters=" << iters << " seed=" << seed
             << " stream=" << stream << "\n";

    // Initialize system with command line options
    struct trial trial;
    pcg32 RNG(seed, stream);
    std::function<void(struct trial &)> initStates = NULL;
    std::function<void(struct trial &)> initNaturalFreqs = NULL;
    if (ic == "uniform") {
        printf("\n\nINITIALIZING UNIFORM\n\n");
        initStates = [](struct trial &t){
            initUniform(t, 0);
        };
        num_waves = 0;
    } else if (ic.substr(0, 4) == "wave") {
        printf("\n\nINITIALIZING WAVES\n\n");
        num_waves = stoi(ic.substr(4));
        initStates = [num_waves](struct trial &t){
            initWave(t, num_waves, false);
        };
    } else {
        std::cout << "Initial configuration " << ic << " not understood!\n"
            << "Try one of:\n"
            << "uniform\n"
            << "wave[n] - where n is the number of waves\n";
        exit(EXIT_FAILURE);
    }

    initNaturalFreqs = [gmean, gstddev](struct trial &t) {
        initGaussianNaturalFreqs(t, gmean, gstddev);
    };
    initTrial(trial, optN, optK, opta, RNG, initNaturalFreqs, initStates);
    printf("Trial initialized: %d %d %d\n", trial.N0, trial.N1, trial.N2);

    printf("Memory usage: %.2f MB\n", sizeof(trial)/1e3);
    printf("Initializing with:\n"
           "%s"
           "log-phases %d\n"
           "log-interval %d\n"
           "log-both %d\n"
           "ic %s\n",
           metadata.str().c_str(),
           log_phases,
           log_interval,
           log_both,
           ic.c_str()
           );

    
    std::map<int, int> hist; // count states that transitioned

    int interval = 0;
    struct timespec start, finish;
    if (log_both > 0) {
        char suffix[16];
        sprintf(suffix, "a-%.6f", opta);
        std::string fn1 = defaultFilePath(optN, optK, "trial", suffix);
        FILE *file1 = fopen(fn1.c_str(), "w");
        fprintf(file1, "%s", metadata.str().c_str());
        fprintf(file1, "r,psi,cycles,N0,N1,N2,t\n");

        std::string fn2 = defaultFilePath(optN, optK, "phasetrial", suffix);
        FILE *file2 = fopen(fn2.c_str(), "w");
        fprintf(file2, "%s", metadata.str().c_str());
        fprintf(file2, "phi1,...,phiN,t\n");

        clock_gettime(CLOCK_MONOTONIC, &start);
        for (int i=0; i<iters; i++) {
            interval++;
            if (interval == log_interval) {
                interval = 0;
                int window_size = (num_waves > 0) ? trial.N/(3*num_waves)/2 : 0;
                fprintf(
                        file1,
                        "%f,%f,%d,%d,%d,%d,%f\n",
                        kuramotoOP(trial),
                        psiOP(trial),
                        cycles(trial, window_size),
                        trial.N0, trial.N1, trial.N2,
                        trial.t);
                compressAndLogPhases(file2, trial, log_both);
            }

            ++hist[update(trial)];  //Update dynamics and log which site transitioned
        }
        clock_gettime(CLOCK_MONOTONIC, &finish);

        fclose(file1);
        fclose(file2);
        std::cout << "Written to\n" << fn1 << "\n" << fn2 << "\n";
    } else if (log_phases > 0) {
        FILE *file = fopen(filename.c_str(), "w");
        fprintf(file, "%s", metadata.str().c_str());
        fprintf(file, "phases\n");
        clock_gettime(CLOCK_MONOTONIC, &start);
        for (int i=0; i<iters; i++) {
            interval++;
            if (interval == log_interval) {
                interval = 0;
                compressAndLogPhases(file, trial, log_phases);
            }

            ++hist[update(trial)];  //Update dynamics and log which site transitioned
        }
        clock_gettime(CLOCK_MONOTONIC, &finish);

        fclose(file);
        std::cout << "Written to\n" << filename << "\n";
    } else {
        FILE *file = fopen(filename.c_str(), "w");
        fprintf(file, "%s", metadata.str().c_str());
        fprintf(file, "r,psi,cycles,N0,N1,N2,t\n");
        clock_gettime(CLOCK_MONOTONIC, &start);
        for (int i=0; i<iters; i++) {
            interval++;
            if (interval == log_interval) {
                interval = 0;
                int window_size = (num_waves > 0) ? trial.N/(3*num_waves)/2 : 0;
                fprintf(file, "%f,%f,%d,%d,%d,%d,%f\n",
                        kuramotoOP(trial),
                        psiOP(trial),
                        cycles(trial, window_size),
                        trial.N0, trial.N1, trial.N2,
                        trial.t);
            }

            ++hist[update(trial)];  //Update dynamics and log which site transitioned
        }
        clock_gettime(CLOCK_MONOTONIC, &finish);

        fclose(file);
        std::cout << "Written to\n" << filename << "\n";
    }

    if (trial.N<50) {
        for (auto p : hist) {
            std::cout << std::fixed << std::setprecision(1) << std::setw(2)
                      << p.first << ' '
                      << std::string(p.second/(1+4*iters/1000), '*') << '\n';
        }
    }
    FILE *diag = fopen("diagnostic.dat", "w");
    fprintf(diag, "freq\n");
    for (auto p : hist) fprintf(diag, "%d\n", p.second);
    fclose(diag);

    double runtime = finish.tv_sec - start.tv_sec;
    runtime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout << "Run time: " << runtime << " s\n\n";

    return 0;
}
