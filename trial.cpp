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
#include "dynamics.hpp"

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    /*Format a string with C-like syntax*/
    size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    std::unique_ptr<char[]> buf( new char[ size ] );
    snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

std::string defaultFileName(int N, int K, double a, int log_phases)
{
    int v = 0;
    const char* prefix = (log_phases > 0) ? "phasetrial" : "trial";
    std::string fname = string_format(
            "%s-N-%06dK-%05da-%.6f_v%d.dat",
            prefix, N, K, a, v
            );
    while (std::ifstream(fname)) {
        v++;
        fname = string_format(
                "%s-N-%06dK-%05da-%.6f_v%d.dat",
                prefix, N, K, a, v
                );
    }
    return fname;
}

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
    int optN = 1000;
    int optK = 300;
    int iters = 50000;
    int seed = 23;
    int stream = 42;
    double opta = 1.6;
    int log_phases = 0;
    int log_both = 0;
    int log_interval = 0;
    std::string filename = "";
    std::string IC = "uniform";

    char opt;
    int opt_idx = 0;
    std::string sopt = "n:k:a:i:s:r:o:c:";
    static struct option lopt[] = {
        {"size",          required_argument, 0, 'n'},
        {"range",         required_argument, 0, 'k'},
        {"coupling",      required_argument, 0, 'a'},
        {"iters",         required_argument, 0, 'i'},
        {"log-phases",    required_argument, 0,  1 },
        {"log-both",      required_argument, 0,  2 },
        {"log-interval",  required_argument, 0,  3 },
        {"random-seed",   required_argument, 0, 's'},
        {"random-stream", required_argument, 0, 'r'},
        {"outfile",       required_argument, 0, 'o'},
        {"IC",            required_argument, 0, 'c'},
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
            IC = std::string(optarg);
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
        filename = defaultFileName(optN, optK, opta, log_phases);
    }
    std::ostringstream metadata;
    metadata << "N=" << optN << " K=" << optK << " a=" << opta
             << " ic=" << IC.c_str() << " log-interval=" << log_interval
             << " iters=" << iters << " seed=" << seed
             << " stream=" << stream << "\n";

    // Initialize system with command line options
    struct trial trial;
    if (IC == "uniform")
        initUniform(trial, optN, optK, opta, seed, stream);
    else if (IC == "wave")
        initWave(trial, optN, optK, opta, seed, stream, 6);

    printf("Memory usage: %.2f MB\n", sizeof(trial)/1e3);
    printf("Initializing with:\n"
           "%s"
           "log-phases %d\n"
           "log-interval %d\n"
           "ic %s\n",
           metadata.str().c_str(),
           log_phases,
           log_interval,
           IC.c_str()
           );

    
    std::map<int, int> hist; // count states that transitioned

    int interval=1;
    struct timespec start, finish;
    if (log_both > 0) {
        std::string fn1 = defaultFileName(optN, optK, opta, 0);
        FILE *file1 = fopen(fn1.c_str(), "w");
        fprintf(file1, "%s", metadata.str().c_str());
        fprintf(file1, "r,psi,N0,N1,N2,t\n");

        std::string fn2 = defaultFileName(optN, optK, opta, 1);
        FILE *file2 = fopen(fn2.c_str(), "w");
        fprintf(file2, "%s", metadata.str().c_str());
        fprintf(file2, "phases\n");

        clock_gettime(CLOCK_MONOTONIC, &start);
        for (int i=0; i<iters; i++) {
            if (interval > log_interval) {
                interval = 1;
                fprintf(
                        file1,
                        "%f,%f,%d,%d,%d,%f\n",
                        kuramotoOP(trial),
                        psiOP(trial),
                        trial.N0, trial.N1, trial.N2,
                        trial.t);
                compressAndLogPhases(file2, trial, log_both);
            }
            interval++;

            int index = update(trial); // update the dynamics
            ++hist[index];             // log which site transitioned
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
            if (interval > log_interval) {
                interval = 1;
                compressAndLogPhases(file, trial, log_phases);
            }
            interval++;

            int index = update(trial); // update the dynamics
            ++hist[index];             // log which site transitioned
        }
        clock_gettime(CLOCK_MONOTONIC, &finish);

        fclose(file);
    } else {
        FILE *file = fopen(filename.c_str(), "w");
        fprintf(file, "%s", metadata.str().c_str());
        fprintf(file, "r,psi,N0,N1,N2,t\n");
        clock_gettime(CLOCK_MONOTONIC, &start);
        for (int i=0; i<iters; i++) {
            if (interval > log_interval) {
                interval = 1;
                fprintf(file, "%f,%f,%d,%d,%d,%f\n",
                        kuramotoOP(trial),
                        psiOP(trial),
                        trial.N0, trial.N1, trial.N2,
                        trial.t);
            }
            interval++;

            int index = update(trial); // update the dynamics
            ++hist[index];             // log which site transitioned
        }
        clock_gettime(CLOCK_MONOTONIC, &finish);

        fclose(file);
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
    std::cout << "Run time: " << runtime << " s\n";

    return 0;
}