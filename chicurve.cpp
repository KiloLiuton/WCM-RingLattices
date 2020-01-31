#include <string>
#include <getopt.h>
#include <iomanip>
#include <omp.h>
#include <vector>
#include "dynamics.hpp"
#include "logging.hpp"

int main(int argc, char *argv[])
{
    int optN = 1000;
    int optK = 200;
    double a = 1.0;
    double A = 4.0;
    int na = 20;
    int maxthrds = 8;
    int num_trials = 50;
    int iters = 50000;
    int burn = 50000;
    int seed = 23;
    int log_interval = 5;
    std::string filename = "";
    std::string IC = "uniform";

    char opt;
    int opt_idx = 0;
    char shortopt[] = {"n:k:t:i:b:a:A:o:c:"};
    static struct option longopt[] = {
        {"size",          required_argument, 0, 'n'},
        {"range",         required_argument, 0, 'k'},
        {"amin",          required_argument, 0, 'a'},
        {"amax",          required_argument, 0, 'A'},
        {"na",            required_argument, 0, '3'},
        {"maxthrds",      required_argument, 0, '4'},
        {"trials",        required_argument, 0, 't'},
        {"iters",         required_argument, 0, 'i'},
        {"burn",          required_argument, 0, 'b'},
        {"log-interval",  required_argument, 0,  1 },
        {"random-seed",   required_argument, 0,  2 },
        {"outfile",       required_argument, 0, 'o'},
        {"IC",            required_argument, 0, 'c'},
        {0,               0,                 0,  0 }
    };
    while ((opt=getopt_long(argc, argv, shortopt, longopt, &opt_idx)) != -1) {
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
        case 't':
            num_trials = atoi(optarg);
            break;
        case 'i':
            iters = atoi(optarg);
            break;
        case 'b':
            burn = atoi(optarg);
            break;
        case 'a':
            a = atof(optarg);
            break;
        case 'A':
            A = atof(optarg);
            break;
        case 1:
            log_interval = abs(atoi(optarg));
            break;
        case 2:
            seed = atoi(optarg);
            break;
        case 3:
            na = atoi(optarg);
            break;
        case 4:
            maxthrds = atoi(optarg);
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
        const char* prefix = "chicurve";
        filename = defaultFilePath(optN, optK, prefix, "");
    }

    std::ostringstream metadata;
    metadata << "N=" << optN << " K=" << optK
             << " ic=" << IC.c_str() << " log-interval=" << log_interval
             << " iters=" << iters << " burn=" << burn
             << " seed=" << seed << "\n";
    printf("Initializing with:\n"
           "%s"
           "log-interval %d\n"
           "ic %s\n",
           metadata.str().c_str(),
           log_interval,
           IC.c_str()
           );
    omp_set_num_threads(maxthrds);

    std::vector<double> couplings;
    for (int i=0; i<na; i++) couplings.push_back(a+i*(A-a)/(na-1));
    double trial_avg_r   = 0;
    double trial_avg_psi = 0;
    for (auto a: couplings) {
        #pragma omp parallel default(none) shared(std::cout,optN,optK,a,iters,burn,IC,seed,num_trials) reduction(+:trial_avg_r,trial_avg_psi)
        {
            #pragma omp single nowait
            {
                int nthrds = omp_get_num_threads();
                printf("Memory usage: %.2f MB\n", nthrds*sizeof(struct trial)/1e3);
            }
            #pragma omp for
            for (int i=0; i<num_trials; i++) {
                std::cout << "Hi, thrd " << omp_get_thread_num() << std::endl;
                struct trial trial;
                std::function<void(struct trial &)> initStates = NULL;
                std::function<void(struct trial &)> initNaturalFreqs = NULL;
                if (IC == "uniform") {
                    initStates = [](struct trial &t){ initUniform(t, 0); };
                } else if (IC.substr(0, 4) == "wave") {
                    int num_waves = stoi(IC.substr(4));
                    initStates = [num_waves](struct trial &t){ initWave(t, num_waves, false); };
                }
                pcg32 RNG(seed, i);
                initTrial(trial, optN, optK, a, RNG, initStates, initNaturalFreqs);
                double r = 0, psi = 0;
                for (int j=0; j<burn; j++) update(trial);
                for (int j=0; j<iters; j++) {
                    r   += kuramotoOP(trial) / trial.rates_sum;
                    psi += psiOP(trial)      / trial.rates_sum;
                    update(trial);
                }
                trial_avg_r += r;
                trial_avg_psi += psi;
            }
        }
    }

    return 0;
}
