#include <string>
#include <getopt.h>
#include <iomanip>
#include <omp.h>
#include <vector>
#include "dynamics.hpp"
#include "logging.hpp"

struct trialData {
    double r, r2, psi, psi2, duration;
};

struct trialData accumulateTrial(struct trial &trial, int iters, int burn, int log_interval)
{
    for (int i=0; i<burn; i++) {
        update(trial);
    }

    double r = 0;
    double r2 = 0;
    double psi = 0;
    double psi2 = 0;
    double duration = 0;
    double dt = 0;
    int interval = 0;
    if (log_interval <= 0) log_interval = 1;
    for (int i=0; i<iters; i++) {
        update(trial);
        dt += trial.dt;
        interval++;
        if (interval == log_interval) {
            double tmp = kuramotoOP(trial) * dt;
            r += tmp;
            r2 += tmp*tmp;
            tmp = psiOP(trial) * dt;
            psi += tmp;
            psi2 += tmp*tmp;

            interval = 0;
            dt = 0;
        }
    }
    duration        = trial.t;
    struct trialData tdat = {
        r / duration,
        r2 / duration,
        psi / duration,
        psi2 / duration,
        duration
    };
    return tdat;
}

int main(int argc, char *argv[])
{
    int optN = 500;
    int optK = 50;
    double a = 1.0;
    double A = 4.0;
    int na = 20;
    int maxthrds = 6;
    int num_trials = 100;
    int iters = 50000;
    int burn = 30000;
    int seed = 23;
    int log_interval = 1;
    std::string filename = "";
    std::string IC = "uniform";

    char opt;
    int opt_idx = 0;
    char shortopt[] = {"n:k:t:i:b:a:A:o:c:"};
    static struct option longopt[] = {
        {"size",          required_argument, 0, 'n'},
        {"range",         required_argument, 0, 'k'},
        {"trials",        required_argument, 0, 't'},
        {"iters",         required_argument, 0, 'i'},
        {"burn",          required_argument, 0, 'b'},
        {"amin",          required_argument, 0, 'a'},
        {"amax",          required_argument, 0, 'A'},
        {"na",            required_argument, 0,  3 },
        {"maxthrds",      required_argument, 0,  4 },
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
    metadata << "N=" << optN << " K=" << optK << " a=" << a << " A=" << A
             << " na=" << na
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

    FILE *file = fopen(filename.c_str(), "w");
    fprintf(file, "%s", metadata.str().c_str());
    fprintf(file, "a,r,psi\n");
    for (int i=0; i<na; i++) {
        double coupl = (na > 1) ? a + (A-a)*i/(na-1) : a;
        printf("[%02d/%02d] Running a=%.4f",i+1, na, coupl);
        double trial_avg_r   = 0;
        double trial_avg_psi = 0;
#pragma omp parallel num_threads(maxthrds) default(none) \
        shared(optN,optK,coupl,iters,burn,IC,seed,num_trials,log_interval) \
        reduction(+:trial_avg_r,trial_avg_psi)
        {
#pragma omp single
            {
                printf(" [%d threads]", omp_get_num_threads());
            }
#pragma omp for 
            for (int i=0; i<num_trials; i++) {
                struct trial trial;
                struct trialData tdata;

                pcg32 RNG(seed, i);
                std::function<void(struct trial &)> initStates = NULL;
                std::function<void(struct trial &)> initNaturalFreqs = NULL;
                if (IC == "uniform") {                                                       
                    initStates = [](struct trial &t){                                        
                        initUniform(t, 0);                                                   
                    };                                                                       
                } else if (IC.substr(0, 4) == "wave") {                                      
                    int num_waves = stoi(IC.substr(4));                                      
                    initStates = [num_waves](struct trial &t){                               
                        initWave(t, num_waves, false);                                       
                    };                                                                       
                }
                initNaturalFreqs = [](struct trial &t) {
                    initGaussianNaturalFreqs(t, 1.0, 0.1);
                };

                initTrial(trial, optN, optK, coupl, RNG, initStates, initNaturalFreqs);
                tdata = accumulateTrial(trial, iters, burn, log_interval);
                trial_avg_r   += tdata.r;
                trial_avg_psi += tdata.psi;
            }
        }
        trial_avg_r   /= num_trials;
        trial_avg_psi /= num_trials;
        fprintf(file, "%f,%f,%f\n", coupl, trial_avg_r, trial_avg_psi);
        printf(" Done!  <r> = %.4f, <psi> = %.4f\n", trial_avg_r, trial_avg_psi);
    }
    printf("Done, %s\n", filename.c_str());
    fclose(file);

    return 0;
}
