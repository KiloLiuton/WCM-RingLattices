#include <string>
#include <getopt.h>
#include <iomanip>
#include <omp.h>
#include <vector>
#include "dynamics.hpp"
#include "logging.hpp"

struct trialData {
    double r, r2, psi, psi2, cycl, duration;
};

struct trialData accumulateTrial(struct trial &trial, int iters, int burn, int log_interval, int num_waves)
{
    for (int i=0; i<burn; i++) {
        update(trial);
    }

    double r = 0;
    double r2 = 0;
    double psi = 0;
    double psi2 = 0;
    double cycl = 0;
    double duration = 0;
    double dt = 0;
    int interval = 0;
    if (log_interval < 1) log_interval = 1;
    for (int i=0; i<iters; i++) {
        update(trial);
        dt += trial.dt;
        interval++;
        if (interval == log_interval) {
            double tmp = kuramotoOP(trial) * dt;
            r += tmp;
            r2 += tmp*tmp;
            tmp = psiOP(trial);
            psi += tmp * dt;
            psi2 += tmp * tmp * dt;
            int window_size = (num_waves > 0) ? trial.N/(3*num_waves)/2 : 0;
            cycl += (double) cycles(trial, window_size) * dt;
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
        cycl / duration,
        duration
    };
    return tdat;
}

int main(int argc, char *argv[])
{
    int optN = 200;
    int optK = 99;
    double a = 1.0;
    double A = 3.6;
    int na = 30;
    double gmean = 1.0;
    double gstddev = 0.1;
    int maxthrds = 10;
    int num_trials = 100;
    int iters = 9000;
    int burn = 2000;
    int seed = 23;
    int log_interval = 1;
    std::string filename = "";
    std::string ic = "wave2";
    int num_waves = 2;

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
        {"gmean",         required_argument, 0,  5 },
        {"gstddev",       required_argument, 0,  6 },
        {"maxthrds",      required_argument, 0,  4 },
        {"log-interval",  required_argument, 0,  1 },
        {"random-seed",   required_argument, 0,  2 },
        {"outfile",       required_argument, 0, 'o'},
        {"ic",            required_argument, 0, 'c'},
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
        case 5:
            gmean = atof(optarg);
            break;
        case 6:
            gstddev = atof(optarg);
            break;
        case 4:
            maxthrds = atoi(optarg);
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
        const char* prefix = "chicurve";
        filename = defaultFilePath(optN, optK, prefix, "");
    }

    std::ostringstream metadata;
    metadata << "N=" << optN << " K=" << optK << " a=" << a << " A=" << A
             << " na=" << na << " gmean=" << gmean << " gstddev=" << gstddev
             << " num_trials=" << num_trials
             << " ic=" << ic.c_str() << " log-interval=" << log_interval
             << " iters=" << iters << " burn=" << burn
             << " seed=" << seed << "\n";
    printf("Initializing with:\n"
           "%s"
           "log-interval %d\n"
           "ic %s\n",
           metadata.str().c_str(),
           log_interval,
           ic.c_str()
           );

    FILE *file = fopen(filename.c_str(), "w");
    fprintf(file, "%s", metadata.str().c_str());
    fprintf(file, "a,r,psi,cycles,chir,chipsi,duration\n");

    // Define function to initialize natural frequencies distribution w ~ g
    std::function<void(struct trial &)> g = NULL;
    g = [gmean, gstddev](struct trial &t) {
        initGaussianNaturalFreqs(t, gmean, gstddev);
    };

    // Define function to initialize initial configuration
    std::function<void(struct trial &)> initStates = NULL;
    if (ic == "uniform") {                                                       
        initStates = [](struct trial &t){                                        
            initUniform(t, 0);                                                   
        };                                                                       
        num_waves = 0;
    } else if (ic.substr(0, 4) == "wave") {                                      
        num_waves = stoi(ic.substr(4));                                      
        initStates = [num_waves](struct trial &t){                               
            initWave(t, num_waves, false);                                       
        };                                                                       
    }
    for (int i=0; i<na; i++) {
        double coupl = (na > 1) ? a + (A-a)*i/(na-1) : a;
        double trial_avg_r      = 0;
        double trial_avg_r2     = 0;
        double trial_avg_psi    = 0;
        double trial_avg_psi2   = 0;
        double trial_avg_cycles = 0;
        double trial_avg_t      = 0;
#pragma omp parallel num_threads(maxthrds) default(none) \
        shared(optN,optK,coupl,iters,burn,seed,num_trials,log_interval,g,initStates,num_waves) \
        reduction(+:trial_avg_r,trial_avg_psi,trial_avg_r2,trial_avg_psi2,trial_avg_cycles,trial_avg_t)
        {
#pragma omp single
            {
                printf("[%d threads]", omp_get_num_threads());
            }
#pragma omp for 
            for (int i=0; i<num_trials; i++) {
                struct trial t;
                struct trialData tdat;

                pcg32 RNG(seed, i);

                initTrial(t, optN, optK, coupl, RNG, g, initStates);
                tdat = accumulateTrial(t, iters, burn, log_interval, num_waves);
                trial_avg_r      += tdat.r;
                trial_avg_r2     += tdat.r*tdat.r;
                trial_avg_psi    += tdat.psi;
                trial_avg_psi2   += tdat.psi*tdat.psi;
                trial_avg_cycles += tdat.cycl;
                trial_avg_t      += tdat.duration;
            }
        }
        trial_avg_r      /= num_trials;
        trial_avg_r2     /= num_trials;
        trial_avg_psi    /= num_trials;
        trial_avg_psi2   /= num_trials;
        trial_avg_cycles /= num_trials;
        trial_avg_t      /= num_trials;
        double chir = optN * (trial_avg_r2 - trial_avg_r*trial_avg_r);
        double chipsi = optN * (trial_avg_psi2 - trial_avg_psi*trial_avg_psi);
        fprintf(
                file,
                "%f,%f,%f,%f,%f,%f,%f\n",
                coupl, trial_avg_r, trial_avg_psi, trial_avg_cycles, chir, chipsi, trial_avg_t
                );
        printf(" a=%.4f done!  <r> =%.4f <psi> =%.4f [%02d/%02d]\n",coupl, trial_avg_r, trial_avg_psi, i+1, na);
    }
    printf("Done, %s\n", filename.c_str());
    fclose(file);

    return 0;
}
