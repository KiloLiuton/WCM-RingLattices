TARGETS = trial chi test
TGTS = $(patsubst %, sim_%, $(TARGETS))

MAXN ?= 25001

all: $(TARGETS)

trial: trial.cpp dynamics.hpp logging.hpp
	g++ trial.cpp -o sim_trial -Wall -Ipcg_random -lm -O3 -march=native -DMAXN=$(MAXN)

chi: chicurve.cpp dynamics.hpp logging.hpp
	g++ chicurve.cpp -o sim_chi -Wall -Ipcg_random -lm -fopenmp -O3 -march=native -DMAXN=$(MAXN)

test: test.cpp dynamics.hpp logging.hpp
	g++ test.cpp -o sim_test -Wall -Ipcg_random -lm -O3 -march=native -DMAXN=$(MAXN) -fopenmp

clean:
	rm -f $(TGTS)
