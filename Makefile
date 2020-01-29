TARGETS = trial chi test
TGTS = $(patsubst %, %$(TARGETS), {sim_, sim_, sim_})

all: $(TARGETS)

trial: trial.cpp dynamics.hpp logging.hpp
	g++ trial.cpp -o sim_trial -Wall -Ipcg_random -lm -O3 -march=native

chi: chicurve.cpp dynamics.hpp logging.hpp
	g++ chicurve.cpp -o sim_chi -Wall -Ipcg_random -lm -fopenmp -O3 -march=native

test: test.cpp dynamics.hpp logging.hpp
	g++ test.cpp -o sim_test -Wall -Ipcg_random -lm -O3 -march=native

clean:
	rm -f $(TGTS)
