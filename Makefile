TARGETS = trial chi

all: $(TARGETS)

trial: trial.cpp dynamics.hpp
	g++ trial.cpp -o sim_trial -Wall -Ipcg_random -lm -O3 -march=native

chi: chicurve.cpp dynamics.hpp
	g++ chicurve.cpp -o sim_chi -Wall -Ipcg_random -lm -fopenmp -O3 -march=native

clean:
	rm -f $(TARGETS)
