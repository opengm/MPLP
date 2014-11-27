CC=g++
CFLAGS=-g -O3 #-Wall #-Wextra
LDFLAGS=

MPLP_CYCLE_ALG_TRIPLET=muldim_arr.o read_model_file.o mplp_alg.o cycle_tighten_main.o

EXECUTABLES=solver

.cpp.o:
	g++ $(CFLAGS) -c $<

all: solver

solver: $(MPLP_CYCLE_ALG_TRIPLET)
	g++ $(CFLAGS) $(MPLP_CYCLE_ALG_TRIPLET) -o $@

clean:
	rm -rf *.o $(EXECUTABLES)

cycle_tighten_main.o: cycle.h

muldim_arr.o: muldim_arr.h

read_model_file.o: read_model_file.h

mplp_alg.o: mplp_alg.h