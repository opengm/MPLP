CC=g++
CFLAGS=-g -O3 #-Wall #-Wextra
LDFLAGS=
INCLUDES := -I./include

MPLP_CYCLE_ALG_TRIPLET=src/muldim_arr.o src/read_model_file.o src/mplp_alg.o src/cycle_tighten_main.o
MPLP_CYCLE_ALG_TRIPLET2=muldim_arr.o read_model_file.o mplp_alg.o cycle_tighten_main.o

EXECUTABLES=solver

.cpp.o:
	g++ $(CFLAGS) ${INCLUDES} -c $<

all: solver

solver: $(MPLP_CYCLE_ALG_TRIPLET)
	g++ $(CFLAGS) $(MPLP_CYCLE_ALG_TRIPLET2) -o $@

clean:
	rm -rf *.o $(EXECUTABLES)

src/cycle_tighten_main.o: ./include/MPLP/cycle.h

src/muldim_arr.o: ./include/MPLP/muldim_arr.h

src/read_model_file.o: ./include/MPLP/read_model_file.h

src/mplp_alg.o: ./include/MPLP/mplp_alg.h