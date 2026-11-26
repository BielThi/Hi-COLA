# Build scripts/pofk_snapshot/pofk_snapshot (P(k) of a Gadget snapshot with Hi-COLA's estimator).
# Same compiler, FML library and FFTW as Hi-COLA's nbody (HiCOLA/FML/COLASolver/Makefile), from the
# hicola conda env:
#   conda activate hicola && make -C scripts/pofk_snapshot

SHELL := /bin/bash

# MPI only (nbody also uses OpenMP, but the jobs run with OMP_NUM_THREADS=1); double precision grids
# like nbody (no SINGLE_PRECISION_FFTW)
CC          = mpicxx -std=c++1z -O3 -Wall -Wextra -march=native
FML_INCLUDE = $(abspath $(CURDIR)/../../HiCOLA)
OPTIONS     = -DUSE_MPI -DUSE_FFTW
INC         = -I$(FML_INCLUDE) -I$(CONDA_PREFIX)/include
LIB         = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
LINK        = -lfftw3_mpi -lfftw3

VPATH := $(FML_INCLUDE)/FML/Global/:$(FML_INCLUDE)/FML/GadgetUtils/:$(FML_INCLUDE)/FML/FileUtils/
OBJS  = Main.o Global.o GadgetUtils.o FileUtils.o

all: pofk_snapshot
.PHONY: all clean

pofk_snapshot: $(OBJS)
	$(CC) -o $@ $^ $(OPTIONS) $(LIB) $(LINK)

%.o: %.cpp Makefile
	$(CC) -c -o $@ $< $(OPTIONS) $(INC)

clean:
	rm -f pofk_snapshot *.o
