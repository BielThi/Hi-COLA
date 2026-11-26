# Build scripts/pofk_snapshot/pofk_snapshot (P(k) of a Gadget snapshot with Hi-COLA's estimator).
# Same compiler, FML library and FFTW as Hi-COLA's nbody (HiCOLA/FML/COLASolver/Makefile), from the
# hicola conda env:
#   conda activate hicola && make -C scripts/pofk_snapshot
#
# Two programs, as the particle ID width of the GADGET files is fixed when compiling FML's reader
# (GADGET_LONG_INT_IDS) and the snapshots have either:
#   pofk_snapshot          4-byte IDs (nbody built without USE_LONG_IDS, snapshots before 2026-11)
#   pofk_snapshot_longids  8-byte IDs (nbody built with USE_LONG_IDS = true)
# jobs/pofk_snapshot.sh reads the width from the snapshot and runs the matching one.

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

all: pofk_snapshot pofk_snapshot_longids
.PHONY: all clean

pofk_snapshot: $(addprefix build/ids4/,$(OBJS))
	$(CC) -o $@ $^ $(OPTIONS) $(LIB) $(LINK)

pofk_snapshot_longids: $(addprefix build/ids8/,$(OBJS))
	$(CC) -o $@ $^ $(OPTIONS) -DGADGET_LONG_INT_IDS $(LIB) $(LINK)

build/ids4/%.o: %.cpp Makefile | build/ids4
	$(CC) -c -o $@ $< $(OPTIONS) $(INC)

build/ids8/%.o: %.cpp Makefile | build/ids8
	$(CC) -c -o $@ $< $(OPTIONS) -DGADGET_LONG_INT_IDS $(INC)

build/ids4 build/ids8:
	mkdir -p $@

clean:
	rm -rf pofk_snapshot pofk_snapshot_longids build *.o
