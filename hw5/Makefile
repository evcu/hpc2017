EXECS= multigrid_1d multigrid_1d_mpi
CC=gcc
MPICC=mpicc

FLAGS=-O3 -fopenmp -Wall -lm -g

all: ${EXECS}

multigrid_1d: multigrid_1d.c
	${CC} ${FLAGS} $^ -o multigrid_1d

multigrid_1d_mpi: multigrid_1d_mpi.c
	${MPICC} ${FLAGS} -o multigrid_1d_mpi $^


clean:
	rm -rf ${EXECS} $(wildcard *.dSYM)
