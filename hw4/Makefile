CC=mpicc
FLAGS=-O3 -lm
EXECS= ssort jacobi-mpi2D jacobi-mpi2D-nonBlocking mpi_solved1 mpi_solved2 mpi_solved3 mpi_solved4 mpi_solved5 mpi_solved6 mpi_solved7 

all: ${EXECS}

ssort: ssort.c
	${CC} ${FLAGS} $^ -o $@

jacobi-mpi2D: jacobi-mpi2D.c
	${CC} ${FLAGS} $^ -o $@

jacobi-mpi2D-nonBlocking: jacobi-mpi2D-nonBlocking.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved1: mpi_bug1.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved2: mpi_bug2.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved3: mpi_bug3.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved4: mpi_bug4.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved5: mpi_bug5.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved6: mpi_bug6.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved7: mpi_bug7.c
	${CC} ${FLAGS} $^ -o $@

clean:
	rm -rf $(wildcard mpi_solved*)
	rm -rf jacobi-mpi2D jacobi-mpi2D.dSYM
	rm -rf ssort ssort.dSYM
	rm -rf jacobi-mpi2D-nonBlocking jacobi-mpi2D-nonBlocking.dSYM
	rm -rf output*.txt