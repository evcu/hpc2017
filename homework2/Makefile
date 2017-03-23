CC=gcc-6
FLAGS=-O3 -fopenmp -Wall -g -lm
FLAGS2=-O3 -Wall -g -lm
EXECS=omp_solved2 omp_solved3 omp_solved4 omp_solved5 omp_solved6 jacobi2D-omp gs2D-omp
all: ${EXECS}

no-omp: jacobi2D-no-omp gs2D-no-omp

omp_solved2: omp_bug2.c
	${CC} ${FLAGS} $^ -o $@

omp_solved3: omp_bug3.c
	${CC} ${FLAGS} $^ -o $@

omp_solved4: omp_bug4.c
	${CC} ${FLAGS} $^ -o $@

omp_solved5: omp_bug5.c
	${CC} ${FLAGS} $^ -o $@

omp_solved6: omp_bug6.c
	${CC} ${FLAGS} $^ -o $@

jacobi2D-omp: jacobi2D-omp.c
	${CC} ${FLAGS} $^ -o $@

gs2D-omp: gs2D-omp.c
	${CC} ${FLAGS} $^ -o $@

jacobi2D-no-omp: jacobi2D-omp.c
	${CC} ${FLAGS2} $^ -o $@

gs2D-no-omp: gs2D-omp.c
	${CC} ${FLAGS2} $^ -o $@

clean:
	rm -rf $(wildcard omp_solved*)
	rm -rf jacobi2D-omp jacobi2D-omp.dSYM
	rm -rf gs2D-omp gs2D-omp.dSYM
	rm -rf gs2D-no-omp gs2D-no-omp.dSYM
	rm -rf jacobi2D-no-omp jacobi2D-no-omp.dSYM

