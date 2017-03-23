#!/bin/bash
echo `date `

N_array=( 100 1000 2000 )
ITER=1000
for var in "$@"
do
	for i in "${N_array[@]}"
	do
	    export OMP_NUM_THREADS=$var
		echo "$var $i"
	    ../jacobi2D-omp $i 1000 | grep Time
	done
done