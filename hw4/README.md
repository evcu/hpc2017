## HW4 - HPC 2017 - Utku Evci
There are bugs and MPI-jacobi and sample-sort algorithm. To test the algorithms run. `mpi_solved5` runs forever by nature.
```
git clone https://github.com/evcu/hpc2017.git
cd hpc2017/hw4
module load intel-2016
make
mpirun -np 4 ./mpi_solved1
mpirun -np 4 ./mpi_solved2
mpirun -np 4 ./mpi_solved3
mpirun -np 4 ./mpi_solved4
mpirun -np 4 ./mpi_solved5
mpirun -np 4 ./mpi_solved6
mpirun -np 4 ./mpi_solved7
mpirun -np 16 ./jacobi-mpi2D 100 100
mpirun -np 16 ./jacobi-mpi2D-nonBlocking 100 100
mpirun -np 10 ./ssort 100
python test.py 10
```

Here the last line is a script that tests the ouput files whether they are sorted or not with a scan through output files and lines

I've used following scripts on stampede
```
sbatch q2
sbatch q3
```

The output files of these jobs are coppied into `console_output/` folder.

### Results
Results/Report is [here](https://github.com/evcu/hpc2017/blob/master/hw4/result/Results.ipynb)
