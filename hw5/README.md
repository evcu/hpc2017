## HW5 - HPC 2017 - Utku Evci
HW-5 Mpi implementation of 1d_multi_grid. There was a fundamental difference between the previous Jacobi implementations and multi-grid implementations. The intersecting border between neigbors was two numbers for previous implementations and it is chosen to be 1 for the serial sample implementation. Therefore a special handling needed to be done. 
```
git clone https://github.com/evcu/hpc2017-hw5.git
cd hpc2017-hw5
module load intel-2016
make
mpirun -np 4 ./multigrid_1d_mpi $((2**20)) 100 //5.250111 s
./multigrid_1d $((2**20)) 100 //7.033151 s
make clean
```

Here the last line is a script that tests the ouput files whether they are sorted or not with a scan through output files and lines

I've used following scripts on stampede
```
sbatch strongScaling.job
sbatch weakScaling.job
```

The output files of these jobs are coppied into `results/` folder.

### Results
Results/Report is [here](https://github.com/evcu/hpc2017-hw5/blob/master/results/Results.ipynb)
