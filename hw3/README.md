Note: I ran all the commands from crunchy4 hosting them on crunchy3 and crunchy1. As a summary:
    - 1) 0.0725s communication latency
    - 2) 49mb/s bandwith

I've run the following command to time the code and estimate the latency
> time mpirun -np 4 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 10000 > /dev/null

I run the code three times and got the following result:
```
$ time mpirun -np 4 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 10000 > /dev/null
real	0m3.358s
user	0m0.083s
sys	    0m0.565s
time mpirun -np 4 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 10000 > /dev/null
real	0m3.333s
user	0m0.061s
sys	    0m0.471s
time mpirun -np 4 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 10000 > /dev/null
real	0m3.179s
user	0m0.058s
sys	    0m0.348s
```

I assume is the communication latency is only related to the real time. The average is 3.290s and there is 4 communication happening per round. There is 10000 round. Therefore 3.290/40000=0.0725 s is the time required for each communication.

When I run

> time mpirun -np 40 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 1000 > /dev/null

such that the overall number of communication is same, I got a decrease in the speed and in this scenario it took 6-7 seconds even though the total number of communications is same.

## PART2
Enable the `#DEFINE PART2ONLY` flag at the beginning of the code.

I ran the the following command three times

```
$ time mpirun -np 6 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 20 > /dev/null

real	0m4.766s
user	0m0.028s
sys	    0m0.080s

time mpirun -np 6 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 20 > /dev/null

real	0m4.816s
user	0m0.032s
sys	    0m0.073s

time mpirun -np 6 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 20 > /dev/null

real	0m4.930s
user	0m0.026s
sys	    0m0.079s 
```

Total number of bytes transferred is `(2000000*120.0)`. The the bandwith is therefore `(2000000*120.0)/((766+816+930)/3000.0+4)=49,614,112` bytes/second, which is almost 50mb/s.

As observed in the first part I got a decrease in speed of factor 2 when I increased the number of processes per machine.

```
$ time mpirun -np 60 --hosts crunchy1,crunchy3 -perhost 1 ./int_ring 2 > /dev/null

real	0m8.172s
user	0m0.036s
sys	    0m0.120s
```
