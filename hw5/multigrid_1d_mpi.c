/* Multigrid for solving -u''=f for x in (0,1)
 * Usage: ./multigrid_1d < Nfine > < iter > [s-steps]
 * NFINE: number of intervals on finest level, must be power of 2
 * ITER: max number of V-cycle iterations
 * S-STEPS: number of Jacobi smoothing steps; optional
 * Author: Georg Stadler, April 2017
 */
#include <stdio.h>
#include <math.h>
#include "util.h"
#include <string.h>
#include <mpi.h>
int mpirank, p;

/* compuate norm of residual */
double compute_norm(double *u, int N)
{
  int i;
  double norm = 0.0, allnorm = 0.0;
  for (i = 0; i <= N; i++)
    norm += u[i] * u[i];

  /* use allreduce for convenience; a reduce would also be sufficient */
  MPI_Allreduce(&norm, &allnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sqrt(allnorm);
}

/* set vector to zero */
void set_zero (double *u, int N) {
  int i;
  for (i = 0; i <= N; i++)
    u[i] = 0.0;
}

/* debug function */
void output_to_screen (double *u, int N,char *c) {
  int i;
  printf("%d: %s ",mpirank, c);
  for (i = 0; i < N+1; 
    i++)
    printf("%f ", u[i]);
  printf("\n");
}

/* coarsen uf from length N+1 to lenght N/2+1
   assuming N = 2^l
*/
void coarsen(double *uf, double *uc, int N) {
  int ic;
  double ll,rr;
  MPI_Status status, status1;
  for (ic = 1; ic < N/2; ++ic)
    uc[ic] = 0.5 * uf[2*ic] + 0.25 * (uf[2*ic-1]+uf[2*ic+1]);

    if (mpirank < p - 1) {
      /* If not the last process, send/recv bdry values to the right */
      MPI_Send(&uf[N-1], 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD);
      MPI_Recv(&rr, 1, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &status);
    }
    if (mpirank > 0) {
      /* If not the first process, send/recv bdry values to the left */
      MPI_Send(&uf[1], 1, MPI_DOUBLE, mpirank-1, 123, MPI_COMM_WORLD);
      MPI_Recv(&ll, 1, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &status1);
    }

  if (mpirank > 0) {
    //get uf[-1]
    // printf("ll: %f",ll);
    uc[0] = 0.5 * uf[0] + 0.25 * (ll+uf[1]);
  }
   if (mpirank  < p - 1) {
    //get -2 uf[N+1]
    // printf("rr: %f",rr);
    uc[N/2] = 0.5 * uf[N] + 0.25 * (uf[N-1]+rr);
  }
}


/* refine u from length N+1 to lenght 2*N+1
   assuming N = 2^l, and add to existing uf
*/
void refine_and_add(double *u, double *uf, int N)
{
  MPI_Status status;
  int i;
  uf[1] += 0.5 * (u[0] + u[1]);
  for (i = 1; i < N; ++i) {
    uf[2*i] += u[i];
    uf[2*i+1] += 0.5 * (u[i] + u[i+1]);
  }
  if (mpirank < p - 1) {
  /* If not the last process, send/recv bdry values to the right */
  uf[2*N] += u[N];
  MPI_Send(&(uf[2*N]), 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD);
  }
  if (mpirank > 0) {
    /* If not the first process, send/recv bdry values to the left */
    MPI_Recv(&(uf[0]), 1, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &status);
  }
}

/* compute residual vector */
void compute_residual(double *u, double *rhs, double *res, int N, double invhsq)
{
  int i;
   double ll,rr;
  MPI_Status status, status1;

  if (mpirank < p - 1) {
      /* If not the last process, send/recv bdry values to the right */
      MPI_Send(&u[N-1], 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD);
      MPI_Recv(&rr, 1, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &status);
    }
  if (mpirank > 0) {
      /* If not the first process, send/recv bdry values to the left */
      MPI_Send(&u[1], 1, MPI_DOUBLE, mpirank-1, 123, MPI_COMM_WORLD);
      MPI_Recv(&ll, 1, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &status1);
    }

  for (i = 1; i < N; i++)
    res[i] = (rhs[i] - (2.*u[i] - u[i-1] - u[i+1]) * invhsq);

   if (mpirank > 0) {
    //get uf[-1]
    res[0] = (rhs[0] - (2.*u[0] - ll - u[1]) * invhsq);
  }
   if (mpirank  < p - 1) {
    //get -2 uf[N+1]
    res[N] = (rhs[N] - (2.*u[N] - u[N-1] - rr) * invhsq);
    
  }

}


/* compute residual and coarsen */
void compute_and_coarsen_residual(double *u, double *rhs, double *resc,
          int N, double invhsq)
{
  double *resf = calloc(sizeof(double), N+1);
  compute_residual(u, rhs, resf, N, invhsq);
  coarsen(resf, resc, N);
  free(resf);
}

/* Perform Jacobi iterations on u */
void jacobi(double *lu, double *rhs, int lN, double hsq, int ssteps)
{
  int i, j;
  double ll,rr;
  MPI_Status status;
  MPI_Request request_out1, request_in1;
  MPI_Request request_out2, request_in2;

  /* Jacobi damping parameter -- plays an important role in MG */
  double omega = 2./3.;
  double *lunew = calloc(sizeof(double), lN+2);
  for (j = 0; j < ssteps; ++j) {

    if (mpirank < p - 1) {
      /* If not the last process, send/recv bdry values to the right */
      MPI_Irecv(&rr, 1, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &request_in1);
      MPI_Isend(&(lu[lN-1]), 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD, &request_out1);
    }
    if (mpirank > 0) {
      /* If not the first process, send/recv bdry values to the left */
      MPI_Irecv(&ll, 1, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &request_in2);
      MPI_Isend(&(lu[1]), 1, MPI_DOUBLE, mpirank-1, 123, MPI_COMM_WORLD, &request_out2);
    }

    /* Jacobi step for all the inner points */
    for (i = 1; i < lN; i++){
      lunew[i]  = lu[i] +  omega * 0.5 * (hsq*rhs[i] + lu[i-1] + lu[i+1] - 2*lu[i]);
    }

    /* check if Isend/Irecv are done */
    if (mpirank < p - 1) {
      MPI_Wait(&request_out1, &status);
      MPI_Wait(&request_in1, &status);
      lunew[lN]  = lu[lN] +  omega * 0.5 * (hsq*rhs[lN] + lu[lN-1] + rr - 2*lu[lN]);
    }
    if (mpirank > 0) {
      MPI_Wait(&request_out2, &status);
      MPI_Wait(&request_in2, &status);
      lunew[0]  = lu[0] +  omega * 0.5 * (hsq*rhs[0] + ll + lu[1] - 2*lu[0]);
    }
     memcpy(lu, lunew, (lN+1)*sizeof(double));
  }
//  output_to_screen(lu,lN);
  free (lunew);
}

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
    /* get name of host running MPI process */
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);

  int i, Nfine, l, pN, iter, max_iters, levels, ssteps = 3;

  if ((argc < 3 || argc > 4) && mpirank == 0 ) {
    fprintf(stderr, "Usage: ./multigrid_1d Nfine maxiter [s-steps]\n");
    fprintf(stderr, "Nfine: # of intervals, Nfine/mpi-process_mpi must be power of two number\n");
    fprintf(stderr, "s-steps: # jacobi smoothing steps (optional, default is 3)\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  sscanf(argv[1], "%d", &Nfine);
  sscanf(argv[2], "%d", &max_iters);
  if (argc > 3)
    sscanf(argv[3], "%d", &ssteps);

    /* compute number of unknowns handled by each process */
  if ((Nfine % p != 0) && mpirank == 0 ) {
    printf("Nfine: %d, p: %d\n", Nfine, p);
    printf("Exiting. Nfine must be a multiple of p\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  pN = Nfine / p;
  /* compute number of multigrid levels */
  levels = floor(log2(pN));
  if (mpirank == 0 ) {
    printf("Multigrid Solve using V-cycles for -u'' = f on (0,1)\n");
    printf("Rank %d: Number of intervals = %d, max_iters = %d\n", mpirank,pN, max_iters);
    printf("Rank %d: Number of MG levels: %d \n", mpirank,levels);
  }
  /* timing */
  timestamp_type time1, time2;
  MPI_Barrier(MPI_COMM_WORLD);
  get_timestamp(&time1);

  /* Allocation of vectors, including left and right bdry points */
  double *u[levels], *rhs[levels];
  /* N, h*h and 1/(h*h) on each level */
  int *N = (int*) calloc(sizeof(int), levels);
  double *invhsq = (double* ) calloc(sizeof(double), levels);
  double *hsq = (double* ) calloc(sizeof(double), levels);
  double *res = (double *) calloc(sizeof(double), pN+1);
  for (l = 0; l < levels; ++l) {
    N[l] = pN / (int) pow(2,l);
    double h = 1.0 / (N[l]*p);
      if (mpirank == 0 ) {
    printf("l: %d: h=%f\n",l,h);
  }
    hsq[l] = h * h;
    if (mpirank == 0 ) {
    printf("MG level %2d, N = %8d\n", l, N[l]);
     }
    invhsq[l] = 1.0 / hsq[l];
    u[l]    = (double *) calloc(sizeof(double), N[l]+1);
    rhs[l] = (double *) calloc(sizeof(double), N[l]+1);
  }
  /* rhs on finest mesh */
  for (i = 0; i <= N[0]; ++i) {
    rhs[0][i] = 1.0;
  }
  /* set boundary values (unnecessary if calloc is used) */
  u[0][0] = u[0][N[0]] = 0.0;
  // if (mpirank == p - 1) u[0][N[0]] = 0.0;
  // if (mpirank == 0) u[0][0]= 0.0;
  double res_norm, res0_norm, tol = 1e-6;

  /* initial residual norm */
  compute_residual(u[0], rhs[0], res, N[0], invhsq[0]);
  res_norm = res0_norm = compute_norm(res, N[0]);
  if (mpirank == 0 ) {
  printf("Initial Residual: %f\n", res0_norm); 
  }
  for (iter = 0; iter < max_iters && res_norm/res0_norm > tol; iter++) {
    /* V-cycle: Coarsening */
    for (l = 0; l < levels-1; ++l) {
      /* pre-smoothing and coarsen */
      jacobi(u[l], rhs[l], N[l], hsq[l], ssteps);
      // output_to_screen(u[l],N[l],"jac");
      compute_and_coarsen_residual(u[l], rhs[l], rhs[l+1], N[l], invhsq[l]);
      // output_to_screen(rhs[l+1],N[l+1],"res");
      /* initialize correction for solution with zero */
      set_zero(u[l+1],N[l+1]);
    }
    /* V-cycle: Solve on coarsest grid using many smoothing steps */
    jacobi(u[levels-1], rhs[levels-1], N[levels-1], hsq[levels-1], 50);
    // output_to_screen(u[levels-1],N[levels-1],"jacDown");
    /* V-cycle: Refine and correct */
    for (l = levels-1; l > 0; --l) {
      /* refine and add to u */
            // output_to_screen(u[l-1],N[l-1],"refAddbefore");

      refine_and_add(u[l], u[l-1], N[l]);
      /* post-smoothing steps */
      // output_to_screen(u[l-1],N[l-1],"refAdd");
      // output_to_screen(rhs[l-1],N[l-1],"resBefore");
      jacobi(u[l-1], rhs[l-1], N[l-1], hsq[l-1], ssteps);
      
      // output_to_screen(u[l-1],N[l-1],"jacUp");
    }

    if (0 == (iter % 1)) {
      compute_residual(u[0], rhs[0], res, N[0], invhsq[0]);
      res_norm = compute_norm(res, N[0]);
      if (mpirank == 0 ) {
      printf("[Iter %d] Residual norm: %2.8f\n", iter, res_norm);
      }
    }
  }

  /* Clean up */
  free (hsq);
  free (invhsq);
  free (N);
  free(res);
  for (l = levels-1; l >= 0; --l) {
    free(u[l]);
    free(rhs[l]);
  }

  /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  if (mpirank == 0 ) {
    printf("Time elapsed is %f seconds.\n", elapsed);
  }
  MPI_Finalize();
  return 0;
}
