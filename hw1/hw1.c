/* multiply vector components, write into a vector,
 *  and compute the inner product  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

int main (int argc, char **argv)
{
  long  i, n;
  double *u,*f,*tempv;
  double prod,h2;
  long p, passes;
  double tolerance=1e4,threshold;
  if (argc != 3) {
    fprintf(stderr, "Function needs vector size and number of passes as input arguments!\n");
    abort();
  }
  n = atol(argv[1]);
  passes = atol(argv[2]);
  h2 = pow(1.0/(n+1),2.0);
  printf("h2 is %f\n",h2);
  u = (double *) malloc(sizeof(double) * n);
  f = (double *) malloc(sizeof(double) * n);
  tempv = (double *) malloc(sizeof(double) * n);

  /* fill vectors */
  for (i = 0; i < n; ++i) {
    u[i] = 0;
    f[i] = 1*h2;
  }

  timestamp_type time1, time2;
  get_timestamp(&time1);

  for (p = 0; p < passes; ++p)
    {
    	for (i = 0; i < n; ++i) {
    		tempv[i] = u[i];
    	}
      	for (i = 1; i < n-1; ++i) {
	      	u[i] = (f[i] + (tempv[i+1]+tempv[i-1]))/2;
      	}
		u[0] = (f[0] +tempv[1])/2;
		u[n-1] = (f[n-1] +tempv[n-2])/2;
	
      /*eval*/
      prod = 0;
      for (i = 0; i < n; ++i) {
  		prod += pow((f[i] - (2*u[i]-u[i+1]-u[i-1])),2.0);
    	}
    prod += pow((f[0] - (2*u[0]-u[1])),2.0);
	prod += pow((f[n-1] - (2*u[n-1]-u[n-2])),2.0);
    prod = sqrt(prod);
    printf("Iteration %ld: norm of residual is %.20f\n",p+1,prod);
    if (p == 0){
    	threshold = prod/tolerance;
    }
    else if (threshold>prod){
    	printf("The tolerance value exceeded at iteration: %ld\n",p+1);
    	goto early_stop;
    }
 }
printf("Residual decrese factor is %f\n",(threshold*tolerance)/prod);

early_stop:
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("Time elapsed is %f seconds.\n", elapsed);

  printf("%f GB/s\n", 4*n*sizeof(double)*passes/1e9/elapsed);
  printf("%f GFlops/s\n", 2*n*passes/1e9/elapsed);

  free(u);
  free(f);
  free(tempv);
  return 0;
}
