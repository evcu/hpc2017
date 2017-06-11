/* Communication ping-pong:
 * Exchange between messages between mpirank
 * 0 <-> 1, 2 <-> 3, ....
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#define PART2ONLY

int part1(int round){
  int rank,size, tag, origin, destination,i;
  MPI_Status status;

  char hostname[1024];
  gethostname(hostname, 1024);

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int message_out = 0;
  int message_in = 0;
  tag = 99;
  
  for(i=0;i<round;i++){
    destination = (rank + 1)%size;
    origin = (rank-1+size)%size ;
    if(rank==0){
      message_out = message_in;
      MPI_Send(&message_out, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
      MPI_Recv(&message_in,  1, MPI_INT, origin,      tag, MPI_COMM_WORLD, &status);
      printf("round %05d, rank %d hosted on %s received from %d the message %d\n",i+1, rank, hostname, origin, message_in);
    }
    else{
      MPI_Recv(&message_in,  1, MPI_INT, origin,      tag, MPI_COMM_WORLD, &status);
      printf("round %05d, rank %d hosted on %s received from %d the message %d\n",i, rank, hostname, origin, message_in);
      message_out = message_in + rank;
      MPI_Send(&message_out, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
    }
    
  }
  MPI_Finalize();
  return 0;
}

int part2(int round){
	int rank,size, tag, origin, destination,i;
  MPI_Status status;

  char hostname[1024];
  gethostname(hostname, 1024);

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int N = 500000;

  long *message_out = calloc(N, sizeof(long));
  long *message_in = calloc(N, sizeof(long));
  tag = 99;
  
  for(i=0;i<round;i++){
    destination = (rank + 1)%size;
    origin = (rank-1+size)%size ;
    if(rank==0){
      MPI_Send(message_out, N, MPI_LONG, destination, tag, MPI_COMM_WORLD);
      MPI_Recv(message_in,  N, MPI_LONG, origin,      tag, MPI_COMM_WORLD, &status);
      printf("round %05d, rank %d hosted on %s received from %d the array of size %d bytes\n",i+1, rank, hostname, origin, N*4);
    }
    else{
      MPI_Recv(message_in,  N, MPI_LONG, origin,      tag, MPI_COMM_WORLD, &status);
      printf("round %05d, rank %d hosted on %s received from %d the array of size %d bytes\n",i, rank, hostname, origin, N*4);
      MPI_Send(message_out, N, MPI_LONG, destination, tag, MPI_COMM_WORLD);
    }
    
  }
  free( message_out );
  free( message_in );
  MPI_Finalize();
  return 0;
}

int main( int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
    int round = 0;
  sscanf(argv[1], "%d", &round);
  #ifdef PART2ONLY
    return part2(round);
  #else
    return part1(round);
  #endif
  

 
}
