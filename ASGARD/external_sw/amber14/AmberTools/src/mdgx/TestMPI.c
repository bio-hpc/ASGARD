#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
 
int main(int argc, char* argv[])
{
  MPI_Comm dup_comm_world, world_comm;
  MPI_Group world_group;
  int world_rank, world_size, rank;
 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  MPI_Comm_dup( MPI_COMM_WORLD, &dup_comm_world );

  /* Exercise Comm_create by creating an equivalent */
  /* to dup_comm_world (sans attributes)            */
  MPI_Comm_group( dup_comm_world, &world_group );
  MPI_Comm_create( dup_comm_world, world_group, &world_comm );

  printf("Created com %12d on process %d\n", (int)world_comm, world_rank);
  system("sleep 0.2");

  MPI_Comm_rank( world_comm, &rank );
  if (rank != world_rank) {
    printf( "incorrect rank in world comm: %d\n", rank );
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 3001);
  }

  int i, j;
  int* myints;
  myints = (int*)malloc(world_size*sizeof(int));
  for (i = 0; i < world_size; i++) {
    myints[i] = world_rank;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for (i = 0; i < world_size; i++) {
    if (i == world_rank) {
      printf("myints on %d: [ ", i);
      for (j = 0; j < world_size; j++) {
	printf("%4d", myints[j]);
      }
      printf("];\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  system("sleep 0.2");
  if (world_rank == 0) {
    printf("Now, let's try a 3 -> 1 message.\n");
  }

  /*** Post a message from 3 -> 1 to copy some integers ***/
  MPI_Status req[2];
  if (world_rank == 1) {
    MPI_Recv(myints, world_size, MPI_INT, 3, 475, world_comm, &req[0]);
  }
  if (world_rank == 3) {
    MPI_Send(myints, world_size, MPI_INT, 1, 475, world_comm);
  }

  for (i = 0; i < world_size; i++) {
    if (i == world_rank) {
      printf("myints on %d: [ ", i);
      for (j = 0; j < world_size; j++) {
	printf("%4d", myints[j]);
      }
      printf("];\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  system("sleep 0.2");
  if (world_rank == 0) {
    printf("Now, let's try a 2 -> 0 message.\n");
  }

  /*** Post a message from 2 -> 0 to copy some integers ***/
  if (world_rank == 0) {
    MPI_Recv(myints, world_size, MPI_INT, 2, 475, world_comm, &req[0]);
  }
  if (world_rank == 2) {
    MPI_Send(myints, world_size, MPI_INT, 0, 475, world_comm);
  }

  for (i = 0; i < world_size; i++) {
    if (i == world_rank) {
      printf("myints on %d: [ ", i);
      for (j = 0; j < world_size; j++) {
	printf("%4d", myints[j]);
      }
      printf("];\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  system("sleep 0.2");

  MPI_Finalize();
  return 0;
}
 

