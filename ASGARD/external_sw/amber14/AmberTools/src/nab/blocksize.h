/*
 * Define the number of contiguous loop indices mapped to each
 * OpenMP thread or MPI task.  This parameter is the log (base 2)
 * of the Chunk Size as defined for the OpenMP schedule clause.
 * For example, a value of 3 equates to a Chunk Size of 8.  This
 * parameter is used in nblist, nbond, egb and egb2.
 */

#define BLOCKSIZE (32)

