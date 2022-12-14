#include "defreal.h"
void	nrerror( char [] );
REAL_T  *vector( size_t, size_t );
int	*ivector( int, int );
int	*ipvector( int, int );
REAL_T	**matrix( int, int, int, int );
int	**imatrix( int, int, int, int );
void    free_vector( REAL_T *, size_t, size_t );
void	free_ivector( int *, int, int );
void	free_matrix( REAL_T **, int, int, int, int );
void	free_imatrix( int **, int, int, int, int );
