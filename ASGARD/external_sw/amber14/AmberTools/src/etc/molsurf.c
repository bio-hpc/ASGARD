#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;

static MOLECULE_T *m;

static REAL_T area;

int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static REAL_T __ft0001__;
if( argc != 3 ){
fprintf( stderr, "Usage: molsurf <pqr-file> <probe-radius>\n" );
exit( 1 );
}
m = getpdb( argv[2 - 1], "-pqr" );
area = molsurf(  &m, NULL, FTEMP( __ft0001__, atof( argv[3 - 1] ) ) );
printf( "surface area = %10.3f\n", area );


	exit( 0 );
}
