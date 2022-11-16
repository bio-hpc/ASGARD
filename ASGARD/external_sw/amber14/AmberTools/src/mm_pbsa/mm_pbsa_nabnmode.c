#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;

































 struct xmin_opt{INT_T mol_struct_opt;INT_T maxiter;REAL_T grms_tol;INT_T method;INT_T numdiff;INT_T m_lbfgs;INT_T iter;REAL_T xmin_time;INT_T ls_method;INT_T ls_maxiter;REAL_T ls_maxatmov;REAL_T beta_armijo;REAL_T c_armijo;REAL_T mu_armijo;REAL_T ftol_wolfe;REAL_T gtol_wolfe;INT_T ls_iter;INT_T print_level;INT_T error_flag;};





static  struct xmin_opt xo;

static MOLECULE_T *m;

static INT_T ier, natm;

static REAL_T *xyz,  *grad;

static REAL_T energy, grms;










int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static INT_T __gdab0001__;
static INT_T __gdab0002__;
static INT_T __it0001__;
static INT_T __it0002__;
static INT_T __it0003__;
static INT_T __it0004__;
static REAL_T __ft0001__;
static REAL_T __ft0002__;
m = getpdb( argv[2 - 1], NULL );
readparm( m, argv[3 - 1] );
natm =  *( NAB_mri( m, "natoms" ) );

__gdab0001__ = 3 * natm;DA_ALLOC( xyz = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "xyz" );
__gdab0002__ = 3 * natm;DA_ALLOC( grad = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "grad" );
setxyz_from_mol(  &m, NULL, xyz );


mm_options( argv[4 - 1] );
mme_init( m, NULL, "::ZZZZ", xyz, NULL );


xmin_opt_init(  &xo );
xo . maxiter = atoi( argv[5 - 1] );
xo . grms_tol = atof( argv[6 - 1] );
xo . method = 3;
xo . numdiff = 1;
xo . m_lbfgs = 3;
xo . ls_method = 2;
xo . ls_maxiter = 20;
xo . ls_maxatmov = 1.500000E-01;
xo . print_level = 1;

energy = mme( xyz, grad, ITEMP( __it0001__, 0 ) );
energy = xmin( mme,  &natm, xyz, grad,  &energy,  &grms,  &xo );
if( grms > xo . grms_tol )return(  - 3 );







ier = nmode( xyz, ITEMP( __it0001__, 3 * natm ), mme2, ITEMP( __it0002__, 0 ), ITEMP( __it0003__, 0 ), FTEMP( __ft0001__, 0.000000E+00 ), FTEMP( __ft0002__, 0.000000E+00 ), ITEMP( __it0004__, 0 ) );
return( ier );



	exit( 0 );
}
