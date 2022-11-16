#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;






























 struct AmberNetcdf{REAL_T temp0;REAL_T restartTime;INT_T isNCrestart;INT_T ncid;INT_T frameDID;INT_T ncframe;INT_T currentFrame;INT_T atomDID;INT_T ncatom;INT_T ncatom3;INT_T coordVID;INT_T velocityVID;INT_T cellAngleVID;INT_T cellLengthVID;INT_T spatialDID;INT_T labelDID;INT_T cell_spatialDID;INT_T cell_angularDID;INT_T spatialVID;INT_T timeVID;INT_T cell_spatialVID;INT_T cell_angularVID;INT_T TempVID;};



























 struct xmin_opt{INT_T mol_struct_opt;INT_T maxiter;REAL_T grms_tol;INT_T method;INT_T numdiff;INT_T m_lbfgs;INT_T iter;REAL_T xmin_time;INT_T ls_method;INT_T ls_maxiter;REAL_T ls_maxatmov;REAL_T beta_armijo;REAL_T c_armijo;REAL_T mu_armijo;REAL_T ftol_wolfe;REAL_T gtol_wolfe;INT_T ls_iter;INT_T print_level;INT_T error_flag;};







static  struct xmin_opt xo;

static  struct AmberNetcdf nc_ascii;

static MOLECULE_T *mol;

static INT_T natm, nmode_results, iter;

static INT_T framenum, traj_complete;

static INT_T i, j;

static REAL_T *xyz,  *grad;

static REAL_T energy, min_energy;

static FILE_T *asciitraj;






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
if( argc != 7 ){
printf( "    Bad command-line arguments\n" );
printf( "Usage: %s {pdb} {prmtop} {maxiter} {drms} {'string of MM options'} {nc/trj}\n", argv[1 - 1] );
exit( 1 );
}



mol = getpdb( argv[2 - 1], NULL );
readparm( mol, argv[3 - 1] );
natm =  *( NAB_mri( mol, "natoms" ) );
__gdab0001__ = 3 * natm;DA_ALLOC( xyz = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "xyz" );
__gdab0002__ = 3 * natm;DA_ALLOC( grad = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "grad" );
setxyz_from_mol(  &mol, NULL, xyz );


mm_options( argv[6 - 1] );
mme_init( mol, NULL, "::Z", xyz, NULL );







xmin_opt_init(  &xo );
xo . maxiter = atoi( argv[4 - 1] );
xo . grms_tol = atof( argv[5 - 1] );
xo . method = 3;
xo . numdiff = 1;
xo . m_lbfgs = 3;
xo . ls_method = 2;
xo . ls_maxiter = 20;
xo . ls_maxatmov = 1.500000E-01;
xo . print_level = 2;



framenum = 0;
if( ( netcdfLoad(  &nc_ascii, argv[7 - 1] ) ) == 0 ){
netcdfLoad(  &nc_ascii, argv[7 - 1] );
printf( "\n Processing NETCDF traj\n\n" );
while( netcdfGetNextFrame(  &nc_ascii, xyz, NULL, NULL ) ){
energy = mme( xyz, grad,  &framenum );
energy = xmin( mme,  &natm, xyz, grad,  &energy,  &min_energy,  &xo );
if( min_energy < ( atof( argv[5 - 1] ) ) ){
printf( "     ----Convergence Satisfied---- \n\n" );
nmode( xyz, ITEMP( __it0001__, 3 * natm ), mme2, ITEMP( __it0002__, 0 ), ITEMP( __it0003__, 1 ), FTEMP( __ft0001__, 0.000000E+00 ), FTEMP( __ft0002__, 0.000000E+00 ), ITEMP( __it0004__, 0 ) );
framenum ++ ;
}
else{
printf( " \n    |----Not minimized properly!----|\n " );
printf( "   |---- Entropy not Calculated---|\n\n" );
}
}
netcdfClose(  &nc_ascii );
}
else{
printf( "\n Processing ASCII traj\n\n" );
asciitraj = fopen( argv[7 - 1], "r" );
if( asciitraj == NULL ){
printf( "\n Unable to open mdcrd (%s) !\n\n", argv[7 - 1] );
exit( 1 );
}
traj_complete = 0;
NAB_getline( asciitraj );
for( i = 1;;i ++  ){
for( j = 1;j <= 3 *  *( NAB_mri( mol, "natoms" ) );j ++  ){
if( ( fscanf( asciitraj, "%lf",  &xyz[j - 1] ) ) < 1 ){
traj_complete = 1;
break;
}
}
if( traj_complete )break;
energy = mme( xyz, grad,  &framenum );
energy = xmin( mme,  &natm, xyz, grad,  &energy,  &min_energy,  &xo );
if( min_energy < ( atof( argv[5 - 1] ) ) ){
printf( "     ----Convergence Satisfied---- \n\n" );
nmode( xyz, ITEMP( __it0001__, 3 * natm ), mme2, ITEMP( __it0002__, 0 ), ITEMP( __it0003__, 1 ), FTEMP( __ft0001__, 0.000000E+00 ), FTEMP( __ft0002__, 0.000000E+00 ), ITEMP( __it0004__, 0 ) );
framenum ++ ;
}
else{
printf( " \n    |----Not minimized properly!----|\n " );
printf( "   |---- Entropy not Calculated---|\n\n" );
}
}
}



	exit( 0 );
}
