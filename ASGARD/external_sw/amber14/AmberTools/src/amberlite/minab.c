#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;





















static MOLECULE_T *mol;

static STRING_T *outpdb = NULL;

static STRING_T *diel = NULL,  *gb = NULL,  *pdb = NULL,  *prmtop = NULL,  *aexp_move = NULL,  *aexp_restrain = NULL,  *mmopt = NULL;

static REAL_T *mol_xyz,  *reference_xyz;

static REAL_T *gradient,  *velocity;

static REAL_T cut, fret, en0, rgbmax, wcons, rmsgrad;

static INT_T nsnb, niter, ntpr, imin, ncycle;


int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static INT_T __gdab0001__;
static INT_T __gdab0002__;
static INT_T __gdab0003__;
static INT_T __gdab0004__;
static INT_T __it0001__;
static REAL_T __ft0001__;
static STRING_T *__st0001__ = NULL;
if(  !argv[2 - 1] ||  !argv[3 - 1] ||  !argv[4 - 1] ||  !argv[5 - 1] ||  !argv[6 - 1] ){
printf( "\n---------------------------------------------------------\n" );
printf( " minab version 0.7 (May 2010)\n" );
printf( "---------------------------------------------------------\n" );
printf( "version with no non-bonded list update, cutoff = 100 A, rgbmax = 15 A \n" );
printf( "usage: minab pdb prm outpdb gbflag niter ['restraints' resforce cutoff]\n" );
printf( "where: pdb        = PDB file with initial coordinates\n" );
printf( "       prm        = parameter-topology file name\n" );
printf( "       outpdb     = file name for refined coordinates (PDB format)\n" );
printf( "       gbflag     = integer (1, 2, 5, 7, or 8 for GB ON, else GB OFF)\n" );
printf( "       niter      = integer (number of iterations)\n" );
printf( "       restraints = atom expression for restrained atoms (':residues:atoms')\n" );
printf( "                    (must be included in 'quotes')!\n" );
printf( "       resforce   = force constant for restraints (kcal/mol/A2)\n" );
printf( "                    (must be given when restraints are specified!)\n" );
exit( 0 );
}
NAB_strcpy(  &pdb, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[2 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &prmtop, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[3 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &outpdb, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[4 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &gb, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[5 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );


niter = atoi( argv[6 - 1] );


aexp_move = NULL;

if(  !argv[7 - 1] ){
NAB_strcpy(  &aexp_restrain, "::ZZZZ" );
wcons = 0.000000E+00;
}
else{
NAB_strcpy(  &aexp_restrain, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[7 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
if(  !argv[8 - 1] ){
printf( "\nyou have specified atoms to be restrained...\n" );
printf( "the restraint force constant must also be given\n\n" );
exit( 0 );
}
wcons = atof( argv[8 - 1] );
}
mol = getpdb( pdb, NULL );readparm( mol, prmtop );
__gdab0001__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( mol_xyz = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "mol_xyz" );
__gdab0002__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( reference_xyz = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "reference_xyz" );
__gdab0003__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( gradient = ( REAL_T * )malloc( __gdab0003__ * ( sizeof( REAL_T ) ) ), "main", "gradient" );
setxyz_from_mol(  &mol, NULL, mol_xyz );
setxyz_from_mol(  &mol, NULL, reference_xyz );

if( ( EQ( gb, "1" ) ) || ( EQ( gb, "2" ) ) || ( EQ( gb, "5" ) ) || ( EQ( gb, "7" ) ) || ( EQ( gb, "8" ) ) ){NAB_strcpy(  &diel, "C" );}
else{NAB_strcpy(  &gb, "0" );NAB_strcpy(  &diel, "R" );}

cut = 100;rgbmax = 1.500000E+01;
nsnb = niter + 1;

NAB_strcpy(  &mmopt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "cut=%lf,nsnb=%d,diel=%s,gb=%s,rgbmax=%lf,wcons=%lf,ntpr=10", cut, nsnb, diel, gb, rgbmax, wcons ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( mmopt );
mme_init( mol, aexp_move, aexp_restrain, reference_xyz, NULL );
en0 = mme( mol_xyz, gradient, ITEMP( __it0001__, 1 ) );
rmsgrad = 1.000000E-01;
imin = conjgrad( mol_xyz, ITEMP( __it0001__, 3 *  *( NAB_mri( mol, "natoms" ) ) ),  &fret, mme,  &rmsgrad, FTEMP( __ft0001__, 1.000000E+01 ),  &niter );
printf( "------------------------------\ninitial energy: %8.3lf kcal/mol\n", en0 );
printf( "final   energy: %8.3lf kcal/mol\n", fret );
if( imin ==  - 1 ){printf( "minimizer stopped because of bad line search\n" );}
else if( imin ==  - 2 ){printf( "minimizer stopped because search direction was uphill\n" );}
else if( imin ==  - 3 ){printf( "minimizer stopped because number of iterations was exceeded\n" );}
else if( imin ==  - 4 ){printf( "minimizer stopped because function value could not be reduced further\n" );}
else printf( "minimizer finished after %d iterations\n", imin );
setmol_from_xyz(  &mol, NULL, mol_xyz );
putpdb( outpdb, mol, NULL );
printf( "refined coordinates written to %s\n", outpdb );
printf( "------------------------------\n" );
exit( 0 );


	exit( 0 );
}
