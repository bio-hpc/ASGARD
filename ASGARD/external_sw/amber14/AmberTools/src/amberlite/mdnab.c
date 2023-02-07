#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;





















static MOLECULE_T *mol;

static FILE_T *trajectory;

static STRING_T *diel = NULL,  *pdb = NULL,  *prmtop = NULL,  *traj = NULL;

static STRING_T *aexp_move = NULL,  *aexp_restrain = NULL;

static STRING_T *mmoptglobal = NULL,  *mmoptoutput = NULL,  *mmopttemp = NULL,  *mmopt = NULL,  *swcons = NULL;

static STRING_T *gb = NULL;

static REAL_T *mol_xyz,  *reference_xyz;

static REAL_T *gradient,  *velocity;

static REAL_T cut, dt, epsext, gamma_ln, tautp, temp0, tempi, rgbmax, wcons;

static INT_T nsnb, nsteps, ntpr, ntpr_md, ntwx, rattle, zerov, mdsteps;

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
static STRING_T *__st0001__ = NULL;
if(  !argv[2 - 1] ||  !argv[3 - 1] ||  !argv[4 - 1] ||  !argv[5 - 1] ||  !argv[6 - 1] ){
printf( "\n---------------------------------------------------------\n" );
printf( " mdnab version 1.0 (October 2009)\n" );
printf( "---------------------------------------------------------\n" );
printf( "usage: mdnab pdb prm traj gbflag picosecs ['restraints' resforce]\n\n" );
printf( "where: pdb        = PDB file name\n" );
printf( "       prm        = parameter-topology file name\n" );
printf( "       traj       = file name for trajectory (binpos format)\n" );
printf( "                    (the extension binpos is automatically added)\n" );
printf( "       gbflag     = integer (0 for GB OFF, 1, 2, 5, 7, or 8 for GB ON)\n" );
printf( "       picosecs   = integer (time of production phase)\n" );
printf( "                    (mdsteps = picosecs * 1000/2, because rattle is used)\n" );
printf( "       restraints = atom expression for restrained atoms ':residues:atoms'\n" );
printf( "                    (the expression must be included in 'quotes')!\n" );
printf( "       resforce   = force constant for restraints (kcal/mol/A2)\n" );
printf( "                    (must be given when restraints are specified!)\n" );
exit( 0 );
}
NAB_strcpy(  &pdb, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[2 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &prmtop, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[3 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &traj, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s.binpos", argv[4 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &gb, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[5 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mdsteps = 500 * ( atoi( argv[6 - 1] ) );
aexp_move = NULL;
if(  !argv[7 - 1] ){
NAB_strcpy(  &aexp_restrain, "::ZZZZ" );
wcons = 0.000000E+00;
}
else{
NAB_strcpy(  &aexp_restrain, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[7 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
if(  !argv[8 - 1] ){
printf( "\nyou have specified atoms to be restrained...\n" );exit( 0 );
printf( "the restraint force constant must also be given\n" );exit( 0 );
}
wcons = atof( argv[8 - 1] );
}
mol = getpdb( pdb, NULL );readparm( mol, prmtop );
__gdab0001__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( mol_xyz = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "mol_xyz" );
__gdab0002__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( reference_xyz = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "reference_xyz" );
__gdab0003__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( gradient = ( REAL_T * )malloc( __gdab0003__ * ( sizeof( REAL_T ) ) ), "main", "gradient" );
__gdab0004__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( velocity = ( REAL_T * )malloc( __gdab0004__ * ( sizeof( REAL_T ) ) ), "main", "velocity" );
setxyz_from_mol(  &mol, NULL, mol_xyz );
setxyz_from_mol(  &mol, NULL, reference_xyz );
trajectory = fopen( traj, "w" );
cut = 1.200000E+01;
rgbmax = cut;
nsnb = 25;
zerov = 0;
rattle = 1;dt = 2.000000E-03;
if( ( EQ( gb, "1" ) ) || ( EQ( gb, "2" ) ) || ( EQ( gb, "5" ) ) || ( EQ( gb, "7" ) ) || ( EQ( gb, "8" ) ) ){NAB_strcpy(  &diel, "C" );}
else{NAB_strcpy(  &gb, "0" );NAB_strcpy(  &diel, "R" );}



NAB_strcpy(  &mmoptglobal, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "cut=%lf,nsnb=%d,diel=%s,gb=%s,rgbmax=%lf,rattle=%d,dt=%lf", cut, nsnb, diel, gb, rgbmax, rattle, dt ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );





tempi = 50;temp0 = 100;gamma_ln = 20;
nsteps = 100;ntpr = nsteps + 1;ntpr_md = nsteps / 10;ntwx = 0;

NAB_strcpy(  &mmoptoutput, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ntpr=%d,ntpr_md=%d,ntwx=%d", ntpr, ntpr_md, ntwx ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

NAB_strcpy(  &mmopttemp, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "zerov=%d, tempi=%lf,temp0=%lf,gamma_ln=%lf,wcons=%lf", zerov, tempi, temp0, gamma_ln, wcons ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &mmopt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%s,%s", mmoptglobal, mmoptoutput, mmopttemp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( mmopt );
mme_init( mol, aexp_move, aexp_restrain, reference_xyz, NULL );
md( 3 *  *( NAB_mri( mol, "natoms" ) ), nsteps, mol_xyz, gradient, velocity, mme );


tempi = 0;temp0 = 150;gamma_ln = 20;
nsteps = 300;ntpr = nsteps + 1;ntpr_md = nsteps / 10;ntwx = 0;

NAB_strcpy(  &mmoptoutput, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ntpr=%d,ntpr_md=%d,ntwx=%d", ntpr, ntpr_md, ntwx ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

NAB_strcpy(  &mmopttemp, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "zerov=%d, tempi=%lf,temp0=%lf,gamma_ln=%lf,wcons=%lf", zerov, tempi, temp0, gamma_ln, wcons ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &mmopt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%s,%s", mmoptglobal, mmoptoutput, mmopttemp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( mmopt );
mme_init( mol, aexp_move, aexp_restrain, reference_xyz, NULL );
md( 3 *  *( NAB_mri( mol, "natoms" ) ), nsteps, mol_xyz, gradient, velocity, mme );


tempi = 0;temp0 = 200;gamma_ln = 10;
nsteps = 600;ntpr = nsteps + 1;ntpr_md = nsteps / 10;ntwx = 0;

NAB_strcpy(  &mmoptoutput, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ntpr=%d,ntpr_md=%d,ntwx=%d", ntpr, ntpr_md, ntwx ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

NAB_strcpy(  &mmopttemp, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "zerov=%d, tempi=%lf,temp0=%lf,gamma_ln=%lf,wcons=%lf", zerov, tempi, temp0, gamma_ln, wcons ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &mmopt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%s,%s", mmoptglobal, mmoptoutput, mmopttemp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( mmopt );
mme_init( mol, aexp_move, aexp_restrain, reference_xyz, NULL );
md( 3 *  *( NAB_mri( mol, "natoms" ) ), nsteps, mol_xyz, gradient, velocity, mme );


tempi = 0;temp0 = 250;gamma_ln = 5;
nsteps = 1000;ntpr = nsteps + 1;ntpr_md = nsteps / 10;ntwx = 0;

NAB_strcpy(  &mmoptoutput, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ntpr=%d,ntpr_md=%d,ntwx=%d", ntpr, ntpr_md, ntwx ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

NAB_strcpy(  &mmopttemp, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "zerov=%d, tempi=%lf,temp0=%lf,gamma_ln=%lf,wcons=%lf", zerov, tempi, temp0, gamma_ln, wcons ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &mmopt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%s,%s", mmoptglobal, mmoptoutput, mmopttemp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( mmopt );
mme_init( mol, aexp_move, aexp_restrain, reference_xyz, NULL );
md( 3 *  *( NAB_mri( mol, "natoms" ) ), nsteps, mol_xyz, gradient, velocity, mme );


tempi = 0;temp0 = 300;gamma_ln = 2;
nsteps = 3000;ntpr = nsteps + 1;ntpr_md = nsteps / 10;ntwx = 0;

NAB_strcpy(  &mmoptoutput, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ntpr=%d,ntpr_md=%d,ntwx=%d", ntpr, ntpr_md, ntwx ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

NAB_strcpy(  &mmopttemp, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "zerov=%d, tempi=%lf,temp0=%lf,gamma_ln=%lf,wcons=%lf", zerov, tempi, temp0, gamma_ln, wcons ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &mmopt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%s,%s", mmoptglobal, mmoptoutput, mmopttemp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( mmopt );
mme_init( mol, aexp_move, aexp_restrain, reference_xyz, NULL );
md( 3 *  *( NAB_mri( mol, "natoms" ) ), nsteps, mol_xyz, gradient, velocity, mme );


nsteps = 5000;ntpr = nsteps + 1;ntpr_md = nsteps / 10;ntwx = 0;

NAB_strcpy(  &mmoptoutput, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ntpr=%d,ntpr_md=%d,ntwx=%d", ntpr, ntpr_md, ntwx ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

NAB_strcpy(  &mmopttemp, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "zerov=%d, tempi=%lf,temp0=%lf,gamma_ln=%lf,wcons=%lf", zerov, tempi, temp0, gamma_ln, wcons ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &mmopt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%s,%s", mmoptglobal, mmoptoutput, mmopttemp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( mmopt );
mme_init( mol, aexp_move, aexp_restrain, reference_xyz, NULL );
md( 3 *  *( NAB_mri( mol, "natoms" ) ), nsteps, mol_xyz, gradient, velocity, mme );


nsteps = mdsteps;ntpr = nsteps + 1;ntpr_md = 500;ntwx = 500;

NAB_strcpy(  &mmoptoutput, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ntpr=%d,ntpr_md=%d,ntwx=%d", ntpr, ntpr_md, ntwx ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

NAB_strcpy(  &mmopttemp, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "t=0., zerov=%d, tempi=%lf,temp0=%lf,gamma_ln=%lf,wcons=%lf", zerov, tempi, temp0, gamma_ln, wcons ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &mmopt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%s,%s", mmoptglobal, mmoptoutput, mmopttemp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( mmopt );
mme_init( mol, aexp_move, aexp_restrain, reference_xyz, trajectory );
md( 3 *  *( NAB_mri( mol, "natoms" ) ), nsteps, mol_xyz, gradient, velocity, mme );
printf( "\ntrajectory with %s picoseconds was written to %s...\n\n", argv[6 - 1], traj );


	exit( 0 );
}
