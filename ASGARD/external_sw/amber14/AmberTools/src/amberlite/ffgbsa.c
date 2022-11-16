#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;





















static ATOM_T *a;
static MOLECULE_T *mol;

static REAL_T epot,  *gradient,  *mol_xyz, dummy[3], sasa, proberad;

static STRING_T *molmec_opt = NULL,  *pdb_file = NULL,  *prmtop_file = NULL,  *gb = NULL,  *rgbmax = NULL,  *cutoff = NULL,  *diel = NULL,  *sa = NULL;

int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static INT_T __gdab0001__;
static INT_T __gdab0002__;
static INT_T __it0001__;
static STRING_T *__st0001__ = NULL;
if(  !argv[2 - 1] ||  !argv[3 - 1] ||  !argv[4 - 1] ||  !argv[5 - 1] ){
printf( "\n---------------------------------------------------------\n" );
printf( " ffgbsa version 1.1 (September 2009)\n" );
printf( "---------------------------------------------------------\n" );
printf( "usage: ffgbsa pdb prm gbflag saflag\n" );
printf( "where: pdb        = PDB file name\n" );
printf( "       prm        = parameter-topology file name\n" );
printf( "       gbflag     = integer (1, 2, 5, 7 or 8 for GB ON, else OFF)\n" );
printf( "       saflag     = integer (0 for SA OFF, 1 for SA ON)\n\n" );
exit( 0 );
}
NAB_strcpy(  &pdb_file, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[2 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );NAB_strcpy(  &prmtop_file, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[3 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &gb, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[4 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );NAB_strcpy(  &sa, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s", argv[5 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
if( ( EQ( gb, "1" ) ) || ( EQ( gb, "2" ) ) || ( EQ( gb, "5" ) ) || ( EQ( gb, "7" ) ) || ( EQ( gb, "8" ) ) ){NAB_strcpy(  &diel, "C" );}else{NAB_strcpy(  &gb, "0" );NAB_strcpy(  &diel, "R" );}
NAB_strcpy(  &molmec_opt, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "cut=100, rgbmax=100 ,diel=%s, gb=%s", diel, gb ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mol = getpdb( pdb_file, NULL );readparm( mol, prmtop_file );
__gdab0002__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( mol_xyz = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "mol_xyz" );__gdab0001__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( gradient = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "gradient" );
setxyz_from_mol(  &mol, NULL, mol_xyz );
mm_options( molmec_opt );mme_init( mol, "::", "::ZZZ", dummy, NULL );
epot = mme( mol_xyz, gradient, ITEMP( __it0001__, 0 ) );
if( NE( sa, "0" ) ){

for( a = NULL;a = NAB_mnext( mol, a ); ){
if( EQ( substr( a->a_atomname, 1, 1 ), "H" ) ){a->a_radius = 1.200000E+00 + 1.400000E+00;}
else if( EQ( substr( a->a_atomname, 1, 1 ), "N" ) ){a->a_radius = 1.550000E+00 + 1.400000E+00;}
else if( EQ( substr( a->a_atomname, 1, 1 ), "O" ) ){a->a_radius = 1.500000E+00 + 1.400000E+00;}
else if( EQ( substr( a->a_atomname, 1, 1 ), "C" ) ){a->a_radius = 1.700000E+00 + 1.400000E+00;}
else if( EQ( substr( a->a_atomname, 1, 1 ), "S" ) ){a->a_radius = 1.800000E+00 + 1.400000E+00;}
else if( EQ( substr( a->a_atomname, 1, 1 ), "F" ) ){a->a_radius = 1.470000E+00 + 1.400000E+00;}
else if( EQ( substr( a->a_atomname, 1, 1 ), "P" ) ){a->a_radius = 1.800000E+00 + 1.400000E+00;}
else if( ( EQ( substr( a->a_atomname, 1, 2 ), "Br" ) ) || ( EQ( substr( a->a_atomname, 1, 2 ), "BR" ) ) )
{a->a_radius = 1.850000E+00 + 1.400000E+00;}
else if( ( EQ( substr( a->a_atomname, 1, 2 ), "Cl" ) ) || ( EQ( substr( a->a_atomname, 1, 2 ), "CL" ) ) )
{a->a_radius = 1.750000E+00 + 1.400000E+00;}
else if( EQ( substr( a->a_atomname, 1, 1 ), "I" ) ){a->a_radius = 1.980000E+00 + 1.400000E+00;}

else{a->a_radius = 1.500000E+00 + 1.400000E+00;}
}

proberad = 0.000000E+00;
sasa = molsurf(  &mol, STEMP( __st0001__, "::" ),  &proberad );printf( "sasa: %10.2lf\n", sasa );
printf( "Esasa = 0.0072 * sasa = %10.2lf\n", sasa * 7.200000E-03 );
}


	exit( 0 );
}
