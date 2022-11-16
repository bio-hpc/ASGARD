#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;

































 struct AmberNetcdf{REAL_T temp0;REAL_T restartTime;INT_T isNCrestart;INT_T ncid;INT_T frameDID;INT_T ncframe;INT_T currentFrame;INT_T atomDID;INT_T ncatom;INT_T ncatom3;INT_T coordVID;INT_T velocityVID;INT_T cellAngleVID;INT_T cellLengthVID;INT_T spatialDID;INT_T labelDID;INT_T cell_spatialDID;INT_T cell_angularDID;INT_T spatialVID;INT_T timeVID;INT_T cell_spatialVID;INT_T cell_angularVID;INT_T TempVID;};




static STRING_T *name = NULL;



 struct mdOptions{STRING_T *pdb,  *prmtop,  *rst,  *traj;};

















 struct rismOptions{STRING_T *xvv,  *guv,  *cuv,  *huv,  *uuv,  *asymp,  *quv,  *chgdist,  *exchem,  *solvene,  *entropy,  *exchemGF,  *solveneGF,  *entropyGF,  *exchemUC,  *solveneUC,  *entropyUC,  *potUV,  *volfmt;INT_T closureOrder, asympcorr;REAL_T buffer, solvcut;REAL_T grdspcX, grdspcY, grdspcZ, solvboxX, solvboxY, solvboxZ;INT_T ngX, ngY, ngZ;REAL_T mdiis_del, mdiis_restart;INT_T mdiis_nvec, mdiis_method;INT_T maxstep, npropagate;INT_T centering, zerofrc, apply_rism_force, polarDecomp, entropicDecomp;INT_T ntwrism, verbose, progress, saveprogress;REAL_T uccoeff1, uccoeff2;};





static  struct mdOptions mdOpt;





static  struct rismOptions rismOpt;


static STRING_T * *closure;

static REAL_T *tolerance;








INT_T error( STRING_T * *message );
INT_T error( STRING_T * *message ){
if( mytaskid == 0 ){
return( fprintf( stderr, "ERROR:  %s\n",  *message ) );
}
return(  - 1 );
}







INT_T setDefaults(  struct mdOptions *mdOpt,  struct rismOptions *rismOpt );
INT_T setDefaults(  struct mdOptions *mdOpt,  struct rismOptions *rismOpt ){
(  *mdOpt ) . pdb = NULL;
(  *mdOpt ) . prmtop = NULL;
(  *mdOpt ) . rst = NULL;
(  *mdOpt ) . traj = NULL;
(  *rismOpt ) . xvv = NULL;
(  *rismOpt ) . guv = NULL;
(  *rismOpt ) . huv = NULL;
(  *rismOpt ) . cuv = NULL;
(  *rismOpt ) . uuv = NULL;
(  *rismOpt ) . asymp = NULL;
(  *rismOpt ) . quv = NULL;
(  *rismOpt ) . chgdist = NULL;
(  *rismOpt ) . exchem = NULL;
(  *rismOpt ) . solvene = NULL;
(  *rismOpt ) . entropy = NULL;
(  *rismOpt ) . exchemGF = NULL;
(  *rismOpt ) . solveneGF = NULL;
(  *rismOpt ) . entropyGF = NULL;
(  *rismOpt ) . exchemUC = NULL;
(  *rismOpt ) . solveneUC = NULL;
(  *rismOpt ) . entropyUC = NULL;
(  *rismOpt ) . potUV = NULL;
NAB_strcpy(  &(  *rismOpt ) . volfmt, "dx" );
NAB_strcpy(  &closure[1 - 1], "kh" );
(  *rismOpt ) . closureOrder = 1;
(  *rismOpt ) . asympcorr = 1;
(  *rismOpt ) . buffer = 14;
(  *rismOpt ) . solvcut =  - 1;
(  *rismOpt ) . grdspcX = 5.000000E-01;
(  *rismOpt ) . grdspcY = 5.000000E-01;
(  *rismOpt ) . grdspcZ = 5.000000E-01;
(  *rismOpt ) . ngX = 0;
(  *rismOpt ) . ngY = 0;
(  *rismOpt ) . ngZ = 0;
(  *rismOpt ) . solvboxX = 0;
(  *rismOpt ) . solvboxY = 0;
(  *rismOpt ) . solvboxZ = 0;
tolerance[1 - 1] = 1.000000E-05;
(  *rismOpt ) . mdiis_del = 7.000000E-01;
(  *rismOpt ) . mdiis_restart = 10;
(  *rismOpt ) . mdiis_nvec = 5;
(  *rismOpt ) . mdiis_method = 2;
(  *rismOpt ) . maxstep = 10000;
(  *rismOpt ) . npropagate = 5;
(  *rismOpt ) . centering = 1;
(  *rismOpt ) . zerofrc = 1;
(  *rismOpt ) . apply_rism_force = 0;
(  *rismOpt ) . polarDecomp = 0;
(  *rismOpt ) . entropicDecomp = 0;
(  *rismOpt ) . ntwrism = 0;
(  *rismOpt ) . verbose = 0;
(  *rismOpt ) . progress = 1;
(  *rismOpt ) . saveprogress = 0;
(  *rismOpt ) . uccoeff1 =  - 3.312000E+00;
(  *rismOpt ) . uccoeff2 = 1.152000E+00;
return( 0 );
}







INT_T printOptions(  struct mdOptions *mdOpt,  struct rismOptions *rismOpt );
INT_T printOptions(  struct mdOptions *mdOpt,  struct rismOptions *rismOpt ){
setDefaults( mdOpt, rismOpt );

fprintf( stderr, "\n" );
fprintf( stderr, "%-20s %s\n", "Key", "Default" );
fprintf( stderr, "%-20s %s\n", "---", "-------" );

fprintf( stderr, "%-20s %s\n", "--pdb", (  *mdOpt ) . pdb );
fprintf( stderr, "%-20s %s\n", "--prmtop", (  *mdOpt ) . prmtop );
fprintf( stderr, "%-20s %s\n", "--rst", (  *mdOpt ) . rst );
fprintf( stderr, "%-20s %s\n", "-y|--traj", (  *mdOpt ) . traj );
fprintf( stderr, "%-20s %s\n", "--xvv", (  *rismOpt ) . xvv );
fprintf( stderr, "%-20s %s\n", "--guv", (  *rismOpt ) . guv );
fprintf( stderr, "%-20s %s\n", "--huv", (  *rismOpt ) . huv );
fprintf( stderr, "%-20s %s\n", "--cuv", (  *rismOpt ) . cuv );
fprintf( stderr, "%-20s %s\n", "--uuv", (  *rismOpt ) . uuv );
fprintf( stderr, "%-20s %s\n", "--asymp", (  *rismOpt ) . asymp );
fprintf( stderr, "%-20s %s\n", "--quv", (  *rismOpt ) . quv );
fprintf( stderr, "%-20s %s\n", "--chgdist", (  *rismOpt ) . chgdist );










fprintf( stderr, "%-20s %s\n", "--volfmt", (  *rismOpt ) . volfmt );
fprintf( stderr, "%-20s %s\n", "--closure", closure[1 - 1] );
fprintf( stderr, "%-20s %d\n", "--closureorder", (  *rismOpt ) . closureOrder );
fprintf( stderr, "%-20s %i\n", "--asympcorr", (  *rismOpt ) . asympcorr );
fprintf( stderr, "%-20s %f\n", "--buffer", (  *rismOpt ) . buffer );
fprintf( stderr, "%-20s %f\n", "--solvcut", (  *rismOpt ) . solvcut );

fprintf( stderr, "%-20s %f,%f,%f\n", "--grdspc", (  *rismOpt ) . grdspcX, (  *rismOpt ) . grdspcY, (  *rismOpt ) . grdspcZ );

fprintf( stderr, "%-20s %d,%d,%d\n", "--ng", (  *rismOpt ) . ngX, (  *rismOpt ) . ngY, (  *rismOpt ) . ngZ );

fprintf( stderr, "%-20s %f,%f,%f\n", "--solvbox", (  *rismOpt ) . solvboxX, (  *rismOpt ) . solvboxY, (  *rismOpt ) . solvboxZ );
fprintf( stderr, "%-20s %f\n", "--tolerance", tolerance[1 - 1] );
fprintf( stderr, "%-20s %f\n", "--mdiis_del", (  *rismOpt ) . mdiis_del );
fprintf( stderr, "%-20s %f\n", "--mdiis_restart", (  *rismOpt ) . mdiis_restart );
fprintf( stderr, "%-20s %d\n", "--mdiis_nvec", (  *rismOpt ) . mdiis_nvec );

fprintf( stderr, "%-20s %d\n", "--maxstep", (  *rismOpt ) . maxstep );
fprintf( stderr, "%-20s %d\n", "--npropagate", (  *rismOpt ) . npropagate );
fprintf( stderr, "%-20s %d\n", "--centering", (  *rismOpt ) . centering );


fprintf( stderr, "%-20s %d\n", "--polarDecomp", (  *rismOpt ) . polarDecomp );


fprintf( stderr, "%-20s %d\n", "--verbose", (  *rismOpt ) . verbose );


return( 0 );
}







INT_T usage(  );
INT_T usage(  ){
STRING_T *empty = NULL;

INT_T i;

if( mytaskid == 0 ){
fprintf( stderr, "USAGE: %s --pdb pdbfile --prmtop prmtopfile [--rst rstfile] [-y|--traj netCDFfile]\n", name );
for( i = 0;i < ( length( name ) );i ++  ){
NAB_strcpy(  &empty, NAB_strcat( empty, " " ) );
}
fprintf( stderr, "       %s --xvv Xvv_filename [--guv Guv_rootname]\n", empty );
fprintf( stderr, "       %s [--cuv Cuv_rootname] [--huv Huv_rootname]\n", empty );
fprintf( stderr, "       %s [--uuv Uuv_rootname] [--asymp asymp_rootname]\n", empty );
fprintf( stderr, "       %s [--quv Quv_rootname] [--chgdist chgdist_rootname]\n", empty );





fprintf( stderr, "       %s [--volfmt volume_format]\n", empty );
fprintf( stderr, "       %s [--closure kh|hnc|pse(1|2|...)[ closure2[ ...]]]\n", empty );
fprintf( stderr, "       %s [--[no]asympcorr] [--buffer distance] [--solvcut distance]\n", empty );
fprintf( stderr, "       %s [--grdspc dx,dy,dz] [--ng nx,ny,nz] [-solvbox lx,ly,lz]\n", empty );
fprintf( stderr, "       %s [--tolerance tol1[ tol2 [...]] [--mdiis_del step_size]\n", empty );
fprintf( stderr, "       %s [--mdiis_restart threshold] [--mdiis_nvec vectors]\n", empty );
fprintf( stderr, "       %s [--maxstep 1000] [--npropagate #_old_solutions]\n", empty );

fprintf( stderr, "       %s [--[no]polarDecomp] [--centering -4..4]\n", empty );


fprintf( stderr, "       %s [--verbose 0|1|2]\n", empty );
printOptions(  &mdOpt,  &rismOpt );
}
exit( 1 );
return( 0 );
}










INT_T testValue( STRING_T * *key, STRING_T * *value );
INT_T testValue( STRING_T * *key, STRING_T * *value ){
STRING_T *__st0001__ = NULL;
STRING_T *__st0002__ = NULL;
if( ( NAB_rematch(  *value, "^-" ) ) &&  !( NAB_rematch(  *value, "^-[0-9]*$" ) ) ){
error( STEMP( __st0002__, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "'%s' must be followed by a value.\n",  *key ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) ) );
usage(  );
}
return( 0 );
}









INT_T numCSV( STRING_T * *str );
INT_T numCSV( STRING_T * *str ){
INT_T pos0, pos1, num;

num = 0;
pos0 = 0;
pos1 = 0;
while( ( pos0 = NAB_index( substr(  *str, pos1 + 1, length(  *str ) ), "," ) ) != 0 ){
num ++ ;
pos1 += pos0;
if( pos1 + 1 > ( length(  *str ) ) ){break;}
}
return(  ++ num );
}







INT_T checkOptions(  struct mdOptions *mdOpt,  struct rismOptions *rismOpt );
INT_T checkOptions(  struct mdOptions *mdOpt,  struct rismOptions *rismOpt ){
STRING_T *__st0001__ = NULL;
if( EQ( (  *mdOpt ) . pdb, NULL ) ){
error( STEMP( __st0001__, "a PDB file is required\n" ) );
usage(  );
}
if( EQ( (  *mdOpt ) . prmtop, NULL ) ){
error( STEMP( __st0001__, "a PRMTOP file is required\n" ) );
usage(  );
}
if( EQ( (  *rismOpt ) . xvv, NULL ) ){
error( STEMP( __st0001__, "an XVV file is required\n" ) );
usage(  );
}
return( 0 );
}




static MOLECULE_T *m;

static REAL_T *p_xyz,  *f_xyz,  *v_xyz;

static  struct AmberNetcdf nc;

static REAL_T time, dgrad, fret;

static INT_T i, j;

static STRING_T *value = NULL,  *key = NULL,  *tmpstr = NULL;

static INT_T tempInt;

static REAL_T ftmp;

static FILE_T *asciiTraj;

static INT_T iframe, ipos, trajDone;

static INT_T nclosure, ntolerance;


int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static INT_T __gdab0001__;
static INT_T __gdab0002__;
static INT_T __gdab0003__;
static INT_T __it0001__;
static STRING_T *__st0001__ = NULL;
static STRING_T *__st0002__ = NULL;
NAB_strcpy(  &name, argv[1 - 1] );


__gdab0001__ = 1;DA_ALLOC( closure = ( STRING_T * * )calloc( __gdab0001__, sizeof( STRING_T * ) ), "main", "closure" );
nclosure = 1;

__gdab0002__ = 1;DA_ALLOC( tolerance = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "tolerance" );
ntolerance = 1;
setDefaults(  &mdOpt,  &rismOpt );





i = 2;
while( i <= argc ){
NAB_strcpy(  &key, argv[i - 1] );
NAB_strcpy(  &value, "-" );
if( EQ( key, "--printDefaults" ) ){
setDefaults(  &mdOpt,  &rismOpt );
printOptions(  &mdOpt,  &rismOpt );
exit( 0 );
}else if( EQ( key, "--pdb" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &mdOpt . pdb, value );
}else if( EQ( key, "--rst" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &mdOpt . rst, value );
}else if( ( EQ( key, "--traj" ) ) || ( EQ( key, "-y" ) ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &mdOpt . traj, value );
}else if( EQ( key, "--prmtop" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &mdOpt . prmtop, value );
}else if( EQ( key, "--xvv" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . xvv, value );
}else if( EQ( key, "--guv" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . guv, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--cuv" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . cuv, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--huv" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . huv, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--uuv" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . uuv, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--asymp" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . asymp, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--quv" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . quv, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--chgdist" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . chgdist, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--exchem" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . exchem, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--solvene" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . solvene, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--entropy" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . entropy, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--exchemGF" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . exchemGF, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--solveneGF" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . solveneGF, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--entropyGF" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . entropyGF, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--exchemUC" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . exchemUC, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--solveneUC" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . solveneUC, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--entropyUC" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . entropyUC, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--potUV" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
NAB_strcpy(  &rismOpt . potUV, value );
rismOpt . ntwrism = 1;
}else if( EQ( key, "--volfmt" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if(  !( NAB_rematch( value, "dx" ) ) &&  !( NAB_rematch( value, "xyzv" ) ) ){
error( STEMP( __st0001__, "--volfmt must be one of dx or xyzv\n" ) );
usage(  );
}
NAB_strcpy(  &rismOpt . volfmt, value );
}else if( EQ( key, "--closure" ) ){

j = i;
while( j + 1 <= argc &&  !( NAB_rematch( argv[j + 1 - 1], "^-" ) ) ){
j ++ ;
}
nclosure = j - i;

free( ( closure ) );
__gdab0001__ = nclosure;DA_ALLOC( closure = ( STRING_T * * )calloc( __gdab0001__, sizeof( STRING_T * ) ), "main", "closure" );
for( j = 1;j <= nclosure;j ++  ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if(  !( NAB_rematch( value, "[kK][hH]" ) ) &&  !( NAB_rematch( value, "[hH][nN][cC]" ) ) &&  !( NAB_rematch( value, "[pP][sS][eE]" ) ) ){
error( STEMP( __st0001__, "--closure must be one of kh, hnc, or pse\n" ) );
usage(  );
}
NAB_strcpy(  &closure[j - 1], value );
}
}else if( EQ( key, "--closureorder" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%i",  &rismOpt . closureOrder ) ) != 1 || rismOpt . closureOrder < 1 ){
error( STEMP( __st0001__, "--closureOrder takes an integer > 0\n" ) );
usage(  );
}
}else if( EQ( key, "--asympcorr" ) ){
rismOpt . asympcorr = 1;
}else if( EQ( key, "--noasympcorr" ) ){
rismOpt . asympcorr = 0;
}else if( EQ( key, "--buffer" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( ( sscanf( value, "%lf",  &rismOpt . buffer ) ) != 1 ){
error( STEMP( __st0001__, "--buffer takes a float\n" ) );
usage(  );
}
}else if( EQ( key, "--solvcut" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%lf",  &rismOpt . solvcut ) ) != 1 || rismOpt . solvcut <= 0 ){
error( STEMP( __st0001__, "--solvcut takes a float > 0\n" ) );
usage(  );
}
}else if( EQ( key, "--grdspc" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
tempInt = numCSV(  &value );
if( tempInt == 3 ){
if( 


( sscanf( value, "%lf,%lf,%lf",  &rismOpt . grdspcX,  &rismOpt . grdspcY,  &rismOpt . grdspcZ ) ) != 3 || rismOpt . grdspcX <= 0 || rismOpt . grdspcY <= 0 || rismOpt . grdspcZ <= 0 ){
error( STEMP( __st0001__, "--grdspc takes floats > 0\n" ) );
usage(  );
}
}else if( tempInt == 1 ){
if( 
( sscanf( value, "%lf",  &rismOpt . grdspcX ) ) != 1 || rismOpt . grdspcX <= 0 ){
error( STEMP( __st0001__, "--grdspc takes floats\n" ) );
usage(  );
}
rismOpt . grdspcY = rismOpt . grdspcX;
rismOpt . grdspcZ = rismOpt . grdspcX;
}else{
error( STEMP( __st0001__, "--grdspc only takes one or three values\n" ) );
usage(  );
}
}else if( EQ( key, "--ng" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
tempInt = numCSV(  &value );
if( tempInt == 3 ){
if( 

( sscanf( value, "%i,%i,%i",  &rismOpt . ngX,  &rismOpt . ngY,  &rismOpt . ngZ ) ) != 3 || rismOpt . ngX <= 0 || rismOpt . ngY <= 0 || rismOpt . ngZ <= 0 ){
error( STEMP( __st0001__, "--ng takes integers > 0\n" ) );
usage(  );
}
}else if( tempInt == 1 ){
if( 
( sscanf( value, "%i",  &rismOpt . ngX ) ) != 1 || rismOpt . ngX <= 0 ){
error( STEMP( __st0001__, "--ng takes integers\n" ) );
usage(  );
}
rismOpt . ngY = rismOpt . ngX;
rismOpt . ngZ = rismOpt . ngX;
}else{
error( STEMP( __st0001__, "--ng only takes one or three values\n" ) );
usage(  );
}
}else if( EQ( key, "--solvbox" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
tempInt = numCSV(  &value );
if( tempInt == 3 ){
if( 


( sscanf( value, "%lf,%lf,%lf",  &rismOpt . solvboxX,  &rismOpt . solvboxY,  &rismOpt . solvboxZ ) ) != 3 || rismOpt . solvboxX <= 0 || rismOpt . solvboxY <= 0 || rismOpt . solvboxZ <= 0 ){
error( STEMP( __st0001__, "--solvbox takes floats > 0\n" ) );
usage(  );
}
}else if( tempInt == 1 ){
if( 
( sscanf( value, "%lf",  &rismOpt . solvboxX ) ) != 1 || rismOpt . solvboxX <= 0 ){
error( STEMP( __st0001__, "--solvbox takes floats > 0\n" ) );
usage(  );
}
rismOpt . solvboxY = rismOpt . solvboxX;
rismOpt . solvboxZ = rismOpt . solvboxX;
}else{
error( STEMP( __st0001__, "--solvbox only takes one or three values\n" ) );
usage(  );
}
}else if( EQ( key, "--tolerance" ) ){

j = i;
while( j + 1 <= argc &&  !( NAB_rematch( argv[j + 1 - 1], "^-" ) ) ){
NAB_strcpy(  &value, argv[j + 1 - 1] );
j ++ ;
}
ntolerance = j - i;

free( ( tolerance ) );
__gdab0002__ = ntolerance;DA_ALLOC( tolerance = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "tolerance" );
for( j = 1;j <= ntolerance;j ++  ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%lf",  &tolerance[j - 1] ) ) != 1 || tolerance[j - 1] <= 0 ){
error( STEMP( __st0001__, "--tolerance takes a float > 0\n" ) );
usage(  );
}
}
}else if( EQ( key, "--mdiis_del" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%lf",  &rismOpt . mdiis_del ) ) != 1 || rismOpt . mdiis_del <= 0 || rismOpt . mdiis_del > 2 ){
error( STEMP( __st0001__, "--mdiis_del takes a float > 0 and <= 2\n" ) );
usage(  );
}
}else if( EQ( key, "--mdiis_restart" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%lf",  &rismOpt . mdiis_restart ) ) != 1 || rismOpt . mdiis_restart <= 0 ){
error( STEMP( __st0001__, "--mdiis_del takes a float > 0\n" ) );
usage(  );
}
}else if( EQ( key, "--mdiis_nvec" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%i",  &rismOpt . mdiis_nvec ) ) != 1 || rismOpt . mdiis_nvec < 1 ){
error( STEMP( __st0001__, "--mdiis_nvec takes an integer > 0\n" ) );
usage(  );
}
}else if( EQ( key, "--mdiis_method" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%i",  &rismOpt . mdiis_method ) ) != 1 || rismOpt . mdiis_method < 0 || rismOpt . mdiis_method > 2 ){
error( STEMP( __st0001__, "--mdiis_method must be 0, 1 or 2\n" ) );
usage(  );
}
}else if( EQ( key, "--maxstep" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%i",  &rismOpt . maxstep ) ) != 1 || rismOpt . maxstep < 1 ){
error( STEMP( __st0001__, "--maxstep takes an integer > 0\n" ) );
usage(  );
}
}else if( EQ( key, "--npropagate" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%i",  &rismOpt . npropagate ) ) != 1 || rismOpt . npropagate < 0 ){
error( STEMP( __st0001__, "--npropagate takes an integer >= 0\n" ) );
usage(  );
}
}else if( EQ( key, "--centering" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%i",  &rismOpt . centering ) ) != 1 || rismOpt . centering <  - 4 || rismOpt . centering > 4 ){
error( STEMP( __st0001__, "--centering must be between -4 and 4\n" ) );
usage(  );
}








}else if( EQ( key, "--polarDecomp" ) ){
rismOpt . polarDecomp = 1;
}else if( EQ( key, "--nopolarDecomp" ) ){
rismOpt . polarDecomp = 0;
}else if( EQ( key, "--entropicDecomp" ) ){
rismOpt . entropicDecomp = 1;
}else if( EQ( key, "--noentropicDecomp" ) ){
rismOpt . entropicDecomp = 0;
}else if( EQ( key, "--verbose" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
if( 
( sscanf( value, "%i",  &rismOpt . verbose ) ) != 1 || rismOpt . verbose < 0 || rismOpt . verbose > 2 ){
error( STEMP( __st0001__, "--verbose must be 0, 1 or 2\n" ) );
usage(  );
}
}else if( EQ( key, "--progress" ) ){
rismOpt . progress = 1;
}else if( EQ( key, "--noprogress" ) ){
rismOpt . progress = 0;
}else if( EQ( key, "--saveprogress" ) ){
rismOpt . saveprogress = 1;
}else if( EQ( key, "--nosaveprogress" ) ){
rismOpt . saveprogress = 0;
}else if( EQ( key, "--uccoeff" ) ){
i ++ ;
if( i <= argc ){
NAB_strcpy(  &value, argv[i - 1] );
}
testValue(  &key,  &value );
tempInt = numCSV(  &value );
if( tempInt == 2 ){
if( ( sscanf( value, "%lf,%lf",  &rismOpt . uccoeff1,  &rismOpt . uccoeff2 ) ) != 2 ){
error( STEMP( __st0001__, "--uccoeff takes floats\n" ) );
usage(  );
}
}else if( tempInt == 1 ){
if( ( sscanf( value, "%lf",  &rismOpt . uccoeff1 ) ) != 1 ){
error( STEMP( __st0001__, "--uccoeff takes floats\n" ) );
usage(  );
}
}else{
error( STEMP( __st0001__, "--uccoeff only takes one or two values\n" ) );
usage(  );
}
}else{
error( STEMP( __st0002__, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "unknown option: '%s'\n", key ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) ) );
usage(  );
}
i ++ ;
}

checkOptions(  &mdOpt,  &rismOpt );




m = getpdb( mdOpt . pdb, NULL );
readparm( m, mdOpt . prmtop );
__gdab0001__ = 3 *  *( NAB_mri( m, "natoms" ) );DA_ALLOC( p_xyz = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "p_xyz" );
__gdab0002__ = 3 *  *( NAB_mri( m, "natoms" ) );DA_ALLOC( f_xyz = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "f_xyz" );
__gdab0003__ = 3 *  *( NAB_mri( m, "natoms" ) );DA_ALLOC( v_xyz = ( REAL_T * )malloc( __gdab0003__ * ( sizeof( REAL_T ) ) ), "main", "v_xyz" );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "cut=%e", sqrt( 1.000000E+38 ) ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "rism=1, ntpr_rism=1, xvvfile=%s", rismOpt . xvv ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
if( rismOpt . guv ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "guvfile=%s", rismOpt . guv ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . cuv ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "cuvfile=%s", rismOpt . cuv ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . huv ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "huvfile=%s", rismOpt . huv ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . uuv ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "uuvfile=%s", rismOpt . uuv ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . asymp ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "asympfile=%s", rismOpt . asymp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . quv ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "quvfile=%s", rismOpt . quv ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . chgdist ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "chgdistfile=%s", rismOpt . chgdist ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . exchem ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "exchemfile=%s", rismOpt . exchem ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . solvene ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "solvenefile=%s", rismOpt . solvene ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . entropy ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "entropyfile=%s", rismOpt . entropy ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . exchemGF ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "exchemGFfile=%s", rismOpt . exchemGF ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . solveneGF ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "solveneGFfile=%s", rismOpt . solveneGF ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . entropyGF ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "entropyGFfile=%s", rismOpt . entropyGF ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . exchemUC ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "exchemUCfile=%s", rismOpt . exchemUC ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . solveneUC ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "solveneUCfile=%s", rismOpt . solveneUC ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . entropyUC ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "entropyUCfile=%s", rismOpt . entropyUC ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
if( rismOpt . potUV ){
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "potUVfile=%s", rismOpt . potUV ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "volfmt=%s", rismOpt . volfmt ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &tmpstr, "" );
if( nclosure >= 1 ){
NAB_strcpy(  &tmpstr, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "closure=%s", closure[1 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
for( i = 2;i <= nclosure;i ++  ){
NAB_strcpy(  &tmpstr, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%s", tmpstr, closure[i - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s, closureOrder=%i", tmpstr, rismOpt . closureOrder ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "buffer=%g, solvcut=%g", rismOpt . buffer, rismOpt . solvcut ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "grdspcx=%g, grdspcy=%g, grdspcz=%g", rismOpt . grdspcX, rismOpt . grdspcY, rismOpt . grdspcZ ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "solvboxx=%g, solvboxy=%g, solvboxz=%g", rismOpt . solvboxX, rismOpt . solvboxY, rismOpt . solvboxZ ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ngx=%d, ngy=%d, ngz=%d", rismOpt . ngX, rismOpt . ngY, rismOpt . ngZ ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
NAB_strcpy(  &tmpstr, "" );
if( ntolerance >= 1 ){
NAB_strcpy(  &tmpstr, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "tolerance=%g", tolerance[1 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
for( i = 2;i <= ntolerance;i ++  ){
NAB_strcpy(  &tmpstr, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s,%g", tmpstr, tolerance[i - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s, mdiis_del=%g", tmpstr, rismOpt . mdiis_del ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "mdiis_restart=%g", rismOpt . mdiis_restart ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "mdiis_nvec=%d, mdiis_method=%d", rismOpt . mdiis_nvec, rismOpt . mdiis_method ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "maxstep=%d, npropagate=%d", rismOpt . maxstep, rismOpt . npropagate ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "centering=%d, zerofrc=%d, apply_rism_force=%d", rismOpt . centering, rismOpt . zerofrc, rismOpt . apply_rism_force ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "polarDecomp=%d", rismOpt . polarDecomp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "entropicDecomp=%d", rismOpt . entropicDecomp ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "ntwrism=%d, verbose=%d, progress=%d", rismOpt . ntwrism, rismOpt . verbose, rismOpt . progress ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "asympCorr=%d", rismOpt . asympcorr ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "saveprogress=%d", rismOpt . saveprogress ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "uccoeff=%g,%g", rismOpt . uccoeff1, rismOpt . uccoeff2 ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mme_init( m, NULL, "::Z", p_xyz, NULL );

if( mdOpt . traj ){

if( ( netcdfLoad(  &nc, mdOpt . traj ) ) == 0 ){
netcdfLoad(  &nc, mdOpt . traj );
if( mytaskid == 0 )printf( "\nProcessing NetCDF trajectory: %s\n", mdOpt . traj );
while( netcdfGetNextFrame(  &nc, p_xyz, NULL, NULL ) ){
if( mytaskid == 0 )printf( "\nFrame: %d of %d\n", nc . currentFrame, nc . ncframe );
mme( p_xyz, f_xyz,  &nc . currentFrame );
}
netcdfClose(  &nc );
}else{

if( mytaskid == 0 ){
if( numtasks > 1 ){
error( STEMP( __st0001__, "ASCII trajectories not supported for more than one process" ) );
exit( 1 );
}
printf( "\nProcessing ASCII trajectory: %s\n", mdOpt . traj );
}
trajDone = 0;
asciiTraj = fopen( mdOpt . traj, "r" );
if( asciiTraj == NULL ){
error( STEMP( __st0002__, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "Failed to open '%s'", mdOpt . traj ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) ) );
exit( 1 );
}
NAB_getline( asciiTraj );
for( iframe = 1;;iframe ++  ){
for( ipos = 1;ipos <= 3 *  *( NAB_mri( m, "natoms" ) );ipos ++  ){
if( ( fscanf( asciiTraj, "%lf",  &p_xyz[ipos - 1] ) ) < 1 ){
trajDone = 1;
break;
}
}
if( trajDone )break;
if( mytaskid == 0 ){
printf( "\nFrame: %d\n", iframe );
}
mme( p_xyz, f_xyz,  &iframe );
}
}

}else{
if( mdOpt . rst ){
if( mytaskid == 0 )printf( "\nProcessing restart file: %s\n", mdOpt . rst );
getxv(  &mdOpt . rst, ITEMP( __it0001__,  *( NAB_mri( m, "natoms" ) ) ),  &time, p_xyz, v_xyz );
}else{
if( mytaskid == 0 )printf( "\nProcessing PDB file: %s\n", mdOpt . pdb );
setxyz_from_mol(  &m, NULL, p_xyz );
}
mme( p_xyz, f_xyz, ITEMP( __it0001__, 1 ) );

}

if( mytaskid == 0 ){printf( "\n3D-RISM processing complete.\n" );}

mme_rism_max_memory(  );
mme_timer(  );


	exit( 0 );
}
