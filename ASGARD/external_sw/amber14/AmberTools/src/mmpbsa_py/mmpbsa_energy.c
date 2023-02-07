#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;





































 struct AmberNetcdf{REAL_T temp0;REAL_T restartTime;INT_T isNCrestart;INT_T ncid;INT_T frameDID;INT_T ncframe;INT_T currentFrame;INT_T atomDID;INT_T ncatom;INT_T ncatom3;INT_T coordVID;INT_T velocityVID;INT_T cellAngleVID;INT_T cellLengthVID;INT_T spatialDID;INT_T labelDID;INT_T cell_spatialDID;INT_T cell_angularDID;INT_T spatialVID;INT_T timeVID;INT_T cell_spatialVID;INT_T cell_angularVID;INT_T TempVID;};







 struct CLOptions{STRING_T *mdin,  *prmtop,  *inpcrd,  *traj,  *mdout;};





 struct GBOptions{INT_T igb, gbsa, rungb;REAL_T extdiel, saltcon, surften, rgbmax;};







 struct PBOptions{INT_T runpb, inp, smoothopt, radiopt, npbopt, solvopt, maxitn;INT_T nfocus, bcopt, eneopt, fscale, dbfopt;REAL_T epsin, epsout, istrng, dprob, iprob, accept, fillratio;REAL_T space, cutnb, sprob, cavity_surften, cavity_offset;};



static  struct GBOptions gbopt;

static  struct PBOptions pbopt;

static  struct CLOptions clopt;



INT_T ReadGBin( STRING_T * *filename );
INT_T ReadGBin( STRING_T * *filename ){
STRING_T *option = NULL,  *line = NULL,  *value = NULL;

FILE_T *gbfile;


gbfile = fopen(  *filename, "r" );
if( gbfile == NULL ){
fprintf( stderr, "Error: Cannot open %s for reading!\n",  *filename );
return(  - 1 );
}


NAB_strcpy(  &line, NAB_getline( gbfile ) );
if( NE( line, "GB" ) )
return( 1 );

gbopt . rungb = 1;

while( NAB_strcpy(  &line, NAB_getline( gbfile ) ) ){
sscanf( line, "%s = %s", NAB_readstring(  &option ), NAB_readstring(  &value ) );
if( EQ( option, "igb" ) )gbopt . igb = atoi( value );
else if( EQ( option, "gbsa" ) )gbopt . gbsa = atoi( value );
else if( EQ( option, "extdiel" ) )gbopt . extdiel = atof( value );
else if( EQ( option, "saltcon" ) )gbopt . saltcon = atof( value );
else if( EQ( option, "surften" ) )gbopt . surften = atof( value );
else if( EQ( option, "rgbmax" ) )gbopt . rgbmax = atof( value );
else{
fprintf( stderr, "Error: Unknown GB option %s!", option );
fclose( gbfile );
return(  - 1 );
}
}

return( fclose( gbfile ) );

}



INT_T ReadPBin( STRING_T * *filename );
INT_T ReadPBin( STRING_T * *filename ){
STRING_T *option = NULL,  *line = NULL,  *value = NULL;

FILE_T *pbfile;


pbfile = fopen(  *filename, "r" );
if( pbfile == NULL ){
fprintf( stderr, "Error: Cannot open %s for reading!\n",  *filename );
return(  - 1 );
}


NAB_strcpy(  &line, NAB_getline( pbfile ) );
if( NE( line, "PB" ) ){

fprintf( stderr, "Error: Cannot determine type of input file (%s)!",  *filename );
return( 1 );
}

pbopt . runpb = 1;

while( NAB_strcpy(  &line, NAB_getline( pbfile ) ) ){
sscanf( line, "%s = %s", NAB_readstring(  &option ), NAB_readstring(  &value ) );
if( EQ( option, "inp" ) )pbopt . inp = atoi( value );
else if( EQ( option, "smoothopt" ) )pbopt . smoothopt = atoi( value );
else if( EQ( option, "radiopt" ) )pbopt . radiopt = atoi( value );
else if( EQ( option, "npbopt" ) )pbopt . npbopt = atoi( value );
else if( EQ( option, "solvopt" ) )pbopt . solvopt = atoi( value );
else if( EQ( option, "maxitn" ) )pbopt . maxitn = atoi( value );
else if( EQ( option, "nfocus" ) )pbopt . nfocus = atoi( value );
else if( EQ( option, "bcopt" ) )pbopt . bcopt = atoi( value );
else if( EQ( option, "eneopt" ) )pbopt . eneopt = atoi( value );
else if( EQ( option, "fscale" ) )pbopt . fscale = atoi( value );
else if( EQ( option, "dbfopt" ) )pbopt . dbfopt = atoi( value );
else if( EQ( option, "epsin" ) )pbopt . epsin = atof( value );
else if( EQ( option, "epsout" ) )pbopt . epsout = atof( value );
else if( EQ( option, "istrng" ) )pbopt . istrng = atof( value );
else if( EQ( option, "dprob" ) )pbopt . dprob = atof( value );
else if( EQ( option, "iprob" ) )pbopt . iprob = atof( value );
else if( EQ( option, "accept" ) )pbopt . accept = atof( value );
else if( EQ( option, "fillratio" ) )pbopt . fillratio = atof( value );
else if( EQ( option, "space" ) )pbopt . space = atof( value );
else if( EQ( option, "cutnb" ) )pbopt . cutnb = atof( value );
else if( EQ( option, "sprob" ) )pbopt . sprob = atof( value );
else if( EQ( option, "cavity_surften" ) )pbopt . cavity_surften = atof( value );
else if( EQ( option, "cavity_offset" ) )pbopt . cavity_offset = atof( value );
else{
fprintf( stderr, "Error: Unknown PB option %s!", option );
fclose( pbfile );
return(  - 1 );
}

}

return( ( fclose( pbfile ) ) );

}


INT_T SetDefaults(  );
INT_T SetDefaults(  ){


NAB_strcpy(  &clopt . mdin, "mdin" );
NAB_strcpy(  &clopt . prmtop, "prmtop" );
NAB_strcpy(  &clopt . inpcrd, "pdb" );
NAB_strcpy(  &clopt . traj, "mdcrd" );
NAB_strcpy(  &clopt . mdout, "mdout" );


gbopt . rungb = 0;
gbopt . gbsa = 0;
gbopt . igb = 1;
gbopt . extdiel = 7.850000E+01;
gbopt . surften = 7.200000E-03;
gbopt . saltcon = 0.000000E+00;
gbopt . rgbmax = 9.990000E+02;


pbopt . runpb = 0;
pbopt . inp = 2;
pbopt . smoothopt = 1;
pbopt . radiopt = 1;
pbopt . npbopt = 0;
pbopt . solvopt = 1;
pbopt . maxitn = 100;
pbopt . nfocus = 2;
pbopt . fscale = 8;
pbopt . dbfopt = 1;
pbopt . epsin = 1.000000E+00;
pbopt . epsout = 8.000000E+01;
pbopt . istrng = 0.000000E+00;
pbopt . dprob = 1.400000E+00;
pbopt . iprob = 2.000000E+00;
pbopt . accept = 1.000000E-03;
pbopt . fillratio = 2.000000E+00;
pbopt . space = 5.000000E-01;
pbopt . bcopt = 5;
pbopt . eneopt = 2;
pbopt . cutnb = 0.000000E+00;
pbopt . sprob = 5.570000E-01;
pbopt . cavity_surften = 3.780000E-02;
pbopt . cavity_offset =  - 5.692000E-01;

return( 0 );
}


INT_T printusage( STRING_T * *prog_name );
INT_T printusage( STRING_T * *prog_name ){
printf( " Usage: %s [-O] -i mdin -o mdout -p prmtop -c pdb -y mdcrd\n",  *prog_name );
return( 0 );
}



INT_T CheckValues(  );
INT_T CheckValues(  ){
INT_T isinerr;

isinerr = 0;
if( gbopt . rungb ){
if( gbopt . igb != 1 && gbopt . igb != 2 && gbopt . igb != 5 && gbopt . igb != 7 && gbopt . igb != 8 ){
fprintf( stderr, "Error: IGB must be 1, 2, 5, 7, or 8!\n" );
isinerr = 1;
}if( gbopt . extdiel <= 0 ){
fprintf( stderr, "Error: EXTDIEL must be positive!\n" );
isinerr = 1;
}if( gbopt . surften < 0 ){
fprintf( stderr, "Error: SURFTEN must be non-negative!\n" );
isinerr = 1;
}if( gbopt . saltcon < 0 ){
fprintf( stderr, "Error: SALTCON must be non-negative!\n" );
isinerr = 1;
}if( gbopt . rgbmax <= 0 ){
fprintf( stderr, "Error: RGBMAX must be positive!\n" );
isinerr = 1;
}if( gbopt . rgbmax < 20 && gbopt . rgbmax > 0 ){
fprintf( stderr, "Warning: Low value for RGBMAX. Consider using default.\n" );
}
}else{
if( pbopt . inp != 0 && pbopt . inp != 1 && pbopt . inp != 2 ){
fprintf( stderr, "Error: INP must be 0, 1, or 2!\n" );
isinerr = 1;
}if( pbopt . inp == 1 ){
fprintf( stderr, "Warning: inp=1 was old default\n" );
if( ( floor( 1.000000E+04 * pbopt . cavity_surften ) ) == ( floor( 1.000000E+04 * 3.780000E-02 ) ) ){
fprintf( stderr, "Warning: cavity_surften=0.0378 not recommended for inp=1, switching to inp=1 default value: 0.0050\n" );
pbopt . cavity_surften = 5.000000E-03;
}
if( ( floor( 1.000000E+04 * pbopt . cavity_offset ) ) == ( floor(  - 5.692000E-01 * 1.000000E+04 ) ) ){
fprintf( stderr, "Warning: cavity_offset=-0.5692 not recommended for inp=1, switching to inp=1 default value: 0.000\n" );
pbopt . cavity_offset = 0.000000E+00;
}
if( ( floor( pbopt . sprob * 1.000000E+04 ) ) == ( floor( 5.570000E-01 * 1.000000E+04 ) ) ){
fprintf( stderr, "Warning: sprob=.557 not recommended for inp=1, switching to inp=1 default value: 1.400\n" );
pbopt . sprob = 1.400000E+00;
}
if( pbopt . radiopt != 0 && pbopt . inp == 1 ){
fprintf( stderr, "Warning: radiopt should be set to 0 for inp=1\n" );
pbopt . radiopt = 0;
}
}if( pbopt . smoothopt != 0 && pbopt . smoothopt != 1 && pbopt . smoothopt != 2 ){
fprintf( stderr, "Error: SMOOTHOPT must be 0, 1, or 2!\n" );
isinerr = 1;
}if( pbopt . radiopt != 1 && pbopt . radiopt != 0 ){
fprintf( stderr, "Error: RADIOPT must be 0, 1, or 2!\n" );
isinerr = 1;
}if( pbopt . npbopt != 0 && pbopt . npbopt != 1 ){
fprintf( stderr, "Error: NPBOPT must be 0 or 1!\n" );
isinerr = 1;
}if( pbopt . solvopt < 1 || pbopt . solvopt == 7 || pbopt . solvopt > 8 ){
fprintf( stderr, "Error: SOLVOPT must be 1, 2, 3, 4, 5, 6, or 8!\n" );
isinerr = 1;
}if( pbopt . maxitn < 1 ){
fprintf( stderr, "Error: MAXITN must be a positive integer!\n" );
isinerr = 1;
}if( pbopt . nfocus != 1 && pbopt . nfocus != 2 ){
fprintf( stderr, "Error: NFOCUS must be 1 or 2!\n" );
isinerr = 1;
}if( pbopt . fscale <= 0 ){
fprintf( stderr, "Error: NFOCUS must be non-negative!\n" );
isinerr = 1;
}if( pbopt . epsin < 0 || pbopt . epsout < 0 ){
fprintf( stderr, "Error: EPSIN/OUT must be non-negative!\n" );
isinerr = 1;
}if( pbopt . istrng < 0 ){
fprintf( stderr, "Error: ISTRNG must be non-negative!\n" );
isinerr = 1;
}if( pbopt . dprob < 0 ){
fprintf( stderr, "Error: DPROB must be non-negative!\n" );
isinerr = 1;
}if( pbopt . iprob < 0 ){
fprintf( stderr, "Error: IPROB must be non-negative!\n" );
isinerr = 1;
}if( pbopt . accept <= 0 ){
fprintf( stderr, "Error: ACCEPT must be positive!\n" );
isinerr = 1;
}if( pbopt . fillratio <= 0 ){
fprintf( stderr, "Error: FILLRATIO must be positive!\n" );
isinerr = 1;
}if( pbopt . space <= 0 ){
fprintf( stderr, "Error: SPACE must be positive!\n" );
isinerr = 1;
}if( pbopt . bcopt != 1 && pbopt . bcopt != 5 && pbopt . bcopt != 6 && pbopt . bcopt != 10 ){
fprintf( stderr, "Error: BCOPT must be 1, 5, 6, or 8!\n" );
isinerr = 1;
}if( pbopt . eneopt != 1 && pbopt . eneopt != 2 ){
fprintf( stderr, "Error: ENEOPT must be 1 or 2!\n" );
isinerr = 1;
}if( pbopt . cutnb < 0 ){
fprintf( stderr, "Error: CUTNB must be non-negative!\n" );
isinerr = 1;
}if( pbopt . sprob < 0 ){
fprintf( stderr, "Error: SPROB must be non-negative!\n" );
isinerr = 1;
}if( pbopt . cavity_surften < 0 ){
fprintf( stderr, "Error: CAVITY_SURFTEN must be non-negative!\n" );
isinerr = 1;
}
}
return( ( isinerr ) );
}







static INT_T i, j, fr_result, frame, trajDone;

static REAL_T kappa, cut;

static REAL_T *coords,  *grad;

static  struct AmberNetcdf nctraj;

static FILE_T *output,  *asciitraj;

static STRING_T *line = NULL,  *prog_name = NULL;

static MOLECULE_T *mol;






















int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static INT_T __gdab0001__;
static INT_T __gdab0002__;
static STRING_T *__st0001__ = NULL;
SetDefaults(  );

NAB_strcpy(  &prog_name, argv[1 - 1] );
i = 2;
while( i <= argc ){
if( EQ( argv[i - 1], "-i" ) )
NAB_strcpy(  &clopt . mdin, argv[ ++ i - 1] );
else if( EQ( argv[i - 1], "-p" ) )
NAB_strcpy(  &clopt . prmtop, argv[ ++ i - 1] );
else if( EQ( argv[i - 1], "-c" ) )
NAB_strcpy(  &clopt . inpcrd, argv[ ++ i - 1] );
else if( EQ( argv[i - 1], "-o" ) )
NAB_strcpy(  &clopt . mdout, argv[ ++ i - 1] );
else if( EQ( argv[i - 1], "-y" ) )
NAB_strcpy(  &clopt . traj, argv[ ++ i - 1] );
else if( EQ( argv[i - 1], "-r" ) )
i ++ ;
else if( ( EQ( argv[i - 1], "-h" ) ) || ( EQ( argv[i - 1], "--help" ) ) || ( EQ( argv[i - 1], "--h" ) ) ){
printusage(  &prog_name );
exit( 0 );
}else if( NE( argv[i - 1], "-O" ) ){
fprintf( stderr, "Error: Bad flag %s!\n", argv[i - 1] );
printusage(  &prog_name );
exit( 1 );
}

i ++ ;
}






fr_result = ReadGBin(  &clopt . mdin );

if( fr_result ==  - 1 )
exit(  - 1 );
else if( fr_result == 1 ){
if( ( ReadPBin(  &clopt . mdin ) ) != 0 )
exit(  - 1 );
}


if( ( CheckValues(  ) ) == 1 )exit( 1 );


output = fopen( clopt . mdout, "w" );
nabout = output;


mol = getpdb( clopt . inpcrd, NULL );
readparm( mol, clopt . prmtop );
__gdab0001__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( coords = ( REAL_T * )malloc( __gdab0001__ * ( sizeof( REAL_T ) ) ), "main", "coords" );
__gdab0002__ = 3 *  *( NAB_mri( mol, "natoms" ) );DA_ALLOC( grad = ( REAL_T * )malloc( __gdab0002__ * ( sizeof( REAL_T ) ) ), "main", "grad" );


if( gbopt . rungb == 1 ){
kappa = sqrt( 1.080600E-01 * gbopt . saltcon );

cut = 9.990000E+02;
if( gbopt . rgbmax > cut )cut = gbopt . rgbmax;


mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "e_debug=2, gb=%d, rgbmax=%lf, gbsa=%d, surften=%lf, cut=%lf", gbopt . igb, gbopt . rgbmax, gbopt . gbsa, gbopt . surften, cut ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "epsext=%lf, kappa=%lf", gbopt . extdiel, kappa ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}else{

mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "e_debug=3, ipb=2, inp=%d, epsin=%lf, epsout=%lf, smoothopt=%d", pbopt . inp, pbopt . epsin, pbopt . epsout, pbopt . smoothopt ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "istrng=%lf, radiopt=%d, dprob=%lf, iprob=%lf, npbopt=%d", pbopt . istrng, pbopt . radiopt, pbopt . dprob, pbopt . iprob, pbopt . npbopt ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "solvopt=%d, accept=%lf, maxitn=%d, fillratio=%lf, space=%lf", pbopt . solvopt, pbopt . accept, pbopt . maxitn, pbopt . fillratio, pbopt . space ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );

mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "nfocus=%d, fscale=%d, bcopt=%d, eneopt=%d, cutnb=%lf, sprob=%lf", pbopt . nfocus, pbopt . fscale, pbopt . bcopt, pbopt . eneopt, pbopt . cutnb, pbopt . sprob ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
mm_options( SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "cavity_surften=%lf, cavity_offset=%lf", pbopt . cavity_surften, pbopt . cavity_offset ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
}

mme_init( mol, NULL, "::Z", coords, NULL );


frame = 0;
if( ( netcdfLoad(  &nctraj, clopt . traj ) ) == 0 ){
netcdfLoad(  &nctraj, clopt . traj );
fprintf( nabout, "Processing NetCDF trajectory (%s)\n\n", clopt . traj );
while( netcdfGetNextFrame(  &nctraj, coords, NULL, NULL ) ){
fprintf( nabout, "Processing frame %d\n", nctraj . currentFrame );
mme( coords, grad,  &frame );
fprintf( nabout, "\n" );
frame ++ ;
}
netcdfClose(  &nctraj );
}else{
fprintf( nabout, "Processing ASCII trajectory (%s)\n\n", clopt . traj );
asciitraj = fopen( clopt . traj, "r" );

if( asciitraj == NULL ){
fprintf( stderr, "Error: Could not open trajectory file (%s) for reading!", clopt . traj );
exit(  - 1 );
}

trajDone = 0;
NAB_getline( asciitraj );

for( i = 1;;i ++  ){
for( j = 1;j <= 3 *  *( NAB_mri( mol, "natoms" ) );j ++  ){
if( ( fscanf( asciitraj, "%lf",  &coords[j - 1] ) ) < 1 ){
trajDone = 1;
break;
}
}
if( trajDone )break;

fprintf( nabout, "Processing frame %d\n", i );
mme( coords, grad,  &frame );
fprintf( nabout, "\n" );
frame ++ ;
}
}

if( gbopt . rungb )
fprintf( nabout, "MM/GBSA processing done!\n" );
else
fprintf( nabout, "MM/PBSA processing done!\n" );


exit( fclose( output ) );



	exit( 0 );
}
