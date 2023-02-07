#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

/* change USE_MKSTEMP in this define to dont_USE_MKSTEMP if  */
/* function tmpnam is desired for temporary file creation. */
#define USE_MKSTEMP absolutely
#ifndef USE_MKSTEMP
#include <fcntl.h>  /* O_CREAT, O_ etc. */
#endif

typedef	unsigned long	NAB_SIZE_T;

static	void	n2c( int, int, char *, int, int, int, char *[], char [], char * );
static	void	cc( int, int, int, int, char [], int, char *[], char [] );
static	int	split( char [], char *[], char * );
static	int	execute( char **, int, int, int );

int main( int argc, char *argv[] )
{
	int	ac;
	int	aopt = 0;
	int	copt = 0;
	int	cgdopt = 0;
	char	*cgdval = NULL;
	int	noassert = 0;
	int	nodebug = 0;
	int	oopt = 0;
	char	cppstring[ 1024 ];
	char	*nfmask;
	char	*dotp;
	char	ofname[ 256 ];
	
	cppstring[ 0 ] = '\0';
	if( argc == 1 ){
		fprintf( stderr,
"usage: %s [-avs] [-c] [-Dstring] [-noassert] [-nodebug] [-o file] [-v] file(s)\n",
			argv[ 0 ] ); 
		exit( 1 );
	}

	if( ( nfmask = ( char * )malloc( (argc + 1)*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "%s: can't allocate arg mask\n", argv[ 0 ] );
		exit( 1 );
	}
	
	/* get nab options:	*/
	strcpy( ofname, "a.out" );
	for( ac = 1; ac < argc; ac++ ){
		nfmask[ ac ] = 0;
		if( strcmp( argv[ ac ], "-avs" ) == 0 ){
			aopt = 1;
			argv[ ac ] = NULL;
		}else if( strcmp( argv[ ac ], "-c" ) == 0 ){
			copt = 1;
		}else if( strncmp( argv[ ac ], "-cgdebug", 8 ) == 0 ){
			cgdopt = 1;
			cgdval = argv[ ac ]; 
			argv[ ac ] = NULL;
		}else if( strncmp( argv[ ac ], "-D", 2 ) == 0 ){
			strcat( cppstring, argv[ ac ] );
			strcat( cppstring, " ");
		}else if( strncmp( argv[ ac ], "-I", 2 ) == 0 ){
			strcat( cppstring, argv[ ac ] );
			strcat( cppstring, " ");
		}else if( strncmp( argv[ ac ], "-noassert", 9 ) == 0 ){
			noassert = 1;
			argv[ ac ] = NULL;
		}else if( strncmp( argv[ ac ], "-nodebug", 8 ) == 0 ){
			nodebug = 1;
			argv[ ac ] = NULL;
		}else if( strcmp( argv[ ac ], "-o" ) == 0 ){
			oopt = 1;
			ac++;
			if( ac == argc ){
				fprintf( stderr, "%s: -o requires file name\n",
					argv[ 0 ] );
				exit( 1 );
			}else{
				strcpy( ofname, argv[ ac ] );
				nfmask[ ac ] = 0;
			}
		}else if( strcmp( argv[ ac ], "-v" ) == 0 ){
			cgdopt = 1;
			cgdval = "";
			argv[ ac ] = NULL;
		}else if( *argv[ ac ] != '-' ){
			if( (dotp = strrchr( argv[ ac ], '.' )) ){
				if( strcmp( dotp, ".nab" ) == 0 )
					nfmask[ ac ] = 1;
			}
		}
	}

	n2c( aopt, cgdopt, cgdval, noassert, nodebug,
		argc, argv, nfmask, cppstring );

	cc( aopt, copt, cgdopt, oopt, ofname, argc, argv, nfmask );

	exit( 0 );
}

static	void	n2c( int aopt, int cgdopt, char *cgdval, int noassert, int nodebug,
	int argc, char *argv[], char nfmask[], char *cppstring )
/*
int	aopt;
int	cgdopt;
char	*cgdval;
int	noassert;
int	nodebug;
int	argc;
char	*argv[];
char	nfmask[];
char	*cppstring;
*/
{
	int	ac;
	char	n2c_ofname[ 256 ];
	char	cpp_ofname[ 256 ];

	char	cmd[ 1024 ];
	char	*fields[ 256 ];
	int	nfields;
	int	cpp_ofd, n2c_ofd;
	int	status;
    char *amberhome;

#ifndef USE_MKSTEMP
	tmpnam( n2c_ofname );
	tmpnam( cpp_ofname );
	if( (n2c_ofd = open( n2c_ofname, O_WRONLY|O_CREAT, 0644 )) < 0 ){
		perror( n2c_ofname ); exit(1);
	}
#else
	strcpy( n2c_ofname, "/tmp/n2c_ofname_XXXXXX" );
	if( ( n2c_ofd = mkstemp( n2c_ofname ) ) < 0 ){
		perror( n2c_ofname );
		exit(1);
	}
#endif  /* USE_MKSTEMP */

	for( ac = 1; ac < argc; ac++ ){
		if( nfmask[ ac ] ){
			if( access( argv[ ac ], F_OK ) ){
				fprintf( stderr,
					"%s: %s: no such file (arg # %d)\n",
					argv[ 0 ], argv[ ac ], ac );
				unlink( n2c_ofname );
				exit( 1 );
			}

           /*  get temp file for cpp output, and run CPP :  */

            amberhome = (char * ) getenv("AMBERHOME");
            if( amberhome == NULL ){
               fprintf( stderr, "AMBERHOME is not set!\n" );
               exit(1);
            }
			sprintf( cmd, "%s/bin/%s %s -I%s/include %s ",
				amberhome, CPP, cppstring, amberhome,
				argv[ ac ] ? argv[ ac ] : "" );
			if( cgdopt ) fprintf( stderr, "cpp cmd: %s\n", cmd );
			nfields = split( cmd, fields, " " );

#ifndef USE_MKSTEMP
			if( (cpp_ofd =
				open( cpp_ofname, O_RDWR|O_CREAT, 0644 )) < 0 ){
				perror( cpp_ofname ); exit(1);
			}
#else
			strcpy( cpp_ofname, "/tmp/cpp_ofname_XXXXXX" );
			if( (cpp_ofd = mkstemp( cpp_ofname ) ) < 0 ){
				perror( cpp_ofname );
				exit(1);
			}
#endif /* USE_MKSTEMP */

			status = execute( fields, 0, cpp_ofd, 2 );
			if( status != 0 ){
				unlink( cpp_ofname );
				fprintf( stderr, "cpp failed!\n" ); exit(1);
			}
			lseek( cpp_ofd, 0L, L_SET );

            /*  next run nab2c:  */

#ifdef MPI
			sprintf( cmd, "%s/bin/mpinab2c %s %s %s %s -nfname %s",
#else
			sprintf( cmd, "%s/bin/nab2c %s %s %s %s -nfname %s",
#endif
				amberhome,
				aopt ? "-avs" : "",
				cgdopt ? cgdval : "",
				noassert ? "-noassert" : "",
				nodebug ? "-nodebug" : "",
				argv[ ac ] ? argv[ ac ] : "" );
			if( cgdopt ) fprintf( stderr, "nab2c cmd: %s\n", cmd );
			nfields = split( cmd, fields, " " );

			status = execute( fields, cpp_ofd, n2c_ofd, 2 );
			unlink( cpp_ofname );
			if( status != 0 ){
				unlink( n2c_ofname );
				fprintf( stderr, "nab2c failed!\n" ); exit(1);
			}
		}
	}
	unlink( n2c_ofname );
}

static	void	cc( int aopt, int copt, int cgdopt, int oopt, char ofname[],
	int argc, char *argv[], char nfmask[] )
/*
int	aopt;
int	copt;
int	cgdopt;
int	oopt;
char	ofname[];
int	argc;
char	*argv[];
char	nfmask[];
*/
{
	int	ac;
	char	*dotp, *amberhome;
	char	cmd[ 1024 ], word[ 1024 ];
	char	*fields[256];
	int 	status;
	int 	nfields;

    amberhome = (char *) getenv("AMBERHOME");
    if( amberhome == NULL ){
       fprintf( stderr, "AMBERHOME is not set!\n" );
       exit(1);
    }
	sprintf( cmd, "%s -I%s/include", CC, amberhome );
	if( aopt ){
#ifdef AVSDIR
		sprintf( word, " -I%s/include", AVSDIR );
		strcat( cmd, word );
#else
		fprintf( stderr, "-avs not supported.\n" );
#endif
	}
	for( ac = 1; ac < argc; ac++ ){
		if( nfmask[ ac ] ){
			dotp = strrchr( argv[ ac ], '.' );
			strcpy( &dotp[ 1 ], "c" );
			sprintf( word, " %s", argv[ ac ] );
			strcat( cmd, word );
		}else if( argv[ ac ] ){
			sprintf( word, " %s", argv[ ac ] );
			strcat( cmd, word );
		}
	}
	if( !copt ){
#ifdef AVSDIR
		if( aopt ){
			sprintf( word, " -L%s/lib -lflow_c -lgeom -lutil",
				AVSDIR );
			strcat( cmd, word );
		}
#endif
		sprintf( word, " -L%s/lib -lnab -lcifparse", amberhome );
		strcat( cmd, word );
		sprintf( word, " %s ", FLIBS );
		strcat( cmd, word );
		strcat( cmd, " -lm" );
	}
	if( cgdopt ) fprintf( stderr, "cc cmd: %s\n", cmd );
	nfields = split(cmd,fields," ");
	status = execute(fields, 0, 1, 2); 
	if (status != 0) {
		fprintf(stderr,"cc failed!\n"); exit (1);
	}
}

static	int	split( char str[], char *fields[], char *fsep )
{
	int	nf, flen;
	char	*sp, *fp, *efp, *nfp;

	if( !str )
		return( 0 );

	/* Use fsep of white space is special */ 
	if( strspn( fsep, " \t\n" ) ){
		for( nf = 0, sp = str; ; ){
			fp = sp + strspn( sp, fsep );
			if( !*fp )
				return( nf );
			if( (efp = strpbrk( fp, fsep )) ){
				if( !( flen = efp - fp ) )
					return( nf );
				nfp = (char *)malloc( (flen + 1) * sizeof(char) );
				strncpy( nfp, fp, flen );
				nfp[ flen ] = '\0';
				fields[ nf ] = nfp;
				sp = efp;
				fields[++nf] = NULL;
			}else{
				flen = strlen( fp );
				nfp = (char *)malloc( (flen + 1) * sizeof(char) );
				strcpy( nfp, fp );
				fields[ nf ] = nfp;
				fields[++nf]=NULL;
				return( nf );
			}
		}
	}else{
		for( nf = 0, sp = str; ; ){
			if( (fp = strchr( sp, *fsep )) ){
				flen = fp - sp;
				nfp = (char *)malloc( (flen + 1) * sizeof(char) );
				strncpy( nfp, sp, flen );
				nfp[ flen ] = '\0';
				fields[ nf ] = nfp;
				fields[++nf]=NULL;
				sp = fp + 1;
			}else{
				flen = strlen( sp );
				nfp = (char *)malloc( (flen + 1) * sizeof(char) );
				strcpy( nfp, sp );
				fields[ nf ] = nfp;
				fields[++nf] = NULL;
				return( nf );
			}
		}
	}
}

static	int	execute( char **args, int sin, int sout, int serr )
/*                   from David A Curry, 
                     "Using C on the UNIX System" pp. 105-106    */
/*
char	**args;
int	sin, sout, serr;
*/
{
	int status;
    pid_t pid;

	if ((pid = fork()) < 0) {
		perror("fork");
		exit( 1 );
	}

	if( pid == 0 ) {
		if( sin != 0 ) {
			close( 0 );  dup( sin );
		}
		if( sout != 1 ) {
			close( 1 );  dup( sout );
		}
		if( serr != 2 ) {
			close( 2 );  dup( serr );
		}

		execvp( *args, args );
		perror( *args );
		exit( 1 );
	}

	while( wait(&status) != pid ) ;
	return( status );
}
