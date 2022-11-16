#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
// Need conditionals for Linux, Cygwyn, Solaris, etc. so
// getopt_long and getopt_long_only prototypes can be found.
#include "nabcode.h"
#include "get_the_opts.h"

void
usage()
{
    fprintf(stderr, "Usage: peptide [-o|-outf string] [-l|-lib string] %s",
        "-t|-stype string -s|-seq string [-?|--help] [--usage]\n");
    fprintf(stderr, "\t-t|-stype string -s|-seq string are required\n");
}

static STRING_T *buf;
static STRING_T *s_type;
static STRING_T *sequence;
static STRING_T *ofile;

INT_T
get_the_opts(INT_T *ac, STRING_T * const *av, struct opt_t *o)
{
    static struct option longopts[] = {
	{ "outf",  required_argument, NULL, 'o' },
	{ "lib",   required_argument, NULL, 'l' },
	{ "stype", required_argument, NULL, 't' },
	{ "seq",   required_argument, NULL, 's' },
	{ "help",  no_argument,       NULL, 'h' },
	{ "usage", no_argument,       NULL, 'u' },
	{ NULL,    0,                 NULL,  0  }
    };
    int ch;

    buf      = malloc( MAXSTRINGLENGTH * sizeof(STRING_T) );
	if( buf == NULL ){
		fprintf( stderr, "bad malloc\n" );
		exit(1);
	}
    s_type   = malloc( MAXSTRINGLENGTH * sizeof(STRING_T) );
    sequence = malloc( MAXSTRINGLENGTH * sizeof(STRING_T) );
    ofile    = malloc( MAXSTRINGLENGTH * sizeof(STRING_T) );

    // This is the default conf.lib file name
    strcpy(buf, getenv("AMBERHOME"));
	fprintf( stderr, "AMBERHOME is %s\n", buf );
    strcat(buf, "/dat/reslib/conf.lib");
    // Default output is to stdout.
    // (Write to /dev/stdout if no file name is given.)
    // (May need to mount fdescfs(5) if /dev/stdout is not available.
    //  Maybe only for BSD derivatives.)
    strcpy(ofile, "/dev/stdout");

    strcpy(sequence, "");
    strcpy(s_type, "");

    while ((ch = getopt_long_only(*ac, av, "o:l:t:s:", longopts, NULL)) != -1)
	switch (ch) {
	case 'o':
	    strcpy(ofile, optarg);
	    break;
	case 'l':
	    strcpy(buf, optarg);
	    break;
	case 't':
	    strcpy(s_type, optarg);
	    break;
	case 's':
	    strcpy(sequence, optarg);
	    break;
	case 'h':
	case 'u':
	default:
	    usage();
	    exit(1);
	}
    *ac -= optind;
     av += optind;

    if (strlen(sequence) == 0 || strlen(s_type) == 0) {
	usage();
	exit(1);
    }

#if defined(DEBUG)
    fprintf(stderr, "buf = %s\n", buf);
    fprintf(stderr, "s_type = %s\n", s_type);
    fprintf(stderr, "sequence = %s\n", sequence);
    fprintf(stderr, "ofile = %s\n", ofile);
    exit(0);
#endif
    o->conf_file   = buf;
    o->struct_type = s_type;
    o->seq         = sequence;
    o->outfile     = ofile;
    return 0;
}
