#ifndef lint
static char yysccsid[] = "@(#)yaccpar	1.9 (Berkeley) 02/21/93";
#endif
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define yyclearin (yychar=(-1))
#define yyerrok (yyerrflag=0)
#define YYRECOVERING (yyerrflag!=0)
#define YYPREFIX "yy"
#line 75 "parser.y"
#include	<unistd.h>
#include        "basics.h"

#include        "classes.h"

#include        "dictionary.h"
#include        "parmLib.h"

#include        "commands.h"
#include	"block.h"
#include	"parser.h"

#include        "leap.h"
#include        "block.h"

#include        "help.h"

#define         MESSAGEFILTER   MESSPARSER




#define         NULLSTR         "null"

/*
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        GLOBAL VARIABLES

*/

                /* The global variable GfCurrentInput defines the file */
                /* from which input is currently read.  When the file */
                /* is empty, then switch back to stdin */

#define	MAXINPUT	1000
#define	MAXINPUTFILES	10		/* Maximum 10 input files */
					/* can be open at once */


char            GsInputLine[MAXINPUT] = "";
BOOL		GbLastLine = FALSE;
BOOL		bCmdDeleteObj;
int             GiInputPos = 0;
PARMLIB		GplAllParameters;
RESULTt		GrMainResult;
BLOCK		GbCommand = NULL;
BLOCK		GbExecute = NULL;
int		GiClipPrompts = 0;
BOOL		GbGraphicalEnvironment;
STRING		GsProgramName;

extern int	iMemDebug;

static	STRING	*SbFirstSourceFiles = NULL;
static	int	iFirstSource = 0;
static	BOOL	SbUseStartup = TRUE;



/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *
 *	The following is used by the parser to maintain a stack of
 *	files where input is received from.  If the file is NULL
 *	then input is received from the main program in the form
 *	of BLOCKS.  The main program is then responsible for reading
 *	the stdin ( in the command line interface ) or for gathering
 *	keypress events from X-Windows ( in the graphical interface ).
 *
 */


int             GiInputFileStackPos = 0;
FILE*           GfaInputFileStack[MAXINPUTFILES];



/*
 *----------------------------------------------------------------
 *
 *	Not quite GLOBAL variables used by the parser.
 */

                                /* Arguments to routines are passed through */
                                /* an array */
#define MAXARGS         10
#define MAXLISTNEXT     10


ATOM            aDummy;
ASSOC           aaArgs[MAXARGS];
int             iArgCount, i;

                                /* List stuff is used for input of nested */
                                /* lists */
#define MAXLISTNEST     10
ASSOC           aaLists[MAXLISTNEST];
int             iCurrentList = -1;
#define PUSHLIST()      iCurrentList++
#define POPLIST()       iCurrentList--
#define CURRENTLIST     aaLists[iCurrentList]


OBJEKT          o0;
double          dTemp;
ASSOC           aAssoc;
STRING          sTemp;
BOOL		bQuit = FALSE;
BOOL		bCommandFound = FALSE;

                /* There seems to be a problem with YACC not properly */
                /* declaring yylval and yyval */
typedef struct  {
	ASSOC		aVal;
	double		dVal;
	STRING		sVal;
	FUNCTION	fCallback;
} YYSTYPEt;

#define YYSTYPE YYSTYPEt


extern  OBJEKT  oGetObject();           /* ( STRING ) */
extern  int     yyparse();
#line 138 "y.tab.c"
#define LVARIABLE 257
#define LSTRING 258
#define LNUMBER 259
#define LASSIGN 260
#define LENDOFCOMMAND 261
#define LOPENLIST 262
#define LCLOSELIST 263
#define LOPENPAREN 264
#define LCLOSEPAREN 265
#define LQUIT 266
#define LCOMMA 267
#define LDOT 268
#define LCOMMAND 269
#define LDUMMY 270
#define LNULL 271
#define LNOTSINGLECHAR 272
#define YYERRCODE 256
short yylhs[] = {                                        -1,
    0,    1,    1,    1,    2,    2,    3,    3,    5,    5,
    8,    6,    6,    6,    6,    6,    6,   11,    4,    7,
    7,    9,   10,   10,   10,   12,
};
short yylen[] = {                                         2,
    1,    1,    1,    2,    2,    2,    3,    2,    1,    1,
    0,    4,    1,    1,    1,    1,    1,    0,    3,    0,
    2,    1,    0,    1,    2,    1,
};
short yydefred[] = {                                      0,
    0,    0,    2,   22,    0,    1,    3,    0,    0,   18,
    4,    0,    5,    6,    0,   15,   14,   13,   11,   16,
   17,   10,    7,    9,   26,    0,   24,   20,   25,    0,
   12,   21,
};
short yydgoto[] = {                                       5,
    6,    7,    8,    9,   23,   25,   30,   28,   10,   26,
   15,   27,
};
short yysindex[] = {                                   -252,
 -254, -250,    0,    0,    0,    0,    0, -253, -249,    0,
    0, -256,    0,    0, -231,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0, -231,    0,    0,    0, -238,
    0,    0,
};
short yyrindex[] = {                                      0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0, -245,    0,    0, -239,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0, -232,    0,    0,    0,    0,
    0,    0,
};
short yygindex[] = {                                      0,
    0,    0,    0,   -1,    0,  -12,    0,    0,    0,    0,
    0,   -3,
};
#define YYTABLESIZE 40
short yytable[] = {                                      24,
   16,   17,   18,    1,    2,   19,   11,   13,    3,   12,
   22,   14,    4,   20,   21,    8,    4,   32,   16,   17,
   18,   23,   29,   19,   31,   16,   17,   18,   19,    0,
   19,   20,   21,    0,    0,    0,    0,    0,   20,   21,
};
short yycheck[] = {                                      12,
  257,  258,  259,  256,  257,  262,  261,  261,  261,  260,
   12,  261,  269,  270,  271,  261,  269,   30,  257,  258,
  259,  261,   26,  262,  263,  257,  258,  259,  261,   -1,
  262,  270,  271,   -1,   -1,   -1,   -1,   -1,  270,  271,
};
#define YYFINAL 5
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 272
#if YYDEBUG
char *yyname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"LVARIABLE","LSTRING","LNUMBER",
"LASSIGN","LENDOFCOMMAND","LOPENLIST","LCLOSELIST","LOPENPAREN","LCLOSEPAREN",
"LQUIT","LCOMMA","LDOT","LCOMMAND","LDUMMY","LNULL","LNOTSINGLECHAR",
};
char *yyrule[] = {
"$accept : input",
"input : line",
"line : LENDOFCOMMAND",
"line : instruct",
"line : error LENDOFCOMMAND",
"instruct : assign LENDOFCOMMAND",
"instruct : function LENDOFCOMMAND",
"assign : LVARIABLE LASSIGN express",
"assign : LVARIABLE LASSIGN",
"express : rawexp",
"express : function",
"$$1 :",
"rawexp : LOPENLIST $$1 elements LCLOSELIST",
"rawexp : LNUMBER",
"rawexp : LSTRING",
"rawexp : LVARIABLE",
"rawexp : LDUMMY",
"rawexp : LNULL",
"$$2 :",
"function : cmdname $$2 args",
"elements :",
"elements : elements rawexp",
"cmdname : LCOMMAND",
"args :",
"args : arg",
"args : args arg",
"arg : rawexp",
};
#endif
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 500
#define YYMAXDEPTH 500
#endif
#endif
int yydebug;
int yynerrs;
int yyerrflag;
int yychar;
short *yyssp;
YYSTYPE *yyvsp;
YYSTYPE yyval;
YYSTYPE yylval;
short yyss[YYSTACKSIZE];
YYSTYPE yyvs[YYSTACKSIZE];
#define yystacksize YYSTACKSIZE
#line 403 "parser.y"
/*------------------------------------------------------------

        ROUTINES

*/

static  BOOL    SbGotUngetc = FALSE;
static  char    ScUngetc;


/*
 *      yyerror
 *
 *      Respond to errors.
 */
int
yyerror( char *sStr )
{
    VP0(( "ERROR: %s\n", sStr ));
    return 1;
}

FILE *
fINPUTFILE()            
{
	if ( GiInputFileStackPos < 0 )
		return(NULL);
	return( GfaInputFileStack[GiInputFileStackPos] );
}



/*
 *	zbGetLine
 *
 *	Get the next line of input from the current input source.
 *	This may be GbCommand or GbExecute depending if
 *	the input is coming from the user or from an execute file.
 *	Return FALSE if there is no more lines to be received from 
 *	the BLOCK.
 */
BOOL
zbGetLine( char *sLine, BOOL *bPFromExecute )
{
BOOL		bGotBlock;
char		c;

    if ( fINPUTFILE() == NULL ) {
	*bPFromExecute = FALSE;
	return(bBlockReadLine( GbCommand, sLine ));
    }

    *bPFromExecute = TRUE;

		/* If there is no line in the execute BLOCK */
		/* then read in another block from the current file */

    if ( bBlockEndOfRead(GbExecute) ) {
	    BlockEmpty( GbExecute );
	    bGotBlock = FALSE;
	    while ( !feof(fINPUTFILE()) && !bGotBlock ) {
		c = fgetc(fINPUTFILE());
		if ( feof(fINPUTFILE()) ) break;
		bGotBlock = bBlockAddChar( GbExecute, c );
	    }

		/* If a complete BLOCK was not read then */
		/* append a final '\n' character, if that doesn't */
		/* make this a complete BLOCK then there is an error */
		/* in the input, also we are at the end of the file */
		/* so pop the file off the stack and say that we have */
		/* a complete BLOCK */

	    if ( !bGotBlock ) {
	    	bGotBlock = bBlockAddChar( GbExecute, '\n' );
	    	if ( fINPUTFILE() != NULL )
	    		fclose( fINPUTFILE() );
	    	INPUTPOPFILE();
	    }
    }
    return(bBlockReadLine( GbExecute, sLine ));
}


		
	


/*
 *      zcGetChar
 *
 *      Get the next character from the line buffer.
 *	If there are no more characters in the line buffer and
 *	GbLastLine is TRUE then return '\0', otherwise if the
 *	line buffer is empty, fill it and return the next character
 *	in it.
 */
char
zcGetChar()
{
char            c;
BOOL		bFromExecute;

                /* If there is a pushed character then return it */

    if ( SbGotUngetc ) {
        SbGotUngetc = FALSE;
        c = ScUngetc;
        goto DONE;
    }

		/* Now if the input line is empty then fill it */
    if ( GsInputLine[GiInputPos] == '\0' ) {
	if ( GbLastLine ) {
	    c = '\0';
	    GbLastLine = FALSE;
	    goto DONE;
	}
	GbLastLine = zbGetLine( GsInputLine, &bFromExecute );
	GiInputPos = 0;
	if ( bFromExecute ) {
	    for ( i=GiClipPrompts; i<=iINPUTSTACKDEPTH(); i++ ) VP2(( ">" ));
	    VP2(( " " ));
	    VP2(( "%s", GsInputLine ));
	} else {
	    VPLOG(( "> %s", GsInputLine ));
	}
    }

    c = GsInputLine[GiInputPos++];

DONE:
    return(c);

}



/*
 *      zUngetc
 *
 *      Push one character back to the file.
 */
void
zUngetc( char c )
{
    SbGotUngetc = TRUE;
    ScUngetc = c;
}






/*
 *      cGetChar
 *
 *      Return the next character, skipping over comments
 *      which start with '#' and end with '\n'.
 */
char
cGetChar()
{
char    c;

	while (1) {
        	c = zcGetChar();
        	if ( c == '#' )
            		while ( ( c=zcGetChar() ) != '\n' )  /* Nothing */ ;
		else
			break;
        } 
	return(c);
}

 
int
iOneCharToken( int c )
{
	switch(c) {
		case ';':
        		return(LENDOFCOMMAND);
		case '=':
        		return(LASSIGN);
		case '*':
        		return(LDUMMY);
		case '(':
        		return(LOPENPAREN);
		case ')':
        		return(LCLOSEPAREN);
		case '{':
        		return(LOPENLIST);
		case '}':
        		return(LCLOSELIST);
		default:
        		return(LNOTSINGLECHAR);
	}
}

/*
 *      yylex
 *
 *      Lexical analyzer for the parser.
 *      Read characters from stdin and return the token types
 *      read and place the value read into the global UNION
 *      yylval.
 *
 *      The things that it recognizes are:
 *              LDOUBLE         [-]###.###E## or ###
 *              LSTRING         "xx xx xxx" or '$'everything up to ' ' ',' ';'
 *              commands        xxxxxxx
 *              LVARIABLE       xxxxxxx which are not commands
 *              LTERMINATOR     ;
 *              LASSIGN         =
 *
 *	Modified 17 November 1992 - David A. Rivkin
 *		Added checking the alias table for command matches.
 *	Total rewrite October 1993 - Bill Ross
 *
 */
int
yylex()
{
STRING          sStr;
int             j, iMax, tok;
BOOL            bGotExp, bGotDot;
char            c;
STRING		sCmd;
STRING		sPossibleCmd;

                /* Skip over blanks, tabs, end of lines etc */

    while ( (c=cGetChar())==' '  ||  c=='\t'  ||  c=='\n'  ||  c == ',' );

    if ( c == '\0' ) 
	return(LENDOFCOMMAND);

    /*
     *  Check the 1-character possibilities: , ; = * ( ) { }
     */
    tok = iOneCharToken( c );
    if ( tok != LNOTSINGLECHAR ) {
        MESSAGE(( "Parsed /%c/\n", c ));
	return(tok);
    }

    /*
     *  it isn't a 1-char thing; read in the rest 
     *	and push back the terminating char
     */
    sStr[0] = c;
    for (j=1;;j++) {

	if ( j >= sizeof(STRING) )
	    DFATAL(( "string too long" ));

	c = cGetChar();
	/*
	 *  NULL terminates anything (?)
	 */
	if ( c == '\0' ) {
	    	sStr[j] = '\0';
	    	break;
	}
	/*
	 *  allow anything inside quotes; chop closing quote
	 */
	if ( sStr[0] == '"' ) {
		if ( c == '"' ) {
			sStr[j] = '\0';
			break;
		}
		sStr[j] = c;
		continue;
	}
	/*
	 *  whitespace is a delimiter outside of quotes
	 */
	if ( c == ' '  ||  c == '\t'  ||  c == '\n'  ||  c == ',' ) {
	    	sStr[j] = '\0';
	    	break;
	}
	/*
	 *  special case for $-type strings: allow embedded single-char
	 *	tokens, except ';'
	 */
	if ( sStr[0] == '$' ) {
		if ( c == ';' ) {
			zUngetc( c );
	    		sStr[j] = '\0';
	    		break;
		}
		sStr[j] = c;
		continue;
	}
	tok = iOneCharToken( c );
	if ( tok != LNOTSINGLECHAR ) {
		zUngetc( c );
	    	sStr[j] = '\0';
	    	break;
	}
	sStr[j] = c;
    }

    /*
     *  see if it's a number
     */
    bGotExp = FALSE;
    bGotDot = FALSE;
    if ( isdigit(sStr[0]) || sStr[0] == '-' || sStr[0] == '+' || 
				( sStr[0] == '.' && isdigit(sStr[1]) ) ) {

        for ( j=0; j<sizeof(STRING); j++ ) {
            MESSAGE(( "Thinking NUMBER got: %c\n", sStr[j] ));
	    switch ( sStr[j] ) {
		case '\0':
        	    if ( sscanf( sStr, "%lf", &yylval.dVal ) != 1 ) {
			VP0(( " Couldn't scan NUMBER from (%s)\n", sStr ));
			return(LNULL);
		    }
        	    MESSAGE(( "Parsed a number: %lf\n", yylval.dVal ));
        	    return(LNUMBER);
		case '.':
		    if ( bGotDot ) {
			VP0(( "(Multiple '.' in NUMBER-like thing (%s))\n", 
								sStr ));
			goto notnum;
		    }
		    if ( bGotExp ) {
			VP0(( 
			 "('.' follows exponent in NUMBER-like thing (%s))\n",
								sStr ));
			goto notnum;
		    }
        	    bGotDot = TRUE;
		    break;
		case 'e':
		case 'E':
            	    if ( bGotExp ) {
			VP0(( "(Multiple 'e' in NUMBER-like thing (%s))\n", 
								sStr ));
			goto notnum;
		    }
                    bGotExp = TRUE;
                    break;
		case '+':
		case '-':
                    break;
		default:
            	    if ( !isdigit(sStr[j]) ) {
			goto notnum;
		    }
            	    break;
	    }
	}
    }
notnum:
    /* 
     *  see if it's a string in quotes
     */
    if ( sStr[0] == '"' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE(( "Parsed a STRING: %s\n", sStr ));
        return(LSTRING);
    }

    /* 
     *  see if it's a string prefixed w/ '$'
     */
    if ( sStr[0] == '$' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE(( "Parsed a STRING: %s\n", sStr ));
        return(LSTRING);
    }

                /* LASTLY!!!!!!!! */
    /* 
     *  see if it's a variable/command 
     */
    strcpy( yylval.sVal, sStr );
    strcpy( sPossibleCmd, sStr );
    StringLower( sPossibleCmd );

    		/* Check if there is an alias that is an exact match */
    if ( (iMax = iVarArrayElementCount( GvaAlias )) ) {
	ALIAS		aAlias;
	aAlias = PVAI( GvaAlias, ALIASt, 0 );
	for ( i=0; i<iMax; i++, aAlias++ ) {
	    if ( strcmp( aAlias->sName, sPossibleCmd ) == 0 ) {
	    	strcpy( sPossibleCmd, aAlias->sCommand );
	    }
        }
    }
                /* Check if there is an exact match of the command */
                /* If a command has already been found for this input
                	line, then do not consider the string a command
                	but rather as a STRING variable */
                	
    if ( !bCommandFound ) {
	for ( j=0; strlen(cCommands[j].sName) != 0; j++ ) {
	    strcpy( sCmd, cCommands[j].sName );
	    StringLower( sCmd );
            if ( strcmp( sCmd, sPossibleCmd ) == 0 ) {
		yylval.fCallback = cCommands[j].fCallback;
		MESSAGE(( "Parsed a command: %s\n", sStr ));
		bCommandFound = TRUE;
		return(LCOMMAND);
	    }
        }
    }


                /* If the variable name is null then return LNULL */

    if ( strcmp( sStr, NULLSTR ) == 0 ) 
	return(LNULL);
    
                /* Return the variable name */

    strcpy( yylval.sVal, sStr );
    MESSAGE(( "Parsed a variable: %s\n", sStr ));
    return(LVARIABLE);

}




/*
 *      oGetObject
 *
 *      If the string is a variable then return the OBJEKT that
 *      it is attached to, otherwise if the string is
 *      a string like: 'unit.mol.res.atom' then parse the
 *      individual names and search the CONTAINERS for
 *      the subcontainers.
 */
OBJEKT
oGetObject( char *sName )
{
CONTAINER       cCont[5];
int             j, k, iSeq;
OBJEKT          oObj;
STRING          sLine, sHead;
BOOL		bDot, bAt, bPdbSeq;
STRING		sGroup;
LIST		lGroup;
OSTRING		osString;
LOOP		lRes;
RESIDUE		rRes;

    oObj = oVariable( sName );
    if ( oObj != NULL ) return(oObj);

        /* Now try to parse the name */
    strcpy( sLine, sName );
    bDot = FALSE;
    bAt = FALSE;
    bPdbSeq = FALSE;
    for ( k=0; k<strlen(sName); k++ ) {
	if ( sLine[k] == '.' ) {
	    sLine[k]=' ';
	    bDot = TRUE;
/*fprintf(stderr, "GOTDOT\n"); */
	} else if ( sLine[k] == '@' ) {
	    sLine[k] = ' ';
	    bAt = TRUE;
	} else if ( sLine[k] == '%' ) {
	    if ( !bDot ) {
	        sLine[k] = ' ';
	        bPdbSeq = TRUE;
	    }
	}
    }


    sRemoveFirstString( sLine, sHead );
    cCont[0] = (CONTAINER)oVariable(sHead);
    if ( cCont[0] == NULL ) {
/* fprintf(stderr, "STRING %s\n", sName); */
    	/* It is not an object variable so...
    	   return the whole string as a OSTRING */
	goto String;
    }
    
/*fprintf(stderr, "VARIABLE %c %c\n", bAt, bPdbSeq);*/
    if ( bPdbSeq ) {
	sRemoveLeadingSpaces( sLine );
	sRemoveFirstString( sLine, sHead );
	if ( strlen(sHead) == 0 ) {
	    goto String;
	}
	if ( isdigit(sHead[0]) ) {

			/* Make sure the rest are digits */
			/* If not return NULL */
	    for ( j=1; j<strlen(sHead); j++ ) {
		if ( !isdigit(sHead[j]) ) {
			    /* It is not an object variable so...
			       return the whole string as a OSTRING */
		    goto String;
		}
	    }
	
			/* Find the PDB sequence number */

	    cCont[1] = NULL;	
	    iSeq = atoi(sHead);
	    lRes = lLoop( (OBJEKT)cCont[0], RESIDUES );
	    while ( (rRes = (RESIDUE)oNext(&lRes)) ) {
		if ( iResiduePdbSequence(rRes)==iSeq ) {
		    cCont[1] = (CONTAINER) rRes;
		}
	    }

	    if ( !cCont[1] ) {
		goto String;
	    }

	    if ( bDot ) {
		sRemoveLeadingSpaces( sLine );
		sRemoveFirstString( sLine, sHead );
		if ( strlen(sHead) == 0 ) {
		    goto String;
		}
		if ( isdigit(sHead[0]) ) {

				/* Make sure the rest are digits */
				/* If not return NULL */
		    for ( j=1; j<strlen(sHead); j++ ) {
			if ( !isdigit(sHead[j]) ) {
				    /* It is not an object variable so...
				       return the whole string as a OSTRING */
			    goto String;
			}
		    }
			
		    iSeq = atoi(sHead);
		    cCont[2] = cContainerFindSequence( cCont[1], 
						DIRECTCONTENTS, iSeq );
		} else {
		    cCont[2] = 
			cContainerFindName( cCont[1], DIRECTCONTENTS, sHead );
		}
		return((OBJEKT)cCont[2]);
	    } else {
		return((OBJEKT)cCont[1]);
	    }
	} else {
	    goto String;
	}
    }


    if ( bDot ) {
	k = 0;
	do {
	    sRemoveLeadingSpaces( sLine );
	    sRemoveFirstString( sLine, sHead );
	    if ( strlen(sHead) == 0 ) break;
	    k++;
	    if ( cCont[k-1]->oHeader.cObjType == PARMSETid ) {
		/*  semi-HACK - parmsets don't have contents in this sense */
		cCont[k] = NULL;
	    } else if ( isdigit(sHead[0]) ) {

			    /* Make sure the rest are digits */
			    /* If not return NULL */
		for ( j=1; j<strlen(sHead); j++ ) {
		    if ( !isdigit(sHead[j]) ) {
   				/* It is not an object variable so...
			    	   return the whole string as a OSTRING */
			goto String;
    		    }
		}
		iSeq = atoi(sHead);
		cCont[k] = cContainerFindSequence( cCont[k-1], 
					    DIRECTCONTENTS, iSeq );
	    } else {
		cCont[k] = 
			cContainerFindName( cCont[k-1], DIRECTCONTENTS, sHead );
	    }
	    if ( cCont[k] == NULL ) break;
	} while ( strlen(sHead) != 0 ) ;
        if ( cCont[k] == NULL ) {
    		/* It is not an object variable so...
    		   return the whole string as a OSTRING */
		goto String;
    	}
	return((OBJEKT)cCont[k]);
    }

			/* If group notation then return the group */

    if ( bAt ) {
	sRemoveLeadingSpaces( sLine );
	sRemoveFirstString( sLine, sGroup );
	if ( iObjectType(cCont[0]) != UNITid ) {
    		/* It is not an object variable so...
    		   return the whole string as a OSTRING */
		goto String;
	}

	lGroup = lUnitGroup( (UNIT)cCont[0], sGroup );
	return((OBJEKT)lGroup);
    }

String:

   /* 
    *  It is not an object variable so...
    *	   return the whole string as a OSTRING
    *	   -- need to set refs to 0 so that
    *	      it will be freed later.. HACK
    */

    osString = (OSTRING)oCreate(OSTRINGid);
    OStringDefine( osString, sName );
    ((OBJEKT)(osString))->iReferences = 0;
    return((OBJEKT)osString );
}





    
/*
 *================================================================
 *
 *	Public routines
 */




/*
 *	ParseArguments
 *
 *	Parse the arguments that the user passes to LEaP from
 *	the command line arguments.
 */
void
ParseArguments( int argc, char *argv[] )
{
char		c;
extern	char	*optarg;

    while ( (c = getopt( argc, argv, "hsI:f:" )) != (char)(EOF) ) {
	switch (c) {
	    case 'h':
		printf( "Usage: %s [options]\n", argv[0] );
		printf( "Options:\n" );
		printf( " -h         Generate this message.\n" );
		printf( " -s         Ignore %s startup file.\n", LEAPRC );
		printf( " -I {dir}   Add {dir} to search path.\n" );
		printf( " -f {file}  Source {file}.\n" );
		exit(1);
	    case 's':
		printf( "-s: Ignoring startup file: %s\n", LEAPRC );
		SbUseStartup = FALSE;
		break;
	    case 'I':
		printf( "-I: Adding %s to search path.\n", optarg );
		BasicsAddDirectory( optarg, 1 );
		break;
	    case 'f':
		printf( "-f: Source %s.\n", optarg );
		if ( iFirstSource == 0 ) {
			MALLOC( SbFirstSourceFiles, STRING *, sizeof(STRING) );
			iFirstSource = 1;
		} else {
			iFirstSource++;
			REALLOC( SbFirstSourceFiles, STRING *, SbFirstSourceFiles,
					iFirstSource * sizeof(STRING));
		}
		strcpy( SbFirstSourceFiles[iFirstSource-1], optarg );
		break;
	}
    }
}






/*
 *	ParseInit
 *
 *	Initialize the parser.
 *	If SbStartup is TRUE the execute the LEAPRC script.
 */
void
ParseInit( RESULTt *rPResult )
{
FILE	*fStartup;
int	iFile;

    VP0(( "\nWelcome to LEaP!\n" ));

#ifdef  DEBUG
    VP0(( "LEaP is running in DEBUG mode!\n" ));
#endif

		/* Initialize memory manager debugging */
    INITMEMORYDEBUG();

    HelpInitialize();

                /* Initialize the first file in the stack to be stdin */
                
    GfaInputFileStack[0] = NULL;
    GbExecute = bBlockCreate();

                /* Create a few OBJEKTs that will be used by the parser */

    aDummy = (ATOM)oCreate(ATOMid);
    ContainerSetName( aDummy, "DUMMY" );
    GplAllParameters = plParmLibCreate();

    VariablesInit();
    GrMainResult.iCommand = CNONE;
    rPResult->iCommand = CNONE;
    
                /* Parse the LEAPRC file if bUseStartup is TRUE */

    if ( SbUseStartup ) {
        fStartup = FOPENNOCOMPLAIN( LEAPRC, "r" );
        if ( fStartup == NULL ) {
	    VP0(( "(no %s in search path)\n", LEAPRC ));
	} else {
	    /*
	     *  source the leaprc
	     */
	    VP0(( "Sourcing %s: %s\n", LEAPRC, GsBasicsFullName ));
	    INPUTPUSHFILE( fStartup );
	    GiClipPrompts = 1;
	    while ( fINPUTFILE() != NULL ) {
	        yyparse();
		if ( GrMainResult.iCommand == CQUIT ) {
	    	    if ( fINPUTFILE() != NULL )
	    	    	fclose( fINPUTFILE() );
	    	    INPUTPOPFILE();
		}
	    }
	    *rPResult = GrMainResult;
	    GiClipPrompts = 0;
	}
    }

		/* Parse the first source file specified */
		/* on the command line using the -f option */

    for (iFile=0; iFile<iFirstSource; iFile++) {
	fStartup = FOPENCOMPLAIN( SbFirstSourceFiles[iFile], "r" );
	if ( fStartup != NULL ) {
	    VP0(( "Sourcing: %s\n", GsBasicsFullName ));
	    INPUTPUSHFILE( fStartup );
	    GiClipPrompts = 1;
	    while ( fINPUTFILE() != NULL ) {
		yyparse();
		if ( GrMainResult.iCommand == CQUIT ) {
	    	    if ( fINPUTFILE() != NULL )
	    	    	fclose( fINPUTFILE() );
	    	    INPUTPOPFILE();
		}
	    }
	    *rPResult = GrMainResult;
	    GiClipPrompts = 0;
	}
    }
    if ( SbFirstSourceFiles )
	FREE( SbFirstSourceFiles );
}




/*
 *	ParseBlock
 *
 *	Parse a BLOCK containing one complete command.
 *	Return in rResult the result of the command.
 */
void
ParseBlock( BLOCK bBlock, RESULTt *rPResult )
{
		/* Set up the BLOCK from which to read the command */

    MESSAGE(( "Parsing block: %s\n", sBlockText(bBlock) ));
    GbCommand = bBlock;
    BlockResetRead( GbCommand );

    GrMainResult.iCommand = CNONE;

		/* Parse the BLOCK */
		/* Keep parsing as long as the 'execute' command */
		/* keeps setting the GLOBAL variable GbContinueParsing */

    do {
	yyparse();
	if ( GrMainResult.iCommand == CQUIT ) {
	    if ( fINPUTFILE() != NULL )
	    	fclose( fINPUTFILE() );
	    INPUTPOPFILE();
	}
    } while ( fINPUTFILE() != NULL );

	/* Reset the bCommandFound variable as a new command may be
	   available in a new input */
    bCommandFound = FALSE;
    
	/* Copy the RESULT from the global result variable */

    *rPResult = GrMainResult;

}





/*
 *	ParseShutdown
 *
 *	Shutdown the parser, release all variables setup in ParseInit
 *	Only needed if debugging memory mgt, since prog mem is all
 *	freed when process exits anyway.
 */
void
ParseShutdown()
{
	if ( !iMemDebug )
		return;

	Destroy( (OBJEKT *)&aDummy );
	VariablesDestroy();
	ParmLibDestroy( &GplAllParameters );

	BlockDestroy( &GbExecute );

	HelpShutdown();
}
#line 1120 "y.tab.c"
#define YYABORT goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR goto yyerrlab
int
yyparse()
{
    register int yym, yyn, yystate;
#if YYDEBUG
    register char *yys;
    extern char *getenv();

    if (yys = getenv("YYDEBUG"))
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = (-1);

    yyssp = yyss;
    yyvsp = yyvs;
    *yyssp = yystate = 0;

yyloop:
    if ((yyn = yydefred[yystate]) != 0) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) != 0 && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yyssp >= yyss + yystacksize - 1)
        {
            goto yyoverflow;
        }
        *++yyssp = yystate = yytable[yyn];
        *++yyvsp = yylval;
        yychar = (-1);
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) != 0 && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;
#ifdef lint
    goto yynewerror;
#endif
yynewerror:
    yyerror("syntax error");
#ifdef lint
    goto yyerrlab;
#endif
yyerrlab:
    ++yynerrs;
yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yyssp]) != 0 && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yyssp, yytable[yyn]);
#endif
                if (yyssp >= yyss + yystacksize - 1)
                {
                    goto yyoverflow;
                }
                *++yyssp = yystate = yytable[yyn];
                *++yyvsp = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yyssp);
#endif
                if (yyssp <= yyss) goto yyabort;
                --yyssp;
                --yyvsp;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = (-1);
        goto yyloop;
    }
yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    yyval = yyvsp[1-yym];
    switch (yyn)
    {
case 1:
#line 224 "parser.y"
{
			return 0;
		}
break;
case 3:
#line 231 "parser.y"
{
                        bCommandFound = FALSE;
                        }
break;
case 4:
#line 235 "parser.y"
{
                            VP0(( "\n" ));
                            yyerrok;
                        }
break;
case 7:
#line 246 "parser.y"
{
                                /* Set the value of the variable */
                            aAssoc = yyvsp[0].aVal;
			    if ( aAssoc != NULL ) {
                                MESSAGE(( "Assigning a value to %s\n", 
                                                        yyvsp[-2].sVal ));
                                VariableSet( yyvsp[-2].sVal, oAssocObject(aAssoc) );
				MESSAGE(( "DEREF (assign) - %s\n",
							sAssocName(aAssoc) ));
                                DEREF( aAssoc );
			    } else {
				MESSAGE(("Not assigning value to %s - rmving\n",
							yyvsp[-2].sVal ));
				VariableRemove( yyvsp[-2].sVal );
			    }

                        }
break;
case 8:
#line 264 "parser.y"
{
                            MESSAGE(( "Removing variable %s\n", yyvsp[-1].sVal ));
                            VariableRemove( yyvsp[-1].sVal );
                        }
break;
case 11:
#line 275 "parser.y"
{
                            PUSHLIST();
                                /* Create an ASSOC for the list */
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, oCreate(LISTid) );
                            CURRENTLIST = aAssoc;
                        }
break;
case 12:
#line 284 "parser.y"
{
                            yyval.aVal = yyvsp[-1].aVal;
                            POPLIST();
                        }
break;
case 13:
#line 289 "parser.y"
{
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            o0 = oCreate(ODOUBLEid);
                            ODoubleSet( o0, yyvsp[0].dVal );
                            AssocSetObject( aAssoc, o0 );
                            yyval.aVal = aAssoc;
                        }
break;
case 14:
#line 298 "parser.y"
{
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            o0 = oCreate(OSTRINGid);
                            OStringDefine( (OSTRING) o0, yyvsp[0].sVal );
                            AssocSetObject( aAssoc, o0 );
			    DEREF( o0 );	/* keeps count = 1 */
                            yyval.aVal = aAssoc;
                        }
break;
case 15:
#line 308 "parser.y"
{
			    OBJEKT	o = oGetObject(yyvsp[0].sVal);
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, yyvsp[0].sVal );
                            AssocSetObject( aAssoc, o ); /* REF's o */
                            yyval.aVal = aAssoc;
                        }
break;
case 16:
#line 316 "parser.y"
{
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, aDummy );
                            yyval.aVal = aAssoc;
                        }
break;
case 17:
#line 323 "parser.y"
{
                            MESSAGE(( "Parsed a null\n" ));
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, NULL );
                            yyval.aVal = aAssoc;
                        }
break;
case 18:
#line 334 "parser.y"
{ iArgCount = 0;
                        }
break;
case 19:
#line 337 "parser.y"
{
                                /* Execute the command */
			    MESSAGE(( "executing function\n"));
			    bCmdDeleteObj = FALSE;
                            o0 = yyvsp[-2].fCallback( iArgCount, aaArgs );
			    if ( o0 != NULL ) {
                            	aAssoc = (ASSOC)oCreate(ASSOCid);
                            	AssocSetObject( aAssoc, o0 );
                            	yyval.aVal = aAssoc;
			    } else {
				MESSAGE(( "func == NULL---\n"));
				yyval.aVal = NULL;
			    }

                                /* DEREF each of the arguments */

                            for ( i=0; i<iArgCount; i++ ) {
				if ( bCmdDeleteObj ) {
					MESSAGE(( "bCmdDeleteObj---\n" ));
					DEREF( oAssocObject( aaArgs[i] ) );
				}
				MESSAGE(( "DEREF (function) - %s\n",
						sAssocName(aaArgs[i])));
                                DEREF( aaArgs[i] );
                            }
                        }
break;
case 21:
#line 368 "parser.y"
{
                                /* Get the element and add it to the list */
                            MESSAGE(( "Adding to list:\n" ));
                            ListAddToEnd( (LIST)oAssocObject(CURRENTLIST), 
							(OBJEKT) yyvsp[0].aVal );
                            yyval.aVal = CURRENTLIST;
                        }
break;
case 26:
#line 386 "parser.y"
{
                            aaArgs[iArgCount++] = yyvsp[0].aVal;
                        }
break;
#line 1427 "y.tab.c"
    }
    yyssp -= yym;
    yystate = *yyssp;
    yyvsp -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yyssp = YYFINAL;
        *++yyvsp = yyval;
        if (yychar < 0)
        {
            if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yyssp, yystate);
#endif
    if (yyssp >= yyss + yystacksize - 1)
    {
        goto yyoverflow;
    }
    *++yyssp = yystate;
    *++yyvsp = yyval;
    goto yyloop;
yyoverflow:
    yyerror("yacc stack overflow");
yyabort:
    return (1);
yyaccept:
    return (0);
}
