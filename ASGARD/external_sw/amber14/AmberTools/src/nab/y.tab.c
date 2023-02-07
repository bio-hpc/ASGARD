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
#line 2 "nabgrm.y"
#include <stdio.h>
#include "nab.h"
#include "cgen.h"
#include "errormsg.h"

extern	VALUE_T	val;
static	VALUE_T	v_type;

int yyerror();

typedef	union	{
	int	ival;
	NODE_T	*npval;
} YYSTYPE;

# define YYSTYPE_IS_DECLARED 1

#line 30 "y.tab.c"
#define SYM_ADDRESS 257
#define SYM_ALLOCATE 258
#define SYM_AND 259
#define SYM_ASSERT 260
#define SYM_ASSIGN 261
#define SYM_ATOM 262
#define SYM_ATSIGN 263
#define SYM_ATTRIBUTE 264
#define SYM_BOUNDS 265
#define SYM_BREAK 266
#define SYM_CALL 267
#define SYM_COMMA 268
#define SYM_CONTINUE 269
#define SYM_DEALLOCATE 270
#define SYM_DEBUG 271
#define SYM_DECL 272
#define SYM_DELETE 273
#define SYM_DONT_MATCH 274
#define SYM_DYNAMIC 275
#define SYM_ELSE 276
#define SYM_EQUAL 277
#define SYM_ERROR 278
#define SYM_FILE 279
#define SYM_FLOAT 280
#define SYM_FLOAT_LIT 281
#define SYM_FOR 282
#define SYM_FOREACH 283
#define SYM_GREATER 284
#define SYM_GREATER_EQUAL 285
#define SYM_HASHED 286
#define SYM_IDENT 287
#define SYM_IF 288
#define SYM_IN 289
#define SYM_INDEX 290
#define SYM_INDIRECT 291
#define SYM_INT 292
#define SYM_INT_LIT 293
#define SYM_LBRACE 294
#define SYM_LBRACK 295
#define SYM_LESS 296
#define SYM_LESS_EQUAL 297
#define SYM_LIST 298
#define SYM_LPAREN 299
#define SYM_MATCH 300
#define SYM_MATRIX 301
#define SYM_MINUS 302
#define SYM_MINUS_ASSIGN 303
#define SYM_MINUS_MINUS 304
#define SYM_MODULUS 305
#define SYM_MODULUS_ASSIGN 306
#define SYM_MOLECULE 307
#define SYM_NEGATE 308
#define SYM_NOT 309
#define SYM_NOT_EQUAL 310
#define SYM_OR 311
#define SYM_PARM 312
#define SYM_PERIOD 313
#define SYM_PLUS 314
#define SYM_PLUS_ASSIGN 315
#define SYM_PLUS_PLUS 316
#define SYM_POINT 317
#define SYM_POINTS_TO 318
#define SYM_RBRACE 319
#define SYM_RBRACK 320
#define SYM_RESIDUE 321
#define SYM_RETURN 322
#define SYM_RPAREN 323
#define SYM_SEMICOLON 324
#define SYM_SIZE_T 325
#define SYM_SLASH 326
#define SYM_SLASH_ASSIGN 327
#define SYM_STAR 328
#define SYM_STAR_ASSIGN 329
#define SYM_STMTLIST 330
#define SYM_STRING 331
#define SYM_STRING_LIT 332
#define SYM_STRUCT 333
#define SYM_TEST 334
#define SYM_TYPE 335
#define SYM_UPARROW 336
#define SYM_UPARROW_ASSIGN 337
#define SYM_WHILE 338
#define YYERRCODE 256
short yylhs[] = {                                        -1,
    0,   22,   22,   23,   23,   65,   65,   66,   66,   21,
   21,   21,   21,   71,   74,   70,   70,   63,   63,   63,
   63,   63,   63,   63,   63,   63,   63,   63,   68,   68,
   31,   31,   30,   49,   49,   75,   75,   73,   73,   10,
   10,    9,    9,    8,    8,   45,   45,   78,   46,   79,
   47,   38,   38,   42,   42,   41,   80,   81,   82,   44,
   32,   32,   56,   56,   43,   43,   64,   64,   64,   64,
   64,   64,   64,   64,   64,   64,   64,   64,   64,    4,
   11,   14,   83,   15,   17,   19,   20,   24,   28,   51,
   84,   51,   39,   85,   86,   61,   77,   87,   88,   89,
   50,   90,   91,   36,   34,   34,   37,   92,   93,   33,
   35,   35,   40,   40,   94,   95,   96,   76,   18,   18,
   26,   26,   27,   27,   55,   55,   55,    7,   13,   62,
   62,   25,   25,   16,   16,    3,    3,   69,   69,   29,
   29,   59,   59,   59,   59,   59,   59,   59,   53,   53,
    1,    1,    6,    6,    5,   52,   52,   12,   12,   12,
   12,   12,   12,   12,   60,   60,   60,   60,   60,   60,
   60,   60,   60,    2,    2,   57,   57,   57,   57,   54,
   54,   72,   72,   48,   58,   58,   67,
};
short yylen[] = {                                         2,
    2,    1,    0,    1,    2,    1,    0,    1,    2,    1,
    1,    1,    1,    2,    3,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    2,    5,
    1,    2,    3,    1,    3,    1,    3,    1,    4,    1,
    1,    1,    3,    1,    1,    2,    3,    0,    3,    0,
    6,    1,    0,    1,    3,    2,    0,    0,    0,    8,
    1,    0,    1,    2,    1,    0,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    3,
    3,    2,    0,    4,    2,    3,    3,    3,    2,    2,
    0,    5,    2,    0,    0,    5,    2,    0,    0,    0,
    7,    0,    0,    6,    1,    1,    3,    0,    0,    7,
    1,    0,    1,    0,    0,    0,    0,    7,    3,    1,
    1,    3,    1,    3,    1,    1,    1,    4,    3,    1,
    3,    1,    3,    1,    3,    1,    3,    1,    3,    1,
    3,    1,    1,    1,    1,    2,    4,    3,    2,    2,
    1,    0,    1,    3,    1,    1,    3,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,
};
short yydefred[] = {                                      0,
   26,   24,   22,   20,   18,   25,   28,   23,   27,   19,
   21,    0,    0,    0,    0,    2,   12,   13,    0,   16,
    0,    0,   10,   11,  184,    0,    5,    0,    0,    0,
    0,    0,    0,    0,  186,  102,   98,  185,   83,    0,
  182,  181,  183,  180,   94,  187,  115,    0,   68,  126,
   69,  127,   70,   72,    0,   71,   73,   74,   75,    0,
    0,   67,    0,    0,   76,    0,    0,   77,  145,    0,
    0,  143,    0,   78,  123,    8,    1,    0,  144,    0,
    0,    0,   79,   46,    0,    0,   14,    0,    0,    0,
    0,    0,    0,   82,   85,    0,    0,    0,  120,    0,
    0,    0,    0,    0,    0,    0,    0,  172,  167,  170,
  169,  173,  165,  166,  171,  168,    0,    0,    0,   89,
  179,  178,  177,  176,    0,   93,    0,    0,  125,    0,
  158,    0,  160,  163,    0,  159,  162,  161,  164,    0,
  150,    0,    9,  175,  174,    0,    0,  146,   97,   47,
   57,   49,    0,    0,    0,   15,    0,    0,    0,   80,
   81,   86,    0,    0,   87,    0,   88,  103,   99,    0,
  148,   95,  116,  135,  133,  131,  139,    0,    0,  151,
  155,   91,    0,    0,  129,  124,  141,  137,    0,   45,
   41,    0,   40,    0,   44,    0,    0,   37,   32,   30,
    0,    0,  119,  122,    0,    0,   84,    0,    0,  147,
    0,    0,    0,  128,   58,    0,   17,    0,   63,    0,
   39,    0,    0,   52,    0,    0,   33,  111,  106,    0,
    0,  105,    0,  100,   96,  117,  154,   92,  157,    0,
   64,   43,   51,    0,   56,   35,  104,  108,    0,    0,
    0,    0,    0,   55,    0,  107,  101,  118,   59,  113,
    0,    0,  109,   60,    0,  110,
};
short yydgoto[] = {                                      13,
  178,  146,   48,   49,  179,  180,   50,  192,  193,  194,
   51,  140,   52,   53,   54,   55,   56,   98,   57,   58,
   14,   15,   16,   59,   60,   99,   61,   62,   63,  157,
  158,  215,  229,  230,  231,   64,  232,  222,   65,  261,
  223,  224,  252,  152,   17,   18,   19,   66,  202,   67,
   68,  184,   69,   70,   71,  216,  125,   72,   73,  117,
   74,   75,   20,   76,   77,   78,   79,  217,   80,   22,
   23,   81,   89,   24,   90,   82,   83,   86,  154,  189,
  240,  262,  104,  212,  106,  208,  103,  206,  250,  102,
  205,  255,  265,  107,  209,  251,
};
short yysindex[] = {                                    -84,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0, -265,    0,  -84, -132,    0,    0,    0, -270,    0,
 -285, -265,    0,    0,    0, -246,    0,  -18,  -18, -268,
 -263,  -18,  139,  -18,    0,    0,    0,    0,    0,  -18,
    0,    0,    0,    0,    0,    0,    0,   22,    0,    0,
    0,    0,    0,    0, -197,    0,    0,    0,    0, -225,
 -230,    0, -233, -132,    0, -221, -132,    0,    0, -265,
 -204,    0, -240,    0,    0,    0,    0, -132,    0, -282,
  -18, -132,    0,    0, -211, -179,    0, -177, -146, -200,
  490, -194, -193,    0,    0, -192,  -18, -189,    0, -125,
 -180, -152, -142, -132, -172,  -18, -140,    0,    0,    0,
    0,    0,    0,    0,    0,    0,  -18,  -18,  -18,    0,
    0,    0,    0,    0,  -18,    0,  -18, -123,    0, -269,
    0,  -18,    0,    0, -127,    0,    0,    0,    0,  -18,
    0,  -18,    0,    0,    0,  -18, -206,    0,    0,    0,
    0,    0,  -59, -136, -265,    0,  490, -155, -265,    0,
    0,    0, -157, -247,    0,  -18,    0,    0,    0, -235,
    0,    0,    0,    0,    0,    0,    0, -154, -100,    0,
    0,    0,  -97, -147,    0,    0,    0,    0,  -84,    0,
    0,  -94,    0, -145,    0,  -84, -177,    0,    0,    0,
  -92, -138,    0,    0,  -18,  -18,    0, -135,  -18,    0,
  -18, -132,  -18,    0,    0,  -84,    0, -265,    0,  586,
    0, -141,  -77,    0, -265, -265,    0,    0,    0, -131,
 -120,    0, -170,    0,    0,    0,    0,    0,    0, -132,
    0,    0,    0,  -84,    0,    0,    0,    0, -265, -130,
 -126, -114, -132,    0,  -18,    0,    0,    0,    0,    0,
 -117, -112,    0,    0,  -18,    0,
};
short yyrindex[] = {                                     55,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,   92,  188,    0,    0,    0,  -81,    0,
  -73,    0,    0,    0,    0, -258,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0, -241,    0,    0,
    0,    0,    0,    0, -203,    0,    0,    0,    0, -178,
    0,    0,  516,    0,    0,  238,    0,    0,    0,    0,
  367,    0,  457,    0,    0,    0,    0,  210,    0,  560,
    0,    0,    0,    0,    0,    0,    0, -249, -109,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0, -281,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0, -103,    1,    0,  412,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,  367,    0,    0,    0,
    0,    0,    0,    0,    0,    0,  -98,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,  -99,    0,
    0,    0,  -95,    0,    0,    0,    0,    0,  140,    0,
    0,  -89,    0,    0,    0,  -91, -219,    0,    0,    0,
  -88,    0,    0,    0,  -86,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,  194,    0,    0,    0,    0,
    0,    0,  -79,    0,    0,    0,    0,    0,    0,    0,
    0,    0,  296,    0,    0,    0,    0,    0,    0,  -71,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,  -68,    0,  -78,    0,    0,    0,    0,    0,
    0,    0,    0,    0,  -70,    0,
};
short yygindex[] = {                                      0,
    0,    0, -101,    0,    0,   24,    0,    0,   32,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,  240,    0,  137,  -93,  -26,    0,  114,    0,
  101,    0,    0,    0,   -5,    0,    0,    0,    0,    0,
    0,   18,    0,    0,    0,    0,    0,   -7,   38,    0,
    0,   52,    0,  -62,   84,    0,    0,    0,  185,    0,
    0,  149,  -80,  -54,    0, -104,    0,   37,  151, -156,
    0,    0,   53, -148,  122,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,
};
#define YYTABLESIZE 918
short yytable[] = {                                     170,
   90,   92,   93,  163,   26,   96,  100,  101,  141,  126,
  159,   85,  128,  105,   88,  174,   25,  134,   38,  144,
  166,   25,   28,  143,   29,  132,  134,  149,   29,  121,
   30,  145,  218,   31,   32,   33,   21,   34,   87,  225,
  219,  121,  121,  135,  188,   35,   36,   91,   38,   50,
   21,   25,   37,   84,    3,   94,  131,   38,   39,  218,
   95,  118,  129,   40,  132,   29,   41,  241,   42,  134,
  164,  122,  204,   43,   38,  171,  159,  127,  134,  172,
   44,  134,  134,  207,  141,  119,   45,  225,  132,  130,
  132,    4,  123,  120,  124,  142,   46,   42,  133,   42,
  181,  134,   47,   38,   38,  183,  135,  132,  135,   44,
  136,   44,  150,  186,  151,  143,  132,  153,  249,  132,
  132,  155,  137,  156,  138,   28,  195,   29,  127,  160,
  161,  162,  139,   30,  165,  253,   31,   32,   33,  100,
   34,  130,  166,  167,  130,  130,  168,  197,   35,   36,
  171,  201,  182,  130,   25,   37,  169,  238,  173,  185,
   38,   39,  196,  200,  147,  203,   40,  211,  210,   41,
  213,   42,  214,  220,  221,  226,   43,    1,  228,  234,
    2,  243,  236,   44,  181,  227,  183,    7,  235,   45,
  244,  247,  257,  195,    3,    4,  258,  233,  143,   46,
  147,  147,  147,  248,  259,   47,  263,    5,  147,    6,
  197,  264,   48,   17,   36,  190,    6,  197,  201,  152,
   31,   35,    7,  153,  156,  147,  191,   25,  260,  147,
   42,   53,    8,   38,  237,   34,    9,  112,  228,   40,
   10,  256,   41,   54,   42,  114,   11,   66,   12,   43,
   65,  242,  112,   27,  175,  187,   44,  199,   90,  266,
   90,  254,   35,  246,  239,  148,   90,  176,   25,   90,
   90,   90,   46,   90,   38,  177,  198,  245,    0,    0,
   40,   90,   90,   41,    0,   42,    0,   90,   90,    0,
   43,    0,    0,   90,   90,  108,    0,   44,  109,   90,
    0,    0,   90,    0,   90,  110,  111,    0,    0,   90,
  112,    0,    3,   46,    3,    0,   90,  113,  114,   90,
    3,  115,   90,    3,    3,    3,    0,    3,    0,    0,
    0,  116,   90,    0,    0,    3,    3,    0,   90,    0,
    0,    3,    3,    0,    0,    0,    0,    3,    3,    4,
    0,    4,    0,    3,    0,    0,    3,    4,    3,    0,
    4,    4,    4,    3,    4,    0,    0,    0,    0,    0,
    3,    0,    4,    4,    0,    0,    3,    0,    4,    4,
    0,    0,    0,    0,    4,    4,    3,    0,    0,    0,
    4,    0,    3,    4,    0,    4,    0,   62,    0,   62,
    4,    0,    0,    0,    0,   62,    0,    4,   62,   62,
   62,    0,   62,    4,    0,    0,    0,    0,    0,   35,
   62,   62,    0,    4,    0,   25,   62,   62,    0,    4,
    0,   38,   62,   62,    0,    0,    0,   97,   62,    0,
   41,   62,   42,   62,    0,    0,    0,   43,   62,    0,
    0,   61,    0,   61,   44,   62,    0,    0,   62,   61,
    0,   62,   61,   61,   61,    0,   61,    0,    0,    0,
   46,   62,    0,    0,   61,   61,    0,   62,    0,    0,
   61,   61,    0,    0,    0,    0,   61,   61,    0,    0,
    0,    0,   61,    0,    0,   61,  125,   61,  125,    0,
  125,    0,   61,    0,    0,  125,    0,    0,    0,   61,
    0,  125,   61,    0,  125,   61,    0,    0,    0,    0,
    0,  125,  125,    0,    0,   61,  125,    0,    0,    0,
    0,   61,  125,  125,  125,    0,    0,  125,    0,  125,
  125,  125,  125,  125,    0,    0,    0,  125,  125,    0,
  125,  125,  125,  125,  125,    0,  125,  125,  125,    0,
  125,  125,    0,  125,  125,  125,  125,    0,    0,  125,
    0,    0,  125,  125,  125,    0,    0,    0,    0,  125,
  125,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  125,  125,  125,    0,    0,  125,    0,  125,  125,  125,
  125,  125,    0,    0,    0,  125,  125,    0,  125,  125,
  125,  125,    0,    0,    0,    0,    0,    0,    0,  125,
    0,  125,  125,  125,  125,  142,    0,    0,    0,  142,
    0,  125,  125,    0,  142,    0,    0,    0,    0,    0,
  142,    0,    0,  142,    0,    0,    0,    0,    0,    0,
  142,  142,    0,    0,    0,  142,    0,    0,    0,    0,
    0,    0,  142,  142,    0,    0,  142,    0,  142,    0,
  149,  142,    0,    0,  149,    0,  142,  142,    0,  149,
  142,    0,    0,    0,    0,  149,  142,    0,  149,  142,
  142,    0,  142,    0,  142,  149,  149,    0,    0,    0,
  149,    0,  142,    0,    0,    0,    0,  149,  149,    0,
    0,  149,    0,  149,    0,  140,  149,    0,    0,  140,
    0,  149,  149,    0,  140,  149,    0,    0,    0,    0,
  140,  149,    0,  140,  149,  149,    0,  149,    0,  149,
  140,  140,    0,    0,    0,  140,    0,  149,    0,    0,
    0,    1,  140,  140,    2,    0,  140,    0,  140,    0,
    0,  140,    0,    0,    0,    0,  140,  140,    3,    4,
  140,    0,    0,    0,  138,    0,  140,    0,    0,  140,
  140,    5,  140,  138,  140,    0,    0,    0,    0,  138,
    6,    0,  138,    0,    0,    0,    7,    0,    0,  138,
  138,    0,    0,    0,  138,    0,    8,    0,    0,    0,
    9,  138,  138,    0,   10,  138,    0,  138,  136,    0,
   11,    0,    0,    0,    0,  138,  138,  136,    0,  138,
    0,    0,    0,  136,    0,  138,  136,    0,  138,  138,
    0,    0,    0,  136,  136,    0,    0,    0,  136,    0,
    0,    0,    0,    0,    0,  136,  136,    0,    0,  136,
  190,    0,    0,    0,    0,    0,   35,    0,    0,  136,
  136,    0,   25,    0,    0,    0,    0,    0,   38,  136,
    0,    0,  136,  136,   40,    0,    0,   41,    0,   42,
    0,    0,    0,    0,   43,    0,    0,    0,    0,    0,
    0,   44,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,   46,
};
short yycheck[] = {                                     104,
    0,   28,   29,   97,   12,   32,   33,   34,   71,   64,
   91,   19,   67,   40,   22,  117,  287,  259,  268,  302,
  268,  287,  258,   78,  260,  295,  268,   82,  287,  263,
  266,  314,  189,  269,  270,  271,    0,  273,  324,  196,
  189,  323,  324,  313,  146,  281,  282,  294,  268,  299,
   14,  287,  288,  324,    0,  324,  261,  293,  294,  216,
  324,  259,   70,  299,  268,  324,  302,  216,  304,  311,
   97,  305,  166,  309,  324,  323,  157,  299,  320,  106,
  316,  323,  324,  319,  147,  311,  322,  244,  295,  268,
  295,    0,  326,  324,  328,  336,  332,  304,  303,  304,
  127,  306,  338,  323,  324,  132,  313,  311,  313,  316,
  315,  316,  324,  140,  294,  170,  320,  295,  289,  323,
  324,  268,  327,  324,  329,  258,  153,  260,  299,  324,
  324,  324,  337,  266,  324,  240,  269,  270,  271,  166,
  273,  320,  268,  324,  323,  324,  299,  155,  281,  282,
  323,  159,  276,   70,  287,  288,  299,  212,  299,  287,
  293,  294,  299,  319,   81,  323,  299,  268,  323,  302,
  268,  304,  320,  268,  320,  268,  309,  262,  205,  206,
  265,  323,  209,  316,  211,  324,  213,    0,  324,  322,
  268,  323,  323,  220,  279,  280,  323,  205,  253,  332,
  117,  118,  119,  324,  319,  338,  324,  292,  125,    0,
  218,  324,  294,  287,  324,  275,  301,  225,  226,  323,
  319,  281,  307,  323,  320,  142,  286,  287,  255,  146,
  320,  323,  317,  293,  211,  324,  321,  324,  265,  299,
  325,  249,  302,  323,  304,  324,  331,  319,  333,  309,
  319,  220,  323,   14,  118,  142,  316,  157,  258,  265,
  260,  244,  281,  226,  213,   81,  266,  119,  287,  269,
  270,  271,  332,  273,  293,  125,  155,  225,   -1,   -1,
  299,  281,  282,  302,   -1,  304,   -1,  287,  288,   -1,
  309,   -1,   -1,  293,  294,  274,   -1,  316,  277,  299,
   -1,   -1,  302,   -1,  304,  284,  285,   -1,   -1,  309,
  289,   -1,  258,  332,  260,   -1,  316,  296,  297,  319,
  266,  300,  322,  269,  270,  271,   -1,  273,   -1,   -1,
   -1,  310,  332,   -1,   -1,  281,  282,   -1,  338,   -1,
   -1,  287,  288,   -1,   -1,   -1,   -1,  293,  294,  258,
   -1,  260,   -1,  299,   -1,   -1,  302,  266,  304,   -1,
  269,  270,  271,  309,  273,   -1,   -1,   -1,   -1,   -1,
  316,   -1,  281,  282,   -1,   -1,  322,   -1,  287,  288,
   -1,   -1,   -1,   -1,  293,  294,  332,   -1,   -1,   -1,
  299,   -1,  338,  302,   -1,  304,   -1,  258,   -1,  260,
  309,   -1,   -1,   -1,   -1,  266,   -1,  316,  269,  270,
  271,   -1,  273,  322,   -1,   -1,   -1,   -1,   -1,  281,
  281,  282,   -1,  332,   -1,  287,  287,  288,   -1,  338,
   -1,  293,  293,  294,   -1,   -1,   -1,  299,  299,   -1,
  302,  302,  304,  304,   -1,   -1,   -1,  309,  309,   -1,
   -1,  258,   -1,  260,  316,  316,   -1,   -1,  319,  266,
   -1,  322,  269,  270,  271,   -1,  273,   -1,   -1,   -1,
  332,  332,   -1,   -1,  281,  282,   -1,  338,   -1,   -1,
  287,  288,   -1,   -1,   -1,   -1,  293,  294,   -1,   -1,
   -1,   -1,  299,   -1,   -1,  302,  259,  304,  261,   -1,
  263,   -1,  309,   -1,   -1,  268,   -1,   -1,   -1,  316,
   -1,  274,  319,   -1,  277,  322,   -1,   -1,   -1,   -1,
   -1,  284,  285,   -1,   -1,  332,  289,   -1,   -1,   -1,
   -1,  338,  295,  296,  297,   -1,   -1,  300,   -1,  302,
  303,  304,  305,  306,   -1,   -1,   -1,  310,  311,   -1,
  313,  314,  315,  316,  259,   -1,  261,  320,  263,   -1,
  323,  324,   -1,  326,  327,  328,  329,   -1,   -1,  274,
   -1,   -1,  277,  336,  337,   -1,   -1,   -1,   -1,  284,
  285,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
  295,  296,  297,   -1,   -1,  300,   -1,  302,  303,  304,
  305,  306,   -1,   -1,   -1,  310,  311,   -1,  313,  314,
  315,  316,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  324,
   -1,  326,  327,  328,  329,  259,   -1,   -1,   -1,  263,
   -1,  336,  337,   -1,  268,   -1,   -1,   -1,   -1,   -1,
  274,   -1,   -1,  277,   -1,   -1,   -1,   -1,   -1,   -1,
  284,  285,   -1,   -1,   -1,  289,   -1,   -1,   -1,   -1,
   -1,   -1,  296,  297,   -1,   -1,  300,   -1,  302,   -1,
  259,  305,   -1,   -1,  263,   -1,  310,  311,   -1,  268,
  314,   -1,   -1,   -1,   -1,  274,  320,   -1,  277,  323,
  324,   -1,  326,   -1,  328,  284,  285,   -1,   -1,   -1,
  289,   -1,  336,   -1,   -1,   -1,   -1,  296,  297,   -1,
   -1,  300,   -1,  302,   -1,  259,  305,   -1,   -1,  263,
   -1,  310,  311,   -1,  268,  314,   -1,   -1,   -1,   -1,
  274,  320,   -1,  277,  323,  324,   -1,  326,   -1,  328,
  284,  285,   -1,   -1,   -1,  289,   -1,  336,   -1,   -1,
   -1,  262,  296,  297,  265,   -1,  300,   -1,  302,   -1,
   -1,  305,   -1,   -1,   -1,   -1,  310,  311,  279,  280,
  314,   -1,   -1,   -1,  259,   -1,  320,   -1,   -1,  323,
  324,  292,  326,  268,  328,   -1,   -1,   -1,   -1,  274,
  301,   -1,  277,   -1,   -1,   -1,  307,   -1,   -1,  284,
  285,   -1,   -1,   -1,  289,   -1,  317,   -1,   -1,   -1,
  321,  296,  297,   -1,  325,  300,   -1,  302,  259,   -1,
  331,   -1,   -1,   -1,   -1,  310,  311,  268,   -1,  314,
   -1,   -1,   -1,  274,   -1,  320,  277,   -1,  323,  324,
   -1,   -1,   -1,  284,  285,   -1,   -1,   -1,  289,   -1,
   -1,   -1,   -1,   -1,   -1,  296,  297,   -1,   -1,  300,
  275,   -1,   -1,   -1,   -1,   -1,  281,   -1,   -1,  310,
  311,   -1,  287,   -1,   -1,   -1,   -1,   -1,  293,  320,
   -1,   -1,  323,  324,  299,   -1,   -1,  302,   -1,  304,
   -1,   -1,   -1,   -1,  309,   -1,   -1,   -1,   -1,   -1,
   -1,  316,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,  332,
};
#define YYFINAL 13
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 338
#if YYDEBUG
char *yyname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"SYM_ADDRESS","SYM_ALLOCATE",
"SYM_AND","SYM_ASSERT","SYM_ASSIGN","SYM_ATOM","SYM_ATSIGN","SYM_ATTRIBUTE",
"SYM_BOUNDS","SYM_BREAK","SYM_CALL","SYM_COMMA","SYM_CONTINUE","SYM_DEALLOCATE",
"SYM_DEBUG","SYM_DECL","SYM_DELETE","SYM_DONT_MATCH","SYM_DYNAMIC","SYM_ELSE",
"SYM_EQUAL","SYM_ERROR","SYM_FILE","SYM_FLOAT","SYM_FLOAT_LIT","SYM_FOR",
"SYM_FOREACH","SYM_GREATER","SYM_GREATER_EQUAL","SYM_HASHED","SYM_IDENT",
"SYM_IF","SYM_IN","SYM_INDEX","SYM_INDIRECT","SYM_INT","SYM_INT_LIT",
"SYM_LBRACE","SYM_LBRACK","SYM_LESS","SYM_LESS_EQUAL","SYM_LIST","SYM_LPAREN",
"SYM_MATCH","SYM_MATRIX","SYM_MINUS","SYM_MINUS_ASSIGN","SYM_MINUS_MINUS",
"SYM_MODULUS","SYM_MODULUS_ASSIGN","SYM_MOLECULE","SYM_NEGATE","SYM_NOT",
"SYM_NOT_EQUAL","SYM_OR","SYM_PARM","SYM_PERIOD","SYM_PLUS","SYM_PLUS_ASSIGN",
"SYM_PLUS_PLUS","SYM_POINT","SYM_POINTS_TO","SYM_RBRACE","SYM_RBRACK",
"SYM_RESIDUE","SYM_RETURN","SYM_RPAREN","SYM_SEMICOLON","SYM_SIZE_T",
"SYM_SLASH","SYM_SLASH_ASSIGN","SYM_STAR","SYM_STAR_ASSIGN","SYM_STMTLIST",
"SYM_STRING","SYM_STRING_LIT","SYM_STRUCT","SYM_TEST","SYM_TYPE","SYM_UPARROW",
"SYM_UPARROW_ASSIGN","SYM_WHILE",
};
char *yyrule[] = {
"$accept : program",
"program : defpart stmtpart",
"defpart : defs",
"defpart :",
"defs : def",
"defs : def defs",
"stmtpart : stmts",
"stmtpart :",
"stmts : stmt",
"stmts : stmts stmt",
"def : type_decl",
"def : var_decl",
"def : func_decl",
"def : func_def",
"type_decl : struct_type SYM_SEMICOLON",
"var_decl : type var_list SYM_SEMICOLON",
"type : simple_type",
"type : struct_type",
"simple_type : SYM_INT",
"simple_type : SYM_SIZE_T",
"simple_type : SYM_FLOAT",
"simple_type : SYM_STRING",
"simple_type : SYM_FILE",
"simple_type : SYM_POINT",
"simple_type : SYM_BOUNDS",
"simple_type : SYM_MATRIX",
"simple_type : SYM_ATOM",
"simple_type : SYM_RESIDUE",
"simple_type : SYM_MOLECULE",
"struct_type : SYM_STRUCT id",
"struct_type : SYM_STRUCT id SYM_LBRACE field_list SYM_RBRACE",
"field_list : field",
"field_list : field field_list",
"field : simple_type id_list SYM_SEMICOLON",
"id_list : id",
"id_list : id SYM_COMMA id_list",
"var_list : var",
"var_list : var SYM_COMMA var_list",
"var : id",
"var : id SYM_LBRACK aspec SYM_RBRACK",
"aspec : as_list",
"aspec : SYM_HASHED",
"as_list : asize",
"as_list : asize SYM_COMMA as_list",
"asize : expr",
"asize : SYM_DYNAMIC",
"func_decl : func_hdr SYM_SEMICOLON",
"func_decl : func_hdr id SYM_SEMICOLON",
"$$1 :",
"func_def : func_hdr $$1 func_body",
"$$2 :",
"func_hdr : type id $$2 SYM_LPAREN formals SYM_RPAREN",
"formals : fp_list",
"formals :",
"fp_list : f_parm",
"fp_list : f_parm SYM_COMMA fp_list",
"f_parm : type var",
"$$3 :",
"$$4 :",
"$$5 :",
"func_body : SYM_LBRACE $$3 f_defpart $$4 f_stmtpart SYM_RBRACE $$5 SYM_SEMICOLON",
"f_defpart : lv_decls",
"f_defpart :",
"lv_decls : var_decl",
"lv_decls : lv_decls var_decl",
"f_stmtpart : stmts",
"f_stmtpart :",
"stmt : expr_stmt",
"stmt : alloc_stmt",
"stmt : assert_stmt",
"stmt : break_stmt",
"stmt : continue_stmt",
"stmt : cmpd_stmt",
"stmt : dealloc_stmt",
"stmt : debug_stmt",
"stmt : delete_stmt",
"stmt : for_stmt",
"stmt : if_stmt",
"stmt : return_stmt",
"stmt : while_stmt",
"alloc_stmt : SYM_ALLOCATE expr SYM_SEMICOLON",
"assert_stmt : SYM_ASSERT expr SYM_SEMICOLON",
"break_stmt : SYM_BREAK SYM_SEMICOLON",
"$$6 :",
"cmpd_stmt : SYM_LBRACE $$6 stmts SYM_RBRACE",
"continue_stmt : SYM_CONTINUE SYM_SEMICOLON",
"dealloc_stmt : SYM_DEALLOCATE expr SYM_SEMICOLON",
"debug_stmt : SYM_DEBUG dbg_list SYM_SEMICOLON",
"delete_stmt : SYM_DELETE expr SYM_SEMICOLON",
"expr_stmt : expr SYM_SEMICOLON",
"if_stmt : if_hdr stmt",
"$$7 :",
"if_stmt : if_hdr stmt SYM_ELSE $$7 stmt",
"for_stmt : for_hdr stmt",
"$$8 :",
"$$9 :",
"return_stmt : SYM_RETURN $$8 expr $$9 SYM_SEMICOLON",
"while_stmt : while_hdr stmt",
"$$10 :",
"$$11 :",
"$$12 :",
"if_hdr : SYM_IF $$10 SYM_LPAREN $$11 expr $$12 SYM_RPAREN",
"$$13 :",
"$$14 :",
"for_hdr : SYM_FOR $$13 SYM_LPAREN $$14 for_ctrl SYM_RPAREN",
"for_ctrl : for_in",
"for_ctrl : for_count",
"for_in : id SYM_IN id",
"$$15 :",
"$$16 :",
"for_count : for_expr SYM_SEMICOLON $$15 for_test_expr SYM_SEMICOLON $$16 for_expr",
"for_expr : expr",
"for_expr :",
"for_test_expr : expr",
"for_test_expr :",
"$$17 :",
"$$18 :",
"$$19 :",
"while_hdr : SYM_WHILE $$17 SYM_LPAREN $$18 expr $$19 SYM_RPAREN",
"dbg_list : SYM_LPAREN e_list SYM_RPAREN",
"dbg_list : e_list",
"e_list : expr",
"e_list : expr SYM_COMMA e_list",
"expr : rval",
"expr : lval assignop expr",
"lval : id",
"lval : ar_lval",
"lval : at_lval",
"ar_lval : lval SYM_LBRACK i_list SYM_RBRACK",
"at_lval : lval SYM_PERIOD SYM_IDENT",
"rval : disj",
"rval : disj SYM_OR rval",
"disj : conj",
"disj : conj SYM_AND disj",
"conj : a_expr",
"conj : a_expr relop a_expr",
"a_expr : term",
"a_expr : term addop a_expr",
"term : factor",
"term : factor mulop term",
"factor : primary",
"factor : primary SYM_UPARROW factor",
"primary : lval",
"primary : num",
"primary : string",
"primary : incr",
"primary : unop primary",
"primary : id SYM_LPAREN actuals SYM_RPAREN",
"primary : SYM_LPAREN expr SYM_RPAREN",
"incr : incrop lval",
"incr : lval incrop",
"actuals : ap_list",
"actuals :",
"ap_list : a_parm",
"ap_list : a_parm SYM_COMMA ap_list",
"a_parm : expr",
"i_list : expr",
"i_list : expr SYM_COMMA i_list",
"assignop : SYM_ASSIGN",
"assignop : SYM_PLUS_ASSIGN",
"assignop : SYM_MINUS_ASSIGN",
"assignop : SYM_STAR_ASSIGN",
"assignop : SYM_SLASH_ASSIGN",
"assignop : SYM_MODULUS_ASSIGN",
"assignop : SYM_UPARROW_ASSIGN",
"relop : SYM_LESS",
"relop : SYM_LESS_EQUAL",
"relop : SYM_EQUAL",
"relop : SYM_NOT_EQUAL",
"relop : SYM_GREATER_EQUAL",
"relop : SYM_GREATER",
"relop : SYM_MATCH",
"relop : SYM_DONT_MATCH",
"relop : SYM_IN",
"addop : SYM_PLUS",
"addop : SYM_MINUS",
"mulop : SYM_STAR",
"mulop : SYM_SLASH",
"mulop : SYM_MODULUS",
"mulop : SYM_ATSIGN",
"incrop : SYM_PLUS_PLUS",
"incrop : SYM_MINUS_MINUS",
"unop : SYM_MINUS",
"unop : SYM_NOT",
"id : SYM_IDENT",
"num : SYM_INT_LIT",
"num : SYM_FLOAT_LIT",
"string : SYM_STRING_LIT",
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
#line 507 "nabgrm.y"

#include "lex.yy.c"
#line 702 "y.tab.c"
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
#line 184 "nabgrm.y"
{ CG_genend(); }
break;
case 3:
#line 186 "nabgrm.y"
{ yyval.npval = NULL; }
break;
case 7:
#line 190 "nabgrm.y"
{ yyval.npval = NULL; }
break;
case 14:
#line 200 "nabgrm.y"
{ yyval.npval = node( SYM_DECL, 0, yyvsp[-1].npval, 0 );
				  CG_genvardecl( yyval.npval, 0, 0, 0 ); }
break;
case 15:
#line 203 "nabgrm.y"
{ yyval.npval = node( SYM_DECL, 0, yyvsp[-2].npval, yyvsp[-1].npval );
				  CG_genvardecl( yyval.npval, 0, 0, 0 ); }
break;
case 16:
#line 205 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 17:
#line 206 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 18:
#line 207 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_INT;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 19:
#line 211 "nabgrm.y"
{
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_SIZE_T;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 20:
#line 215 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_FLOAT;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 21:
#line 219 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_STRING;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 22:
#line 223 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_FILE;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 23:
#line 227 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_POINT;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 24:
#line 231 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_BOUNDS;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 25:
#line 235 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_MATRIX;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 26:
#line 239 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_ATOM;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 27:
#line 243 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_RESIDUE;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 28:
#line 247 "nabgrm.y"
{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_MOLECULE;
			yyval.npval = node( SYM_TYPE, &v_type, 0, 0 ); }
break;
case 29:
#line 251 "nabgrm.y"
{
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_USER;
			yyval.npval = node( SYM_TYPE, &v_type, node( SYM_STRUCT, 0, yyvsp[0].npval, 0 ), 0 ); }
break;
case 30:
#line 255 "nabgrm.y"
{
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_USER;
			yyval.npval = node( SYM_TYPE, &v_type, node( SYM_STRUCT, 0, yyvsp[-3].npval, yyvsp[-1].npval ), 0 ); }
break;
case 31:
#line 259 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[0].npval, 0 ); }
break;
case 32:
#line 261 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[-1].npval, yyvsp[0].npval ); }
break;
case 33:
#line 263 "nabgrm.y"
{ yyval.npval = node( SYM_DECL, 0, yyvsp[-2].npval, yyvsp[-1].npval ); }
break;
case 34:
#line 264 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[0].npval, 0 ); }
break;
case 35:
#line 266 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 36:
#line 267 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[0].npval, 0 ); }
break;
case 37:
#line 268 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 38:
#line 269 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 39:
#line 271 "nabgrm.y"
{ yyval.npval = node( SYM_LBRACK, 0, yyvsp[-3].npval, yyvsp[-1].npval ); }
break;
case 40:
#line 272 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 41:
#line 273 "nabgrm.y"
{ yyval.npval = node( SYM_HASHED, 0, 0, 0 ); }
break;
case 42:
#line 274 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[0].npval, 0 ); }
break;
case 43:
#line 276 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 44:
#line 277 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 45:
#line 278 "nabgrm.y"
{ yyval.npval = node( SYM_DYNAMIC, 0, 0, 0 ); }
break;
case 46:
#line 281 "nabgrm.y"
{ CG_genop( NULL, SYM_SEMICOLON);
				  CG_genfend( NULL ); }
break;
case 47:
#line 284 "nabgrm.y"
{ CG_genop( NULL, SYM_SEMICOLON );
				  CG_genfend( yyvsp[-1].npval ); }
break;
case 48:
#line 287 "nabgrm.y"
{ CG_genfstart(); }
break;
case 49:
#line 288 "nabgrm.y"
{ CG_genfend( NULL ); }
break;
case 50:
#line 289 "nabgrm.y"
{ CG_genfhdr( yyvsp[-1].npval, yyvsp[0].npval ); }
break;
case 51:
#line 291 "nabgrm.y"
{ CG_genplist( yyvsp[-1].npval ); }
break;
case 52:
#line 292 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 53:
#line 293 "nabgrm.y"
{ yyval.npval = NULL; }
break;
case 54:
#line 294 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[0].npval, 0 ); }
break;
case 55:
#line 296 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 56:
#line 298 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[0].npval, NULL ); 
			  	  yyval.npval = node( SYM_DECL, 0, yyvsp[-1].npval, yyval.npval ); }
break;
case 57:
#line 300 "nabgrm.y"
{ CG_genpdecls();
				  CG_genop( NULL, SYM_LBRACE ); }
break;
case 58:
#line 302 "nabgrm.y"
{ CG_genedefs( TRUE ); }
break;
case 59:
#line 304 "nabgrm.y"
{ CG_genestmts( TRUE );
				  CG_genop( NULL, SYM_RBRACE ); }
break;
case 60:
#line 306 "nabgrm.y"
{ yyval.npval=NULL; }
break;
case 62:
#line 308 "nabgrm.y"
{ yyval.npval = NULL; }
break;
case 66:
#line 312 "nabgrm.y"
{ yyval.npval = NULL; }
break;
case 80:
#line 329 "nabgrm.y"
{ CG_genmain();
				  yyval.npval = node( SYM_ALLOCATE, 0, 0, yyvsp[-1].npval );
				  CG_genexpr( yyval.npval );
				  CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 81:
#line 334 "nabgrm.y"
{ CG_genmain();
				  yyval.npval = node( SYM_ASSERT, 0, 0, yyvsp[-1].npval );
				  CG_genassert( yyval.npval ); }
break;
case 82:
#line 338 "nabgrm.y"
{ CG_genmain();
				  CG_genrword( SYM_BREAK );
				  CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 83:
#line 341 "nabgrm.y"
{ CG_genmain();
				  CG_genop( NULL, SYM_LBRACE ); }
break;
case 84:
#line 344 "nabgrm.y"
{ CG_genop( NULL, SYM_RBRACE ); }
break;
case 85:
#line 346 "nabgrm.y"
{ CG_genmain();
				  CG_genrword( SYM_CONTINUE );
				  CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 86:
#line 350 "nabgrm.y"
{ CG_genmain();
				  yyval.npval = node( SYM_DEALLOCATE, 0, 0, yyvsp[-1].npval );
			 	  CG_genexpr(yyval.npval);
				  CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 87:
#line 355 "nabgrm.y"
{ CG_genmain();
				  yyval.npval = node( SYM_DEBUG, 0, 0, yyvsp[-1].npval );
				  CG_gendebug( yyval.npval ); }
break;
case 88:
#line 359 "nabgrm.y"
{ CG_genmain();
				  yyval.npval = node( SYM_DELETE, 0, 0, yyvsp[-1].npval );
			 	  CG_genexpr(yyval.npval);
				  CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 89:
#line 364 "nabgrm.y"
{ CG_genmain(); CG_genexpr( yyvsp[-1].npval );
				  CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 91:
#line 368 "nabgrm.y"
{ CG_genrword( SYM_ELSE ); }
break;
case 94:
#line 370 "nabgrm.y"
{ CG_genmain();
				  CG_genrword(SYM_RETURN);
				  CG_genop( NULL, SYM_LPAREN ); }
break;
case 95:
#line 373 "nabgrm.y"
{ CG_genexpr( yyvsp[0].npval ); }
break;
case 96:
#line 374 "nabgrm.y"
{ CG_genop( NULL, SYM_RPAREN );
				  CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 98:
#line 378 "nabgrm.y"
{ CG_genmain(); CG_genrword( SYM_IF ); }
break;
case 99:
#line 379 "nabgrm.y"
{ CG_genop( NULL, SYM_LPAREN ); }
break;
case 100:
#line 380 "nabgrm.y"
{ yyval.npval = node( SYM_TEST, 0, 0, yyvsp[0].npval );
				  CG_genexpr( yyval.npval ); }
break;
case 101:
#line 382 "nabgrm.y"
{ CG_genop( NULL, SYM_RPAREN ); }
break;
case 102:
#line 384 "nabgrm.y"
{ CG_genmain(); CG_genrword( SYM_FOR ); }
break;
case 103:
#line 385 "nabgrm.y"
{ CG_genop( NULL, SYM_LPAREN ); }
break;
case 104:
#line 387 "nabgrm.y"
{ CG_genop( NULL, SYM_RPAREN ); }
break;
case 107:
#line 390 "nabgrm.y"
{ yyval.npval = node( SYM_FOREACH, 0, yyvsp[-2].npval, yyvsp[0].npval );
				  CG_genexpr( yyval.npval ); }
break;
case 108:
#line 393 "nabgrm.y"
{ CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 109:
#line 395 "nabgrm.y"
{ CG_genop( NULL, SYM_SEMICOLON ); }
break;
case 111:
#line 396 "nabgrm.y"
{ CG_genexpr( yyvsp[0].npval ); }
break;
case 112:
#line 397 "nabgrm.y"
{ yyval.npval = NULL; }
break;
case 113:
#line 398 "nabgrm.y"
{ yyval.npval = node( SYM_TEST, 0, 0, yyvsp[0].npval );
				  CG_genexpr( yyval.npval ); }
break;
case 114:
#line 400 "nabgrm.y"
{ yyval.npval = NULL; }
break;
case 115:
#line 402 "nabgrm.y"
{ CG_genmain(); CG_genrword( SYM_WHILE ); }
break;
case 116:
#line 403 "nabgrm.y"
{ CG_genop( NULL, SYM_LPAREN ); }
break;
case 117:
#line 404 "nabgrm.y"
{ yyval.npval = node( SYM_TEST, 0, 0, yyvsp[0].npval );
				  CG_genexpr( yyval.npval ); }
break;
case 118:
#line 406 "nabgrm.y"
{ CG_genop( NULL, SYM_RPAREN ); }
break;
case 119:
#line 409 "nabgrm.y"
{ yyval.npval = yyvsp[-1].npval; }
break;
case 120:
#line 410 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 121:
#line 411 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[0].npval, 0 ); }
break;
case 122:
#line 413 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 123:
#line 415 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 124:
#line 417 "nabgrm.y"
{ yyval.npval = node( yyvsp[-1].ival, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 125:
#line 418 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 126:
#line 419 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 127:
#line 420 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 128:
#line 422 "nabgrm.y"
{ yyval.npval = node( SYM_LBRACK, 0, yyvsp[-3].npval, yyvsp[-1].npval ); }
break;
case 129:
#line 424 "nabgrm.y"
{ yyval.npval = node( SYM_PERIOD, 0, yyvsp[-2].npval,
				  node( SYM_ATTRIBUTE, &val, 0, 0 ) ); }
break;
case 130:
#line 426 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 131:
#line 428 "nabgrm.y"
{ yyval.npval = node( SYM_OR, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 132:
#line 429 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 133:
#line 431 "nabgrm.y"
{ yyval.npval = node( SYM_AND, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 134:
#line 432 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 135:
#line 434 "nabgrm.y"
{ yyval.npval = node( yyvsp[-1].ival, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 136:
#line 435 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 137:
#line 437 "nabgrm.y"
{ yyval.npval = node( yyvsp[-1].ival, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 138:
#line 438 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 139:
#line 440 "nabgrm.y"
{ yyval.npval = node( yyvsp[-1].ival, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 140:
#line 441 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 141:
#line 443 "nabgrm.y"
{ yyval.npval = node( SYM_UPARROW, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 142:
#line 444 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 143:
#line 445 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 144:
#line 446 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 145:
#line 447 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 146:
#line 448 "nabgrm.y"
{ yyval.npval = node( yyvsp[-1].ival, 0, NULL, yyvsp[0].npval ); }
break;
case 147:
#line 450 "nabgrm.y"
{ yyval.npval = node( SYM_CALL, 0, yyvsp[-3].npval, yyvsp[-1].npval ); }
break;
case 148:
#line 452 "nabgrm.y"
{ yyval.npval = node( SYM_LPAREN, 0, NULL, yyvsp[-1].npval ); }
break;
case 149:
#line 453 "nabgrm.y"
{ yyval.npval = node( yyvsp[-1].ival, 0, 0, yyvsp[0].npval ); }
break;
case 150:
#line 454 "nabgrm.y"
{ yyval.npval = node( yyvsp[0].ival, 0, yyvsp[-1].npval, 0 ); }
break;
case 151:
#line 455 "nabgrm.y"
{ yyval.npval = yyvsp[0].npval; }
break;
case 152:
#line 456 "nabgrm.y"
{ yyval.npval = NULL; }
break;
case 153:
#line 457 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[0].npval, 0 ); }
break;
case 154:
#line 459 "nabgrm.y"
{ yyval.npval = node( SYM_LIST, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 155:
#line 460 "nabgrm.y"
{ yyval.npval = node( SYM_PARM, 0, 0, yyvsp[0].npval ); }
break;
case 156:
#line 461 "nabgrm.y"
{ yyval.npval = node( SYM_INDEX, 0, yyvsp[0].npval, 0 ); }
break;
case 157:
#line 463 "nabgrm.y"
{ yyval.npval = node( SYM_INDEX, 0, yyvsp[-2].npval, yyvsp[0].npval ); }
break;
case 158:
#line 464 "nabgrm.y"
{ yyval.ival = SYM_ASSIGN; }
break;
case 159:
#line 466 "nabgrm.y"
{ yyval.ival = SYM_PLUS_ASSIGN; }
break;
case 160:
#line 468 "nabgrm.y"
{ yyval.ival = SYM_MINUS_ASSIGN; }
break;
case 161:
#line 470 "nabgrm.y"
{ yyval.ival = SYM_STAR_ASSIGN; }
break;
case 162:
#line 472 "nabgrm.y"
{ yyval.ival = SYM_SLASH_ASSIGN; }
break;
case 163:
#line 474 "nabgrm.y"
{ yyval.ival = SYM_MODULUS_ASSIGN; }
break;
case 164:
#line 476 "nabgrm.y"
{ yyval.ival = SYM_UPARROW_ASSIGN; }
break;
case 165:
#line 477 "nabgrm.y"
{ yyval.ival = SYM_LESS; }
break;
case 166:
#line 479 "nabgrm.y"
{ yyval.ival = SYM_LESS_EQUAL; }
break;
case 167:
#line 480 "nabgrm.y"
{ yyval.ival = SYM_EQUAL; }
break;
case 168:
#line 481 "nabgrm.y"
{ yyval.ival = SYM_NOT_EQUAL; }
break;
case 169:
#line 483 "nabgrm.y"
{ yyval.ival = SYM_GREATER_EQUAL; }
break;
case 170:
#line 484 "nabgrm.y"
{ yyval.ival = SYM_GREATER; }
break;
case 171:
#line 485 "nabgrm.y"
{ yyval.ival = SYM_MATCH; }
break;
case 172:
#line 487 "nabgrm.y"
{ yyval.ival = SYM_DONT_MATCH; }
break;
case 173:
#line 488 "nabgrm.y"
{ yyval.ival = SYM_IN; }
break;
case 174:
#line 489 "nabgrm.y"
{ yyval.ival = SYM_PLUS; }
break;
case 175:
#line 490 "nabgrm.y"
{ yyval.ival = SYM_MINUS; }
break;
case 176:
#line 491 "nabgrm.y"
{ yyval.ival = SYM_STAR; }
break;
case 177:
#line 492 "nabgrm.y"
{ yyval.ival = SYM_SLASH; }
break;
case 178:
#line 493 "nabgrm.y"
{ yyval.ival = SYM_MODULUS; }
break;
case 179:
#line 494 "nabgrm.y"
{ yyval.ival = SYM_ATSIGN; }
break;
case 180:
#line 495 "nabgrm.y"
{ yyval.ival = SYM_PLUS_PLUS; }
break;
case 181:
#line 497 "nabgrm.y"
{ yyval.ival = SYM_MINUS_MINUS; }
break;
case 182:
#line 498 "nabgrm.y"
{ yyval.ival = SYM_NEGATE; }
break;
case 183:
#line 499 "nabgrm.y"
{ yyval.ival = SYM_NOT; }
break;
case 184:
#line 501 "nabgrm.y"
{ yyval.npval = node( SYM_IDENT, &val, 0, 0 ); }
break;
case 185:
#line 502 "nabgrm.y"
{ yyval.npval = node( SYM_INT_LIT, &val, 0, 0 ); }
break;
case 186:
#line 503 "nabgrm.y"
{ yyval.npval = node( SYM_FLOAT_LIT, &val, 0, 0 ); }
break;
case 187:
#line 505 "nabgrm.y"
{ yyval.npval = node( SYM_STRING_LIT, &val, 0, 0 ); }
break;
#line 1528 "y.tab.c"
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
