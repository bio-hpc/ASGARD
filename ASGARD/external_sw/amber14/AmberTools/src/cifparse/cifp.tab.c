#ifndef lint
static char cifpsccsid[] = "@(#)yaccpar	1.9 (Berkeley) 02/21/93";
#endif
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define cifpclearin (cifpchar=(-1))
#define cifperrok (cifperrflag=0)
#define YYRECOVERING (cifperrflag!=0)
#define YYPREFIX "cifp"
#line 46 "cifparse.y"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cifparse.h"
int cifplex();
void cifperror();

static int curItemNo, curValueNo, itemIndex;
static int  *fieldList = NULL;

#line 58 "cifparse.y"
typedef union {
	char TempBuffer[MAXVALUELENGTH+1];
} YYSTYPE;
#line 28 "y.tab.c"
#define ITEMNAME 257
#define VALUE 258
#define LOOP 259
#define DATABLOCK 260
#define UNKNOWN 261
#define MISSING 262
#define YYERRCODE 256
short cifplhs[] = {                                        -1,
    0,    0,    4,    3,    3,    6,    6,    7,    7,    8,
    8,    1,    9,    2,    2,    2,    5,
};
short cifplen[] = {                                         2,
    1,    2,    2,    0,    2,    2,    2,    2,    2,    1,
    2,    1,    1,    1,    1,    1,    1,
};
short cifpdefred[] = {                                      4,
    0,    0,   17,    2,    4,   12,   13,    0,    5,    0,
    0,    0,   14,   15,   16,    6,    9,   10,    0,    8,
   11,
};
short cifpdgoto[] = {                                       1,
    8,   16,    2,    4,    5,    9,   10,   19,   11,
};
short cifpsindex[] = {                                      0,
 -252, -240,    0,    0,    0,    0,    0, -246,    0, -251,
 -244, -240,    0,    0,    0,    0,    0,    0, -246,    0,
    0,
};
short cifprindex[] = {                                      0,
    0,    2,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    3,    0,    0,    0,    0,    0,    0,    1,    0,
    0,
};
short cifpgindex[] = {                                      0,
   -6,  -10,    9,    0,    0,    0,    0,    0,    0,
};
#define YYTABLESIZE 263
short cifptable[] = {                                      18,
    7,    1,    3,   17,   20,    6,   13,    3,   21,   14,
   15,   13,    6,   12,   14,   15,    6,    0,    7,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    7,    0,    7,
    7,    1,    3,
};
short cifpcheck[] = {                                      10,
    0,    0,    0,   10,   11,  257,  258,  260,   19,  261,
  262,  258,  257,    5,  261,  262,  257,   -1,  259,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,  257,   -1,  259,
  260,  260,  260,
};
#define YYFINAL 1
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 262
#if YYDEBUG
char *cifpname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"ITEMNAME","VALUE","LOOP",
"DATABLOCK","UNKNOWN","MISSING",
};
char *cifprule[] = {
"$accept : Datablocks",
"Datablocks : Lines",
"Datablocks : Datablocks Datablock",
"Datablock : DatablockName Lines",
"Lines :",
"Lines : Lines Line",
"Line : ItemName Value",
"Line : ItemNameList ValueList",
"ItemNameList : Loop ItemName",
"ItemNameList : ItemNameList ItemName",
"ValueList : Value",
"ValueList : ValueList Value",
"ItemName : ITEMNAME",
"Loop : LOOP",
"Value : VALUE",
"Value : UNKNOWN",
"Value : MISSING",
"DatablockName : DATABLOCK",
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
int cifpdebug;
int cifpnerrs;
int cifperrflag;
int cifpchar;
short *cifpssp;
YYSTYPE *cifpvsp;
YYSTYPE cifpval;
YYSTYPE cifplval;
short cifpss[YYSTACKSIZE];
YYSTYPE cifpvs[YYSTACKSIZE];
#define cifpstacksize YYSTACKSIZE
#line 155 "cifparse.y"

void cifperror(s)
char *s;
/* ---------------------------------------------------------------------------
 * Purpose:  cifperror()
 * 
             Report the location of a parsing error...
 -----------------------------------------------------------------------------*/
 
{
  fprintf( stderr,"%s near line %d\n", s, lineNo);
}

void ndb_cif_process_item_name_list()
/* ----------------------------------------------------------------------
   Purpose: ndb_cif_process_item_name_list()

            Registers the item keyword for the the current item in the 
	    current category.  Maintains an index array of "valid" keyword 
	    names in fieldList[].  This array is used as an indirectly 
	    reference between keywords and values ...  
 * ----------------------------------------------------------------------*/
{
  int prevCategoryId, categoryId;
  char itemKeyword[MxNameLen], categoryName[MxNameLen];

/*
  fprintf( stderr,"Processing item name list at line %d keyword %s\n", lineNo, TempKeyword);
*/
  prevCategoryId = ndb_cif_current_category();
  categoryId = ndb_cif_get_category_name_from_item_name(categoryName, TempKeyword);


  fieldList = (int *) realloc(fieldList, (curItemNo+1) * sizeof(int));
  if (categoryId != 0 && categoryId == prevCategoryId) {
    ndb_cif_get_item_keyword_from_item_name(itemKeyword, TempKeyword);
    ndb_cif_put_item_keyword(itemKeyword);
    itemIndex++;
    fieldList[curItemNo] = itemIndex;
  }
  else {
    fprintf( stderr,"Syntax error line %d with %s\n", lineNo, TempKeyword);
    fieldList[curItemNo] = -1;
  }
  curItemNo++;

}

void ndb_cif_process_value_list()
/* ----------------------------------------------------------------------
     Purpose:  ndb_cif_process_value_list()

               Add the current value to the appropriate column in the 
               the current row.  Start a new row if necessary.
 * ----------------------------------------------------------------------*/
{
  int fieldNo;

/*
  fprintf( stderr,"Processing value list at line %d value %s\n", lineNo, TempValue);
*/
  if (fieldList[curValueNo] != -1) {
    fieldNo = fieldList[curValueNo];
    if (fieldNo == 1) ndb_cif_new_row();
    if (fieldNo != 0) ndb_cif_put_item_value(fieldList[curValueNo], TempValue);
  }
  curValueNo++;

  if (curValueNo == curItemNo) curValueNo = 0;
}

void ndb_cif_process_item_name_value_pair()
/* ----------------------------------------------------------------------
      Purpose: ndb_cif_process_item_name_value_pair()

               Assign the current value to its associated item name.
 * ----------------------------------------------------------------------*/
{
  int categoryId, tmpCategoryId;
  char categoryName[MxNameLen], itemKeyword[MxNameLen];

/*
  fprintf( stderr,"Processing item name value pair at line %d value %s\n", lineNo, TempKeyword);
*/

  tmpCategoryId = ndb_cif_current_category();
  categoryId = ndb_cif_get_category_name_from_item_name(categoryName,TempKeyword);
  curItemNo  = 1;
  curValueNo = 0;
  itemIndex  = 1;
  if (categoryId == 0) {
    if (strcmp(categoryName,"")) {
      ndb_cif_new_category(categoryName);
    }
    else {
      fprintf( stderr, "Missing category name component in %s at line %d\n", 
	     TempKeyword, lineNo);
      return ;
    }
  }
  else if (tmpCategoryId != categoryId) {
    fprintf( stderr,"Category conflict between %s and %s at line %d\n",  
		    cifFiles.datablocks[0].categories[tmpCategoryId].categoryName,
		    categoryName,
		    lineNo);
    return;
  }

  ndb_cif_get_item_keyword_from_item_name(itemKeyword, TempKeyword);
  ndb_cif_put_item_keyword(itemKeyword);
  fieldList = (int *) realloc(fieldList, sizeof(int));
  fieldList[0] = ndb_cif_current_col();
  ndb_cif_process_value_list();
}


void ndb_cif_process_loop_declaration()
/* ----------------------------------------------------------------------
     Purpose: ndb_cif_process_loop_declaration()

              Handles initialization for a new loop, by creating a new 
              category and adding the current item name to this category.
 * ---------------------------------------------------------------------- */
{
  char  categoryName[MxNameLen];
  int   categoryId;

/*
  fprintf( stderr,"Processing loop declaration at line %d value %s\n", lineNo, TempKeyword);
*/

  categoryId = ndb_cif_get_category_name_from_item_name(categoryName, TempKeyword);
  if (categoryId == 0) {
    ndb_cif_new_category(categoryName);
    ndb_cif_process_item_name_list();
  }
  else {
    fprintf( stderr,"Duplicate category name %s at line %d\n",categoryName, lineNo);
  }
}

#if 0
void cifpprint(FILE *file, int type, YYSTYPE value)
{
  fprintf(file,"\nType  = %d\n",type);
  fprintf(file,"Value = %s\n",value.TempBuffer);
}
#endif

#line 333 "y.tab.c"
#define YYABORT goto cifpabort
#define YYREJECT goto cifpabort
#define YYACCEPT goto cifpaccept
#define YYERROR goto cifperrlab
int
cifpparse()
{
    register int cifpm, cifpn, cifpstate;
#if YYDEBUG
    register char *cifps;
    extern char *getenv();

    if (cifps = getenv("YYDEBUG"))
    {
        cifpn = *cifps;
        if (cifpn >= '0' && cifpn <= '9')
            cifpdebug = cifpn - '0';
    }
#endif

    cifpnerrs = 0;
    cifperrflag = 0;
    cifpchar = (-1);

    cifpssp = cifpss;
    cifpvsp = cifpvs;
    *cifpssp = cifpstate = 0;

cifploop:
    if ((cifpn = cifpdefred[cifpstate]) != 0) goto cifpreduce;
    if (cifpchar < 0)
    {
        if ((cifpchar = cifplex()) < 0) cifpchar = 0;
#if YYDEBUG
        if (cifpdebug)
        {
            cifps = 0;
            if (cifpchar <= YYMAXTOKEN) cifps = cifpname[cifpchar];
            if (!cifps) cifps = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, cifpstate, cifpchar, cifps);
        }
#endif
    }
    if ((cifpn = cifpsindex[cifpstate]) != 0 && (cifpn += cifpchar) >= 0 &&
            cifpn <= YYTABLESIZE && cifpcheck[cifpn] == cifpchar)
    {
#if YYDEBUG
        if (cifpdebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, cifpstate, cifptable[cifpn]);
#endif
        if (cifpssp >= cifpss + cifpstacksize - 1)
        {
            goto cifpoverflow;
        }
        *++cifpssp = cifpstate = cifptable[cifpn];
        *++cifpvsp = cifplval;
        cifpchar = (-1);
        if (cifperrflag > 0)  --cifperrflag;
        goto cifploop;
    }
    if ((cifpn = cifprindex[cifpstate]) != 0 && (cifpn += cifpchar) >= 0 &&
            cifpn <= YYTABLESIZE && cifpcheck[cifpn] == cifpchar)
    {
        cifpn = cifptable[cifpn];
        goto cifpreduce;
    }
    if (cifperrflag) goto cifpinrecovery;
#ifdef lint
    goto cifpnewerror;
#endif
cifpnewerror:
    cifperror("syntax error");
#ifdef lint
    goto cifperrlab;
#endif
cifperrlab:
    ++cifpnerrs;
cifpinrecovery:
    if (cifperrflag < 3)
    {
        cifperrflag = 3;
        for (;;)
        {
            if ((cifpn = cifpsindex[*cifpssp]) != 0 && (cifpn += YYERRCODE) >= 0 &&
                    cifpn <= YYTABLESIZE && cifpcheck[cifpn] == YYERRCODE)
            {
#if YYDEBUG
                if (cifpdebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *cifpssp, cifptable[cifpn]);
#endif
                if (cifpssp >= cifpss + cifpstacksize - 1)
                {
                    goto cifpoverflow;
                }
                *++cifpssp = cifpstate = cifptable[cifpn];
                *++cifpvsp = cifplval;
                goto cifploop;
            }
            else
            {
#if YYDEBUG
                if (cifpdebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *cifpssp);
#endif
                if (cifpssp <= cifpss) goto cifpabort;
                --cifpssp;
                --cifpvsp;
            }
        }
    }
    else
    {
        if (cifpchar == 0) goto cifpabort;
#if YYDEBUG
        if (cifpdebug)
        {
            cifps = 0;
            if (cifpchar <= YYMAXTOKEN) cifps = cifpname[cifpchar];
            if (!cifps) cifps = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, cifpstate, cifpchar, cifps);
        }
#endif
        cifpchar = (-1);
        goto cifploop;
    }
cifpreduce:
#if YYDEBUG
    if (cifpdebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, cifpstate, cifpn, cifprule[cifpn]);
#endif
    cifpm = cifplen[cifpn];
    cifpval = cifpvsp[1-cifpm];
    switch (cifpn)
    {
case 6:
#line 81 "cifparse.y"
{
  ndb_cif_process_item_name_value_pair();
}
break;
case 8:
#line 91 "cifparse.y"
{
  ndb_cif_process_loop_declaration();
}
break;
case 9:
#line 96 "cifparse.y"
{
  ndb_cif_process_item_name_list();
}
break;
case 10:
#line 103 "cifparse.y"
{
  ndb_cif_rewind_column();
  ndb_cif_process_value_list();
}
break;
case 11:
#line 108 "cifparse.y"
{
  ndb_cif_process_value_list();
}
break;
case 12:
#line 114 "cifparse.y"
{
  /*  sprintf(TempKeyword,"%s",$1); */
  strncpy(TempKeyword,cifpvsp[0].TempBuffer,MxNameLen);
}
break;
case 13:
#line 121 "cifparse.y"
{
  curItemNo = 0;  curValueNo = 0;  itemIndex = 0;
}
break;
case 14:
#line 127 "cifparse.y"
{
  /*  sprintf(TempValue,"%s",$1); */
  strncpy(TempValue,cifpvsp[0].TempBuffer,MAXVALUELENGTH);
}
break;
case 15:
#line 132 "cifparse.y"
{
  strcpy(TempValue,"");
}
break;
case 16:
#line 136 "cifparse.y"
{
  strcpy(TempValue,"");
}
break;
case 17:
#line 142 "cifparse.y"
{
  int idatablock;
  idatablock = ndb_cif_get_datablock_id(&(cifpvsp[0].TempBuffer)[5]);
  if (idatablock == 0)
    ndb_cif_new_datablock(&(cifpvsp[0].TempBuffer)[5]);
  else {
    ndb_cif_move_datablock(&(cifpvsp[0].TempBuffer)[5]);
    ndb_cif_reset_datablock();
  }
}
break;
#line 550 "y.tab.c"
    }
    cifpssp -= cifpm;
    cifpstate = *cifpssp;
    cifpvsp -= cifpm;
    cifpm = cifplhs[cifpn];
    if (cifpstate == 0 && cifpm == 0)
    {
#if YYDEBUG
        if (cifpdebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        cifpstate = YYFINAL;
        *++cifpssp = YYFINAL;
        *++cifpvsp = cifpval;
        if (cifpchar < 0)
        {
            if ((cifpchar = cifplex()) < 0) cifpchar = 0;
#if YYDEBUG
            if (cifpdebug)
            {
                cifps = 0;
                if (cifpchar <= YYMAXTOKEN) cifps = cifpname[cifpchar];
                if (!cifps) cifps = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, cifpchar, cifps);
            }
#endif
        }
        if (cifpchar == 0) goto cifpaccept;
        goto cifploop;
    }
    if ((cifpn = cifpgindex[cifpm]) && (cifpn += cifpstate) >= 0 &&
            cifpn <= YYTABLESIZE && cifpcheck[cifpn] == cifpstate)
        cifpstate = cifptable[cifpn];
    else
        cifpstate = cifpdgoto[cifpm];
#if YYDEBUG
    if (cifpdebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *cifpssp, cifpstate);
#endif
    if (cifpssp >= cifpss + cifpstacksize - 1)
    {
        goto cifpoverflow;
    }
    *++cifpssp = cifpstate;
    *++cifpvsp = cifpval;
    goto cifploop;
cifpoverflow:
    cifperror("yacc stack overflow");
cifpabort:
    return (1);
cifpaccept:
    return (0);
}
