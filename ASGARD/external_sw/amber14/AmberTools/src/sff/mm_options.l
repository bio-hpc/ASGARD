%{
extern int mytaskid;

#ifdef flex

#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r = mmo_yyinput(b, ms))
int mmo_mmoinput( char*, int );

#else

#undef input
#undef unput
#define input() (*gopts++)
#define unput(c) (*--gopts = c)

#endif

#undef ECHO
#define ECHO yyout=stdout;if(mytaskid==0) fprintf(yyout, "\tmm_options:  %s\n",yytext)

char* setVarStr ( char *str);
int setVarInt ( char *str);
double setVarDbl ( char *str);
void setArrayStr ( char *str, char *a, int n, int len);
void setArrayInt ( char *str, int *a, int n);
void setArrayDbl ( char *str, double *a, int n);
%}
Float [+-]?((([0-9]+\.?[0-9]*)|([0-9]*\.?[0-9]+))([eE][-+]?[0-9]+)?)
Int [+-]?[0-9]+
PosInt [+]?[1-9][0-9]*
NegInt [-]?[1-9][0-9]*
NonNegInt [+]?[0-9]+
Str [^,\ =]+
Assign (\ *=\ *|\ +)
ListSep (\ *,\ *)
%%

cut{Assign}{Float}		         { ECHO; cut = setVarDbl(yytext); }
cutnp{Assign}{Float}	                 { ECHO; cutnp = setVarDbl(yytext); }
scnb{Assign}{Float}		         { ECHO; scnb = setVarDbl(yytext); }
scee{Assign}{Float}		         { ECHO; scee = setVarDbl(yytext); }
ntpr{Assign}[1-9][0-9]*		         { ECHO; ntpr = setVarInt(yytext); 
                                           if( ntpr==0) ntpr = 9999999; }
hcp{Assign}[012]			 { ECHO; hcp = setVarInt(yytext); }
hcp_h1{Assign}{Float}                    { ECHO; hcp_h1 = setVarDbl(yytext); }
hcp_h2{Assign}{Float}                    { ECHO; hcp_h2 = setVarDbl(yytext); }
hcp_h3{Assign}{Float}                    { ECHO; hcp_h3 = setVarDbl(yytext); }

nscm{Assign}{NonNegInt}                  { ECHO; nscm = setVarInt(yytext); }
nsnb{Assign}{NonNegInt}			 { ECHO; nsnb = setVarInt(yytext); }
nsnp{Assign}{NonNegInt}			 { ECHO; nsnp = setVarInt(yytext); }
ntwx{Assign}{NonNegInt}			 { ECHO; ntwx = setVarInt(yytext); }
nchk{Assign}{NonNegInt}			 { ECHO; nchk = setVarInt(yytext); }
nchk2{Assign}{NonNegInt}		 { ECHO; nchk2 = setVarInt(yytext); }

dt{Assign}{Float}		         { ECHO; dt = setVarDbl(yytext); }
t{Assign}{Float}		         { ECHO; t = setVarDbl(yytext); }
tautp{Assign}{Float}	                 { ECHO; tautp = setVarDbl(yytext); }
temp0{Assign}{Float}	                 { ECHO; temp0 = setVarDbl(yytext); }
tempi{Assign}{Float}	                 { ECHO; tempi = setVarDbl(yytext); }
vlimit{Assign}{Float}	                 { ECHO; vlimit = setVarDbl(yytext); }
gamma_ln{Assign}{Float}	                 { ECHO; gamma_ln = setVarDbl(yytext); }
genmass{Assign}{Float}	                 { ECHO; genmass = setVarDbl(yytext); }
wcons{Assign}{Float}	                 { ECHO; wcons = setVarDbl(yytext); }
ntpr_md{Assign}[1-9][0-9]*		 { ECHO; ntpr_md = setVarInt(yytext); 
                                           if( ntpr_md==0) ntpr_md = 9999999; }
zerov{Assign}{NonNegInt}		 { ECHO; zerov = setVarInt(yytext); }
k4d{Assign}{Float}		         { ECHO; k4d = setVarDbl(yytext); }

ips{Assign}[0123]			 { ECHO; ips = setVarInt(yytext); }

gb{Assign}[012578]			 { ECHO; gb = setVarInt(yytext); }
gb_debug{Assign}[01]		         { ECHO; gb_debug = setVarInt(yytext); }
nr_debug{Assign}[01]		         { ECHO; nr_debug = setVarInt(yytext); }
gb2_debug{Assign}[01]		         { ECHO; gb2_debug = setVarInt(yytext); }
gbsa_debug{Assign}[01]		         { ECHO; gbsa_debug = setVarInt(yytext); }
e_debug{Assign}[0123]		         { ECHO; e_debug = setVarInt(yytext); }
epsext{Assign}{Float}	                 { ECHO; epsext = setVarDbl(yytext); }
kappa{Assign}{Float}	                 { ECHO; kappa = setVarDbl(yytext); }
rgbmax{Assign}{Float}	                 { ECHO; rgbmax = setVarDbl(yytext); }

surften{Assign}{Float}	                 { ECHO; surften = setVarDbl(yytext); }
dradius{Assign}{Float}	                 { ECHO; dradius = setVarDbl(yytext); }
deltar{Assign}{Float}	                 { ECHO; deltar = setVarDbl(yytext); }
rwater{Assign}{Float}	                 { ECHO; rwater = setVarDbl(yytext); }
alphanp{Assign}{Float}	                 { ECHO; alphanp = setVarDbl(yytext); }
kappanp{Assign}{Float}	                 { ECHO; kappanp = setVarDbl(yytext); }
min_volume{Assign}{Float}	         { ECHO; min_volume = setVarDbl(yytext); }
max_set_size{Assign}{NonNegInt}		 { ECHO; max_set_size = setVarInt(yytext); }
MPI_min_tasks{Assign}{NonNegInt}	 { ECHO; MPI_min_tasks = setVarInt(yytext); }
cull_np_lists{Assign}[01]  	         { ECHO; cull_np_lists = setVarInt(yytext); }
use_lower_tri{Assign}[01]  	         { ECHO; use_lower_tri = setVarInt(yytext); }
dynamic_loops{Assign}[01]  	         { ECHO; dynamic_loops = setVarInt(yytext); }
static_arrays{Assign}[01]  	         { ECHO; static_arrays = setVarInt(yytext); }
blocksize{Assign}{NonNegInt}		 { ECHO; blocksize = setVarInt(yytext); }
gbsa{Assign}[01]	                 { ECHO; gbsa = setVarInt(yytext); }
dim{Assign}[34]                          { ECHO; dim = setVarInt(yytext); }

rism{Assign}[01]                         { ECHO; rismData.rism=setVarInt(yytext);}
xvvfile{Assign}{Str}                     { ECHO; xvvfile=setVarStr(yytext);}
guvfile{Assign}{Str}                     { ECHO; guvfile=setVarStr(yytext); }
huvfile{Assign}{Str}                     { ECHO; huvfile=setVarStr(yytext); }
cuvfile{Assign}{Str}                     { ECHO; cuvfile=setVarStr(yytext); }
uuvfile{Assign}{Str}                     { ECHO; uuvfile=setVarStr(yytext); }
asympfile{Assign}{Str}                   { ECHO; asympfile=setVarStr(yytext); }
quvfile{Assign}{Str}                     { ECHO; quvfile=setVarStr(yytext); }
chgdistfile{Assign}{Str}                 { ECHO; chgdistfile=setVarStr(yytext); }
exchemfile{Assign}{Str}                  { ECHO; exchemfile=setVarStr(yytext); }
solvenefile{Assign}{Str}                 { ECHO; solvenefile=setVarStr(yytext); }
entropyfile{Assign}{Str}                 { ECHO; entropyfile=setVarStr(yytext); }
exchemGFfile{Assign}{Str}                { ECHO; exchemGFfile=setVarStr(yytext); }
solveneGFfile{Assign}{Str}               { ECHO; solveneGFfile=setVarStr(yytext); }
entropyGFfile{Assign}{Str}               { ECHO; entropyGFfile=setVarStr(yytext); }
exchemUCfile{Assign}{Str}                { ECHO; exchemUCfile=setVarStr(yytext); }
solveneUCfile{Assign}{Str}               { ECHO; solveneUCfile=setVarStr(yytext); }
entropyUCfile{Assign}{Str}               { ECHO; entropyUCfile=setVarStr(yytext); }
potUVfile{Assign}{Str}                   { ECHO; potUVfile=setVarStr(yytext); }
volfmt{Assign}{Str}                      { ECHO; volfmt=setVarStr(yytext); }
 /*
closure{Assign}0                         { ECHO; closure="hnc"; }
closure{Assign}1                         { ECHO; closure="kh"; }
closure{Assign}2                         { ECHO; closure="pse"; }
 */

 /* read in an array of closure names (strings). */
 /* This is tricky because we need to avoid eating trailing keywords */
 /* in the same mm_options call. To do this we first enforce the rule */
 /* that there must be a comma before a trailing keyword. This differs */
 /* from the rest of mm_options input where either whitespace or a comma */
 /* may be used. The we have two rules:*/
 /* 1) More keywords follow the list. Reads everything in list, */
 /* including the comma, with the requirement of being followed by a */
 /* keyword. */
closure{Assign}({Str}{ListSep}){1,10}/{Str}{Assign} {
               /* find the trailing comma and replace it with a null character */
               char * comma = strrchr(yytext,',');
               *comma='\0';
               ECHO; 
               setArrayStr(yytext,closure[0],CLOSURELEN,NCLOSURE); 
}
 /* 2) Last keyword in the input line.  Should read the rest of the */
 /* line. By enforcing comma separated structure we avoid eating up */
 /* key-value pairs that are separated by white space only*/
closure{Assign}({Str}\ *,\ *){0,9}{Str} { ECHO; setArrayStr(yytext,closure[0],CLOSURELEN,NCLOSURE); }

closureOrder{Assign}{PosInt}             { ECHO; rismData.closureOrder=setVarInt(yytext); }
asympCorr{Assign}[01]                    { ECHO; rismData.asympCorr=setVarInt(yytext); }
solvcut{Assign}{Float}     	         { ECHO; rismData.solvcut=setVarDbl(yytext); }
buffer{Assign}{Float}      	         { ECHO; rismData.buffer=setVarDbl(yytext); }
grdspcx{Assign}{Float}	                 { ECHO; rismData.grdspc[0]=setVarDbl(yytext); }
grdspcy{Assign}{Float}	                 { ECHO; rismData.grdspc[1]=setVarDbl(yytext); }
grdspcz{Assign}{Float}	                 { ECHO; rismData.grdspc[2]=setVarDbl(yytext); }
grdspc{Assign}{Float}(\ *,\ *{Float}){2,2}      { ECHO; setArrayDbl(yytext,rismData.grdspc,3);}
ngx{Assign}{Int}                       	 { ECHO; rismData.ng3[0]=setVarInt(yytext); }
ngy{Assign}{Int}        	         { ECHO; rismData.ng3[1]=setVarInt(yytext); }
ngz{Assign}{Int}        	         { ECHO; rismData.ng3[2]=setVarInt(yytext); }
ng{Assign}{Int}(\ *,\ *{Int}){2,2}       { ECHO; setArrayInt(yytext,rismData.ng3,3);}
solvboxx{Assign}{Float}	                 { ECHO; rismData.solvbox[0]=setVarDbl(yytext); }
solvboxy{Assign}{Float}	                 { ECHO; rismData.solvbox[1]=setVarDbl(yytext); }
solvboxz{Assign}{Float}	                 { ECHO; rismData.solvbox[2]=setVarDbl(yytext); }
solvbox{Assign}{Float}(\ *,\ *{Float}){2,2} { ECHO; setArrayDbl(yytext,rismData.solvbox,3);}
tolerance{Assign}{Float}({ListSep}{Float})* { ECHO; setArrayDbl(yytext,rismtol,NCLOSURE); }
mdiis_del{Assign}{Float}	         { ECHO; rismData.mdiis_del=setVarDbl(yytext); }
mdiis_restart{Assign}{Float}	         { ECHO; rismData.mdiis_restart=setVarDbl(yytext); }
mdiis_nvec{Assign}{PosInt}	       	 { ECHO; rismData.mdiis_nvec=setVarInt(yytext); }
mdiis_method{Assign}[012]                { ECHO; rismData.mdiis_method=setVarInt(yytext); }
maxstep{Assign}{PosInt}	       	         { ECHO; rismData.maxstep=setVarInt(yytext); }
npropagate{Assign}[0-5]                  { ECHO; rismData.npropagate=setVarInt(yytext); }
centering{Assign}[0-4\-\+]+              { ECHO; rismData.centering=setVarInt(yytext); }
zerofrc{Assign}[01]                      { ECHO; rismData.zerofrc=setVarInt(yytext); }
apply_rism_force{Assign}[01]             { ECHO; rismData.apply_rism_force=setVarInt(yytext); }
polarDecomp{Assign}[01]                  { ECHO; rismData.polarDecomp=setVarInt(yytext); }
entropicDecomp{Assign}[01]               { ECHO; rismData.entropicDecomp=setVarInt(yytext); }
rismnrespa{Assign}{PosInt}	       	 { ECHO; rismData.rismnrespa=setVarInt(yytext); }
fcestride{Assign}{NonNegInt}	         { ECHO; rismData.fcestride=setVarInt(yytext); }
fcecut{Assign}{Float}	                 { ECHO; rismData.fcecut=setVarDbl(yytext); }
fcenbasis{Assign}{PosInt}	       	 { ECHO; rismData.fcenbasis=setVarInt(yytext); }
fcecrd{Assign}[012]                      { ECHO; rismData.fcecrd=setVarInt(yytext); }
saveprogress{Assign}{NonNegInt}	         { ECHO; rismData.saveprogress=setVarInt(yytext); }
ntwrism{Assign}{Int}	       	         { ECHO; rismData.ntwrism=setVarInt(yytext); }
verbose{Assign}[012]                     { ECHO; rismData.verbose=setVarInt(yytext); }
progress{Assign}{NonNegInt}	         { ECHO; rismData.progress=setVarInt(yytext); }
ntpr_rism{Assign}{PosInt}                { ECHO; ntpr_rism=setVarInt(yytext); }
uccoeff{Assign}{Float}({ListSep}{Float})* { ECHO; setArrayDbl(yytext,rismData.uccoeff,2); }

ipb{Assign}[012]			 { ECHO; pbsa = setVarInt(yytext); }
inp{Assign}[012]			 { ECHO; inp = setVarInt(yytext); }
epsin{Assign}{Float}	                 { ECHO; epsin = setVarDbl(yytext); }
epsout{Assign}{Float}	                 { ECHO; epsout = setVarDbl(yytext); }
smoothopt{Assign}[012]		         { ECHO; smoothopt = setVarInt(yytext); }
istrng{Assign}{Float}	                 { ECHO; istrng = setVarDbl(yytext); }
ivalence{Assign}{Float}	                 { ECHO; ivalence = setVarDbl(yytext); }
radiopt{Assign}[01]		         { ECHO; radiopt = setVarInt(yytext); }
dprob{Assign}{Float}	                 { ECHO; dprob = setVarDbl(yytext); }
iprob{Assign}{Float}	                 { ECHO; iprob = setVarDbl(yytext); }
npbopt{Assign}[01]			 { ECHO; npbopt = setVarInt(yytext); }
solvopt{Assign}[1-68]		         { ECHO; solvopt = setVarInt(yytext); }
accept{Assign}{Float}	                 { ECHO; accept = setVarDbl(yytext); }
maxitn{Assign}{NonNegInt}		 { ECHO; maxitn = setVarInt(yytext); }
fillratio{Assign}{Float}	         { ECHO; fillratio = setVarDbl(yytext); }
space{Assign}{Float}	                 { ECHO; space = setVarDbl(yytext); }
nfocus{Assign}[12]			 { ECHO; nfocus = setVarInt(yytext); }
fscale{Assign}[1-9]+		         { ECHO; fscale = setVarInt(yytext); }
bcopt{Assign}([156]|10)			 { ECHO; bcopt = setVarInt(yytext); }
eneopt{Assign}[12]			 { ECHO; eneopt = setVarInt(yytext); }
dbfopt{Assign}[12]			 { ECHO; dbfopt = setVarInt(yytext); }
frcopt{Assign}[0-3]		         { ECHO; frcopt = setVarInt(yytext); }
cutnb{Assign}{Float}	                 { ECHO; cutnb = setVarDbl(yytext); }
sprob{Assign}{Float}	                 { ECHO; sprob = setVarDbl(yytext); }
npbverb{Assign}[0-1]		         { ECHO; npbverb = setVarInt(yytext); }
arcres{Assign}{Float}	                 { ECHO; arcres = setVarDbl(yytext); }
cavity_surften{Assign}{Float}	         { ECHO; cavity_surften = setVarDbl(yytext); } 
cavity_offset{Assign}{Float}	         { ECHO; cavity_offset = setVarDbl(yytext); }
nsnba{Assign}{NonNegInt}		 { ECHO; nsnba = setVarInt(yytext); }
npbgrid{Assign}{NonNegInt}		 { ECHO; npbgrid = setVarInt(yytext); }
maxarcdot{Assign}{NonNegInt}		 { ECHO; maxarcdot = setVarInt(yytext); }

diel{Assign}Yammp			 { ECHO; dield = -5; }
diel{Assign}R94			         { ECHO; dield = -4; }
diel{Assign}RL94			 { ECHO; dield = -3; }
diel{Assign}RL			         { ECHO; dield = -2; }
diel{Assign}R			         { ECHO; dield =  0; }
diel{Assign}C			         { ECHO; dield =  1; }
dielc{Assign}{Float}	                 { ECHO; dielc = setVarDbl(yytext); }

rattle{Assign}[012]			 { ECHO; irattle = setVarInt(yytext); }

[\ \n\t,]+						;

.						{	
#ifdef flex
  if(mytaskid==0){
    char c;
    fprintf( stderr, "Error: mm_option: unknown keyword or invalid value: %s", yytext );
    for(c=input();c != ' ' && c != '=' && c != EOF;c=input())
      fprintf( stderr, "%c",c);
    fprintf( stderr,"\n");
  }
#else
  if(mytaskid==0)fprintf( stderr, "Error: mm_option: unknown keyword or invalid value: %s\n", --gopts );
#endif
  
  exit( 1 );
}
%%

int yywrap(){
	return( 1 );
}

#ifdef flex


int  mmo_yyinput( buf, max_size )

char *buf;
int  max_size;
{
	int n;

	if( max_size > &gopts[mmoinputlim] - mmoinputptr ) 
		n =  &gopts[mmoinputlim] - mmoinputptr;
	else
		n = max_size;

	if( n>0 ){
		memcpy( buf, mmoinputptr, n );
		mmoinputptr += n;
	}
	return n;
}
#endif

/**
 * Extract a single value from an option statement
 * Takes a string of the form
 * <keyword>=<val1>
 * with an unlimited amount of white space.
 *IN:
 *   str : string 
 *OUT:
 *   returns the string after the keyword
 */
char* setVarStr ( char *str){
  char *copy,*token;
  copy = strdup(str);
  /* the first token is the variable name */
  token=strtok(copy,"= ");
  /* the remaining characters should be the string */
  token = strtok(NULL," ,=");
  return token;
}

/**
 * Extract a single value from an option statement
 * Takes a string of the form
 * <keyword>=<val1>
 * with an unlimited amount of white space.
 *IN:
 *   str : string 
 *OUT:
 *   returns the double after the keyword
 */
double setVarDbl ( char *str){
  double v[1];
  setArrayDbl(str,v,1);
  return v[0];
}

/**
 * Extract a single value from an option statement
 * Takes a string of the form
 * <keyword>=<val1>
 * with an unlimited amount of white space.
 *IN:
 *   str : string 
 *OUT:
 *   returns the integer after the keyword
 */
int setVarInt ( char *str){
  int v[1];
  setArrayInt(str,v,1);
  return v[0];
}

/**
 * Extract mulitple values from an option statement
 * Takes a string of the form
 * <keyword>=<val1>,<val2>,...
 * with an unlimited amount of white space.  Up to n Values are stored
 * in array a.  A 1D character array is used and treated as a 2D array
 * that is [len][n] and must be preallocated.  Each string in the
 * array is set to empty on entry.
 *IN:
 *   str : string 
 *   a   : array of char (len,n)
 *   len : length of 
 *   n   : number of elements of a
 *SIDE EFFECTS:
 *   Values from string are stored in a
 */
void setArrayStr ( char *str, char *a, int len, int n){
  int i;
  char *copy, *token;
  copy = malloc(strlen(str)+1);
  if(copy == NULL){
    fprintf(stderr,"Error: mm_option: failed to allocate memory parsing:\n");
    fprintf(stderr,"'%s'\n",str);
    exit(1);
  }
  strcpy(copy,str);
  /* set all strings in the array to empty */
  for(i=0; i<n; i++){
     *(a+len*i)='\0';
   }

  /* the first token is the variable name */
  token = strtok(copy,"= ");
  /* the remaining tokens should all be floating point numbers */
  for(i = 0; i< n && (token = strtok(NULL," ,=")) != NULL; i++){
        fflush(stderr);
    strncpy((a+len*i), token, len);
  }
  free(copy);  
}

/**
 * Extract mulitple values from an option statement
 * Takes a string of the form
 * <keyword>=<val1>,<val2>,...
 * with an unlimited amount of white space.  Up to n Values are stored in array a.
 *IN:
 *   str : string 
 *   a   : array of doubles
 *   n   : length of a
 *SIDE EFFECTS:
 *   Values from string are stored in a
 */
void setArrayDbl ( char *str, double *a, int n){
  int i;
  char *copy, *token;

  copy = malloc(strlen(str)+1);
  if(copy == NULL){
    fprintf(stderr,"Error: mm_option: failed to allocate memory parsing:\n");
    fprintf(stderr,"'%s'\n",str);
    exit(1);
  }
  strcpy(copy,str);
  /* the first token is the variable name */
  token = strtok(copy,"= ");
  /* the remaining tokens should all be floating point numbers */
  for(i = 0; i< n && (token = strtok(NULL," ,=")) != NULL; i++){
    sscanf(token,"%lg",&(a[i]));
  }
  free(copy);  
}

/**
 * Extract mulitple values from an option statement
 * Takes a string of the form
 * <keyword>=<val1>,<val2>,...
 * with an unlimited amount of white space.  Up to n Values are stored in array a.
 *IN:
 *   str : string 
 *   a   : array of doubles
 *   n   : length of a
 *SIDE EFFECTS:
 *   Values from string are stored in a
 */
void setArrayInt ( char *str, int *a, int n){
  int i;
  char *copy, *token;
  copy = malloc(strlen(str)+1);
  if(copy == NULL){
    fprintf(stderr,"Error: mm_option: failed to allocate memory parsing:\n");
    fprintf(stderr,"'%s'\n",str);
    exit(1);
  }
  strcpy(copy,str);
  /* the first token is the variable name */
  token = strtok(copy,"= ");
  /* the remaining tokens should all be integer numbers */
  for(i = 0; i< n && (token = strtok(NULL," ,=")) != NULL; i++){
    sscanf(token,"%i",&(a[i]));
  }
  free(copy);
}
