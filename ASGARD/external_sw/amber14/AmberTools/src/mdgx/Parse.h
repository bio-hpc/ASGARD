#ifndef ParseHeadings
#define ParseHeadings

#include "MatrixDS.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "ParamFitDS.h"

int CountWords(char* line);

char ToUpper(char c);

double RealXpYf(char* word, int X, int Y);

int WordIsNumber(char* word);

int WordIsAtomType(char* word);

cmat ParseWords(char* line);

void RemoveWhiteSpace(char* a, int asize);

void EqualSpace(char* line);

void RemoveComments(char* line);

void NixCommaCarriage(char* line);

int AdvanceToSegment(FILE *inp, char* segname, int scan0);

int DetectNamelistEnd(char* line, char* errmsg);

int ReadNamelistLine(char* line, cmat *lwords, char* callfunc, FILE *inp);

void SeekString(cmat L, char* val, char* sname, char* salias);

void SeekStringInc(cmat L, char* val, char* sname, char* salias, int *counter);

void SeekSSR(cmat L, char* val1, char* val2, double *val3, char* sname,
             char* salias, int *counter);

void SeekS3R(cmat L, char* val1, char* val2, char* val3, double *val4,
	     char* sname, char* salias, int *counter);

void SeekStringPlusVal(cmat L, char* val1, double *val2, char* sname,
		       char* salias, int *counter);

void SeekRecord(cmat L, cmat *C, char* sname, char* salias, int *counter);

cmat SeekNString(cmat L, cmat* val, int* fspec, char* sname, char* salias);

void SeekReal(cmat L, double *val, char* sname, char* salias);

void SeekNReal(cmat L, double* val, char* sname, char* salias, int maxidx);

void SeekInt(cmat L, int *val, char* sname, char* salias);

void SeekLLInt(cmat L, long long int *val, char* sname, char* salias);

int* ParseAmbMask(char* maskstr, prmtop *tp, coord *crd);

FILE* FOpenSafe(char* fname, int ovrwrt);

void SeekTorsionID(cmat L, prmset *mp, char* sname, char* salias,
                   int *maxhadj);

void SeekRecast(cmat L, prmset *mp, char* sname, char* salias, int *maxhold,
                int specinst);

long long int ReadNumericalShorthand(char* numstr);

#endif
