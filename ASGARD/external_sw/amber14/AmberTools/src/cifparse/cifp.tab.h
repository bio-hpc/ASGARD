#define ITEMNAME 257
#define VALUE 258
#define LOOP 259
#define DATABLOCK 260
#define UNKNOWN 261
#define MISSING 262
typedef union {
	char TempBuffer[MAXVALUELENGTH+1];
} YYSTYPE;
extern YYSTYPE cifplval;
