// This from $AMBERHOME/src/nab/stringutil.c
#define MAXSTRINGLENGTH 1069

struct opt_t {
    STRING_T *conf_file,  *struct_type,  *seq,  *outfile;
};

INT_T get_the_opts(INT_T *ac, STRING_T * const *av, struct opt_t *o);
