// Process command line options

struct opt_t {
    string conf_file, struct_type, seq, outfile;
} opts;

int get_the_opts(int ac, string av[1], struct opt_t o);
