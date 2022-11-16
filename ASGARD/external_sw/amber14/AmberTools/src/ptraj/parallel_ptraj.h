
void printMPIerr(int err, char *actionName);
int parallel_open_file_read(coordinateInfo *C, char *filename);
int parallel_open_file_write(coordinateInfo *C, char *filename);
int parallel_close_file(coordinateInfo *C);
char *parallel_fgets(char *buffer, int num, coordinateInfo *C);
int parallel_fseek(coordinateInfo *C, int frame);
int parallel_rewind(coordinateInfo *C);
int parallel_fseek_end(coordinateInfo *C);
int parallel_get_position(coordinateInfo *C, long int *offset);
int parallel_fread(coordinateInfo *C);

