// Wrappers for MPI file routines
#include "ptraj.h"

/* ================= MPI WRAPPERS ========================== */
#ifdef MPI
#include "mpi.h"
/*
 * printMPIerr()
 * Wrapper for MPI_Error string.
 */
void printMPIerr(int err, char *actionName) {
  int len,eclass,i;
  char buffer[1024];

  MPI_Error_string(err,buffer,&len);
  MPI_Error_class(err,&eclass);
  // Remove newlines from MPI error string
  for (i=0; i<len; i++)
    if (buffer[i]=='\n') buffer[i]=':';
  fprintf(stdout,"[%i] MPI ERROR %d: %s: [%s]\n",worldrank,eclass,actionName,buffer);

  return;
}

/*
 * parallel_check_error()
 * All ranks pass in error value. Compute the sum. return non-zero if error on
 * any nodes.
 */
int parallel_check_error(int err) {
  int errtotal;

  errtotal=0;
  MPI_Allreduce(&err,&errtotal,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  return errtotal;
}
#endif

/*
 * parallel_open_file_read()
 * Use MPI to open a coordinate file for reading.
 * Return MPI file handle.
 */
int parallel_open_file_read(coordinateInfo *C, char *filename) {
#ifdef MPI
  int err;
  MPI_File *mfp;

  C->mfp=NULL;
  if (prnlev>0) fprintf(stdout,"[%i] parallel_open_file_read: Opening input file %s\n",worldrank,filename);
  mfp = (MPI_File*) malloc(sizeof(MPI_File));
  err=MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, mfp);
  if (err!=MPI_SUCCESS) printMPIerr(err,"parallel_open_file_read()");
  // Check that everyone opened the file. 
  if (parallel_check_error(err)!=0) {
    free(mfp); 
    return 1;
  } else {
    C->mfp=mfp;
    return 0;
  }
#endif
  return 1;
}

/*
 * parallel_open_file_write()
 * Use MPI to open a coordinate file for writing, delete if present.
 * Return MPI file handle.
 */
int parallel_open_file_write(coordinateInfo *C, char *filename) {
#ifdef MPI
  int err,errtotal;
  MPI_File *mfp;
  
  C->mfp=NULL;
  // Remove file if present
  MPI_File_delete(filename,MPI_INFO_NULL);

  if (prnlev>0) fprintf(stdout,"[%i] parallel_open_file_write: Opening output file %s\n",worldrank,filename);
  mfp = (MPI_File*) malloc(sizeof(MPI_File));
  err=MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, mfp);
  if (err!=MPI_SUCCESS) printMPIerr(err,"parallel_open_file_write()");
  // Check that everyone opened the file. 
  if (parallel_check_error(err)!=0) {
    free(mfp);  
    return 1;
  } else {
    C->mfp=mfp;
    return 0;
  }
#endif
  return 1;
}

/*
 * parallel_close_file()
 * Close MPI file.
 */
int parallel_close_file(coordinateInfo *C) {
#ifdef MPI
  int err;
  if (C->mfp==NULL) return 0;
  //MPI_Barrier(MPI_COMM_WORLD);
  if (prnlev>0) fprintf(stdout,"[%i] parallel_close_file: Closing file.\n",worldrank);
  err=MPI_File_close((MPI_File *) C->mfp);
  if (err!=MPI_SUCCESS) printMPIerr(err,"parallel_close_file");
  // DAN ROE - Free the memory used by the file pointer
  free(C->mfp);
  C->mfp=NULL;
#endif 
 return 0;
}

/*
 * parallel_fgets()
 * Like fgets, use mpi file routines to get all chars up to and including 
 * null or newline. Returns buffer, or NULL on error. 
 */
char *parallel_fgets(char *buffer, int num, coordinateInfo *C) {
#ifdef MPI
  int i,err;

  for (i=0; i<num-1; i++) {
    err=MPI_File_read(*(C->mfp),buffer+i,1,MPI_CHAR,MPI_STATUS_IGNORE);
    if (err!=MPI_SUCCESS) {
      printMPIerr(err,"parallel_fgets");
      return NULL;
    }
    if (buffer[i]=='\n' || buffer[i]=='\0') {i++; break;} // Always have i be 1 more char ahead   
  }

  if (i==num && buffer[i-1]=='\n')
    buffer[i-1]='\0';
  else
    buffer[i]='\0';

  return buffer;
#endif
  return NULL;
}

/*
 * parallel_fseek()
 */
int parallel_fseek(coordinateInfo *C, int frame) {
#ifdef MPI
  int err;

  err=MPI_File_seek( *(C->mfp), C->titleSize+(frame*C->frameSize), MPI_SEEK_SET);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err,"trajFile_fseek");
    return 1;
  } else 
    return 0;
#endif
  return 1;
}

/*
 * parallel_rewind()
 */
int parallel_rewind(coordinateInfo *C) {
#ifdef MPI
  int err;

  err=MPI_File_seek( *(C->mfp), 0L, MPI_SEEK_SET);
  return err;
#endif
  return 1;
}

/* parallel_fseek_end()
 */
int parallel_fseek_end(coordinateInfo *C) {
#ifdef MPI
  int err;

  err=MPI_File_seek( *(C->mfp),0L,MPI_SEEK_END);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err,"trajFile_fseek_end:");
    return 1;
  } else
    return 0;
#endif
  return 1;
}

/* parallel_get_position()
 * Wrapper for get position/ftell 
 * NOTE: Careful about long conversions. MPI_Offset is long long I think 
 */
int parallel_get_position(coordinateInfo *C, long int *offset) {
# ifdef MPI
  int err;
  MPI_Offset moffset;

  err=MPI_File_get_position(*(C->mfp), &moffset);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err,"trajFile_get_position");
    *offset=-1;
    return -1;
  }
  *offset=(long int) moffset;

  return 0;
#endif
  return 1;
}


/* parallel_fread()    
 * Wrapper for fread using MPI routines. Put a frame into the frame buffer.
 * To be consistent with older ptraj routines return 0 on error, 1 on success. 
 */
int parallel_fread(coordinateInfo *C) {
#ifdef MPI
  int err;
  MPI_Status status;

  err=MPI_File_read( *(C->mfp),C->buffer,C->frameSize,MPI_CHAR,&status);
  if (err!=MPI_SUCCESS) {
    printMPIerr(err,"trajFile_fread");
    return 0;
  }

  return 1;
#endif
  return 0;
}
    
