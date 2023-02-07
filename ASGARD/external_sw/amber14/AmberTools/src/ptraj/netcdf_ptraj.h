// netcdf_ptraj.h
#ifdef BINTRAJ
#  ifdef MPI
#    include "../../include/pnetcdf.h"
#    define nc_strerror ncmpi_strerror
#  else
#    include "netcdf.h"
#  endif 
#endif

void dan_netcdf_debug(int ncid);
int NETCDF_info_debug(netcdfTrajectoryInfo *N, char *filename);
int NETCDF_open(coordinateInfo *trajInfo);
int NETCDF_close(coordinateInfo *trajInfo);
int NETCDF_setup(coordinateInfo *trajInfo, int *actualAtoms);
int NETCDF_setupOutput(coordinateInfo *trajInfo, int atoms);
#ifdef BINTRAJ
int netcdfDefineVariable( int ncid, char *name, nc_type xtype, int ndims, int dimids[], int *varidp );
void netcdfDefineDimension( int ncid, char *name, int length, int *dimidp );
char *netcdfGetAttributeText( int ncid, int vid, char *attribute );
void netcdfPutAttributeText( int ncid, int vid, char *attribute, char *text );
int netcdfGetDimensionInfo(int ncid, char *attribute, int *length);
#endif
