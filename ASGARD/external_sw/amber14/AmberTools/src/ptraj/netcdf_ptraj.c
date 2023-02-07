// NETCDF file routines
#include "ptraj.h"

// Defines for converting serial netcdf routines to parallel netcdf routines.
#ifdef MPI
#define nc_inq ncmpi_inq
#define nc_inq_varname ncmpi_inq_varname
#define nc_close ncmpi_close
#define nc_inq_varid ncmpi_inq_varid
#define nc_inq_varnatts ncmpi_inq_varnatts
#define nc_inq_attlen ncmpi_inq_attlen
#define nc_get_att_text ncmpi_get_att_text
#define nc_get_att_double ncmpi_get_att_double
#define nc_set_fill ncmpi_set_fill
#define nc_enddef ncmpi_enddef
#define nc_put_vara_text ncmpi_put_vara_text
#define nc_def_var ncmpi_def_var
#define nc_def_dim ncmpi_def_dim
#define nc_put_att_text ncmpi_put_att_text
#define nc_inq_dimid ncmpi_inq_dimid
#define nc_inq_dimlen ncmpi_inq_dimlen
#endif

/* DAN ROE 
 * dan_netcdf_debug()
 * For use in printing various attributes of a previously opened netcdf file.
 */
void dan_netcdf_debug(int ncid) {
#ifdef BINTRAJ
  int ndimsp, nvarsp, ngattsp,unlimdimidp;
  int err,i;
  char *name;

  /* ncid:    NetCDF ID, from a previous call to nc open or nc create.
   * ndimsp:  Pointer to location for returned number of dimensions defined for 
   *         this netCDF dataset.
   * nvarsp:  Pointer to location for returned number of variables defined for 
   *         this netCDF dataset.
   * ngattsp: Pointer to location for returned number of global attributes 
   *         defined for this netCDF dataset.
   * unlimdimidp: 
   *  Pointer to location for returned ID of the unlimited dimension, if 
   *  there is one for this netCDF dataset. If no unlimited length 
   *  dimension has been defined, -1 is returned.
   */
  name=(char*) safe_malloc(1024*sizeof(char));
  fprintf(stdout,"========== BEG. NETCDF DEBUG ==========\n");
  err=nc_inq(ncid,&ndimsp,&nvarsp,&ngattsp,&unlimdimidp);
  fprintf(stdout,"nc_inq returned %i\n",err);
  if (err==NC_NOERR)
    fprintf(stdout,"ndimsp=%i  nvarsp=%i  ngattsp=%i  unlimdimidp=%i\n",
            ndimsp,nvarsp,ngattsp,unlimdimidp);
  else
    fprintf(stdout,"NETCDF Error occurred.\n");
  /* Print name of each variable defined in netcdf file */ 
  fprintf(stdout,"NC VARIABLES:\n");
  for (i=0; i<nvarsp; i++) {
    err=nc_inq_varname(ncid,i,name);
    fprintf(stdout,"  Var %i - ",i);
    if (err==NC_NOERR)
      fprintf(stdout,"%s\n",name);
    else
      fprintf(stdout,"NETCDF Error occured.\n");
  }  
    
  fprintf(stdout,"==========  END NETCDF DEBUG ==========\n");

  safe_free(name);
#endif
  return;
}

/* NETCDF_info_debug()
 * Print info about the NCInfo structure 
 */
int NETCDF_info_debug(netcdfTrajectoryInfo *N, char *filename) {

  fprintf(stdout,"*********** NETCDF INFO: %s ***********\n",filename);
  if (N==NULL) {
    fprintf(stdout,"  Structure is NULL.\n\n");
    return 0;
  }

  fprintf(stdout,"           ncid %i     currentFrame %i\n",
          N->ncid,N->currentFrame);
  fprintf(stdout,"       frameDID %i       spatialDID %i          atomDID %i\n",
          N->frameDID,N->spatialDID,N->atomDID);
  fprintf(stdout,"    velocityVID %i  cell_spatialDID %i  cell_angularDID %i\n",
          N->velocityVID,N->cell_spatialDID,N->cell_angularDID);
  fprintf(stdout,"       labelDID %i       spatialVID %i          timeVID %i\n",
          N->labelDID,N->spatialVID,N->timeVID);
  fprintf(stdout,"  coordinateVID %i  cell_spatialVID %i  cell_angularVID %i\n",
          N->coordinateVID,N->cell_spatialVID,N->cell_angularVID);
  fprintf(stdout,"      TempVarID %i    cellLengthVID %i     cellAngleVID %i\n",
          N->TempVarID,N->cellLengthVID,N->cellAngleVID);
/*  double velocityScale;
  char *timeUnits;
  char *coordinateUnits;
  char *cellLengthUnits;
  char *cellAngleUnits;
  char *velocityUnits;
  char *Conventions;
  char *ConventionVersion;*/
  if (N->R!=NULL) 
    fprintf(stdout,"  COORDINATE ARRAY IS NULL\n");
  else
    fprintf(stdout,"  COORDINATE ARRAY HAS BEEN ALLOCATED.\n");
  fprintf(stdout,"\n");
  return 0;
}


/*
 * NETCDF_open()
 * Open the trajectory specified by the filename and accessMode in trajInfo as 
 * a NETCDF traj.
 * Return 0 on success, 1 on failure
 */
int NETCDF_open(coordinateInfo *trajInfo) {
#ifdef BINTRAJ
  int err,ncid;

// NOTE: Put in a check, only open if coord is unknown or netcdf
  if (prnlev>0) fprintf(stdout,"[%i] NETCDF_open(): Opening %s\n",
                        worldrank,trajInfo->filename);

  switch (trajInfo->accessMode) {
    case 0: // Read 
#     ifdef MPI
      err = ncmpi_open(MPI_COMM_WORLD, trajInfo->filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
      /* This next line is a test. Apparently it puts the netcdf file in an
       * independent I/O mode. Not sure if it is bad to always put here. 
       * Originally this call was only made from ptrajPreprocess...
       */
      if (err == NC_NOERR)
        err = ncmpi_begin_indep_data(ncid);
#     else
      err = nc_open(trajInfo->filename, NC_NOWRITE, &ncid);
#     endif
    break;
    case 1: // Write
      //omode=NC_WRITE; 
#     ifdef MPI
      err = ncmpi_create(MPI_COMM_WORLD, trajInfo->filename, NC_64BIT_OFFSET, MPI_INFO_NULL, &ncid);
      if (err == NC_NOERR)
        ncmpi_begin_indep_data(ncid);
#     else
      err = nc_create(trajInfo->filename, NC_64BIT_OFFSET, &ncid);
#     endif
    break;
    case 2: // Append
      printfone("Appending of NETCDF files is not supported.\n");
      return 1;
      break;
  }

  /* If opening succeeded and memory hasnt been allocated already
   *  initialize necessary data structure.
   * NOTE: Should this be in NETCDF_setup? If so ncid would have to 
   *       be its own variable in coordinateInfo.
   * NOTE: If this is an output file trajInfo->type has already been set.
   *       Not a huge problem but is a bit circular. TRAJOUT should eventually
   *       only set trajInfo->isNetcdf.
   */
  if (err == NC_NOERR) {
    trajInfo->type = COORD_AMBER_NETCDF;
    if (trajInfo->NCInfo==NULL) {
      trajInfo->NCInfo = (netcdfTrajectoryInfo *) safe_malloc(sizeof(netcdfTrajectoryInfo));
      INITIALIZE_netcdfTrajectoryInfo( trajInfo->NCInfo );
    }
    trajInfo->NCInfo->currentFrame = worldrank;
    // Always set NCID since it can change depending on when file is opened
    trajInfo->NCInfo->ncid = ncid;
    if (prnlev>0) fprintf(stdout,"NETCDF_open(): %s has been assigned ncid of %i\n",
                          trajInfo->filename,ncid);
    return 0;
  }

  // If we are here an error occured. Print the error message before exiting.
  fprintf(stdout,"Error: NETCDF_open(): Could not open %s with accessMode %i\n",
          trajInfo->filename,trajInfo->accessMode);
  fprintf(stdout,"%s\n",nc_strerror(err));
#endif
  // If no BINTRAJ always fail
  return 1;
}

/*
 * NETCDF_close()
 * Close the trajectory specified by the NCInfo structure in trajInfo.
 * Return 0 on success, 1 on failure.
 */
int NETCDF_close(coordinateInfo *trajInfo) {
#ifdef BINTRAJ
  int err;

  if (prnlev>0) fprintf(stdout,"NETCDF_close(): Closing %s\n",trajInfo->filename);
  err = nc_close(trajInfo->NCInfo->ncid);
  if (err==NC_NOERR) return 0;

  printfone("Error closing NetCDF file %s (nc_close), error: %s\n",
            trajInfo->filename, nc_strerror(err));
#endif
  return 1;
}

/*
 * NETCDF_setup()
 * Given a trajinfo structure for which an NCinfo struct has successfully been
 * allocated, setup the netcdf file for reading. Put the number of atoms in the 
 * NETCDF file in actualAtoms.
 */
#undef  ROUTINE
#define ROUTINE "NETCDF_setup()"
int NETCDF_setup(coordinateInfo *trajInfo, int *actualAtoms) {
#ifdef BINTRAJ
  int err,spatial,i;
  netcdfTrajectoryInfo *NCInfo;
  char *filename;
#  ifdef MPI
  MPI_Offset ist;
#  else
  size_t ist;
#  endif

  if (prnlev>0) fprintf(stdout,"NETCDF_setup(): Setting up %s\n",trajInfo->filename);

  NCInfo=trajInfo->NCInfo;
  if (NCInfo == NULL) return 1;

  filename=trajInfo->filename;

  /* DAN ROE: For DEBUG of netcdf files, prints what vars are stored */
  if (prnlev>3) dan_netcdf_debug(NCInfo->ncid);

  /*
   *  Get global attributes (after initializing a structure to hold the data), 
   */
  trajInfo->title =           netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "title");
  trajInfo->application =     netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "application");
  trajInfo->program =         netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "program");
  trajInfo->version =         netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "programVersion");
  NCInfo->Conventions =       netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "Conventions");
  NCInfo->ConventionVersion = netcdfGetAttributeText(NCInfo->ncid, NC_GLOBAL, "ConventionVersion");
  if (strstr(NCInfo->Conventions, "AMBER") == NULL) {
    printfone("WARNING: NetCDF file has Conventions that do not include the string \"AMBER\"\n");
  }
  if (strcmp(NCInfo->ConventionVersion, "1.0") != 0) {
    printfone("WARNING: NetCDF file has ConventionVersion differing from \"1.0\"\n");
  }

  /*
   *  get the NetCDF dimension ID's and sizes for the frames, spatial and atoms
   */
    
  NCInfo->frameDID   = netcdfGetDimensionInfo(NCInfo->ncid, AMBER_NETCDF_FRAME, &(trajInfo->stop));
  NCInfo->spatialDID = netcdfGetDimensionInfo(NCInfo->ncid, AMBER_NETCDF_SPATIAL, &spatial);
  NCInfo->atomDID    = netcdfGetDimensionInfo(NCInfo->ncid, AMBER_NETCDF_ATOM, actualAtoms);

  if (spatial != 3) {
    printfone("ptraj cannot handle NetCDF files with other than 3 dims\n");
    return 1;
  }

  /*
   *  perform a sanity check on time variable and units.
   * DAN ROE: Should errors here really terminate ptraj?
   */
  err =    nc_inq_varid(NCInfo->ncid, AMBER_NETCDF_TIME, &NCInfo->timeVID);
  if (err != NC_NOERR) {
    printfone("Error: NetCDF time variable: %s", nc_strerror(err));
    return 1;
  }
  //  error(ROUTINE, "NetCDF time variable, error: %s", nc_strerror(err));

  err =    nc_inq_varnatts(NCInfo->ncid, NCInfo->timeVID, &i);
  if ( err != NC_NOERR ) {
    printfone("Error: Getting number of time attributes in NetCDF file %s\n", filename);
    return 1;
  }
  //  error(ROUTINE, "Getting number of time attributes in NetCDF file %s\n", filename);
  if ( i != 1 ) {
    printfone("Error: Only one time attribute is expected in NetCDF file %s\n", filename);
    return 1;
  }
  //  error(ROUTINE, "Only one time attribute is expected in NetCDF file %s\n", filename);

  err =    nc_inq_attlen(NCInfo->ncid, NCInfo->timeVID, "units", &ist);
  if (err != NC_NOERR) {
    printfone("Error: Grabbing time units attribute length in NetCDF file: %s\n", 
              nc_strerror(err));
    return 1;
  }
  NCInfo->timeUnits = (char *) safe_malloc(sizeof(char) * (ist+1));

  err =    nc_get_att_text(NCInfo->ncid, NCInfo->timeVID, "units", NCInfo->timeUnits);
  if (err != NC_NOERR) {
    printfone("Error: Could not get time units from NetCDF file %s: %s", 
              filename, nc_strerror(err));
    safe_free(NCInfo->timeUnits);
    return 1;
  }
  if (strcmp("picosecond",NCInfo->timeUnits) != 0)
    printfone("WARNING: Expecting time units in picoseconds, got -%s-\n", NCInfo->timeUnits);

  err =    nc_inq_varid(NCInfo->ncid, AMBER_NETCDF_SPATIAL, &NCInfo->spatialVID);
  if (err != NC_NOERR) {
    printfone("Error: Getting spatialVID in the NetCDF file %s: %s\n", 
            filename, nc_strerror(err));
    return 1;
  }

  /*
   *  NETCDF: check to see what trajectory data is within the NetCDF file
   *  and perform sanity checks...
   *
   *  ARE COORDINATES PRESENT?
   */ 
    err =    nc_inq_varid(NCInfo->ncid, AMBER_NETCDF_COORDS, &NCInfo->coordinateVID);
    if (err != NC_NOERR) {
      printfone("Error: No coordinates are present in the NetCDF file %s\n", filename);
      return 1;
    } 

    err =    nc_inq_varnatts(NCInfo->ncid, NCInfo->coordinateVID, &i);
    if ( err != NC_NOERR ) {
      printfone("Error: Getting number of coordinate attributes in NetCDF file %s\n", filename);
      return 1;
    }
    if ( i != 1 ) {
      printfone("Error: Only a single coordinate attribute is expected in NetCDF file %s\n", 
                filename);
      return 1;
    }

    err =    nc_inq_attlen(NCInfo->ncid, NCInfo->coordinateVID, "units", &ist);
    if (err != NC_NOERR) {
      printfone("Error: Getting coordinateVID attribute length. %s\n",nc_strerror(err));
      return 1;
    }    
    NCInfo->coordinateUnits = (char *) safe_malloc(sizeof(char) * (ist+1));

    err =    nc_get_att_text(NCInfo->ncid, NCInfo->coordinateVID, "units", NCInfo->coordinateUnits);
    if (err != NC_NOERR) {
      printfone("Error: Could not get coordinate units from NetCDF file %s", filename);
      safe_free(NCInfo->coordinateUnits);
      return 1;
    }
    if (strcmp("angstrom",NCInfo->coordinateUnits) != 0)
      printfone("WARNING: Expecting coordinate units in angstroms, got %s\n", 
                NCInfo->coordinateUnits);

  /*
   *  ARE CELL_LENGTHS (i.e. box info) PRESENT?
   */ 
    err =    nc_inq_varid(NCInfo->ncid, "cell_angles", &NCInfo->cellAngleVID);
    if ((err != NC_NOERR)&&(prnlev>0)) 
      printfone("Warning: NetCDF cell angle variable ID: %s\n", nc_strerror(err));
    err =    nc_inq_varid(NCInfo->ncid, "cell_lengths", &NCInfo->cellLengthVID);
    if ((err != NC_NOERR)&&(prnlev>0)) 
      printfone("Warning: NetCDF cell length variable ID: %s\n", nc_strerror(err));

    // Set up box information
    if (err == NC_NOERR) {
      // Angle information
      err =    nc_inq_varnatts(NCInfo->ncid, NCInfo->cellAngleVID, &i);
      if ( i > 0 ) {
        err =    nc_inq_attlen(NCInfo->ncid, NCInfo->cellAngleVID, "units", &ist);
        if (err == NC_NOERR) {
          NCInfo->cellAngleUnits = (char *) safe_malloc(sizeof(char) * (ist+1));
          err=nc_get_att_text(NCInfo->ncid, NCInfo->cellAngleVID, "units", NCInfo->cellAngleUnits);
        }
      }

      // Cell length information
      err =    nc_inq_varnatts(NCInfo->ncid, NCInfo->cellLengthVID, &i);
      if ( i > 0 ) {
        err =    nc_inq_attlen(NCInfo->ncid, NCInfo->cellLengthVID, "units", &ist);
        if (err == NC_NOERR) {
          NCInfo->cellLengthUnits = (char *) safe_malloc(sizeof(char) * (ist+1));
          err=nc_get_att_text(NCInfo->ncid, NCInfo->cellLengthVID, "units", 
                              NCInfo->cellLengthUnits);
        }
      }
      trajInfo->isBox = 1;
    } // End box information setup

  /*
   *  ARE VELOCITIES PRESENT?
   */ 

    err =    nc_inq_varid(NCInfo->ncid, "velocities", &NCInfo->velocityVID);
    if (err == NC_NOERR) {
      trajInfo->isVelocity = 1;
      err =    nc_inq_varnatts(NCInfo->ncid, NCInfo->velocityVID, &i);
      if ( i > 1 ) {
        err =    nc_inq_attlen(NCInfo->ncid, NCInfo->velocityVID, "units", &ist);
        if (err == NC_NOERR) {
          NCInfo->velocityUnits = (char *) safe_malloc(sizeof(char) * (ist+1));
          err = nc_get_att_text(NCInfo->ncid, NCInfo->velocityVID, "units", NCInfo->velocityUnits);
        }
        err = nc_get_att_double(NCInfo->ncid, NCInfo->velocityVID, "scale_factor", 
                                &NCInfo->velocityScale);
      }
    }

    /*
     * Are replica temperatures present?
     */
    err=   nc_inq_varid(NCInfo->ncid,"temp0",&NCInfo->TempVarID);
    if (err == NC_NOERR) {
      if (prnlev>0) printfone("\nNetCDF file has replica temperatures.\n");
    } else {
      if (prnlev>0) printfone("\nNetCDF file does not have replica temperatures.\n");
      NCInfo->TempVarID=-1;
    }

  return 0;
#endif //BINTRAJ

  return 1;
}

/*
 * NETCDF_setupOutput()
 * Given a coordinateInfo structure which has been opened for writing NETCDF
 * (NCInfo is allocated), set up the NCInfo structure.
 * Return 0 on success, 1 on failure.
 */
int NETCDF_setupOutput(coordinateInfo *trajInfo, int atoms) {
#ifdef BINTRAJ
  netcdfTrajectoryInfo *NCInfo;
  int dimensionID[NC_MAX_VAR_DIMS],err,oldMode;
#  ifdef MPI
  MPI_Offset start[3], count[3];
#  else
  size_t start[3], count[3];
#  endif
  char xyz[3];
  char abc[15] = { 'a', 'l', 'p', 'h', 'a',
                   'b', 'e', 't', 'a', ' ',
                   'g', 'a', 'm', 'm', 'a' };

  if (trajInfo->NCInfo==NULL) return 1;
  NCInfo = trajInfo->NCInfo;
  
  /*
   *  define the global dimensions of the NetCDF file
   */
  netcdfDefineDimension(NCInfo->ncid, AMBER_NETCDF_FRAME, NC_UNLIMITED, &NCInfo->frameDID);
  netcdfDefineDimension(NCInfo->ncid, AMBER_NETCDF_SPATIAL, 3, &NCInfo->spatialDID);
  netcdfDefineDimension(NCInfo->ncid, AMBER_NETCDF_ATOM, atoms, &NCInfo->atomDID);
  netcdfDefineDimension(NCInfo->ncid, AMBER_NETCDF_LABEL, AMBER_NETCDF_LABELLEN, &NCInfo->labelDID);
  netcdfDefineDimension(NCInfo->ncid, AMBER_NETCDF_CELL_SPATIAL, 3, &NCInfo->cell_spatialDID);
  netcdfDefineDimension(NCInfo->ncid, AMBER_NETCDF_CELL_ANGULAR, 3, &NCInfo->cell_angularDID);

  /*
   *  put global attributes
   */
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "title", trajInfo->title);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "application", trajInfo->application);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "program", trajInfo->program);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "programVersion", trajInfo->version);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "Conventions", "AMBER");
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "ConventionVersion", "1.0");

  /*
   *  handle definition of non-optional variables
   */
  dimensionID[0] = NCInfo->spatialDID;
  netcdfDefineVariable(NCInfo->ncid, AMBER_NETCDF_SPATIAL, NC_CHAR, 1, dimensionID, &NCInfo->spatialVID);

  dimensionID[0] = NCInfo->frameDID;
  netcdfDefineVariable(NCInfo->ncid, AMBER_NETCDF_TIME, NC_FLOAT, 1, dimensionID, &NCInfo->timeVID);
  netcdfPutAttributeText(NCInfo->ncid, NCInfo->timeVID, "units", "picosecond");

  dimensionID[0] = NCInfo->frameDID;
  dimensionID[1] = NCInfo->atomDID;
  dimensionID[2] = NCInfo->spatialDID;
  netcdfDefineVariable(NCInfo->ncid, AMBER_NETCDF_COORDS, NC_FLOAT, 3, dimensionID, &NCInfo->coordinateVID);
  netcdfPutAttributeText(NCInfo->ncid, NCInfo->coordinateVID, "units", "angstrom");

  dimensionID[0] = NCInfo->cell_spatialDID;
  netcdfDefineVariable(NCInfo->ncid, AMBER_NETCDF_CELL_SPATIAL, NC_CHAR, 1, dimensionID, &NCInfo->cell_spatialVID);

  dimensionID[0] = NCInfo->cell_angularDID;
  dimensionID[1] = NCInfo->labelDID;
  netcdfDefineVariable(NCInfo->ncid, AMBER_NETCDF_CELL_ANGULAR, NC_CHAR, 2, dimensionID, &NCInfo->cell_angularVID);

  // Set up Box coords
  if (trajInfo->isBox) {

    dimensionID[0] = NCInfo->frameDID;
    dimensionID[1] = NCInfo->cell_spatialDID;
    netcdfDefineVariable(NCInfo->ncid, "cell_lengths", NC_DOUBLE, 2, dimensionID, &NCInfo->cellLengthVID);
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->cellLengthVID, "units", "angstrom");

    dimensionID[1] = NCInfo->cell_angularDID;
    netcdfDefineVariable(NCInfo->ncid, "cell_angles", NC_DOUBLE, 2, dimensionID, &NCInfo->cellAngleVID);
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->cellAngleVID, "units", "degree");

  }

  /* DAN ROE: Replica temperature */
  if (trajInfo->isREMDTRAJ) {
    fprintf(stdout,"NETCDF: Defining replica temperature in output trajectory.\n");
    dimensionID[0] = NCInfo->frameDID;
    netcdfDefineVariable(NCInfo->ncid, "temp0", NC_DOUBLE, 1, dimensionID, &NCInfo->TempVarID);
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->TempVarID,"units","kelvin");
  }

  /*
   *  set fill mode
   */
  err = nc_set_fill(NCInfo->ncid, NC_NOFILL, &oldMode);
  if (err != NC_NOERR) {
    printfone("NetCDF setting fill value: %s\n", nc_strerror(err));
    return 1;
  }

  /*
   *  end of NetCDF definitions
   */
  err = nc_enddef(NCInfo->ncid);
  if (err != NC_NOERR) {
    printfone("NetCDF error on ending definitions: %s\n", nc_strerror(err));
    return 1;
  }

  /*
   *  specify spatial dimension labels
   */
  start[0] = 0;
  count[0] = 3;

  xyz[0] = 'x';
  xyz[1] = 'y';
  xyz[2] = 'z';
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->spatialVID, start, count, xyz);
  if (err != NC_NOERR) {
    printfone("Error on NetCDF output of spatial VID 'x', 'y' and 'z': %s\n",
              nc_strerror(err));
    return 1;
  }

  xyz[0] = 'a';
  xyz[1] = 'b';
  xyz[2] = 'c';
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->cell_spatialVID, start, count, xyz);
  if (err != NC_NOERR) {
    printfone("Error on NetCDF output of spatial VID 'x', 'y' and 'z': %s\n",
              nc_strerror(err));
    return 1;
  }

  start[0] = 0;
  start[1] = 0;
  count[0] = 3;
  count[1] = 5;
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->cell_angularVID, start, count, abc);
  if (err != NC_NOERR) {
    printfone("Error on NetCDF output of angular VID 'alpha', 'beta ' and 'gamma': %s\n",
              nc_strerror(err));
    return 1;
  }

  //outInfo->NCInfo = NCInfo;
  return 0;

#endif
  return 1;
}

/*
 * The functions below should eventually only be required within this file.
 */

#ifdef BINTRAJ
/*
 * netcdfDefineVariable()
 *  PF - Multiptraj
 *  New method called in ptraj.c used to write a netcdf file.
 *  Checks to see if using MPI to call correct method
 * return 0 on success, 1 on failure.
 */
int netcdfDefineVariable( int ncid, char *name, nc_type xtype, int ndims, int dimids[], int *varidp ) {
  int err;

  err = nc_def_var(ncid, name, xtype, ndims, dimids, varidp);

  if (err == NC_NOERR) 
    return 0;
  else
    fprintf(stdout, "netcdfDefineVariable: Error defining variable (%s): %s\n",
            name, nc_strerror(err));
  return 1;
}

/*
 *  PF - Multiptraj
 *  New method called in ptraj.c used to write a netcdf file.
 *  Checks to see if using MPI, and if so use ncmpi instead
 *  from library parallel-netcdf
 */

   void
netcdfDefineDimension( int ncid, char *name, int length, int *dimidp )
{
  int err;

  err = nc_def_dim(ncid, name, length, dimidp);
  
  if (err != NC_NOERR) {
    fprintf(stdout, "netcdfDefineDimension: Error defining dimension (%s): %s\n",
            name, nc_strerror(err));
  }
}


   char *
netcdfGetAttributeText( int ncid, int vid, char *attribute )
{
  int err;
  char *returnText;
#  ifdef MPI
  MPI_Offset ist;
#  else
  size_t ist;
#  endif

  err = nc_inq_attlen(ncid, vid, attribute, &ist);
  if (err == NC_NOERR) {
    returnText = (char *) safe_malloc(sizeof(char) * (ist+1));
    err = nc_get_att_text(ncid, vid, attribute, returnText);
    if (err != NC_NOERR) {
      fprintf(stdout, "netcdfGetAttributeText: Could not get attribute (%s) from NetCDF file", attribute);
      safe_free(returnText);
      returnText = NULL;
    }
  } else
    fprintf(stdout, "netcdfGetAttributeText: Error grabbing attribute (%s): %s\n", 
            attribute, nc_strerror(err));

  return returnText;

}

/*
 *  PF - Multiptraj
 *  Checks to see if using MPI to call correct method
 */

   void
netcdfPutAttributeText( int ncid, int vid, char *attribute, char *text )
{
  int err;

  err = nc_put_att_text(ncid, vid, attribute, strlen(text), text);

  if (err != NC_NOERR) {
    fprintf(stdout, "netcdfPutAttributeText: Error putting attribute (%s): %s\n", 
            attribute, nc_strerror(err));
  }
}

   int
netcdfGetDimensionInfo(int ncid, char *attribute, int *length)
{
  int err, dimID;
#  ifdef MPI
  MPI_Offset slength;
#  else
  size_t slength;
#  endif

  *length = 0;

  err = nc_inq_dimid(ncid, attribute, &dimID);
  if (err != NC_NOERR) {
    fprintf(stdout, "Yikes! netcdfGetDimensionInfo on %s, error: %s\n",
             attribute, nc_strerror(err));
  } else {

    err = nc_inq_dimlen(ncid, dimID, &slength);
    if (err != NC_NOERR) {
      fprintf(stdout, "Do'h!!! netcdfGetDimensionInfo length on %s, error: %s\n",
              attribute, nc_strerror(err));
      slength = 0;
    }
  }
  *length = (int) slength;
  return dimID;
}


#endif
