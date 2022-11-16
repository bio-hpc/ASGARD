#include "config.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream> 

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

using namespace std;
#include "emil.h"

const char *usageAndLicencingMessage=
"##################################################################################\n"
"#\n"
"# EMIL module for calculation of absolute free energies of biomolecules.\n"
"# By Joshua T Berryman, University of Luxembourg, 2013.\n"
"#\n"
"# This code is covered by the GPL, except when distributed with non-gpl parts of\n"
"# AMBER (when it is covered by the same licence as the rest of that distribution).\n"
"# See licence.txt in the source directory.\n"
"#\n"
"# If you use this module for research, please cite:\n"
"#                         Joshua T Berryman and Tanja Schilling,\n"
"#                         J. Chem. Theory Comput., 2013, 9 (1), 679-686.\n"
"#\n"
"#\n"
"##################################################################################\n";

#if 0
//NAB functions need this fancy declare
extern "C" {
  extern MOLECULE_T *loadAmberTopCrd(STRING_T **crdFileName, STRING_T **topFileName, STRING_T **pdbFileName );
  extern MOLECULE_T *copymolecule(MOLECULE_T *m);
  extern int         setxyz_from_mol( MOLECULE_T **m, char **aex, POINT_T xyz[] );
  extern void        initFrameUtils();
  extern void        cenmas_noImg( int, int [], REAL_T [], REAL_T *, REAL_T *, REAL_T *, REAL_T * );
}
#endif


void hsc::setupList(){

  int     icell;
  double  invxCell[3], rc;
  
  box->update();

  rc = rcutLiquid; //set neighbourlist cell size == interaction cutoff distance
  
  //the same cell list serves for soft force and well swapping
  if ( rc < SF_MAX_R ) rc = SF_MAX_R;

  // dimensions
  for (int i=0;i<3;i++){
    box->nxCell[i] = int(box->x[i]/rc);
    box->xCell[i] = box->x[i]/box->nxCell[i];
  }
  box->nCells = box->nxCell[0]*box->nxCell[1]*box->nxCell[2];

  if ((box->nxCell[0] < 3)||(box->nxCell[1] < 3) ||(box->nxCell[2] < 3)){
    MERR << "EMIL Error: at least one box dimension is less than 3 cells" << endl;
    MERR << "Box size: " << box->nxCell[0] << "x" 
                         << box->nxCell[1] << "x"
                         << box->nxCell[2]<< " cells." << endl;
    exit(8);
  }

  cells = new Cell[box->nCells];
  for (int i=0; i<box->nCells; i++){
      cells[i].myId = i;
  }
  
  for(int i = 0; i < 3; i++){
     invxCell[i] = 1.0 / box->xCell[i];
  }
  
  // put particles into cells
  for (int i=0; i<N; i++){

    icell = selectCell(part[i].R, invxCell );

    if (icell >= 0 && icell < box->nCells) {
      part[i].insertToCell(cells[icell]); 
    }
    else if (icell == box->nCells) {
      part[i].insertToCell(cells[icell-1]);
    }
    else {
	MERR << "EMIL Error: problem with particle->cell assignment!" << endl;
	MERR << "Number of cell " << icell << " of " << box->nCells  << endl;
	MERR << i << " Coordinates: " << part[i].R[0] << " " << part[i].R[1] 
	     << " " << part[i].R[2] << endl; 
	exit(8);
    }
  }

  // find neighbouring cells
  setUpNeighbours();
}

void hsc::setUpNeighbours(){

  int x,y,z,count,cl,cln;
   
  for (int k=0; k<box->nxCell[2]; k++){
    for (int j=0; j<box->nxCell[1]; j++){
      for (int i=0; i<box->nxCell[0]; i++){

	count = 0;
	for (int zs = -1; zs <= 1; zs++){

	  // Cheap hack to get positive values only from the modulo operator 
	  z = (k+zs + 100*box->nxCell[2])%box->nxCell[2]; 

	  for (int ys = -1; ys <= 1; ys++){
	    y = (j+ys + 100*box->nxCell[1])%box->nxCell[1];

	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell
	      if ( xs!=0 || ys!=0 || zs!=0 ){

		x = (i+xs + 100*box->nxCell[0])%box->nxCell[0];

		cl  = i + box->nxCell[0]*j + box->nxCell[1]*box->nxCell[0]*k;
		cln = z*box->nxCell[1]*box->nxCell[0] + y*box->nxCell[0] + x;
		cells[cl].neighbours[count] = &cells[cln];
		count ++;

	      }
	    }
	  }
	}

      }
    }
  } 
}


// ATTENTION: box needs to fit three rcut in each dimension!

void hsc::setupListW(){

  int icell;
  double rc = rcutLiquid; //(Here: cut off radius of wells)
  double invxCellW[3];

  //the same cell list serves for soft force and well swapping
  if ( rc < SF_MAX_R ) rc = SF_MAX_R;

  // dimensions
  for (int i=0;i<3;i++){
    box->nxCellW[i] = int(box->x[i]/rc);
    box->xCellW[i]  = box->x[i]/box->nxCellW[i];
    invxCellW[i]    = 1.0 / box->xCellW[i];
  }
  box->nCellsW = box->nxCellW[0]*box->nxCellW[1]*box->nxCellW[2];

  if ((box->nxCellW[0] < 2)||(box->nxCellW[1] < 2)||(box->nxCellW[2] < 2)){
    MERR << "EMIL Error: Box dimension is less than 2 times rcut" << endl;
    MERR << "rcut is: " << rc << endl;
    MERR << "box is: " << box->x[0] << "x" << box->x[1] << "x" << box->x[2] << endl;
    exit(8);
  }

  cellsW = new Cell[box->nCellsW];
  for (int i=0; i<box->nCellsW; i++){
      cellsW[i].myId = i;
  }

  // put particles into cells
  for (int i=0; i<N; i++){
    
    icell = selectCellW( part[i].R, invxCellW );
    if ( icell < box->nCellsW && icell >= 0 ) part[i].insertToCellW(cellsW[icell]);
    else if (icell == box->nCellsW) part[i].insertToCellW(cellsW[icell-1]);
    else {
	MERR << "EMIL Error: problem with particle->cell assignment!" << endl;
	MERR << "Number of well-cell " << icell << " of " 
	     << box->nCellsW  << endl;
	MERR << "Coordinates: " << part[i].R[0] << " " << part[i].R[1] 
	     << " " << part[i].R[2] << endl; 
	exit(8);
    }
  }

   // assign cells to wells and wells to cells
  for (int i=0; i<N; i++){
    
    icell = selectCellW( wells[i]->R, invxCellW );
    
    if (icell < box->nCellsW){
      wells[i]->cell = &cellsW[icell];
      wells[i]->next          = cellsW[icell].firstWell;
      cellsW[icell].firstWell = wells[i];
    }
    else if (icell == box->nCellsW) {
      wells[i]->cell              = &cellsW[icell - 1];
      wells[i]->next              = cellsW[icell - 1].firstWell;
      cellsW[icell - 1].firstWell = wells[i];
    }
    else {
	MERR << "EMIL Error: problem with well->cell assignment!" << endl;
	MERR << "Number of well-cell " << icell << " of " 
	     << box->nCellsW  << endl;
	MERR << "Coordinates: " << wells[i]->R[0] << " " << wells[i]->R[1] 
	     << " " << wells[i]->R[2] << endl; 
	exit(8);
    }
  }
  // find neighbouring cells
  setUpNeighboursW();

}

void hsc::setUpNeighboursW(){

  int x,y,z,count,cl,cln;
   
  for (int k=0; k<box->nxCellW[2]; k++){
    for (int j=0; j<box->nxCellW[1]; j++){
      for (int i=0; i<box->nxCellW[0]; i++){
	count = 0;
	for (int zs = -1; zs <= 1; zs++){
	  // Cheap hack to get positive values only from the modulo operator 
	  z = (k+zs + 100*box->nxCellW[2])%box->nxCellW[2]; 
	  for (int ys = -1; ys <= 1; ys++){
	    y = (j+ys + 100*box->nxCellW[1])%box->nxCellW[1];
	    for (int xs = -1; xs <= 1; xs++){
	      // don't count the middle cell
	      if (xs!=0 || ys!=0 || zs!=0){
		x = (i+xs + 100*box->nxCellW[0])%box->nxCellW[0];
		cl = i + box->nxCellW[0]*j + box->nxCellW[1]*box->nxCellW[0]*k;
		cln = z*box->nxCellW[1]*box->nxCellW[0] + y*box->nxCellW[0] + x;
		cellsW[cl].neighbours[count] = &cellsW[cln];
		count ++;
	      }
	    }
	  }
	}
      }
    }
  } 
}


void hsc::setUpNew(double *crd, double *frc, double *vel, double *boxVec){

  ifstream infile;
  int      p, i;
  bool     haveBox; 
    
#ifndef SSTI_BUFFER_LOG
  //open the logfile
  if( logFileName.size() > 0 ){
    logFile = new ofstream(logFileName.c_str());
    if( logFile->bad() ){
      MERR << "Error, could not open logfile, name: " << logFileName << endl;
    }else{
      MERR << "Directing output to file: " << logFileName << endl;
    }
  }else{
      MERR << "Directing output to stdout." << endl;
      logFile = &cout;
  }
#else
  logFile    = new stringstream();
#ifdef USE_MPI 
  if( myTaskId == 0 ){ //assume 0==master.
#endif
    logFile_fp = fopen( logFileName.c_str(), "w");
    if( !logFile_fp ){
      MERR << "Failed to open logfile: " << logFileName << endl;
      exit( 8 );
    }
#ifdef USE_MPI 
  }
#endif
#endif
 *logFile << usageAndLicencingMessage;
  flushLog();  


  //setup the box
  box = new Box;
  haveBox  = 0;
  //periodic = 0;
  periodic = 1;
  if( boxVec ){
    if( boxVec[0] > 0.0 && boxVec[1] > 0.0 && boxVec[2] > 0.0 ){
        *logFile << "# Setting up emil, with boxinfo: " << boxVec[0] << " " << boxVec[1] << " " << boxVec[2] << endl;
        MERR << "# Setting up emil, with boxinfo: " << boxVec[0] << " " << boxVec[1] << " " << boxVec[2] << endl;
        box->x[0] = boxVec[0];
        box->x[1] = boxVec[1];
        box->x[2] = boxVec[2];
        box->update();
        haveBox  = 1;
        periodic = 1;
    }
  }

  //setup the structures which hold the particle info
  part           = new Particle[N];

  // setup the array of floats which holds coordinates of each atom.
  restraintCoords = new double[3*N];//don't need restraint minima for every single atom, but its easier for indexing.
  particleCoords  = crd;
  forces          = frc;
  velocities      = vel;
  i = 0;
  for(p = 0; p < N; p++){
    part[p].R    = &particleCoords[i];
    part[p].myId = p;
    i           += 3;
  }
  

  //Rather than write separate code for no PBCs, just invent a box
  //into which the molecule will fit. Note that can have box but not be periodic!
  if( haveBox == 0 ){
    double maxR[3], minR[3];
    maxR[0] = particleCoords[0];
    maxR[1] = particleCoords[1];
    maxR[2] = particleCoords[2];
    minR[0] = particleCoords[0];
    minR[1] = particleCoords[1];
    minR[2] = particleCoords[2];
    for(i = 3; i < 3*N; ){
        if( particleCoords[i] > maxR[0] ) maxR[0] = particleCoords[i];
        if( particleCoords[i] < minR[0] ) minR[0] = particleCoords[i];
        i++;
        if( particleCoords[i] > maxR[1] ) maxR[1] = particleCoords[i];
        if( particleCoords[i] < minR[1] ) minR[1] = particleCoords[i];
        i++;
        if( particleCoords[i] > maxR[2] ) maxR[2] = particleCoords[i];
        if( particleCoords[i] < minR[2] ) minR[2] = particleCoords[i];
        i++;
    }


    /* this is not a very futureproof way of finding the longest of the 
    ** various cutoffs in the system */
    rcut_longest  = rcutLiquid;
    if( rcutSolid > rcut_longest)  rcut_longest = rcutLiquid;
    if( SF_MAX_R  > rcut_longest)  rcut_longest = SF_MAX_R;
    rcut_longest2 = rcut_longest * rcut_longest;

    box->x[0] = maxR[0] - minR[0] + 2*rcut_longest;
    box->x[1] = maxR[1] - minR[1] + 2*rcut_longest;
    box->x[2] = maxR[2] - minR[2] + 2*rcut_longest;

    if( box->x[0] < 3.0*rcut_longest) box->x[0] = 3.0*rcut_longest;
    if( box->x[1] < 3.0*rcut_longest) box->x[1] = 3.0*rcut_longest;
    if( box->x[2] < 3.0*rcut_longest) box->x[2] = 3.0*rcut_longest;


   *logFile << "# Setting up emil, with fake boxinfo: " << box->x[0] << " " << box->x[1] << " " << box->x[2] << endl;
   MERR << "# Setting up emil, with fake boxinfo: " << box->x[0] << " " << box->x[1] << " " << box->x[2] << endl;
  }


  iran = &seed;
 
 *logFile<< "# read parameters:" << endl;
 *logFile<< "# N " << N << endl;
 *logFile<< "# seed " << seed << " Displ " << maxDisplace << endl;
  
  flushLog();

  initialize();
  
  //we need this mask for rigid-body fitting constraint sets to solute molecules.
  //initially we support only one solute molecule... later we will have to have >1 values allowed for the mask.
  soluteMask       = new int[N];
  topSoluteAtom    = 0;
  for( i = 0; i < N; i++ ){
    if( part[i].isRoot ){
      if(part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL ){
        soluteMask[i] = 1; 
        topSoluteAtom = i;
      }else{
        soluteMask[i] = 0;
      }
    }
  }
  
  #if 0
    initFrameUtils();
  #endif
  
  //cenmas_noImg( topSoluteAtom, soluteMask, particleCoords, &x, &y, &z, &mtot );
  
  //image particles s.t. the solute geometric centre is at the origin
  //for(p = 0; p < N; p++){
  //  particleCoords[3 * p]     -= x;
  //  particleCoords[3 * p + 1] -= y;
  //  particleCoords[3 * p + 2] -= z;
  //}
  
 *logFile<< "# done setup.\n#" << endl;
  flushLog();
}

  
void hsc::setUpFromFile() {

  ifstream infile2;
  char line[200];
  char quantity[200];
  float number;
  int check, i, p;

  // read configuration from ConfigFile
  int allSet = 0;
  int count = 0;
  double Rin[3];

  infile2.open(configFile.c_str());
  
  if (infile2.bad()) {
      MERR << "EMIL Error: " << configFile << " not found\n";
      exit(8);
  }
  box = new Box;
 
  while (infile2.peek() != EOF) {
      infile2.getline(line,sizeof(line));
      if(allSet == 5){

	check = sscanf(line, "%le %le %le", &Rin[0], &Rin[1], &Rin[2]);

	if (check == 3) { 
	  if (count < N){
	    for (int j=0; j<3; j++){
	      part[count].R[j] = Rin[j];
	    }
	    count++;
	  }
	  else {
	    MERR << "EMIL Error: Too many lines in " << configFile << endl;
	    exit(8);
	  }
	}
	
	else { 
	  MERR << "EMIL Error: wrong line format in " << configFile << endl;
      MERR << "Line was: \n\"" << line << "\"" << endl;
	  exit(8);
	}
      }	

      else {
	check = sscanf(line, "%s%f", quantity, &number);
	if (check != 2) {
	  MERR << "EMIL Error: wrong line format in " << configFile << endl;
      MERR << "EMIL Line was: \n\"" << line << "\"" << endl;
      MERR << "EMIL Counted: " << check << " fields when expected 1 string and 1 float." << endl;
	  exit(8);
	}
	if (strcmp(quantity,"#Number")==0) {
	  N=int(number); 	 
	  part           = new Particle[N];
	  particleCoords = new double[3*N];
	  i = 0;
	  for(p = 0; p < N; p++){
	    part[p].R = &particleCoords[i+=3];
	  }
	  
	  allSet++;}
	else if (strcmp(quantity,"#Density")==0) {allSet++;} 
	else if (strcmp(quantity,"#Boxx")==0) {box->x[0]=number; allSet++;}
	else if (strcmp(quantity,"#Boxy")==0) {box->x[1]=number; allSet++;}
	else if (strcmp(quantity,"#Boxz")==0) {box->x[2]=number; allSet++;}
      }
  }
  
  infile2.close();

  if (count==0) {
    MERR << "EMIL Error: parameter missing in " << configFile << endl;
    exit(8);
  }

 *logFile<< "# read file: " << configFile << endl;

  iran = &seed;

  initialize();

 *logFile<< "# read parameters and configuration" << endl;
 *logFile<< "# N " << N << endl;
 *logFile<< "# Eq " << eqSteps << " EqEval " << eqEvalSteps << " EqSnap " 
       << eqSnapSteps << endl;
 *logFile<< "# Mc " << mcSteps << " McEval " << mcEvalSteps << " McSnap " 
       << mcSnapSteps << endl;
 *logFile<< "# seed " << seed << " Displ " << maxDisplace << endl;

}




