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

#include "emil.h"

using namespace std;

void hsc::readParamFile(){

  ifstream infile;
  char line[200];
  char quantity[200];
  char value[200], *p;
  int check;
  int lineCount;
  
  infile.open(parameterFile.c_str());
  if (!infile.is_open()) {
    MERR << "Error: parameterFile '" << parameterFile << "' not opened/found\n";
    exit(8);
  }
  MERR << "# Param file: " << parameterFile << endl;

  //init this here.
  numMaterials     = 0;
  materialType     = new HSC_POTENTIAL_T[64];
  materialAtomMask = new string[64];
  
  //set defaults
  bx              =   0.0;
  seed            =   1;
  assignmentMethod = HSC_ASSIGNMENT_MC;
  assignmentEvery =   1;
  mcSteps         =   0;
  mcEvalSteps     =   0;
  mcSnapSteps     =   0;
  eqSteps         =   0;
  eqEvalSteps     =   0;
  eqSnapSteps     =   0;
  maxDisplace     =   0.2;
  epsilonSolid    =   5.0;
  epsilonTrap     =   5.0;
  wingForce       =   1.111;
  rcutSolid       =   5.0;
  rcutLiquid      =   5.0;
  req             =   0.03;
  swapTriesPerChain  = 0.0;
  relocTriesPerChain = 0.0;
  targetSwapRate  =   0.01;
  targetRelocRate =   0.001; 
  softcore_scale  =   0.8;  
  saveWellsEvery  =   1000;
  printEvery      =     10;
  stepLambdaBackwards = false;
  fixLambda           = false;
  
  lineCount = 0;
  while (infile.peek() != EOF) {
    infile.getline(line,sizeof(line));
    check = sscanf(line, "%s%s", quantity, value);
  
    MERR << "read: " << line << endl;
    if (check != 2 || line[0] == '#') {
      //pass through any comments or empty lines, otherwise complain.
      if( check != -1 && line[0] != '#' ){
        MERR << check << " Warning! Wrong line format "<< line << " in file: " << parameterFile << endl;
      }
      continue;
    }
    lineCount++;

    //set the key word lowercase
    for (p=quantity; *p; ++p) *p = tolower(*p);
 
    //do some pattern matching
    if (strcmp(quantity,"number")==0) {
      N=atoi(value);
    }
    else if (strcmp(quantity,"assignmentmethod")==0) {
        assignmentMethod=(HSC_ASSIGNMENT_T)atoi(value);
        if( assignmentMethod < 0 || assignmentMethod >= HSC_ASSIGNMENT_NUMASSIGNMENTS){
      MERR << "EMIL Error! Assignment algorithm not recognised:\n" << line << endl;
      MERR << "Exiting.";
      exit(8);
        }
    }

    else if (strcmp(quantity,"assignmentevery")==0) assignmentEvery=atoi(value);
    else if (strcmp(quantity,"boxside")==0) bx=atof(value);
    else if (strcmp(quantity,"seed")==0) seed=atoi(value); 
    else if (strcmp(quantity,"mcsteps")==0) mcSteps=atoi(value);
    else if (strcmp(quantity,"mcevalsteps")==0) mcEvalSteps=atoi(value);
    else if (strcmp(quantity,"mcsnapsteps")==0) mcSnapSteps=atoi(value);
    else if (strcmp(quantity,"eqsteps")==0) eqSteps=atoi(value); 
    else if (strcmp(quantity,"eqevalsteps")==0) eqEvalSteps=atoi(value);
    else if (strcmp(quantity,"eqsnapsteps")==0) eqSnapSteps=atoi(value);  
    else if (strcmp(quantity,"maxdisplacement")==0) maxDisplace=atof(value);
    else if (strcmp(quantity,"targetswaprate")==0)  targetSwapRate = atof(value);
    else if (strcmp(quantity,"targetrelocrate")==0)  targetRelocRate = atof(value);
    else if (strcmp(quantity,"epsilonwell")==0) epsilonSolid=atof(value);
    else if (strcmp(quantity,"rwell")==0)       rcutSolid=atof(value);
    else if (strcmp(quantity,"epsilontrap")==0) epsilonTrap=atof(value);
    else if (strcmp(quantity,"rtrap")==0)       rcutLiquid=atof(value);
    else if (strcmp(quantity,"wingforce")==0)   wingForce=atof(value);
    else if (strcmp(quantity,"reqtrap")==0)     req=atof(value);
    else if (strcmp(quantity,"pdbfilename")==0) pdbFileName.assign(value); 
    else if (strcmp(quantity,"topfilename")==0) topFileName.assign(value);
    else if (strcmp(quantity,"crdfilename")==0) crdFileName.assign(value);
    else if (strcmp(quantity,"logfilename")==0) logFileName.assign(value);
    else if (strcmp(quantity,"savewellsevery")==0) saveWellsEvery=atoi(value);
    else if (strcmp(quantity,"printevery")==0)     printEvery=atoi(value);
    else if (strcmp(quantity,"lambda")==0) {
      amberLambda = atof(value);
      restLambda  = atof(value);
      fixLambda   = true;
    }
    else if (strcmp(quantity,"steplambdabackwards")==0) stepLambdaBackwards=(bool)atoi(value);
    else if (strcmp(quantity,"trajoutfilename")==0) trajOutFileName.assign(value);
    else if (strcmp(quantity,"crdoutfilename")==0) crdOutFileName.assign(value);
    else if (strcmp(quantity,"swaptriesperchain")==0) swapTriesPerChain=atof(value);
    else if (strcmp(quantity,"reloctriesperchain")==0) relocTriesPerChain=atof(value);
    else if (strcmp(quantity,"softcorescale")==0) softcore_scale=atof(value);
    else if (strcmp(quantity,"solidres")==0)  {
      materialType[numMaterials] = HSC_POTENTIAL_EINSTEIN_CRYSTAL;
      materialAtomMask[numMaterials].assign(value);
      numMaterials++;
    }
    else if (strcmp(quantity,"liquidres")==0){
      materialType[numMaterials] = HSC_POTENTIAL_WINGWELL_LIQUID;
      materialAtomMask[numMaterials].assign(value);
      numMaterials++;
    }
    else if (strcmp(quantity,"assignedliquidres")==0){
      materialType[numMaterials] = HSC_POTENTIAL_HARMONIC_LIQUID;
      materialAtomMask[numMaterials].assign(value);
      numMaterials++;
    }
    else if (strcmp(quantity,"wobbliumres")==0){
      materialType[numMaterials] = HSC_POTENTIAL_WOBBLIUM;
      materialAtomMask[numMaterials].assign(value);
      numMaterials++;
    }
    else{
      MERR << "EMIL Error! Line not recognised:\n" << line << endl;
      MERR << "Exiting.";
      exit(8);
    }
  }

  MERR << "read " << lineCount << " key-value pairs from " << parameterFile << endl;

  infile.close();



}

void hsc::writeOutConfig(char name[16]){

  char temp[16];
  char allname[50];

  strcpy(allname, name);
  strcat(allname, "_");
  //  sprintf(temp,"%.*lf", 3, N/box->V);
  //  strcat(allname, temp);
  //  strcat(allname, "_");
  strcat(allname, temp);

  ofstream outfile(allname, ios::out);
  outfile.precision(10);

  outfile << "#Number " << N << endl;
  outfile << "#Density " << N/box->V << endl;
  outfile << "#Boxx " << box->x[0] << endl;
  outfile << "#Boxy " << box->x[1] << endl;
  outfile << "#Boxz " << box->x[2] << endl;

  for (int i=0; i<N; i++){
    outfile << part[i].R[0] << " " << part[i].R[1]  << " "
            << part[i].R[2] <<  endl;
  }

  outfile.close();
}


void hsc::writeOutWells(){

  ofstream outfile( wellSaveFileName.c_str(), ios::out );
  
  outfile.precision(10);

  outfile << "#Number " << N << endl;
  outfile << "#Density " << N/box->V << endl;
  outfile << "#Boxx " << box->x[0] << endl;
  outfile << "#Boxy " << box->x[1] << endl;
  outfile << "#Boxz " << box->x[2] << endl;

  for (int i=0; i<N; i++){
    
    outfile << part[i].theWell->R[0] << " " << part[i].theWell->R[1]  << " "
            << part[i].theWell->R[2] <<  endl;
  }

  outfile.close();
}


void hsc::writeOutBinding(char name[16]){

  char temp[16];
  char allname[50];

  strcpy(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 3, N/box->V);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, epsilonTrap );
  strcat(allname, temp);

  ofstream outfile(allname, ios::out);
  outfile.precision(10);

  outfile << "#Number " << N << endl;
  outfile << "#Density " << N/box->V << endl;
  outfile << "#Boxx " << box->x[0] << endl;
  outfile << "#Boxy " << box->x[1] << endl;
  outfile << "#Boxz " << box->x[2] << endl;
  outfile << "# Particle xyz Well xyz Diff " << endl;

  for (int i=0; i<N; i++){   
    outfile << part[i].R[0] << " " << part[i].R[1] << " " << part[i].R[2] 
      << " " << part[i].theWell->R[0] << " " 
      << part[i].theWell->R[1]  << " "
      << part[i].theWell->R[2] << " " << sqrt(part[i].rijsqW) << endl;
  }
  
  outfile.close();
}

void hsc::writeOutBindingToWells(char name[16]){

  double dist;
  char temp[16];
  char allname[50];

  strcpy(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 3, N/box->V);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, epsilonTrap );
  strcat(allname, temp);

  ofstream outfile(allname, ios::out);
  outfile.precision(10);

  outfile << "#Number " << N << endl;
  outfile << "#Density " << N/box->V << endl;
  outfile << "#Boxx " << box->x[0] << endl;
  outfile << "#Boxy " << box->x[1] << endl;
  outfile << "#Boxz " << box->x[2] << endl;
  outfile << "# Well xyz Particlexyz Diff " << endl;

  for (int i=0; i<N; i++){

    dist = computeDistance( *(wells[i]->theParticle) );
    
    outfile << wells[i]->R[0] << " " << wells[i]->R[1] << " " << wells[i]->R[2] 
    << " " << wells[i]->theParticle->R[0] << " " 
    << wells[i]->theParticle->R[1]  << " "
    << wells[i]->theParticle->R[2] << " " << sqrt(dist) << endl;
  }
  
  outfile.close();
}

void hsc::writeOutSimData(char name[16]){

  char temp[16];
  char allname[50];

  strcpy(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 3, N/box->V);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, epsilonTrap);
  strcat(allname, temp);

  ofstream outfile(allname, ios::out);

  outfile << "# N: " << N << endl;
  outfile << "# Box x, y, z: " << box->x[0] << " " << box->x[1] 
          << " " << box->x[2] << endl;
  outfile << "# Epsilon trap: " << epsilonTrap << endl;
  outfile << "# Model: linear" << endl;
  outfile << "# Cutoff radius: " << rcutSolid << endl;
  outfile << "# Translation acceptance rate: " 
          << float(transAccept)/float(transMoves) << endl;
  outfile << "# Swap acceptance rate: " 
          << float(swapAccept)/float(swapMoves) << endl;
  outfile << "# Reloc acceptance rate: " 
          << float(relocAccept)/float(relocMoves) << endl;
  outfile << "# <Phi/N> = " << meanPhi/float(measureCount * N) << endl; 
  outfile << "# <Distance from well> = "  
          << meanDistanceFromWell/float(measureCount) << endl; 
  outfile << "# Ratio of particles outside wells = " 
          << meanRatioOutside/float(measureCount) << endl;
  outfile.close();

}

void hsc::flushLog(){
  
#ifdef SSTI_BUFFER_LOG //flush the output buffer
    if( myTaskId == 0 ){
      fprintf(logFile_fp, "%s", (char *)logFile->str().c_str() );
    }
    delete logFile;
    logFile = new stringstream();
#endif

}







