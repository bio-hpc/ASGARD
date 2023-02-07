/*!
 \file gaussianParser.cpp
 \brief Parses gaussian files
 \author Martin Peters
 
 Reads and writes guassian files
 
 $Date: 2010/08/19 13:48:48 $
 $Revision: 1.20 $
 
 ----------------------------------------------------------------------------
 
 MTK++ - C++ package of modeling libraries.
 
 Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)
 
 This file is part of MTK++.
 
 MTK++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MTK++ is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lessser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 ----------------------------------------------------------------------------
 */
#include <sstream>

#include "gaussianParser.h"
#include "zmatParser.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"
#include "Molecule/bond.h"
#include "Molecule/angle.h"
#include "Molecule/torsion.h"
#include "Molecule/improper.h"
#include "Utils/vector3d.h"

#include "Statistics/sheet.h"
#include "Statistics/table.h"

#include "StringManip.h"
#include "Utils/index.h"

// #include "Utils/diagonalize.h"

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;

#include "Utils/diagonalize_eigen.h"

#include "Log/errorHandler.h"
#include "Diagnostics/MTKException.h"

namespace MTKpp
{
	
	// ============================================================
	// Function : gaussianParser()
	// ------------------------------------------------------------
	// Constructor for the class.
	// ============================================================
	gaussianParser::gaussianParser():baseParser()
	{
		this->itsChkPtFile = "";
		
		this->bChkPt = 0;
		
		this->itsMem = "";
		
		this->bMem = 0;
		
		this->bNProc = 0;
		
		this->bWriteInternalCoords = true;
		this->bWriteCartCoords = false;
		
		this->itsTheory = "HF";
		
		this->itsBasisSet = "6-31G*";
		
		this->itsBasisSetFile = "";
		
		this->itsCharge = 0;
		
		this->itsMultiplicity = 1;
		
		this->itsVerbosity = "N";
		
		this->bZMatrixGenerated = false;
		
		this->bModRedundant = false;

		this->g09corhigh = false;
		
		pGaussEVectors = 0;
	}
	
	// ============================================================
	// Function : ~gaussianParser()
	// ------------------------------------------------------------
	// Destructor for the class
	// All data is destroyed.
	// ============================================================
	gaussianParser::~gaussianParser()
	{
		if (this->bZMatrixGenerated) {
			delete pZmatParser;
		}
	}
	
	// ============================================================
	// Function : Read
	// ------------------------------------------------------------
	// parsers a gaussian output file
	// ------------------------------------------------------------
	// Format:
	// ============================================================
	void gaussianParser::Read(const std::string &gaussfile, sheet* pSheet)
	{
		std::ifstream igau;
		igau.open(gaussfile.c_str());
		
		int end   = gaussfile.length();
		int slash = gaussfile.find_last_of("/");
		std::string file_name = gaussfile.substr(slash+1,(end-slash-5));
		
		if (!igau) {
			//std::cout << "\nUNABLE TO OPEN GAUSSIAN LOG FILE"
			//          << "\nFILENAME = " << gaussfile
			//          << "\nEXITING...\n" << std::endl;
			//exit(1); // make this return an integer...
			std::stringstream ss;
			ss << "\nUNABLE TO OPEN GAUSSIAN LOG FILE"
			<< "\nFILENAME = " << gaussfile
			<< "\nEXITING...\n" << std::endl;
			std::cout << ss.str();
			throw MTKException(ss.str()); // make this return an integer...
						      //return;
		}
		
		std::string fileline = "";
		std::vector<std::string> words;
		std::vector<std::string> words2;
		
		int nAtoms = 0;
		
		bool bReadZMatrix = false;
		int intIndex = 1;
		
		bool bGFreq = false;
		
		double scaleFreq = 0.8929;
		
		// - Read the file - //
		while (igau) {
			getline(igau,fileline);
			
			if (containsSubStr(fileline, "ModRedundant")) {
				bModRedundant = true;
			}

			if (containsSubStr(fileline, " Redundant internal coordinates found in file.")) {
				g09corhigh = true;
			}

			if (fileline.substr(0,31) == " Redundant internal coordinates") {
				bModRedundant = true;
				if (g09corhigh == false) {
					getline(igau,fileline); // chk file
					getline(igau,fileline); // charge, multiplicity
				}
				getline(igau,fileline);
				while (fileline.substr(0,37) != " Recover connectivity data from disk.") {
					std::vector<std::string> v;
					splitString(fileline, ",", v, 0);
					atomLabels.push_back(v[0]);
					vector3d* pCoords = new vector3d(0.0);
					pCoords->setX(string2Double(v[2]));
					pCoords->setY(string2Double(v[3]));
					pCoords->setZ(string2Double(v[4]));
					molCoords.push_back(pCoords);
					nAtoms++;
					getline(igau,fileline);
				}
				
				int nEVs = 0;
				if (nAtoms == 2) {
					nEVs = 1;
				}
				else {
					nEVs = (3*nAtoms)-6;
				}
				//eigenvectors.resize((3*nAtoms)-6, (3*nAtoms));
				eigenvectors.resize(nEVs, (3*nAtoms));
				
				if (pSheet) {
					pGaussEValues = pSheet->addTable();
					if (!pGaussEValues) exit(1);
					pGaussEValues->setName("Gaussian EigenValues");
					
					pGaussEVectors = pSheet->addTable();
					if (!pGaussEVectors) exit(1);
					pGaussEVectors->setName("Gaussian EigenVectors");
					pGaussEVectors->setSizes((3*nAtoms)-6, (3*nAtoms));
					
					pGaussAtoms = pSheet->addStringTable();
					if (!pGaussAtoms) exit(1);
					pGaussAtoms->setName("Gaussian Atoms");
					
					pGaussCoords = pSheet->addTable();
					if (!pGaussCoords) exit(1);
					pGaussCoords->setName("Gaussian Coordinates");
					
					pGaussAtoms->setSizes(nAtoms, 1);
					
					pGaussEValues->setSizes(nEVs, 1);
					pGaussEVectors->setSizes(nEVs, (3*nAtoms));
					pGaussCoords->setSizes(1, (3*nAtoms));
					
					pGaussAtoms->initialize("");
					pGaussEValues->initialize(0.0);
					pGaussEVectors->initialize(0.0);
					pGaussCoords->initialize(0.0);
					
					for (int a = 0; a < nAtoms; a++) {
						pGaussAtoms->setCellValue(a, 0, atomLabels[a]);
						for (unsigned int a2 = 0; a2 < 3; a2++) {
							pGaussCoords->setCellValue(0, (a*3) + a2, (*molCoords[a])[a2]);
						}
					}
				}
			}
			if (bModRedundant) {
				if (containsSubStr(fileline, "calculate D2E/DX2 analytically")) {
					fileline = replaceCharacter(fileline, '!', ' ');
					splitString(fileline, " ", words, 0);
					std::string intType = words[1].substr(0,1);
					std::string intIndices = words[1].substr(2,words[1].size()-3);
					splitString(intIndices, ",", words2, 0);
					if (intType == "R") {
						int at1 = atoi(words2[0].c_str());
						int at2 = atoi(words2[1].c_str());
						//std::cout << intIndex << " BOND " << at1 << " " << at2 << std::endl;
						modRedBonds[indexAB(at1, at2, MAXATOMS)] = words[0];
						zmatDataType[words[0]] = 1;
					}
					else if (intType == "A") {
						int at1 = atoi(words2[0].c_str());
						int at2 = atoi(words2[1].c_str());
						int at3 = atoi(words2[2].c_str());
						//std::cout << intIndex <<  " ANGLE " << at1 << " " << at2 << " " << at3 << std::endl;
						
						try {
							ULONG_KIND angleIndex = indexABC_ULL(at1, at2, at3,
											     MAXATOMS,MAXATOMS);
							//modRedAngles[indexABC(at1, at2, at3, MAXATOMS, MAXATOMS)] = words[0];
							modRedAngles[angleIndex] = words[0];
						}
						catch (MTKException& e) {
							std::cout << "Error in gaussianParser: " << e.message << std::endl;
						}
						
						
						zmatDataType[words[0]] = 2;
					}
					else if (intType == "D") {
						int at1 = atoi(words2[0].c_str());
						int at2 = atoi(words2[1].c_str());
						int at3 = atoi(words2[2].c_str());
						int at4 = atoi(words2[3].c_str());
						//std::cout << intIndex << " DIHEDRAL " << at1 << " " << at2 << " " << at3 << " " << at4 << std::endl;
						modRedDihedrals[indexABCD(at1, at2, at3, at4, MAXATOMS, MAXATOMS, MAXATOMS)] = words[0];
						zmatDataType[words[0]] = 3;
					}
					zmatDataNames.push_back(words[0]);
					zmatData[words[0]] = strtod(words[2].c_str(),0);
					words.clear();
					words2.clear();
					intIndex++;
				}
			}
			
			if (fileline.substr(0,41) == " Z-Matrix taken from the checkpoint file:") {
				getline(igau,fileline);
				fileline = " Symbolic Z-matrix:";
			}
			if (fileline.substr(0,19) == " Symbolic Z-matrix:") {
				getline(igau,fileline);
				splitString(fileline, " ", words, 0);
				this->itsCharge = atoi(words[2].c_str());
				this->itsMultiplicity = atoi(words[5].c_str());
				words.clear();
				
				getline(igau,fileline);
				fileline = replaceCharacter(fileline, ',', ' ');
				splitString(fileline, " ", words, 0);
				if (words.size() == 4) continue;
				int nAt = 0;
				while (words[0] != "Variables:") {
					
					zmatrix.push_back(words);
					nAt++;
					nAtoms++;
					if (nAt > 1) {
						std::vector<int> bd;
						bd.push_back(nAt);
						bd.push_back(atoi(words[1].c_str()));
						zmatrixBonds.push_back(bd);
					}
					if (nAt > 2) {
						std::vector<int> ang;
						ang.push_back(nAt);
						ang.push_back(atoi(words[1].c_str()));
						ang.push_back(atoi(words[3].c_str()));
						zmatrixAngles.push_back(ang);
					}
					if (nAt > 3) {
						std::vector<int> tor;
						tor.push_back(nAt);
						tor.push_back(atoi(words[1].c_str()));
						tor.push_back(atoi(words[3].c_str()));
						tor.push_back(atoi(words[5].c_str()));
						zmatrixTorsions.push_back(tor);
					}
					
					if (words.size() > 7) {
						zmatDataType[words[2]] = 1;
						zmatDataNames.push_back(words[2]);
						zmatDataType[words[4]] = 2;
						zmatDataAngles.push_back(words[4]);
						zmatDataType[words[6]] = 3;
						zmatDataTorsions.push_back(words[6]);
					}
					else if (words.size() == 5) {
						zmatDataType[words[2]] = 1;
						zmatDataNames.push_back(words[2]);
						zmatDataType[words[4]] = 2;
						zmatDataAngles.push_back(words[4]);
					}
					else if (words.size() == 3) {
						zmatDataType[words[2]] = 1;
						zmatDataNames.push_back(words[2]);
					}
					
					words.clear();
					getline(igau,fileline);
					fileline = replaceCharacter(fileline, ',', ' ');
					splitString(fileline, " ", words, 0);
				}
				for (unsigned int x = 0; x < zmatDataAngles.size(); x++) {
					zmatDataNames.push_back(zmatDataAngles[x]);
				}
				for (unsigned int x = 0; x < zmatDataTorsions.size(); x++) {
					zmatDataNames.push_back(zmatDataTorsions[x]);
				}
				
				words.clear();
				getline(igau,fileline);
				fileline = replaceCharacter(fileline, '=', ' ');
				splitString(fileline, " ", words, 0);
				while (words.size() == 2) {
					zmatData[words[0]] = strtod(words[1].c_str(),0);
					words.clear();
					getline(igau,fileline);
					fileline = replaceCharacter(fileline, '=', ' ');
					splitString(fileline, " ", words, 0);
				}
				bReadZMatrix = true;
				
				int nEVs = 0;
				if (nAtoms == 2) {
					nEVs = 1;
				}
				else {
					nEVs = (3*nAtoms)-6;
				}
				//eigenvectors.resize((3*nAtoms)-6, (3*nAtoms));
				eigenvectors.resize(nEVs, (3*nAtoms));
				
				
				if (pSheet) {
					pGaussEValues = pSheet->addTable();
					if (!pGaussEValues) exit(1);
					pGaussEValues->setName("Gaussian EigenValues");
					
					pGaussEVectors = pSheet->addTable();
					if (!pGaussEVectors) exit(1);
					pGaussEVectors->setName("Gaussian EigenVectors");
					pGaussEVectors->setSizes((3*nAtoms)-6, (3*nAtoms));
					
					pGaussAtoms = pSheet->addStringTable();
					if (!pGaussAtoms) exit(1);
					pGaussAtoms->setName("Gaussian Atoms");
					
					//pGaussCoords = pSheet->addTable();
					//if (!pGaussCoords) exit(1);
					//pGaussCoords->setName("Gaussian Coordinates");
					
					pGaussAtoms->setSizes(nAtoms, 1);
					
					pGaussEValues->setSizes(nEVs, 1);
					pGaussEVectors->setSizes(nEVs, (3*nAtoms));
					
					//pGaussEValues->setSizes((3*nAtoms)-6, 1);
					//pGaussEVectors->setSizes((3*nAtoms)-6, (3*nAtoms));
					//pGaussCoords->setSizes(1, (3*nAtoms));
					
					pGaussAtoms->initialize("");
					pGaussEValues->initialize(0.0);
					pGaussEVectors->initialize(0.0);
					//pGaussCoords->initialize(0.0);
					
					for (int a = 0; a < nAtoms; a++) {
						pGaussAtoms->setCellValue(a, 0, zmatrix[a][0]);
						//for (unsigned int a2 = 0; a2 < 3; a2++) {
						//  pGaussCoords->setCellValue(0, (a*3) + a2, molCoords[a][a2]);
						//}
					}
					
				}
			}
			
			// Update zmatrix data with optimized values
			if (fileline.substr(0,51) == "                       !   Optimized Parameters   !") {
				// Skip next 4 lines
				for (int p = 0; p < 5; p++) getline(igau, fileline);
				while (fileline.substr(0,10) != " ---------") {
					words.clear();
					splitString(fileline, " ", words, 0);
					zmatData[words[1]] = strtod(words[2].c_str(),0);
					getline(igau,fileline);
				}
			}
			
			int counter = 0;
			unsigned int s = 0;
			double fc = 0.0;
			if (fileline.substr(0,26) == " Internal force constants:") {
				getline(igau,fileline); // skip line
				for (unsigned int i = 0; i < zmatData.size(); i++) {
					words.clear();
					getline(igau,fileline);
					splitString(fileline, " ", words, 0);
					s = words.size();
					words[s-1] = replaceCharacter(words[s-1], 'D', 'e');
					fc = atof(words[s-1].c_str());
					//if (fc > 1) {
					//fc /= 10.0;
					//std::cout << i+1 << " " << zmatDataNames[i] << " " << fc << std::endl;
					//}
					forceConstants[zmatDataNames[i]] = fc;
					counter++;
					if (s == 6) {
						int skip = zmatData.size() - counter;
						for (int sk = 0; sk < skip; sk++) {
							getline(igau,fileline);
						}
						getline(igau, fileline); // skip line
					}
				}
			}
			
			if (fileline.substr(0,16) == " Frequencies -- " and !bGFreq) {
				words.clear();
				splitString(fileline, " ", words, 0);
				s = words.size();
				for (unsigned int f = 2; f < s; f++) {
					frequencies.push_back(string2Double(words[f]));
					std::cout << words[f] << std::endl;
					
					int lFreq = static_cast<int>(frequencies.size());
					pGaussEValues->setCellValue(lFreq-1, 0, string2Double(words[f]) * scaleFreq);
				}
				int lastFreq = static_cast<int>(frequencies.size());
				
				// skip 4 lines
				getline(igau, fileline);
				getline(igau, fileline);
				getline(igau, fileline);
				getline(igau, fileline);
				for (int a = 0; a < nAtoms; a++) {
					getline(igau, fileline);
					words.clear();
					splitString(fileline, " ", words, 0);
					s = words.size();
					if (s == 11) {
						std::cout << words[2] << " " << words[3] << " " << words[4] << " ";
						std::cout << string2Double(words[2]) << " " << string2Double(words[3]) << " " << string2Double(words[4]) << std::endl;
						eigenvectors(lastFreq-3, a*3  ) = string2Double(words[2]);
						eigenvectors(lastFreq-3, a*3+1) = string2Double(words[3]);
						eigenvectors(lastFreq-3, a*3+2) = string2Double(words[4]);
						
						eigenvectors(lastFreq-2, a*3  ) = string2Double(words[5]);
						eigenvectors(lastFreq-2, a*3+1) = string2Double(words[6]);
						eigenvectors(lastFreq-2, a*3+2) = string2Double(words[7]);
						
						eigenvectors(lastFreq-1, a*3  ) = string2Double(words[8]);
						eigenvectors(lastFreq-1, a*3+1) = string2Double(words[9]);
						eigenvectors(lastFreq-1, a*3+2) = string2Double(words[10]);
						
						pGaussEVectors->setCellValue(lastFreq-3, a*3  , string2Double(words[2]));
						pGaussEVectors->setCellValue(lastFreq-3, a*3+1, string2Double(words[3]));
						pGaussEVectors->setCellValue(lastFreq-3, a*3+2, string2Double(words[4]));
						
						pGaussEVectors->setCellValue(lastFreq-2, a*3  , string2Double(words[5]));
						pGaussEVectors->setCellValue(lastFreq-2, a*3+1, string2Double(words[6]));
						pGaussEVectors->setCellValue(lastFreq-2, a*3+2, string2Double(words[7]));
						
						pGaussEVectors->setCellValue(lastFreq-1, a*3  , string2Double(words[8]));
						pGaussEVectors->setCellValue(lastFreq-1, a*3+1, string2Double(words[9]));
						pGaussEVectors->setCellValue(lastFreq-1, a*3+2, string2Double(words[10]));
					}
					else if (s == 8) {
						eigenvectors(lastFreq-2, a*3  ) = string2Double(words[2]);
						eigenvectors(lastFreq-2, a*3+1) = string2Double(words[3]);
						eigenvectors(lastFreq-2, a*3+2) = string2Double(words[4]);
						
						eigenvectors(lastFreq-1, a*3  ) = string2Double(words[5]);
						eigenvectors(lastFreq-1, a*3+1) = string2Double(words[6]);
						eigenvectors(lastFreq-1, a*3+2) = string2Double(words[7]);
						
						pGaussEVectors->setCellValue(lastFreq-2, a*3  , string2Double(words[2]));
						pGaussEVectors->setCellValue(lastFreq-2, a*3+1, string2Double(words[3]));
						pGaussEVectors->setCellValue(lastFreq-2, a*3+2, string2Double(words[4]));
						
						pGaussEVectors->setCellValue(lastFreq-1, a*3  , string2Double(words[5]));
						pGaussEVectors->setCellValue(lastFreq-1, a*3+1, string2Double(words[6]));
						pGaussEVectors->setCellValue(lastFreq-1, a*3+2, string2Double(words[7]));
					}
					else if (s == 5) {
						eigenvectors(lastFreq-1, a*3  ) = string2Double(words[2]);
						eigenvectors(lastFreq-1, a*3+1) = string2Double(words[3]);
						eigenvectors(lastFreq-1, a*3+2) = string2Double(words[4]);
						
						pGaussEVectors->setCellValue(lastFreq-1, a*3  , string2Double(words[2]));
						pGaussEVectors->setCellValue(lastFreq-1, a*3+1, string2Double(words[3]));
						pGaussEVectors->setCellValue(lastFreq-1, a*3+2, string2Double(words[4]));
					}
					else {
						std::cout << " ERROR " << std::endl;
					}
				}
			}
			
			if (containsSubStr(fileline, " Frequencies --- ")) {
				bGFreq = true;
				words.clear();
				splitString(fileline, " ", words, 0);
				s = words.size();
				int lastFreq = static_cast<int>(frequencies.size());
				for (unsigned int f = 2; f < s; f++) {
					pGaussEValues->setCellValue(lastFreq, 0, string2Double(words[f]) * scaleFreq);
					frequencies.push_back(string2Double(words[f]));
					//std::cout << words[f] << std::endl;
					lastFreq = static_cast<int>(frequencies.size());
				}
				
				// skip 4 lines
				getline(igau, fileline);
				getline(igau, fileline);
				getline(igau, fileline);
				getline(igau, fileline);
				for (int a = 0; a < nAtoms*3; a++) {
					getline(igau, fileline);
					words.clear();
					splitString(fileline, " ", words, 0);
					s = words.size();
					int nCurFreqs = s-3;
					for (unsigned int s2 = 3; s2 < s; s2++) {
						eigenvectors(lastFreq - nCurFreqs, a) = string2Double(words[s2]);
						pGaussEVectors->setCellValue(lastFreq - nCurFreqs, a, string2Double(words[s2]));
						nCurFreqs--;
					}
				}
			}
		}
		
#ifdef DEBUG
		if (!bModRedundant) {
			std::cout << " gaussianParser::zmatrix " << std::endl;
			for (unsigned int i = 0; i < zmatrix.size(); i++) {
				for (unsigned int j = 0; j < zmatrix[i].size(); j++) {
					std::cout << zmatrix[i][j] << " ";
				}
				std::cout << " " << std::endl;
			}
		}
		
		typedef std::map<std::string, double>::iterator zmatDataIterator;
		for (zmatDataIterator z = zmatData.begin(); z != zmatData.end(); z++) {
			std::cout << z->first << " " << std::showpoint << z->second << std::endl;
		}
		std::cout << " " << std::endl;
		
		int c = 0;
		typedef std::map<std::string, double>::iterator fcIterator;
		for (fcIterator z = forceConstants.begin(); z != forceConstants.end(); z++) {
			if (zmatDataType[z->first] == 1) {
				std::cout << c+1 << " " << z->first << " " << zmatDataType[z->first] << " " << std::showpoint
				<< (z->second*HB2TOKCALMOLA2)/2 << std::endl;
			}
			else {
				std::cout << c+1 << " " << z->first << " " << std::showpoint << (z->second*H2KCALMOL)/2 << std::endl;
			}
			
			c++;
		}
		std::cout << " " << std::endl;
#endif
		
	}
	
	// ============================================================
	// Function : Write
	// ------------------------------------------------------------
	// Write a gaussian input file.
	// ------------------------------------------------------------
	// Format:
	// ============================================================
	void gaussianParser::Write(const std::string &gaussfile, molecule* pMolecule,
				   std::vector< vector3d > &coordinates)
	{
		std::string routecard;
		std::ofstream ogauss;
		ogauss.open(gaussfile.c_str());
		
		if (pMolecule == 0) {
			std::cout << "\nERROR!"
			<< "\nMolecule is missing or empty."
			<< "\nGaussian file generation has been stopped." << std::endl;
			return;
		}
		
		if (!ogauss) {
			std::cout << "\nERROR!"
			<< "\nGaussian job file \"" << gaussfile << "\" could not be created or opened for writing."
			<< "\nGaussian file generation has been stopped." << std::endl;
			return;
		}
		
		ogauss << "$RunGauss" << std::endl;
		
		if (bChkPt) {
			ogauss << "%Chk=" << this->itsChkPtFile << std::endl;
		}
		if (bMem) {
			ogauss << "%Mem=" << this->itsMem << std::endl;
		}
		if (bNProc) {
			ogauss << "%NProcShared=" << this->itsNProc << std::endl;
		}
		
		routecard = "#" + this->itsVerbosity + " " + this->itsTheory + "/" + this->itsBasisSet;

		if (!itsCommandOptions.empty()) {
			for (mapIterator p = itsCommandOptions.begin(); p != itsCommandOptions.end(); p++) {
				std::string c =  p->first;
				std::vector<std::string> o = p->second;
				std::string opts = "";
				if (o.size() == 1) {
					if ( p->second[0] != "NONE") {
						opts = "=" + p->second[0];
					}
				}
				else {
					opts = "(";
					for (unsigned int i = 0; i < o.size(); i++) {
						opts = opts + o[i];
						if (i != o.size()-1) {
							opts = opts + ",";
						}
					}
					opts = opts + ")";
				}
				
				if (routecard.length() + c.length() + opts.length() < 72) {
					routecard = routecard + " " + c + opts;
				} else {
					if (routecard.length() != 0) {
						ogauss << routecard << std::endl;
					}
					routecard = c + opts;
				}
			}
		}
		for (unsigned int i = 0; i < this->iops.size(); i++) {
			if (routecard.length() + this->iops[i].length() < 72) {
				routecard = routecard + " " + this->iops[i];
			} else {
				if (routecard.length() != 0) {
					ogauss << routecard << std::endl;
				}
				routecard = this->iops[i];
			}
		}
		
		// Flush command block and follow with a blank line
		ogauss << routecard << "\n" << std::endl;
		
		// Write name (if required) and follow with a blank line
		if (bWriteMoleculeName) {
			ogauss << pMolecule->getName() << "\n"
			<< std::endl;
		}
		
		// Write charge and multiplicity (if required). Do not follow with a blank line.
		if (bWriteChargeAndMult) {
			ogauss << this->itsCharge << " " << this->itsMultiplicity << std::endl;
		}
		
		if (bWriteInternalCoords) {
			int failure = 0;
			if (!bZMatrixGenerated) {
				//zmatParser* pZmatParser = new zmatParser();
				//int failure = pZmatParser->genZmatrix(pMolecule, zmatrix, zmatData, atomMap);
				failure = this->generateZMatrix(pMolecule);
			}
			
			if (!failure) {
				for (unsigned int q = 0; q < zmatrix.size(); q++) {
					for (unsigned int w = 0; w < zmatrix[q].size(); w++) {
						ogauss << zmatrix[q][w] << " ";
					}
					ogauss << " " << std::endl;
				}
				ogauss << " " << std::endl;
				
				typedef std::map<std::string, double>::iterator zmatDataIterator;
				for (zmatDataIterator z = zmatData.begin(); z != zmatData.end(); z++) {
					ogauss << z->first << " " << std::showpoint << z->second << std::endl;
				}
				ogauss << " " << std::endl;
			}
			else {
				return;
			}
		}
		else if (bWriteCartCoords) {
			std::vector<atom*> atomList = pMolecule->getAtomList();
			for (unsigned int a = 0; a < atomList.size(); a++) {
				ogauss << atomList[a]->getElementSymbol() << " "
				<< std::showpoint << atomList[a]->getX() << " "
				<< std::showpoint << atomList[a]->getY() << " "
				<< std::showpoint << atomList[a]->getZ() << std::endl;
			}
			// End molecule specification
			ogauss << std::endl;
		}
		
		std::vector<std::string> optOptions = this->getCommandOption("Opt");
		if (optOptions.size() > 0) {
			bool bMD = 0;
			for (unsigned int x = 0; x < optOptions.size(); x++) {
				if (optOptions[x] == "ModRedundant") {
					bMD = 1;
				}
			}
			if (bMD == 1) {
				std::ifstream iMDFile;
				iMDFile.open(this->itsModRedundantFile.c_str());
				
				if (!iMDFile) {
					ogauss << std::endl;
				}
				else {
					std::string mdFileline = "";
					while (iMDFile){
						getline(iMDFile, mdFileline);
						ogauss << mdFileline << std::endl;
					}
					iMDFile.close();
				}
			}
		}
		
		if (this->itsBasisSet.substr(0,3) == "Gen") {
			std::ifstream iBSFile;
			iBSFile.open(this->itsBasisSetFile.c_str());
			
			if (!iBSFile) {
				std::cout << "\nERROR!"
				<< "\nBasis set \"Gen\" requires a custom file."
				<< "\nFile \"" << this->itsBasisSetFile << "\" is missing or unreadable."
				<< "\nPlease choose another file. For instructions, consult your program's"
				<< "\ndocumentation."
				<< "\nGaussian file generation has been stopped.\n" << std::endl;
				return;
			}
			
			std::string bsFileline = "";
			while (iBSFile){
				getline(iBSFile, bsFileline);
				ogauss << bsFileline << std::endl;
			}
			iBSFile.close();
		}
		
		// +BPR
		if (this->getCommandOption("Pseudo").size() > 0) {
			std::ifstream iPseudoFile;
			iPseudoFile.open(this->itsPseudoPotentialFile.c_str());
			
			if (!iPseudoFile) {
				std::cout << "\nERROR!"
				<< "\nPseudo=Read or Pseudo=Cards requires a custom file."
				<< "\nFile \"" << this->itsPseudoPotentialFile << "\" is missing or unreadable."
				<< "\nPlease choose another file. For instructions, consult your program's"
				<< "\ndocumentation."
				<< "\nGaussian file generation has been stopped.\n" << std::endl;
				return;
			}
			
			std::string pseudoFileline = "";
			while (iPseudoFile){
				getline(iPseudoFile, pseudoFileline);
				ogauss << pseudoFileline << std::endl;
			}
			iPseudoFile.close();
		}
		// -BPR
		
		if (this->getCommandOption("Pop").size() > 0) {
			if (!itsMKRadii.empty()) {
				for (dMapIterator p = itsMKRadii.begin(); p != itsMKRadii.end(); p++) {
					ogauss << p->first << "=" << p->second << "\n";
				}
				// End MK Radii spec.
				ogauss << " " << std::endl;
			}
		}
		
		ogauss.close();
	}
	
	// ============================================================
	// Function : Write
	// ------------------------------------------------------------
	// Write a gaussian input file.
	// ============================================================
	void gaussianParser::Write(const std::string &gaussfile, molecule* pMolecule)
	{
		if (pMolecule != 0) {
			std::vector< vector3d > coordinates;
			pMolecule->getCoordinates(coordinates);
			this->Write(gaussfile, pMolecule, coordinates);
		}
	}
	
	// ============================================================
	// Function : Write
	// ------------------------------------------------------------
	// Write a gaussian input file.
	// ============================================================
	void gaussianParser::Write(const std::string &gaussfile,collection* pCollection, const int &molId)
	{
		molecule* pMolecule = pCollection->getMolecule(molId);
		
		if (pMolecule != 0) {
			std::vector< vector3d > coordinates;
			pMolecule->getCoordinates(coordinates);
			this->Write(gaussfile, pMolecule, coordinates);
		}
	}
	
	// ============================================================
	// Function : generateZMatrix
	// ------------------------------------------------------------
	// Generate z-matrix
	// ============================================================
	int gaussianParser::generateZMatrix(molecule* pMolecule)
	{
		pZmatParser = new zmatParser();
		
		int failure = pZmatParser->genZmatrix(pMolecule, zmatrix, zmatData, atomMap);
		bZMatrixGenerated = true;
		return failure;
	}
	
	// ============================================================
	// Function : readMappingFile
	// ------------------------------------------------------------
	// Read mapping file name
	// ============================================================
	void gaussianParser::readMappingFile(std::string mapFile)
	{
		std::ifstream imap;
		imap.open(mapFile.c_str());
		
		if (!imap) {
			std::cout << "\nUNABLE TO OPEN MAP FILE"
			<< "\nFILENAME = " << mapFile << std::endl;
			return;
		}
		std::string fileline;
		
		while (imap) {
			getline(imap,fileline);
			if (fileline == "") break;
			std::vector<std::string> words;
			splitString(fileline, " ", words, 0);
			if (words.size() == 2) {
				int at1 = atoi(words[0].c_str());
				int at2 = atoi(words[1].c_str());
				atomMap[at1] = at2;
			}
			else {
				std::cout << fileline << std::endl;
				std::cout << " Error in gaussianParser::readMappingFile ... exiting " << std::endl;
				throw MTKException(" Error in gaussianParser::readMappingFile ... exiting ");
			}
		}
		imap.close();
	}
	
	// ============================================================
	// Function : writeMappingFile
	// ------------------------------------------------------------
	// Write mapping file name
	// ============================================================
	void gaussianParser::writeMappingFile(std::string mapFile)
	{
		if (!bZMatrixGenerated) {
			return;
		}
		
		std::ofstream omap;
		omap.open(mapFile.c_str());
		
		if (!omap) {
			std::cout << "\nUNABLE TO OPEN MAP FILE"
			<< "\nFILENAME = " << mapFile << std::endl;
			return;
		}
		for (iMapIterator ii = atomMap.begin(); ii != atomMap.end(); ii++) {
			omap << ii->first << " "  << ii->second << std::endl;
		}
		omap.close();
	}
	
	// ============================================================
	// Function : addCommandOption
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::addCommandOption(const std::string &c)
	{
		std::vector<std::string> popOption;
		popOption.push_back("NONE");
		this->itsCommandOptions[c] = popOption;
	}
	
	// ============================================================
	// Function : addCommandOption
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::addCommandOption(const std::string &c, const std::string &o)
	{
		std::vector<std::string> popOption;
		popOption.push_back(o);
		this->itsCommandOptions[c] = popOption;
	}
	
	// ============================================================
	// Function : addCommandOption
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::addCommandOption(const std::string &c, std::vector<std::string> &o)
	{
		this->itsCommandOptions[c] = o;
	}
	
	// ============================================================
	// Function : getCommandOption
	// ------------------------------------------------------------
	//
	// ============================================================
	std::vector<std::string> gaussianParser::getCommandOption(const std::string &c)
	{
		mapIterator p = itsCommandOptions.find(c);
		
		if (p != itsCommandOptions.end()){
			return p->second;
		}
		std::vector<std::string> blank;
		return blank;
	}
	
	// ============================================================
	// Function : removeCommandOption
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::removeCommandOption(const std::string &c)
	{
		mapIterator p = itsCommandOptions.find(c);
		
		if (p != itsCommandOptions.end()) {
			itsCommandOptions.erase(p);
		}
	}
	
	// ============================================================
	// Function : addIop
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::addIop(const std::string &i)
	{
		this->iops.push_back(i);
	}
	
	// ============================================================
	// Function : clearIop
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::clearIop()
	{
		this->iops.clear();
	}
	
	// ============================================================
	// Function : setChkPt
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setChkPt(const std::string &c)
	{
		this->bChkPt = 1;
		this->itsChkPtFile = c;
	}
	
	// ============================================================
	// Function : setMem
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setMem(const std::string &m)
	{
		this->bMem = 1;
		this->itsMem = m;
	}
	
	// ============================================================
	// Function : setNProc
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setNProc(const std::string &n)
	{
		this->bNProc = 1;
		this->itsNProc = n;
	}
	
	// ============================================================
	// Function : setTheory
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setTheory(const std::string &t)
	{
		this->itsTheory = t;
	}
	
	// ============================================================
	// Function : setBasisSet
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setBasisSet(const std::string &b)
	{
		this->itsBasisSet = b;
	}
	
	// ============================================================
	// Function : setBasisSetFile
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setBasisSetFile(const std::string &b)
	{
		this->itsBasisSetFile = b;
	}
	
	// ============================================================
	// Function : setPseudoPotentialFile
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setPseudoPotentialFile(const std::string &b)
	{
		this->itsPseudoPotentialFile = b;
	}
	
	// ============================================================
	// Function : setModRedundantFile
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setModRedundantFile(const std::string &b)
	{
		this->itsModRedundantFile = b;
	}
	
	// ============================================================
	// Function : setNoCoords
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setNoCoords()
	{
		this->bWriteInternalCoords = false;
		this->bWriteCartCoords = false;
	}
	
	// ============================================================
	// Function : setWriteMoleculeName
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setWriteMoleculeName(int i)
	{
		if (i) {
			this->bWriteMoleculeName = true;
		} else {
			this->bWriteMoleculeName = false;
		}
	}
	
	// ============================================================
	// Function : setWriteMoleculeName
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setWriteChargeAndMult(int i)
	{
		if (i) {
			this->bWriteChargeAndMult = true;
		} else {
			this->bWriteChargeAndMult = false;
		}
	}
	
	// ============================================================
	// Function : setVerbosity
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setVerbosity(const std::string &i)
	{
		this->itsVerbosity = i;
	}
	
	// ============================================================
	// Function : setCharge
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setCharge(const int &i)
	{
		this->itsCharge = i;
	}
	
	// ============================================================
	// Function : setMultiplicity
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setMultiplicity(const int &i)
	{
		this->itsMultiplicity = i;
	}
	
	// ============================================================
	// Function : setCartesian
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setCartesian(int i)
	{
		if (i) {
			this->bWriteCartCoords = true;
			this->bWriteInternalCoords = false;
		}
		else {
			this->bWriteCartCoords = false;
			this->bWriteInternalCoords = true;
		}
	}
	
	// ============================================================
	// Function : setInternal
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setInternal(int i)
	{
		if (i) {
			this->bWriteInternalCoords = true;
			this->bWriteCartCoords = false;
		}
		else {
			this->bWriteInternalCoords = false;
			this->bWriteCartCoords = true;
		}
	}
	
	// ============================================================
	// Function : setMKRadii
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::setMKRadii(std::string e, double d)
	{
		this->itsMKRadii[e] = d;
	}
	
	// ============================================================
	// Function : readFormattedChkPtFile
	// ------------------------------------------------------------
	//
	// ============================================================
	void gaussianParser::readFormattedChkPtFile(std::string fchkPtFile, sheet* pSheet)
	{
		
		std::ifstream ifchk;
		ifchk.open(fchkPtFile.c_str());
		
		if (!ifchk) {
			std::stringstream ss;
			ss << "\nUNABLE TO OPEN GAUSSIAN FCHK FILE"
			<< "\nFILENAME = " << fchkPtFile
			<< "\nEXITING...\n" << std::endl;
			std::cout << ss.str();
			throw MTKException(ss.str());
		}
		
		std::string fileline = "";
		
		// 1st line is name
		getline(ifchk,fileline);
		pSheet->setName(fileline);
		int nAtoms = 0;
		int nLowerTriangle = 0;
		
		//int indexI = 0;
		//int indexJ = 0;
		
		// - Read the file - //
		while (ifchk) {
			getline(ifchk,fileline);
			std::vector<std::string> words;
			
			if (fileline.substr(0,15) == "Number of atoms") {
				splitString(fileline, " ", words, 0);
				nAtoms = atoi(words[4].c_str());
			}
			else if (fileline.substr(0,29) == "Current cartesian coordinates") {
				table<double>* cartesianCoords = pSheet->addTable();
				cartesianCoords->setName("Current cartesian coordinates");
				cartesianCoords->setSizes(nAtoms,3);
				
				//ublas::matrix<double> &coords = cartesianCoords->getMatrix();
				Eigen::Matrix<double, Dynamic, Dynamic> &coords = cartesianCoords->getMatrix();
				
				std::vector<std::string> coord;
				int counter = 0;
				for (int ii = 0; ii < nAtoms; ii++) {
					for (int iii = 0; iii < 3; iii++) {
						if (counter == 5) {
							counter = 0;
						}
						if (counter == 0) {
							coord.clear();
							getline(ifchk,fileline);
							splitString(fileline, " ", coord, 0);
						}
						coords(ii, iii) = strtod(coord[counter].c_str(), 0) / BOHR2ANG;
						counter++;
					}
				}
			}
			else if (fileline.substr(0,25) == "Cartesian Force Constants") {
				splitString(fileline, " ", words, 0);
				table<double>* cartForceConstants = pSheet->addTable();
				cartForceConstants->setName("Cartesian Force Constants");
				cartForceConstants->setSizes(nAtoms*3,nAtoms*3);
				
				//ublas::matrix<double> &cFCs = cartForceConstants->getMatrix();
				Eigen::Matrix<double, Dynamic, Dynamic> &cFCs = cartForceConstants->getMatrix();
				
				nLowerTriangle = atoi(words[words.size()-1].c_str());
				//indexI = 0;
				//indexJ = 0;
				int counter = 0;
				
				std::vector<std::string> fc;
				for (int ii = 0; ii < nAtoms*3; ii++) {
					for (int iii = 0; iii < ii+1; iii++) {
						if (counter == 5) {
							counter = 0;
						}
						if (counter == 0) {
							fc.clear();
							getline(ifchk,fileline);
							splitString(fileline, " ", fc, 0);
						}
						cFCs(ii, iii) = strtod(fc[counter].c_str(), 0);
						cFCs(iii, ii) = cFCs(ii, iii);
						counter++;
					}
				}
			}
		}
	}
	
	// ============================================================
	// Function : getForceConstant
	// ------------------------------------------------------------
	// Seminario. Calculation of intramolecular force fields from second-derivative tensors.
	// Int. J. Quantum Chem (1996) vol. 60 (7) pp. 59-65
	// ============================================================
	int gaussianParser::getForceConstant(sheet* pSheet, int atomA, int atomB, double& d, double& fc)
	{
		table<double>* cartForceConsts = pSheet->getTable("Cartesian Force Constants");
		if (!cartForceConsts) {
			std::stringstream ss;
			ss << " Can't find Cartesian Force Constants in formatted checkpoint file ... exiting " << std::endl;
			std::cout << ss.str();
			throw MTKException(ss.str());
		}
		
		//ublas::matrix<double> &cartForceConstsMatrix = cartForceConsts->getMatrix();
		Eigen::Matrix<double, Dynamic, Dynamic> &cartForceConstsMatrix = cartForceConsts->getMatrix();
		
		table<double>* coords = pSheet->getTable("Current cartesian coordinates");
		if (!coords) {
			std::stringstream ss;
			ss << " Can't find Current cartesian coordinates in formatted checkpoint file ... exiting " << std::endl;
			std::cout << ss.str();
			throw MTKException(ss.str());
		}
		
		//ublas::matrix<double> &coordsMatrix = coords->getMatrix();
		//unsigned int coordMatrixRows = coordsMatrix.size1();
		//unsigned int cartForceConstsMatrixRows = cartForceConstsMatrix.size1();
		
		Eigen::Matrix<double, Dynamic, Dynamic> &coordsMatrix = coords->getMatrix();
		unsigned int coordMatrixRows = coordsMatrix.rows();
		unsigned int cartForceConstsMatrixRows = cartForceConstsMatrix.rows();
		
		if (static_cast<unsigned int>(atomA) > coordMatrixRows or
		    static_cast<unsigned int>(atomB) > coordMatrixRows) {
			return 0;
		}
		if (static_cast<unsigned int>(atomA*3) > cartForceConstsMatrixRows or
		    static_cast<unsigned int>(atomB*3) > cartForceConstsMatrixRows) {
			return 0;
		}
		
		
		vector3d vA;
		vector3d vB;
		for (int ati = 0; ati < 3; ati++) {
			vA[ati] = coordsMatrix(atomA,ati);// * BOHR2ANG;
			vB[ati] = coordsMatrix(atomB,ati);// * BOHR2ANG;
		}
		
		d = vA.dist(vB);
		
		vector3d Uab1 = vA - vB;
		vector3d Uab = Uab1.unit();
		
		int atomAtemp = atomB * 3;
		int atomBtemp = atomA * 3;
		
		//ublas::matrix<double, ublas::column_major> AB(3,3);
		//ublas::vector<double> eigenValuesAB(3);
		//for (int e = 0; e < 3; e++) {
		//  eigenValuesAB(e) = 0.0;
		//}
		
		Eigen::Matrix<double, Dynamic, Dynamic> AB(3,3);
		
		for (int ati = 0; ati < 3; ati++) {
			for (int atj = 0; atj < 3; atj++) {
				AB(ati, atj) = - cartForceConstsMatrix(atomAtemp, atomBtemp+atj);
			}
			atomAtemp++;
		}
		
		SelfAdjointEigenSolver<MatrixXd> eigensolver(AB);
		MatrixXd evectors = eigensolver.eigenvectors();
		VectorXd eigenValuesAB = eigensolver.eigenvalues();
		eigenValueSort(evectors, eigenValuesAB, 1);
		
		double fc1 = 0.0;
		double fc2 = 0.0;
		
		for (int e = 0; e < 3; e++) {
			vector3d ev(evectors(0,e), evectors(1,e), evectors(2,e));
			//std::cout << " EIGENVALUE " << eigenValuesAB[e] << " " << ev * Uab << std::endl;
			
			if (std::abs(ev * Uab) > 0.95) {
				fc1 = eigenValuesAB[e];
				//std::cout << " a) Dist A-B: " << vA.dist(vB) * ANG2BOHR << " " << ev * Uab << std::endl;
			}
			fc2 += (eigenValuesAB[e] * std::abs(Uab * ev));
		}
		
		/*
		 // boost/eigen
		 
		 int result = diagonalize(AB, eigenValuesAB);
		 if (result != 0) {
		 std::stringstream ss;
		 ss << " Error in diagonalization ... exiting " << std::endl;
		 std::cout << ss.str();
		 throw MTKException(ss.str());
		 }
		 
		 eigenValueSort(AB, eigenValuesAB, 1);
		 //std::cout << " EIGENVALUES: " << eigenValuesAB[0] << " " << eigenValuesAB[1] << " "
		 //          << eigenValuesAB[2] << std::endl;
		 
		 for (int e = 0; e < 3; e++) {
		 vector3d ev(AB(0,e), AB(1,e), AB(2,e));
		 //std::cout << " EIGENVALUE " << eigenValuesAB[e] << " " << ev * Uab << std::endl;
		 
		 if (std::abs(ev * Uab) > 0.95) {
		 fc1 = eigenValuesAB[e];
		 //std::cout << " a) Dist A-B: " << vA.dist(vB) * ANG2BOHR << " " << ev * Uab << std::endl;
		 }
		 fc2 += (eigenValuesAB[e] * std::abs(Uab * ev));
		 }
		 */
		
		//std::cout << " fc1 = " << fc1 << " fc2 = " << fc2 << " Dist = " << vA.dist(vB) / BOHR2ANG;
		//std::cout << " b) Dist A-B: " << vA.dist(vB) * ANG2BOHR << std::endl;
		if (fc1 > fc2) {
			fc = fc1;
		}
		else {
			fc = fc2;
		}
		fc = fc * HB2TOKCALMOLA2;
		//std::cout << " --> FORCE CONSTANT = " << fc << std::endl;
		return 0;
	}
	
	// ============================================================
	// Function : getForceConstant
	// ------------------------------------------------------------
	// Seminario. Calculation of intramolecular force fields from second-derivative tensors.
	// Int. J. Quantum Chem (1996) vol. 60 (7) pp. 59-65
	// ============================================================
	int gaussianParser::getForceConstant(sheet* pSheet, int atomA, int atomB, int atomC, double& a, double& fc)
	{
		table<double>* cartForceConsts = pSheet->getTable("Cartesian Force Constants");
		//ublas::matrix<double> &cartForceConstsMatrix = cartForceConsts->getMatrix();
		Eigen::Matrix<double, Dynamic, Dynamic> &cartForceConstsMatrix = cartForceConsts->getMatrix();
		
		table<double>* coords = pSheet->getTable("Current cartesian coordinates");
		//ublas::matrix<double> &coordsMatrix = coords->getMatrix();
		Eigen::Matrix<double, Dynamic, Dynamic> &coordsMatrix = coords->getMatrix();
		
		//unsigned int coordMatrixRows = coordsMatrix.size1();
		unsigned int coordMatrixRows = coordsMatrix.rows();
		
		if (static_cast<unsigned int>(atomA) > coordMatrixRows or
		    static_cast<unsigned int>(atomB) > coordMatrixRows or
		    static_cast<unsigned int>(atomC) > coordMatrixRows) {
			return 0;
		}
		
		//unsigned int cartForceConstsMatrixRows = cartForceConstsMatrix.size1();
		unsigned int cartForceConstsMatrixRows = cartForceConstsMatrix.rows();
		
		if (static_cast<unsigned int>(atomA*3) > cartForceConstsMatrixRows or
		    static_cast<unsigned int>(atomB*3) > cartForceConstsMatrixRows or
		    static_cast<unsigned int>(atomC*3) > cartForceConstsMatrixRows) {
			return 0;
		}
		
		vector3d vA;
		vector3d vB;
		vector3d vC;
		//ublas::matrix<double, ublas::column_major> eigenvectorsAB(3,3);
		//ublas::matrix<double, ublas::column_major> eigenvectorsCB(3,3);
		
		Eigen::Matrix<double, Dynamic, Dynamic> AB(3,3);
		Eigen::Matrix<double, Dynamic, Dynamic> CB(3,3);
		
		//ublas::vector<double> eigenValuesAB(3);
		//for (int e = 0; e < 3; e++) {
		//  eigenValuesAB(e) = 0.0;
		//}
		
		//ublas::vector<double> eigenValuesCB(3);
		//for (int e = 0; e < 3; e++) {
		//  eigenValuesCB(e) = 0.0;
		//}
		
		for (int ati = 0; ati < 3; ati++) {
			vA[ati] = coordsMatrix(atomA,ati) * BOHR2ANG;
			vB[ati] = coordsMatrix(atomB,ati) * BOHR2ANG;
			vC[ati] = coordsMatrix(atomC,ati) * BOHR2ANG;
		}
		
		a = angle(vA, vB, vC) * RAD2DEG;
		
		vector3d Uab1 = vA - vB;
		vector3d Uab = Uab1.unit();
		
		vector3d Ucb1 = vC - vB;
		vector3d Ucb = Ucb1.unit();
		
		vector3d u_N1 = cross(Ucb, Uab);
		vector3d u_N = u_N1.unit();
		
		vector3d u_PA = cross(u_N, Uab);
		vector3d u_PC = cross(Ucb, u_N);
		
		//std::cout << " Dist1: " << vA.dist(vB) * ANG2BOHR << std::endl;
		//std::cout << " Dist2: " << vC.dist(vB) * ANG2BOHR << std::endl;
		
		double dist2AB = vA.dist2(vB);
		double dist2CB = vC.dist2(vB);
		
		int atomAtemp = atomB * 3;
		int atomBtemp = atomA * 3;
		
		for (int ati = 0; ati < 3; ati++) {
			for (int atj = 0; atj < 3; atj++) {
				AB(ati, atj) = - cartForceConstsMatrix(atomAtemp, atomBtemp+atj);
			}
			atomAtemp++;
		}
		
		SelfAdjointEigenSolver<MatrixXd> eigensolver(AB);
		MatrixXd eigenvectorsAB = eigensolver.eigenvectors();
		VectorXd eigenValuesAB = eigensolver.eigenvalues();
		eigenValueSort(eigenvectorsAB, eigenValuesAB, 1);
		
		/*
		 int result = diagonalize(eigenvectorsAB, eigenValuesAB);
		 if (result != 0) {
		 std::stringstream ss;
		 ss << " Error in diagonalization ... exiting " << std::endl;
		 std::cout << ss.str();
		 throw MTKException(ss.str());
		 //exit(0);
		 }
		 eigenValueSort(eigenvectorsAB, eigenValuesAB, 1);
		 */
		
		atomAtemp = atomB * 3;
		atomBtemp = atomC * 3;
		
		for (int ati = 0; ati < 3; ati++) {
			for (int atj = 0; atj < 3; atj++) {
				CB(ati, atj) = - cartForceConstsMatrix(atomAtemp, atomBtemp+atj);
			}
			atomAtemp++;
		}
		
		SelfAdjointEigenSolver<MatrixXd> eigensolverCB(CB);
		MatrixXd eigenvectorsCB = eigensolverCB.eigenvectors();
		VectorXd eigenValuesCB = eigensolverCB.eigenvalues();
		eigenValueSort(eigenvectorsCB, eigenValuesCB, 1);
		
		/*
		 result = diagonalize(eigenvectorsCB, eigenValuesCB);
		 if (result != 0) {
		 std::stringstream ss;
		 ss << " Error in diagonalization ... exiting " << std::endl;
		 std::cout << ss.str();
		 throw MTKException(ss.str());
		 //exit(0);
		 }
		 eigenValueSort(eigenvectorsCB, eigenValuesCB, 1);
		 */
		double abContribution = 0.0;
		for (int e = 0; e < 3; e++) {
			vector3d ev(eigenvectorsAB(0,e), eigenvectorsAB(1,e), eigenvectorsAB(2,e));
			abContribution += (eigenValuesAB[e] * std::abs(u_PA * ev));
		}
		abContribution = 1.0 / (abContribution * dist2AB);
		
		double cbContribution = 0.0;
		for (int e = 0; e < 3; e++) {
			vector3d ev(eigenvectorsCB(0,e), eigenvectorsCB(1,e), eigenvectorsCB(2,e));
			cbContribution += (eigenValuesCB[e] * std::abs(u_PC * ev));
		}
		cbContribution = 1.0 / (cbContribution * dist2CB);
		
		fc = (1.0 / (abContribution + cbContribution)) * H2KCALMOL;
		//std::cout << " angle " << angle(vA, vB, vC) * RAD2DEG
		//          << " angleForceConstant = " << fc << std::endl;
		return 0;
	}
	
	// ============================================================
	// Function : getForceConstantZMAT
	// ------------------------------------------------------------
	//
	// ============================================================
	int gaussianParser::getForceConstantZMAT(int atomA, int atomB, double& r, double& fc)
	{
		fc = 0.0;
		std::string bdName = "";
		std::string eMess = "atom a:" + int2String(atomA) + " atom b:" + int2String(atomB) + "\n";
		
		if (this->bModRedundant) {
			bdName = modRedBonds[indexAB(atomA, atomB, MAXATOMS)];
		}
		else {
			if (atomMap[atomA] > atomMap[atomB]) {
				bdName = zmatrix[atomMap[atomA]-1][2];
			}
			else {
				bdName = zmatrix[atomMap[atomB]-1][2];
			}
		}
		if (bdName != "") {
			//std::cout << " " << zmatData[bdName] << " " << (forceConstants[bdName] * HB2TOKCALMOLA2) << std::endl;
			eMess = " Bond: " + bdName + " req = " + double2String(zmatData[bdName]) + " keq = " + double2String(forceConstants[bdName] * HB2TOKCALMOLA2);
			r = zmatData[bdName];
			fc = (forceConstants[bdName] * HB2TOKCALMOLA2);
		}
		else {
			int nAtoms = static_cast<int>(molCoords.size());
			if (atomA > nAtoms or atomB > nAtoms) {
				r = 0.0;
				fc = -1.0;
				eMess += "nAtoms:" + int2String(nAtoms ) + " index out of range";
				errorLogger.throwError("gaussianParser::getForceConstantZMAT", eMess, MTK_ERROR);
				return 1;
			}
			else {
				eMess = " Bond not defined, setting k to 0 ";
				//std::cout << " BOND NOT DEFINED " << std::endl;
				r = molCoords[atomA-1]->dist(*(molCoords[atomB-1]));
				fc = 0.0;
				errorLogger.throwError("gaussianParser::getForceConstantZMAT", eMess, INFO);
				return 1;
			}
		}
		errorLogger.throwError("gaussianParser::getForceConstantZMAT", eMess, INFO);
		return 0;
	}
	
	// ============================================================
	// Function : getForceConstantZMAT
	// ------------------------------------------------------------
	//
	// ============================================================
	int gaussianParser::getForceConstantZMAT(int atomA, int atomB, int atomC, double& t, double& fc)
	{
		fc = 0.0;
		std::string agName = "";
		std::string eMess = "atom a:" + int2String(atomA) + " atom b:" + int2String(atomB) 
		+ " atom c:" + int2String(atomC) + "\n";
		
		if (this->bModRedundant) {
			
			//agName = modRedAngles[indexABC(atomA, atomB, atomC, MAXATOMS, MAXATOMS)];
			
			try {
				ULONG_KIND angleIndex = indexABC_ULL(atomA, atomB, atomC,
								     MAXATOMS,MAXATOMS);
				agName = modRedAngles[angleIndex];
			}
			catch (MTKException& e) {
				std::cout << "Error in getForceConstantZMAT: " << e.message << std::endl;
			}
		}
		else {
			for (int i = 0; i < static_cast<int>(zmatrix.size()); i++) {
				if (static_cast<int>(zmatrix[i].size()) < 5) continue;
				if (atomMap[atomA] == i+1) {
					if ((atoi(zmatrix[i][1].c_str()) == atomMap[atomB]) and
					    (atoi(zmatrix[i][3].c_str()) == atomMap[atomC])) {
						agName = zmatrix[i][4];
						break;
					}
				}
				else if (atomMap[atomC] == i+1) {
					if ((atoi(zmatrix[i][1].c_str()) == atomMap[atomB]) and
					    (atoi(zmatrix[i][3].c_str()) == atomMap[atomA])) {
						agName = zmatrix[i][4];
						break;
					}
				}
			}
		}
		if (agName != "") {
			//std::cout << " " << zmatData[agName] << " " << (forceConstants[agName] * H2KCALMOL) << std::endl;
			eMess = " Angle: " + agName + " req = " + double2String(zmatData[agName]) + " keq = " + double2String(forceConstants[agName] * H2KCALMOL);
			t = zmatData[agName];
			fc = (forceConstants[agName] * H2KCALMOL);
		}
		else {
			int nAtoms = static_cast<int>(molCoords.size());
			if (atomA > nAtoms or atomB > nAtoms or atomC > nAtoms) {
				t = 0.0;
				fc = -1.0;
				eMess += "nAtoms:" + int2String(nAtoms ) + " index out of range";
				errorLogger.throwError("gaussianParser::getForceConstantZMAT", eMess, MTK_ERROR);
				return 1;
			}
			else {
				eMess = " Angle not defined, setting k to 0 ";
				//std::cout << " ANGLE NOT DEFINED, SETTING K TO 0 " << std::endl;
				t = RAD2DEG * angle(*(molCoords[atomA-1]), *(molCoords[atomB-1]), *(molCoords[atomC-1]));
				fc = 0.0;
				errorLogger.throwError("gaussianParser::getForceConstantZMAT", eMess, INFO);
				return 1;
			}
		}
		errorLogger.throwError("gaussianParser::getForceConstantZMAT", eMess, INFO);
		return 0;
	}
	
	// ============================================================
	// Function : getFrequencies
	// ------------------------------------------------------------
	//
	// ============================================================
	std::vector<double> gaussianParser::getFrequencies()
	{
		return this->frequencies;
	}
	
} // MTKpp namespace
