/*!
 \file MCPB.cpp

 \brief Developes Metal Cluster Force Field Parameters.

 \author Martin B. Peters

 $Date: 2010/05/04 21:43:32 $
 $Revision: 1.1 $

 ----------------------------------------------------------------------------

 MTK++ - C++ package of modeling libraries.

 Copyright (C) 2005-2009  (see AUTHORS file for a list of contributors)

 ----------------------------------------------------------------------------
 */
#include "Utils/printHeader.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/selection.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/angle.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"
#include "Molecule/atomTyper.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/stdGroup.h"
#include "Molecule/stdFrag.h"
#include "Molecule/parameters.h"
#include "Molecule/protonate.h"
#include "Molecule/metalCenter.h"

// - UTILS
#include "Utils/vector3d.h"
#include "Utils/constants.h"
#include "Utils/diagonalize.h"
#include "Utils/printHeader.h"

// - STATS
#include "Statistics/sheet.h"
#include "Statistics/table.h"

// - AMBER
#include "MM/amber.h"

// - PARSERS
#include "Parsers/baseParser.h"
#include "Parsers/elementParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/molParser.h"
#include "Parsers/sdfParser.h"
#include "Parsers/gaussianParser.h"
#include "Parsers/paramParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/prepParser.h"
#include "Parsers/frcmodParser.h"
#include "Parsers/acParser.h"
#include "Parsers/amberParser.h"
#include "Parsers/inputParser.h"
#include "Parsers/mtkppParser.h"
#include "Parsers/dMParser.h"
#include "Parsers/commLineOptions.h"
#include "Parsers/StringManip.h"

// - Log
#include "Log/errorHandler.h"

// - Package name
#include "config.h"

#include <math.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <vector>

using namespace MTKpp;

std::string getEnvVar(std::string const& key)
{
    char const* val = getenv(key.c_str()); 
    return val == NULL ? std::string() : std::string(val);
}

void printFunctionList()
{
	std::string functionList = "\n MCPB Functions:\n";
	functionList += "      addBondAndAngleParameters\n";
	functionList += "      Add bond and angle parameters if missing\n";
	functionList += "      syntax: addBondAndAngleParameters /col/mol groupName\n";
	functionList += "\n";
	functionList += "      addFragment\n";
	functionList += "      Add Fragment to atom\n";
	functionList += "      syntax: addFragment 6MemRings/6CH bd /col/Mol//34 ag /col/Mol//27 tr /col/Mol//9 180.0\n";
	functionList += "\n";
	functionList += "      addHs\n";
	functionList += "      Add Hydrogens to the cluster cuCYM4\n";
	functionList += "      syntax: addHs" ;
	functionList += "      syntax: addHs /COL/MOL\n";
	functionList += "\n";
	functionList += "      addStdMol\n";
	functionList += "      Add standard molecule to lib file\n";
	functionList += "      syntax: addStdMol /COL/MOL Group\n";
	functionList += "\n";
	functionList += "      addToResidue\n";
	functionList += "      Add atoms in HOH-772 to ON2\n";
	functionList += "           bb_heavy  == backbone [ca, n, c, o]\n";
	functionList += "           bb  == backbone [ca, h, ha, n, nh, c, o]\n";
	functionList += "           bbb == backbone [ca, h, ha, n, nh, c, o, cb]\n";
	functionList += "      syntax: addToResidue /1FWJ/CLR/ON2-22 /1FWJ/7/HOH-772\n";
	functionList += "\n";
	functionList += "      appendResidue\n";
	functionList += "      Add CYS-12@CB of molecule to residue CY1 of cluster cuCYM4\n";
	functionList += "      syntax: appendResidue /1FEE/cuCYM4/CY1 /1FEE/1/CYS-12/.CB.\n";
	functionList += "\n";
	functionList += "      assignConnectivity\n";
	functionList += "      Assigns all bonds, angles, torsions and impropers\n";
	functionList += "      syntax: assignConnectivity\n";
	functionList += "\n";
	functionList += "      assignDisulfideBonds\n";
	functionList += "      Assigns all disulfide bonds (Needs to be carried out before atom typing)\n";
	functionList += "      syntax: assignDisulfideBonds\n";
	functionList += "\n";
	functionList += "      assignParameters\n";
	functionList += "      Assigns bond/angle/torsion/improper parameters\n";
	functionList += "      syntax: assignParameters /COL/MOL\n";
	functionList += "\n";
	functionList += "      assignStdFeatures\n";
	functionList += "      Assigns std features to molecule\n";
	functionList += "      syntax: assignStdFeatures /col/mol\n";
	functionList += "\n";
	functionList += "      atomType\n";
	functionList += "      Assigns atom types in the collection\n";
	functionList += "      syntax: atomType\n";
	functionList += "\n";
	functionList += "      basisSet\n";
	functionList += "      Set Gaussian Basis Set\n";
	functionList += "      syntax: basisSet 6-31G*\n";
	functionList += "      syntax: basisSet GEN bs.txt\n";
	functionList += "      syntax: basisSet GEN.6D.7F bs.txt\n";
	functionList += "\n";
	functionList += "      capResidue\n";
	functionList += "      Cap residue R1 with NME and ACE\n";
	functionList += "      syntax: capResidue /1FEE/cuCYM4/CY1/.N.. ACE\n";
	functionList += "\n";
	functionList += "      clusterCharge\n";
	functionList += "      Set Gaussian Charge\n";
	functionList += "      syntax: clusterCharge cuCYM4 -3\n";
	functionList += "\n";
	functionList += "      clusterSpin\n";
	functionList += "      Set Gaussian Spin\n";
	functionList += "      syntax: clusterSpin 0\n";
	functionList += "\n";
	functionList += "      copyAtomType\n";
	functionList += "      Copy atom type named NB in parm94 to atom type named NX in 1CA2\n";
	functionList += "      syntax: copyAtomType parm94/NB 1CA2/NX\n";
	functionList += "\n";
	functionList += "      copyStdResidue\n";
	functionList += "      Copy residue named HIS to a residue named HS1\n";
	functionList += "      syntax: copyStdResidue aminoAcids94/HIS myLib/HS1\n";
	functionList += "\n";
	functionList += "      createBond\n";
	functionList += "      Create bond\n";
	functionList += "      syntax: createBond /1FEE/cuCYM4//CU.. /1FEE/cuCYM4//.SG.\n";
	functionList += "\n";
	functionList += "      createMolecule\n";
	functionList += "      Create molecule named cuCYM4\n";
	functionList += "      syntax: createMolecule cuCYM4\n";
	functionList += "\n";
	functionList += "      createResidue";
	functionList += "      Create residue named HS1 in molecule named MOL\n";
	functionList += "      syntax: createResidue HS1 in MOL\n";
	functionList += "\n";
	functionList += "      createStdGroup\n";
	functionList += "      Create standard group in current directory\n";
	functionList += "      syntax: createStdGroup myLib\n";
	functionList += "\n";
	functionList += "      findMetalCenters\n";
	functionList += "      Find all metal centers in the collection\n";
	functionList += "      syntax: findMetalCenters\n";
	functionList += "\n";
	functionList += "      gaussianCharges\n";
	functionList += "      Generate a Gaussian input file for partial-charge computation\n";
	functionList += "      Optional: set the Gaussian job file name\n";
	functionList += "      syntax: gaussianCharges //cuCYM4 cuCYM4.com\n";
	functionList += "\n";
	functionList += "      gaussianMem\n";
	functionList += "      Set amount of memory requested for Gaussian\n";
	functionList += "      syntax: gaussianMem 3600MB\n";
	functionList += "\n";
	functionList += "      gaussianMoldenFormat\n";
	functionList += "      Request Gaussian log files formatted for viewing in Molden\n";
	functionList += "      (adds the GFInput and IOp(6/7=3) keywords)\n";
	functionList += "      syntax: gaussianMoldenFormat\n";
	functionList += "\n";
	functionList += "      gaussianNProc\n";
	functionList += "      Set number of processors requested for Gaussian\n";
	functionList += "      syntax: gaussianNProc 2\n";
	functionList += "\n";
	functionList += "      gaussianOptAndFC\n";
	functionList += "      Generate Gaussian input files for optimisation and force constants\n";
	functionList += "      Optional: set the Gaussian job file name\n";
	functionList += "      syntax: gaussianOptAndFC //cuCYM4 cuCYM4.com\n";
	functionList += "\n";
	functionList += "      gaussianVerbosity\n";
	functionList += "      Set the verbosity of Gaussian output ([T]erse, [N]ormal, [P]rolix)\n";
	functionList += "      syntax: gaussianVerbosity [T|N|P]\n";
	functionList += "\n";
	functionList += "      levelOfTheory\n";
	functionList += "      Set Gaussian Theory Level\n";
	functionList += "      syntax: levelOfTheory HF\n";
	functionList += "\n";
	functionList += "      listFragments\n";
	functionList += "      List available fragments in a particular library\n";
	functionList += "      syntax: listFragments terminal\n";
	functionList += "\n";
	functionList += "      loadLib\n";
	functionList += "      Loads AMBER library files into MCPB\n";
	functionList += "      syntax: loadLib ~/MTKpp/data/amino94.xml\n";
	functionList += "\n";
	functionList += "      loadParam\n";
	functionList += "      Loads AMBER Parameters into MCPB\n";
	functionList += "      syntax: loadParam ~/MTKpp/data/parm94.xml\n";
	functionList += "\n";
	functionList += "      modRedundant\n";
	functionList += "      Add modReduntant definitions to gaussian input file\n";
	functionList += "      syntax: modRedundant md.txt\n";
	functionList += "\n";
	functionList += "      nmodeMatch\n";
	functionList += "      Compare nMode and Gaussian Normal Modes\n";
	functionList += "      syntax: nmodeMatch filename.xml\n";
	functionList += "\n";
	functionList += "      optimizePolarHs\n";
	functionList += "      Optimize polar hydrogens in a collection\n";
	functionList += "      syntax: optimizePolarHs\n";
	functionList += "\n";
	functionList += "      print\n";
	functionList += "      Print to screen details of structure\n";
	functionList += "      syntax: print cuCYM4\n";
	functionList += "\n";
	functionList += "      printAtomTypes\n";
	functionList += "      Print Available atom types\n";
	functionList += "      syntax: printAtomTypes\n";
	functionList += "\n";
	functionList += "      pseudoPotentials\n";
	functionList += "      Add pseudo potential definitions to gaussian input file\n";
	functionList += "      syntax: pseudoPotentials filename.txt\n";
	functionList += "\n";
	functionList += "      quit\n";
	functionList += "      Exits program\n";
	functionList += "      syntax: quit\n";
	functionList += "\n";
	functionList += "      readFormattedChkPtFile\n";
	functionList += "      Read Formatted Checkpoint file\n";
	functionList += "      syntax: readFormattedChkPtFile file.fchk\n";
	functionList += "\n";
	functionList += "      readGaussianOutput\n";
	functionList += "      Read Gaussian Output\n";
	functionList += "      syntax: readGaussianOutput cuCYM4.log\n";
	functionList += "\n";
	functionList += "      readMolZmatMapping\n";
	functionList += "      Read Molecule <--> Z-Matrix mapping file\n";
	functionList += "      syntax: readMolZmatMapping file.map\n";
	functionList += "\n";
	functionList += "      readNMode\n";
	functionList += "      Read frequencies from nmode\n";
	functionList += "      syntax: readNMode /Col/mol file_nmd2.out\n";
	functionList += "\n";
	functionList += "      readNModeVectors\n";
	functionList += "      Read eigenvalues and eigenvectors from nmode vecs file\n";
	functionList += "      syntax: readNModeVectors /Col/mol file_nmd2.vecs file_nmd2.molden\n";
	functionList += "\n";
	functionList += "      readPdb\n";
	functionList += "      Reads a PDB file\n";
	functionList += "      syntax: readPdb 1FEE 1FEE.pdb\n";
	functionList += "\n";
	functionList += "      readRespCharges\n";
	functionList += "      Read RESP charges into molecule\n";
	functionList += "      syntax: readRespCharges /COL/MOL file.resp2.chg\n";
	functionList += "\n";
	functionList += "      readSdf\n";
	functionList += "      Read sd file\n";
	functionList += "      syntax: readSdf 1ABC CLR.sdf\n";
	functionList += "\n";
	functionList += "      renumber\n";
	functionList += "      Renumbers atoms and residues in the collection\n";
	functionList += "      syntax: renumnber\n";
	functionList += "\n";
	functionList += "      respgenAdditions\n";
	functionList += "       Add info to respgen files\n";
	functionList += "       syntax: respgenAdditions groupName fileName bb\n";
	functionList += "\n";
	functionList += "            bb Definitions:\n";
	functionList += "            - 0 No restraints\n";
	functionList += "            - 1 Heavy Atoms in Backbone (bb_heavy, [ca, n, c, o]) set to parm94 values\n";
	functionList += "            - 2 Atoms in Backbone (bb, [ca, h, ha, n, nh, c, o]) set to parm94 values\n";
	functionList += "            - 3 Atoms in Backbone plus CB (bbb, [ca, h, ha, n, nh, c, o, cb]) set to parm94 values\n";
	functionList += "\n";
	functionList += "      set\n";
	functionList += "      Set variable\n";
	functionList += "      syntax: set variable_name value\n";
	functionList += "\n";
	functionList += "      setAtomName\n";
	functionList += "      Set atom name\n";
	functionList += "      syntax: setAtomName /1CA2/znCLR/ACE-1/.CA. to .CH3\n";
	functionList += "\n";
	functionList += "      setAtomType\n";
	functionList += "      Set atom type named of HS1@NE2 to NX from 1CA2\n";
	functionList += "      syntax: setAtomType myLib/HS1/.NE2 1CA2/NX\n";
	functionList += "\n";
	functionList += "      setFormalCharge\n";
	functionList += "      Set Formal Charge on atom\n";
	functionList += "      syntax: setFormalCharge /1FEE/cuCYM4//CU.. 1\n";
	functionList += "\n";
	functionList += "      setLoggingLevel\n";
	functionList += "      Set the verbosity of error/warning/info messages\n";
	functionList += "      syntax: setLoggingLevel 1\n";
	functionList += "           Values:\n";
	functionList += "             1 - Error\n";
	functionList += "             2 - Warning\n";
	functionList += "             3 - Debug\n";
	functionList += "             4 - Info\n";
	functionList += "\n";
	functionList += "      setMaxFileID\n";
	functionList += "      Set the id of a new atom in the collection\n";
	functionList += "      syntax: setMaxFileID /1L6J/znCLR /1L6J/1\n";
	functionList += "\n";
	functionList += "      setMKRadii\n";
	functionList += "      Set Merz-Kollman radii for element\n";
	functionList += "      syntax: setMKRadii cu 0.91\n";
	functionList += "\n";
	functionList += "      setResidueName\n";
	functionList += "      Set residue name\n";
	functionList += "      syntax: setResidueName /1CA2/1/HIS-119 to HIE\n";
	functionList += "\n";
	functionList += "      source\n";
	functionList += "      Sources a global file\n";
	functionList += "      syntax: source file_name\n";
	functionList += "\n";
	functionList += "      updateForceConstants\n";
	functionList += "      Determine force constants and add them to parm file\n";
	functionList += "      syntax: updateForceConstants /COL/MOL Group X Y\n";
	functionList += "             X Values:\n";
	functionList += "              - 0 Do not update bonds and angle equilibrium values\n";
	functionList += "              - 1 Do update bonds and angle equilibrium values (req)\n";
	functionList += "             Y Values:\n";
	functionList += "              - 0 Seminario Method\n";
	functionList += "              - 1 Z-matrix Method\n";
	functionList += "\n";
	functionList += "      updateFrequencies\n";
	functionList += "      Scale frequencies\n";
	functionList += "      syntax: updateFrequencies /Col/mol 0.9806\n";
	functionList += "\n";
	functionList += "      updateRespCharges\n";
	functionList += "      Add RESP charge to lib file\n";
	functionList += "      syntax: updateRespCharges /COL/MOL Group\n";
	functionList += "\n";
	functionList += "      writeData\n";
	functionList += "      Write contents\n";
	functionList += "      syntax: writeData f.xml\n";
	functionList += "\n";
	functionList += "      writeFrcmodFile\n";
	functionList += "      Writes AMBER Parameters\n";
	functionList += "      syntax: writeFrcmodFile ZnCCCC.frcmod 1A5T\n";
	functionList += "\n";
	functionList += "      writeLeap\n";
	functionList += "      Write the metal center bonding info for leap\n";
	functionList += "      syntax: writeLeap name pdbFile\n";
	functionList += "\n";
	functionList += "      writeLib\n";
	functionList += "      Write standard library\n";
	functionList += "      syntax: writeLib groupName fileName.xml\n";
	functionList += "\n";
	functionList += "      writeMol\n";
	functionList += "      Write mol file\n";
	functionList += "      syntax: writeMol /1ABC/clr cuCYS.mol\n";
	functionList += "\n";
	functionList += "      writeParams\n";
	functionList += "      Write all new parameters\n";
	functionList += "      syntax: writeParams groupName fileName.xml\n";
	functionList += "\n";
	functionList += "      writePdb\n";
	functionList += "      Write pdb file\n";
	functionList += "      syntax: writePdb /1ABC/clr cuCYS.pdb\n";
	functionList += "\n";
	functionList += "      writePrepFile\n";
	functionList += "      Writes AMBER prep file\n";
	functionList += "      syntax: writePrepFile ZnCCCC.prep 1A5T\n";
	functionList += "\n";
	functionList += "      writePrmtop\n";
	functionList += "      Write prmtop and coordinate files for sander\n";
	functionList += "      syntax: writePrmtop inpcrd prmtop\n";
	functionList += "\n";
	functionList += "      writeSdf\n";
	functionList += "      Write sd file\n";
	functionList += "      syntax: writeSdf /1ABC/clr CLR.sdf\n";
	functionList += "\n";
	functionList += "      writeState\n";
	functionList += "      Write state xml file\n";
	functionList += "      syntax: writeState file.mtk\n";
	std::cout << functionList << std::endl;
}

/*!
 \brief Developes metal cluster force field parameters
 \param argc
 \param argv
 \return success
 */
int main (int argc, char **argv)
{
	std::string prog_name = "MCPB";
	std::vector<std::string> authors;
	std::string author = "Martin B. Peters";
	authors.push_back(author);

	// PARSE COMMAND LINE

	// 1. CREATE AN OBJECT
	commLineOptions *clo = new commLineOptions();

	// 2. SET PREFERENCES
	clo->noUsage();

	// 3. SET THE USAGE/HELP
	clo->addUsage( "  MCPB: Semi-automated tool for metalloprotein         " );
	clo->addUsage( "          parameterization                           \n" );
	clo->addUsage( "    usage:  MCPB [flags] [options]                   \n" );
	clo->addUsage( "  options:                                             " );
	clo->addUsage( "          -i script file                               " );
	clo->addUsage( "          -l log file                                  " );
	clo->addUsage( "    flags:                                             " );
	clo->addUsage( "          -h help                                      " );
	clo->addUsage( "          -f function list                             " );

	// 4. SET THE OPTION STRINGS/CHARACTERS
	clo->setOption("inFile", 'i');
	clo->setOption("logFile", 'l');
	clo->setFlag("funcs", 'f');

	// 5. PROVIDE THE COMMANDLINE
	clo->processCommandArgs( argc, argv );

	clo->usageOn();

	// 6. GET THE VALUES
	if ( clo->getFlag( "help" ) || clo->getFlag( 'h' ) ) {
		printHeader(std::cout, prog_name, authors);
		clo->printUsage();
		return 0;
	}

	if ( clo->getFlag( "funcs" ) || clo->getFlag( 'f' ) ) {
		printHeader(std::cout, prog_name, authors);
		printFunctionList();
		return 0;
	}

	std::string inputFile = "";
	std::string logFile = "MCPB.log";

	if ( clo->getValue( "i" ) != 0 ) {
		inputFile = std::string(clo->getValue("i"));
	}
	else if ( clo->getValue( "inFile" ) != 0 ) {
		inputFile =  std::string(clo->getValue("inFile"));
	}
	else {
		printHeader(std::cout, prog_name, authors);
		clo->printUsage();
		std::cout << " Please provide an MCPB file " << std::endl;
		return 0;
	}

	if ( clo->getValue("l") != 0) {
		logFile = std::string(clo->getValue("l"));
	}
	else if ( clo->getValue("logFile") != 0) {
		logFile =  std::string(clo->getValue("logFile"));
	}

	if (inputFile == logFile) {
		printHeader(std::cout, prog_name, authors);
		clo->printUsage();

		std::cout << " The input and log files are the same ... exiting " << std::endl;
		exit(1);
	}

	// 6. OPEN LOG FILE
	std::ofstream oLog;
	oLog.open(logFile.c_str());

	if (!oLog) {
		printHeader(std::cout, prog_name, authors);
		clo->printUsage();
		std::cout << "\nUNABLE TO OPEN LOG FILE"
		<< "\nFILENAME = " << logFile << std::endl;
		exit(1);
	}

	// Set errorLog stream to the log file
	MTKpp::errorLogger.setStream(&oLog);

	// Print MTK++ copyright message
	printHeader(oLog, prog_name, authors);

	// 8. DONE
	delete clo;

	// --------------------------------------------------------------------- //

	std::string errorMessage = "";

	std::vector<std::vector<std::string> > inputFileContents;

	int failure = readInputFile(inputFile, inputFileContents);

	if (failure) {
		printHeader(std::cout, prog_name, authors);
		clo->printUsage();

		MTKpp::errorLogger.throwError("MCPB", "Failed to open MCPB input file ... exiting", MTK_ERROR);
		exit(1);
	}

	std::string AMBERHOME = getenv("AMBERHOME");
	collection* pCollection = new collection();
	pCollection->addParameters();
	pCollection->addStdLibrary();

	molecule* pMolecule = 0;
	submolecule* pSubMolecule = 0;

	// READ ELEMENT DATA
	elementParser* pElementParser = new elementParser(pCollection->pElements);
	std::string elementXmlFile = AMBERHOME + "/dat/mtkpp/elements.xml";
	pElementParser->Read(elementXmlFile);
	/*
	 std::vector<std::string> bb_heavy;
	 bb_heavy.push_back(" CA ");
	 bb_heavy.push_back(" N  ");
	 bb_heavy.push_back(" C  ");
	 bb_heavy.push_back(" O  ");

	 std::vector<std::string> bb_heavy_types;
	 bb_heavy_types.push_back("CT");
	 bb_heavy_types.push_back("N");
	 bb_heavy_types.push_back("C");
	 bb_heavy_types.push_back("O");
	 */
	std::vector<std::string> bb_heavy;
	bb_heavy.push_back(" CA :CT");
	bb_heavy.push_back(" N  :N");
	bb_heavy.push_back(" C  :C");
	bb_heavy.push_back(" O  :O");
	/*
	 std::vector<std::string> bb;
	 bb.push_back(" CA ");
	 bb.push_back(" H  ");
	 bb.push_back(" HA ");
	 bb.push_back(" N  ");
	 bb.push_back(" HN ");
	 bb.push_back(" C  ");
	 bb.push_back(" O  ");

	 std::vector<std::string> bb_types;
	 bb_types.push_back("CT");
	 bb_types.push_back("H");
	 bb_types.push_back("H1");
	 bb_types.push_back("N");
	 bb_types.push_back("H");
	 bb_types.push_back("C");
	 bb_types.push_back("O");
	 */

	std::vector<std::string> bb;
	bb.push_back(" CA :CT");
	bb.push_back(" H  :H");
	bb.push_back(" HA :H1");
	bb.push_back(" N  :N");
	bb.push_back(" HN :H");
	bb.push_back(" C  :C");
	bb.push_back(" O  :O");
	/*
	 std::vector<std::string> bbb;
	 bbb.push_back(" CA ");
	 bbb.push_back(" H  ");
	 bbb.push_back(" HA ");
	 bbb.push_back(" N  ");
	 bbb.push_back(" HN ");
	 bbb.push_back(" C  ");
	 bbb.push_back(" O  ");
	 bbb.push_back(" CB ");

	 std::vector<std::string> bbb_types;
	 bbb_types.push_back("CT");
	 bbb_types.push_back("H");
	 bbb_types.push_back("H1");
	 bbb_types.push_back("N");
	 bbb_types.push_back("H");
	 bbb_types.push_back("C");
	 bbb_types.push_back("O");
	 bbb_types.push_back("CT");
	 */

	std::vector<std::string> bbb;
	bbb.push_back(" CA :CT");
	bbb.push_back(" H  :H");
	bbb.push_back(" HA :H1");
	bbb.push_back(" N  :N");
	bbb.push_back(" HN :H");
	bbb.push_back(" C  :C");
	bbb.push_back(" O  :O");
	bbb.push_back(" CB :CT");

	// CREATE PARSERS
	pdbParser* pPdbParser = 0;
	pPdbParser = new pdbParser();

	gaussianParser* pGParser = 0;
	pGParser = new gaussianParser();

	// Store some matrices
	sheet* pSheet = new sheet();
	pSheet->setName("MCPB Data");
	dMParser* pDMParser = new dMParser();

	typedef std::map<int, Bond*>::iterator BondMapIterator;
	typedef std::map<ULONG_KIND, Angle*>::iterator AngleMapIterator;
	typedef std::vector<atom*>::iterator atIterator;
	typedef std::map<ULONG_KIND, Torsion*>::iterator TorsionMapIterator;
	typedef std::map<int, Improper*>::iterator   ImproperMapIterator;

	std::map<std::string, std::string> variableMap;
	typedef std::vector<std::vector<std::string> >::iterator myIterator;
	bool sourceFound = true;
	while (sourceFound) {
		sourceFound = false;
		int myIndex = 0;
		int sourceIndex = 0;
		for (myIterator i = inputFileContents.begin(); i != inputFileContents.end(); i++) {
			std::vector<std::string> cur = *i;

			if (cur[0] == "source") {
				/*!
				 @ingroup MCPB_commands
				 \code
				 Function: source

				 Description: Source a global file

				 syntax: source file_name
				 \endcode
				 */
				sourceFound = true;
				if ((cur.size() == 2)) {
					sourceIndex = myIndex;
					if (fileExists(cur[1])) {
						std::vector<std::vector<std::string> > sourceInputFileContents;
						failure = readInputFile(cur[1], sourceInputFileContents);
						if (failure) {
							errorMessage =  " Error reading file " + cur[1];
							MTKpp::errorLogger.throwError("MCPB::source", errorMessage, 1);
							exit(1);
						}
						inputFileContents.insert(i+1, sourceInputFileContents.begin(),
									 sourceInputFileContents.end());
					}
					else {
						errorMessage =  " Error reading file " + cur[1];
						MTKpp::errorLogger.throwError("MCPB::source", errorMessage, 1);
						exit(1);
					}
				}
				else {
					errorMessage =  " Improper use of the source command ";
					MTKpp::errorLogger.throwError("MCPB::source", errorMessage, 1);
					exit(1);
				}
			}
			if (sourceFound) {
				inputFileContents[sourceIndex][0] = "# source_read";
				break;
			}
			myIndex++;
		}
	}

	for (unsigned int i = 0; i < inputFileContents.size(); i++) {
		if (inputFileContents[i][0] == "# source_read") continue;

		MTKpp::errorLogger.throwError("MCPB", inputFileContents[i][0], INFO);

		if (inputFileContents[i][0] == "quit") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: quit

			 Description: Exits program

			 syntax: quit
			 \endcode
			 */
			MTKpp::errorLogger.throwError("MCPB", " MCPB Exited Normally", INFO);
			exit(0);
		}

		else if (inputFileContents[i][0] == "setLoggingLevel") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: setLoggingLevel

			 Description: Set the verbosity of error/warning/info messages

			 syntax: setLoggingLevel 1

			 Values:
			 1 - Error
			 2 - Warning
			 3 - Debug
			 4 - Info
			 \endcode
			 */
			if ((inputFileContents[i].size() == 2)) {
				int l = atoi(inputFileContents[i][1].c_str());
				MTKpp::errorLogger.setLevel(l);
			}
		}

		else if (inputFileContents[i][0] == "set") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: set

			 Description: Set variable

			 syntax: set variable_name value
			 \endcode
			 */
			if ((inputFileContents[i].size() == 3)) {
                //std::cout << "SET COMMAND " << inputFileContents[i][2] << std::endl;
                if (inputFileContents[i][2].compare(0,1,"$") == 0) {   
                  std::string envVariable = "";
                  std::string tVar = inputFileContents[i][2];
                  envVariable = getEnvVar(tVar.substr(1,tVar.size()-1));
                  if (envVariable == "") {
                      std::cout << " Environment variable " << tVar.substr(1,tVar.size()-1)
                                << " is unset " << std::endl;
                    exit(1);
                  }
                  //std::cout << "Environment variable = " << envVariable << std::endl;
                  inputFileContents[i][2] = envVariable;
			    }

				variableMap[inputFileContents[i][1]] = inputFileContents[i][2];
				for (unsigned int j = i+1; j < inputFileContents.size(); j++) {
					bool bReplaced = false;
					for (unsigned int k = 0; k < inputFileContents[j].size(); k++) {
						if (containsSubStr(inputFileContents[j][k], inputFileContents[i][1])) {
							inputFileContents[j][k] = replaceSubStr(inputFileContents[j][k],
												inputFileContents[i][1], inputFileContents[i][2]);
							bReplaced = true;
						}
					}
				}
			}
		}

		else if (inputFileContents[i][0] == "loadParam") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: loadParam

			 Description: Loads AMBER Parameters into MCPB

			 syntax: loadParam ~/MTKpp/data/parm94.xml
			 \endcode
			 */
			if ((inputFileContents[i].size() != 2) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::loadParam", " Incorrect use of loadParam ", MTK_ERROR);
				exit(1);
			}
			else {
				paramParser* pParamParser = new paramParser(pCollection->getParameters());
				if (pParamParser) {
					pParamParser->Read(inputFileContents[i][1]);
				}
				delete pParamParser;
			}
		}

		else if (inputFileContents[i][0] == "loadLib") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: loadLib

			 Description: Loads AMBER library files into MCPB

			 syntax: loadLib ~/MTKpp/data/amino94.xml
			 \endcode
			 */
			if ((inputFileContents[i].size() != 2) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::loadLib", " Incorrect use of loadLib ... exiting", MTK_ERROR);
				exit(1);
			}
			else {
				stdLibrary* pStdLibrary = pCollection->getStdLibrary();
				if (pStdLibrary) {
					parameters* pParms = pCollection->getParameters();
					if (!pParms) {
						MTKpp::errorLogger.throwError("MCPB::loadLib",
									      " Please read in parameters first ... exiting", MTK_ERROR);
						exit(1);
					}
					stdLibParser* pStdLibParser = new stdLibParser(pCollection, pStdLibrary, pParms);
					if (pStdLibParser) {
						pStdLibParser->Read(inputFileContents[i][1]);
						delete pStdLibParser;
					}
				}
			}
		}

		else if (inputFileContents[i][0] == "printGroupCharge") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: printFrag

			 Description: Print fragment details

			 syntax: printGroupCharge Zn-CCCC

			 \endcode
			 */
			if ((inputFileContents[i].size() != 2) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::printGroupCharge", " Incorrect use of printGroupCharge ... exiting", MTK_ERROR);
				exit(1);
			}
			else {
				stdLibrary* pStdLibrary = pCollection->getStdLibrary();
				if (pStdLibrary) {
					stdGroup* pStdGroup =  pStdLibrary->getStdGroup(inputFileContents[i][1]);
					std::cout << " Group " << pStdGroup->getName() << "\n";
					if (pStdGroup) {
						molecule* pM = pStdGroup->getStdMolecule();
						if (pM) {
							double groupCharge = 0.0;
							std::vector<submolecule*> pSM = pM->getSubMoleculeList();
							//std::cout << "  Number of submolecules = " << pSM.size() << "\n";
							for (unsigned int x = 0; x < pSM.size(); x++) {
								stdFrag* pSF = pStdGroup->getStdFrag(pSM[x]->getName());
								if (pSF) {
									std::cout << " " << pSF->getSymbol() << "\n";
									double fCharge = 0.0;
									std::vector<atom*> pA = pSM[x]->getAtomList();
									for (unsigned int y = 0; y < pA.size(); y++) {
										stdAtom* pSA = pA[y]->getStdAtom();
										if (pSA) {
											fCharge += pSA->atmCharge;
											std::cout << "       |" << pSA->identity << "| " << pSA->atmCharge << "\n";
										}
										else {
											std::cout << " Error find std atom " << std::endl;
											exit(1);
										}
									}
									groupCharge += fCharge;
									std::cout << "     Total Fragment Charge " << pSF->getCharge() << "\n";
								}
							}
							/*std::vector<stdFrag*> pL = pStdGroup->getStdFragList();
							for (unsigned int x = 0; x < pL.size(); x++) {
								std::cout << "   Fragment " << pL[x]->getSymbol() << "\n";
								std::vector<stdAtom*> pL2 = pL[x]->getStdAtomList();
								double fCharge = 0.0;
								for (unsigned int y = 0; y < pL2.size(); y++) {
									fCharge += pL2[y]->atmCharge;
									std::cout << "       Atom " << pL2[y]->identity << " " << pL2[y]->atmCharge << "\n";
								}
								groupCharge += fCharge;
								//groupCharge += pL[x]->getCharge();
								std::cout << "   Fragment Charge " << pL[x]->getCharge() << "\n";
							}
							*/

							std::cout   << "  Total Charge = " << groupCharge << "\n" << std::endl;
						}
					}
				}
			}
		}

		else if (inputFileContents[i][0] == "printFrag") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: printFrag

			 Description: Print fragment details

			 syntax: printFrag Zn-CCCC CY1

			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::printFrag", " Incorrect use of printFrag ... exiting", MTK_ERROR);
				exit(1);
			}
			else {
				stdLibrary* pStdLibrary = pCollection->getStdLibrary();
				if (pStdLibrary) {
                  stdGroup* pStdGroup =  pStdLibrary->getStdGroup(inputFileContents[i][1]);

                  double groupCharge = 0.0;
                  std::vector<stdFrag*> pL = pStdGroup->getStdFragList();
                  for (unsigned int x = 0; x < pL.size(); x++) {
                    groupCharge += pL[x]->getCharge();

                  }
                  std::cout << " Group " << inputFileContents[i][1] << " charge = "
                            << groupCharge << std::endl;
                  
                  stdFrag* pStdFrag = pStdGroup->getStdFrag(inputFileContents[i][2]);
                  if (pStdFrag) {
                    std::cout << "   Fragment " << inputFileContents[i][2] << " charge = "
                              << pStdFrag->getCharge() << std::endl;
                  }
		        }
		   }
		}
	
		else if (inputFileContents[i][0] == "readPdb") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: readPdb

			 Description: Reads a PDB file

			 syntax: readPdb 1FEE 1FEE.pdb
			 \endcode
			 */
			if (inputFileContents[i].size() != 3) {
				MTKpp::errorLogger.throwError("MCPB::readPdb", " Incorrect use of readPdb ... exiting", MTK_ERROR);
				exit(1);
			}
			else {
				pPdbParser->Read(inputFileContents[i][2], pCollection);
				pCollection->setName(inputFileContents[i][1]);
			}
		}

		else if (inputFileContents[i][0] == "writeFrcmodFile") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writeFrcmodFile

			 Description: Writes AMBER Parameters

			 syntax: writeFrcmodFile ZnCCCC.frcmod 1A5T
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::writeFrcmodFile", " Incorrect use of writeFrcmodFile ", MTK_ERROR);
				exit(1);
			}
			else {
				frcmodParser* pFrcmodParser = new frcmodParser(pCollection->getParameters(), inputFileContents[i][2]);
				if (pFrcmodParser) {
					pFrcmodParser->Write(inputFileContents[i][1], inputFileContents[i][2]);
				}
				delete pFrcmodParser;
			}
		}

		else if (inputFileContents[i][0] == "writePrepFile") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writePrepFile

			 Description: Writes AMBER prep file

			 syntax: writePrepFile ZnCCCC.prep 1A5T
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::writePrepFile", " Incorrect use of writePrepFile ", MTK_ERROR);
				exit(1);
			}
			else {
				stdLibrary* pStdLibrary = pCollection->getStdLibrary();
				if (pStdLibrary) {
					parameters* pParms = pCollection->getParameters();
					if (!pParms) {
						MTKpp::errorLogger.throwError("MCPB::writePrepFile",
									      " Please read in parameters first ", MTK_ERROR);
						exit(1);
					}
					prepParser* pPrepParser = new prepParser();
					if (pPrepParser) {
						stdGroup* pStdGroup = pStdLibrary->getStdGroup(inputFileContents[i][2]);
						pPrepParser->Write(inputFileContents[i][1], pStdGroup);
						delete pPrepParser;
					}
				}
			}
		}

		else if (inputFileContents[i][0] == "printAtomTypes") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: printAtomTypes

			 Description: Print Available atom types

			 syntax: printAtomTypes
			 \endcode
			 */
			if (!pCollection) {
				MTKpp::errorLogger.throwError("MCPB::printAtomTypes", " Incorrect use of printAtomTypes ... exiting", MTK_ERROR);
				exit(1);
			}
			else {
				pCollection->getParameters()->printAtomTypes();
			}
		}

		else if (inputFileContents[i][0] == "assignDisulfideBonds") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function assignDisulfideBonds

			 Description: Assigns all disulfide bonds (Needs to be carried out before atomtyping)

			 syntax assignDisulfideBonds
			 \endcode
			 */
			if (pCollection) {
				connections* pConnections = new connections(pCollection);
				pConnections->assignDisulfideBonds();
				delete pConnections;
			}
		}

		else if (inputFileContents[i][0] == "atomType") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: atomType

			 Description: Assigns atom types in the collection

			 syntax: atomType
			 \endcode
			 */
			atomTyper* pAtomTyper = new atomTyper();

			if (inputFileContents[i].size() == 1) {
				pAtomTyper->atomTypeByLib(pCollection);
			}
			else if (inputFileContents[i].size() == 2) {
				selection* pSeln = new selection(pCollection);
				failure = pSeln->parse(inputFileContents[i][1]);
				if (failure) {
					exit(1);
				}

				molecule* pMol = pSeln->getMol();
				if (!pMol) {
					MTKpp::errorLogger.throwError("MCPB::atomType",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				else {
					pAtomTyper->atomTypeByLib(pMol);
				}
			}
			delete pAtomTyper;
		}

		else if (inputFileContents[i][0] == "atomType2") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: atomType ---- fix this in the future

			 Description: Assigns atom types in the collection

			 syntax: atomType
			 \endcode
			 */
			atomTyper* pAtomTyper = new atomTyper(0); // preceive histidines
			pAtomTyper->atomTypeByLib(pCollection);
			delete pAtomTyper;
		}

		else if (inputFileContents[i][0] == "optimizePolarHs") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: optimizePolarHs

			 Description: Optimize polar hydrogens in a collection

			 syntax: optimizePolarHs
			 \endcode
			 */
			protonate* pProt = new protonate(pCollection);
			pProt->optimizePolarHs();
			delete pProt;
		}

		else if (inputFileContents[i][0] == "assignConnectivity") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: assignConnectivity

			 Description: Assigns all bonds, angles, torsions and impropers

			 syntax: assignConnectivity
			 \endcode
			 */
			if (pCollection) {
				connections* pConnections = new connections(pCollection);
				pConnections->run();
				delete pConnections;
			}
		}

		else if (inputFileContents[i][0] == "assignParameters") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: assignParameters

			 Description: Assigns bond/angle/torsion/improper parameters

			 syntax: assignParameters /COL/MOL
			 \endcode
			 */
			connections* pConnections = new connections(pCollection);

			if (inputFileContents[i].size() == 1) {
				pConnections->assignStd();
			}
			else if (inputFileContents[i].size() == 2) {
				selection* pSeln = new selection(pCollection);
				failure = pSeln->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::assignParameters",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				molecule* pMol = pSeln->getMol();
				if (!pMol) {
					MTKpp::errorLogger.throwError("MCPB::assignParameters",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				else {
					pConnections->assignStd(pMol);
				}
			}
			delete pConnections;
		}

		else if (inputFileContents[i][0] == "addBondAndAngleParameters") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: addBondAndAngleParameters

			 Description: Add bond and angle parameters if missing

			 syntax: addBondAndAngleParameters /col/mol groupName
			 \endcode
			 */
			if ((!pCollection) or (inputFileContents[i].size() != 3)) {
				MTKpp::errorLogger.throwError("MCPB::addBondAndAngleParameters",
							      " Error in addBondAndAngleParameters ... exiting ", MTK_ERROR);
				exit(1);
			}

			selection* pSeln = new selection(pCollection);
			failure = pSeln->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::addBondAndAngleParameters",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			molecule* pMol = pSeln->getMol();
			if (!pMol) {
				MTKpp::errorLogger.throwError("MCPB::addBondAndAngleParameters",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}

			parameters* pParams = pCollection->getParameters();

			atom* pAtom = 0;
			std::vector<atom*> molAtoms = pMol->getAtomList();
			for (atIterator a = molAtoms.begin(); a != molAtoms.end(); a++) {
				pAtom = *a;
				if (!pAtom->getStdAtom()) {
					errorMessage = " Missing atom parameters for |" + pAtom->getName() + "| ... exiting";
					MTKpp::errorLogger.throwError("MCPB::addBondAndAngleParameters",
								      errorMessage, MTK_ERROR);
					exit(1);
				}
			}

			Bond* pBond = 0;
			std::map<int, Bond*> molBonds = pMol->getBondMap();

			if (!molBonds.empty()) {
				for (BondMapIterator b = molBonds.begin(); b != molBonds.end(); b++) {
					pBond = b->second;
					if (!pBond->pBondParam) {
						atom* pAt1 = pBond->atom1;
						atom* pAt2 = pBond->atom2;
						if (pAt1 and pAt2) {
							if (pAt1->getParent()->getName() != "CH3" and
							    pAt2->getParent()->getName() != "CH3") {
								stdAtom* pStdAtom1 = pAt1->getStdAtom();
								stdAtom* pStdAtom2 = pAt2->getStdAtom();
								if (pStdAtom1 and pStdAtom2) {
									if (!pParams->getBondParam(pStdAtom1->type, pStdAtom2->type)) {
										oLog << " Missing bond parameters for |"
										<< pBond->atom1->getName() << "|-|"
										<< pBond->atom2->getName() << "|  |"
										<< pStdAtom1->type << "|-|"
										<< pStdAtom2->type << "|"
										//<< pBond->atom1->getParent()->getName() << ":"
										//<< pBond->atom1->getParent()->getSubMolId() << "--"
										//<< pBond->atom2->getParent()->getName() << ":"
										//<< pBond->atom2->getParent()->getSubMolId()
										<< std::endl;
										bondParam* pBP = pParams->addBondParam();
										pBP->atomType1 = pStdAtom1->type;
										pBP->atomType2 = pStdAtom2->type;
										pBP->keq = 0.0;
										pBP->req = pBond->size;
										pBP->groupName = inputFileContents[i][2];
										pBP->optimize = "t";
									}
								}
							}
						}
					}
				}
			}

			Angle* pAngle = 0;
			std::map<ULONG_KIND, Angle*> molAngles = pMol->getAngleMap();

			if (!molAngles.empty()) {
				for (AngleMapIterator a = molAngles.begin(); a != molAngles.end(); a++) {
					pAngle = a->second;
					if (!pAngle->pAngleParam) {
						atom* pAt1 = pAngle->atom1;
						atom* pAt2 = pAngle->atom2;
						atom* pAt3 = pAngle->atom3;
						if (pAt1 and pAt2 and pAt3) {
							if (pAt1->getParent()->getName() != "CH3" and
							    pAt2->getParent()->getName() != "CH3" and
							    pAt3->getParent()->getName() != "CH3") {
								stdAtom* pStdAtom1 = pAt1->getStdAtom();
								stdAtom* pStdAtom2 = pAt2->getStdAtom();
								stdAtom* pStdAtom3 = pAt3->getStdAtom();
								if (pStdAtom1 and pStdAtom2 and pStdAtom3) {
									if (!pParams->getAngleParam(pStdAtom1->type, pStdAtom2->type, pStdAtom3->type)) {
										oLog << " Missing angle parameters for |"
										<< pAngle->atom1->getName() << "|-|"
										<< pAngle->atom2->getName() << "|-|"
										<< pAngle->atom3->getName() << "|   |"
										<< pStdAtom1->type << "|-|"
										<< pStdAtom2->type << "|-|"
										<< pStdAtom3->type << "|"
										<< std::endl;
										angleParam* pAP = pParams->addAngleParam();
										pAP->atomType1 = pStdAtom1->type;
										pAP->atomType2 = pStdAtom2->type;
										pAP->atomType3 = pStdAtom3->type;
										pAP->keq = 0.0;
										pAP->req = pAngle->size;
										pAP->groupName = inputFileContents[i][2];
										pAP->optimize = "t";
									}
								}
							}
						}
					}
				}
			}

			Torsion* pTorsion = 0;
			std::map<ULONG_KIND, Torsion*> molTorsions = pMol->getTorsionMap();

			if (!molTorsions.empty()) {
				for (TorsionMapIterator t = molTorsions.begin(); t != molTorsions.end(); t++) {
					pTorsion = t->second;
					if (!pTorsion->parametersAssigned) {
						atom* pAt1 = pTorsion->atom1;
						atom* pAt2 = pTorsion->atom2;
						atom* pAt3 = pTorsion->atom3;
						atom* pAt4 = pTorsion->atom4;

						if (pAt1 and pAt2 and pAt3 and pAt4) {
							if (pAt1->getParent()->getName() != "CH3" and
							    pAt2->getParent()->getName() != "CH3" and
							    pAt3->getParent()->getName() != "CH3" and
							    pAt4->getParent()->getName() != "CH3") {

								stdAtom* pStdAtom1 = pAt1->getStdAtom();
								stdAtom* pStdAtom2 = pAt2->getStdAtom();
								stdAtom* pStdAtom3 = pAt3->getStdAtom();
								stdAtom* pStdAtom4 = pAt4->getStdAtom();

								if (pStdAtom1 and pStdAtom2 and
								    pStdAtom3 and pStdAtom4) {

									if (pParams->getTorsionParamList(pStdAtom1->type, pStdAtom2->type,
													 pStdAtom3->type, pStdAtom4->type).size() == 0) {
										oLog << " Missing torsion parameters for |"
										<< pAt1->getName() << "|-|"
										<< pAt2->getName() << "|-|"
										<< pAt3->getName() << "|-|"
										<< pAt4->getName() << "|   |"
										<< pStdAtom1->type << "|-|"
										<< pStdAtom2->type << "|-|"
										<< pStdAtom3->type << "|-|"
										<< pStdAtom4->type << "|"
										<< std::endl;

										torsionParam* pTP = pParams->addTorsionParam();
										pTP->atomType1 = pStdAtom1->type;
										pTP->atomType2 = pStdAtom2->type;
										pTP->atomType3 = pStdAtom3->type;
										pTP->atomType4 = pStdAtom4->type;

										// Hoops et al --> torsions are zero
										pTP->Nt = 3.0;
										pTP->Vn = 0.0;
										pTP->gamma = 0.0;
										pTP->npth = 3;

										pTP->groupName = inputFileContents[i][2];
										//pTP->optimize = "f";
									}
								}
							}
						}
					}
				}
			}


		}

		else if (inputFileContents[i][0] == "createStdGroup") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: createStdGroup

			 Description: Create standard group in current directory

			 syntax: createStdGroup myLib
			 \endcode
			 */
			if (inputFileContents[i].size() != 2) {
				MTKpp::errorLogger.throwError("MCPB::createStdGroup",
							      " Incorrect use of createStdGroup ", MTK_ERROR);
				exit(1);
			}
			else {
				stdLibrary* pStdLib = pCollection->getStdLibrary();
				if (!pStdLib) {
					MTKpp::errorLogger.throwError("MCPB::createStdGroup",
								      " Incorrect use of createStdGroup ", MTK_ERROR);
					exit(1);
				}
				stdGroup* pStdGroup = pStdLib->addStdGroup();
				if (!pStdGroup) {
					MTKpp::errorLogger.throwError("MCPB::createStdGroup",
								      " Incorrect use of createStdGroup ", MTK_ERROR);
					exit(1);
				}
				pStdGroup->setName(inputFileContents[i][1]);
			}
		}

		else if (inputFileContents[i][0] == "copyStdResidue") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: copyStdResidue

			 Description: Copy residue named HIS to a residue named HS1

			 syntax: copyStdResidue aminoAcids94/HIS myLib/HS1
			 \endcode
			 */
			if (inputFileContents[i].size() != 3) {
				MTKpp::errorLogger.throwError("MCPB::copyStdResidue",
							      " Incorrect use of copyStdResidue ", MTK_ERROR);
				exit(1);
			}
			else {
				std::vector<std::string> words1;
				splitString(inputFileContents[i][1], "/", words1, 0);
				if (words1.size() != 2) {
					MTKpp::errorLogger.throwError("MCPB::copyStdResidue",
								      " Incorrect use of copyStdResidue ", MTK_ERROR);
					exit(1);
				}
				words1[1] = replaceCharacter(words1[1], '.', ' ');

				std::vector<std::string> words2;
				splitString(inputFileContents[i][2], "/", words2, 0);
				if (words2.size() != 2) {
					MTKpp::errorLogger.throwError("MCPB::copyStdResidue",
								      " Incorrect use of copyStdResidue ", MTK_ERROR);
					exit(1);
				}
				words2[1] = replaceCharacter(words2[1], '.', ' ');

				stdLibrary* pStdLib = pCollection->getStdLibrary();

				stdGroup* pStdGroup1 = pStdLib->getStdGroup(words1[0]);

				if (!pStdGroup1) {
					MTKpp::errorLogger.throwError("MCPB::copyStdResidue",
								      " Incorrect use of copyStdResidue ", MTK_ERROR);
					exit(1);
				}

				stdFrag* pStdFrag1 = pStdGroup1->getStdFrag(words1[1]);
				if (!pStdFrag1) {
					MTKpp::errorLogger.throwError("MCPB::copyStdResidue",
								      " Incorrect use of copyStdResidue ", MTK_ERROR);
					exit(1);
				}

				stdGroup* pStdGroup2 = pStdLib->getStdGroup(words2[0]);
				if (!pStdGroup2) {
					MTKpp::errorLogger.throwError("MCPB::copyStdResidue",
								      " Incorrect use of copyStdResidue ", MTK_ERROR);
					exit(1);
				}
				stdFrag* pStdFrag2 = pStdGroup2->addStdFrag(pStdFrag1);
				pStdFrag2->setSymbol(words2[1]);
				std::string frag8LCode = pStdFrag2->getCode();
				std::string new8LCode = "";

				if (frag8LCode.size() == 8) {
					new8LCode = frag8LCode.substr(0,5) + words2[1];
				}
				else {
					new8LCode = "     " + words2[1];
				}
				pStdFrag2->setCode(new8LCode);
			}
		}

		else if (inputFileContents[i][0] == "copyAtomType") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: copyAtomType

			 Description: Copy atom type named NB in parm94 to atom type named NX in 1CA2

			 syntax: copyAtomType parm94/NB 1CA2/NX
			 \endcode
			 */
			if (inputFileContents[i].size() != 3) {
				MTKpp::errorLogger.throwError("MCPB::copyAtomType",
							      " Incorrect use of copyAtomType ", MTK_ERROR);
				exit(1);
			}
			else {
				std::vector<std::string> words1;
				splitString(inputFileContents[i][1], "/", words1, 0);
				if (words1.size() != 2) {
					MTKpp::errorLogger.throwError("MCPB::copyAtomType",
								      " Incorrect use of copyAtomType ", MTK_ERROR);
					exit(1);
				}

				std::vector<std::string> words2;
				splitString(inputFileContents[i][2], "/", words2, 0);
				if (words2.size() != 2) {
					MTKpp::errorLogger.throwError("MCPB::copyAtomType",
								      " Incorrect use of copyAtomType ", MTK_ERROR);
					exit(1);
				}

				words1[1] = replaceCharacter(words1[1], '.', ' ');
				words2[1] = replaceCharacter(words2[1], '.', ' ');
				parameters* pParams = pCollection->getParameters();
				if (!pParams) {
					MTKpp::errorLogger.throwError("MCPB::copyAtomType",
								      " Can't find parameters ", MTK_ERROR);
					exit(1);
				}

				atomType* pAtomType1 = pParams->getAtomType(words1[1]);

				if (pAtomType1) {
					if (pAtomType1->groupName != words1[0]) {
						errorMessage = pAtomType1->groupName + " != " + words1[0];
						MTKpp::errorLogger.throwError("MCPB::copyAtomType",
									      errorMessage, MTK_ERROR);
						exit(1);
					}
					//atomType* pNewAtomType = pParams->addAtomType(pAtomType1, words2[1], words2[0]);
					pParams->addAtomType(pAtomType1, words2[1], words2[0]);
					/*
					 if (!pNewAtomType) {
					 MTKpp::errorLogger.throwError("MCPB::copyAtomType",
					 " Incorrect use of copyAtomType ", MTK_ERROR);
					 exit(1);
					 }
					 */
				}
				else {
					errorMessage = " Error in copyAtomType can't find " + words1[1] + " in library " + words1[0];
					MTKpp::errorLogger.throwError("MCPB::copyAtomType",
								      errorMessage, MTK_ERROR);
					exit(1);
				}
			}
		}

		else if (inputFileContents[i][0] == "setAtomType") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: setAtomType

			 Description: Set atom type named of HS1@NE2 to NX from 1CA2

			 syntax: setAtomType myLib/HS1/.NE2 1CA2/NX
			 \endcode
			 */
			if (inputFileContents[i].size() != 3) {
				MTKpp::errorLogger.throwError("MCPB::setAtomType",
							      " Incorrect use of setAtomType ", MTK_ERROR);
				exit(1);
			}
			else {
				std::vector<std::string> words1;
				splitString(inputFileContents[i][1], "/", words1, 0);
				if (words1.size() != 3) {
					MTKpp::errorLogger.throwError("MCPB::setAtomType",
								      " Incorrect use of setAtomType ", MTK_ERROR);
					exit(1);
				}
				words1[2] = replaceCharacter(words1[2], '.', ' ');

				if (words1[2].size() != 4) {
					MTKpp::errorLogger.throwError("MCPB::setAtomType",
								      " Error in setAtomType -> atom name must be of size 4 ... exiting ", MTK_ERROR);
					exit(1);
				}

				stdLibrary* pStdLib = pCollection->getStdLibrary();

				stdGroup* pStdGroup = pStdLib->getStdGroup(words1[0]);
				if (!pStdGroup) {
					errorMessage = " Incorrect use of setAtomType: Can't find stdGroup: " + words1[0];
					MTKpp::errorLogger.throwError("MCPB::setAtomType",
								      errorMessage, MTK_ERROR);
					exit(1);
				}

				stdFrag* pStdFrag = pStdGroup->getStdFrag(words1[1]);
				if (!pStdFrag) {
					errorMessage = " Incorrect use of setAtomType: Can't find stdFrag: " + words1[1];
					MTKpp::errorLogger.throwError("MCPB::setAtomType",
								      errorMessage, MTK_ERROR);
					exit(1);
				}

				stdAtom* pStdAtom = pStdFrag->getStdAtom(words1[2]);
				if (!pStdAtom) {
					errorMessage = " Incorrect use of setAtomType: Can't find stdAtom: " + words1[2];
					MTKpp::errorLogger.throwError("MCPB::setAtomType",
								      errorMessage, MTK_ERROR);
					exit(1);
				}

				std::vector<std::string> words2;
				splitString(inputFileContents[i][2], "/", words2, 0);
				if (words2.size() != 2) {
					MTKpp::errorLogger.throwError("MCPB::setAtomType",
								      " Incorrect use of setAtomType ", MTK_ERROR);
					exit(1);
				}

				if (words2[1].size() != 2) {
					MTKpp::errorLogger.throwError("MCPB::setAtomType",
								      " Error in setAtomType -> atom type must be of size 2 ... exiting ", MTK_ERROR);
					exit(1);
				}
				pStdAtom->type = words2[1];
			}
		}

		else if (inputFileContents[i][0] == "createMolecule") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: createMolecule

			 Description: Create molecule named cuCYM4

			 syntax: createMolecule cuCYM4
			 \endcode
			 */
			if (inputFileContents[i].size() != 2) {
				MTKpp::errorLogger.throwError("MCPB::createMolecule",
							      " Incorrect use of createMolecule ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				pMolecule = pCollection->addMolecule();
				pMolecule->setName(inputFileContents[i][1]);
			}
		}

		else if (inputFileContents[i][0] == "setMaxFileID") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: setMaxFileID

			 Description:

			 syntax: setMaxFileID /1L6J/znCLR /1L6J/1
			 \endcode
			 */
			if ((!pCollection) or (inputFileContents[i].size() != 3)) {
				MTKpp::errorLogger.throwError("MCPB::setMaxFileID",
							      " Incorrect use of setMaxFileID ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				selection* pSeln1 = new selection(pCollection);
				failure = pSeln1->parse(inputFileContents[i][1]);
				errorMessage = " Error in selection parsing " +
				inputFileContents[i][1] + " ... exiting ";
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::setMaxFileID",
								      errorMessage, MTK_ERROR);
					exit(1);
				}
				molecule* pMol1 = pSeln1->getMol();
				if (!pMol1) {
					MTKpp::errorLogger.throwError("MCPB::setMaxFileID",
								      errorMessage, MTK_ERROR);
					exit(1);
				}

				selection* pSeln2 = new selection(pCollection);
				failure = pSeln2->parse(inputFileContents[i][2]);
				errorMessage = " Error in selection parsing " +
				inputFileContents[i][2] + " ... exiting ";
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::setMaxFileID",
								      errorMessage, MTK_ERROR);
					exit(1);
				}
				molecule* pMol2 = pSeln2->getMol();
				if (!pMol2) {
					MTKpp::errorLogger.throwError("MCPB::setMaxFileID",
								      errorMessage, MTK_ERROR);
					exit(1);
				}
				pMol1->setMaxFileID(pMol2->getMaxFileID()+1);
			}
		}

		else if (inputFileContents[i][0] == "createResidue") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: createResidue

			 Description: Create residue named HS1 in molecule named MOL

			 syntax: createResidue HS1 in MOL
			 \endcode
			 */
			if (inputFileContents[i].size() != 4) {
				MTKpp::errorLogger.throwError("MCPB::createResidue",
							      " Incorrect use of createResidue ", MTK_ERROR);
				exit(1);
			}
			else {
				if (inputFileContents[i][2] == "in") {
					pMolecule = pCollection->getMolecule(inputFileContents[i][3]);
					pSubMolecule = pMolecule->addSubMolecule();
					pSubMolecule->setName(inputFileContents[i][1]);
					pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());
				}
			}
		}

		else if (inputFileContents[i][0] == "addToResidue") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: addToResidue

			 Description: Add atoms in CYS-12 to CY1

			 bb_heavy  == backbone [ca, n, c, o]
			 bb  == backbone [ca, h, ha, n, nh, c, o]
			 bbb == backbone [ca, h, ha, n, nh, c, o, cb]

			 syntax: addToResidue ///CY1 //1/CYS-12
			 syntax: addToResidue ///CY1 //1/CYS-12 not bb
			 syntax: addToResidue ///GLY //1/PHE-10 and bb
			 ^            ^
			 |            +-- selection-expression
			 +-- residue name
			 \endcode
			 */
			std::vector<std::string>::iterator result;

			if ((!pCollection) or (inputFileContents[i].size() < 3)) {
				MTKpp::errorLogger.throwError("MCPB::addToResidue",
							      " Incorrect use of addToResidue ", MTK_ERROR);
				exit(1);
			}

			selection* pSeln1 = new selection(pCollection);
			failure = pSeln1->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::addToResidue",
							      " Error in selection parsing ", MTK_ERROR);
				exit(1);
			}

			submolecule* pSMol = pSeln1->getSMol();
			if (!pSMol) {
				MTKpp::errorLogger.throwError("MCPB::addToResidue",
							      " Error in selection parsing ", MTK_ERROR);
				exit(1);
			}

			selection* pSeln2 = new selection(pCollection);
			failure = pSeln2->parse(inputFileContents[i][2]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::addToResidue",
							      " Error in selection parsing ", MTK_ERROR);
				exit(1);
			}

			std::vector<atom*> selnAtoms = pSeln2->getAtoms();
			std::string dd = "";

			if (inputFileContents[i].size() == 5) {
				if (inputFileContents[i][3] == "not") {
					if (inputFileContents[i][4] == "bb") {
						for (unsigned int p = 0; p < selnAtoms.size(); p++) {
							stdAtom* pStdSelAtom = selnAtoms[p]->getStdAtom();
							dd = "";
							if (pStdSelAtom) {
								dd = pStdSelAtom->identity + ":" + pStdSelAtom->type;
								//result = std::find(bb.begin(), bb.end(), selnAtoms[p]->getName());
								result = std::find(bb.begin(), bb.end(), dd);
								if (result == bb.end()) {
									pSMol->addAtom(selnAtoms[p]);
								}
							}
						}
					}
					else if (inputFileContents[i][4] == "bbb") {
						for (unsigned int p = 0; p < selnAtoms.size(); p++) {
							stdAtom* pStdSelAtom = selnAtoms[p]->getStdAtom();
							dd = "";
							if (pStdSelAtom) {
								dd = pStdSelAtom->identity + ":" + pStdSelAtom->type;
								//result = std::find(bbb.begin(), bbb.end(), selnAtoms[p]->getName());
								result = std::find(bbb.begin(), bbb.end(), dd);
								if (result == bbb.end()) {
									pSMol->addAtom(selnAtoms[p]);
								}
							}
						}
					}
					else if (inputFileContents[i][4] == "bb_heavy") {
						for (unsigned int p = 0; p < selnAtoms.size(); p++) {
							stdAtom* pStdSelAtom = selnAtoms[p]->getStdAtom();
							dd = "";
							if (pStdSelAtom) {
								dd = pStdSelAtom->identity + ":" + pStdSelAtom->type;
								//result = std::find(bb_heavy.begin(), bb_heavy.end(), selnAtoms[p]->getName());
								result = std::find(bb_heavy.begin(), bb_heavy.end(), dd);
								if (result == bb_heavy.end()) {
									pSMol->addAtom(selnAtoms[p]);
								}
							}
						}
					}
				}
				else if (inputFileContents[i][3] == "and") {
					if (inputFileContents[i][4] == "bb") {
						for (unsigned int p = 0; p < selnAtoms.size(); p++) {
							stdAtom* pStdSelAtom = selnAtoms[p]->getStdAtom();
							dd = "";
							if (pStdSelAtom) {
								dd = pStdSelAtom->identity + ":" + pStdSelAtom->type;
								result = std::find(bb.begin(), bb.end(), dd);
								//result = std::find(bb.begin(), bb.end(), selnAtoms[p]->getName());
								if (result != bb.end()) {
									pSMol->addAtom(selnAtoms[p]);
								}
							}
						}
					}
					else if (inputFileContents[i][4] == "bbb") {
						for (unsigned int p = 0; p < selnAtoms.size(); p++) {
							stdAtom* pStdSelAtom = selnAtoms[p]->getStdAtom();
							dd = "";
							if (pStdSelAtom) {
								dd = pStdSelAtom->identity + ":" + pStdSelAtom->type;
								result = std::find(bbb.begin(), bbb.end(), dd);
								//result = std::find(bbb.begin(), bbb.end(), selnAtoms[p]->getName());
								if (result != bbb.end()) {
									pSMol->addAtom(selnAtoms[p]);
								}
							}
						}
					}
					else if (inputFileContents[i][4] == "bb_heavy") {
						for (unsigned int p = 0; p < selnAtoms.size(); p++) {
							stdAtom* pStdSelAtom = selnAtoms[p]->getStdAtom();
							dd = "";
							if (pStdSelAtom) {
								dd = pStdSelAtom->identity + ":" + pStdSelAtom->type;
								result = std::find(bb_heavy.begin(), bb_heavy.end(), dd);
								//result = std::find(bb_heavy.begin(), bb_heavy.end(), selnAtoms[p]->getName());
								if (result != bb_heavy.end()) {
									pSMol->addAtom(selnAtoms[p]);
								}
							}
						}
					}
				}
				else {
					MTKpp::errorLogger.throwError("MCPB::addToResidue",
								      " Unknown argument " + inputFileContents[i][3] + " ... exiting ", MTK_ERROR);
					exit(1);
				}
			}
			else {
				for (unsigned int p = 0; p < selnAtoms.size(); p++) {
					pSMol->addAtom(selnAtoms[p]);
					MTKpp::errorLogger.throwError("MCPB::addToResidue",
								      " Adding atom " + selnAtoms[p]->getName()
								      + " " + int2String(selnAtoms[p]->getFileID()), INFO);
				}
			}

			molecule* pMol1 = pSeln1->getMol();
			molecule* pMol2 = pSeln2->getMol();
			std::vector<atom*> clusterResAtoms = pMol1->getAtomList();
			if (clusterResAtoms.size() == 0) {
				MTKpp::errorLogger.throwError("MCPB::addToResidue",
							      " exiting ", MTK_ERROR);
				exit(1);
			}

			/*
			 std::cout << "MCPB::addToResidue Molecule name: " << pMol2->getName() << std::endl;

			 if (pMol2->getName() == "MOH") {
			 std::vector<atom*> mohAtoms = pMol2->getAtomList();
			 for (unsigned int m = 0; m < mohAtoms.size(); m++) {
			 std::cout << mohAtoms[m]->getFileID() << " " << mohAtoms[m]->getName() << std::endl;
			 }

			 if (mohAtoms.size() == 2) {
			 Bond* pBond1 = pMol2->getBond(mohAtoms[0], mohAtoms[1]);
			 if (pBond1) {
			 std::cout << " Bond exists " << std::endl;
			 }
			 }
			 }

			 */

			// Copy over bonds
			for (unsigned int a1 = 0; a1 < clusterResAtoms.size(); a1++) {
				atom* pAt1 = pMol2->getAtom(clusterResAtoms[a1]->getFileID(), 0, 1);
				if (!pAt1) {
					//std::cout << " Error finding atom " << clusterResAtoms[a1]->getFileID() << std::endl;
					continue;
				}

				for (unsigned int a2 = a1+1; a2 < clusterResAtoms.size(); a2++) {
					atom* pAt2 = pMol2->getAtom(clusterResAtoms[a2]->getFileID(), 0, 1);
					if (!pAt2) {
						//std::cout << " no atom " << clusterResAtoms[a2]->getFileID() << std::endl;
						continue;
					}

					if (pMol2->hasBond(pAt1, pAt2)) {
						Bond* pBond1 = pMol2->getBond(pAt1, pAt2);

						if (!pBond1) {
							//std::cout << " no bond " << pAt1->getFileID() << "-" << pAt2->getFileID() << std::endl;
							continue;
						}

						if (!pMol1->hasBond(clusterResAtoms[a1], clusterResAtoms[a2])) {
							Bond* pBond = pMol1->addBond(clusterResAtoms[a1], clusterResAtoms[a2],
										     pBond1->type, pBond1->stereo, pBond1->topology, 0.0);
							/*
							 MTKpp::errorLogger.throwError("MCPB::addToResidue",
							 " Adding bond " + pBond->atom1->getName()
							 + " " + pBond->atom2->getName(), INFO);
							 */
							if (pBond) {
								clusterResAtoms[a1]->addBondedAtom(clusterResAtoms[a2]);
								clusterResAtoms[a2]->addBondedAtom(clusterResAtoms[a1]);
							}
							else {
								std::cout << " MCPB::addToResidue Error adding bond " << std::endl;
							}
						}
						else {
							Bond* pBond = pMol1->getBond(clusterResAtoms[a1], clusterResAtoms[a2]);
							if (pBond) {
								MTKpp::errorLogger.throwError("MCPB::addToResidue",
											      " Bond already exists: " + pBond->atom1->getName()
											      + " " + pBond->atom2->getName(), INFO);
							}
						}
					}
				}
			}

			// Assign angles and torsions
			connections* pConnections = new connections(pCollection);
			pConnections->assignAngles(pMol1);
			pConnections->assignTorsions(pMol1);
			pConnections->assignImpropers(pMol1);
			delete pConnections;

			delete pSeln1;
			delete pSeln2;
		}

		else if (inputFileContents[i][0] == "appendResidue") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: appendResidue

			 Description: Add CYS-12@CB of molecule to residue CY1 of cluster cuCYM4

			 syntax: appendResidue /1FEE/cuCYM4/CY1 /1FEE/1/CYS-12/.CB.
			 \endcode
			 */
			if (inputFileContents[i].size() != 3) exit(1);

			selection* pSeln1 = new selection(pCollection);
			failure = pSeln1->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::appendResidue",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			molecule* pMol1 = pSeln1->getMol();
			if (!pMol1) {
				MTKpp::errorLogger.throwError("MCPB::appendResidue",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			submolecule* pSMol1 = pSeln1->getSMol();
			if (!pSMol1) {
				MTKpp::errorLogger.throwError("MCPB::appendResidue",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}

			selection* pSeln2 = new selection(pCollection);
			failure = pSeln2->parse(inputFileContents[i][2]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::appendResidue",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			molecule* pMol2 = pSeln2->getMol();
			if (!pMol2) {
				MTKpp::errorLogger.throwError("MCPB::appendResidue",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}

			std::vector<atom*> seln2Atoms = pSeln2->getAtoms();
			for (unsigned int p = 0; p < seln2Atoms.size(); p++) {
				pSMol1->addAtom(seln2Atoms[p]);
			}

			std::vector<atom*> clusterResAtoms = pSMol1->getAtomList();

			// Copy over bonds and assign angles and torsions
			for (unsigned int a1 = 0; a1 < clusterResAtoms.size(); a1++) {
				atom* pAt1 = pMol2->getAtom(clusterResAtoms[a1]->getFileID(), 0, 1);
				if (!pAt1) continue;
				for (unsigned int a2 = a1+1; a2 < clusterResAtoms.size(); a2++) {
					atom* pAt2 = pMol2->getAtom(clusterResAtoms[a2]->getFileID(), 0, 1);
					if (!pAt2) continue;
					if (pMol2->hasBond(pAt1, pAt2)) {
						Bond* pBond1 = pMol2->getBond(pAt1, pAt2);

						if (!pBond1) continue;
						if (!pMol1->hasBond(clusterResAtoms[a1], clusterResAtoms[a2])) {
							Bond* pBond = pMol1->addBond(clusterResAtoms[a1], clusterResAtoms[a2],
										     pBond1->type, pBond1->stereo, pBond1->topology, 0.0);
							if (pBond) {
								clusterResAtoms[a1]->addBondedAtom(clusterResAtoms[a2]);
								clusterResAtoms[a2]->addBondedAtom(clusterResAtoms[a1]);
							}
						}
					}
				}
			}
			connections* pConnections = new connections(pCollection);  // do this for just the added submolecule
			pConnections->assignAngles();
			pConnections->assignTorsions();
			delete pConnections;
			delete pSeln1;
			delete pSeln2;
		}

		else if (inputFileContents[i][0] == "addHs") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: addHs

			 Description: Add Hydrogens to the cluster cuCYM4

			 syntax: addHs
			 syntax: addHs /COL/MOL
			 \endcode
			 */
			std::string eMessage = "";

			if (inputFileContents[i].size() == 1) {
				eMessage += "Adding Hydrogens to collection\n";
				std::vector<molecule*> molList = pCollection->getMoleculeList();
				for (unsigned int m = 0; m < molList.size(); m++) {
					int nAtsBefore = molList[m]->getNumAtoms();
					if (molList[m]->getNumAtoms() < 100) {
						molList[m]->determineRings();
					}
					int maxID = 0;
					for (unsigned int m2 = 0; m2 < molList.size(); m2++) {
						if (maxID < molList[m2]->getMaxFileID()) {
							maxID = molList[m2]->getMaxFileID()+1;
						}
					}
					molList[m]->setMaxFileID(maxID);
					molList[m]->addHydrogens();
					int nAtsAfter = molList[m]->getNumAtoms();
					eMessage += "      " + molList[m]->getName() + ": " + int2String(nAtsAfter-nAtsBefore)
					+ " Hydrogen Atom/s added.\n";

				}
			}
			else if (inputFileContents[i].size() == 2) {
				selection* pSeln = new selection(pCollection);
				failure = pSeln->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::addHs",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				if (pSeln->getSelectionType() == 1) { // type: 0 == collection, 1 == molecule, 2 == submolecule
					molecule* pMol = pSeln->getMol();
					if (!pMol) {
						MTKpp::errorLogger.throwError("MCPB::addHs",
									      " Error in selection parsing ... exiting ", MTK_ERROR);
						exit(1);
					}
					else {
						if (pMol->getNumAtoms() < 100) {
							pMol->determineRings();
						}
						int nAtsBefore = pMol->getNumAtoms();
						pMol->addHydrogens();
						int nAtsAfter = pMol->getNumAtoms();
						eMessage += " Adding Hydrogens to molecule: " + pMol->getName() + ": "
						+ int2String(nAtsAfter-nAtsBefore) + " Hydrogen Atom/s added.\n";
					}
				}
				else if (pSeln->getSelectionType() == 2) {
					eMessage += "Adding Hydrogens to submolecule: \n";

				}
				else {
					eMessage += "Adding Hydrogens to collection \n";
					protonate* pProt = new protonate(pCollection);
					pProt->run();
					delete pProt;
				}
			}
			MTKpp::errorLogger.throwError("MCPB::addHs", eMessage, INFO);
		}

		else if (inputFileContents[i][0] == "setResidueName") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: setResidueName

			 Description: Set residue name

			 syntax: setResidueName /1CA2/1/HIS-119 to HIE
			 \endcode
			 */
			if (inputFileContents[i].size() != 4) {
				MTKpp::errorLogger.throwError("MCPB::setResidueName",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			selection* pSeln = new selection(pCollection);
			failure = pSeln->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::setResidueName",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			if (pSeln->getSelectionType() == 1) {
				molecule* pSelMol = pSeln->getMol();
				if (pSelMol) {
					pSelMol->setName(inputFileContents[i][3]);
				}
			}
			else if (pSeln->getSelectionType() == 2) {
				submolecule* pSelSubMol = pSeln->getSMol();
				if (std::string(inputFileContents[i][3]).size() == 3 and pSelSubMol) {
					pSelSubMol->setName(inputFileContents[i][3]);
				}
			}
		}

		else if (inputFileContents[i][0] == "setMoleculeName") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: setMoleculeName

			 Description: Set molecule name

			 syntax: setResidueName /1AMP//HOH-935 to MOH
			 \endcode
			 */
			if (inputFileContents[i].size() != 4) {
				MTKpp::errorLogger.throwError("MCPB::setMoleculeName",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			selection* pSeln = new selection(pCollection);
			failure = pSeln->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::setMoleculeName",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			if (pSeln->getSelectionType() == 1) {
				molecule* pSelMol = pSeln->getMol();
				if (pSelMol) {
					pSelMol->setName(inputFileContents[i][3]);
				}
			}
			else if (pSeln->getSelectionType() == 2) {
				submolecule* pSelSubMol = pSeln->getSMol();
				if (pSelSubMol) {
					molecule* pSelMol = pSelSubMol->getParent();
					if (pSelMol) {
						pSelMol->setName(inputFileContents[i][3]);
					}
				}
			}
		}

		else if (inputFileContents[i][0] == "setAtomName") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: setAtomName

			 Description: Set atom name

			 syntax: setAtomName /1CA2/znCLR/ACE-1/.CA. to .CH3
			 \endcode
			 */
			if (inputFileContents[i].size() != 4) {
				MTKpp::errorLogger.throwError("MCPB::setAtomName",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			selection* pSeln = new selection(pCollection);
			failure = pSeln->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::setAtomName",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			if (pSeln->getSelectionType() == 3) {
				atom* pAtom = pSeln->getAtom();
				if ((std::string(inputFileContents[i][3]).size() == 4) && (pAtom)) {
					std::string newName = replaceCharacter(inputFileContents[i][3], '.', ' ');
					pAtom->setName(newName);
				}
			}
		}

		else if (inputFileContents[i][0] == "createBond") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: createBond

			 Description: Create bond

			 syntax: createBond /1FEE/cuCYM4//CU.. /1FEE/cuCYM4//.SG.
			 \endcode
			 */
			if (inputFileContents[i].size() != 3) {
				MTKpp::errorLogger.throwError("MCPB::createBond",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			selection* pSeln1 = new selection(pCollection);
			failure = pSeln1->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::createBond",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			std::vector<atom*> selnAtoms1 = pSeln1->getAtoms();
			if (selnAtoms1.size() == 0) {
				MTKpp::errorLogger.throwError("MCPB::createBond",
							      " No atoms found in selection 1 ... exiting ", MTK_ERROR);
				exit(1);
			}

			selection* pSeln2 = new selection(pCollection);
			failure = pSeln2->parse(inputFileContents[i][2]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::createBond",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			std::vector<atom*> selnAtoms2 = pSeln2->getAtoms();
			if (selnAtoms2.size() == 0) {
				MTKpp::errorLogger.throwError("MCPB::createBond",
							      " No atoms found in selection 2 ... exiting ", MTK_ERROR);
				exit(1);
			}

			for (unsigned int a1 = 0; a1 < selnAtoms1.size(); a1++) {
				molecule* pLocalMol1 = selnAtoms1[a1]->getParent()->getParent();
				for (unsigned int a2 = 0; a2 < selnAtoms2.size(); a2++) {
					molecule* pLocalMol2 = selnAtoms2[a2]->getParent()->getParent();
					if (pLocalMol1 == pLocalMol2) {
						if (!pLocalMol1->hasBond(selnAtoms1[a1], selnAtoms2[a2])) {
							vector3d* coord1 = selnAtoms1[a1]->getCoords();
							vector3d* coord2 = selnAtoms2[a2]->getCoords();
							double bondDistance = coord1->dist(*coord2);
							Bond* pBond = pLocalMol2->addBond(selnAtoms1[a1], selnAtoms2[a2], 1, 2, 0, bondDistance);
							if (pBond) {
								selnAtoms1[a1]->addBondedAtom(selnAtoms2[a2]);
								selnAtoms2[a2]->addBondedAtom(selnAtoms1[a1]);
							}
						}
					}
					else {
						MTKpp::errorLogger.throwError("MCPB::createBond",
									      " Trying to Add a Bond between two different molecules ... exiting ", MTK_ERROR);
						exit(1);
					}
				}
			}
			connections* pConnections = new connections(pCollection);
			molecule* selMol = pSeln1->getMol();
			pConnections->assignAngles(selMol);
			pConnections->assignTorsions(selMol);
			delete pConnections;
			delete pSeln1;
			delete pSeln2;
		}

		else if (inputFileContents[i][0] == "print") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: print

			 Description: Print to screen details of structure

			 syntax: print cuCYM4
			 \endcode
			 */
			selection* pSel = new selection(pCollection);
			failure = pSel->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::print",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}

			std::cout << "---------------------------------" << std::endl;

			if (pSel->getSelectionType() == 0) {
				std::vector<molecule*> molList = pCollection->getMoleculeList();
				for (unsigned int k = 0; k < molList.size(); k++) {
					molecule* selMol = molList[k];
					std::cout << " m:" << k+1 << " " << selMol->getName() << std::endl;

					std::vector<submolecule*> residues = selMol->getSubMoleculeList();

					double charge = 0.0;
					for (unsigned int h = 0; h < residues.size(); h++) {
						std::cout << "  r: " << residues[h]->getName() << " " << residues[h]->getSubMolId() << std::endl;
						std::vector<atom*> residueAtoms = residues[h]->getAtomList();
						std::cout <<   "       atomName FileID Type Valence Hybridization   (Charge)" << std::endl;
						double resCharge = 0.0;
						for (unsigned int y = 0; y < residueAtoms.size(); y++) {
							stdAtom* pLStdAtom = residueAtoms[y]->getStdAtom();
							std::cout << "       "							
							 << residueAtoms[y]->getName() << " " << residueAtoms[y]->getFileID()
							 << " " << residueAtoms[y]->getType() << " " << residueAtoms[y]->getValence()
							 << " " << residueAtoms[y]->getHybridization() << "   ";
							if (pLStdAtom) {
								std::cout << "(" << pLStdAtom->atmCharge << ")";
								charge += pLStdAtom->atmCharge;
								resCharge += pLStdAtom->atmCharge;
							}
							std::cout << " " << std::endl;
						}
						std::cout << "   Residue Charge = " << charge << std::endl;
					}
					std::cout << "   Molecular Charge = " << charge << std::endl;

					std::map<int, Bond*> clusterBonds = selMol->getBondMap();
					if (clusterBonds.size() > 0) {
						std::cout <<   "       at1-at2 Type Topology " << std::endl;
						for (BondMapIterator b = clusterBonds.begin(); b != clusterBonds.end(); b++) {
							Bond* pBond = b->second;
							std::cout << "       " << pBond->atom1->getFileID() << "-" << pBond->atom2->getFileID()
							<< " " << pBond->type << " " << pBond->topology << std::endl;
						}
					}
				}
			}
			if (pSel->getSelectionType() == 1) {
				molecule* selMol = pSel->getMol();
				std::cout << " cluster: " << selMol->getName() << std::endl;

				std::vector<submolecule*> residues = selMol->getSubMoleculeList();

				for (unsigned int h = 0; h < residues.size(); h++) {
					std::cout << "  residue: " << residues[h]->getName() << std::endl;
					std::vector<atom*> residueAtoms = residues[h]->getAtomList();
					std::cout <<   "       atomName FileID Type Valence Hybridization " << std::endl;
					for (unsigned int y = 0; y < residueAtoms.size(); y++) {
						std::cout << "       " << residueAtoms[y]->getName() << " " << residueAtoms[y]->getFileID()
						<< " " << residueAtoms[y]->getType() << " " << residueAtoms[y]->getValence()
						<< " " << residueAtoms[y]->getHybridization()
						<< std::endl;
					}
				}

				std::map<int, Bond*> clusterBonds = selMol->getBondMap();
				if (clusterBonds.size() > 0) {
					std::cout <<   "       at1-at2 Type Topology " << std::endl;
					for (BondMapIterator b = clusterBonds.begin(); b != clusterBonds.end(); b++) {
						Bond* pBond = b->second;
						std::cout << "       " << pBond->atom1->getFileID() << "-" << pBond->atom2->getFileID()
						<< " " << pBond->type << " " << pBond->topology << std::endl;
					}
				}
			}
			std::cout << "---------------------------------" << std::endl;
		}

		else if (inputFileContents[i][0] == "capResidue") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: capResidue

			 Description: Cap residue R1 with NME and ACE

			 syntax: capResidue /1FEE/cuCYM4/CY1/.N.. ACE

			 syntax: capResidue /1FEE/cuCYM4/CY1/.C.. NME

			 syntax: capResidue /1FEE/cuCYM4/CY1/.SG. CH3
			 \endcode
			 */
			if (inputFileContents[i].size() != 3) {
				MTKpp::errorLogger.throwError("MCPB::capResidue",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}

			selection* pSeln = new selection(pCollection);
			failure = pSeln->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::capResidue",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			std::vector<atom*> selnAtoms = pSeln->getAtoms();

			stdLibrary* pStdLibrary = pCollection->getStdLibrary();
			stdFrag* pStdFrag = pStdLibrary->getStdFrag(inputFileContents[i][2]);
			if (!pStdFrag) {
				MTKpp::errorLogger.throwError("MCPB::capResidue",
							      " Cannot find standard fragment: " + inputFileContents[i][2], MTK_ERROR);
				exit(1);
			}

			// Get 3 other atoms bonded 1-2, 1-3, and 1-4 to selected atom
			atom* agAtom = 0;
			atom* trAtom = 0;

			std::vector<atom*> angleAtoms = selnAtoms[0]->getBondedAtoms();
			std::vector<atom*> trAtoms;

			if (angleAtoms.size() > 0) {
				agAtom = angleAtoms[0];
				trAtoms = agAtom->getBondedAtoms();
			}
			if (trAtoms.size() > 0) {
				for (unsigned int p = 0; p < trAtoms.size(); p++) {
					if ((trAtoms[p] != selnAtoms[0]) and (trAtoms[p] != agAtom) ) {
						trAtom = trAtoms[p];
						break;
					}
				}
			}
			if (!agAtom or !trAtom) {
				MTKpp::errorLogger.throwError("MCPB::capResidue",
							      " Error (make sure bonds are assigned) ...  exiting ", MTK_ERROR);
				exit(1);
			}

			std::vector<vector3d*> fragCoords;
			if (selnAtoms[0]->getName() == " N  ") { // if .N.. go backwards through the standard residue
				pStdFrag->generateCoordinates(selnAtoms[0]->getCoords(), agAtom->getCoords(), trAtom->getCoords(), 0);
				fragCoords = pStdFrag->getCoordinates();
			}
			else { // if .C.. or any other atom go forwards through the standard residue
				pStdFrag->generateCoordinates(selnAtoms[0]->getCoords(), agAtom->getCoords(), trAtom->getCoords(), 1);
				fragCoords = pStdFrag->getCoordinates();
			}

			molecule* clr = pSeln->getMol();
			if (!clr) {
				MTKpp::errorLogger.throwError("MCPB::capResidue",
							      " Error ...  exiting ", MTK_ERROR);
				exit(1);
			}
			submolecule* pSubMoleculeNew = clr->addSubMolecule();
			pSubMoleculeNew->setName(inputFileContents[i][2]);
			pSubMoleculeNew->setSubMolId(clr->getNumSubMolecules());

			// stdAtom's
			atom* pAtom = 0;
			std::vector<stdAtom*> atomList = pStdFrag->getStdAtomList();
			for (unsigned int k = 0; k < atomList.size(); k++) {
				pAtom = pSubMoleculeNew->addAtom();
				pAtom->setElement(pCollection->pElements->getElement(atomList[k]->atSymbol));
				pAtom->setFileID(k+1);
				pAtom->setCoords(fragCoords[k]->getX(),fragCoords[k]->getY(),fragCoords[k]->getZ());
				pAtom->setName(atomList[k]->identity);
			}
		}

		else if (inputFileContents[i][0] == "setFormalCharge") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: setFormalCharge

			 Description: Set Formal Charge on atom

			 syntax: setFormalCharge /1FEE/cuCYM4//CU.. 1
			 \endcode
			 */
			if (inputFileContents[i].size() != 3) exit(1);
			selection* pSeln = new selection(pCollection);
			failure = pSeln->parse(inputFileContents[i][1]);
			if (failure) {
				MTKpp::errorLogger.throwError("MCPB::setFormalCharge",
							      " Error in selection parsing ... exiting ", MTK_ERROR);
				exit(1);
			}
			std::vector<atom*> selnAtoms = pSeln->getAtoms();

			int chg = atoi(inputFileContents[i][2].c_str());
			for (unsigned int j = 0; j < selnAtoms.size(); j++) {
				selnAtoms[j]->setFormalCharge(chg);
			}
		}

		else if (inputFileContents[i][0] == "writePdb") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writePdb

			 Description: Write pdb file

			 syntax: writePdb /1ABC/clr cuCYS.pdb
			 \endcode
			 */
			if (inputFileContents[i].size() < 2) {
				MTKpp::errorLogger.throwError("MCPB::writePdb",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else if (inputFileContents[i].size() == 2) {
				pPdbParser->Write(inputFileContents[i][1], pCollection);
			}
			else if (inputFileContents[i].size() == 3) {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::writePdb",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}

				if (pSel->getSelectionType() == 1) {
					molecule* selMol = pSel->getMol();
					pPdbParser->Write(inputFileContents[i][2], selMol);
				}
				else {
				    MTKpp::errorLogger.throwError("MCPB::writePdb",
								      " Error molecule was not selected ... exiting ", MTK_ERROR);
			        exit(1);
				}
			}
		}

		else if (inputFileContents[i][0] == "writeMol") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writeMol

			 Description: Write mol file

			 syntax: writeMol /1ABC/clr cuCYS.mol
			 \endcode
			 */
			if (inputFileContents[i].size() < 3) {
				MTKpp::errorLogger.throwError("MCPB::writeMol",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else if (inputFileContents[i].size() == 3) {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::writeMol",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}

				if (pSel->getSelectionType() == 1) {
					molecule* selMol = pSel->getMol();
					if (selMol) {
						molParser* pMolParser = new molParser();
						pMolParser->Write(inputFileContents[i][2], selMol);
						delete pMolParser;
					}
				}
			}
		}

		else if (inputFileContents[i][0] == "writeSdf") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writeSdf

			 Description: Write sd file

			 syntax: writeSdf /1ABC/clr CLR.sdf
			 \endcode
			 */
			if (inputFileContents[i].size() < 3) {
				MTKpp::errorLogger.throwError("MCPB::writeSdf",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else if (inputFileContents[i].size() == 3) {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::writeSdf",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}

				if (pSel->getSelectionType() == 0) {
					sdfParser* pSdfParser = new sdfParser();
					pSdfParser->Write(inputFileContents[i][2], pCollection);
					delete pSdfParser;
				}

				if (pSel->getSelectionType() == 1) {
					molecule* selMol = pSel->getMol();
					if (selMol) {
						baseParser* pBaseParser = new baseParser();
						sdfParser* pSdfParser = new sdfParser();

						// remove all funcGroups
						typedef std::map<std::string, std::string>::iterator PropertyMapIterator;
						std::map<std::string, std::string> molProperties = selMol->getProperties();

						std::vector<std::string> toBeDeleted;
						for (PropertyMapIterator m = molProperties.begin(); m != molProperties.end(); m++) {
							std::string functionalGroup = m->first;
							std::string::size_type loc = functionalGroup.find("FuncGroup:", 0);
							if (loc != std::string::npos) toBeDeleted.push_back(m->first);
						}

						for (unsigned int f = 0; f < toBeDeleted.size(); f++) {
							selMol->delProperty(toBeDeleted[f]);
						}

						std::vector<std::string> fGroups;
						std::vector<submolecule*> submols = selMol->getSubMoleculeList();
						for (unsigned int s = 0; s < submols.size(); s++) {
							if (submols[s]->hasStdFrag()) {
								stdFrag* sFrag = submols[s]->getStdFrag();
								if (sFrag) {
									std::string fragName = sFrag->getSymbol();
									fGroups.push_back(fragName);
									int numFuncGroup = 0;
									for (unsigned int p = 0; p < fGroups.size(); p++) {
										if (fGroups[p] == fragName) numFuncGroup++;
									}

									std::stringstream funcGroupNumber;
									funcGroupNumber << numFuncGroup;
									std::string fName = "FuncGroup:" + fragName
									+ ":" + funcGroupNumber.str();
									std::string fValue = "";

									std::vector<atom*> fragAtomList = submols[s]->getAtomList();
									for (unsigned int x = 0; x < fragAtomList.size(); x++) {
										std::stringstream stdAtomIndex;
										std::stringstream atomIndex;

										stdAtom* pStdAtom = fragAtomList[x]->getStdAtom();
										if (pStdAtom) {
											stdAtomIndex << pStdAtom->index;
											atomIndex << fragAtomList[x]->getIndex();

											fValue += stdAtomIndex.str() + ";" +
											atomIndex.str()    + "|";
										}
									}
									selMol->addProperty(fName, fValue);
								}
							}
						}

						pSdfParser->Write(pBaseParser->OpenFile(inputFileContents[i][2]), selMol);
						delete pSdfParser;
					}
				}
			}
		}

		else if (inputFileContents[i][0] == "readSdf") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: readSdf

			 Description: Read sd file

			 syntax: readSdf 1ABC CLR.sdf
			 \endcode
			 */
			if (inputFileContents[i].size() < 3) {
				MTKpp::errorLogger.throwError("MCPB::readSdf",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else if (inputFileContents[i].size() == 3) {
				sdfParser* pSdfParser = new sdfParser();
				pSdfParser->Read(inputFileContents[i][2], pCollection);
				delete pSdfParser;
			}
			if (pCollection->getNumberMolecules() > 1) {
				MTKpp::errorLogger.throwError("MCPB::readSdf",
							      " Expecting just one molecule ... exiting ", MTK_ERROR);
				exit(1);
			}
			molecule* pMol = pCollection->getMolecule(1);
			if (!pMol) exit(1);

			pMol->setName(inputFileContents[i][1]);
			typedef std::map<std::string, std::string>::iterator PropertyMapIterator;
			std::map<std::string, std::string> molProperties = pMol->getProperties();
			for (PropertyMapIterator m = molProperties.begin(); m != molProperties.end(); m++) {
				std::string functionalGroup = m->first;
				std::string::size_type loc = functionalGroup.find("FuncGroup:", 0);
				if (loc != std::string::npos) {
					std::vector<std::string> vAtomPairs;
					splitString(m->second, "|", vAtomPairs, 0);

					std::vector<std::string> vfName;
					splitString(functionalGroup, ":", vfName, 0);
					if (vfName.size() == 3) {
						stdLibrary* pStdLib = pCollection->getStdLibrary();
						if (pStdLib) {
							stdFrag* pStdFg = pStdLib->getStdFrag(vfName[1]);
							if (pStdFg) {
								for (unsigned int a = 0; a < vAtomPairs.size(); a++) {
									std::vector<std::string> vAtomPair;
									splitString(vAtomPairs[a], ";", vAtomPair, 0);
									stdAtom* pStdAt = pStdFg->getStdAtom(atoi(vAtomPair[0].c_str()));
									atom* pAt = pMol->getAtom(atoi(vAtomPair[1].c_str()), 0, 1);
									if (pStdAt and pAt) {
										pAt->setName(pStdAt->identity);
										pAt->setStdAtom(pStdAt);
									}
								}
							}
						}
					}
				}
			}
			connections* pConnections = new connections(pCollection);
			pConnections->assignAngles(pMol);
			pConnections->assignStdBondsAngles(pMol);
			// uncomment when we extend Hoops et. al.
			//pConnections->assignTorsions(pMol);
			//pConnections->assignImpropers(pMol);
			//pConnections->assignStd(pMol);
			delete pConnections;
		}

		//////////////////////////////////////
		////////////// GAUSSIAN //////////////
		//////////////////////////////////////

		else if (inputFileContents[i][0] == "levelOfTheory") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: levelOfTheory

			 Description: Set Gaussian Theory Level

			 syntax: levelOfTheory HF
			 \endcode
			 */
			pGParser->setTheory(inputFileContents[i][1]);
		}

		else if (inputFileContents[i][0] == "basisSet") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: basisSet

			 Description: Set Gaussian Basis Set

			 syntax: basisSet 6-31G*

			 syntax: basisSet GEN bs.txt
			 syntax: basisSet GEN.6D.7F bs.txt
			 \endcode
			 */
			if (inputFileContents[i].size() == 2) {
				pGParser->setBasisSet(inputFileContents[i][1]);
			}
			else if (inputFileContents[i].size() == 3) {
				inputFileContents[i][1] = replaceCharacter(inputFileContents[i][1], '.', ' ');
				pGParser->setBasisSet(inputFileContents[i][1]);
				pGParser->setBasisSetFile(inputFileContents[i][2]);
			}
			else {
				MTKpp::errorLogger.throwError("MCPB::basisSet",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}

		else if (inputFileContents[i][0] == "pseudoPotentials") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: pseudoPotentials

			 Description: Specify pseudopotentials

			 syntax: pseudoPotential pseudo.txt
			 \endcode
			 */
			if (inputFileContents[i].size() == 2) {
				pGParser->addCommandOption("Pseudo","Read");
				pGParser->setPseudoPotentialFile(inputFileContents[i][1]);
			} else {
				MTKpp::errorLogger.throwError("MCPB::pseudoPotentials",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}

		}

		else if (inputFileContents[i][0] == "modRedundant") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: modRedundant

			 Description: Specify redundant coordinates

			 syntax: modRedundant modred.txt
			 \endcode
			 */
			if (inputFileContents[i].size() == 2) {
				pGParser->addCommandOption("Opt","ModRedundant");
				pGParser->setModRedundantFile(inputFileContents[i][1]);
			} else {
				MTKpp::errorLogger.throwError("MCPB::modRedundant",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}

		else if (inputFileContents[i][0] == "clusterCharge") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: clusterCharge

			 Description: Set Gaussian Charge

			 syntax: clusterCharge cuCYM4 -3
			 \endcode
			 */
			if (inputFileContents[i].size() == 3) {
				pGParser->setCharge(atoi(inputFileContents[i][2].c_str()));
			} else {
				MTKpp::errorLogger.throwError("MCPB::clusterCharge",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "clusterSpin") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: clusterSpin
			 
			 Description: Set Gaussian Spin
			 
			 syntax: clusterSpin 0
			 \endcode
			 */
			if (inputFileContents[i].size() == 2) {
				pGParser->setMultiplicity(atoi(inputFileContents[i][1].c_str()));
			} else {
				MTKpp::errorLogger.throwError("MCPB::clusterSpin",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "gaussianMoldenFormat") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: gaussianMoldenFormat
			 
			 Description: Use Molden formatted output in log file. Print out
			 details of the basis set and the molecular orbitals.
			 
			 syntax: gaussianMoldenFormat (bare word)
			 \endcode
			 */
			pGParser->addCommandOption("GFInput");
			pGParser->addIop("IOp(6/7=3)");
		}
		
		else if (inputFileContents[i][0] == "gaussianVerbosity") {
			if (inputFileContents[i].size() == 2) {
				pGParser->setVerbosity(inputFileContents[i][1]);
			} else {
				MTKpp::errorLogger.throwError("MCPB::gaussianVerbosity",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "gaussianMem") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: gaussianMem
			 
			 Description: Set amount of memory requested for Gaussian
			 
			 syntax: gaussianMem 3600MB
			 \endcode
			 */
			if (inputFileContents[i].size() == 2) {
				pGParser->setMem(inputFileContents[i][1]);
			} else {
				MTKpp::errorLogger.throwError("MCPB::gaussianMem",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "gaussianNProc") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: gaussianNProc
			 
			 Description: Set number of processors requested for Gaussian
			 
			 syntax: gaussianNProc 2
			 \endcode
			 */
			if (inputFileContents[i].size() == 2) {
				pGParser->setNProc(inputFileContents[i][1]);
			} else {
				MTKpp::errorLogger.throwError("MCPB::gaussianNProc",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "gaussianOptAndFC") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: gaussianOptAndFC
			 
			 Description: Set Gaussian input name
			 
			 syntax: gaussianOptAndFC //cuCYM4 cuCYM4.com
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::gaussianOptAndFC",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::gaussianOptAndFC",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				if (pSel->getSelectionType() == 1) {
					molecule* pSelMol = pSel->getMol();
					std::string fileNameBase = baseName(inputFileContents[i][2]);
					pGParser->setChkPt(fileNameBase+"_opt.chk");
					
					// optimization
					pGParser->setCartesian(1);
					pGParser->setWriteMoleculeName(1);
					pGParser->setWriteChargeAndMult(1);
					pGParser->addCommandOption("Integral", "(Grid=UltraFine)");
					
					std::vector<std::string> optOptions = pGParser->getCommandOption("Opt");
					if (optOptions.size() == 0) {
						pGParser->addCommandOption("Opt");
					}
					
					pGParser->addCommandOption("SCF", "XQC");
					pGParser->addCommandOption("Geom", "PrintInputOrient");
					pGParser->Write(fileNameBase+"_opt.com", pSelMol);
					pGParser->removeCommandOption("Opt");
					
					// k's
					// Do not write coordinates, name, charge or multiplicity,
					// since we are using allcheckpoint
					pGParser->setNoCoords();
					pGParser->setWriteMoleculeName(0);
					pGParser->setWriteChargeAndMult(0);
					pGParser->addCommandOption("Freq", "NoRaman");
					pGParser->addCommandOption("Geom", "AllCheckpoint");
					pGParser->addCommandOption("Guess", "Read");
					pGParser->addIop("IOp(7/33=1)");
					pGParser->Write(fileNameBase+"_fc.com", pSelMol);
					pGParser->clearIop();
					
					/*
					 // optimization
					 pGParser->addCommandOption("Opt", "Z-Matrix");
					 pGParser->generateZMatrix(pSelMol);
					 // [molecule index] = zmatrix index;
					 pGParser->writeMappingFile(fileNameBase+"_zmat.map");
					 pGParser->addCommandOption("Geom", "PrintInputOrient");
					 pGParser->Write(fileNameBase+"_opt.com", pSelMol);
					 pGParser->removeCommandOption("Opt");
					 */
				}
			}
		}
		
		else if (inputFileContents[i][0] == "gaussianCharges") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: gaussianCharges
			 
			 Description: Set Gaussian input name
			 
			 syntax: gaussianCharges //cuCYM4 cuCYM4.com
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::gaussianCharges",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::gaussianCharges",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				if (pSel->getSelectionType() == 1) {
					molecule* pSelMol = pSel->getMol();
					std::string fileNameBase = baseName(inputFileContents[i][2]);
					pGParser->setChkPt(fileNameBase+"_mk.chk");
					pGParser->setWriteMoleculeName(1);
					pGParser->setWriteChargeAndMult(1);
					pGParser->setCartesian(1);
					
					// q's
					std::vector<std::string> popOptions;
					popOptions.push_back("MK");
					popOptions.push_back("ReadRadii");
					pGParser->addCommandOption("Pop", popOptions);
					pGParser->removeCommandOption("Freq");
					pGParser->addIop("IOp(6/33=2)");
					
					pGParser->addCommandOption("Integral", "(Grid=UltraFine)");
					pGParser->addCommandOption("SCF", "XQC");
					
					pGParser->Write(fileNameBase+"_mk.com", pSelMol);
					
					// Create espgen script
					// espgen -i fileNameBase+"_mk.log" -o fileNameBase+"_mk.esp"
					
					// Create resp files
					// 1st line: title
					// Create resp script
					
					// Create script to update charges in lib file
				}
			}
		}
		
		else if (inputFileContents[i][0] == "respgenAdditions") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: respgenAdditions
			 
			 Description: Add info to respgen files
			 
			 syntax: respgenAdditions groupName fileName bb
			 
			 bb Definitions:
			 - 0 No restraints
			 - 1 Heavy Atoms in Backbone (bb_heavy, [ca, n, c, o]) set to parm94 values
			 - 2 Atoms in Backbone (bb, [ca, h, ha, n, nh, c, o]) set to parm94 values
			 - 3 Atoms in Backbone plus CB (bbb, [ca, h, ha, n, nh, c, o, cb]) set to parm94 values
			 \endcode
			 */
			if ((inputFileContents[i].size() < 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::respgenAdditions",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				int iBB = 0;
				if (inputFileContents[i].size() == 4) {
					iBB = string2Int(inputFileContents[i][3]);
				}
				
				std::vector<std::string>::iterator result;
				acParser* pAcParser = new acParser();
				pAcParser->Write(inputFileContents[i][2]+".ac", pCollection);
				delete pAcParser;
				std::string respAddFile = inputFileContents[i][2] + "_respAdds";
				std::ofstream orespAdd;
				orespAdd.open(respAddFile.c_str());
				
				if (!orespAdd) {
					MTKpp::errorLogger.throwError("MCPB::respgenAdditions",
								      " Unable to open RESP additions file ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				std::vector<molecule*> molList = pCollection->getMoleculeList();
				
				char temp[100];
				for (unsigned int m = 0; m < molList.size(); m++) {
					molecule* pMol = molList[m];
					std::vector<submolecule*> submols = pMol->getSubMoleculeList();
					for (unsigned int s = 0; s < submols.size(); s++) {
						if (submols[s]->hasStdFrag()) {
							stdFrag* pStdFrag = submols[s]->getStdFrag();
							if (pStdFrag->getParent()->getName() != inputFileContents[i][1]) {
								sprintf(temp,"%5d%10.5f",
									pStdFrag->numStdAtoms(), 0.0);
								orespAdd << temp << std::endl;
								std::vector<atom*> atomList = submols[s]->getAtomList();
								char temp2[80];
								int counter = 0;
								for (unsigned int a = 0; a < atomList.size(); a++) {
									if (counter > 7) { // was 6
										orespAdd << "" << std::endl;
										counter = 0;
									}
									sprintf(temp2,"%5d%5d", 1, atomList[a]->getIndex());
									orespAdd << temp2;
									counter++;
								}
								orespAdd << "" << std::endl;
							}
						}
					}
				}
				if (iBB) {
					std::string freezingAtoms = " \n Atoms frozen: ";
					for (unsigned int m = 0; m < molList.size(); m++) {
						molecule* pMol = molList[m];
						std::vector<submolecule*> submols = pMol->getSubMoleculeList();
						for (unsigned int s = 0; s < submols.size(); s++) {
							if (submols[s]->hasStdFrag()) {
								stdFrag* pStdFrag = submols[s]->getStdFrag();
								if (pStdFrag->getParent()->getName() != inputFileContents[i][1]) continue;
								std::vector<atom*> subAtoms = submols[s]->getAtomList();
								for (unsigned int a = 0; a < subAtoms.size(); a++) {
									stdAtom* pStdAtom_a = subAtoms[a]->getStdAtom();
									
									if (iBB == 1 and pStdAtom_a) {
										std::string dd = pStdAtom_a->identity + ":" + pStdAtom_a->type;
										
										//result = std::find(bb_heavy.begin(), bb_heavy.end(), subAtoms[a]->getName());
										//result = std::find(bb_heavy.begin(), bb_heavy.end(), pStdAtom_a->identity);
										result = std::find(bb_heavy.begin(), bb_heavy.end(), dd);
										
										if (result != bb_heavy.end()) {
											sprintf(temp,"%5d%10.5f", 1, pStdAtom_a->atmCharge);
											freezingAtoms += "\n BB Heavy "  + submols[s]->getName() + " " + int2String(submols[s]->getSubMolId()) + " " + dd;
											
											//std::cout << dd << " " << subAtoms[a]->getName() << " ---> " << subAtoms[a]->getStdAtom()->type << std::endl;
											orespAdd << temp << std::endl;
											sprintf(temp,"%5d%5d", 1, subAtoms[a]->getIndex());
											orespAdd << temp << std::endl;
										}
									}
									else if (iBB == 2 and pStdAtom_a) {
										std::string dd = pStdAtom_a->identity + ":" + pStdAtom_a->type;
										
										//result = std::find(bb.begin(), bb.end(), subAtoms[a]->getName());
										result = std::find(bb.begin(), bb.end(), dd);
										
										if (result != bb.end()) {
											freezingAtoms += "\n BB "  + submols[s]->getName() + " " + int2String(submols[s]->getSubMolId()) + " " + dd;
											sprintf(temp,"%5d%10.5f", 1, pStdAtom_a->atmCharge);
											orespAdd << temp << std::endl;
											sprintf(temp,"%5d%5d", 1, subAtoms[a]->getIndex());
											orespAdd << temp << std::endl;
										}
									}
									else if (iBB == 3 and pStdAtom_a) {
										std::string dd = pStdAtom_a->identity + ":" + pStdAtom_a->type;
										
										//result = std::find(bbb.begin(), bbb.end(), subAtoms[a]->getName());
										result = std::find(bbb.begin(), bbb.end(), dd);
										
										if (result != bbb.end()) {
											freezingAtoms += "\n BBB "  + submols[s]->getName() + " "
											              + int2String(submols[s]->getSubMolId()) + " " + dd;
											sprintf(temp,"%5d%10.5f", 1, pStdAtom_a->atmCharge);
											orespAdd << temp << std::endl;
											sprintf(temp,"%5d%5d", 1, subAtoms[a]->getIndex());
											orespAdd << temp << std::endl;
										}
									}
									else {
										MTKpp::errorLogger.throwError("MCPB::respgenAdditions",
													      " Unknown option ... exiting ", MTK_ERROR);
									}
								}
							}
						}
					}
					MTKpp::errorLogger.throwError("MCPB::respgenAdditions",
								      freezingAtoms, INFO);
				}
				orespAdd << "\n" << std::endl;
				orespAdd.close();
			}
		}
		
		else if (inputFileContents[i][0] == "setMKRadii") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: setMKRadii
			 
			 Description: Set Merz-Kollman radii for element
			 
			 syntax: setMKRadii cu 0.91
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pGParser)) {
				MTKpp::errorLogger.throwError("MCPB::setMKRadii",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				double radii = strtod(inputFileContents[i][2].c_str(), 0);
				pGParser->setMKRadii(inputFileContents[i][1], radii);
			}
		}

		else if (inputFileContents[i][0] == "readGaussianOutput") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: readGaussianOutput
			 
			 Description: Read Gaussian Output
			 
			 syntax: readGaussianOutput cuCYM4.log
			 \endcode
			 */
			if ((inputFileContents[i].size() == 2) or (pGParser)) {
				pGParser->Read(inputFileContents[i][1], pSheet);
			}
			else {
				MTKpp::errorLogger.throwError("MCPB::readGaussianOutput",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "readMolZmatMapping") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: readMolZmatMapping
			 
			 Description: Read Molecule <--> Z-Matrix mapping file
			 
			 syntax: readMolZmatMapping file.map
			 \endcode
			 */
			if ((inputFileContents[i].size() == 2) and (pGParser)) {
				pGParser->readMappingFile(inputFileContents[i][1]);
			}
			else {
				MTKpp::errorLogger.throwError("MCPB::readMolZmatMapping",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "readFormattedChkPtFile") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: readFormattedChkPtFile
			 
			 Description:
			 
			 syntax: readFormattedChkPtFile file.fchk
			 \endcode
			 */
			if ((inputFileContents[i].size() == 2) or (pGParser)) {
				pGParser->readFormattedChkPtFile(inputFileContents[i][1], pSheet);
			}
			else {
				MTKpp::errorLogger.throwError("MCPB::readFormattedChkPtFile",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "writeLib") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writeLib
			 
			 Description: Write standard library
			 
			 syntax: writeLib groupName fileName.xml
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::writeLib",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				parameters* pParms = pCollection->getParameters();
				if (!pParms) {
					MTKpp::errorLogger.throwError("MCPB::writeLib",
								      " Please read in parameters first ... exiting ", MTK_ERROR);
					exit(1);
				}
				stdLibrary* pStdLib = pCollection->getStdLibrary();
				stdLibParser* pStdLibParser = new stdLibParser(pCollection, pStdLib, pParms);
				pStdLibParser->Write(inputFileContents[i][2], inputFileContents[i][1]);
				delete pStdLibParser;
			}
		}
		
		else if (inputFileContents[i][0] == "writeParams") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writeParams
			 
			 Description: Write all new parameters
			 
			 syntax: writeParams groupName fileName.xml
			 \endcode
			 */
			if ((!pCollection) or (inputFileContents[i].size() != 3)) {
				MTKpp::errorLogger.throwError("MCPB::writeParams",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			paramParser* pParamParser = new paramParser(pCollection->getParameters());
			if (pParamParser) {
				pParamParser->Write(inputFileContents[i][2], inputFileContents[i][1]);
				delete pParamParser;
			}
		}
		
		else if (inputFileContents[i][0] == "listFragments") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: listFragments
			 
			 Description: List available fragments in a particular library
			 
			 syntax: listFragments terminal
			 \endcode
			 */
			if ((inputFileContents[i].size() != 2) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::listFragments",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				stdLibrary* pStdLibrary = pCollection->getStdLibrary();
				if (pStdLibrary) {
					stdGroup* pStdGroup = pStdLibrary->getStdGroup(inputFileContents[i][1]);
					if (pStdGroup) {
						std::vector<stdFrag*> fragList = pStdGroup->getStdFragList();
						for (unsigned int f = 0; f < fragList.size(); f++) {
							std::cout << " 3L: " << fragList[f]->getSymbol() << "| Name: "
							<< fragList[f]->getName() << std::endl;
						}
					}
				}
			}
		}
		
		else if (inputFileContents[i][0] == "addFragment") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: addFragment
			 
			 Description: Add Fragment to atom
			 
			 syntax: addFragment 6MemRings/6CH bd /col/Mol//34 ag /col/Mol//27 tr /col/Mol//9 180.0
			 \endcode
			 */
			if (inputFileContents[i].size() < 9) {
				MTKpp::errorLogger.throwError("MCPB::addFragment",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else if (inputFileContents[i].size() == 9) {
				
				std::vector<std::string> words;
				splitString(inputFileContents[i][1], "/", words, 0);
				if (words.size() != 2) {
					MTKpp::errorLogger.throwError("MCPB::addFragment",
								      " Incorrect use ... exiting ", MTK_ERROR);
					exit(1);
				}
				words[1] = replaceCharacter(words[1], '.', ' ');
				
				stdLibrary* pStdLib = pCollection->getStdLibrary();
				
				stdGroup* pStdGroup = pStdLib->getStdGroup(words[0]);
				
				if (!pStdGroup) {
					MTKpp::errorLogger.throwError("MCPB::addFragment",
								      " Incorrect use ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				stdFrag* pStdFrag = pStdGroup->getStdFrag(words[1]);
				if (!pStdFrag) {
					MTKpp::errorLogger.throwError("MCPB::addFragment",
								      " Incorrect use ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				selection* pSel1 = new selection(pCollection);
				failure = pSel1->parse(inputFileContents[i][3]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::addFragment",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				selection* pSel2 = new selection(pCollection);
				failure = pSel2->parse(inputFileContents[i][5]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::addFragment",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				selection* pSel3 = new selection(pCollection);
				failure = pSel3->parse(inputFileContents[i][7]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::addFragment",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				atom* bdAtom = 0;
				atom* agAtom = 0;
				atom* trAtom = 0;
				
				if (pSel1->getSelectionType() == 3) {
					bdAtom = pSel1->getAtom();
					pSubMolecule = bdAtom->getParent();
					pMolecule = pSubMolecule->getParent();
				}
				
				if (pSel2->getSelectionType() == 3) {
					agAtom = pSel2->getAtom();
				}
				
				if (pSel3->getSelectionType() == 3) {
					trAtom = pSel3->getAtom();
				}
				
				// Check that there is an available site
				// right now assume there is
				
				if (bdAtom and agAtom and trAtom and pSubMolecule) {
					stdAtom* pStdAtom = pStdFrag->getStdAtom(1);
					double oldTorsion = pStdAtom->bondTorsion;
					pStdAtom->bondTorsion = strtod(inputFileContents[i][8].c_str(), 0);
					pStdFrag->generateCoordinates(bdAtom->getCoords(), agAtom->getCoords(), trAtom->getCoords(), 1);
					pStdAtom->bondTorsion = oldTorsion;
					std::vector<vector3d*> stdFragCoords = pStdFrag->getCoordinates();
					
					int startPt = pMolecule->getAtomIndex();
					// stdAtoms
					std::vector<stdAtom*> stdFragAtoms = pStdFrag->getStdAtomList();
					pSubMolecule = pMolecule->addSubMolecule();
					pSubMolecule->setName(pStdFrag->getSymbol());
					pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());
					
					for (unsigned int k = 0; k < stdFragAtoms.size(); k++) {
						atom* pAtom = pSubMolecule->addAtom();
						pAtom->setElement(pCollection->pElements->getElement(stdFragAtoms[k]->atSymbol));
						pAtom->setFileID(pMolecule->getNumAtoms()+1);
						pAtom->setCoords(stdFragCoords[k]->getX(),stdFragCoords[k]->getY(),stdFragCoords[k]->getZ());
						pAtom->setName(stdFragAtoms[k]->identity);
					}
					
					for (unsigned int k = 0; k < stdFragAtoms.size(); k++) {
						if (stdFragAtoms[k]->chain == "M") {
							int at1 = startPt + k;
							atom* pBondAtom1 = pMolecule->getAtom(at1, 1, 0);
							
							Bond* pBond = pMolecule->addBond(pBondAtom1, bdAtom, 1, 0, 0, 0.0);
							if (!pBond) {
								MTKpp::errorLogger.throwError("MCPB::addFragment",
											      " Error ... exiting ", MTK_ERROR);
								exit(1);
							}
							pBondAtom1->addBondedAtom(bdAtom);
							bdAtom->addBondedAtom(pBondAtom1);
							break;
						}
					}
					
					// stdBonds
					std::vector<stdBond*> stdBondList = pStdFrag->getStdBondList();
					for (unsigned int k = 0; k < stdBondList.size(); k++) {
						if ((stdBondList[k]->atom1 < 0) or (stdBondList[k]->atom2 < 0)) continue;
						int at1 = startPt + stdBondList[k]->atom1 - 1;
						int at2 = startPt + stdBondList[k]->atom2 - 1;
						int bondType = stdBondList[k]->type;
						int bondStereo = stdBondList[k]->stereo;
						int bondTopology = stdBondList[k]->topology;
						
						atom* pBondAtom1 = pMolecule->getAtom(at1, 1, 0);
						if (!pBondAtom1) {
							MTKpp::errorLogger.throwError("MCPB::addFragment",
										      " Error ... exiting ", MTK_ERROR);
							exit(1);
						}
						
						atom* pBondAtom2 = pMolecule->getAtom(at2, 1, 0);
						if (!pBondAtom2) {
							MTKpp::errorLogger.throwError("MCPB::addFragment",
										      " Error ... exiting ", MTK_ERROR);
							exit(1);
						}
						
						Bond* pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);
						if (!pBond) {
							MTKpp::errorLogger.throwError("MCPB::addFragment",
										      " Error ... exiting ", MTK_ERROR);
							exit(1);
						}
						pBondAtom1->addBondedAtom(pBondAtom2);
						pBondAtom2->addBondedAtom(pBondAtom1);
					}
					
					// stdLoop's
					std::vector<stdLoop*> loopList = pStdFrag->getStdLoopList();
					for (unsigned int k = 0; k < loopList.size(); k++) {
						if ((loopList[k]->atom1 < 0) or (loopList[k]->atom2 < 0)) continue;
						int at1 = startPt + loopList[k]->atom1 - 1;
						int at2 = startPt + loopList[k]->atom2 - 1;
						int bondType = loopList[k]->type;
						int bondStereo = loopList[k]->stereo;
						int bondTopology = 1;
						
						atom* pBondAtom1 = pMolecule->getAtom(at1, 1, 0);
						if (!pBondAtom1) {
							MTKpp::errorLogger.throwError("MCPB::addFragment",
										      " Error ... exiting ", MTK_ERROR);
							exit(1);
						}
						
						atom* pBondAtom2 = pMolecule->getAtom(at2, 1, 0);
						if (!pBondAtom2) {
							MTKpp::errorLogger.throwError("MCPB::addFragment",
										      " Error ... exiting ", MTK_ERROR);
							exit(1);
						}
						
						Bond* pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);
						if (!pBond) {
							MTKpp::errorLogger.throwError("MCPB::addFragment",
										      " Error ... exiting ", MTK_ERROR);
							exit(1);
						}
						pBondAtom1->addBondedAtom(pBondAtom2);
						pBondAtom2->addBondedAtom(pBondAtom1);
					}
				}
			}
		}
		
		else if (inputFileContents[i][0] == "tmpWriteZMatrix") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: tmpWriteZMatrix
			 
			 Description:
			 
			 syntax: tmpWriteZMatrix /1A5T/1A5T file.zmat
			 \endcode
			 */
			if (inputFileContents[i].size() == 3) {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					std::cout << " Error in selection parsing ... exiting " << std::endl;
					exit(1);
				}
				
				if (pSel->getSelectionType() == 1) {
					molecule* pSelMol = pSel->getMol();
					std::vector<atom*> myAtoms = pSelMol->getAtomList();
					table<double>* molCoords = pSheet->getTable("Current cartesian coordinates");
					//ublas::matrix<double> &molCoordsMatrix = molCoords->getMatrix();
					Eigen::Matrix<double, Dynamic, Dynamic> &molCoordsMatrix = molCoords->getMatrix();
					for (unsigned int z = 0; z < myAtoms.size(); z++) {
						myAtoms[z]->setCoords(molCoordsMatrix(z,0),
								      molCoordsMatrix(z,1),
								      molCoordsMatrix(z,2));
					}
					pGParser->generateZMatrix(pSelMol);
					// k's
					pGParser->setCartesian(0);
					pGParser->setWriteMoleculeName(1);
					pGParser->setWriteChargeAndMult(1);
					pGParser->addCommandOption("Freq", "NoRaman");
					pGParser->addIop("IOp(7/33=1)");
					pGParser->Write(inputFileContents[i][2], pSelMol);
					pGParser->writeMappingFile(inputFileContents[i][2]+".map");
				}
			}
		}
		
		else if (inputFileContents[i][0] == "updateForceConstants") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: updateForceConstants
			 
			 Description:
			 
			 syntax: updateForceConstants /COL/MOL Group X Y
			 X Values:
			 - 0 Do not update bonds and angle equilibrium values
			 - 1 Do update bonds and angle equilibrium values (req)
			 Y Values:
			 - 0 Seminario Method
			 - 1 Z-matrix Method
			 \endcode
			 */
			if ((inputFileContents[i].size() != 5) or (!pCollection)) {
				errorMessage = " Incorrect use of updateForceConstants";
				errorMessage += "   syntax: updateForceConstants /COL/MOL Group 1 1";
				MTKpp::errorLogger.throwError("MCPB::updateForceConstants", errorMessage, MTK_ERROR);
				exit(1);
			}
			else {
				int updateRvalues = string2Int(inputFileContents[i][3]);
				int fcMethod = string2Int(inputFileContents[i][4]);
				
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::updateForceConstants",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				Bond* pBond = 0;
				Angle* pAngle = 0;
				std::string eM = "\n";
				
				if (pSel->getSelectionType() == 1) {
					molecule* pSelMol = pSel->getMol();
					std::map<int, Bond*> molBonds = pSelMol->getBondMap();
					
					if (!molBonds.empty()) {
						for (BondMapIterator b = molBonds.begin(); b != molBonds.end(); b++) {
							pBond = b->second;
							if (pBond->pBondParam) {
								if (pBond->pBondParam->optimize) {
									std::vector<Bond*> allBondsSameType;
									for (BondMapIterator b2 = molBonds.begin(); b2 != molBonds.end(); b2++) {
										Bond* pBond2 = b2->second;
										if (pBond->pBondParam == pBond2->pBondParam) {
											allBondsSameType.push_back(pBond2);
										}
									}
									double averageReq = 0.0;
									double averageKeq = 0.0;
									for (unsigned int d = 0; d < allBondsSameType.size(); d++) {
										atom* pAt1 = allBondsSameType[d]->atom1;
										atom* pAt2 = allBondsSameType[d]->atom2;
										if (pAt1 and pAt2) {
											double bondForceConstant = 0.0;
											double bondSize = 0.0;
											
											int f = -1;
											if (fcMethod == 0 and pSheet) { // Seminario
												f = pGParser->getForceConstant(pSheet, pAt1->getFileID()-1, pAt2->getFileID()-1,
															       bondSize, bondForceConstant);
											}
											else if (fcMethod == 1){ // Z-Matrix
												f = pGParser->getForceConstantZMAT(pAt1->getFileID(),
																   pAt2->getFileID(),
																   bondSize,
																   bondForceConstant);
											}
											else {
												MTKpp::errorLogger.throwError("MCPB::updateForceConstants",
															      " Incorrect FC Method Defined ", MTK_ERROR);
												exit(1);
											}
											
											if (!f) {
												averageKeq += bondForceConstant / 2;
												averageReq += bondSize;
												bool lessThanZero = false;
												if (bondForceConstant < 0.0) {
													allBondsSameType[d]->pBondParam->keq = 0.0;
													bondForceConstant = 0.0;
													lessThanZero = true;
												}
												eM += "Bond: |" + pAt1->getStdAtom()->type + "|-|" + pAt2->getStdAtom()->type + "| "
												+ double2String(bondSize) + " " + double2String(bondForceConstant) + " " + int2String(lessThanZero)
												+ " " + int2String(allBondsSameType.size()) + "\n";
											}
											else {
												averageReq += bondSize;
												eM += "Bond: |" + pAt1->getStdAtom()->type + "|-|" + pAt2->getStdAtom()->type + "| Undefined \n";
											}
										}
									}
									averageReq = averageReq / double(allBondsSameType.size());
									averageKeq = averageKeq / double(allBondsSameType.size());
									for (unsigned int d = 0; d < allBondsSameType.size(); d++) {
										allBondsSameType[d]->pBondParam->keq = averageKeq;
										allBondsSameType[d]->pBondParam->optimize = false;
										if (updateRvalues) {
											allBondsSameType[d]->pBondParam->req = averageReq;
										}
									}
									eM += "  Average Bond: " + double2String(averageReq) + " " + double2String(averageKeq * 2) + "\n";
								}
							}
						}
					}
					
					std::map<ULONG_KIND, Angle*> molAngles = pSelMol->getAngleMap();
					if (!molAngles.empty()) {
						for (AngleMapIterator b = molAngles.begin(); b != molAngles.end(); b++) {
							pAngle = b->second;
							if (pAngle->pAngleParam) {
								if (pAngle->pAngleParam->optimize) {
									std::vector<Angle*> allAnglesSameType;
									for (AngleMapIterator b2 = molAngles.begin(); b2 != molAngles.end(); b2++) {
										Angle* pAngle2 = b2->second;
										if (pAngle->pAngleParam == pAngle2->pAngleParam) {
											allAnglesSameType.push_back(pAngle2);
										}
									}
									double averageReq = 0.0;
									double averageKeq = 0.0;
									for (unsigned int d = 0; d < allAnglesSameType.size(); d++) {
										atom* pAt1 = allAnglesSameType[d]->atom1;
										atom* pAt2 = allAnglesSameType[d]->atom2;
										atom* pAt3 = allAnglesSameType[d]->atom3;
										
										if (pAt1 and pAt2 and pAt3) {
											double angleForceConstant = 0.0;
											double angleSize = 0.0;
											
											int f = -1;
											if (fcMethod == 0 and pSheet) { // Seminario
												f = pGParser->getForceConstant(pSheet, pAt1->getFileID()-1, pAt2->getFileID()-1, pAt3->getFileID()-1,
															       angleSize, angleForceConstant);
											}
											else if (fcMethod == 1){ // Z-Matrix
												f = pGParser->getForceConstantZMAT(pAt1->getFileID(),
																   pAt2->getFileID(), pAt3->getFileID(), angleSize, angleForceConstant);
											}
											else {
												MTKpp::errorLogger.throwError("MCPB::updateForceConstants",
															      " Incorrect FC Method Defined ", MTK_ERROR);
												exit(1);
											}
											
											if (!f) {
												double aT = angleForceConstant;
												angleForceConstant = angleForceConstant / 2;
												bool lessThanZero = false;
												if (angleForceConstant < 0.0) {
													angleForceConstant = 0.0;
													lessThanZero = true;
												}
												
												eM += "Angle: |" + pAt1->getStdAtom()->type + "|-|" + pAt2->getStdAtom()->type + "|-|"
												+ pAt3->getStdAtom()->type + "| "
												+ double2String(angleSize) + " " + double2String(aT) + " " + int2String(lessThanZero)
												+ " " + int2String(allAnglesSameType.size()) + "\n";
												
												averageReq += angleSize;
												averageKeq += angleForceConstant;
											}
											else {
												averageReq += angleSize;
												eM += "Angle: |" + pAt1->getStdAtom()->type + "|-|" + pAt2->getStdAtom()->type + "|-|"
												+ pAt3->getStdAtom()->type + "| Undefined\n";
											}
										}
									}
									
									averageReq = averageReq / double(allAnglesSameType.size());
									averageKeq = averageKeq / double(allAnglesSameType.size());
									for (unsigned int d = 0; d < allAnglesSameType.size(); d++) {
										allAnglesSameType[d]->pAngleParam->optimize = false;
										allAnglesSameType[d]->pAngleParam->keq = averageKeq;
										if (updateRvalues) {
											allAnglesSameType[d]->pAngleParam->req = averageReq * DEG2RAD;
										}
									}
									eM += "  Average Angle: " + double2String(averageReq) + " " + double2String(averageKeq * 2) + "\n";
								}
							}
						}
					}
					
					MTKpp::errorLogger.throwError("MCPB::updateForceConstants", eM, INFO);
				}
			}
		}
		
		else if (inputFileContents[i][0] == "readRespCharges") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: readRespCharges
			 
			 Description:
			 
			 syntax: readRespCharges /COL/MOL file.resp2.chg
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::readRespCharges",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::readRespCharges",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				std::vector<double> molCharges;
				std::ifstream iRespChargeFile;
				iRespChargeFile.open(inputFileContents[i][2].c_str());
				
				if (!iRespChargeFile) {
					MTKpp::errorLogger.throwError("MCPB::readRespCharges",
								      " Unable to open RESP2 chg file ... exiting ", MTK_ERROR);
					exit(1);
				}
				std::string fileline = "";
				// - Read the file - //
				while (iRespChargeFile) {
					getline(iRespChargeFile,fileline);
					for (unsigned int d = 0; d < fileline.size(); d+=10) {
						std::string sChg = fileline.substr(d,10);
						molCharges.push_back(string2Double(sChg));
					}
					/*
					 std::vector<std::string> words;
					 splitString(fileline, " ", words, 0);
					 if (words.size() > 0) {
					 for (unsigned int h = 0; h < words.size(); h++) {
					 molCharges.push_back(string2Double(words[h]));
					 std::cout << " " << string2Double(words[h]);
					 }
					 words.clear();
					 }
					 //std::cout << " " <<std::endl;*/
				}
				
				if (pSel->getSelectionType() == 1) {
					molecule* pSelMol = pSel->getMol();
					if (!pSelMol) exit(1);
					std::vector<atom*> lAtomList = pSelMol->getAtomList();
					if (molCharges.size() != lAtomList.size()) {
						//std::cout << " N Charges = " << molCharges.size() << " N Atoms = " << lAtomList.size();
						MTKpp::errorLogger.throwError("MCPB::readRespCharges",
									      " Error in readRespCharges, charge file does not have the same\n number of atoms as the pdb file ... exiting ", MTK_ERROR);
					}
					for (unsigned int k = 0; k < lAtomList.size(); k++) {
						lAtomList[k]->setZcharge(molCharges[k]);
					}
				}
			}
		}
		
		else if (inputFileContents[i][0] == "updateRespCharges") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: updateRespCharges
			 
			 Description:
			 
			 syntax: updateRespCharges /COL/MOL Group
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::updateRespCharges",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::updateRespCharges",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				stdLibrary* pStdLibrary = pCollection->getStdLibrary();
				stdGroup* pStdGroup = 0;
				if (pStdLibrary) {
					pStdGroup = pStdLibrary->getStdGroup(inputFileContents[i][2]);
					if (!pStdGroup) {
						MTKpp::errorLogger.throwError("MCPB::updateRespCharges",
									      " Can't find stdGroup ... exiting ", MTK_ERROR);
						exit(1);
					}
				}
				else {
					MTKpp::errorLogger.throwError("MCPB::updateRespCharges",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
				}
				
				/*
				 molecule* pStdMolecule = 0;
				 bool bStdMolecule = false;
				 if (!pStdGroup->hasStdMolecule()) {
				 pStdMolecule = pCollection->addMolecule();
				 pStdMolecule->setName("Reference");
				 pStdGroup->setStdMolecule(pStdMolecule);
				 bStdMolecule = true;
				 }
				 else {
				 pStdMolecule = pStdGroup->getStdMolecule();
				 }
				 */
				std::string updateRespChargesMessage = "\n Group Name: " + pStdGroup->getName();
				
				double groupCharge = 0.0;
				if (pSel->getSelectionType() == 1) {
					molecule* pSelMol = pSel->getMol();
					if (!pSelMol) {
						exit(1);
					}
					std::vector<submolecule*> submolList = pSelMol->getSubMoleculeList();
					for (unsigned int t = 0; t < submolList.size(); t++) {
						stdFrag* pStdFrag = submolList[t]->getStdFrag();
						if (pStdFrag) {
							bool doIt = pStdGroup->hasStdFrag(pStdFrag->getSymbol());
							if (doIt) {
								//std::cout << " " << pStdFrag->getSymbol() << "\n";
								//submolecule* pStdSubmol = pStdMolecule->addSubMolecule();
								//pStdSubmol->copy(submolList[t]);
								std::vector<atom*> lAtomList = submolList[t]->getAtomList();
								for (unsigned int r = 0; r < lAtomList.size(); r++) {
									stdAtom* pStdAtom = pStdFrag->getStdAtom(lAtomList[r]->getName());
									if (pStdAtom) {
										//std::cout << "   |" << lAtomList[r]->getName() << " | " <<  lAtomList[r]->getZcharge() << "\n";
										pStdAtom->atmCharge = lAtomList[r]->getZcharge();
										//groupCharge += lAtomList[r]->getZcharge();
										groupCharge += pStdAtom->atmCharge;
									}
								}
							}
						}
					}
					/*
					 if (bStdMolecule) {
					 std::map<int, Bond*> bondMap = pSelMol->getBondMap();
					 if (!bondMap.empty()) {
					 for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
					 Bond* pBond = b->second;
					 atom* pBdAtom1 = pBond->atom1;
					 atom* pBdAtom2 = pBond->atom2;
					 
					 atom* pNewAtom1 = pStdMolecule->getAtom(pBdAtom1->getFileID(), 0, 1, 0);
					 atom* pNewAtom2 = pStdMolecule->getAtom(pBdAtom2->getFileID(), 0, 1, 0);
					 if (pNewAtom1 and pNewAtom2) {
					 pStdMolecule->addBond(pNewAtom1, pNewAtom2, pBond->type,
					 pBond->stereo, pBond->topology, pBond->size);
					 }
					 }
					 }
					 }
					 */
				}
				updateRespChargesMessage += " Charge: " + double2String(groupCharge, 2);
				//std::cout << " Group charge = " << double2String(groupCharge, 2) << "\n";
				MTKpp::errorLogger.throwError("MCPB::updateRespCharges",
							      updateRespChargesMessage, INFO);
			}
		}

		else if (inputFileContents[i][0] == "respgen") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: respgen

			 Description: Create resp1 and resp2 files (recreates functionality of respgen)

			 syntax: respgen /COL/MOL Group fileName
			 \endcode
			 */

			if ((inputFileContents[i].size() != 4) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::respgen",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::respgen",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}

				stdLibrary* pStdLibrary = pCollection->getStdLibrary();
				stdGroup* pStdGroup = 0;
				if (pStdLibrary) {
					pStdGroup = pStdLibrary->getStdGroup(inputFileContents[i][2]);
					if (!pStdGroup) {
						MTKpp::errorLogger.throwError("MCPB::respgen",
						 " Can't find stdGroup ... exiting ", MTK_ERROR);
						exit(1);
					}
				}
				else {
					MTKpp::errorLogger.throwError("MCPB::respgen",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
				}

				if (pSel->getSelectionType() == 1) {
				  molecule* pSelMol = pSel->getMol();
					if (!pSelMol) {
					  std::cout << " Error in respgen selection \n";
					  exit(1);
					}

				std::vector<std::string>::iterator result;

				acParser* pAcParser = new acParser();
				pAcParser->Write(inputFileContents[i][3]+".ac", pCollection);
				delete pAcParser;

				std::string bb0FileResp1 = inputFileContents[i][3] + "_bb0.resp1";
				std::string bb0FileResp2 = inputFileContents[i][3] + "_bb0.resp2";
				std::ofstream obb0FileResp1;
				std::ofstream obb0FileResp2;
				obb0FileResp1.open(bb0FileResp1.c_str());
				obb0FileResp2.open(bb0FileResp2.c_str());

				std::string bb1FileResp1 = inputFileContents[i][3] + "_bb1.resp1";
				std::string bb1FileResp2 = inputFileContents[i][3] + "_bb1.resp2";
				std::ofstream obb1FileResp1;
				std::ofstream obb1FileResp2;
				obb1FileResp1.open(bb1FileResp1.c_str());
				obb1FileResp2.open(bb1FileResp2.c_str());

				std::string bb2FileResp1 = inputFileContents[i][3] + "_bb2.resp1";
				std::string bb2FileResp2 = inputFileContents[i][3] + "_bb2.resp2";
				std::ofstream obb2FileResp1;
				std::ofstream obb2FileResp2;
				obb2FileResp1.open(bb2FileResp1.c_str());
				obb2FileResp2.open(bb2FileResp2.c_str());

				std::string bb3FileResp1 = inputFileContents[i][3] + "_bb3.resp1";
				std::string bb3FileResp2 = inputFileContents[i][3] + "_bb3.resp2";
				std::ofstream obb3FileResp1;
				std::ofstream obb3FileResp2;
				obb3FileResp1.open(bb3FileResp1.c_str());
				obb3FileResp2.open(bb3FileResp2.c_str());

				if (!obb0FileResp1 or !obb0FileResp2 or !obb1FileResp1 or !obb1FileResp2 or
				    !obb2FileResp1 or !obb2FileResp2 or !obb3FileResp1 or !obb3FileResp2) {
					MTKpp::errorLogger.throwError("MCPB::respgen",
								      " Unable to open RESP files ... exiting ", MTK_ERROR);
					exit(1);
				}

				std::vector<molecule*> molList = pCollection->getMoleculeList();

				std::string terminalFrags = "";
				char temp[100];
				bool bFirst = true;
				for (unsigned int m = 0; m < molList.size(); m++) {
					molecule* pMol = molList[m];
					std::vector<submolecule*> submols = pMol->getSubMoleculeList();
					for (unsigned int s = 0; s < submols.size(); s++) {
						if (submols[s]->hasStdFrag()) {
							stdFrag* pStdFrag = submols[s]->getStdFrag();
							if (pStdFrag->getParent()->getName() != inputFileContents[i][2]) {
								if (bFirst) {
									sprintf(temp,"%5d%10.5f\n", pStdFrag->numStdAtoms(), 0.0);
									bFirst = false;
								}
								else {
									sprintf(temp,"\n%5d%10.5f\n", pStdFrag->numStdAtoms(), 0.0);
								}
								
								terminalFrags += std::string(temp);
								std::vector<atom*> atomList = submols[s]->getAtomList();
								char temp2[80];
								int counter = 0;
								for (unsigned int a = 0; a < atomList.size(); a++) {
									if (counter > 7) { // was 6
										terminalFrags += "\n";
										counter = 0;
									}
									sprintf(temp2,"%5d%5d", 1, atomList[a]->getIndex());
									terminalFrags += std::string(temp2);
									counter++;
								}
								//terminalFrags += "\n";
							}
						}
					}
				}
				//std::cout << terminalFrags << std::endl;

				std::string bb0 = terminalFrags + "\n";
				std::string bb1 = terminalFrags + "\n";
				std::string bb2 = terminalFrags + "\n";
				std::string bb3 = terminalFrags + "\n";

				int lInd = 0;
				for (unsigned int m = 0; m < molList.size(); m++) {
					molecule* pMol = molList[m];

					if (pMol->getName() == "Reference") continue;

					std::vector<submolecule*> submols = pMol->getSubMoleculeList();
					for (unsigned int s = 0; s < submols.size(); s++) {
						if (submols[s]->hasStdFrag()) {
							stdFrag* pStdFrag = submols[s]->getStdFrag();
							//std::cout << " respgen " <<  pStdFrag->getParent()->getName() << " " << inputFileContents[i][2] << " \n ";
							if (pStdFrag->getParent()->getName() != inputFileContents[i][2]) continue;
							std::vector<atom*> subAtoms = submols[s]->getAtomList();
							for (unsigned int a = 0; a < subAtoms.size(); a++) {
								stdAtom* pStdAtom_a = subAtoms[a]->getStdAtom();
								lInd++;
								// iBB == 1
								if (pStdAtom_a) {
									std::string dd = pStdAtom_a->identity + ":" + pStdAtom_a->type;
		
									//result = std::find(bb_heavy.begin(), bb_heavy.end(), subAtoms[a]->getName());
									//result = std::find(bb_heavy.begin(), bb_heavy.end(), pStdAtom_a->identity);
									result = std::find(bb_heavy.begin(), bb_heavy.end(), dd);

									if (result != bb_heavy.end()) {
										sprintf(temp,"%5d%10.5f\n", 1, pStdAtom_a->atmCharge);
										//freezingAtoms += "\n BB Heavy "  + submols[s]->getName() + " " + int2String(submols[s]->getSubMolId()) + " " + dd;
										//std::cout << "\n BB Heavy "  + submols[s]->getName() + " " + int2String(submols[s]->getSubMolId()) + " " + dd;
										//std::cout << dd << " " << subAtoms[a]->getName() << " ---> " << subAtoms[a]->getStdAtom()->type << std::endl;
										//orespAdd << temp << std::endl;
										bb1 += std::string(temp);
										sprintf(temp,"%5d%5d\n", 1, subAtoms[a]->getIndex());
										//orespAdd << temp << std::endl;
										bb1 += std::string(temp);
									}

								// iBB == 2
									dd = pStdAtom_a->identity + ":" + pStdAtom_a->type;

									//result = std::find(bb.begin(), bb.end(), subAtoms[a]->getName());
									result = std::find(bb.begin(), bb.end(), dd);

									if (result != bb.end()) {
										//freezingAtoms += "\n BB "  + submols[s]->getName() + " " + int2String(submols[s]->getSubMolId()) + " " + dd;
										sprintf(temp,"%5d%10.5f\n", 1, pStdAtom_a->atmCharge);
										//orespAdd << temp << std::endl;
										bb2 += std::string(temp);
										sprintf(temp,"%5d%5d\n", 1, subAtoms[a]->getIndex());
										//orespAdd << temp << std::endl;
										bb2 += std::string(temp);
									}
								// iBB == 3 
									dd = pStdAtom_a->identity + ":" + pStdAtom_a->type;

									//result = std::find(bbb.begin(), bbb.end(), subAtoms[a]->getName());
									result = std::find(bbb.begin(), bbb.end(), dd);

									if (result != bbb.end()) {
										//freezingAtoms += "\n BBB "  + submols[s]->getName() + " "
										//              + int2String(submols[s]->getSubMolId()) + " " + dd;
										sprintf(temp,"%5d%10.5f\n", 1, pStdAtom_a->atmCharge);
										bb3 += std::string(temp);
										//orespAdd << temp << std::endl;
										sprintf(temp,"%5d%5d\n", 1, subAtoms[a]->getIndex());
										bb3 += std::string(temp);
										//orespAdd << temp << std::endl;
									}
								}
								else {
									MTKpp::errorLogger.throwError("MCPB::respgen",
										" Unknown option ... exiting ", MTK_ERROR);
								}
							}
						}
					}
				}

				std::string resp1 = "";
				resp1 += "Resp charges for organic molecule\n";
				resp1 += "\n";
				resp1 += " &cntrl\n";
				resp1 += "\n";
				resp1 += " nmol = 1,\n";
				resp1 += " ihfree = 1,\n";
				resp1 += " ioutopt = 1,\n";
				resp1 += "\n";
				resp1 += " &end\n";
				resp1 += "    1.0\n";
				resp1 += "Resp charges for organic molecule\n";

				//std::cout << pStdGroup->getCharge() << " " << pSelMol->getNumAtoms() << std::endl;
				sprintf(temp,"%5d%5d\n", int(pStdGroup->getCharge()),  pSelMol->getNumAtoms());
				resp1 += std::string(temp);

				std::vector<submolecule*> sList = pSelMol->getSubMoleculeList();
				typedef std::map<std::string, int>::iterator mapIterator;

				std::map<std::string, int> lMap;
				int ind = 0;
				for (unsigned int p = 0; p < sList.size(); p++) {
					std::vector<atom*> aList = sList[p]->getAtomList();
					for (unsigned int a = 0; a < aList.size(); a++) {
						if (!aList[a]->getStdAtom()) {
							std::cout << " Error ... exiting()\n";
							exit(1);
						}
						int elNum = aList[a]->getAtomicNum();
						std::string key = "";
						if (elNum == 1) {
							std::vector<atom*> hBondedAtoms = aList[a]->getBondedAtoms();
							if (hBondedAtoms.size() == 1) {
								key = aList[a]->getParent()->getName() + ":" + hBondedAtoms[0]->getName() + ":" + aList[a]->getStdAtom()->type;
							}
							else {
								key = aList[a]->getParent()->getName() + ":" + aList[a]->getStdAtom()->type;
							}
						}
						else {
							key = aList[a]->getParent()->getName() + ":" + aList[a]->getName();
						}

						mapIterator mi = lMap.find(key);

						int fg = 0;
						if (mi != lMap.end()){
							fg = lMap[key]+1;
						}
						else {
							lMap[key] = ind;
						}
						//sprintf(temp,"%4d %-11s%5d%5d\n", ind+1, key.c_str(), aList[a]->getAtomicNum(), fg);
						sprintf(temp,"%5d%5d\n", aList[a]->getAtomicNum(), fg);
						resp1 += std::string(temp);
						ind++;
					}
				}

				std::string resp2 = "";

				resp2 += "Resp charges for organic molecule\n";
				resp2 += "\n";
				resp2 +=  " &cntrl\n";
				resp2 +=  "\n";
				resp2 +=  " nmol = 1,\n";
				resp2 +=  " ihfree = 1,\n";
				resp2 +=  " ioutopt = 1,\n";
				resp2 +=  " iqopt = 2,\n";
				resp2 +=  " qwt = 0.001,\n";
				resp2 +=  "\n";
				resp2 +=  " &end\n";
				resp2 +=  "    1.0\n";
				resp2 +=  "Resp charges for organic molecule\n";

				sprintf(temp,"%5d%5d\n", int(pStdGroup->getCharge()),  pSelMol->getNumAtoms());
				resp2 += std::string(temp);

				ind = 0;
				std::map<std::string, int> lMap2;

				for (unsigned int p = 0; p < sList.size(); p++) {
					std::vector<atom*> aList = sList[p]->getAtomList();
					for (unsigned int a = 0; a < aList.size(); a++) {

						if (!aList[a]->getStdAtom()) {
							std::cout << " Error ... exiting()\n";
							exit(1);
						}
						int elNum = aList[a]->getAtomicNum();
						std::string key = "";
						if (elNum == 1) {
							std::vector<atom*> hBondedAtoms = aList[a]->getBondedAtoms();
							if (hBondedAtoms.size() == 1) {
								key = aList[a]->getParent()->getName() + ":" + hBondedAtoms[0]->getName() + ":" + 
								 aList[a]->getStdAtom()->type;
							}
 							else {
								key = aList[a]->getParent()->getName() + ":" + aList[a]->getStdAtom()->type;
							}
						}
						else {
							key = aList[a]->getParent()->getName() + ":" + aList[a]->getName();
						}

						mapIterator mi = lMap2.find(key);

						int fg = 0;
						if (mi != lMap2.end()){
							fg = lMap2[key]+1;
						}
						else {
							lMap2[key] = ind;
						}
						std::string atName = aList[a]->getName();
						std::string atType = aList[a]->getStdAtom()->type;
						int atNum = aList[a]->getAtomicNum();

						if (atName == " CH3" or atName == " CB " or containsSubStr(key, " CB :H1") or
							containsSubStr(key, " CB :HC") or 
							containsSubStr(key, " CH3:HC") or containsSubStr(key, " CH3:H1") or
							containsSubStr(key, "GLY: CA ")) {// or (atNum > 16)) {
							//if (atName == " CH3" or atName == " CB " or  atType == "H1" or  atType == "HC") {
							//sprintf(temp,"%4d %-11s%5d%5d\n", ind+1, key.c_str(), aList[a]->getAtomicNum(), fg);
							sprintf(temp,"%5d%5d\n", aList[a]->getAtomicNum(), fg);
						}
						else {
							//sprintf(temp,"%4d %-11s%5d%5d\n", ind+1, key.c_str(), aList[a]->getAtomicNum(), -99);
							sprintf(temp,"%5d%5d\n", aList[a]->getAtomicNum(), -99);
						}
						resp2 += std::string(temp);

						ind++;
					}
				}

				std::string twoBlankLine = "\n\n";
				bb0 += twoBlankLine;
				bb1 += twoBlankLine;
				bb2 += twoBlankLine;
				bb3 += twoBlankLine;

				obb0FileResp1 << resp1;
				obb0FileResp1 << bb0;
				obb0FileResp2 << resp2;
				obb0FileResp2 << bb0;

				obb1FileResp1 << resp1;
				obb1FileResp1 << bb1;
				obb1FileResp2 << resp2;
				obb1FileResp2 << bb1;

				obb2FileResp1 << resp1;
				obb2FileResp1 << bb2;
				obb2FileResp2 << resp2;
				obb2FileResp2 << bb2;

				obb3FileResp1 << resp1;
				obb3FileResp1 << bb3;
				obb3FileResp2 << resp2;
				obb3FileResp2 << bb3;

				obb0FileResp1.close();
				obb0FileResp2.close();
				obb1FileResp1.close();
				obb1FileResp2.close();
				obb2FileResp1.close();
				obb2FileResp2.close();
				obb3FileResp1.close();
				obb3FileResp2.close();
				}
			}
		}

		else if (inputFileContents[i][0] == "addStdMol") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: addStdMol
			 
			 Description:
			 
			 syntax: addStdMol /COL/MOL Group
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::addStdMol",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::addStdMol",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				stdLibrary* pStdLibrary = pCollection->getStdLibrary();
				stdGroup* pStdGroup = 0;
				if (pStdLibrary) {
					pStdGroup = pStdLibrary->getStdGroup(inputFileContents[i][2]);
					if (!pStdGroup) {
						MTKpp::errorLogger.throwError("MCPB::addStdMol",
									      " Can't find stdGroup ... exiting ", MTK_ERROR);
						exit(1);
					}
				}
				else {
					MTKpp::errorLogger.throwError("MCPB::addStdMol",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
				}
				
				if (pSel->getSelectionType() == 1) {
					molecule* pSelMol = pSel->getMol();
					if (!pSelMol) {
						exit(1);
					}
					pStdGroup->setStdMolecule(pSelMol);
				}
			}
		}
		
		else if (inputFileContents[i][0] == "assignStdFeatures") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: assignStdFeatures
			 
			 Description: Assigns std features to molecule
			 
			 syntax: assignStdFeatures /col/mol
			 \endcode
			 */
			if (inputFileContents[i].size() == 2) {
				
				selection* pSeln = new selection(pCollection);
				failure = pSeln->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::assignStdFeatures",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				if (pSeln->getSelectionType() == 1) {
					molecule* pLig = pSeln->getMol();
					if (!pLig) {
						MTKpp::errorLogger.throwError("MCPB::assignStdFeatures",
									      " Incorrect use ... exiting ", MTK_ERROR);
						exit(1);
					}
					atomTyper* pAtomTyper = new atomTyper();
					pAtomTyper->assignStdProperties(pLig);
					delete pAtomTyper;
				}
			}
		}
		
		else if (inputFileContents[i][0] == "findMetalCenters") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: findMetalCenters
			 
			 Description: Find all metal centers in the collection
			 
			 syntax: findMetalCenters
			 \endcode
			 */
			if ((inputFileContents[i].size() != 1) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::findMetalCenters",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			MTKpp::errorLogger.throwError("MCPB::findMetalCenters", " Start ", INFO);
			if (pCollection->hasMetal()) {
				pCollection->findMetals();
				pCollection->determineMetalEnvironments();
			}
			else {
				MTKpp::errorLogger.throwError("MCPB::findMetalCenters",
							      " Collection doesn't contain any metals ", INFO);
			}
		}
		
		else if (inputFileContents[i][0] == "writeLeap") {
			/*!
			 \code
			 Function: writeLeap
			 
			 Description: Write the metal center bonding info for leap
			 
			 syntax: writeLeap name pdbFile
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::writeLeap",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				std::string leapFile = inputFileContents[i][1] + ".leaprc";
				std::ofstream oleap;
				oleap.open(leapFile.c_str());
				
				if (!oleap) {
					std::cout << "\nUNABLE TO OPEN LEAP FILE"
					<< "\nFILENAME = " << leapFile << std::endl;
					exit(1);
				}
				
				std::string pmlFile = inputFileContents[i][1] + ".pml";
				std::ofstream opml;
				opml.open(pmlFile.c_str());
				
				if (!opml) {
					std::cout << "\nUNABLE TO OPEN PML FILE"
					<< "\nFILENAME = " << pmlFile << std::endl;
					exit(1);
				}
				
				int numberOfStdResidues = 0;
				std::vector<molecule*> mols = pCollection->getMoleculeList();
				for (unsigned int m = 0; m < mols.size(); m++) {
					if (mols[m]->getName() == "Reference") {
						numberOfStdResidues += mols[m]->getNumSubMolecules();
					}
				}
				
				oleap << "source leaprc.ff94" << std::endl;
				oleap << "frcmodMetals = loadamberparams metals.frcmod" << std::endl;
				oleap << "frcmod = loadamberparams " + inputFileContents[i][1] + ".frcmod" << std::endl;
				oleap << "loadamberprep " + inputFileContents[i][1] + ".prep" << std::endl;
				oleap << "model = loadpdb " + inputFileContents[i][2] << std::endl;
				
				opml << "load " << inputFileContents[i][2] << "\n";
				opml << "set valence=0" << "\n";
				opml << "hide lines" << "\n";
				opml << "set sphere_scale=0.4" << "\n";
				opml << "colorByElement" << "\n";
				opml << "bg white" << "\n";
				
				// Disulfide bonds
				for (unsigned int m = 0; m < mols.size(); m++) {
					std::map<int, Bond*> bondMap = mols[m]->getBondMap();
					for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
						Bond* pBond = b->second;
						if (pBond->kind == 3) {
							std::string name1 = pBond->atom1->getName();
							std::string name1Strip = stripString(name1, " ");
							std::string name2 = pBond->atom2->getName();
							std::string name2Strip = stripString(name2, " ");
							
							int colIndex1 = pBond->atom1->getParent()->getColIndex() - numberOfStdResidues;
							int colIndex2 = pBond->atom2->getParent()->getColIndex() - numberOfStdResidues;
							
							oleap << "bond model." + int2String(colIndex1)
							<< "." << name1Strip
							<< " model." + int2String(colIndex2)
							<< "." << name2Strip << std::endl;
						}
					}
				}
				
				std::vector<metalCenter*> mcs = pCollection->getMetalCenters();
				for (unsigned int m = 0; m < mcs.size(); m++) {
					atom* pMetalAtom = mcs[m]->getMetalAtom();
					std::string name = pMetalAtom->getName();
					std::string nameStrip = stripString(name, " ");
					
					int colIndex1 = pMetalAtom->getParent()->getColIndex() - numberOfStdResidues;
					
					std::string s = "bond model." + int2String(colIndex1) + "." + nameStrip + " ";
					
					//opml << "create metal" << m+1 << ", resi " << colIndex1 << "\n";
					
					opml << "create metal" << m+1 << ", resi " << pMetalAtom->getParent()->getSubMolId() << "\n";
					opml << "set sphere_scale=0.25, elem " << pMetalAtom->getElementSymbol() << "\n";
					opml << "show spheres, metal" << m+1 << "\n";
					
					std::vector<atom*> primaryBondedAtoms;
					mcs[m]->getPrimaryShellAtoms(primaryBondedAtoms);
					for (unsigned int n = 0; n < primaryBondedAtoms.size(); n++) {
						std::string name2 = primaryBondedAtoms[n]->getName();
						std::string name2Strip = stripString(name2, " ");
						
						int colIndex2 = primaryBondedAtoms[n]->getParent()->getColIndex() - numberOfStdResidues;
						
						oleap << s << "model." + int2String(colIndex2)
						<< "." << name2Strip << std::endl;
						
						//opml << "create r" << m+1 << n+1 << ", resi " << colIndex2 << "\n";
						opml << "create r" << m+1 << n+1 << ", resi " << primaryBondedAtoms[n]->getParent()->getSubMolId() << "\n";
						opml << "show sticks, r" << m+1 << n+1 << "\n";
						
						opml << "myLine2 bd" << m << n << ", metal" << m+1 << ", r" << m+1 << n+1 << " and name "
						<< primaryBondedAtoms[n]->getName() << ", r=0.0, g=0.0, b=0.0 \n";
					}
				}
				
				oleap << "solvateoct model TIP3PBOX 8.0" << std::endl;
				oleap << "saveamberparm model " + inputFileContents[i][1] + "_mcpb.top "
				+ inputFileContents[i][1] + "_mcpb.crd" << std::endl;
				oleap << "savePdb model " + inputFileContents[i][1] + "_mcpb_h2o.pdb" << std::endl;
				oleap << "quit" << std::endl;
				oleap.close();
				opml.close();
			}
		}
		
		else if (inputFileContents[i][0] == "writePrmtop") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writePrmtop
			 
			 Description: Write prmtop and coordinate files for sander
			 
			 syntax: writePrmtop inpcrd prmtop
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::writePrmtop",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				amberParser* pPrmTop = new amberParser();
				pPrmTop->Write(inputFileContents[i][1], inputFileContents[i][2], pCollection);
				delete pPrmTop;
			}
		}
		
		else if (inputFileContents[i][0] == "readNMode") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: readNMode
			 
			 Description: Read frequencies from nmode
			 
			 syntax: readNMode /Col/mol file_nmd2.out
			 \endcode
			 */
			if ((inputFileContents[i].size() != 3) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::readNMode",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				std::ifstream inmde;
				inmde.open(inputFileContents[i][2].c_str());
				
				if (!inmde) {
					MTKpp::errorLogger.throwError("MCPB::readNMode",
								      " Can't open nmode output file ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				std::vector<std::string> words;
				std::vector<double> frequencies;
				std::string freqString;
				std::string nmde_fileline;
				while (inmde) {
					getline(inmde, nmde_fileline);
					if (containsSubStr(nmde_fileline, "vibrational  ")) {
						while (nmde_fileline != "") {
							getline(inmde, nmde_fileline);
							if (nmde_fileline.size() == 0) {
								break;
							}
							words.clear();
							splitString(nmde_fileline, " ", words, 0);
							if (string2Int(words[0]) > 6) {
								frequencies.push_back(string2Double(words[1]));
								freqString += (words[1] + "|");
							}
						}
					}
				}
				std::cout << " " << frequencies.size() << " were read. " << std::endl;
				
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::readNMode",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				if (pSel->getSelectionType() == 1) {
					molecule* selMol = pSel->getMol();
					if (selMol) {
						int nAtoms = selMol->getNumAtoms();
						int nFreqs = 3*nAtoms - 6;
						if (static_cast<int>(frequencies.size()) != nFreqs) {
							errorMessage = " Number of nmode frequencies = " + int2String(frequencies.size())
							+ " doesn't match 3N-6 " + int2String(nFreqs);
							MTKpp::errorLogger.throwError("MCPB::readNMode",
										      errorMessage, MTK_ERROR);
							exit(1);
						}
						std::cout << freqString << std::endl;
						selMol->addProperty("frequencies", freqString);
					}
				}
				else {
					MTKpp::errorLogger.throwError("MCPB::readNMode",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
				}
			}
		}
		
		else if (inputFileContents[i][0] == "readNModeVectors") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: readNModeVectors
			 
			 Description: Read eigenvalues and eigenvectors from nmode vecs file
			 
			 syntax: readNModeVectors /Col/mol file_nmd2.vecs file_nmd2.molden
			 \endcode
			 */
			if ((inputFileContents[i].size() != 4) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::readNModeVectors",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				/* //boost/eigen
				 table<double>* pNModeEValues = pSheet->addTable();
				 if (!pNModeEValues) exit(1);
				 pNModeEValues->setName("NMode EigenValues");
				 
				 table<double>* pNModeEVectors = pSheet->addTable();
				 if (!pNModeEVectors) exit(1);
				 pNModeEVectors->setName("NMode EigenVectors");
				 
				 table<double>* pNModeCoords = pSheet->addTable();
				 if (!pNModeCoords) exit(1);
				 pNModeCoords->setName("NMode Coordinates");
				 
				 ublas::vector<double> coordinates;
				 ublas::vector<double> eigenvalues;
				 ublas::matrix<double, ublas::column_major> eigenvectors;
				 
				 std::ifstream inmde;
				 inmde.open(inputFileContents[i][2].c_str());
				 
				 if (!inmde) {
				 MTKpp::errorLogger.throwError("MCPB::readNModeVectors",
				 " Can't open nmode output file ... exiting ", MTK_ERROR);
				 exit(1);
				 }
				 
				 std::vector<std::string> words;
				 std::string nmde_fileline;
				 int nCs = 0;
				 int nFs = 0;
				 int nLs = 0;
				 while (inmde) {
				 getline(inmde, nmde_fileline);
				 if (containsSubStr(nmde_fileline, "NORMAL COORDINATE FILE")) {
				 getline(inmde, nmde_fileline);
				 words.clear();
				 splitString(nmde_fileline, " ", words, 0);
				 nCs = string2Int(words[0]);
				 nFs = nCs - 6;
				 if (nFs == 0) nFs++;
				 nLs = ceil(double(nCs / 7.0));
				 if (nLs == 0) nLs++;
				 
				 //std::cout << nCs << " " << nLs << std::endl;
				 coordinates.resize(nCs);
				 eigenvalues.resize(nFs);
				 eigenvectors.resize(nFs, nCs);
				 
				 pNModeEValues->setSizes(nFs, 1);
				 pNModeEVectors->setSizes(nFs, nCs);
				 pNModeCoords->setSizes(1,nCs);
				 
				 int coo = 0;
				 for (int c = 0; c < nLs; c++) {
				 getline(inmde, nmde_fileline);
				 words.clear();
				 splitString(nmde_fileline, " ", words, 0);
				 for (unsigned int w = 0; w < words.size(); w++) {
				 coordinates(coo) = string2Double(words[w]);
				 coo++;
				 }
				 }
				 }
				 
				 if (containsSubStr(nmde_fileline, "****")) {
				 getline(inmde, nmde_fileline);
				 words.clear();
				 splitString(nmde_fileline, " ", words, 0);
				 int cFreq = string2Int(words[0]);
				 if (cFreq <= nCs) {
				 if (cFreq > 6) {
				 double freq = string2Double(words[1]);
				 eigenvalues(cFreq-7) = freq;
				 pNModeEValues->setCellValue(cFreq-7, 0, freq);
				 
				 //std::cout << nmde_fileline << "| |" << freq << "| |" << eigenvalues(cFreq-7)<< std::endl;
				 }
				 }
				 
				 int freqCounter = 0;
				 for (int c = 0; c < nLs; c++) {
				 getline(inmde, nmde_fileline);
				 if (cFreq <= nCs) {
				 if (cFreq > 6) {
				 words.clear();
				 splitString(nmde_fileline, " ", words, 0);
				 for (unsigned int w = 0; w < words.size(); w++) {
				 double evec = string2Double(words[w]);
				 eigenvectors(cFreq-7, freqCounter) = evec;
				 pNModeEVectors->setCellValue(cFreq-7, freqCounter, evec);
				 freqCounter++;
				 }
				 }
				 }
				 }
				 }
				 }
				 
				 selection* pSel = new selection(pCollection);
				 failure = pSel->parse(inputFileContents[i][1]);
				 if (failure) {
				 MTKpp::errorLogger.throwError("MCPB::readNMode",
				 " Error in selection parsing ... exiting ", MTK_ERROR);
				 exit(1);
				 }
				 
				 std::ofstream omolden;
				 if (pSel->getSelectionType() == 1) {
				 molecule* selMol = pSel->getMol();
				 if (selMol) {
				 
				 omolden.open(inputFileContents[i][3].c_str());
				 if (!omolden) {
				 std::cout << "\nUNABLE TO OPEN MOLDEN FILE"
				 << "\nFILENAME = " << inputFileContents[i][3] << std::endl;
				 exit(1);
				 }
				 
				 omolden << "[Molden Format]" << std::endl;
				 omolden << "[GEOMETRIES] XYZ" << std::endl;
				 std::vector<atom*> as = selMol->getAtomList();
				 for (unsigned int sa = 0; sa < as.size(); sa++) {
				 omolden << as[sa]->getElementSymbol() << " ";
				 omolden << std::showpoint << as[sa]->getX() << " ";
				 omolden << as[sa]->getY() << " ";
				 omolden << as[sa]->getZ() << std::endl;
				 }
				 
				 omolden << "[FR-COORD]" << std::endl;
				 int m = 0;
				 for (unsigned int sa = 0; sa < as.size(); sa++) {
				 omolden << as[sa]->getElementSymbol() << " ";
				 for (int m2 = m; m2 < m+3; m2++) {
				 omolden << std::showpoint << coordinates(m2)/ANG2BOHR << " ";
				 }
				 m+=3;
				 omolden << " " << std::endl;
				 }
				 }
				 }
				 else {
				 MTKpp::errorLogger.throwError("MCPB::readNMode",
				 " Error in selection parsing ... exiting ", MTK_ERROR);
				 }
				 
				 omolden << "[FREQ]" << std::endl;
				 for (int f = 0; f < nFs; f++) {
				 omolden << eigenvalues(f) << std::endl;
				 }
				 
				 omolden << "[FR-NORM-COORD]" << std::endl;
				 for (int f = 0; f < nFs; f++) {
				 omolden << "vibration " << f+1 << std::endl;
				 int co = 0;
				 for (int f2 = 0; f2 < nCs; f2++) {
				 co++;
				 omolden << std::showpoint << eigenvectors(f,f2)/ANG2BOHR << " ";
				 if (co == 3) {
				 omolden << " " << std::endl;
				 co = 0;
				 }
				 }
				 }
				 omolden.close();
				 */
			}
		}
		
		else if (inputFileContents[i][0] == "updateFrequencies") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: updateFrequencies
			 
			 Description:
			 
			 syntax: updateFrequencies /Col/mol 0.9806
			 \endcode
			 */
			if ((inputFileContents[i].size() >= 2) and (pGParser)) {
				double scaleFactor = 1.0;
				if (inputFileContents[i].size() == 3) {
					scaleFactor = strtod(inputFileContents[i][2].c_str(), 0);
				}
				std::vector<double> freqs = pGParser->getFrequencies();
				
				selection* pSel = new selection(pCollection);
				failure = pSel->parse(inputFileContents[i][1]);
				if (failure) {
					MTKpp::errorLogger.throwError("MCPB::updateFrequencies",
								      " Error in selection parsing ... exiting ", MTK_ERROR);
					exit(1);
				}
				
				if (pSel->getSelectionType() == 1) {
					molecule* selMol = pSel->getMol();
					if (selMol) {
						std::string freqString;
						for (unsigned int f = 0; f < freqs.size(); f++) {
							freqString += (double2String(freqs[f] * scaleFactor) + "|");
							std::cout << freqs[f] << std::endl;
						}
						selMol->addProperty("frequencies", freqString);
					}
				}
			}
			else {
				MTKpp::errorLogger.throwError("MCPB::updateFrequencies",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
		}
		
		else if (inputFileContents[i][0] == "writeState") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writeState
			 
			 Description: Write state xml file
			 
			 syntax: writeState file.mtk
			 \endcode
			 */
			if ((inputFileContents[i].size() != 2) or (!pCollection)) {
				MTKpp::errorLogger.throwError("MCPB::writeState",
							      " Incorrect use ... exiting ", MTK_ERROR);
				exit(1);
			}
			else {
				mtkppParser* pMtkppParser = 0;
				pMtkppParser = new mtkppParser();
				pMtkppParser->Write(inputFileContents[i][1], pCollection);
				delete pMtkppParser;
			}
		}
		
		else if (inputFileContents[i][0] == "writeData") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: writeData
			 
			 Description: Write contents
			 
			 syntax: writeData f.xml
			 \endcode
			 */
			if ((inputFileContents[i].size() >= 2) and (pSheet)) {
				pDMParser->write(pSheet, inputFileContents[i][1]);
			}
		}
		
		else if (inputFileContents[i][0] == "nmodeMatch") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: nmodeMatch
			 
			 Description: Compare nMode and Gaussian Normal Modes
			 
			 syntax: nmodeMatch filename.xml
			 \endcode
			 */
			if ((inputFileContents[i].size() >= 2)) {
				pDMParser->read(pSheet, inputFileContents[i][1]);
				/* //boost/eigen
				 table<double>* nmodeEValues = pSheet->getTable("NMode EigenValues");
				 ublas::matrix<double> &nmodeEValuesMatrix = nmodeEValues->getMatrix();
				 
				 table<double>* nmodeEVectors = pSheet->getTable("NMode EigenVectors");
				 ublas::matrix<double> &nmodeEVectorsMatrixT = nmodeEVectors->getMatrix();
				 
				 table<double>* gaussianEValues = pSheet->getTable("Gaussian EigenValues");
				 ublas::matrix<double> &gaussianEValuesMatrix = gaussianEValues->getMatrix();
				 
				 table<double>* gaussianEVectors = pSheet->getTable("Gaussian EigenVectors");
				 ublas::matrix<double> &gaussianEVectorsMatrixT = gaussianEVectors->getMatrix();
				 
				 int nRowsNMode = nmodeEValues->getNumRows();
				 int nRowsGaussian = gaussianEValues->getNumRows();
				 
				 int nColsNMode = nmodeEVectors->getNumColumns();
				 int nColsGaussian = gaussianEVectors->getNumColumns();
				 
				 if (nRowsNMode != nRowsGaussian and nColsNMode != nColsGaussian) {
				 std::cout << " PROBLEM " << std::endl;
				 }
				 else {
				 std::cout << " EIGENVALUES " << std::endl;
				 for (int z = 0; z < nRowsNMode; z++) {
				 std::cout << nmodeEValuesMatrix(z, 0) << " " << gaussianEValuesMatrix(z, 0) << std::endl;
				 }
				 
				 //std::cout << " NMODE " << std::endl;
				 //nmodeEVectors->printMatrix();
				 
				 //std::cout << " Gaussian " << std::endl;
				 //gaussianEVectors->printMatrix();
				 //std::cout << " \n\n " << std::endl;
				 
				 // Normalize Eigenvector matrices
				 for (int zi = 0; zi < nRowsNMode; zi++) {
				 double normNMode = 0.0;
				 double normGaussian = 0.0;
				 for (int zk = 0; zk < nColsNMode; zk++) {
				 normNMode += (nmodeEVectorsMatrixT(zi, zk) * nmodeEVectorsMatrixT(zi, zk));
				 normGaussian += (gaussianEVectorsMatrixT(zi, zk) * gaussianEVectorsMatrixT(zi, zk));
				 }
				 normNMode = sqrt(normNMode);
				 normGaussian = sqrt(normGaussian);
				 for (int zk = 0; zk < nColsNMode; zk++) {
				 //std::cout << gaussianEVectorsMatrixT(zi, zk) << " " << normGaussian << " "
				 //          << gaussianEVectorsMatrixT(zi, zk) / normGaussian << " ";
				 double newValue1 = nmodeEVectorsMatrixT(zi, zk) / normNMode;
				 double newValue2 = gaussianEVectorsMatrixT(zi, zk) / normGaussian;
				 nmodeEVectors->setCellValue(zi, zk, newValue1);
				 gaussianEVectors->setCellValue(zi, zk, newValue2);
				 //std::cout << gaussianEVectors->getCellValue(zi, zk) << std::endl;
				 }
				 }
				 
				 ublas::matrix<double> &nmodeEVectorsMatrix = nmodeEVectors->getMatrix();
				 ublas::matrix<double> &gaussianEVectorsMatrix = gaussianEVectors->getMatrix();
				 
				 //std::cout << " NMODE " << std::endl;
				 //nmodeEVectors->printMatrix();
				 
				 //std::cout << " Gaussian " << std::endl;
				 //gaussianEVectors->printMatrix();
				 //std::cout << " \n\n " << std::endl;
				 
				 // dot product
				 double dotProduct = 0.0;
				 double dotProductMax = 0.0;
				 int dotProductMaxIndex = 0;
				 double eN = 0.0;
				 double eN2 = 0.0;
				 double dff = 0.0;
				 for (int zi = 0; zi < nRowsNMode; zi++) { // NMode
				 for (int zj = 0; zj < nRowsNMode; zj++) { // Gaussian
				 for (int zk = 0; zk < nColsNMode; zk++) { // Loop over coordinates in Eigenvector
				 //std::cout << nmodeEVectorsMatrix(zi, zk) << " " <<  gaussianEVectorsMatrix(zj, zk) << " ";
				 eN += (nmodeEVectorsMatrix(zi, zk) * nmodeEVectorsMatrix(zi, zk));
				 eN2 += (gaussianEVectorsMatrix(zj, zk) * gaussianEVectorsMatrix(zj, zk));
				 dotProduct += (nmodeEVectorsMatrix(zi, zk) * gaussianEVectorsMatrix(zj, zk));
				 }
				 //std::cout << " eN = " << sqrt(eN) << " eN2 = " << sqrt(eN2) << " ";
				 std::cout <<  zi << " - " << zj << " " << dotProduct << " " << dotProductMax << std::endl;
				 
				 eN = 0.0;
				 eN2 = 0.0;
				 dotProduct = std::abs(dotProduct);
				 if (dotProduct > dotProductMax) {
				 dotProductMax = dotProduct;
				 dotProductMaxIndex = zj;
				 }
				 dotProduct = 0.0;
				 }
				 double dffL = nmodeEValuesMatrix(zi, 0) - gaussianEValuesMatrix(dotProductMaxIndex, 0);
				 dffL *= dffL;
				 std::cout << "        nmode index = " << zi << " BestMatch = " << dotProductMaxIndex
				 << " => dotProduct =" << dotProductMax << " EVALUES: "
				 << " " << nmodeEValuesMatrix(zi, 0) << " " << gaussianEValuesMatrix(dotProductMaxIndex, 0) << "  "
				 << dffL
				 << std::endl;
				 
				 dff += dffL;
				 dotProductMax = 0.0;
				 dotProductMaxIndex = 0;
				 }
				 std::cout << " FINAL DFF = " << sqrt(dff/nRowsNMode) << std::endl;
				 }
				 */
			}
		}
		
		else if (inputFileContents[i][0] == "renumber") {
			/*!
			 @ingroup MCPB_commands
			 \code
			 Function: renumber
			 
			 Description: Renumbers collection
			 
			 syntax: renumber
			 \endcode
			 */
			if (pCollection) {
				pCollection->renumber();
			}
		}
		
		else {
			std::string unknownCommand = " Unknown command: \"" + inputFileContents[i][0] + "\"";
			MTKpp::errorLogger.throwError("MCPB", unknownCommand, MTK_ERROR);
		}
	}
	
	// - Clean up - //
	delete pCollection;
	delete pPdbParser;
	delete pSheet;
	return 0;
}
