/*!
   \file DivConInterface.h
   \brief A bunch of routines to deal with DivCon
   \author Martin B. Peters
   \author Duane Williams

   \todo include dcParser ...

   $Date: 2010/03/29 20:33:22 $
   $Revision: 1.7 $

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

#ifndef DIVCONINTERFACE_H
#define DIVCONINTERFACE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>

namespace MTKpp
{

// ============================================================
// Class : DivConInterface()
// ------------------------------------------------------------
/*! 
   \class DivConInterface
   \brief A class to interface with DivCon
   \author Martin Peters
   \author Duane Williams
   \version 0.1
   \date 2005
*/
// ============================================================
class DivConInterface
{
public:

    /*!
       \brief DivConInterface Constructor
    */
    inline DivConInterface() {};

    /*!
       \brief run a DivCon calculation
    */
    inline int runDivCon(std::string divConInput) {
      std::string DivCon = getenv("DIVCON");

      int iDivCon = -1;
      std::string runDivCon = DivCon + " -i " + divConInput;
      iDivCon = system(runDivCon.c_str());
      if (iDivCon != 0 ) {
        std::cout << " DivConInterface:Error executing DivCon for " << divConInput << std::endl;
        return 1;
      }
      else {
        std::cout << " DivConInterface:Divcon Successful: " <<  divConInput << std::endl;
      }
      return 0;
    }

    /*!
       \brief Minimize all molecule in a collection
    */
    inline int Minimize(std::string hamiltonian, std::string chargeModel,
                 std::string optScheme, int maxOpt, collection* pCollection) {
      molecule* pMolecule = 0;
      dcParser* pDcParser = new dcParser();
      pDcParser->setHamiltonian(hamiltonian);
      pDcParser->setChargeModel(chargeModel);
      pDcParser->setOptimizer(optScheme);
      pDcParser->setMaxOpt(maxOpt);
      pDcParser->setDirect(1);

      std::vector<molecule*> molList = pCollection->getMoleculeList();
      std::string dcFileName = "";
      std::string dcOutFileName = "";

      for (unsigned int i = 0; i < molList.size(); i++) {
        pMolecule = molList[i];

        if (molList.size() > 1) {
          std::stringstream ss1;
          ss1 << i+1;
          std::string num = ss1.str().c_str(); 
          dcFileName = "dcMin_" + num + ".in";
          dcOutFileName = "dcMin_" + num + ".out";
        }
        else {
          dcFileName = "dcMin_1.in";
        }

        int nAtoms = pMolecule->getNumAtoms();
        if (nAtoms == 1) {
          pDcParser->setOptimizer("");
          pDcParser->setMaxOpt(0);
        }
        else {
          pDcParser->setOptimizer(optScheme);
          pDcParser->setMaxOpt(maxOpt);
        }

        pMolecule->setTotalCharge(pMolecule->getFormalCharge()); // pMolecule->setTotalCharge(0);
        bool success = true;
        pDcParser->Write(dcFileName, pMolecule, success);
        if (!success) return 1;

        int failure = runDivCon(dcFileName);
        if (failure) {
          delete pDcParser;
          return 1;
        }
      }

      // - Clean up - //
      delete pDcParser;
      return 0;
    }

    /*!
       \brief Minimize all molecule in a collection
    */
    inline int Minimize(std::string hamiltonian, std::string chargeModel,
                 std::string optScheme, int maxOpt, molecule* pMolecule, int num = 1) {
      dcParser* pDcParser = new dcParser();
      pDcParser->setHamiltonian(hamiltonian);
      pDcParser->setChargeModel(chargeModel);
      pDcParser->setOptimizer(optScheme);
      pDcParser->setMaxOpt(maxOpt);
      pDcParser->setDirect(1);

      std::string dcFileName = "";
      std::string dcOutFileName = "";

      std::stringstream ss1;
      ss1 << num;
      std::string numStr = ss1.str().c_str(); 
      dcFileName = "dcMin_" + numStr + ".in";
      dcOutFileName = "dcMin_" + numStr + ".out";

      int nAtoms = pMolecule->getNumAtoms();
      if (nAtoms == 1) {
        pDcParser->setOptimizer("");
        pDcParser->setMaxOpt(0);
      }
      else {
        pDcParser->setOptimizer(optScheme);
        pDcParser->setMaxOpt(maxOpt);
      }

      pMolecule->setTotalCharge(pMolecule->getFormalCharge()); // pMolecule->setTotalCharge(0);

      bool success = true;
      pDcParser->Write(dcFileName, pMolecule, success);
      if (!success) return 1;

      int failure = runDivCon(dcFileName);
      if (failure) {
        delete pDcParser;
        return 1;
      }

      // - Clean up - //
      delete pDcParser;
      return 0;
    }

    /*!
       \brief Run nmr on all molecules in a collection
    */
    inline int runNMR(std::string hamiltonian, int nuclei, collection* pCollection) {
      molecule* pMolecule = 0;
      dcParser* pDcParser = new dcParser();
      pDcParser->setHamiltonian(hamiltonian);
      pDcParser->setDirect(1);
      pDcParser->setCalNuc(nuclei);

      std::string dcFileName = "";

      std::vector<molecule*> molList = pCollection->getMoleculeList();

      for (unsigned int i = 0; i < molList.size(); i++) {
        pMolecule = molList[i];
        std::stringstream ss1;
        ss1 << i+1;
        std::string i_str = ss1.str().c_str(); 
        dcFileName = "dcMin_" + i_str + ".out";

        // test if file is there
        std::ifstream inputStream;
        inputStream.open(dcFileName.c_str());
        if (!inputStream) {
          std::cerr << "Error opening DivCon output file, please run minimization first" << std::endl;
          return 1;
        }

        pDcParser->Read(dcFileName, pMolecule);
        std::string dcFileName2 = "dcNmr_" + i_str + ".in";
        std::string dcOutFileName = "dcNmr_" + i_str + ".out";

        bool success = true;
        pDcParser->Write(dcFileName2, pMolecule, success);
        if (!success) return 1;

        int failure = runDivCon(dcFileName2);
        if (failure) {
          delete pDcParser;
          return 1;
        }
        pDcParser->Read(dcOutFileName, pMolecule);
      }

      // - Clean up - //
      delete pDcParser;
      return 0;
    }
    /*!
       \brief Run nmr on all molecules in a collection
    */
    inline int runNMR(dcParser* pDcParser, collection* pCollection) {
      molecule* pMolecule = 0;
      std::string dcFileName = "";

      std::vector<molecule*> molList = pCollection->getMoleculeList();

      for (unsigned int i = 0; i < molList.size(); i++) {
        pMolecule = molList[i];
        std::stringstream ss1;
        ss1 << i+1;
        std::string i_str = ss1.str().c_str(); 
        dcFileName = "dcMin_" + i_str + ".out";

        // test if file is there
        std::ifstream inputStream;
        inputStream.open(dcFileName.c_str());
        if (!inputStream) {
          std::cerr << "Error opening DivCon output file, please run minimization first" << std::endl;
          return 1;
        }

        pDcParser->Read(dcFileName, pMolecule);
        std::string dcFileName2 = "dcNmr_" + i_str + ".in";
        std::string dcOutFileName = "dcNmr_" + i_str + ".out";

        bool success = true;
        pDcParser->Write(dcFileName2, pMolecule, success);
        if (!success) return 1;

        int failure = runDivCon(dcFileName2);
        if (failure) {
          delete pDcParser;
          return 1;
        }
        pDcParser->Read(dcOutFileName, pMolecule);
      }
      return 0;
    }

    /*!
       \brief Run nmr on all molecules in a collection
    */
    inline int runNMR(std::string hamiltonian, int nuclei, molecule* pMolecule, int n) {
      dcParser* pDcParser = new dcParser();
      pDcParser->setHamiltonian(hamiltonian);
      pDcParser->setDirect(1);
      pDcParser->setCalNuc(nuclei);

      std::string dcFileName = "";

      std::stringstream ss1;
      ss1 << n;
      std::string i_str = ss1.str().c_str(); 
      dcFileName = "dcMin_" + i_str + ".out";

      // test if file is there
      std::ifstream inputStream;
      inputStream.open(dcFileName.c_str());
      if (!inputStream) {
        std::cerr << "Error opening DivCon output file, please run minimization first" << std::endl;
        return 1;
      }

      pDcParser->Read(dcFileName, pMolecule);
      std::string dcFileName2 = "dcNmr_" + i_str + ".in";
      std::string dcOutFileName = "dcNmr_" + i_str + ".out";

      bool success = true;
      pDcParser->Write(dcFileName2, pMolecule, success);
      if (!success) return 1;

      int failure = runDivCon(dcFileName2);
      if (failure) {
        delete pDcParser;
        return 1;
      }
      pDcParser->Read(dcOutFileName, pMolecule);

      // - Clean up - //
      delete pDcParser;
      return 0;
    }

protected:

   //!
};

} // MTKpp namespace

#endif // DIVCONINTERFACE_H

