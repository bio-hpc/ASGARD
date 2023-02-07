/*!
   \file tester.cpp

   \brief A source file used for testing new code

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.6 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/printHeader.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"

// - MTK++ INCLUDE
#include "Utils/vector3d.h"
#include "Utils/constants.h"
#include "Utils/diagonalize.h"

#include "Statistics/sheet.h"
#include "Statistics/table.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/paramParser.h"

// - COMMAND LINE OPTIONS
#include "Parsers/commLineOptions.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;

void wait ( int seconds )
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}

/*!
   \brief File for testing new code
*/
int main (int argc, char **argv)
{
    std::string prog_name = "tester";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);
    printHeader(std::cout, prog_name, authors);

    std::string AMBERHOME = getenv("AMBERHOME");

    collection* pCollection = new collection();
    pCollection->addStdLibrary();
    pCollection->addParameters();

    // Read elements xml file
    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME + "/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);

    // Read parameter library
    paramParser* pParamParser = new paramParser(pCollection->getParameters());
    if (pParamParser) {
      std::string paramXmlFile = AMBERHOME + "/dat/mtkpp/parm94.xml";
      pParamParser->Read(paramXmlFile);
      paramXmlFile = AMBERHOME + "/dat/mtkpp/parm_gaff.xml";
      pParamParser->Read(paramXmlFile);
      //pParamParser->Write("asdf_param.xml", "parm94");
      delete pParamParser;
    }

    // Read standard library
    stdLibrary* pStdLibrary = pCollection->getStdLibrary();
    if (pStdLibrary) {
      stdLibParser* pStdLibParser = new stdLibParser(pStdLibrary, pCollection->getParameters());
      if (pStdLibParser) {
        std::string libXmlFile = AMBERHOME + "/dat/mtkpp/amino94.xml";
        pStdLibParser->Read(libXmlFile);
        //pStdLibParser->Write("asdf.xml");
        delete pStdLibParser;
      }
    }

    //
    // ADD CODE HERE
    //

    //
    // END CODE HERE
    //

    // - Clean up - //
    delete pElementParser;
    delete pCollection;
    return 0;
}

