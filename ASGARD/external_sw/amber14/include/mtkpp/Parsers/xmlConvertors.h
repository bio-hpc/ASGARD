/*!
   \file xmlConvertors.h
   \brief Xml handler

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.6 $

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

#ifndef XMLCONVERTORS_H
#define XMLCONVERTORS_H

#include <iostream>
#include <sstream>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNode.hpp>

#if (XERCES_VERSION_MAJOR == 3)
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMLSSerializer.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/impl/DOMLSSerializerImpl.hpp>
#else
#include <xercesc/dom/DOMWriter.hpp>
#include <xercesc/dom/DOM.hpp>
#endif

#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/framework/XMLFormatter.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
//using namespace std;

XERCES_CPP_NAMESPACE_USE

// ============================================================
// Class : StrX()
// ------------------------------------------------------------
/*! 
   \class StrX
   \brief This is a simple class that lets us do easy (though not terribly efficient) 
trancoding of XMLCh data to local code page for display.
   \version 0.1
   \date 2005
*/
// ============================================================
class StrX{
public :

    /*!
       \brief StrX Constructor
    */
    StrX(const XMLCh* const toTranscode){
        // Call the private transcoding method
        fLocalForm = XMLString::transcode(toTranscode);
    }

    //! StrX Destructor
    ~StrX(){
        XMLString::release(&fLocalForm);
    }
    const char* localForm() const{
        return fLocalForm;
    }
protected:

    //! ???
    char*   fLocalForm;
};

#define XC(xch) StrX(xch).localForm()
#define CX(str) XStrC(str).unicodeForm()
#define X(str) XStrC(str).unicodeForm()

// ============================================================
// Class : XStrC()
// ------------------------------------------------------------
/*! 
   \class XStrC
   \brief This is a simple class that lets us do easy (though not terribly efficient) 
trancoding of char* data to XMLCh data.
   \version 0.1
   \date 2005
*/
// ============================================================
class XStrC{
public :

    /*!
       \brief XStrC Constructor
    */
    XStrC(const char* const toTranscode){
        // Call the private transcoding method
        fUnicodeForm = XMLString::transcode(toTranscode);
    }

    //! XStrC destructor
    ~XStrC(){
        XMLString::release(&fUnicodeForm);
    }
    const XMLCh* unicodeForm() const{
        return fUnicodeForm;
    }

protected:
    //! ???
    XMLCh*   fUnicodeForm;
};

template<typename P>
const char* tocchars(P val){
std::stringstream sval;
sval << val << std::ends;
return (sval.str().c_str());
}

#endif // XMLCONVERTORS_H

