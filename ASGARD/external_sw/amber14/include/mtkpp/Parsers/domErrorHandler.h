/*!
   \file domErrorHandler.h
   \brief Dom xml parser error handling

   $Date: 2010/03/29 20:39:34 $
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

#ifndef DOMERRORHANDLER_H
#define DOMERRORHANDLER_H

#include "xmlConvertors.h"
#include <xercesc/util/OutOfMemoryException.hpp>

// ============================================================
// Class : MyDOMErrorHandler()
// ------------------------------------------------------------
/*! 
   \class MyDOMErrorHandler
   \brief Dom xml parser error handling
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class MyDOMErrorHandler : public DOMErrorHandler{
public:

    /*!
       \brief MyDOMErrorHandler Constructor
    */
    MyDOMErrorHandler();

    //! MyDOMErrorHandler Destructor
    ~MyDOMErrorHandler();

    //! ???
    bool getSawErrors() const;

    //! ???
    bool handleError(const DOMError& domError);

    //! ???
    void resetErrors();

protected:

    //! ???
    MyDOMErrorHandler(const MyDOMErrorHandler&);

    //! ???
    void operator=(const MyDOMErrorHandler&);

    //! ???
    bool fSawErrors;
};


//!  DOMCountHandlers: Overrides of the DOM ErrorHandler interface
inline std::ostream& operator<<(std::ostream& target, const StrX& toDump){
    target << toDump.localForm();
    return target;
}

//! ???
inline bool MyDOMErrorHandler::getSawErrors() const{
    return fSawErrors;
}

#endif // DOMERRORHANDLER_H
