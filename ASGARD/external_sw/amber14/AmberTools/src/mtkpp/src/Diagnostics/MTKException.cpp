/* 
 * File:   MTKException.cpp
 * Author:
 * 
 * Created on May 21, 2009, 9:24 PM
 */

#include "MTKException.h"

namespace MTKpp
{

MTKException::MTKException() {}

MTKException::MTKException(std::string message) {
    this->message=message;
}

MTKException::MTKException(const MTKException& orig) {}

MTKException::~MTKException() throw() {}

} // MTKpp namespace
