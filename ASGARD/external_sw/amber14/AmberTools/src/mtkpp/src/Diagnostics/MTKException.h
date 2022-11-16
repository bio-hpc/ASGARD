/* 
 * File:   MTKException.h
 * Author:
 *
 * Created on May 21, 2009, 9:24 PM
 */

#ifndef _MTKEXCEPTION_H
#define _MTKEXCEPTION_H

#include <string>

#include "boost/exception.hpp"

#ifdef __INTEL_COMPILER

// remark #869: parameter was never referenced
#pragma warning(disable:869)

#endif

namespace MTKpp
{
class MTKException : public boost::exception {
public:
    std::string message;
public:
    MTKException();
    MTKException(std::string message);
    MTKException(const MTKException& orig);
    virtual ~MTKException() throw();

protected:

};

} // MTKpp namespace

#endif /* _MTKEXCEPTION_H */
