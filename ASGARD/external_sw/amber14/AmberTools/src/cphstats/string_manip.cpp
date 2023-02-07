// string_manip.cpp: Contains the routines to do string manipulations

#include <sstream>
#include <cstring>
#include "constants.h"
#include "exceptions.h"
#include "string_manip.h"

using namespace std;

/// Splits a string along a specific delimiter
vector<string> split(char* instring, const char* delim) {
   vector<string> result;
   char* pch;
   pch = strtok(instring, delim);
   while (pch != NULL) {
      result.push_back( string(pch) );
      pch = strtok(NULL, delim);
   }
   return result;
}

vector<string> split(string const& instring, const char* delim) {
   char* buf = (char*)instring.c_str();
   return split(buf, delim);
}

vector<string> split(string const& instring, string const& delim) {
   char* buf = (char*)instring.c_str();
   return split(buf, delim.c_str());
}

vector<string> split(string const& instring) {
   char* buf = (char*)instring.c_str();
   return split(buf, " \n\r");
}

vector<string> split(char* instring) {
   return split(instring, " \r\n");
}

string strip(string const& instring) {
   int first = -1;
   int last = 0;
   for (size_t i = 0; i < instring.size(); i++) {
      if (instring[i] != ' ' && instring[i] != '\n' && instring[i] != '\r')
         last = i;
      if (first == -1 && (instring[i] == ' ' || instring[i] == '\n' ||
          instring[i] == '\r')) continue;
      first = i;
   }
   
   return instring.substr(first, last);
}

string strip(const char* instring) {
   return strip(string(instring));
}

string upper(string const& instring) {
   string ret;
   for (size_t i = 0; i < instring.size(); i++)
      ret += (char) toupper(instring[i]);

   return ret;
}

string upper(const char* instring) {
   return upper(string(instring));
}

int StringToInt(string const& s) {
   int result;
   if (!(istringstream(s) >> result))
      throw InvalidInteger("Could not convert [[ " + s + " ]] to an integer!");
   return result;
}

double StringToDouble(string const& s) {
   double result;
   if (!(istringstream(s) >> result))
      throw InvalidDecimal("Could not convert [[ " + s + " ]] to a decimal!");
   return result;
}

float StringToFloat(string const& s) {
   float result;
   if (!(istringstream(s) >> result))
      throw InvalidDecimal("Could not convert [[ " + s + " ]] to a decimal!");
   return result;
}
