/** cpin.cpp: Handles the code for parsing cpin files and defining titratable
  * residues.
  */

#include <iostream>

#include "constants.h"
#include "cpin.h"
#include "exceptions.h"
#include "string_manip.h"

using namespace std;

/// Cpin constructor from a char array
Cpin::Cpin() :
is_valid_(false)
{ }
   
int Cpin::Parse(const char* cpinname) {
   filename_ = string(cpinname);
   if (filename_.size() > FN_LEN)
      throw StringBufferOverflow("File name [[ " + filename_ + 
            " ]] too large. Adjust FN_LEN in cpin.h and parse_cpin.F90 and recompile");
   char *my_fname = (char*) filename_.c_str();
   char resname[TITR_RES_C+1][40];
   int ierr;
   parse_cpin_(&trescnt_, protcnt_, stateinf_, resname, my_fname, &ierr);

   // Error catch
   if (ierr != 0) {
      cerr << "Error: Could not open or parse " << filename_ << " for reading!" << endl;
      return ierr;
   }

   for (int i = 0; i <= trescnt_; i++) {
      resname_.push_back(string(resname[i]));
   }

   is_valid_ = true;

   // Now define all of the residues
   for (int i = 0; i < trescnt_; i++) {
      int nstates = stateinf_[i].num_states;
      int firststate = stateinf_[i].first_state;
      vector<int> res_protcnt;
      for (int j = 0; j < nstates; j++)
         res_protcnt.push_back( protcnt_[firststate+j] );
      residues_.push_back( TitratableResidue(resname_[i+1], res_protcnt) );
   }

   return 0;
}

int Cpin::Parse(string const& cpinname) {
   return Parse(cpinname.c_str());
}

/// Constructor for TitratableResidue
TitratableResidue::TitratableResidue(string const& resname,
                                     vector<int> const& protcnts) :
resid_(0)
{
   protcnts_ = protcnts;

   // Process the resname
   vector<string> words = split(resname);
   resname_ = words[1]; resid_ = StringToInt(words[2]);

   // Now determine which states are "protonated" and which are "deprotonated"
   int max_prots = protcnts_[0];

   for (size_t i = 1; i < protcnts_.size(); i++)
      max_prots = protcnts_[i] > max_prots ? protcnts_[i] : max_prots;

   for (size_t i = 0; i < protcnts_.size(); i++)
      protonated_.push_back(protcnts_[i] == max_prots);
}
