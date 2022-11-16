// conprob.cpp: contains the code to handle conditional probability evaluations
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cctype>

#include "conprob.h"
#include "string_manip.h"

using namespace std;

ConditionalProb::ConditionalProb() :
valid_(UNASSIGNED) {}

ConditionalProb::ConditionalProb(string const& instring) {
   valid_ = UNASSIGNED;
   instring_ = instring;
}

ConditionalProb::ConditionalProb(const char* instring) {
   valid_ = UNASSIGNED;
   instring_ = string(instring);
}

ConditionalProb::RetType ConditionalProb::Set(Cpin const& cpin) {

   if (instring_.empty()) {
      cerr << "Error: Conditional probability string is empty!" << endl;
      return ERR;
   }

   // Build the active state list
   for (Cpin::ResIterator rit = cpin.begin(); rit != cpin.end(); rit++)
      active_states_.push_back( ActiveState(rit->numStates(), false) );

   // For comparisons
   string DEPROT = string("DEPROTONATED");
   string PROT = string("PROTONATED");

   // We have to parse something like; <number>:<state>,<number>:<state>,...
   vector<string> criteria = split(instring_.c_str(), ",");

   for (vector<string>::const_iterator cit = criteria.begin();
        cit != criteria.end(); cit++) {
      vector<string> parts = split(*cit, ":");
      if (parts.size() != 2) {
         cerr << "Error: Conditional probability string [" << instring_ <<
                 "] not recognized!" << endl;
         return ERR;
      }
      int resid = StringToInt(parts[0]);
      if (resid == 0) {
         cerr << "Error: Could not determine residue name from [" << *cit << "] in ["
              << instring_ << "]" << endl;
         return ERR;
      }
      // Find out which residue index this is
      int residx = -1;
      for (int i = 0; i < cpin.getTrescnt(); i++) {
         if (cpin.getResidues()[i].getResnum() == resid) {
            residx = i;
            break;
         }
      }
      if (residx < 0) {
         cerr << "Error: Residue " << resid << " does not appear to be titratable!" << endl;
         return ERR;
      }
      
      int stateidx;
      istringstream iss(parts[1]);
      if (parts[1].find(";") != string::npos) {
         /* We have specified multiple protonation states we want to match for
          * a single residue
          */
         vector<string> selres = split(parts[1], ";");
         for (vector<string>::const_iterator sit = selres.begin();
               sit != selres.end(); sit++) {
            int state = StringToInt(*sit);
            if (state >= cpin.getResidues()[residx].numStates()) {
               cerr << "Error: State index " << state << " out of range for residue "
                    << resid << "!" << endl;
               return ERR;
            }
            active_states_[residx][state] = true;
         }
     }else if (iss >> stateidx) {
         if (stateidx >= cpin.getResidues()[residx].numStates()) {
            cerr << "Error: State index " << stateidx << " out of range for residue "
                 << resid << "!" << endl;
            return ERR;
         }
         active_states_[residx][stateidx] = true;
//       cout << "Residue " << residx << ", state " << stateidx << " is true." << endl; // DEBUG
     }else {
         parts[1] = upper(parts[1]);
         if (PROT.size() >= parts[1].size() && 
             parts[1] == PROT.substr(0, parts[1].size())) {
            // Label all protonated states of this residue
            for (int i = 0; i < cpin.getResidues()[residx].numStates(); i++) {
               active_states_[residx][i] = cpin.getResidues()[residx].isProtonated(i);
//             // DEBUG
//             if (cpin.getResidues()[residx].isProtonated(i))
//                cout << "Residue " << residx << ", state " << i << " is true." << endl;
            }

        }else if (DEPROT.size() >= parts[1].size() &&
                  parts[1] == DEPROT.substr(0, parts[1].size())) {
            for (int i = 0; i < cpin.getResidues()[residx].numStates(); i++) {
               active_states_[residx][i] = !cpin.getResidues()[residx].isProtonated(i);
//             // DEBUG
//             if (!cpin.getResidues()[residx].isProtonated(i))
//                cout << "Residue " << residx << ", state " << i << " is true." << endl;
            }

        }else {
            cerr << "Error: Could not identify conditional criteria [" <<
                    parts[1] << "]" << endl;
            return ERR;
         }
      }
   }

   // Now filter out all of the residues that were not specified. If all states
   // in a given residue are 'false', then this residue was not specified in the
   // conditional probability, so turn them all to 'true'

   for (int i = 0; i < cpin.getTrescnt(); i++) {
      bool all_false = true;
      for (int j = 0; j < cpin.getResidues()[i].numStates(); j++) {
         if (active_states_[i][j]) {
            all_false = false;
            break;
         }
      }

      if (all_false)
         for (int j = 0; j < cpin.getResidues()[i].numStates(); j++)
            active_states_[i][j] = true;
   }

// DEBUG
// cout << setfill('=') << setw(80) << "=" << endl; setfill(' ');

   return OK;
}

bool ConditionalProb::SatisfiedBy(ProtVector invec) const {
   bool satisfied = true;
   for (size_t i = 0; i < invec.size(); i++)
      satisfied = satisfied && active_states_[i][invec[i]];
   return satisfied;
}
