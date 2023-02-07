#ifndef CPIN_H
#define CPIN_H

#include <string>
#include <vector>
#include "constants.h"

// External Fortran functions for parsing the cpin

typedef struct {
   int num_states, first_atom,  num_atoms,  first_state, first_charge;
} StateInfo;

extern "C" {
   void parse_cpin_(int*, int[TITR_STATES_C], StateInfo[TITR_RES_C],
                        char[TITR_RES_C+1][40], char[FN_LEN], int*);
}

/// A titratable residue
class TitratableResidue {
   public:
      // Constructors (only the implemented/used ones are uncommented)
//    TitratableResidue(const char*);
//    TitratableResidue(std::string const&);
      TitratableResidue(std::string const&, std::vector<int> const&);
//    TitratableResidue(const char*, int*);

      // A way to add a list of states
      void defineStates(std::vector<int>);
      void defineStates(int*);

      // In-line Setters
      void setResname(const char* resname) { resname_ = std::string(resname); }
      void setResname(std::string const& resname) { resname_ = resname; }
      void setResnum(int resnum) { resid_ = resnum; }
      
      // In-line Getters
      std::string getResname() const { return resname_; }
      int getResnum() const { return resid_; }

      // Determine how many states are present
      int numStates() const { return protonated_.size(); }

      // Determine if a particular state is protonated or not
      bool isProtonated(int state) const { return protonated_[state]; }
      // Determine how many protons are in a specific state
      int numProtons(int state) const { return protcnts_[state]; }

   private:
      /// How many protons are in each state?
      std::vector<int> protcnts_;

      /// Is each state protonated?
      std::vector<bool> protonated_;

      /// Name of the titratable residue
      std::string resname_;

      /// Number of the residue
      int resid_;
};

/// Wrapper for calling the Fortran cpin
class Cpin {
   public:
      // Constructors
      Cpin();

      // Get the data
      std::vector<TitratableResidue> getResidues() const { return residues_; }

      // Provide an iterator over the data
      typedef std::vector<TitratableResidue>::const_iterator ResIterator;
      ResIterator begin() const { return residues_.begin(); }
      ResIterator end()   const { return residues_.end();   }

      // Parse the cpin
      int Parse(const char*);
      int Parse(std::string const&);

      // Get the number of residues
      int getTrescnt() const { return trescnt_; }

      // Get the name of the file
      std::string getFilename() const { return filename_; }

      bool isRead() const { return is_valid_; }

   private:
      /// Is our cpin file valid?
      bool is_valid_;

      /// CpHMD state info array (analogous to the one used in sander)
      StateInfo stateinf_[TITR_RES_C];

      /// Number of titratable residues defined in the cpin
      int trescnt_;

      /// Number of 'active' protons in each state
      int protcnt_[TITR_STATES_C];

      /// Name of each titratable residue
      std::vector<std::string> resname_;

      /// List of titratable residues
      std::vector<TitratableResidue> residues_;

      /// Name of the file
      std::string filename_;

};
#endif /* CPIN_H */
