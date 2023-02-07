// prottraj.cpp: Set up vectors of all of the protonation states for each snapshot

#include <cmath>   // log10
#include <iomanip> // setw
#include <sstream> // istringstream
#include <fstream> // ofstream
#include <ios>     // fixed

#include "constants.h"
#include "exceptions.h"
#include "prottraj.h"

// Some simple binary sorting macros
#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a > b ? b : a)

using namespace std;

ProtTraj::ProtTraj(Cpin* cpin, float pH, Record const& recin) :
cpin_(NULL),
nres_(0),
pH_(0.0f),
nframes_(0),
time_step_(0)
{
   cpin_ = cpin;
   nres_ = cpin->getTrescnt();
   pH_ = pH;

   // Set the first point
   for (int i = 0; i < nres_; i++)
      last_point_.push_back(0);

   for (vector<RecordPoint>::const_iterator it = recin.points.begin();
            it != recin.points.end(); it++)
      last_point_[it->residue] = it->state;

   statelist_.push_back(last_point_);
}

/// Loads a full cpout file point by point
void ProtTraj::LoadCpout(CpoutFile cpout) {
   if (time_step_ == 0) time_step_ = cpout.StepSize();
   while (!cpout.Done()) {
      Record r;
      try {
         r = cpout.GetRecord();
      } catch ( CpoutFinished &e ) {
         break;
      }
      // Do not add an empty record
      if (r.points.empty()) continue;
      AddPoint(r);
   }
}

/// Loads a single frame into the trajectory
void ProtTraj::AddPoint(Record const& recin) {
   for (vector<RecordPoint>::const_iterator it = recin.points.begin();
            it != recin.points.end(); it++)
      last_point_[it->residue] = it->state;
   statelist_.push_back(last_point_);
   nframes_++;
}

/// Prints the calcpka-style output
void ProtTraj::PrintCalcpka(ostream& fd, const int start) {

   // First calculate the statistics
   vector<long long int> nprot(nres_, 0ll);
   vector<int> transitions(nres_, 0);
   long long int totprot = 0ll;

   // Start with the first point, but skip it in the iteration
   last_point_ = statelist_[start];
   for (int i = start + 1; i < nframes_ + 1; i++) {
      int j = 0;
      for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
         transitions[j] += (int) (rit->isProtonated(last_point_[j]) !=
                                  rit->isProtonated(statelist_[i][j]));
         long long int protadd = (long long int) (rit->isProtonated(statelist_[i][j]));
         nprot[j] += protadd;
         totprot += (long long int) rit->numProtons(statelist_[i][j]);
         j++;
      }
      // Reassign the last point to the current one before incrementing
      last_point_ = statelist_[i];
   }

   // Now do the printing
   fd.precision(3);
   fd << fixed;
   fd << "Solvent pH is " << setw(8) << pH_ << endl;
   int i = 0;
   for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
      double pKa = pH_ - log10( (double) (nframes_-nprot[i]) / (double) nprot[i] );
      float offset = (float) pKa - pH_;
      fd << setw(3) << rit->getResname() << " " << setw(4) << left << rit->getResnum()
         << ": Offset " << right << setw(6) << offset << "  Pred " << setw(6) << pKa
         << "  Frac Prot " << setw(5) << (double) nprot[i] / (double) nframes_
         << " Transitions " << setw(9) << transitions[i] << endl;
      i++;
   }
   fd << endl << "Average total molecular protonation: " << setw(7)
      << (double)totprot / (double)nframes_ << endl;

   return;
}

/// Prints the calcpka-style output for the entire simulation
void ProtTraj::PrintCalcpka(ostream &fd) {
   PrintCalcpka(fd, 0);
}

/// Print a time series of the desired properties in 'chunks' of the simulation
void ProtTraj::PrintChunks(const int window, string const& fname,
                           const bool print_prot, const bool print_pka) {

   // Open up the file and write a header (set fixed precision)
   ofstream fp(fname.c_str());
   if (!fp.is_open())
      throw FileIOError("Could not open " + fname + " for writing!");
   fp << fixed << "#Time step ";
   for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
      stringstream iss;
      iss << rit->getResname() << " " << rit->getResnum();
      fp << setw(8) << iss.str() << " ";
   }
   fp << " Total Avg. Prot." << endl;

   int interval = window / time_step_;
   int first = 0;
   while (first < nframes_+1) {
      vector<long long int> nprot(nres_, 0ll);
      long long int totprot = 0ll;
      int last = MIN(nframes_+1, first + interval);
      for (int i = first; i < last; i++) {
         int j = 0;
         for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
            long long int protadd = (long long int) (rit->isProtonated(statelist_[i][j]));
            nprot[j] += protadd;
            totprot += (long long int) rit->numProtons(statelist_[i][j]);
            j++;
         }
      }
      // Now write the information
      int j = 0;
      fp << setw(10) << time_step_*(first + interval / 2) << " ";
      for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
         double fracprot = (double)nprot[j] / (double)(last - first);
         if (print_prot) {
            // We want the fraction protonated
            fp << setprecision(5) << setw(8) << fracprot << " ";
        }else if (print_pka) {
            // We want the pKa
            double pKa = pH_ - log10( (1.0 - fracprot) / fracprot );
            fp << setprecision(4) << setw(8) << pKa << " ";
        }else {
            // We want the fraction deprotonated
            fp << setprecision(5) << setw(8) << 1.0-fracprot << " ";
         }
         j++;
      }
      double fracprot = (double)totprot/(double)(last - first);
      fp << setprecision(6) << setw(17) << fracprot << endl;
      first += interval;
   }

   // Now that I'm here (and done), close the file
   fp.close();
   return;
}

/// Prints a cumulative running average time series of the desired property
void ProtTraj::PrintCumulative(string const& fname, const int interval,
                               const bool print_prot, const bool print_pka) {
   
   // Open up the file and write a header
   ofstream fp(fname.c_str());
   if (!fp.is_open())
      throw FileIOError("Could not open " + fname + " for writing!");
   fp << fixed << "#Time step ";
   for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
      stringstream iss;
      iss << rit->getResname() << " " << rit->getResnum();
      fp << setw(8) << iss.str() << " ";
   }
   fp << " Total Avg. Prot." << endl;

   // Now go through the trajectory
   vector<long long int> nprot(nres_, 0ll);
   long long int totprot = 0ll;
   int c = 0; // simple counter
   for (int i = 1; i < nframes_ + 1; i++) {
      
      { // scope this
      int j = 0;
      for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
         long long int protadd = (long long int) (rit->isProtonated(statelist_[i][j]));
         nprot[j] += protadd;
         totprot += (long long int) rit->numProtons(statelist_[i][j]);
         j++;
      }
      }

      if (c * time_step_ >= interval) {
         c = 0;
         // Now we print a point
         int j = 0;
         fp << setw(10) << (i-1)*time_step_ << " ";
         for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
            double fracprot = (double)nprot[j] / (double)(i+1);
            if (print_prot) {
               // We want the fraction protonated
               fp << setprecision(5) << setw(8) << fracprot << " ";
           }else if (print_pka) {
               // We want the pKa
               double pKa = pH_ - log10( (1.0 - fracprot) / fracprot );
               fp << setprecision(4) << setw(8) << pKa << " ";
           }else {
               // We want the fraction deprotonated
               fp << setprecision(5) << setw(8) << 1.0-fracprot << " ";
            }
            j++;
         }
         double fracprot = (double)totprot / (double)i;
         fp << setprecision(6) << setw(17) << fracprot << endl;
      }
      c++; // heh
   }

   fp.close();
   return;
}

/// Print the rolling/running average time series with a given window
void ProtTraj::PrintRunningAvg(const int window, const int interval,
                               string const& fname, const bool print_prot,
                               const bool print_pka) {
   const int halfwin = window / (2 * time_step_);

   // Open the file for writing and print a header
   ofstream fp(fname.c_str());
   if (!fp.is_open())
      throw FileIOError("Could not open " + fname + " for writing!");
   fp << fixed << "#Time step ";
   for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
      stringstream iss;
      iss << rit->getResname() << " " << rit->getResnum();
      fp << setw(8) << iss.str() << " ";
   }
   fp << " Total Avg. Prot." << endl;

   for (int i = 1; i < nframes_ + 1; i+= interval/time_step_) {
      vector<long long int> nprot(nres_, 0ll);
      long long int totprot = 0ll;
      // Loop over all frames we should include here
      for (int j = MAX(0, i-halfwin); j < MIN(nframes_, i+halfwin); j++) {
         // Loop over every residue
         int k = 0;
         for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
            long long int protadd = (long long int) (rit->isProtonated(statelist_[j][k]));
            nprot[k] += protadd;
            totprot += (long long int) rit->numProtons(statelist_[j][k]);
            k++;
         }
      }
      // Print this frame to the output file
      int j = 0;
      int n = MIN(nframes_, i+halfwin) - MAX(0, i-halfwin);
      fp << setw(10) << time_step_*i << " ";
      for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
         double fracprot = (double)nprot[j] / (double)n;
         if (print_prot) {
            // We want the fraction protonated
            fp << setprecision(5) << setw(8) << fracprot << " ";
        }else if (print_pka) {
            // We want the pKa
            double pKa = pH_ - log10( (1.0 - fracprot) / fracprot );
            fp << setprecision(4) << setw(8) << pKa << " ";
        }else {
            // We want the fraction deprotonated
            fp << setprecision(5) << setw(8) << 1.0-fracprot << " ";
         }
         j++;
     }
     double fracprot = (double)totprot / (double)n;
     fp << setprecision(6) << setw(17) << fracprot << endl;
   }

   fp.close();
   return;
}

/** Print the population fraction of the entire simulation for each state and
  * each residue
  */
void ProtTraj::PrintProtPop(string const& fname) {
   // Set up the table of populations, initializing everything to zero
   vector<StateCount> pop_counts;
   // Now loop through every state and every residue, building pop_counts
   for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++) {
      StateCount newstate;
      newstate.nstates = 0;
      for (int i = 0; i < rit->numStates(); i++) {
         newstate.state_cnt.push_back(0ll);
         newstate.prot_cnt.push_back(rit->numProtons(i));
         newstate.nstates++;
      }
      pop_counts.push_back(newstate);
   }

   // Now loop through the entire trajectory, building the population list
   for (int i = 1; i < nframes_ + 1; i++) {
      ProtVector curst = statelist_[i];
      for (int j = 0; j < nres_; j++)
         pop_counts[j].state_cnt[ statelist_[i][j] ] += 1ll;
   }

   // Now print out the population counts
   ofstream fp(fname.c_str());
   if (!fp.is_open())
      throw FileIOError("Could not open " + fname + " for writing!");

   // Get the max number of states
   int maxst = cpin_->getResidues()[0].numStates();
   for (Cpin::ResIterator rit = cpin_->begin(); rit != cpin_->end(); rit++)
      maxst = MAX(maxst, rit->numStates());

   // Write the header
   fp << right << fixed << setw(17) << "Residue Number" << " ";
   for (int i = 0; i < maxst; i++)
      fp << setw(9) << "State" << setw(3) << i << ' ';
   fp << endl << setfill('-') << setw(maxst*13+18) << '-' << endl;

   // Now loop through every residue and print out the populations
   for (int i = 0; i < nres_; i++) {
      stringstream iss;
      iss << "Residue: " << cpin_->getResidues()[i].getResname() << " "
          <<  cpin_->getResidues()[i].getResnum();
      fp << setfill(' ') << left << setw(17) << iss.str() << " ";
      for (int j = 0; j < pop_counts[i].nstates; j++) {
         double pop = (double) pop_counts[i].state_cnt[j] / (double) nframes_;
         fp << setprecision(6) << setw(8) << pop << " (" << pop_counts[i].prot_cnt[j]
            << ") ";
      }
      fp << endl;
   }

   fp.close();
}

/// Prints out the conditional probabilities
void ProtTraj::PrintCondProb(string const& fname,
                             vector<ConditionalProb> const& condprobs) {
   // Determine how large the field has to be to print out all of the
   // conditional probabilities
   size_t max_size = condprobs[0].str().size() + (size_t) 1;
   for (vector<ConditionalProb>::const_iterator it = condprobs.begin();
        it != condprobs.end(); it++) {
      max_size = MAX(max_size, it->str().size() + (size_t) 1);
   }
   max_size = MAX(max_size, 25);

   // Open the file, then loop through all conditional probabilities, calculate
   // the conditional probability, then print it
   ofstream fp(fname.c_str());
   if (!fp.is_open())
      throw FileIOError("Could not open " + fname + " for writing!");
   fp << fixed << setw(max_size) << "Conditional Probabilities" << " "
      << setw(10) << "Fraction" << endl;
   fp.precision(6);
   for (vector<ConditionalProb>::const_iterator it = condprobs.begin();
        it != condprobs.end(); it++) {
      // Now loop through the whole trajectory
      long long int counted = 0;
      for (int i = 1; i < nframes_ + 1; i++)
         counted += (long long int) it->SatisfiedBy(statelist_[i]);

      double frac = (double) counted / (double) nframes_;
      fp << setw(max_size) << it->str() << " " << setw(10) << frac << endl;
   }
   fp.close();
}

void ProtTraj::PrintCondTimeseries(string const& fname, const int window,
                                   vector<ConditionalProb> const& condprobs) {
   // Determine how large the field has to be to print out all of the
   // conditional probabilities
   size_t max_size = 10;
   for (vector<ConditionalProb>::const_iterator it = condprobs.begin();
        it != condprobs.end(); it++) {
      max_size = MAX(max_size, it->str().size() + (size_t) 1);
   }

   // Open the file, then loop through all conditional probabilities, calculate
   // the conditional probability, then print it
   ofstream fp(fname.c_str());
   if (!fp.is_open())
      throw FileIOError("Could not open " + fname + " for writing!");
   fp.precision(5);
   fp << fixed << " Time step";
   for (vector<ConditionalProb>::const_iterator it = condprobs.begin();
        it != condprobs.end(); it++)
      fp << " " << setw(max_size) << it->str();

   fp << endl;

   int interval = window / time_step_;
   int first = 0;
   while (first < nframes_+1) {
      vector<long long int> nprot(nres_, 0ll);
      fp << setw(10) << first*time_step_;
      int last = MIN(nframes_+1, first+interval);
      for (vector<ConditionalProb>::const_iterator it = condprobs.begin();
           it != condprobs.end(); it++) {
         long long int counted = 0ll;
         for (int i = first; i < last; i++)
            counted += (long long int) it->SatisfiedBy(statelist_[i]);
         double frac = (double) counted / (double) (last-first);
         fp << " " << setw(max_size) << frac;
      }
      fp << endl;
      first += interval;
   }

   fp.close();
}
