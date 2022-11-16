// main.cpp: Holds the main program

// Standard C++ headers
#include <iostream>
#include <fstream>
#include <iomanip>

// My headers
#include "cpin.h"
#include "cpout.h"
#include "cloptions.h"
#include "prottraj.h"
#include "test.h"
#include "types.h"
#include "utilities.h"

using namespace std;

int main(int argc, char**argv) {

   // Set up the command-line options and parse them
   CLOptions clopt = CLOptions(argc, argv);
   if (clopt.Parse())
      return 1;

   // Check the input
   if (clopt.CheckInput()) {
      cerr << "Error: Input errors detected! See messages above." << endl;
      return 1;
   }

// test_clopt(clopt); /* Uncomment to test command-line parsing */

   int nres = 0; // number of residues (for consistency tests)

   /* Set up the cpin and print some diagnostic information (only necessary if
    * we're not fixing REMD files
    */
   Cpin my_cpin;
   if (clopt.REMDPrefix().empty()) {
      if ( clopt.Cpin().empty() ) {
         cerr << "Error: No cpin file provided!" << endl;
         return 1;
      }
      if ( my_cpin.Parse(clopt.Cpin()) )
         return 1;
   
      nres = (int) my_cpin.getTrescnt();

      if (nres <= 0) {
         cerr << "Error: Did not detect any residues in " << 
                 clopt.Cpin() << endl;
         return 1;
      }
      if (clopt.Debug()) {
         cout << "There are " << nres << " titratable residues defined in " <<
                 my_cpin.getFilename() << ":" << endl;
         cout << "They are:" << endl;
         for (Cpin::ResIterator it = my_cpin.begin(); it != my_cpin.end(); it++) {
            cout << "\t" << setw(3) << it->getResname() << left << " " << 
                    setw(3) << it->getResnum() << " (" << it->numStates() 
                    << " states) [ ";
            for (int j = 0; j < it->numStates(); j++) {
               if (it->isProtonated(j))
                  cout << "P ";
               else
                  cout << "D ";
            }
            cout << "]" << endl;
         }
      }
   } // if clopt.REMDPrefix().empty()

   // Set up the cpouts
   CpoutList cpouts;
   for (CLOptions::cpout_iterator it = clopt.begin(); it != clopt.end(); it++) {
      CpoutFile c = CpoutFile(*it);
      // Skip over invalid cpouts
      if (!c.Valid()) {
         cerr << "Error: Cpout file " << *it << " is invalid! Skipping." << endl;
         continue;
      }
      // For REMD fixing where a cpin is unnecessary, make sure all cpouts have
      // the same number of residues, so set nres to the first cpout's Nres
      if (nres <= 0) nres = c.Nres();
      // Skip over cpouts with a residue mismatch
      if (c.Nres() != nres) {
         cerr << "Error: Cpout file " << *it << " has " << c.Nres() <<
                 " residues. I expected " << my_cpin.getTrescnt() <<
                 ". Skipping." << endl;
         continue;
      }
      if (clopt.Debug())
         cout << "Added [[ " << *it << " ]] to cpout list." << endl;
      cpouts.push_back(c);
   }

   if (clopt.Debug()) 
      cout << "Analyzing " << clopt.Cpouts().size() << " cpouts." << endl;
   if (cpouts.size() != clopt.Cpouts().size()) {
     cerr << "Error: Number of Cpout files " << cpouts.size() << 
             " does not equal number specified: " << clopt.Cpouts().size() << endl;
     return 1;
   }

   // Special-case REMD re-ordering
   if (!clopt.REMDPrefix().empty())
      return sort_remd_files(cpouts, clopt.REMDPrefix(), clopt.Overwrite());

   ProtTraj stats = ProtTraj(&my_cpin, cpouts[0].pH(), cpouts[0].GetRecord());
   for (cpout_iterator it = cpouts.begin(); it != cpouts.end(); it++) {
      // If we are here, then warn if we are about to use a REMD file
      if (!clopt.Expert())
         it->WarnRemd();
      stats.LoadCpout(*it);
   }

   // Do the normal calcpka-style output if requested
   if (clopt.Calcpka()) {
      if (clopt.Output().empty())
         stats.PrintCalcpka(cout);
      else {
         ofstream fd; fd.open(clopt.Output().c_str());
         stats.PrintCalcpka(fd);
         fd.close();
      }
   }

   // Do chunk analysis
   if (clopt.ChunkWindow() > 0)
      stats.PrintChunks(clopt.ChunkWindow(), clopt.ChunkOutput(),
                        clopt.PrintProtonated() && !clopt.pKa(), clopt.pKa());
   
   // Do cumulative analysis
   if (clopt.doCumulative())
      stats.PrintCumulative(clopt.CumulativeOutput(), clopt.Interval(),
                        clopt.PrintProtonated() && !clopt.pKa(), clopt.pKa());

   // Do running averages
   if (clopt.RunningAvgWindow() > 0)
      stats.PrintRunningAvg(clopt.RunningAvgWindow(), clopt.Interval(), 
                        clopt.RunningAvgOutput(),
                        clopt.PrintProtonated() && !clopt.pKa(), clopt.pKa());

   // Do protonation fraction dump
   if (clopt.Population().size() > 0)
      stats.PrintProtPop(clopt.Population());

   // Do conditional probabilities
   if (clopt.CondProbs().size() > 0) {
      for (CLOptions::prob_iterator it = clopt.condbegin(); 
                                    it != clopt.condend(); it++)
         if (it->Set(my_cpin) == ConditionalProb::ERR) {
            cerr << "Quitting due to errors above." << endl;
            return 1;
         }

      stats.PrintCondProb(clopt.ConditionalOutput(), clopt.CondProbs());
   }

   // Do conditional probability time series (we already set the conditional
   // probabilities above, so no need to do it again)
   if (clopt.CondProbs().size() > 0 && !clopt.ConditionalChunkOut().empty())
      stats.PrintCondTimeseries(clopt.ConditionalChunkOut(), clopt.Interval(),
                                clopt.CondProbs());

   if (clopt.Debug()) cout << "All done!" << endl;

   return 0;
}
