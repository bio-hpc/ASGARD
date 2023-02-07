#ifndef PROTTRAJ_H
#define PROTTRAJ_H
#include <vector>
#include <ostream>

#include "conprob.h"
#include "cpin.h"
#include "cpout.h"
#include "types.h"

class ProtTraj {

   public:
      // Constructors
      ProtTraj(Cpin*, float, Record const&);

      // Destructor
//    ~ProtTraj();

      // Adds a point to the trajectory
      void AddPoint(Record const&);

      // Loads a full cpout
      void LoadCpout(CpoutFile);

      // Prints the calcpka-style output
      void PrintCalcpka(std::ostream &);
      void PrintCalcpka(std::ostream &, const int);

      // Prints the desired values for each residue for chunks of the simulation
      // of given size (in time steps)
      void PrintChunks(const int, std::string const&, const bool, const bool);

      // Prints the cumulative averages for each residue during the course of
      // the simulation
      void PrintCumulative(std::string const&, const int,
                           const bool, const bool);

      // Prints the running averages for each residue during the course of the
      // simulation
      void PrintRunningAvg(const int, const int, std::string const&,
                           const bool, const bool);

      // Prints a summary of the population of each residue in each state
      void PrintProtPop(std::string const&);

      // Prints the conditional probability over the course of a trajectory
      void PrintCondProb(std::string const&, std::vector<ConditionalProb> const&);

      // Prints the time series of conditional probabilities
      void PrintCondTimeseries(std::string const&, const int,
                               std::vector<ConditionalProb> const&);

   private:

      // Cpin file
      Cpin *cpin_;

      // Store the last point loaded
      ProtVector last_point_;

      // Store every snapshot that was saved
      std::vector<ProtVector> statelist_;

      // Store the number of residues we have
      int nres_;

      // Solvent pH
      float pH_;

      // Number of frames
      long long int nframes_;

      // Monte carlo time step
      int time_step_;

      typedef struct {
         std::vector<long long int> state_cnt;
         std::vector<int> prot_cnt;
         int nstates;
      } StateCount;

};
#endif /* PROTTRAJ_H */
