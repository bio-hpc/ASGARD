// cloptions.cpp: Manage the command-line options and parsing

// Standard C++ headers
#include <fstream>
#include <iostream>

// My headers
#include "cloptions.h"
#include "string_manip.h"
#include "utilities.h"

using namespace std;

CLOptions::CLOptions(int argc, char**argv) :
verbose_(0),
do_calcpka_(true),
cumulative_(false),
runavgwin_(0),
chunksize_(0),
overwrite_(false),
interval_(1000),
protonated_(true),
time_step_(0.002f),
pKa_(false),
debug_(false),
expert_(false)
{
   
   // Initialize some strings
   runavgout_ = string("running_avgs.dat");
   cumout_ = string("cumulative.dat");
   chunkout_ = string("chunk.dat");
   condprobf_ = string("conditional_prob.dat");

   // Now parse everything
   int i = 1;
   vector<bool> marked(argc, false);

   prog_ = string(argv[0]);
   prog_ = prog_.substr(prog_.find_last_of('/')+1);
   parse_return_ = OK;
   marked[0] = true; // this is the program name

   while (i < argc) {
      string arg = string(argv[i]);
      bool has_another = (i != argc - 1);
      string argnext;
      if (has_another)
         argnext = string(argv[i+1]);

      if (arg == "-h" || arg == "--help") {
         parse_return_ = HELP;
         break;
     }else if (arg == "-O" || arg == "--overwrite") {
         marked[i] = true;
         overwrite_ = true;
     }else if (arg == "-V" || arg == "--version") {
         parse_return_ = VERSION;
         break;
     }else if (arg == "-i" || arg == "--cpin") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         cpin_ = argnext;
     }else if (arg == "-o" || arg == "--calcpka-output") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         output_ = argnext;
     }else if (arg == "-v" || arg == "--verbose") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         verbose_ = StringToInt(argnext);
     }else if (arg == "-R" || arg == "--running-avg-out") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         runavgout_ = argnext;
     }else if (arg == "--chunk-out") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         chunkout_ = argnext;
     }else if (arg == "--cumulative-out") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         cumout_ = argnext;
     }else if (arg == "--calcpka") {
         marked[i] = true;
         do_calcpka_ = true;
     }else if (arg == "--no-calcpka") {
         marked[i] = true;
         do_calcpka_ = false;
     }else if (arg == "-r" || arg == "--running-avg") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         runavgwin_ = StringToInt(argnext);
     }else if (arg == "-t" || arg == "--time-step") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         time_step_ = StringToFloat(argnext);
     }else if (arg == "--chunk") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         chunksize_ = StringToInt(argnext);
     }else if (arg == "--cumulative") {
         marked[i] = true;
         cumulative_ = true;
     }else if (arg == "--fix-remd") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         reorder_prefix_ = argnext;
     }else if (arg == "-n" || arg == "--interval") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         interval_ = StringToInt(argnext);
     }else if (arg == "-p" || arg == "--protonated") {
         marked[i] = true;
         protonated_ = true;
         pKa_ = false;
     }else if (arg == "-d" || arg == "--deprotonated") {
         marked[i] = true;
         protonated_ = false;
         pKa_ = false;
     }else if (arg == "-a" || arg == "--pKa") {
         marked[i] = true;
         pKa_ = true;
     }else if (arg == "--population") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         population_ = argnext;
     }else if (arg == "-c" || arg == "--conditional") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         condprobs_.push_back( ConditionalProb(argnext) );
     }else if (arg == "--conditional-output") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         condprobf_ = argnext;
     }else if (arg == "--chunk-conditional") {
         if (!has_another) {parse_return_ = INSUFARGS; break;}
         marked[i++] = true;
         marked[i] = true;
         condprob_chunk_ = argnext;
     }else if (arg == "--debug") {
         marked[i] = true;
         debug_ = true;
     }else if (arg == "--expert") {
         marked[i] = true;
         expert_ = true;
     }else if (arg == "--novice") {
         marked[i] = true;
         expert_ = false;
     }else if (arg[0] == '-') {
         cerr << "Unrecognized command-line option: " << arg << endl;
         parse_return_ = ERR;
         break;
      }
      i++;
   }

   if (parse_return_ == OK) {
      // Now fill the list of cpout files
      for (int j = 1; j < argc; j++) {
         if (!marked[j]) {
            marked[j] = true;
            cpouts_.push_back(string(argv[j]));
         }
      }
   }

   if (cpouts_.size() == 0 && parse_return_ == OK) parse_return_ = HELP;
}

void CLOptions::Help() {

   Usage();
   cout << endl;
   cout << "General Options:" << endl;
   cout << "    -h, --help     Print this help and exit." << endl;
   cout << "    -V, --version  Print the version number and exit." << endl;
   cout << "    -O, --overwrite" << endl;
   cout << "                   Allow existing outputs to be overwritten." << endl;
   cout << "    --debug        Print out information about the files that are" << endl;
   cout << "                   being read in and used for the calculations." << endl;
   cout << "    --expert       I will consider you an expert user and NOT warn" << endl;
   cout << "                   you if you try to compute statistics from REMD-based" << endl;
   cout << "                   files before using --fix-remd [NOT default behavior]" << endl;
   cout << "    --novice       I will warn you if you try to use REMD-based files" << endl;
   cout << "                   to compute statistics. [Default behavior]" << endl;
   cout << endl;
   cout << "Input Files and Options:" << endl;
   cout << "    -i FILE, --cpin FILE" << endl;
   cout << "                   Input cpin file (from sander) with titrating residue" << endl;
   cout << "                   information." << endl;
   cout << "    -t FLOAT, --time-step FLOAT" << endl;
   cout << "                   This is the time step in ps you used in your simulations." << endl;
   cout << "                   It will be used to print data as a function of time." << endl;
   cout << "                   Default is 2 fs (0.002)" << endl;
   cout << endl;
   cout << "Output Files:" << endl;
   cout << "    -o FILE, --calcpka-output FILE" << endl;
   cout << "                   File to which the standard `calcpka'-type statistics" << endl;
   cout << "                   are written. Default is stdout" << endl;
   cout << "    -R FILE, --running-avg-out FILE" << endl;
   cout << "                   Output file where the running averages of time series" << endl;
   cout << "                   data for each residue is printed (see [Output Options]" << endl;
   cout << "                   below for details). Default is [running_avgs.dat]" << endl;
   cout << "    --chunk-out FILE" << endl;
   cout << "                   Output file where the time series data calculated" << endl;
   cout << "                   over chunks of the simulation are printed (see " << endl;
   cout << "                   [Output Options] below for details)." << endl;
   cout << "                   Default is [chunk.dat]" << endl;
   cout << "    --cumulative-out FILE" << endl;
   cout << "                   Output file where the cumulative time series data" << endl;
   cout << "                   is printed (see [Output Options] below for details)." << endl;
   cout << "                   Default is [cumulative.dat]" << endl;
   cout << "    --population FILE" << endl;
   cout << "                   Output file where protonation state populations are" << endl;
   cout << "                   printed for every state of every residue." << endl;
   cout << "    --conditional-output FILE" << endl;
   cout << "                   Output file with requested conditional probabilities." << endl;
   cout << "                   Default is [conditional_prob.dat]." << endl;
   cout << "    --chunk-conditional FILE" << endl;
   cout << "                   Prints a time series of the conditional probabilities over" << endl;
   cout << "                   a trajectory split up into chunks." << endl;
   cout << endl;
   cout << "Output Options:" << endl;
   cout << "  These options modify how the output files will appear" << endl;
   cout << endl;
   cout << "    -v INT, --verbose INT" << endl;
   cout << "                   Controls how much information is printed to the" << endl;
   cout << "                   calcpka-style output file. Options are:" << endl;
   cout << "                      (0) Just print fraction protonated. [Default]" << endl;
   cout << "                      (1) Print everything calcpka prints." << endl;
   cout << "    -n INT, --interval INT" << endl;
   cout << "                   An interval between which to print out time series data" << endl;
   cout << "                   like `chunks', `cumulative' data, and running averages." << endl;
   cout << "                   It is also used as the 'window' of the conditional" << endl;
   cout << "                   probability time series (--chunk-conditional)." << endl;
   cout << "                   Default [1000]" << endl;
   cout << "    -p, --protonated" << endl;
   cout << "                   Print out protonation fraction instead of deprotonation" << endl;
   cout << "                   fraction in time series data (Default behavior)." << endl;
   cout << "    -d, --deprotonated" << endl;
   cout << "                   Print out deprotonation fraction instead of protonation" << endl;
   cout << "                   fraction in time series data." << endl;
   cout << "    -a, --pKa      Print predicted pKas (via Henderson-Hasselbalch) in place" << endl;
   cout << "                   of fraction (de)protonated. NOT default behavior." << endl;
   cout << endl;
   cout << "Analysis Options:" << endl;
   cout << "  These options control which analyses are done. By default, only" << endl;
   cout << "  the original, calcpka-style analysis is done." << endl;
   cout << endl;
   cout << "    --calcpka      Triggers the calcpka-style output [On by default]" << endl;
   cout << "    --no-calcpka   Turns off the calcpka-style output" << endl;
   cout << "    -r WINDOW, --running-avg WINDOW" << endl;
   cout << "                   Defines a window size for a moving, running average" << endl;
   cout << "                   time series. <WINDOW> is the number of MD steps (NOT" << endl;
   cout << "                   the number of MC exchange attempts)." << endl;
   cout << "    --chunk WINDOW" << endl;
   cout << "                   Computes the time series data over a chunk of the" << endl;
   cout << "                   simulation of size <WINDOW> time steps. See above for" << endl;
   cout << "                   details." << endl;
   cout << "    --cumulative   Computes the cumulative average time series data (see above" << endl;
   cout << "                   for options) over the course of the trajectory." << endl;
   cout << "    --fix-remd PREFIX" << endl;
   cout << "                   This option will trigger " << prog_ << " to reassemble the " << endl;
   cout << "                   titration data into pH-specific ensembles. This" << endl;
   cout << "                   is an exclusive mode of the program---no other" << endl;
   cout << "                   analyses will be done." << endl;
   cout << "    -c CONDITIONAL, --conditional CONDITIONAL" << endl;
   cout << "                   Evaluates conditional probabilities. CONDITIONAL should be a" << endl;
   cout << "                   string of the format:" << endl;
   cout << "                         <resid>:<state>,<resid>:<state>,..." << endl;
   cout << "                     or" << endl;
   cout << "                         <resid>:PROT,<resid>:DEPROT,..." << endl;
   cout << "                     or" << endl;
   cout << "                         <resid>:<state1>;<state2>,<resid>:PROT,..." << endl;
   cout << "                   Where <resid> is the residue number in the prmtop (NOT the" << endl;
   cout << "                   cpin) and <state> is either the state number or (p)rotonated" << endl;
   cout << "                   or (d)eprotonated, case-insensitive" << endl;
   cout << endl;
   cout << "This program analyzes constant pH output files (cpout) from Amber." << endl;
   cout << "These output files can be compressed using gzip compression. The" << endl;
   cout << "compression will be detected automatically by the file name extension." << endl;
   cout << "You must have the gzip headers for this functionality to work." << endl;

   return;
}

void CLOptions::Usage() {
   cout << "Usage: " << prog_ << " [-O] [-V] [-h] [-i <cpin>] [-t] [-o FILE] [-R FILE -r INT]" << endl;
   cout << "             [--chunk INT --chunk-out FILE] [--cumulative --cumulative-out FILE]" << endl;
   cout << "             [-v INT] [-n INT] [-p|-d] [--calcpka|--no-calcpka] [--fix-remd]" << endl;
   cout << "             [--population FILE] [-c CONDITION -c CONDITION -c ...]" << endl;
   cout << "             [--conditional-output FILE] [--chunk-conditional FILE]" << endl;
   cout << "             cpout1 [cpout2 [cpout3 ...] ]" << endl;

   return;
}

void CLOptions::Version() {
   cout << prog_ << ": Version " << VERSION_STR << endl;
}

int CLOptions::Parse() {
   if (parse_return_ == HELP)
      Help();
   else if (parse_return_ == VERSION)
      Version();
   else if (parse_return_ == INSUFARGS) {
      cerr << "Insufficient arguments!" << endl;
      Usage();
  }else if (parse_return_ == ERR) {
      Usage();
   }

   return parse_return_;
}

int CLOptions::CheckInput() {
   int inerr = 0;
   if (runavgwin_ < 0) {
      cerr << "Error: -r/--running-avg window must be non-negative!" << endl;
      inerr = 1;
   }

   if (chunksize_ < 0) {
      cerr << "Error: --chunk window must be non-negative!" << endl;
      inerr = 1;
   }

   if (time_step_ <= 0) {
      cerr << "Error: --time-step must be a positive number!" << endl;
      inerr = 1;
   }

   if (interval_ <= 0) {
      cerr << "Error: -n/--interval must be a positive integer!" << endl;
      inerr = 1;
   }

   // Make sure no files exist that we want to overwrite
   if (runavgwin_ > 0 && !overwrite_)
      if (fexists(runavgout_)) {
         cerr << "Error: " << runavgout_ << " exists; not overwriting." << endl;
         inerr = 1;
      }

   if (!condprob_chunk_.empty() && !overwrite_)
      if (fexists(condprob_chunk_)) {
         cerr << "Error: " << condprob_chunk_ << " exists; not overwriting." << endl;
         inerr = 1;
      }

   if (do_calcpka_ && !overwrite_)
      if (!output_.empty() && fexists(output_)) {
         cerr << "Error: " << output_ << " exists; not overwriting." << endl;
         inerr = 1;
      }

   if (chunksize_ > 0 && !overwrite_)
      if (fexists(chunkout_)) {
         cerr << "Error: " << chunkout_ << " exists; not overwriting." << endl;
         inerr = 1;
      }

   if (cumulative_ && !overwrite_)
      if (fexists(cumout_)) {
         cerr << "Error: " << cumout_ << " exists; not overwriting." << endl;
         inerr = 1;
      }

   if (!population_.empty() && !overwrite_)
      if (fexists(population_)) {
         cerr << "Error: " << population_ << " exists; not overwriting." << endl;
         inerr = 1;
      }

   if (!condprobf_.empty() && !overwrite_)
      if (fexists(condprobf_)) {
         cerr << "Error: " << condprobf_ << " exists; not overwriting." << endl;
         inerr = 1;
      }

   // Make sure we provided a prefix and didn't accidentally just dump a bunch
   // of cpout files on the command line without a prefix
   if (!reorder_prefix_.empty() && fexists(reorder_prefix_)) {
      cerr << "--------| Interpreting " << reorder_prefix_ << " as the new REMD file prefix. " 
           << "If this is" << endl;
      cerr << "Warning | a cpout file, kill this program and re-run with a correct" << endl;
      cerr << "--------| REMD file prefix." << endl;
   }

   return inerr;
}
