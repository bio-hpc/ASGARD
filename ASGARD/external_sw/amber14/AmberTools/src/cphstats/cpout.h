#ifndef CPOUT_H
#define CPOUT_H

#include <cstdio>
#include <vector>
#include <string>

#include <iostream>

// To support compression
#ifdef HASGZ
#  include "zlib.h"
#endif

#include "constants.h"

typedef struct {
   int residue;
   int state;
} RecordPoint;

typedef struct {
   std::vector<RecordPoint> points;
   int time_step;
   float time;
   float pH;
   bool full;
} Record;

class CpoutFile {
   public:
      enum FileType {ASCII=0, BZIP, GZIP};

      // Constructors
      CpoutFile(std::string const&);
      CpoutFile(const char*);
      
      // Getters
      bool Valid() const     { return valid_;           }
      bool Done() const      { return !valid_ || done_; }
      int Nres() const       { return nres_;            }
      int StepSize() const   { return step_size_;       }
      float pH() const       { return orig_ph_;         }
      float startTime() const{ return start_time_;      }
      void WarnRemd() const {
            if (remd_file_)
               std::cerr << "Warning: " << filename_ << " comes from a pH-REMD simulation! Not valid "
                    << "for pKa calculations." << std::endl << std::endl;
      }

      std::string Filename() { return filename_; }

      // Destructor
//    ~CpoutFile();

      // Get the next record
      Record GetRecord();

   private:
      // Auto-dispatch
      int Gets(char* c, int i) { if (type_ == ASCII) return AsciiGets(c, i);
#                                ifdef HASGZ
                                 if (type_ == GZIP) return GzGets(c, i);
#                                endif
                                 return 1;}
      void Rewind() { if (type_ == ASCII) fseek(fp_, 0, SEEK_SET);
#                     ifdef HASGZ
                      if (type_ == GZIP) gzseek(gzfp_, 0, SEEK_SET);
#                     endif
                    }

      void Close() { if (type_ == ASCII) fclose(fp_);
#                    ifdef HASGZ
                     if (type_ == GZIP) gzclose(gzfp_);
#                    endif
                   }
      // Real methods
#     ifdef HASGZ
      int GzGets(char*, int);
#     endif
      int AsciiGets(char*, int);

      // File objects
      FILE *fp_;
#     ifdef HASGZ
      gzFile gzfp_;
#     endif

      // File type (ASCII? Gzip?)
      FileType type_;

      std::string filename_; // Original file name
      bool valid_;       // Is this a valid file?
      bool done_;        // Is this file done reading?
      float orig_ph_;    // pH on the opening full record
      int step_size_;    // Monte carlo step size
      int nres_;         // Number of residues defined in this cpout
      float start_time_; // The time in the first full record
      bool remd_file_;   // Is this file from a REMD run?
};

#endif /* CPOUT_H */
