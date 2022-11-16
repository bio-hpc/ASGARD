/// test.cpp: Contains functions for testing various parts of the program

#include <cstdio>
#include "constants.h"
#include "test.h"

void test_clopt(CLOptions clopt) {
   
   printf("The verbose level is: %d\n", clopt.Verbose());
   printf("The cpouts are:\n");
   for (CLOptions::cpout_iterator it = clopt.begin(); it != clopt.end(); it++)
      printf("\t%s\n", it->c_str());
   printf("The cpin is: %s\n", clopt.Cpin().c_str());
   if (clopt.doCumulative())
      printf("Cumulative file name: %s\n", clopt.CumulativeOutput().c_str());
   else
      printf("NOT doing cumulative protonation calcs.\n");

   if (clopt.RunningAvgWindow() > 0)
      printf("Running average: Window %d; File %s\n", clopt.RunningAvgWindow(),
                      clopt.RunningAvgOutput().c_str());
   else
      printf("NOT doing running averages.\n");

   if (clopt.ChunkWindow() > 0)
      printf("Chunk statistics: Window %d; File %s\n", clopt.ChunkWindow(),
                      clopt.ChunkOutput().c_str());
   else
      printf("NOT doing chunk analysis.\n");

   if (clopt.Calcpka())
      if (clopt.Output().empty())
         printf("Output file: stdout\n");
      else
         printf("Output file: %s\n", clopt.Output().c_str());
   else
      printf("NOT doing calcpka analysis.\n");

   if (clopt.CheckInput())
      printf("Errors in input detected!\n");

   return;
}

void test_cpouts(std::vector<CpoutFile> cpouts) {
   
   std::vector< std::vector<int> > statelist;

   fprintf(stdout, "I have %d cpouts.\n", (int)cpouts.size());

   for (int i = 0; i < cpouts[0].Nres(); i++) {
      std::vector<int> newstates;
      statelist.push_back( newstates );
   }
   int i = 0;
   for (std::vector<CpoutFile>::iterator it = cpouts.begin();
         it != cpouts.end(); it++) {
      
      fprintf(stdout, "Analyzing Cpout file %s.\n", it->Filename().c_str());
      Record myrec;
      do {
         try {
            myrec = it->GetRecord();
         } catch (CpoutFinished &e) {
            break;
         }
         for (size_t j = 0; j < myrec.points.size(); j++)
            statelist[myrec.points[j].residue].push_back(myrec.points[j].state);
      } while (!it->Done());
      i++;
   }

   for (int j = 0; j < cpouts[0].Nres(); j++) {
      for (std::vector<int>::const_iterator it = statelist[j].begin();
               it != statelist[j].end(); it++) {
         fprintf(stdout, "Residue %4d State: %2d\n", j, *it);
      }
   }
   return;
}
