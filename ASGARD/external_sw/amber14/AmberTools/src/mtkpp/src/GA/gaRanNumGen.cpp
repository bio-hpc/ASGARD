//=============================================================
//   Filename : gaRanNumGen.cpp
// Written by : Martin Peters
// Project    : GA
// ------------------------------------------------------------
// description: Class gaRanNumGen
//=============================================================

//#include "gaRanNumGen.h"

// ============================================================
// Function : gaRanNumGen()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
//gaRanNumGen::gaRanNumGen(int seed):itsSeed(seed) {}

// ============================================================
// Function : ~gaRanNumGen()
// ------------------------------------------------------------
// Destructor for the class.
// All data is destroyed.
// ============================================================
//gaRanNumGen::~gaRanNumGen() {}

// ============================================================
// Function : setSeed()
// ------------------------------------------------------------
// Set Random Number Generation Seed
// ============================================================
/*void gaRanNumGen::setSeed(int seed)
{
    
    if (seed < 0) {
      srand(time(NULL));
    }
    else {
      srand(seed);
    }
}*/

// ============================================================
// Function : ranNumBetweenZeroAndOne()
// ------------------------------------------------------------
// Returns a random number (double) between 0 and 1.
// ============================================================
//void gaRanNumGen::ranNumBetweenZeroAndOne()
//{
//    return ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
//}

// ============================================================
// Function : ranNumBetweenZeroAndX()
// ------------------------------------------------------------
// Returns a random number (double) between 0 and X.
// ============================================================
//void gaRanNumGen::ranNumBetweenZeroAndX(int X)
//{
//    return ranNumBetweenZeroAndOne() * X;
//}

/* generate some random numbers */
//  printf ("A number between 0 and 100: %d\n", rand()%100);
//  printf ("A number between 20 and 30: %d\n", rand()%10+20);
// r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
/*
main()
{
// DEFINE VARIABLES 
   	int seed; 
      // Choose any integer value for seed to initialize the random number generator with. 
      //            The same value of seed will always produce the same sequence of random numbers.
      //           To get a different sequence of random numbers next time, choose a different value of seed.

Note: on Microsoft Visual C++ 6.0, 
if seed is negative than a seed value from the system clock is used. 
Also check for a randomize() function on your compiler. 

double r;		// random value in range [0,1) 

long int M;		// user supplied upper boundary 

double x;		// random value in range [0,M) 
int y;		// random integer in range [0,M) if M is an integer then range = [0,M-1] 
int z;		// random integer in range [1,M+1) if M is an integer then range = [1,M] 
int count;		// just a variable we need to count how many random numbers we've made for this example 

seed = 10000;		// choose a seed value 
srand(seed);		// initialize random number generator

M = 1000;// Choose M. Upper bound 
    for (count=1; count<=20; ++count) {
      r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
      // r is a random floating point value in the range [0,1) {including 0, not including 1}. 
      // Note we must convert rand() and/or RAND_MAX+1 to floating point values to avoid integer division. 
      // In addition, Sean Scanlon pointed out the possibility that RAND_MAX may be the largest positive 
      // integer the architecture can represent, so (RAND_MAX+1) may result in an overflow, or more likely 
      // the value will end up being the largest negative integer the architecture can represent, so to avoid 
      // this we convert RAND_MAX and 1 to doubles before adding. 
      x = (r * M);
      // x is a random floating point value in the range [0,M) {including 0, not including M}. 
      y = (int) x;
      // y is a random integer in the range [0,M) {including 0, not including M}. 
      //    If M is an integer then the range is [0,M-1] {inclusive} 
      z = y + 1;
      // z is a random integer in the range [1,M+1) {including 1, not including M+1}
      //  If M is an integer then the range is [1,M] {inclusive} 
      printf("random number %3d %5f %5f %5d %5d\n",count,r,x,y,z);
    }
} 
*/


