#include "config.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <valarray>
#include <iostream>

/* turn on a harness main() (bottom of the file) */
#define _JV_STANDALONE

// used to 'simulate' 2D arrays using linear arrays
#define sqrmat(x,y,a) ((y)*(a) + (x))

#define cpuClock() ((double)clock()/CLOCKS_PER_SEC)

using namespace std;
#include "emil.h"
#include "linearAssignment.h"


float LinearAssignment::getPairCost(int partId, int wellId){

  float     dist2, dx, dy, dz;
  Particle *p;
  Well     *w;

  p =  part[partId];
  w =  well[wellId];

  /* just do the distance calc locally..easier than navigating the 
     now-quite-bloated Hsc class structure, and also saves casting
     between doubles and floats. */
  dx = (p->R[0] - w->R[0]);
  dy = (p->R[1] - w->R[1]);
  dz = (p->R[2] - w->R[2]);
  

  //it is possible to do this more efficiently 
  if( periodic ) {
         while( dx >=   halfBox[0] ){ dx -= box[0]; }
         while( dx < -1*halfBox[0] ){ dx += box[0]; }
         while( dy >=   halfBox[1] ){ dy -= box[1]; }
         while( dy < -1*halfBox[1] ){ dy += box[1]; }
         while( dz >=   halfBox[2] ){ dz -= box[2]; }
         while( dz < -1*halfBox[2] ){ dz += box[2]; }
  } 

  dist2 = dx*dx + dy*dy + dz*dz;


  return( dist2 );
}


//  Adapted from C++ code by Rob Tyka
//  Adapted from PASCAL version in
//  Jonker, R. and Volgenant, A., 1987. A shortest augmenting
//  path algorithm for dense and sparse linear assignment problems. 
//  Computing 38, pp. 325–340

double LinearAssignment::jonkerVolgenantAssign(){

	bool found_unassigned;
	int i, i1, f0=0, lastfree, f, i0, k, f1;
	int j, j1, j2, end, last, low, up;
	int  *pred    = new int[N];	//predecessor-array for shortest path tree
	int  *free    = new int[N];  //unassigned rows (number f0, index f)
	int  *col     = new int[N];  //col array of columns, scanned
	int  *match	  = new int[N];
	int  *d       = new int[N];  // shortest path lengths

	float min, h, u1, u2, tmp1;
	float  *v       = new float[N];  // column cost
	float  *u       = new float[N];  // row cost

    for(i=0; i<N; i++) match[i] = 0;

    /* set some phony values to silence compiler warnings */
    last = -1;
    j2   = -1;
    min  = FLT_MAX;


	//column reduction
	for(j=N-1; j>=0; j--)
	{
		min = cm_lazy(0,j);
		i1 = 0;
		for(i=1; i<N; i++)
			if(cm_lazy(i,j)<min) {
				min = cm_lazy(i,j);
				i1 = i;
			}
			v[j] = min;

			if(++match[i1]==1)
			{
				partToWell[i1] = j;
				wellToPart[j] = i1;
			}else{
				wellToPart[j] = -1;
			}
	}

	// reduction transfer
	for(i=0; i<N; i++){
		if(match[i]==0){			
			free[f0++] = i;
		}else{
			if(match[i]==1)
			{
				j1 = partToWell[i];
				min = FLT_MAX;
				for(j=0; j<N; j++)
					if(j!=j1){
						if(cm_lazy(i,j)-v[j] < min){
							min = cm_lazy(i,j) - v[j];
						}
					}
				v[j1] = v[j1] - min;
			}
		}
	}

	//augmenting row reduction
	int loopcnt = 0;
	do
	{
		loopcnt++;
		k=0;
		lastfree = f0;
		f0 = 0;
		while(k<lastfree)
		{
			i = free[k];
			k++;

			u1 = cm_lazy(i,0) - v[0];
			j1 = 0;
			u2 = FLT_MAX;
			for(j=1; j<N; j++)
			{
				h = cm_lazy(i,j) - v[j];
				if(h<u2){
					if(h>=u1)
					{
						u2 = h;
						j2 = j;
					}
					else
					{
						u2 = u1;
						u1 = h;
						j2 = j1;
						j1 = j;
					}
				}
			}

			i0 = wellToPart[j1];
			if(u1 < u2){
				v[j1] = v[j1] - (u2 - u1);
			}else{
				if(i0>=0)
				{
					j1 = j2;
					i0 = wellToPart[j2];
				}
			}

			partToWell[i]  = j1;
			wellToPart[j1] = i;

			if(i0 >= 0){
				if(u1 < u2){
					free[--k] = i0;
				}else{
					free[f0++] = i0;
				}
			}
		}
	}
	while(loopcnt < 2);  // routine applied twice

	//augmentation
	for(f=0; f<f0; f++)
	{
		f1 = free[f];
		low = 0;
		up = 0;

		for(j=0; j<N; j++)
		{
			d[j] = cm_lazy(f1,j) - v[j];
			pred[j] = f1;
			col[j] = j;
		}

		found_unassigned = false;
		do
		{
			if(up==low)
			{
				last = low - 1;

				min = d[col[up++]];
				for(k=up; k<N; k++)
				{
					j = col[k];
					h = d[j];
					if(h<=min)
					{
						if(h<min)
						{
							up = low;
							min = h;
						}
						col[k] = col[up];
						col[up++] = j;
					}
				}
				for(k=low; k<up; k++){
					if(wellToPart[col[k]] < 0)
					{
						end = col[k];
						found_unassigned = true;
						break;
					}
				}
			}

			if(!found_unassigned)
			{
				j1 = col[low];
				low++;
				i = wellToPart[j1];
				h = cm_lazy(i,j1) - v[j1] - min;

				for(k=up; k<N; k++)
				{
					j = col[k];
					tmp1 = cm_lazy(i,j) - v[j] - h;
					if(tmp1<d[j])
					{
						pred[j] = i;
						if(tmp1==min){
							if(wellToPart[j]<0)
							{
								end = j;
								found_unassigned = true;
								break;
							}
							else
							{
								col[k] = col[up];
								col[up++] = j;
							}
                        }
						d[j] = tmp1;
					}
				}
			}
		}
		while(!found_unassigned);

		for(k=0; k<=last; k++)
		{
			j1 = col[k];
			v[j1] = v[j1] + d[j1] - min;
		}

		do
		{
			i = pred[end];
			wellToPart[end] = i;
			j1 = end;
			end = partToWell[i];
			partToWell[i] = j1;
		}
		while(i!=f1);
	}

	// work out final cost
	double cost = 0.0;
	for(i=0; i<N; i++)
	{
		j     = partToWell[i];
		cost += (double)cm_lazy(i,j);
	}

	delete [] pred;
	delete [] free;
	delete [] col;
	delete [] match;
	delete [] d;
	delete [] u;
	delete [] v;

	return cost;
}

#ifdef JV_STANDALONE
//driver routine to test & bench Jonker-Volgenant algorithm
int main(int argc, char* argv[]){

   const char *filename = argc < 2 ? "Asg5.in" : argv[1];
   FILE   *fptr = fopen(filename, "r");
   double start, elapsed;

   printf("Reading from %s\n", filename);
   if (fptr == NULL)
   {  puts("Failed to open file."); return 127;  }

   initialize(fptr);
   start = cpuClock();

   //do the calc
   lowerLimit = LinearAssignmentJVC(costmatrix, solution, invsolution, size);


   elapsed = cpuClock() - start;


   printf("Lowest cost: %d\n", lowerLimit);
   for( int idx = 0; idx < size; idx++)
      printf ("%3d", solution[idx]);
   putchar('\n');


   printf("%3.3f seconds\n", elapsed);
   return 0;
}
#endif











