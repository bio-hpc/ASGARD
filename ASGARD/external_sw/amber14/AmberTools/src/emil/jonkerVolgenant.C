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

#include "linearAssignment.h"


#define sqrmat(y,x,a) ((y)*(a) + (x))



//  Adapted from PASCAL version in
//  Jonker, R. and Volgenant, A., 1987. A shortest augmenting
//  path algorithm for dense and sparse linear assignment problems. 
//  Computing 38, pp. 325-340

int LinearAssignment::integerJonkerVolgenantAssign(){

        int nlist;
        int *costm, *row2col, *col2row;

        /* function-local aliases for class variables */
        nlist   = N;
        costm   = intCostMatrix;
        row2col = wellToPart;
        col2row = partToWell;


        bool found_unassigned;
        int i, i1, f0=0, lastfree, f, i0, k, f1;
        int j, j1, j2, end, last, low, up;
        int min, h, u1, u2, tmp1;
        int *pred    = new int[nlist];        //predecessor-array for shortest path tree
        int *free    = new int[nlist]; //unassigned rows (number f0, index f)
        int *col     = new int[nlist]; //col array of columns, scanned
        int *match   = new int[nlist];
        int *d       = new int[nlist];      // shortest path lengths
        int *u       = new int[nlist];   // row cost
        int *v       = new int[nlist];   // column cost

        //needless initialisations to avoid compiler warnings
        last = 0;
        j2   = 0;
        min  = 0;

        for(i=0; i<nlist; i++) match[i] = 0;

        //column reduction
        for(j=nlist-1; j>=0; j--)
        {
                min = costm[sqrmat(0,j,nlist)];
                i1 = 0;
                for(i=1; i<nlist; i++)
                        if(costm[sqrmat(i,j,nlist)]<min)
                        {
                                min = costm[sqrmat(i,j,nlist)];
                                i1 = i;
                        }
                        v[j] = min;

                        if(++match[i1]==1)
                        {
                                col2row[i1] = j;
                                row2col[j] = i1;
                        }else{
                                row2col[j] = -1;
                        }
        }

        // reduction transfer
        for(i=0; i<nlist; i++){
                if(match[i]==0){                        
                        free[f0++] = i;
                }else{
                        if(match[i]==1)
                        {
                                j1 = col2row[i];
                                min = INT_MAX;
                                for(j=0; j<nlist; j++)
                                        if(j!=j1){
                                                if(costm[sqrmat(i,j,nlist)]-v[j] < min){
                                                        min = costm[sqrmat(i,j,nlist)] - v[j];
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

                        u1 = costm[sqrmat(i,0,nlist)] - v[0];
                        j1 = 0;
                        u2 = INT_MAX;
                        for(j=1; j<nlist; j++)
                        {
                                h = costm[sqrmat(i,j,nlist)] - v[j];
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

                        i0 = row2col[j1];
                        if(u1 < u2){
                                v[j1] = v[j1] - (u2 - u1);
                        }else{
                                if(i0>=0)
                                {
                                        j1 = j2;
                                        i0 = row2col[j2];
                                }
                        }

                        col2row[i] = j1;
                        row2col[j1] = i;

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

                for(j=0; j<nlist; j++)
                {
                        d[j] = costm[sqrmat(f1,j,nlist)] - v[j];
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
                                for(k=up; k<nlist; k++)
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
                                        if(row2col[col[k]] < 0)
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
                                i = row2col[j1];
                                h = costm[sqrmat(i,j1,nlist)] - v[j1] - min;

                                for(k=up; k<nlist; k++)
                                {
                                        j = col[k];
                                        tmp1 = costm[sqrmat(i,j,nlist)] - v[j] - h;
                                        if(tmp1<d[j])
                                        {
                                                pred[j] = i;
                                                if(tmp1==min){
                                                        if(row2col[j]<0)
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
                        row2col[end] = i;
                        j1 = end;
                        end = col2row[i];
                        col2row[i] = j1;
                }
                while(i!=f1);
        }

        // work out final cost
        int cost = 0;
        for(i=0; i<nlist; i++)
        {
                j = col2row[i];
                u[i] = costm[sqrmat(i,j,nlist)] - v[j];
                cost += costm[sqrmat(i,j,nlist)];
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
