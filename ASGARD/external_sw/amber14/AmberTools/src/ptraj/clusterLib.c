/*  _______________________________________________________________________
 *
 *                        RDPARM/PTRAJ: 2008
 *  _______________________________________________________________________
 *
 *  This file is part of rdparm/ptraj.
 *
 *  rdparm/ptraj is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  rdparm/ptraj is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You can receive a copy of the GNU General Public License from
 *  http://www.gnu.org or by writing to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *  ________________________________________________________________________
 *
 *  CVS tracking:
 *
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/clusterLib.c,v 10.0 2008/04/15 23:24:11 case Exp $
 *
 *  Revision: $Revision: 10.0 $
 *  Date: $Date: 2008/04/15 23:24:11 $
 *  Last checked in by $Author: case $
 *  ________________________________________________________________________
 *
 *
 *  CONTACT INFO: To learn who the code developers are, who to contact for
 *  more information, and to know what version of the code this is, refer
 *  to the CVS information and the include files (contributors.h && version.h)
 *
 */

#include "contributors.h"
#include "version.h"

/*  ________________________________________________________________________
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <time.h>

#define CLUSTERLIB_MODULE
#include "ptraj.h"


/* Exporting magic, for Windows builds: */
#define CLUSTER_EXPORTS


/*
 *  Functions pertaining to the cluster struct follow - these are
 *  the "members" of the Cluster "class". 
 */

  /*
   *  Returns the number of points available for clustering (NOT the 
   *  number of members in this cluster) 
   */
   int
ClusterGetPointCount(Cluster* This)
{
  return This->PointCount;
}

  /*
   *  Returns true if the specified point is in this cluster 
   */
   int
ClusterIsMember(Cluster* This, int PointIndex)
{
  return This->Mask[PointIndex];
}

  /*
   *  Adds the specified point to this cluster 
   */
   void
ClusterAddMember(Cluster* This, int PointIndex)
{
  This->Mask[PointIndex]=1;
}

  /*
   *  Add all points to this cluster 
   */
   void
ClusterAddEverything(Cluster* This)
{
  int PointIndex;
  for (PointIndex=0; PointIndex<This->PointCount; PointIndex++)
  {
    ClusterAddMember(This,PointIndex);
  }
}

  /*
   *  Get the index of the point being used as this cluster's centroid.
   *  (Most subclasses build a new centroid point, and so don't use this method) 
   */
   int
ClusterGetCentroidIndex(Cluster* This)
{
  return This->CentroidIndex;
}

  /*
   *  Remove a point from this cluster 
   */
   void
ClusterRemoveMember(Cluster* This, int PointIndex)
{
  This->Mask[PointIndex]=0;
}

  /*
   *  Allocate a new Cluster object; free with ClusterFree 
   */
   Cluster*
ClusterNew(Clustering* Owner)
{
  Cluster* pCluster;
  pCluster = (Cluster*)malloc(sizeof(Cluster));
  memset(pCluster,0,sizeof(Cluster));
  pCluster->Owner=Owner;
  pCluster->SSEWithin=0;
  pCluster->PointCount = ClusteringGetPointCount(Owner);
  pCluster->Mask = (int*)malloc(sizeof(int) * pCluster->PointCount);
  memset(pCluster->Mask,0,sizeof(int) * pCluster->PointCount);
  pCluster->CentroidIndex=0;
  pCluster->VirtualFunctions = ClusterVirtualFunctions;
  return pCluster;
}

  /*
   *  Free a Cluster object 
   */
   void
BASEClusterFree(Cluster* This)
{
  safe_free(This->Mask);
  safe_free(This);
}

  /* 
   *  Find the centroid of this cluster.  After altering the cluster's 
   *  membership (adding and removing points), you must call ClusterFindCentroid
   *  to keep the centroid up-to-date 
   */
   void
BASEClusterFindCentroid(Cluster* This)
{
    float Closest = 0; /* Initialization value */
    float ThisPointDistance;
    int PointIndex;
    int PointIndexB;
    SymmetricMatrix* PairwiseDistances;
    /**/
    /* For this base-class implementation: Look for the point whose 
       CUMULATIVE DISTANCE (sum of distances to all other points in
       the cluster) is as small as possible.  Iterate over all pairs
       of points.  Store the index in This->CentroidIndex, and the smallest
       cumulative distance so far in Closest.*/
    PairwiseDistances = ClusteringGetPairwiseDistances(This->Owner);
    for (PointIndex=0; PointIndex<This->PointCount; PointIndex++)
    {
        if (!This->Mask[PointIndex])
        {
            continue;
        }
        ThisPointDistance=0;
        for (PointIndexB=0; PointIndexB<This->PointCount; PointIndexB++)
        {
            if (!This->Mask[PointIndexB])
            {
                continue;
            }
            ThisPointDistance += ClusteringPointToPointDistance(This->Owner,\
                PointIndex,PointIndexB);
        }
        if (Closest==0 || ThisPointDistance<Closest)
        {
            Closest = ThisPointDistance;
            This->CentroidIndex = PointIndex;
        }
    }
}

  /*
   *  Return the number of points in this cluster 
   */
   CLUSTER_API int
ClusterCountPointsInCluster(Cluster* This)
{
    int Count = 0;
    int PointIndex;
    for (PointIndex=0; PointIndex<This->PointCount; PointIndex++)
    {
        Count += This->Mask[PointIndex];
    }
    return Count;
}

  /*
   *  !!! Cluster VTable stuff !!!
   */
   void
ClusterFindCentroid(Cluster* This)
{
  This->VirtualFunctions[F_CLUSTER_FIND_CENTROID](This);
}

   void
ClusterFree(Cluster* This)
{
  This->VirtualFunctions[F_CLUSTER_FREE](This);
}


  /*
   *  Functions pertaining to the clustering struct follow - these are
   *  the "members" of the Clustering "class". 
   */

  /*
   *  Allocate a new cluster object for this clustering.  (Afterwards, you still
   *  must call ClusteringAddClustering if you want the Cluster to our children) 
   */
   Cluster*
BASEClusteringGetNewCluster(Clustering* This)
{
  return ClusterNew(This);
}

  /*
   *  Returns the PairwiseDistances (set at initialization time) for a Clustering 
   */
   SymmetricMatrix*
ClusteringGetPairwiseDistances(Clustering* This)
{
  return This->PairwiseDistances;
}

   /*
    *  Returns the number of points being clustered 
    */
   int
ClusteringGetPointCount(Clustering* This)
{
  return This->PointCount;
}

   /*
    *  Allocates a new Clustering object.  Free with ClusteringFree. 
    */
   Clustering*
ClusteringNew(SymmetricMatrix* PairwiseDistances)
{
  Clustering* This;
  This = (Clustering*)malloc(sizeof(Clustering));
  memset(This,0,sizeof(Clustering));
  This->Head=NULL;
  This->Tail=NULL;
  This->PairwiseDistances = PairwiseDistances;
  This->PointCount = PairwiseDistances->Size;
  This->VirtualFunctions = ClusteringVirtualFunctions;
  This->FloatVirtualFunctions = ClusteringFloatVirtualFunctions;
  This->IntVirtualFunctions = ClusteringIntVirtualFunctions;
  This->attributeArray = NULL;
  This->attributeArrayTorsion = NULL;
  return This;
}

  /*
   *  Frees a Clustering object.  Also frees its clusters! 
   */
   void
BASEClusteringFree(Clustering* This)
{
    ClusterNode* Node;
    ClusterNode* Prev = NULL;
    /**/
    for (Node = This->Head; Node; Node = Node->pNext)
    {
        if (Prev)
        {
            ClusterFree(Prev->Cluster);
            safe_free(Prev);
        }
        Prev = Node;
    }
    if (Prev)
    {
        ClusterFree(Prev->Cluster);
        safe_free(Prev);
    }
    safe_free(This);
}

  /*
   *  Adds a Cluster (probably created by ClusteringGetNewCluster) to our children 
   */
   ClusterNode*
ClusteringAddCluster(Clustering* This, Cluster* NewCluster)
{
    ClusterNode* NewNode;
    /**/
    NewNode = (ClusterNode*)malloc(sizeof(ClusterNode));
    memset(NewNode,0,sizeof(ClusterNode));
    NewNode->Cluster = NewCluster;
    NewCluster->Owner = This;
    if (This->Tail)
    {
        This->Tail->pNext = NewNode;
        NewNode->pPrev = This->Tail;
    }
    if (!This->Head)
    {
        This->Head = NewNode;
    }
    This->Tail = NewNode;
    This->ClusterCount++;
    return NewNode;
}


  /*
   *  Removes a cluster from our list of children.  (e.g. if merging two 
   *  clusters into one, get rid of an old one).  Does *not* free the 
   *  old cluster, just removes its node. 
   */
   void
ClusteringRemoveCluster(Clustering* This, Cluster* OldCluster)
{
    /* Frees the cluster node, but does NOT free the cluster*/
    ClusterNode* TempNode;
    /**/
    for (TempNode = This->Head; TempNode; TempNode = TempNode->pNext)
    {
        if (TempNode->Cluster == OldCluster)
        {
            if (TempNode->pPrev)
            {
                TempNode->pPrev->pNext = TempNode->pNext;
            }
            if (TempNode->pNext)
            {
                TempNode->pNext->pPrev = TempNode->pPrev;
            }
            if (TempNode == This->Head)
            {
                This->Head = TempNode->pNext;
            }
            if (TempNode == This->Tail)
            {
                This->Tail = TempNode->pPrev;
            }
            safe_free(TempNode);
            This->ClusterCount--;
            return;
        }
    }
}

  /*
   *  Remove empty clusters. 
   */
   void
ClusteringRemoveEmptyCluster(Clustering* This)
{
	ClusterNode* Node, *TempNode;
    int ClusterIndex;
    
   	for (Node = This->Head,ClusterIndex=0; Node; ClusterIndex++)
    {	
   		if (ClusterCountPointsInCluster(Node->Cluster) == 0 )
   		{
            TempNode = Node->pPrev;
           	ClusteringRemoveCluster(This, Node->Cluster);
            if (TempNode) 
            	Node = TempNode->pNext;
            else 
            	Node = This->Head;
            fprintf(stdout, "Cluster %d does not contain any points. Delete it. \n", ClusterIndex);
   		} else {
        	Node = Node->pNext;
        }
    }
}

  /*
   *  Swap the order of two cluster next to each, used in the bubble sorting 
   */
   void
ClusteringSwapCluster(Clustering* This, ClusterNode* Node)
{
	ClusterNode* NextNode;
    /**/
    NextNode = Node->pNext;
    if (NextNode == NULL) 
    {
    	return;
    }
    if (This->Head == Node) 
    {
    	This->Head = NextNode;
    }
    else
    {
    	Node->pPrev->pNext = NextNode;
    }
    if (This->Tail == NextNode) 
    {
    	This->Tail = Node;
    }
    else
    {
    	NextNode->pNext->pPrev = Node;
    }
    Node->pNext = NextNode->pNext;
    NextNode->pPrev = Node->pPrev;
    Node->pPrev = NextNode;
    NextNode->pNext = Node;
}

  /*
   *  Bubble sorting of the Cluster, largest to smallest 
   */
   void
ClusteringSortCluster(Clustering* This)
{
	ClusterNode* Node;
	ClusterNode* Node2;
	ClusterNode* SmallestNode = NULL;
    int points1, points2;
    /**/
    while(SmallestNode != This->Head->pNext)
    {
	    for (Node = This->Head; Node->pNext && Node->pNext != SmallestNode;)
	    {
	        points1 = ClusterCountPointsInCluster(Node->Cluster);
	        points2 = ClusterCountPointsInCluster(Node->pNext->Cluster);
	        if (points1 < points2) 
	        {
	        	ClusteringSwapCluster(This, Node);
	   	    }
	        else
	        {
	        	Node = Node->pNext;
	        }
	    }
	    SmallestNode = Node;
    }
}

  /*
   *  Iterate over each cluster, and call FindCentroid on each one.  Useful if you
   *  have just preformed major restructuring on them. 
   */
   void
ClusteringFindAllCentroids(Clustering* This)
{
    ClusterNode* Node;
    for (Node = This->Head; Node; Node=Node->pNext)
    {
        ClusterFindCentroid(Node->Cluster);
    }
}

#define BLOCK_SIZE 20480

/*
 *  Returns the length of the first line in a file.  (It's assumed that
 *  all lines have the same length; we call this when dealing with clustering
 *  output).  The line-length should be the same across files, when you 
 *   compare clustering results! 
 */
   int
GetLineLength(FILE* SourceFile)
{
  char Temp[BLOCK_SIZE];
  int Size = 0;
  char* Result;
  while (1)
  {
    Result = fgets(Temp,BLOCK_SIZE,SourceFile);
    if (!Result || strlen(Result)==0)
    {
      break;
    }
    /* Skip comment-lines: */
    if (Result[0]=='#')
    {
      continue;
    }
    Size += strlen(Result);
    if (Result[strlen(Result)-1]=='\n')
    {
      Size -= 1;
      break;
    }
  }
  fseek(SourceFile,0,0);
  return Size;
}

  /*
   *  Read a clustering from a text-file.  Re-reading the clustering won't 
   *  re-read the point definitions; we only get the membership info.  
   *  Useful for comparing the results of clustering algorithms. 
   */
   void
ClusteringReadFromDisk(Clustering* This,char* SourceFilePath)  
{
    char* FileLine;
    int PointCount;
    char* Temp;
    int PointIndex;
    Cluster* NewCluster;
    char* Result;
    FILE* SourceFile;
    int LineSize;
    SourceFile = fopen(SourceFilePath,"r");
    if (!SourceFile) {
      fprintf(stderr, "Cannot open %s to read\n", SourceFilePath);
    }
    PointCount = This->PointCount;
    LineSize = GetLineLength(SourceFile);
    if (LineSize != PointCount)
    {
      fprintf(stdout,"Warning!  Line length %d doesn't match point count %d for file %s\n",
	      LineSize,PointCount,SourceFilePath);
    }
    /*fprintf(stdout,"LineLength %d, point count %d\n",LineSize,PointCount);*/
    FileLine = (char*)SafeMalloc(__FILE__, __LINE__, sizeof(char)*(BLOCK_SIZE));
    while (1)
    {
        Result = fgets(FileLine, BLOCK_SIZE ,SourceFile);
        if (!Result)
        {
            return;
        }
        if (strlen(FileLine)<1)
        {
            break;
        }
        if (Result[0]=='#')
        {
            continue;
        }
	if (FileLine[0] == '\r' || FileLine[0] == '\n' || FileLine[0]=='\0')
        {
            continue;
        }
	/*fprintf(stdout,"Here is a file line:\n%s",Result);*/
        NewCluster = ClusteringGetNewCluster(This);
        PointIndex=0;
        for (Temp = FileLine; *Temp!='\n' && *Temp!='\0' && *Temp!='\r' && PointIndex<PointCount; Temp++)
        {
            if (*Temp=='X')
            {
              ClusterAddMember(NewCluster,PointIndex);
            }
            if (*Temp != 'X' && *Temp != '.')
            {
              fprintf(stdout,"Warning!  Bogus character %c encountered in cluster output file %s\n",
		      *Temp,SourceFilePath);
            }
            PointIndex++;
        }
        if (PointIndex!=1 && PointIndex!=PointCount)
        {
          fprintf(stdout,"Warning!  Line length %d instead of %d found in cluster output file %s\n",
		  PointIndex,PointCount,SourceFilePath);
        }
	/*fprintf(stdout,"Adding one cluster from disk.\n");*/
        ClusteringAddCluster(This,NewCluster);
    }
    fclose(SourceFile);
    free(FileLine);
}

  /*
   *  Return an array whose nth entry is the cluster index that the nth 
   *  point has been assigned to. Used for comparing results.  Caller
   *  frees the array.
   */
   int*
ClusteringGetClusterNumberArray(Clustering* This)
{
    int* Array;
    int PointIndex;
    int PointCount;
    ClusterNode* Node;
    int ClusterIndex=0;
    /**/
    PointCount = ClusteringGetPointCount(This);
    Array = (int*)malloc(sizeof(int)*PointCount);
    for (Node = This->Head; Node; Node = Node->pNext)
    {
        for (PointIndex=0; PointIndex<PointCount; PointIndex++)
        {
              if (ClusterIsMember(Node->Cluster,PointIndex))
            {
                Array[PointIndex]=ClusterIndex;
            }
        }
        ClusterIndex++;
    }
    return Array; /* caller frees it */
}

  /*
   *  Helper function for comparing results: Assuming random distributions 
   *  (that is, using only our cluster sizes), what are the odds that two 
   *  points fall in the same cluster? 
   */
   float
ClusteringGetOddsSameCluster(Clustering* This)
{
    float Odds = 0;
    int Count;
    int PointCount = This->PointCount;
    ClusterNode* Node;
    /* The odds that a point falls in cluster A is equal to the number 
       of points in A divided by the total points (Count/PointCount).  The
       odds that two points fall in cluster A equals (Count/PointCount)**2.
       The odds that two points share any cluster is the sum of such 
       values. */
    for (Node = This->Head; Node; Node = Node->pNext)
    {
        Count = ClusterCountPointsInCluster(Node->Cluster);
        Odds += (Count*Count / (float)(PointCount*PointCount));
    }
    return Odds;
}

/*
 *
 * DATA: Chi-squared table for 1 degree of freedom (2 categories) 
 */

/*float ChiSquaredThreshold[] = {1.32, 1.64, 2.07, 2.71, 3.84, 5.02, 5.41, 6.63, 7.88, 9.14, 10.83, 12.12};
  float ChiSquaredPValue[] = {0.25, 0.20, 0.15, 0.10, 0.05, 0.025, 0.02, 0.01, 0.005, 0.0025, 0.001, 0.0005};*/

float ChiSquaredThreshold[] = {0.001, 0.005, 0.01, 0.02, 0.05, 
                               0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
			       1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
			       11.0, 12.0, 13.0};
float ChiSquaredPValue[] = {0.975, 0.944, 0.920, 0.888, 0.823,
			    0.752, 0.655, 0.584, 0.527, 0.480, 0.439, 0.403, 0.371, 0.343, 0.317,
			    0.221, 0.157, 0.0833, 0.0455, 0.0254, 0.0143, 0.00815, 0.00468,
			    0.0027, 0.00157, 0.00091, 0.00053, 0.00031};

int ChiSquaredCount = sizeof(ChiSquaredPValue) / sizeof(float);
float ChiSquaredMax = 0.999; /* Higher than any p-value in our table */
float ChiSquaredMin = 0.0002; /* Lower than any p-value in our table */

  /*
   *  Measure the distance between two clusters using chi-squared values.  We divide all
   *  pairs frames into two classes: the "agreement" and "disagreement" classes.  A frame-pair falls
   *  in the "agreement" class if (a) both clusterings put the frames in the same cluster, or (b) both
   *  clusterings put them in different clusters.  Otherwise it's a "diagreement" pair. 
   */
   float
ClusteringMeasureChiSquared(Clustering* This, Clustering* Other)
{
    int* NumbersA;
    int* NumbersB;
    int PointIndex;
    int OtherPointIndex;
    int PointCount;
    int MyA;
    int OtherA;
    int MyB;
    int OtherB;
    int Disagreements;
    int Agreements;
    float OddsSameA;
    float OddsSameB;
    float ExpectedDisagreements;
    float Deviation;
    int PairCount = 0;
    float XSquared;
    int ChiIndex;
    /**/
    PointCount = This->PointCount;
    NumbersA = ClusteringGetClusterNumberArray(This);
    NumbersB = ClusteringGetClusterNumberArray(Other);
    Disagreements = 0;
    Agreements = 0;
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        for (OtherPointIndex=PointIndex+1; OtherPointIndex<PointCount; OtherPointIndex++)
        {
            MyA = NumbersA[PointIndex];
            MyB = NumbersA[OtherPointIndex];
            OtherA = NumbersB[PointIndex];
            OtherB = NumbersB[OtherPointIndex];
            if ((MyA==MyB && OtherA!=OtherB) || (MyA!=MyB && OtherA==OtherB))
            {
              Disagreements++;
            }
	    else
            {
	      Agreements++;
	    }
            PairCount++;
        }
    }
    OddsSameA = ClusteringGetOddsSameCluster(This);
    OddsSameB = ClusteringGetOddsSameCluster(Other);
      /*
       *  How many disagreements do we expect to see?  The odds of a 
       *  disagreement are:
       *    P(Same in A)*P(different in B) + P(Different in A)*P(Same in B)
       *  Expected number of disagreements will be this probability times the
       *  number of point-pairs 
       */
    ExpectedDisagreements = PairCount * ((OddsSameA * (1-OddsSameB)) + (OddsSameB * (1-OddsSameA)));
    fprintf(stdout,"OddsSameA:%.4f OddsSameB: %.4f PairCount:%d",OddsSameA, OddsSameB, PairCount);
    Deviation = ExpectedDisagreements - Disagreements;
    Deviation = Deviation*Deviation / ExpectedDisagreements;
    XSquared = abs(Deviation) - 0.5;
    XSquared = max(0,XSquared);
    XSquared *= 2;

    fprintf(stdout,"ExpectedDis %.4f, ActualDis %d\n",ExpectedDisagreements, Disagreements);
    fprintf(stdout,"XSquared: %f\n",XSquared);
    safe_free(NumbersA);
    safe_free(NumbersB);
    /* Now get a p-value */
    if (XSquared < ChiSquaredThreshold[0])
    {
      fprintf(stdout,"ChiSquared p-value: %.4f", ChiSquaredMax);
      return ChiSquaredMax;
    }
    for (ChiIndex=0; ChiIndex<ChiSquaredCount; ChiIndex++)
    {
      if (ChiSquaredThreshold[ChiIndex]>XSquared)
      {	
	fprintf(stdout,"ChiSquared p-value: %.4f", ChiSquaredPValue[ChiIndex-1]);
	return ChiSquaredPValue[ChiIndex-1];
      }
    }
	fprintf(stdout,"ChiSquared p-value: %.4f", ChiSquaredMin);
    return ChiSquaredMin; 
  
}

  /*
   *  Measure the distance between this cluster and another.  Our distance 
   *  metric is defined as follows: Say two clusterings "agree" on a pair 
   *  of points if they  both assign them to the same cluster, or if they 
   *  both assign them to different clusters.  We iterate over all pairs 
   *  of points, and compute the number of times our two clusterings 
   *  disagree.  Then, we take the ratio of this sum with the "expected" 
   *  number of disagreements.  Two random clusterings will have a distance 
   *  near 1.0; a clustering has distance 0.0 from itself; it is possible 
   *  (though not likely) to have distances >1 
   */
   float
ClusteringMeasureDistance(Clustering* This,Clustering* Other)
{
    int* NumbersA;
    int* NumbersB;
    int PointIndex;
    int OtherPointIndex;
    int PointCount;
    int MyA;
    int OtherA;
    int MyB;
    int OtherB;
    int Distance = 0;
    float OddsSameA;
    float OddsSameB;
    float ExpectedDistance;
    int PairCount = 0;
    /**/
    PointCount = This->PointCount;
    NumbersA = ClusteringGetClusterNumberArray(This);
    NumbersB = ClusteringGetClusterNumberArray(Other);
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        for (OtherPointIndex=PointIndex+1; OtherPointIndex<PointCount; OtherPointIndex++)
        {
            MyA = NumbersA[PointIndex];
            MyB = NumbersA[OtherPointIndex];
            OtherA = NumbersB[PointIndex];
            OtherB = NumbersB[OtherPointIndex];
            if ((MyA==MyB && OtherA!=OtherB) || (MyA!=MyB && OtherA==OtherB))
            {
              Distance++; /* This is one disagreement! */
            }
            PairCount++;
        }
    }
    OddsSameA = ClusteringGetOddsSameCluster(This);
    OddsSameB = ClusteringGetOddsSameCluster(Other);
    /* How many disagreements do we expect to see?  The odds of a 
    disagreement are:
    P(Same in A)*P(different in B) + P(Different in A)*P(Same in B)
    Expected number of disagreements will be this probability times the
    number of point-pairs */
    ExpectedDistance = PairCount * ((OddsSameA * (1-OddsSameB)) + (OddsSameB * (1-OddsSameA)));
    fprintf(stdout,"Distance: %d vs %d expected %f return %f\n",Distance,PairCount,
	    ExpectedDistance,Distance/ExpectedDistance);
    safe_free(NumbersA);
    safe_free(NumbersB);
    return Distance / ExpectedDistance;
}

/*
  DATA:

          p = 0.40      0.25       0.10     0.05       0.025    0.01    0.005    0.0005
   df      
*/
float TTestPValue[8] = {0.40, 0.25,   0.10,     0.05,    0.025,    0.01,    0.005,   0.0005};
float TScore[31][8] = {
/*  1  */  {0.324920 ,1.000000 ,3.077684 ,6.313752 ,12.70620 ,31.82052 ,63.65674 ,636.6192},
/*  2  */  {0.288675 ,0.816497 ,1.885618 ,2.919986 ,4.30265 ,6.96456 ,9.92484 ,31.5991},
/*  3  */  {0.276671 ,0.764892 ,1.637744 ,2.353363 ,3.18245 ,4.54070 ,5.84091 ,12.9240},
/*  4  */  {0.270722 ,0.740697 ,1.533206 ,2.131847 ,2.77645 ,3.74695 ,4.60409 ,8.6103},
/*  5  */  {0.267181 ,0.726687 ,1.475884 ,2.015048 ,2.57058 ,3.36493 ,4.03214 ,6.8688},
/*  6  */  {0.264835 ,0.717558 ,1.439756 ,1.943180 ,2.44691 ,3.14267 ,3.70743 ,5.9588},
/*  7  */  {0.263167 ,0.711142 ,1.414924 ,1.894579 ,2.36462 ,2.99795 ,3.49948 ,5.4079},
/*  8  */  {0.261921 ,0.706387 ,1.396815 ,1.859548 ,2.30600 ,2.89646 ,3.35539 ,5.0413},
/*  9  */  {0.260955 ,0.702722 ,1.383029 ,1.833113 ,2.26216 ,2.82144 ,3.24984 ,4.7809},
/* 10  */  {0.260185 ,0.699812 ,1.372184 ,1.812461 ,2.22814 ,2.76377 ,3.16927 ,4.5869},
/* 11  */  {0.259556 ,0.697445 ,1.363430 ,1.795885 ,2.20099 ,2.71808 ,3.10581 ,4.4370},
/* 12  */  {0.259033 ,0.695483 ,1.356217 ,1.782288 ,2.17881 ,2.68100 ,3.05454 ,4.3178},
/* 13  */  {0.258591 ,0.693829 ,1.350171 ,1.770933 ,2.16037 ,2.65031 ,3.01228 ,4.2208},
/* 14  */  {0.258213 ,0.692417 ,1.345030 ,1.761310 ,2.14479 ,2.62449 ,2.97684 ,4.1405},
/* 15  */  {0.257885 ,0.691197 ,1.340606 ,1.753050 ,2.13145 ,2.60248 ,2.94671 ,4.0728},
/* 16  */  {0.257599 ,0.690132 ,1.336757 ,1.745884 ,2.11991 ,2.58349 ,2.92078 ,4.0150},
/* 17  */  {0.257347 ,0.689195 ,1.333379 ,1.739607 ,2.10982 ,2.56693 ,2.89823 ,3.9651},
/* 18  */  {0.257123 ,0.688364 ,1.330391 ,1.734064 ,2.10092 ,2.55238 ,2.87844 ,3.9216},
/* 19  */  {0.256923 ,0.687621 ,1.327728 ,1.729133 ,2.09302 ,2.53948 ,2.86093 ,3.8834},
/* 20  */  {0.256743 ,0.686954 ,1.325341 ,1.724718 ,2.08596 ,2.52798 ,2.84534 ,3.8495},
/* 21  */  {0.256580 ,0.686352 ,1.323188 ,1.720743 ,2.07961 ,2.51765 ,2.83136 ,3.8193},
/* 22  */  {0.256432 ,0.685805 ,1.321237 ,1.717144 ,2.07387 ,2.50832 ,2.81876 ,3.7921},
/* 23  */  {0.256297 ,0.685306 ,1.319460 ,1.713872 ,2.06866 ,2.49987 ,2.80734 ,3.7676},
/* 24  */  {0.256173 ,0.684850 ,1.317836 ,1.710882 ,2.06390 ,2.49216 ,2.79694 ,3.7454},
/* 25  */  {0.256060 ,0.684430 ,1.316345 ,1.708141 ,2.05954 ,2.48511 ,2.78744 ,3.7251},
/* 26  */  {0.255955 ,0.684043 ,1.314972 ,1.705618 ,2.05553 ,2.47863 ,2.77871 ,3.7066},
/* 27  */  {0.255858 ,0.683685 ,1.313703 ,1.703288 ,2.05183 ,2.47266 ,2.77068 ,3.6896},
/* 28  */  {0.255768 ,0.683353 ,1.312527 ,1.701131 ,2.04841 ,2.46714 ,2.76326 ,3.6739},
/* 29  */  {0.255684 ,0.683044 ,1.311434 ,1.699127 ,2.04523 ,2.46202 ,2.75639 ,3.6594},
/* 30  */  {0.255605 ,0.682756 ,1.310415 ,1.697261 ,2.04227 ,2.45726 ,2.75000 ,3.6460},
/*inf  */  {0.253347 ,0.674490 ,1.281552 ,1.644854 ,1.95996 ,2.32635 ,2.57583 ,3.2905}
};

   float
ClusteringTscore2PValue(float tscore, int df)
{
	int i, j, dfIndex;
    
    if (df <1)
    {
    	fprintf(stderr, "df is less than 1.\n");
        return(-1);
    }
    dfIndex = (df > 31) ? 30 : df-1;
    for (i=0; i<8; i++) 
    {
    	if(TScore[dfIndex][i] > tscore) 
        {
            break;
        }
    }
    return(TTestPValue[i-1]*2);
}

void
ClusteringMeasureTScore(Clustering* This, ClusterNode* ClusterNodeA, ClusterNode* ClusterNodeB) 
{
	float SSEA, SSEB;
    int NA, NB;
    float t;
    int dfA, dfB, dfIndex;
    float df2;
    float s2A, s2B;
    float s2AN, s2BN;
    int i;
    float t1, t2, t3;
    float dist;
    Cluster* ClusterA = ClusterNodeA->Cluster;
    Cluster* ClusterB = ClusterNodeB->Cluster;
    
    dist = ClusteringCentroidToCentroidDistance(This, ClusterNodeA, ClusterNodeB);
    SSEA = ClusterA->SSEWithin;
    SSEB = ClusterB->SSEWithin;
    NA = ClusterCountPointsInCluster(ClusterA);
    NB = ClusterCountPointsInCluster(ClusterB);
    dfA = NA - 1;
    dfB = NB - 1;
    dfIndex = dfA + dfB - 1; 
    
    s2A = SSEA / dfA;
    s2B = SSEB / dfB;
    if (dfA == 0) 
    {
    	printf("Only one frame in NodeA\n");
        return;
    }
    if (dfB == 0) 
    {
    	printf("Only one frame in NodeB\n");
        return;
    }
    if (dfA == 0) s2A = 0;
    if (dfB == 0) s2B = 0;
    s2AN = s2A / NA;
    s2BN = s2B / NB;
    df2 = (s2AN + s2BN)*(s2AN + s2BN)/(s2AN*s2AN/dfA + s2BN*s2BN/dfB);
    
    t1 = dist/sqrt(s2A/NA + s2B/NB);
    t2 = dist/sqrt((SSEA+SSEB)/(dfA+dfB)*(1/(NA+0.0)+1/(NB+0.0)));
    t3 = dist/sqrt(This->SSE/(This->PointCount-This->ClusterCount));
    
    printf("SSEA %f SSEB %f  ", SSEA, SSEB);
    printf("s2A %f s2B %f  ", s2A, s2B);
    printf("NA %d NB %d  ", NA, NB);
    printf("s2AN %f s2BN %f  ", s2AN, s2BN);
    printf("df: %d; df2: %f\n", dfA+dfB, df2);
    printf("dist %f t1: %f; t2: %f; t3: %f\n", dist, t1, t2, t3);
    printf("p1: %f; p2: %f; p3: %f\n",ClusteringTscore2PValue(t1,dfA+dfB), 
	   ClusteringTscore2PValue(t2,dfA+dfB), ClusteringTscore2PValue(t3,dfA+dfB));
    
}

  /* 
   *  Print out a clustering to stdout, for debugging purposes.  Prints 
   *  lines of the form "...X..X", where an X in the xth position on 
   *  the yth row means that point x is a member of cluster y 
   */
   void
ClusteringDebugPrintClusterList(Clustering* This)
{
    ClusteringOutputClusterListToFile(This,stdout,1);
}



  /*
   *  Creates a file at the specified path, and calls ClusteringOutputClusterListToFile 
   */
   void
ClusteringOutputClusterList(Clustering* This,char* FilePath)
{
    FILE* OutFile;
    OutFile = fopen(FilePath,"w");
    ClusteringOutputClusterListToFile(This,OutFile,0);
    ClusteringOutputClusterConsensusToFile(This,OutFile,1);
    ClusteringOutputClusterCondenseToFile(This,OutFile);
    fclose(OutFile);
}

  /*
   *  Calls ClusteringOutputClusterListToFile, *appending* to the specified file 
   */
   void
ClusteringAppendClusterList(Clustering* This,char* FilePath)
{
    FILE* OutFile;
    OutFile = fopen(FilePath,"a");
    ClusteringOutputClusterListToFile(This,OutFile,0);
    ClusteringOutputClusterConsensusToFile(This,OutFile,1);
    ClusteringOutputClusterCondenseToFile(This,OutFile);
    fclose(OutFile);
}

  /*
   *  Prints a cluster list to a file.  If IncludeSeparators is true, 
   *  include extra lines (which we can't parse later) for better 
   *  readability 
   */
   void
ClusteringOutputClusterListToFile(Clustering* This, FILE* OutFile,int IncludeSeparators)
{
    ClusterNode* Node;
    int PointCount;
    int PointIndex;
 
    /**/
    if (IncludeSeparators)
    {
        fprintf(OutFile,"Clustering: %d clusters\n",This->ClusterCount);
    }
    fprintf(OutFile,"#Clustering: %d clusters\n",This->ClusterCount);
    PointCount = This->PointCount;
    for (Node=This->Head; Node; Node = Node->pNext)
    {
        for (PointIndex=0; PointIndex<PointCount; PointIndex++)
        {
            if (ClusterIsMember(Node->Cluster,PointIndex))
            {
                fprintf(OutFile,"X");
            }
            else
            {
                fprintf(OutFile,".");
            }
        }
        fprintf(OutFile,"\n");
    }
    
}


  /* */
   void
ClusteringOutputClusterCondenseToFile(Clustering* This, FILE* OutFile)
{
    ClusterNode* Node;
    int PointCount = This->PointCount;
    int PointIndex, ClusterIndex;
	float* freq;
    int partition = 50;
    float WindowLength, fPos;   
    int currentWindow,currentLength;  
    int i, ci, iPos;
    char c;
    
    /**/	
    WindowLength = (PointCount+0.0)/partition;
    if (PointCount < partition) WindowLength = 1;
    freq = (float*)safe_malloc(sizeof(float) * partition * This->ClusterCount);
    memset(freq, 0, sizeof(float) * partition * This->ClusterCount);
    currentWindow = 0;
    fPos = 0;
    iPos = 0;
    PointIndex = 0;
    while(PointIndex < PointCount) 
    {
        currentLength = 0;
    	for (Node=This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
    	{
            for (PointIndex=iPos; PointIndex<fPos+WindowLength&&PointIndex<PointCount; PointIndex++)
    	    {
	            if (ClusterIsMember(Node->Cluster,PointIndex))
        	        freq[ClusterIndex*partition+currentWindow] += 1;
        	}
            freq[ClusterIndex*partition+currentWindow] /= (PointIndex - iPos);
		}
        iPos = PointIndex;
        fPos += WindowLength;
        currentWindow++;
    }
    
    fprintf(OutFile, "############################################################################################\n");
    fprintf(OutFile, "#                                     Condensed Map                                         \n");
    fprintf(OutFile, "#                                                                                                   \n");
    if (WindowLength == (int)WindowLength)
    {	fprintf(OutFile, "#    %d points divided into %d windows, each window contains %d points.                 \n", PointCount, partition, (int)WindowLength);
    }
    else
    {
    	fprintf(OutFile, "#    %d points divided into %d windows, each window contains %d or %d points.          \n", PointCount, partition, (int)WindowLength,(int)WindowLength+1);
    }
    fprintf(OutFile, "#              occurence == 0.0                                                            \n");
    fprintf(OutFile, "#       0.0 <= occurence <  0.1     .                                                      \n");
    fprintf(OutFile, "#       0.1 <= occurence <  0.2     1                                                      \n");
    fprintf(OutFile, "#       0.2 <= occurence <  0.3     2                                                      \n");
    fprintf(OutFile, "#       0.3 <= occurence <  0.4     3                                                      \n");
    fprintf(OutFile, "#       0.4 <= occurence <  0.5     4                                                      \n");
    fprintf(OutFile, "#       0.5 <= occurence <  0.6     5                                                      \n");
    fprintf(OutFile, "#       0.6 <= occurence <  0.7     6                                                      \n");
    fprintf(OutFile, "#       0.7 <= occurence <  0.8     7                                                      \n");
    fprintf(OutFile, "#       0.8 <= occurence <  0.9     8                                                      \n");
    fprintf(OutFile, "#       0.9 <= occurence <  1.0     9                                                      \n");
    fprintf(OutFile, "#       1.0 == occurence            X                                                      \n");
    fprintf(OutFile, "#                                                                                          \n");
    fprintf(OutFile,"#Clustering: divide %d points into %d clusters\n", PointCount, This->ClusterCount);
    for (Node=This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
    {
    	PointIndex = ClusterCountPointsInCluster(Node->Cluster);
	    fprintf(OutFile,"#Cluster %4d: has %5d points, occurence %.3f\n",ClusterIndex, PointIndex, PointIndex/(PointCount+0.0) );
    }
  	for (Node=This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
   	{
        fprintf(OutFile, "#Cluster %4d   ",ClusterIndex); 
        for (i=0; i<partition; i++)
   	    {
            ci = (int)(freq[partition*ClusterIndex+i] * 10);
            c = ci + '0';
            if (ci == 10) c = 'X';
            if (ci == 0) 
            {
	            if (freq[partition*ClusterIndex+i]>0) 
                	c = '.';
            	else
                	c = ' ';
            }
           /* if (ClusterIsMember(Node->Cluster,PointIndex))*/
       	        fprintf(OutFile, "%c", c);
       	}
        fprintf(OutFile, "\n");
	}
    	
}

  /*  */
   void
ClusteringOutputClusterConsensusToFile(Clustering* This, FILE* OutFile,int IncludeSeparators)
{
    ClusterNode* Node;
    int PointCount;
    int PointIndex, ClusterIndex;
	char symbol[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*";
    char* consensus;
    char c;
    int MaxCluster = strlen(symbol); 
    int LineLength = 50;
    int currentLine,currentStart;
    /**/
    
    if (IncludeSeparators)
    {
        fprintf(OutFile,"##################################################################################\n");
    }
    PointCount = This->PointCount;
    fprintf(OutFile,"#Clustering: divide %d points into %d clusters\n", PointCount, This->ClusterCount);
    for (Node=This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
    {
    	PointIndex = ClusterCountPointsInCluster(Node->Cluster);
	    fprintf(OutFile,"#Cluster %4d: has %5d points, occurence %.3f\n",
		    ClusterIndex, PointIndex, PointIndex/(PointCount+0.0) );
    }
    currentLine = 0;
    if (This->ClusterCount >= MaxCluster) 
    {
        fprintf(stdout, "Too many clusters, clusters beyond %d will be represented using *.\n", MaxCluster-1);
    }
    consensus = (char*)safe_malloc(sizeof(char)*LineLength);
    while(currentLine < (int)(PointCount-1.0)/LineLength) 
    {
        currentStart = currentLine*LineLength;
    	for (Node=This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
    	{
            c = (ClusterIndex < MaxCluster)? symbol[ClusterIndex] : symbol[MaxCluster-1];
            fprintf(OutFile, "#Cluster %4d (%c)  %5d   ",ClusterIndex, c, currentStart + 1); 
            for (PointIndex=currentStart; PointIndex<currentStart+LineLength; PointIndex++)
    	    {
	            if (ClusterIsMember(Node->Cluster,PointIndex))
        	    {
                    fprintf(OutFile,"X");
                    consensus[PointIndex - currentStart] = c;
            	}
                else
                	fprintf(OutFile,".");
        	}
	        fprintf(OutFile, "\n"); 
		}
        currentLine++;
        fprintf(OutFile, "#Consensus         %5d   ",currentStart + 1); 
        for (PointIndex=currentStart; PointIndex<currentStart+LineLength; PointIndex++)
	    {
	        fprintf(OutFile, "%c", consensus[PointIndex-currentStart]);
        }
        fprintf(OutFile, "\n\n");
    }
    currentStart = currentLine*LineLength;
   	for (Node=This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
   	{
        c = (ClusterIndex < MaxCluster)? symbol[ClusterIndex] : symbol[MaxCluster-1];
        fprintf(OutFile, "#Cluster %4d (%c)  %5d   ",ClusterIndex, c, currentStart + 1); 
        for (PointIndex=currentStart; PointIndex<PointCount; PointIndex++)
   	    {
            if (ClusterIsMember(Node->Cluster,PointIndex))
       	    {
                fprintf(OutFile,"X");
           		consensus[PointIndex - currentStart] = c;
            }
            else
               	fprintf(OutFile,".");
       	}
        fprintf(OutFile, "\n"); 
           currentLine++;
	}
    fprintf(OutFile, "#Consensus         %5d   ",currentStart + 1); 
    for (PointIndex=currentStart; PointIndex<PointCount; PointIndex++)
    {
	    fprintf(OutFile, "%c", consensus[PointIndex-currentStart]);
    }
    fprintf(OutFile, "\n\n");
}

  /*
   *  For use in the paper: Produce a tab-delimited file with a snapshot of the clustering
   *  as of a certain cycle.  (All our clustering algorithms are cyclical in one way or another!)  
   *  Not every one will be interesting output; the most interesting for hierarchical are the first
   *  two, and the most interesting for linkage are the last two.  (Linkage: May want to highlight 
   *  the two "closest" points that trigger a merge, although it should be basically obvious where 
   *  the merge is happening.  )
   *  (This function is somewhat redundant from the centripetal output function, but they're both
   *  intended to be one-off code for use only in the paper, so I don't care much)
   */
   void
PAPERGenerateClusterSnapshot(Clustering* This, char* AlgName, int CycleIndex)
{
    float** Coords;
    int PointIndex;
    ClusterNode* Node;
    Cluster* TheCluster;
    int NodeIndex;
    int FreePointIndex = 0;
    int* WroteRep;
    int ClusterColumn = 2;
    int RowIndex;
    int ColumnIndex;
    int ColumnCount;
    char FileName[256];
    FILE* TheFile;
    int* CurrentRow;
    /**/
    ColumnCount = This->ClusterCount * 2; 
    Coords = (float**)malloc(sizeof(float*) * ColumnCount);
    CurrentRow = (int*)malloc(sizeof(int) * ColumnCount);
    memset(CurrentRow, 0, sizeof(int) * ColumnCount);

    /* Allocate and initialize the big array of coordinates */		    
    for (ColumnIndex=0; ColumnIndex<ColumnCount; ColumnIndex++)
    {
      Coords[ColumnIndex] = (float*)malloc(sizeof(float) * This->PointCount);
      memset(Coords[ColumnIndex], 0, sizeof(float) * This->PointCount);
    }

    /* Iterate over all the points.  Add their coordinates to the columns corresponding to the current cluster.*/
    for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
    {
      for (Node = This->Head,NodeIndex=0; Node; Node=Node->pNext,NodeIndex++)
      {
	TheCluster = (Cluster*)Node->Cluster;
	if (!ClusterIsMember(TheCluster, PointIndex))
	{
	  continue;
	}
	Coords[NodeIndex*2][CurrentRow[NodeIndex]] = ClusteringGetAttributeValue(This, PointIndex, 0);
	Coords[NodeIndex*2 + 1][CurrentRow[NodeIndex]] = ClusteringGetAttributeValue(This, PointIndex, 1);
	CurrentRow[NodeIndex] += 1;
      }
    }

    /* Output the array of coordinates */   
    sprintf(FileName,"%s%d.txt", AlgName, CycleIndex);
    TheFile = fopen(FileName,"w");
    for (RowIndex=0; RowIndex < This->PointCount; RowIndex++)
    {
      for (ColumnIndex=0; ColumnIndex < ColumnCount; ColumnIndex++)
	{
	  fprintf(TheFile,"%.4f\t",Coords[ColumnIndex][RowIndex]);
	}
      fprintf(TheFile,"\n");
    }

    free(Coords);
    free(CurrentRow);
}

  /*
   *  Compute the distance from the specified point to the centroid of the specified
   *  cluster.  Must override this, if you construct a centroid! 
   */
   float
BASEClusteringDistanceToCentroid(Clustering* This,int PointIndex, Cluster* Cluster)
{
    return GetSymElement(This->PairwiseDistances,PointIndex,ClusterGetCentroidIndex(Cluster));
}

   CLUSTER_API float
ClusteringMinDistanceToCluster(Clustering* This,int PointIndex, Cluster* Cluster)
{
	int PointCount = This->PointCount;
    int Point;
    float Distance = -1;
    float C2CDistance;
    
    for (Point = 0; Point < PointCount; Point++)
    {
    	if (!ClusterIsMember(Cluster, Point))
        {
        	continue;
        }
        C2CDistance = GetSymElement(This->PairwiseDistances,Point,PointIndex);
        if (Distance < 0 || Distance > C2CDistance)
        {
        	Distance = C2CDistance;
        }
    }
    return Distance;
}

   CLUSTER_API float
ClusteringMaxDistanceToCluster(Clustering* This,int PointIndex, Cluster* Cluster)
{
	int PointCount = This->PointCount;
    int Point;
    float Distance = -1;
    float C2CDistance;
    
    for (Point = 0; Point < PointCount; Point++)
    {
    	if (!ClusterIsMember(Cluster, Point))
        {
        	continue;
        }
        C2CDistance = GetSymElement(This->PairwiseDistances,Point,PointIndex);
        if (Distance < 0 || Distance < C2CDistance)
        {
        	Distance = C2CDistance;
        }
    }
    return Distance;
}

#define ASSIGN_ITERATION_COUNT 3

  /*
   *  Perform K-means clustering 
   */
   void
BASEClusteringClusterMeans(Clustering* This,int* SeedPoints, int DesiredClusterCount, int iteration, int mode)
{
    int PointIndex;
    int PointCount, UnprocessedPointCount;
    int* FinishedPoints;
    Cluster* NewCluster;
    float Distance; 
    float Closest;
    Cluster* ClosestCluster;
    ClusterNode *Node, *OldCluster;
    int IterationIndex, ProcessIndex;
    int NoChange;
    int ClosestClusterIndex, OldClusterIndex;
    
    /**/
    
    PointCount = This->PointCount;
    FinishedPoints = (int*)malloc(sizeof(int)*PointCount);
    memset(FinishedPoints,0,sizeof(int)*PointCount);
    /*  randomly picking */
    srand(time(NULL) + clock());
    /* Add the "seed-clusters" */
    for (PointIndex=0; PointIndex<DesiredClusterCount; PointIndex++)
    {
    	if (!SeedPoints) break; /* Use Decoy structure. */
        NewCluster = ClusteringGetNewCluster(This);
        ClusterAddMember(NewCluster,SeedPoints[PointIndex]);
        ClusterFindCentroid(NewCluster);
        ClusteringAddCluster(This,NewCluster);
        FinishedPoints[SeedPoints[PointIndex]]=1; 
    }
    UnprocessedPointCount = PointCount - DesiredClusterCount;

    /* Assign points in 3 passes.  If a point looked like it belonged to cluster A
       at first, but then we added many other points and altered our cluster shapes, it's
       possible that we'll want to reassign it to cluster B. */
    for (IterationIndex = 0; IterationIndex < iteration; IterationIndex++)
    {
      /* Add each point to an existing cluster, and re-compute centroid*/

   	  fprintf(stdout, "Round %d\n", IterationIndex);
      if (IterationIndex!=0) {
          memset(FinishedPoints,0,sizeof(int) * PointCount);
          UnprocessedPointCount = PointCount;
      }
      NoChange = 1;
      for (ProcessIndex=0; ProcessIndex<PointCount; ProcessIndex++)
      {
          if (mode == 1)
            PointIndex = ProcessIndex;
          else if (mode == 2)
            PointIndex = ChooseNextPoint(FinishedPoints,This->PointCount,UnprocessedPointCount-ProcessIndex);
          
          if (PointIndex < 0) continue;
          if (IterationIndex==0 && FinishedPoints[PointIndex])
          {
              continue;
          }
          if (IterationIndex > 0)
          {
            /* Yank this point out of its cluster, recompute the centroid */
            for (Node = This->Head; Node; Node = Node->pNext)
            {
              if (ClusterIsMember(Node->Cluster,PointIndex))
              {
  	        	/* If this point is alone in its cluster, it's in the right place already! */
	        	if (ClusterCountPointsInCluster(Node->Cluster)==1)
				{
				  continue;
				}
                OldCluster = (ClusterNode *) Node->Cluster;
                ClusterRemoveMember(Node->Cluster,PointIndex);
                ClusterFindCentroid(Node->Cluster);
              }
            }
          }
          Closest = 0;
          for (Node = This->Head; Node; Node = Node->pNext)
          {
              Distance = ClusteringDistanceToCentroid(This,PointIndex,Node->Cluster);
              if (Closest==0 || Distance<Closest)
              {
                  Closest = Distance;
                  ClosestCluster = Node->Cluster;
              }
          }
          ClusterAddMember(ClosestCluster,PointIndex);
          ClusterFindCentroid(ClosestCluster);
          if ( (ClusterNode *) Node->Cluster != OldCluster) { 
              NoChange = 0;
          }    
      }
      if (NoChange == 1)
      {
           fprintf(stdout, "No change in round %d. Skip the rest of iterations.\n", IterationIndex);
           break;
      }
    }
cleanup:
    free(FinishedPoints);
}

  /*
   *  NOTE: implementation of LinkageClustering:
   *
   *  We build a list of LinkageNodes; each LinkageNode corresponds to one of the clusters
   *  in our list, and we keep the two lists in synch.  The LinkageNode objects are where
   *  we keep track of who is close to who; that way, we don't have to recompute all our
   *  distances after each merge!
   *  In the linkage node for the first cluster, we keep track of the best distance to any other cluster.
   *  In the linkage node for each cluster, we keep track of the best distance to any cluster LATER IN
   *  THE LIST.  This optimization lets us do fewer distance checks.  It still finds all the closest
   *  links (easiest to see with a circle-and-arrow diagram)
   */

  /*
   *  Helper for Linkage clustering: Find the closest cluster to NodeA (where "closest"
   *  means minimum point-to-point distance between any two cluster members) 
   */
   void
ClusteringLinkageSeekClosest(Clustering* This,ClusterNode* NodeA,LinkageNode* Node)
{
    ClusterNode* NodeB;
    float Distance;
    for (NodeB=NodeA->pNext; NodeB; NodeB = NodeB->pNext)
    {
        Distance = ClusteringCentroidToCentroidDistance(This,NodeA,NodeB);
        /* We set a linkage node's "ClosestDistance" member to the magic value -1 to
           indicate that it is uninitialized. */
        if (Node->ClosestDistance<0 || Distance<Node->ClosestDistance)
        {
            Node->ClosestDistance=Distance;
            Node->pCluster = NodeB;
        }
    }
}

   void
ClusteringEdgeLinkSeekClosest(Clustering* This,ClusterNode* NodeA,LinkageNode* Node)
{
    ClusterNode* NodeB;
    ClusterNode* BestNode = NULL;
    float Distance;
    int PointIndexA;
    int PointIndexB;
    SymmetricMatrix* PairwiseDistances;
    float BestDistance = -1; /* Special not-yet-assigned value */
    /**/
    PairwiseDistances = This->PairwiseDistances;
    for (NodeB=NodeA->pNext; NodeB; NodeB = NodeB->pNext)
    {
      /* Compare all pairs of points */
      for (PointIndexA=0; PointIndexA < This->PointCount; PointIndexA++)
      {
    	if (!ClusterIsMember(NodeA->Cluster, PointIndexA))
	    {
	      continue;
    	}
        for (PointIndexB=PointIndexA + 1; PointIndexB < This->PointCount; PointIndexB++)
	    {
 	      if (!ClusterIsMember(NodeB->Cluster, PointIndexB))
	      {
	    	continue;
	  	  }
          Distance = SymmetricMatrixElement(PairwiseDistances, PointIndexA, PointIndexB);
          if (BestDistance<0 || Distance < BestDistance)
	      {
	    	BestDistance = Distance;
	    	BestNode = NodeB;
	  	  }
		}
      }
    }
    if (!BestNode)
    {
      /* This case should only occur when NodeA is at the end of the list, so that there are
         no distances to look at.  Do some sanity checking: */
      if (NodeA->pNext)
      {
		fprintf(stderr, "Internal error in link-edge finding: BestNode == NULL\n");
      }
      Node->ClosestDistance = -1; 
      Node->pCluster = NULL;
      return;
    }
    Node->ClosestDistance = BestDistance;
    Node->pCluster = BestNode;
}


  /*
   *  Build initial clusters and LinkageNodes - one for each point 
   */
   LinkageNode*
LinkageBuildInitialClusters(Clustering* This)
{
  int PointIndex;
  Cluster* NewCluster;
  LinkageNode* Node;
  LinkageNode* PreviousNode = NULL;
  LinkageNode* LinkHead = NULL;

  for (PointIndex=0; PointIndex<This->PointCount; PointIndex++)
  {
      NewCluster = ClusteringGetNewCluster(This);
      ClusterAddMember(NewCluster,PointIndex);
      ClusterFindCentroid(NewCluster);
      ClusteringAddCluster(This,NewCluster);
      Node = (LinkageNode*)malloc(sizeof(LinkageNode));
      if (PreviousNode)
      {
          PreviousNode->pNext = Node;
      }
      else
      {
          LinkHead = Node;
      }
      Node->pPrev=PreviousNode;
      Node->pNext=NULL;
      Node->pCluster = This->Tail;
      Node->ClosestDistance=-1;
      PreviousNode=Node;
  }
  return LinkHead;
}

  /*
   *  Single-linkage clustering - start with each point in its own singleton-cluster,
   *  and merge the closest clusters until the number of clusters is low enough (or
   *  cluster sizes are at the max allowed limit).  To make the whole process efficient,
   *  we do NOT compute all the cluster-to-cluster distances at every step.  Rather,
   *  we track the current nearest-neighbor for each cluster.  When you merge clusters,
   *  every other cluster just has to decide whether the new merged cluster will 
   *  displace its old nearest-neighbor.
   */
   void
BASEClusteringClusterLinkage(Clustering* This, int DesiredClusterCount,float Epsilon)
{
    int PointIndex;
    int PointCount;
    Cluster* NewCluster;
    float ClosestLength;
    float C2CDistance;
    /* MergeNodeA and MergeNodeB are the clusters that we merge after our search */
    ClusterNode* NodeA;
    ClusterNode* MergeNodeA;
    ClusterNode* NodeB;
    ClusterNode* MergeNodeB;
    int* ClosestIndex;
    int* ClosestCookies;
    int ClusterIndex;
    int OtherClusterIndex;
    float* ClosestSoFar;
    int MergeMe;
    int IndexOfNewCluster;
    LinkageNode* LinkHead;
    LinkageNode* Node;
    LinkageNode* PreviousNode = NULL;
    LinkageNode* LinkNodeA = NULL;
    float Distance;
    int IndexA;
    int IndexB;
    int CycleIndex;
    /**/
    PointCount = This->PointCount;
    LinkHead = LinkageBuildInitialClusters(This);

    /* Find the nearest neighbor for each of the other clusters. */
    for (NodeA = This->Head, Node=LinkHead; NodeA; NodeA=NodeA->pNext,Node=Node->pNext)
    {
        ClusteringLinkageSeekClosest(This,NodeA,Node);
    }
    CycleIndex = 0;
    while (1)
    {
        CycleIndex++;
        /*PAPERGenerateClusterSnapshot(This, "Link", CycleIndex);*/
        /* Find the closest pair, for merging: */
        ClosestLength = -1;
        IndexA=0;
        for (NodeA = This->Head, Node=LinkHead; NodeA->pNext; NodeA=NodeA->pNext,Node=Node->pNext)
        {
            if (ClosestLength < 0 || Node->ClosestDistance < ClosestLength)
            {
                MergeNodeA = NodeA;
                MergeNodeB = Node->pCluster;
                LinkNodeA=Node;
                ClosestLength = Node->ClosestDistance;
            }
            IndexA++;
        }
        if (ClosestLength<0) /* We found nothing to merge, in this case */
        {
            return;
        }
        /* If the pair are too far apart, stop! */
        if (!DesiredClusterCount && ClosestLength>Epsilon)
        {
            return;
        }

        /* Merge!  Remove cluster A, and corresponding linkage node A */
        ClusteringMergeClusters(This,MergeNodeA,MergeNodeB);
        ClusterFindCentroid(MergeNodeB->Cluster);
        if (This->ClusterCount<=DesiredClusterCount)
        {
            return; /* We've merged enough, we're done */
        }
        fprintf(stdout,"%d\n",This->ClusterCount);
        if (LinkNodeA->pNext)
        {
            LinkNodeA->pNext->pPrev = LinkNodeA->pPrev;
        }
        if (LinkNodeA->pPrev)
        {
            LinkNodeA->pPrev->pNext = LinkNodeA->pNext;
        }
        if (LinkHead==LinkNodeA)
        {
            LinkHead = LinkNodeA->pNext;
        }
        free(LinkNodeA);
        /* Fix up distances for all the nodes */
        for (NodeA = This->Head, Node=LinkHead; NodeA; NodeA=NodeA->pNext,Node=Node->pNext)
        {
            if (NodeA==MergeNodeB) 
            {
                /* This is the cluster we just merged.  Re-do its distances! */
                Node->ClosestDistance=-1;
                ClusteringLinkageSeekClosest(This,NodeA,Node);
            }
            else if (Node->pCluster == MergeNodeA || Node->pCluster == MergeNodeB)
            {
                /* This point used to have cluster A or B as a nearest neighbor.  Cluster A is
                   gone, and B is changed ,so recompute the nearest neighbor! */
                Node->ClosestDistance=-1;
                ClusteringLinkageSeekClosest(This,NodeA,Node);
            }
            else
            {
                /* Check just the distance to the merged node */
                Distance = ClusteringCentroidToCentroidDistance(This,NodeA,MergeNodeB);
                if (Distance<Node->ClosestDistance)
                {
                    Node->ClosestDistance=Distance;
                    Node->pCluster = MergeNodeB;
                }
            }
        }
    }
}


  /*
   *  Edge-joining clustering: Join the clusters with shortest point-to-point (not
   *  centroid-to-centroid) distance.
   */
   void
BASEClusteringClusterEdgeLink(Clustering* This, int DesiredClusterCount,float Epsilon)
{
    int PointIndex;
    int PointCount;
    Cluster* NewCluster;
    float ClosestLength;
    float C2CDistance;
    /* MergeNodeA and MergeNodeB are the clusters that we merge after our search */
    ClusterNode* NodeA;
    ClusterNode* MergeNodeA;
    ClusterNode* NodeB;
    ClusterNode* MergeNodeB;
    int* ClosestIndex;
    int* ClosestCookies;
    int ClusterIndex;
    int OtherClusterIndex;
    float* ClosestSoFar;
    int MergeMe;
    int IndexOfNewCluster;
    LinkageNode* LinkHead;
    LinkageNode* Node;
    LinkageNode* PreviousNode = NULL;
    LinkageNode* LinkNodeA = NULL;
    float Distance;
    int IndexA;
    int IndexB;
    int CycleIndex;
    /**/
    PointCount = This->PointCount;
    LinkHead = LinkageBuildInitialClusters(This);

    /* Find the nearest neighbor for each of the other clusters. */
    for (NodeA = This->Head, Node=LinkHead; NodeA; NodeA=NodeA->pNext,Node=Node->pNext)
    {
        ClusteringEdgeLinkSeekClosest(This,NodeA,Node);
    }
    CycleIndex = 0;
    while (1)
    {
        CycleIndex++;
        /*PAPERGenerateClusterSnapshot(This, "EdgeLink", CycleIndex);*/
		ClusteringOutputStatus(This);
        /* Find the closest pair, for merging: */
        ClosestLength = -1;
        IndexA=0;
        for (NodeA = This->Head, Node=LinkHead; NodeA->pNext; NodeA=NodeA->pNext,Node=Node->pNext)
        {
            if (ClosestLength < 0 || Node->ClosestDistance < ClosestLength)
            {
                MergeNodeA = NodeA;
                MergeNodeB = Node->pCluster;
                LinkNodeA=Node;
                ClosestLength = Node->ClosestDistance;
            }
            IndexA++;
        }
        if (ClosestLength<0) /* We found nothing to merge, in this case */
        {
            return;
        }
        /* If the pair are too far apart, stop! */
        if (!DesiredClusterCount && ClosestLength>Epsilon)
        {
            return;
        }

        /* Merge!  Remove cluster A, and corresponding linkage node A */
        ClusteringMergeClusters(This,MergeNodeA,MergeNodeB);
        ClusterFindCentroid(MergeNodeB->Cluster);
        if (This->ClusterCount<=DesiredClusterCount)
        {
            return; /* We've merged enough, we're done */
        }
        fprintf(stdout,"%d\n",This->ClusterCount);
        if (LinkNodeA->pNext)
        {
            LinkNodeA->pNext->pPrev = LinkNodeA->pPrev;
        }
        if (LinkNodeA->pPrev)
        {
            LinkNodeA->pPrev->pNext = LinkNodeA->pNext;
        }
        if (LinkHead==LinkNodeA)
        {
            LinkHead = LinkNodeA->pNext;
        }
        free(LinkNodeA);
        /* Fix up distances for all the nodes */
        for (NodeA = This->Head, Node=LinkHead; NodeA; NodeA=NodeA->pNext,Node=Node->pNext)
        {
            if (NodeA==MergeNodeB) 
            {
  	        /* This is the cluster we just merged.  Re-find its closest neighbors! */
                Node->ClosestDistance=-1;
                ClusteringEdgeLinkSeekClosest(This,NodeA,Node);
            }
            else if (Node->pCluster == MergeNodeA || Node->pCluster == MergeNodeB)
            {
                /* This point used to have cluster A or B as a nearest neighbor.  The merged
		   cluster is still its nearest neighbor, and the distance is unchanged!  */
                Node->pCluster = MergeNodeB; 
            }
	    /* If the cluster didn't have one of those two clusters as its nearest neighbor, then
	       it WILL NOT have the merged cluster as its nearest neighbor either.  (You can be a nearest
	       neighbor only if you have a pair of points that's very close, and the merged cluster has no new points!)
	    */
        }
    }
}


  /*
   *   Merge two clusters into one large cluster containing all the points from the 
   *     original.  (We remove NodeA from our list, and put the new cluster into NodeB) 
   */
   void
BASEClusteringMergeClusters(Clustering* This,ClusterNode* NodeA,ClusterNode* NodeB)
{
    int PointIndex;
    /**/
    for (PointIndex=0; PointIndex<This->PointCount; PointIndex++)
    {
        if (ClusterIsMember(NodeA->Cluster,PointIndex))
        {
            ClusterAddMember(NodeB->Cluster,PointIndex);
        }
    }
    ClusteringRemoveCluster(This,NodeA->Cluster);
}


  /*
   *  Average-linkage clustering 
   */
   void
ClusteringAverageLinkage(Clustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int* ClusterPointCounts; /* Save for the point count for each cluster as in the cluster cookies. */
    int PointIndex, PointIndexA, PointIndexB;
    int PointCount = This->PointCount;
    Cluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    int ClosestNodeACookie;
    int CountA, CountB;
    float DistA, DistB;
    int PointCountA, PointCountB;
    float Distance, NewDist;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    int finish = 0;

    /* */
    ClusterCookies = (int*)malloc(sizeof(int)*PointCount);
    ClusterPointCounts = (int*)malloc(sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = ClusterNew(This);
        ClusterAddMember(NewCluster,PointIndex);
        ClusteringAddCluster(This,NewCluster);
        ClusterCookies[PointIndex] = PointIndex;
        ClusterPointCounts[PointIndex] = 1;
    }

    /* Keep joining until you're finished */
    while (!finish)
    {
        /* If there's 1 cluster left, then we're definitely done: */
		if (This->ClusterCount < 2)
        {
	  		finish = 1;
            continue;
		}
        MergeNodeB = This->Head;
        MergeNodeA = MergeNodeB->pNext;
        ClosestNodeIndexB = 0; 
        ClosestNodeIndexA = 1; 
        ClosestLength = GetSymElement(ClusterDistances,ClusterCookies[ClosestNodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
        {
            for (NodeA = NodeB->pNext,NodeIndexA=NodeIndexB+1; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
            {
                C2CDistance = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance < ClosestLength)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
        if (!DesiredClusterCount && ClosestLength > Epsilon)
        {
            finish = 1;
            continue;
        }
		/* If everything has been merged together into one big mess, return. */
		if (!MergeNodeA)
        {
	  		finish = 1;
            continue;
        }
		/* Merge A and B into node B */
        ClusteringMergeClusters(This,MergeNodeA,MergeNodeB);
        
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update cookies */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA]; /* Remember the old cookie for Node A. */
        for (NodeIndexA=0; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            if (NodeIndexA>=ClosestNodeIndexA)
            {
                ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
            }
        }
        if (ClosestNodeIndexB>ClosestNodeIndexA)
        {
            ClosestNodeIndexB--;
        }
        /* Update cluster-to-cluster distances */
        CountA = ClusterPointCounts[ClosestNodeACookie];
        CountB = ClusterPointCounts[ClusterCookies[ClosestNodeIndexB]];
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            DistA = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
            NewDist = (DistA * CountA + DistB * CountB) / (CountA + CountB);
	        SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],NewDist);
        }
        ClusterPointCounts[ClusterCookies[ClosestNodeIndexB]] = CountA + CountB;
    }
    safe_free(ClusterDistances);
    safe_free(ClusterCookies);
    safe_free(ClusterPointCounts);
}

  /*
   *  Complete-linkage clustering 
   */
   void
ClusteringCompleteLinkage(Clustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex, PointIndexA, PointIndexB;
    int PointCount = This->PointCount;
    Cluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    int ClosestNodeACookie;
    int PointCountA, PointCountB;
    float Distance, NewDist;
    float DistA, DistB;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    int finish = 0;

    /* */
    ClusterCookies = (int*)malloc(sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = ClusterNew(This);
        ClusterAddMember(NewCluster,PointIndex);
        ClusteringAddCluster(This,NewCluster);
        ClusterCookies[PointIndex] = PointIndex;
    }

    /* Keep joining until you're finished */
    while (!finish)
    {
        /* If there's 1 cluster left, then we're definitely done: */
		if (This->ClusterCount < 2)
        {
	  		finish = 1;
            continue;
		}
        MergeNodeB = This->Head;
        MergeNodeA = MergeNodeB->pNext;
        ClosestNodeIndexB = 0; 
        ClosestNodeIndexA = 1; 
        ClosestLength = GetSymElement(ClusterDistances,ClusterCookies[ClosestNodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
        {
            for (NodeA = NodeB->pNext,NodeIndexA=NodeIndexB+1; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
            {
                C2CDistance = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance < ClosestLength)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
        if (!DesiredClusterCount && ClosestLength > Epsilon)
        {
            finish = 1;
            continue;
        }
		/* If everything has been merged together into one big mess, return. */
		if (!MergeNodeA)
        {
	  		finish = 1;
            continue;
        }
		/* Merge A and B into node B */
        ClusteringMergeClusters(This,MergeNodeA,MergeNodeB);
        
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update cookies */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA]; /* Remember the old cookie for Node A. */
        for (NodeIndexA=0; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            if (NodeIndexA>=ClosestNodeIndexA)
            {
                ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
            }
        }
        if (ClosestNodeIndexB>ClosestNodeIndexA)
        {
            ClosestNodeIndexB--;
        }
        /* Update cluster-to-cluster distances */
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            DistA = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
            if (DistA < DistB)
	            SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],DistA);
        }
    }
    safe_free(ClusterDistances);
    safe_free(ClusterCookies);
}

  /*
   *  Edge-linkage clustering 
   */
   void
ClusteringEdgeLinkage(Clustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex, PointIndexA, PointIndexB;
    int PointCount = This->PointCount;
    Cluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    int ClosestNodeACookie;
    int PointCountA, PointCountB;
    float Distance, NewDist;
    float DistA, DistB;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    int finish = 0;

    /* */
    ClusterCookies = (int*)malloc(sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = ClusterNew(This);
        ClusterAddMember(NewCluster,PointIndex);
        ClusteringAddCluster(This,NewCluster);
        ClusterCookies[PointIndex] = PointIndex;
    }

    /* Keep joining until you're finished */
    while (!finish)
    {
        /* If there's 1 cluster left, then we're definitely done: */
		if (This->ClusterCount < 2)
        {
	  		finish = 1;
            continue;
		}
        MergeNodeB = This->Head;
        MergeNodeA = MergeNodeB->pNext;
        ClosestNodeIndexB = 0; 
        ClosestNodeIndexA = 1; 
        ClosestLength = GetSymElement(ClusterDistances,ClusterCookies[ClosestNodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
        {
            for (NodeA = NodeB->pNext,NodeIndexA=NodeIndexB+1; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
            {
                C2CDistance = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance < ClosestLength)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
        if (!DesiredClusterCount && ClosestLength > Epsilon)
        {
            finish = 1;
            continue;
        }
		/* If everything has been merged together into one big mess, return. */
		if (!MergeNodeA)
        {
	  		finish = 1;
            continue;
        }
		/* Merge A and B into node B */
        ClusteringMergeClusters(This,MergeNodeA,MergeNodeB);
        
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update cookies */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA]; /* Remember the old cookie for Node A. */
        for (NodeIndexA=0; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            if (NodeIndexA>=ClosestNodeIndexA)
            {
                ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
            }
        }
        if (ClosestNodeIndexB>ClosestNodeIndexA)
        {
            ClosestNodeIndexB--;
        }
        /* Update cluster-to-cluster distances */
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            DistA = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
            if (DistA < DistB)
	            SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],DistA);
        }
    }
    safe_free(ClusterDistances);
    safe_free(ClusterCookies);
}




  /*
   *  Perform HIERARCHICAL clustering 
   */
   void
ClusteringClusterHierarchical(Clustering* This,int DesiredClusterCount,float Epsilon)
{
    float EccentricEpsilon;
    Cluster* EccentricCluster;
    int EccentricPointA;
    int EccentricPointB;
    Cluster* NewCluster;
    int CycleIndex;
    /**/
    NewCluster = ClusteringGetNewCluster(This);
    ClusterAddEverything(NewCluster);
    ClusteringAddCluster(This,NewCluster);
    CycleIndex = 0;
    while (1)
    {
        CycleIndex++;
        /*PAPERGenerateClusterSnapshot(This, "Hier", CycleIndex);*/
	ClusteringOutputStatus(This);

        ClusteringFindEccentricCluster(This,&EccentricEpsilon,&EccentricCluster,
            &EccentricPointA,&EccentricPointB);
        if (!DesiredClusterCount && EccentricEpsilon<Epsilon)
        {
            break;
        }
        ClusteringSplitCluster(This,EccentricCluster,EccentricPointA,EccentricPointB);
        if (DesiredClusterCount>0 && This->ClusterCount>=DesiredClusterCount)
        {
            return;
        }
        if (DesiredClusterCount<0 && This->ClusterCount >= -DesiredClusterCount)
        {
            return;
        }
    }
}

  /*
   *  Split a cluster.  Each of the two specified points becomes a "seed" for one of
   *  the two children. Each point in ParentCluster is assigned to the child-cluster
   *  whose seed it is nearest to.  
   */
   void
BASEClusteringSplitCluster(Clustering* This,Cluster* ParentCluster, int PointA, int PointB)
{
   /* If you're close to PointA, you'll live in ParentCluster.  If you're closer to point B,
      you'll live in NewCluster. */
   Cluster* NewCluster;
   int PointCount;
   int PointIndex;
   int ReassignCount = 0;
   int CycleIndex;
   float DistA;
   float DistB;
   /**/
   PointCount = This->PointCount;
   NewCluster = ClusteringGetNewCluster(This);
   for (PointIndex=0; PointIndex<PointCount; PointIndex++)
   {
       if (!ClusterIsMember(ParentCluster,PointIndex))
       {
           continue;
       }
       DistA = GetSymElement(This->PairwiseDistances,PointIndex,PointA);
       DistB = GetSymElement(This->PairwiseDistances,PointIndex,PointB);
       if (DistB<DistA)
       {
           ClusterRemoveMember(ParentCluster,PointIndex);
           ClusterAddMember(NewCluster,PointIndex);
       }
   }
   ClusterFindCentroid(ParentCluster);
   ClusterFindCentroid(NewCluster);
   /* Improve the assignments by re-assigning any points that are closer to the other 
      cluster's centroid */
   for (CycleIndex=0; CycleIndex<3; CycleIndex++)
     {
       fprintf(stdout,"Reassignment cycle %d:\n",CycleIndex);
       for (PointIndex=0; PointIndex<PointCount; PointIndex++)
       {
    	   if (ClusterIsMember(ParentCluster,PointIndex) || ClusterIsMember(NewCluster,PointIndex))
    	   {
    	     DistA = ClusteringDistanceToCentroid(This,PointIndex,ParentCluster);
    	     DistB = ClusteringDistanceToCentroid(This,PointIndex,NewCluster);
    	     if (DistA < DistB && ClusterIsMember(NewCluster,PointIndex))
    	       {
    		 ClusterRemoveMember(NewCluster,PointIndex);
    		 ClusterAddMember(ParentCluster,PointIndex);
		 fprintf(stdout,"Reassign %d",PointIndex);
    		 ReassignCount++;
    	       }
    	     if (DistB < DistA && ClusterIsMember(ParentCluster,PointIndex))
    	       {
    		 ClusterRemoveMember(ParentCluster,PointIndex);
    		 ClusterAddMember(NewCluster,PointIndex);
		 fprintf(stdout,"Reassign %d",PointIndex);
    		 ReassignCount++;
    	       }
    	   }
       }   
       ClusterFindCentroid(ParentCluster);
       ClusterFindCentroid(NewCluster);

     }
   fprintf(stdout,"Reassignment count:%d\n",ReassignCount);
   ClusteringAddCluster(This,NewCluster);
}


  /*
   *  Find the most eccentric cluster.  The eccentricity of a cluster is defined as
   *  the maximal point-to-point distance for points in the cluster.  This function 
   *  also tracks the distance (Epsilon), and the indices of the two eccentric 
   *  points.
   */
   void
ClusteringFindEccentricCluster(Clustering* This, float* Epsilon, Cluster** EccentricCluster, 
			       int* PointA, int* PointB)
{
    float CurrentEccentricity;
    ClusterNode* Node;
    int PointIndexA;
    int PointIndexB;
    int PointCount;
    /**/
    *Epsilon = 0;
    PointCount = This->PointCount;
    for (Node = This->Head; Node; Node=Node->pNext)
    {
        /* Iterate over all pairs of points in this cluster: */
        for (PointIndexA=0; PointIndexA<PointCount; PointIndexA++)
        {
            if (!ClusterIsMember(Node->Cluster,PointIndexA))
            {
                continue;
            }
            for (PointIndexB=PointIndexA+1; PointIndexB<PointCount; PointIndexB++)
            {
                if (!ClusterIsMember(Node->Cluster,PointIndexB))
                {
                    continue;
                }
                /* Is this pair farther apart than we've ever seen? */
                CurrentEccentricity = SymmetricMatrixElement(This->PairwiseDistances,PointIndexA,PointIndexB);
                if (CurrentEccentricity > *Epsilon)
                {
                    *Epsilon = CurrentEccentricity;
                    *EccentricCluster = Node->Cluster;
                    *PointA = PointIndexA;
                    *PointB = PointIndexB;
                }
            }
        }
    }
}

  /*
   *  Return the head of our cluster-list 
   */
   ClusterNode*
ClusteringGetHead(Clustering* This)
{
    return This->Head;
}

  /*
   *  Find some seed-points for K-means clustering.  Take the first point as an 
   *  arbitrary first choice.  Then, at each iteration, add the point whose minimal
   *  distance from any seed is as large as possible.  
   */
   int*
BASEClusteringFindKmeansSeeds(Clustering* This, int Seeds)
{
    int* SeedIndices;
    int SeedIndex;
    float BestDistance;
    int CandidatePointIndex;
    int SkipFlag;
    int CheckIndex;
    float CandidateNearestDistance;
    float Distance;
    int BestDistanceIndex;
    /**/
    SeedIndices = (int*)malloc(sizeof(int) * Seeds);
    SeedIndices[0]=1; /* Arbitrary first-choice */
    for (SeedIndex=1; SeedIndex<Seeds; SeedIndex++)
    {
        BestDistance = 0;
        for (CandidatePointIndex=0; CandidatePointIndex<This->PointCount; CandidatePointIndex++)
        {
            /* Make sure this candidate isn't already a seed */
            SkipFlag=0;
            for (CheckIndex=0; CheckIndex<SeedIndex; CheckIndex++)
            {
                if (SeedIndices[CheckIndex]==CandidatePointIndex)
                {
                    SkipFlag=1;
                }
            }
            if (SkipFlag)
            {
                continue;
            }
            /* Get the closest distance from this candidate to a current seed */
            CandidateNearestDistance = 0;
            for (CheckIndex=0; CheckIndex<SeedIndex; CheckIndex++)
            {
                Distance = GetSymElement(This->PairwiseDistances,SeedIndices[CheckIndex],CandidatePointIndex);
                if (Distance<CandidateNearestDistance || CandidateNearestDistance==0)
                {
                    CandidateNearestDistance = Distance;
                }
            }
            /* Is this the best so far? */
            if (CandidateNearestDistance > BestDistance)
            {
                BestDistance = CandidateNearestDistance;
                BestDistanceIndex = CandidatePointIndex;
            }
        }
        SeedIndices[SeedIndex] = BestDistanceIndex;
    }
    return SeedIndices;
}

  /*
   *  Return the distance between two clusters' points 
   */
   float
BASEClusteringPointToPointDistance(Clustering* This, int A, int B)
{	
	if (This->PairwiseDistances)
	    return GetSymElement(This->PairwiseDistances,A,B);
    else
    	fprintf(stderr, "Error in line %d File %s, no PairwiseDistance matrix.\n", __LINE__, __FILE__);

    return 0.0;
}

  /*
   *  Return the distance between two clusters' centroids 
   */
   float
BASEClusteringCentroidToCentroidDistance(Clustering* This, ClusterNode* NodeA,ClusterNode* NodeB)
{
    return GetSymElement(This->PairwiseDistances,ClusterGetCentroidIndex(NodeA->Cluster),ClusterGetCentroidIndex(NodeB->Cluster));
}

  /*
   *  Free memory from a cobweb cluster.  Frees JUST this node, not its children. 
   */
   void
BASEClusteringCobwebClusterFree(Clustering* This, CobwebCluster* DeadGuy) 
{
    safe_free(DeadGuy);
}

  /*
   *  Helper function for cobweb clustering: Merge two nodes.  (Creates a new parent, with these 
   *  two nodes as children)
   */
   CobwebTreeNode*
CobwebClusterMergeNodes(Clustering* This, CobwebTreeNode* First,CobwebTreeNode* Second)
{
    CobwebTreeNode* Merged;

    Merged = CobwebTreeNodeNew();
    Merged->Cluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
    CobwebTreeNodeAddChild(Merged,First);
    CobwebTreeNodeAddChild(Merged,Second);
    return Merged;
}

  /*
   *  Allocate a new (non-leaf) CobwebTreeNode 
   */
   CLUSTER_API CobwebTreeNode*
CobwebTreeNodeNew()
{
    CobwebTreeNode* This;
    This = (CobwebTreeNode*)malloc(sizeof(CobwebTreeNode));
    memset(This,0,sizeof(CobwebTreeNode));
    return This;
}

  /*
   *  Orphan the specified CobewebTreeNode: Disconnect it from its parents (updating their children-count)
   */
   void
CobwebTreeNodeRemoveChild(CobwebTreeNode* Child)
{
    CobwebTreeNode* Temp;
    if (Child->Parent)
    {
        if (Child->Parent->FirstChild==Child)
        {
            Child->Parent->FirstChild = Child->Next;
        }
        if (Child->Parent->LastChild==Child)
        {
            Child->Parent->LastChild = Child->Previous;
        }
        for (Temp=Child->Parent; Temp; Temp=Temp->Parent)
        {
            Temp->TotalFrameCount -= Child->TotalFrameCount;
        }
    }
    if (Child->Next)
    {
        Child->Next->Previous = Child->Previous;
    }
    if (Child->Previous)
    {
        Child->Previous->Next = Child->Next;
    }
    Child->Next=NULL;
    Child->Previous=NULL;
    Child->Parent = NULL;
}

  /*
   *  Append an item to a CTNList.  Returns the new tail.  
   */
   CTNList*
CTNListAppend(CTNList* Tail, CobwebTreeNode* CobwebNode)
{
    CTNList* Current; 
    /**/
    Current = (CTNList*)malloc(sizeof(CTNList));
    memset(Current,0,sizeof(CTNList));
    Current->Node = CobwebNode;
    Current->Previous = Tail;
    if (Tail)
    {
        Tail->Next = Current;
    }
    return Current;
}

  /*
   *  Pass in a CobwebTreeNode that plays the role of parent.  This 
   *  method returns a CTNList containing each of the parent's child 
   *  CobwebTreeNodes as children 
   */
   CTNList*
CobwebTreeNodeListChildren(CobwebTreeNode* Node)
{
    CTNList* Head = NULL;
    CTNList* Prev = NULL;
    CTNList* Current;
    CobwebTreeNode* Child;
    for (Child=Node->FirstChild; Child; Child=Child->Next)
    {
        Current = (CTNList*)malloc(sizeof(CTNList));
        memset(Current,0,sizeof(CTNList));
        Current->Node = Child;
        Current->Previous = Prev;
        if (Prev)
        {
            Prev->Next = Current;
        }
        if (!Head)
        {
            Head=Current;
        }
        Prev=Current;
    }
    return Head;
}

  /*
   * Free a CTNList.  (Doesn't free the underlying CobwebTreeNodes) 
   */
   void
CTNListFree(CTNList* Head)
{
    CTNList* Node;
    CTNList* Prev = NULL;
    for (Node = Head; Node; Node=Node->Next)
    {
        safe_free(Prev);
        Prev=Node;
    }
    safe_free(Prev);
}

  /*
   *  Create a new Cobweb cluster object.  You should override this,
   *  to allocate the right number of means and stddevs.  
   */
   CobwebCluster*
BASEClusteringNewCobwebCluster(Clustering* ThisClustering, int PointIndex)
{
  CobwebCluster* This;
  /**/
  This = (CobwebCluster*)malloc(sizeof(CobwebCluster));
  memset(This,0,sizeof(CobwebCluster));
  This->PointIndex=PointIndex;
  return This;
}

  /*
   *  Free memory from a CobwebTreeNode and its children 
   */
   CLUSTER_API void
CobwebTreeNodeFree(CobwebTreeNode* Head)
{
    CobwebTreeNode* Node;
    CobwebTreeNode* Prev = NULL;
    if (!Head)
    {
        return;
    }
    CobwebTreeNodeRemoveChild(Head);
    for (Node = Head->FirstChild; Node; Node = Node->Next)
    {
        CobwebTreeNodeFree(Prev);
        Prev=Node;
    }
    CobwebTreeNodeFree(Prev);
    safe_free(Head);
}

  /*
   *  Free memory from a CobwebTreeNode, but DO NOT eliminate the children 
   *  (NR = Not Recursive) 
   */
   void
CobwebTreeNodeFreeNR(CobwebTreeNode* Head)
{
    CobwebTreeNode* Node;
    CobwebTreeNode* Prev = NULL;
    if (!Head)
    {
        return;
    }
    CobwebTreeNodeRemoveChild(Head);
    safe_free(Head);
}

  /*
   *  Add a child to a CobwebTreeNode. 
   */
   void
CobwebTreeNodeAddChild(CobwebTreeNode* Parent, CobwebTreeNode* Child)
{
    CobwebTreeNode* Temp;
    if (!Parent->FirstChild)
    {
        Parent->FirstChild = Child;
        Parent->LastChild = Child;
        Child->Next = NULL;
        Child->Previous = NULL;
    }
    else
    {
        Parent->LastChild->Next = Child;
        Child->Previous = Parent->LastChild;
        Child->Next = NULL;
        Parent->LastChild = Child;
    }
    Child->Parent = Parent;
    /* Proceed up the tree, keeping the framecount accurate in ancestors */
    for (Temp=Parent; Temp; Temp=Temp->Parent)
    {
        Temp->TotalFrameCount += Child->TotalFrameCount;
    }
}

  /*
   *  Allocate a leaf CobwebTreeNode corresponding to the given point
   */
   CobwebTreeNode*
CobwebTreeNodeNewLeaf(Clustering* This, int PointIndex)
{
    CobwebTreeNode* CTN;
    CTN = (CobwebTreeNode*)malloc(sizeof(CobwebTreeNode));
    memset(CTN,0,sizeof(CobwebTreeNode));
    CTN->Cluster = ClusteringNewCobwebCluster(This,PointIndex);
    CTN->TotalFrameCount = 1;
    return CTN;
}

  /*
   *  Helper function for COBWEB clustering: Compute means and standard-deviations
   *  for the specified tree node.  This function should be called on a non-leaf 
   *  node whenever there is a change in its chldren.  For permanent changes to the
   *  tree, you should *also* call this method for all of the node's parents. 
   */
   void
BASECobwebClusterComputeMeans(Clustering* This, CobwebTreeNode* Node)
{
  /* Abstract method: Does nothing! */
  return;
}

  /*
   *  Helper function for CobwebTreeNodes: Replaces one node's tree position with another.
   *  NB: We don't recurse total-frame counts, so we'd better swap back later to keep integrity!
   */
   void
CobwebTreeNodeSwap(CobwebTreeNode* Old, CobwebTreeNode* New)
{
    New->Previous = Old->Previous;
    New->Next = Old->Next;
    New->Parent = Old->Parent;
    if (Old->Previous)
    {
        Old->Previous->Next = New;
    }
    if (Old->Next)
    {
        Old->Next->Previous = New;
    }
    if (Old->Parent)
    {
        Old->Parent->TotalFrameCount -= Old->TotalFrameCount;
        Old->Parent->TotalFrameCount += New->TotalFrameCount;
        if (Old->Parent->FirstChild==Old)
        {
            Old->Parent->FirstChild=New;
        }
        if (Old->Parent->LastChild==Old)
        {
            Old->Parent->LastChild=New;
        }

    }
}

#define MAX_LEVEL 50 
#define INDENT_DEPTH 2

  /*
   *  Print a CobwebTreeNode and its descendants to stdout, recursively. 
   *  If level < 0, No limit.
   */
  void
DebugPrintCobwebTree(Clustering* This, CobwebTreeNode* Root, int Indent, int level)
{
    char IndentString[MAX_LEVEL * INDENT_DEPTH];
    int SpaceIndex;
    CobwebTreeNode* Child;
    
    /**/
    if (level == 0) 
    {
    	return;
    }
    memset(IndentString,0,sizeof(char) * MAX_LEVEL * INDENT_DEPTH);
    Indent = min(Indent,MAX_LEVEL * INDENT_DEPTH - 1);
    for (SpaceIndex=0; SpaceIndex<Indent; SpaceIndex++)
    {
        IndentString[SpaceIndex]=' ';
    }
    if (Root->Cluster->PointIndex==COBWEB_NONLEAF)
    {
        fprintf(stdout,"%sNode: holds %d frames\n",IndentString,Root->TotalFrameCount);
    }
    else
    {
        fprintf(stdout,"%sLeaf: Frame #%d\n",IndentString,Root->Cluster->PointIndex);
    }
    for (Child = Root->FirstChild; Child; Child = Child->Next)
    {
        DebugPrintCobwebTree(This,Child,Indent + INDENT_DEPTH, level - 1);
    }
}

   void
DebugWriteCobwebTree(Clustering* This, FILE* OutFile, CobwebTreeNode* Root, int Indent, int level)
{
    char IndentString[MAX_LEVEL * INDENT_DEPTH];
    int SpaceIndex;
    CobwebTreeNode* Child;
    int AttributeIndex;
    int AttributeCount = ClusteringGetAttributeCount(This);
    
    /**/
    if (level == 0)
    {
    	return;
    }
    memset(IndentString,0,sizeof(char) * MAX_LEVEL * INDENT_DEPTH);
    Indent = min(Indent, MAX_LEVEL * INDENT_DEPTH - 1);
    for (SpaceIndex=0; SpaceIndex<Indent; SpaceIndex++)
    {
        IndentString[SpaceIndex]=' ';
    }
    if (Root->Cluster->PointIndex==COBWEB_NONLEAF)
    {
        fprintf(OutFile,"%sNode: holds %d frames, ",IndentString,Root->TotalFrameCount);
        fprintf(OutFile,"CU is %.4f\n",ClusteringCobwebCategoryUtility(This,Root));
		
        fprintf(OutFile,"%sMeans (", IndentString); 
    	for (AttributeIndex = 0; AttributeIndex < min(AttributeCount-1, 6); AttributeIndex++) 
        {
	       fprintf(OutFile,"%.4f ",Root->Cluster->Means[AttributeIndex]);
	    }
        if (AttributeIndex < AttributeCount-1)
        {
	       fprintf(OutFile,"... ");
        }
        fprintf(OutFile,"%.4f)\n", Root->Cluster->Means[AttributeCount-1]);
		
        fprintf(OutFile,"%sStddev (", IndentString); 
    	for (AttributeIndex = 0; AttributeIndex < min(AttributeCount-1, 6); AttributeIndex++) 
        {
	       fprintf(OutFile,"%.4f ",Root->Cluster->Stddevs[AttributeIndex]);
	    }
        if (AttributeIndex < AttributeCount-1)
        {
	       fprintf(OutFile,"... ");
        }
        fprintf(OutFile,"%.4f)\n", Root->Cluster->Stddevs[AttributeCount-1]);
    }
    else
    {
        fprintf(OutFile,"%sLeaf: Frame #%d\n",IndentString,Root->Cluster->PointIndex);
		fprintf(OutFile,"%sAt (", IndentString); 
    	for (AttributeIndex = 0; AttributeIndex < min(AttributeCount-1, 6); AttributeIndex++) 
        {
	       fprintf(OutFile,"%.4f ",ClusteringGetAttributeValue(This, Root->Cluster->PointIndex, AttributeIndex));
	    }
        if (AttributeIndex < AttributeCount-1)
        {
	       fprintf(OutFile,"... ");
        }
        fprintf(OutFile,"%.4f)\n", ClusteringGetAttributeValue(This, Root->Cluster->PointIndex, AttributeCount-1));
    }
    for (Child = Root->FirstChild; Child; Child = Child->Next)
    {
        DebugWriteCobwebTree(This, OutFile, Child, Indent + INDENT_DEPTH, level - 1);
    }
}

   CobwebTreeNode*
ReadCobwebTree(Clustering* This, char* CobwebTreeFile)
{
    char IndentString[MAX_LEVEL * INDENT_DEPTH];
    char FileLine[BLOCK_SIZE];
    char *Result;
    int SpaceIndex;
    FILE* TreeFile;
    CobwebTreeNode *Root;
    CobwebTreeNode *PreviousNode;
    CobwebTreeNode *NewNode;
    CobwebTreeNode *TempNode;
    int ThisLevel, PreviousLevel;
    int FrameIndex, PoundIndex;
    
    /**/
    TreeFile = fopen(CobwebTreeFile, "r");
    if (!TreeFile)
    {
    	return NULL;
    }
    if (TreeFile)
    {
    	fprintf(stdout, "  Read in the existing CobwebPreCoalesce file...");
    }
	if (!CheckCobwebTreeFile(This, TreeFile))
    {
    	fprintf(stdout, "\nOverwrite the CobwebPreCoalesce file.\n");
    	fclose(TreeFile);
        return NULL;
    }
   	Result = fgets(FileLine, BLOCK_SIZE ,TreeFile);
    if (strncmp(FileLine, "Node: holds", 10))
    {
    	fprintf(stdout, "\nCobwebPreCoalesce file is not correct.\n");
    	return NULL;
    }
    ThisLevel = PreviousLevel = 0;
    Root = CobwebTreeNodeNew();
   	Root->Cluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
    /* Ignore the next two lines for means and stddev. */
   	Result = fgets(FileLine, BLOCK_SIZE ,TreeFile); 
  	Result = fgets(FileLine, BLOCK_SIZE ,TreeFile); 
    PreviousNode = Root;
    Result = fgets(FileLine, BLOCK_SIZE ,TreeFile);
    while (Result)
    {
        SpaceIndex = strcspn(FileLine, "NL");
        ThisLevel = SpaceIndex / INDENT_DEPTH;
        if (FileLine[SpaceIndex] == 'N')
        {
		    NewNode = CobwebTreeNodeNew();
        	NewNode->Cluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
            /* Ignore the next two lines for means and stddev. */
	    	Result = fgets(FileLine, BLOCK_SIZE ,TreeFile); 
	    	Result = fgets(FileLine, BLOCK_SIZE ,TreeFile); 
        } else {
		    PoundIndex = strchr(FileLine, '#') - FileLine;
            sscanf(FileLine+PoundIndex+1, "%d", &FrameIndex);
            NewNode = CobwebTreeNodeNewLeaf(This, FrameIndex);
	    	Result = fgets(FileLine, BLOCK_SIZE ,TreeFile); /* Ignore the next line for attribute values output. */
        }
        if (ThisLevel > PreviousLevel) 
        {
        	CobwebTreeNodeAddChild(PreviousNode, NewNode);
        } else if (ThisLevel == PreviousLevel)
        {
        	CobwebTreeNodeAddChild(PreviousNode->Parent, NewNode);
        } else
        {
        	while (PreviousLevel > ThisLevel)
            {
            	PreviousNode = PreviousNode->Parent;
                PreviousLevel--;
            }
        	CobwebTreeNodeAddChild(PreviousNode->Parent, NewNode);
        }
        for (TempNode = NewNode->Parent;TempNode;TempNode = TempNode->Parent)
        {
        	CobwebClusterComputeMeans(This,TempNode);
        }
        PreviousNode = NewNode;
        PreviousLevel = ThisLevel;
    	Result = fgets(FileLine, BLOCK_SIZE ,TreeFile);
    }
   	fprintf(stdout, "DONE.\n");
    fclose(TreeFile);
    return Root;
}

   int
CheckCobwebTreeFile(Clustering* This, FILE* TreeFile)
{
    char FileLine[BLOCK_SIZE];
    char *Result, *cpos;
    float f1, f2;
    int i, pos;
    int AttributeCount = ClusteringGetAttributeCount(This);
    
    Result = fgets(FileLine, BLOCK_SIZE ,TreeFile);
    cpos = strchr(FileLine, ':');
    sscanf(cpos+1, "%d", &i);
    if (i != This->PointCount)
    {
    	fprintf(stdout, "CobwebPreCoalesce file is not correct. %d frames in the file, %d frames in the trajectory. \n", i, This->PointCount);
    	return 0;
    }
    Result = fgets(FileLine, BLOCK_SIZE ,TreeFile);
    cpos = strchr(FileLine, ':');
    sscanf(cpos+1, "%d", &i);
    if (i != AttributeCount)
    {
    	fprintf(stdout, "CobwebPreCoalesce file is not correct. %d attributes in the file, %d attribute in the ptraj input file. \n", i, AttributeCount);
    	return 0;
    }
    Result = fgets(FileLine, BLOCK_SIZE ,TreeFile);
    cpos = strchr(FileLine, ':');
    sscanf(cpos+1, "%f", &f1);
    if (abs(This->Acuity - f1) > 0.000001)
    {
      fprintf(stdout, "CobwebPreCoalesce file is not correct. The acuity value in the file is %f, %f in the clustering. \n", f1, This->Acuity);
    	return 0;
    }
    Result = fgets(FileLine, BLOCK_SIZE ,TreeFile);
    Result = fgets(FileLine, BLOCK_SIZE ,TreeFile);
    sscanf(FileLine, "%f %f", &f1, &f2);
    if (abs(ClusteringGetAttributeValue(This,1,0)-f1) > 0.000001)
    {
    	fprintf(stdout, "CobwebPreCoalesce file is not correct. The first attribute value for frame 1 in the file is %f, %f in the attributeArray. \n", f1, ClusteringGetAttributeValue(This,1,0));
    	return 0;
    }
    if (abs(ClusteringGetAttributeValue(This,1,AttributeCount-1)-f2) > 0.000001)
    {
    	fprintf(stdout, "CobwebPreCoalesce file is not correct. The first attribute value for frame 1 in the file is %f, %f in the attributeArray. \n", f2, ClusteringGetAttributeValue(This,1,AttributeCount-1));
    	return 0;
    }
	return 1;
}

   void
DebugWriteCobwebTreeToFile(Clustering* This, char* FileName, CobwebTreeNode* Root)
{
  FILE* OutputFile;
  int AttributeCount = ClusteringGetAttributeCount(This);
  /**/
  OutputFile = fopen(FileName,"w");
  if (!OutputFile) 
  {
    fprintf(stdout, "Cannot write CobwebPreCoalesce file to %s\n", FileName);
    return;
  }
  /* The next five lines will used for checking */
  fprintf(OutputFile, "Frames: %d \n", This->PointCount);   
  fprintf(OutputFile, "Attributes: %d  \n", AttributeCount);
  fprintf(OutputFile, "Acuity: %f \n", This->Acuity);   
  fprintf(OutputFile, "Next line contains the first and the last attribute value for the frame 1 (the second frame).\n");
  fprintf(OutputFile, "%f %f\n", ClusteringGetAttributeValue(This,1,0), ClusteringGetAttributeValue(This,1,AttributeCount-1));
  
  DebugWriteCobwebTree(This, OutputFile, Root, 0, -1);
  fclose(OutputFile);
}


  /*
   *  Given a Host node, add the NewChild as one of its descendants.  If the host is
   *  currently a leaf, replace it in the tree with a parent of Host and NewChild.  If
   *  the host is not a leaf, consider whether Host or one its current children will 
   *  make the best host for NewChild.  (Recursive; we may add NewChild several layers down)
   *  Assumed: NewChild is a single frame!
   */
   void
CobwebClusterHostNode(Clustering* This,CobwebTreeNode* Host,CobwebTreeNode* NewChild)
{
   CobwebTreeNode* NewNode;
   CobwebTreeNode* OldParent;
   CobwebTreeNode* Parent;
   if (Host->Cluster->PointIndex==COBWEB_NONLEAF)
   {
       ClusteringCobwebAddNewFrame(This,Host,NewChild,1);
   }
   else
   {
       /* The current host is a leaf.  Make a shared-parent! */
       OldParent = Host->Parent;
       CobwebTreeNodeRemoveChild(Host);
       NewNode = CobwebTreeNodeNew();
       NewNode->Cluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
       CobwebTreeNodeAddChild(NewNode,Host);
       CobwebTreeNodeAddChild(NewNode,NewChild);
       CobwebClusterComputeMeans(This,NewNode);
       CobwebTreeNodeAddChild(OldParent,NewNode);
       /* Recursively FIX the means of higher-up nodes in the tree.  (Insertions at the 
        branch ends don't justify screwing up the means) */
       for (Parent = OldParent; Parent; Parent=Parent->Parent)
       {
           CobwebClusterComputeMeans(This,Parent);
       }
   }
}

/* Helper for AddNewFrame.  Consider merging the two best hosts, if that would improve CU */
CobwebTreeNode* CobwebConsiderMerge(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* BestHost, 
                                    CobwebTreeNode* NextBestHost, float StatusQuoCU)
{
   CobwebTreeNode* MergedMaster;
   CobwebCluster* TempMeansCluster = NULL;
   float CU;
   /**/
   fprintf(stdout,"Considering merge...\n");
   TempMeansCluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
   CobwebTreeNodeRemoveChild(BestHost);
   CobwebTreeNodeRemoveChild(NextBestHost);
   MergedMaster = CobwebClusterMergeNodes(This,BestHost, NextBestHost);
   /* Temporarily insert this node, in place of its sources */
   CobwebTreeNodeAddChild(Root,MergedMaster);
   CobwebClusterComputeMeans(This,MergedMaster);
   CobwebClusterCopyMeans(This,Root->Cluster,TempMeansCluster);
   CobwebClusterComputeMeans(This,Root);
   CU = ClusteringCobwebCategoryUtility(This,Root);
   if (CU <= StatusQuoCU)
   {
       /* No - merging will not improve our CU.  Undo the mischief we've done to our tree: */
       CobwebTreeNodeRemoveChild(MergedMaster);
       ClusteringCobwebClusterFree(This,MergedMaster->Cluster);
       CobwebTreeNodeFreeNR(MergedMaster);
       CobwebClusterCopyMeans(This,TempMeansCluster,Root->Cluster);
       CobwebTreeNodeAddChild(Root,BestHost);
       CobwebTreeNodeAddChild(Root,NextBestHost);
       ClusteringCobwebClusterFree(This,TempMeansCluster);
       return NULL;
   }
   /* Yes - merge is good!*/
   ClusteringCobwebClusterFree(This,TempMeansCluster);

   return MergedMaster;
}

/* Helper function for AddNewFrame.  Consider splitting the best host - that is, removing
   the node BestHost, and inserting all its children where it used to be */
int CobwebConsiderSplit(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode*BestHost, float StatusQuoCU)
{
  CobwebTreeNode* PrevChild;
  CobwebTreeNode* Child;
  CTNList* TempList;
  CTNList* TempListItem;
  float CU;
  CobwebCluster* TempMeansCluster = NULL;
  /**/
  fprintf(stdout,"Considering split...\n");
  TempMeansCluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
  TempList = CobwebTreeNodeListChildren(BestHost);
  CobwebTreeNodeRemoveChild(BestHost);
  
  PrevChild = NULL;
  for (Child=BestHost->FirstChild; Child; Child=Child->Next)
  {
      if (PrevChild)
      {
          CobwebTreeNodeRemoveChild(PrevChild);
          CobwebTreeNodeAddChild(Root,PrevChild);
      }
      PrevChild=Child;
  }
  if (PrevChild)
  {
      CobwebTreeNodeRemoveChild(PrevChild);
      CobwebTreeNodeAddChild(Root,PrevChild);
  }

  CobwebClusterCopyMeans(This,Root->Cluster,TempMeansCluster);
  CobwebClusterComputeMeans(This,Root);
  CU = ClusteringCobwebCategoryUtility(This,Root);
  if (CU <= StatusQuoCU)
  {
      /* No - splitting will not improve our CU.  Undo the mischief we've done to our tree: */
      for (TempListItem=TempList; TempListItem; TempListItem=TempListItem->Next)
      {
          CobwebTreeNodeRemoveChild(TempListItem->Node);
          CobwebTreeNodeAddChild(BestHost,TempListItem->Node);
      }
      CTNListFree(TempList);
      CobwebClusterCopyMeans(This,TempMeansCluster,Root->Cluster);
      CobwebTreeNodeAddChild(Root,BestHost);
      ClusteringCobwebClusterFree(This,TempMeansCluster);
      return 0;
  }
  else
  {
      /* Looks like that zany splitting scheme of ours was a good 
         plan, after all! */
      ClusteringCobwebClusterFree(This,BestHost->Cluster);
      CobwebTreeNodeFreeNR(BestHost);
      ClusteringCobwebClusterFree(This,TempMeansCluster);

      return 1;
  }
  
}

/* Helper function for COBWEB clustering: Add the specified point to the
   cluster tree.  Called recursively with a lower-in-the-tree root. */
void BASEClusteringCobwebAddNewFrame(Clustering* This,CobwebTreeNode* Root,CobwebTreeNode* Frame, int level)
{
    CobwebTreeNode* BestHost = NULL;
    float BestHostCU = 0;
    int RootChildCount = 0;
    CobwebTreeNode* NextBestHost = NULL;
    float NextBestHostCU = 0;
    float RootChildCU;
    CobwebTreeNode* PrevChild;
    CobwebTreeNode* Child;
    CobwebTreeNode* TempTree;
    float StatusQuoCU = 0; /* CategoryUtility */
    float CU = 0; /* CategoryUtility */
    CobwebTreeNode* MergedMaster;
    float* OldMeans;
    float* OldStddevs;
    CTNList* TempList;
    CTNList* TempListItem;
    CobwebTreeNode* Parent;
    int Result;
    /* 
       TempMeansCluster is used to hold the Means and Stddevs arrays when
       we're testing some tree manipulation.
    */
    CobwebCluster* TempMeansCluster = NULL;

    /**/
    /* If there are 0 or 1 children, simply add a child! */
    if (!Root->FirstChild || Root->FirstChild==Root->LastChild)
    {
        CobwebTreeNodeAddChild(Root,Frame);
        CobwebClusterComputeMeans(This, Root);
        if (Root->Parent)
        {
            for (Parent = Root->Parent; Parent; Parent=Parent->Parent)
            {
                CobwebClusterComputeMeans(This,Parent);
            }
        }
        return;
    }
    TempMeansCluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
    /* Consider adding a child to the root, as that may be best: */
    BestHost = Root;
    BestHostCU = TryCategoryUtility(This,Root, Root,Frame);
    /* 
    Ok: We have at least 2 children.  Consider the CUs of adding the new node to each of 
    our children!  Keep track of the best host, and the runner-up. 
    */
    for (Child=Root->FirstChild; Child; Child=Child->Next)
    {
        CU = TryCategoryUtility(This, Root, Child, Frame); 
        if (CU > BestHostCU)
        {
            NextBestHost = BestHost;
            NextBestHostCU = BestHostCU;
            BestHost = Child;
            BestHostCU = CU;

        }
        else if (CU > NextBestHostCU)
        {
            NextBestHost = Child;
            NextBestHostCU = CU;
        }
        RootChildCount++;
    }
    if (BestHost!=Root)
    {
        StatusQuoCU = ClusteringCobwebCategoryUtility(This,Root);       
        if (NextBestHost!=Root && RootChildCount>2)
        {
          /* Try merging the best host with the runner-up: */
          MergedMaster = CobwebConsiderMerge(This, Root, BestHost, NextBestHost, StatusQuoCU);
          if (MergedMaster)
          {
            /* We merged! Add the frame as a child of the merged node. */
            ClusteringCobwebAddNewFrame(This,MergedMaster,Frame, level+1);
            ClusteringCobwebClusterFree(This,TempMeansCluster);
            return;
          }
        }
        /*************************************************/
        /* Consider splitting the best host!  It may be a good idea to replace it in 
           the tree with all-its-children */
        if (BestHost->Cluster->PointIndex==COBWEB_NONLEAF)
        {
              Result = CobwebConsiderSplit(This, Root, BestHost, StatusQuoCU);
            if (Result)
            {
              /* We did the split, and it was good.  Start adding the frame again! */
              ClusteringCobwebClusterFree(This,TempMeansCluster);
              ClusteringCobwebAddNewFrame(This,Root,Frame,level);           
              return;
            }
        }
    }
    /**************************************************/
    /* Add the frame node as a child of the best host */
    if (BestHost==Root)
    {
        CobwebTreeNodeAddChild(Root,Frame);
        CobwebClusterComputeMeans(This,Root);
        if (Root->Parent)
        {
            for (Parent = Root->Parent; Parent; Parent=Parent->Parent)
            {
                CobwebClusterComputeMeans(This,Parent);
            }
        }

    }
    else
    {
        CobwebClusterHostNode(This,BestHost,Frame);
    }
    ClusteringCobwebClusterFree(This,TempMeansCluster);

}

/* TryCategoryUtility temporarily remodels the tree, but restores the old state when it's done! 
   Root is the node whose CU we care about.  Host is either Root, or one of Root's immediate children.
   
*/
CLUSTER_API float TryCategoryUtility(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* Host, CobwebTreeNode* Frame)
{
    CobwebTreeNode* NewNode;
    CobwebTreeNode* ChildA;
    CobwebTreeNode* ChildB;
    CobwebCluster* OldParent;
    CobwebTreeNode* OldPrevious;
    CobwebTreeNode* OldNext;
    CobwebCluster* Swap1 = NULL;
    CobwebCluster* Swap2 = NULL;

    float CU;

    if (Host->Cluster->PointIndex==COBWEB_NONLEAF)
    {
        /* We're adding yet-another-child to Host, and checking the category utility of Root */
        CobwebTreeNodeAddChild(Host,Frame);
        /* Temporarily change our view of the world */
        Swap1 = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
        CobwebClusterCopyMeans(This,Host->Cluster,Swap1);
        CobwebClusterComputeMeans(This,Host);
        if (Root!=Host)
        {
            Swap2 = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
            CobwebClusterCopyMeans(This,Root->Cluster,Swap2);            
            CobwebClusterComputeMeans(This,Root);
        }
        CU = ClusteringCobwebCategoryUtility(This,Root);
        /* Swap back the old stuff */
        CobwebClusterCopyMeans(This,Swap1,Host->Cluster);
        if (Root!=Host)
        {
           CobwebClusterCopyMeans(This,Swap2,Root->Cluster);
        }
        CobwebTreeNodeRemoveChild(Frame);

    }
    else
    {
        /* The current host is a leaf.  Make a shared-parent! 
        We do some black magic here: Higher up in the call-stack, we're iterating over the linked list
        of potential hosts.  We do NOT want to throw off that iteration by rearranging the list!  So,
        the way that we remove and re-add Host from its listing is not standard.
        */
        /* Create a parent for these two kids: */
        NewNode = CobwebTreeNodeNew();
        NewNode->Cluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
        ChildA = CobwebTreeNodeNew();
        ChildA->Cluster = Host->Cluster;
        ChildA->TotalFrameCount = Host->TotalFrameCount;
        CobwebTreeNodeAddChild(NewNode,ChildA);
        ChildB = CobwebTreeNodeNew();
        ChildB->Cluster = Frame->Cluster;
        ChildB->TotalFrameCount = Frame->TotalFrameCount;
        CobwebTreeNodeAddChild(NewNode,ChildB);
        CobwebClusterComputeMeans(This,NewNode);
        /* Pluck the host off the root, and put the newnode on, and see how that looks: */
        CobwebTreeNodeSwap(Host,NewNode);
        Swap1 = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
        CobwebClusterCopyMeans(This,Root->Cluster,Swap1);
        
        CobwebClusterComputeMeans(This,Root);
        CU = ClusteringCobwebCategoryUtility(This,Root);
        /* Put things back the way they were: */
        
        CobwebTreeNodeSwap(NewNode,Host);
        CobwebClusterCopyMeans(This,Swap1,Root->Cluster);
        NewNode->Previous=NULL;
        NewNode->Next=NULL;
        NewNode->Parent=NULL;
        ClusteringCobwebClusterFree(This,NewNode->Cluster);
        CobwebTreeNodeFreeNR(NewNode);
        CobwebTreeNodeFreeNR(ChildA);
        CobwebTreeNodeFreeNR(ChildB);
    }
    Frame->Parent = NULL;
    if (Swap1) ClusteringCobwebClusterFree(This, Swap1);
    if (Swap2) ClusteringCobwebClusterFree(This, Swap2);
    return CU;
}

/* CoalescentNode is a struct that we use in order to transform a CobwebTree into a 
   flat collection of disjoint clusters. */
typedef struct CoalescentNode
{
    struct CobwebTreeNode* Node;
    struct CoalescentNode* Previous;
    struct CoalescentNode* Next;
    int PlannedClusters;
    float CU;
    float CUTotal;
} CoalescentNode;

CoalescentNode* CoalescentNodeNew()
{
    CoalescentNode* This;
    This = (CoalescentNode*)malloc(sizeof(CoalescentNode));
    memset(This,0,sizeof(CoalescentNode));
    return This;
}

/* Helper function
   Update the values of Root as if the FrameNode is added. 
   The original values of Root should be save before this function is called so that it can be restored. Or we can reverse the action by another function.
   Only the values of the Root will be changed.
*/
CLUSTER_API void CobwebIncorporate(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* FrameNode)
{
	CobwebTreeNode *Temp;
    int AttributeIndex, AttributeCount;
    
    
}

/* Helper function. 
   Try merge two children of Root, report the new CU for Root. 
   Return the CU.
   Will free the newly merged node and restore the original state of the cobweb tree. 
   If execute == 1, merge.
   FrameNode should contain only one point, and it is a child of Root. 
   If FrameNode is NULL, ignore it. 
*/
CLUSTER_API float CobwebTryMerge(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* Child1, CobwebTreeNode* Child2, int execute)
{
	CobwebTreeNode* MergedMaster;
    float CU;
    /*
    fprintf(stdout, "Considering merging...\n");
    */
    if (!Child1 || !Child2 || !Root) return -1000;
    if (Child1 == Child2) return -1000;
    CobwebTreeNodeRemoveChild(Child1);
    CobwebTreeNodeRemoveChild(Child2);
    MergedMaster = CobwebClusterMergeNodes(This, Child1, Child2);
    CobwebTreeNodeAddChild(Root, MergedMaster);
    
    CobwebClusterComputeMeans(This, MergedMaster);
    /*CobwebClusterComputerMeans(This, Root); Not necessary. */
    if (!execute)
	    CU = ClusteringCobwebCategoryUtility(This, Root);
    if (!execute)
    {
    	CobwebTreeNodeRemoveChild(MergedMaster);
        ClusteringCobwebClusterFree(This, MergedMaster->Cluster);
        CobwebTreeNodeFreeNR(MergedMaster);
        CobwebTreeNodeAddChild(Root, Child1);
        CobwebTreeNodeAddChild(Root, Child2);
    }
	return CU;    
}

/* Helper function.
   Try to split a Node. Promote all its children to replace the parent node.
   Return the CU for Root. 
   Node should be Root's child.
   Will restore the original cobweb tree.
   If execute == 1, do not restore.
*/
CLUSTER_API float CobwebTrySplit(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* Node, int execute)
{
	CobwebTreeNode* Child;
    CobwebTreeNode* PrevChild;
    CTNList* TempList;
    CTNList* TempListItem;
    float CU = 0;
    
    if (!Node || !Root) return -1000;
    if (Node->FirstChild == Node->LastChild) return -1000;
    TempList = CobwebTreeNodeListChildren(Node);
    CobwebTreeNodeRemoveChild(Node);
    
	PrevChild = NULL;
    for (Child = Node->FirstChild; Child; Child = Child->Next)
    {
    	if (PrevChild)
        {
        	CobwebTreeNodeRemoveChild(PrevChild);
            CobwebTreeNodeAddChild(Root, PrevChild);
        }
    	PrevChild = Child;
    }
    if (PrevChild) 
    {
       	CobwebTreeNodeRemoveChild(PrevChild);
        CobwebTreeNodeAddChild(Root, PrevChild);
    }
    
    /*CobwebClusterComputerMeans(This, Root); Not necessary. */
    if (!execute)
	    CU = ClusteringCobwebCategoryUtility(This, Root);
    if (!execute)
    {
    	for (TempListItem = TempList; TempListItem; TempListItem = TempListItem->Next)
        {
        	CobwebTreeNodeRemoveChild(TempListItem->Node);
            CobwebTreeNodeAddChild(Node, TempListItem->Node);
        }
        CobwebTreeNodeAddChild(Root, Node);
    } else {
    	ClusteringCobwebClusterFree(This, Node->Cluster);
        CobwebTreeNodeFreeNR(Node);
    }
	CTNListFree(TempList);
    return CU;
}

/* Helper function. 
   Merge two children of Root, report the new CU for Root. 
   Will try to split one of the two nodes. If CU increases, try to split both.
*/
CLUSTER_API float CobwebMerge(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* Child1, CobwebTreeNode* Child2)
{
	CobwebTreeNode* MergedMaster;
    float CU11, CU12, CU21, CU22;
    /*
    fprintf(stdout, "Considering merging...\n");
    */
    if (!Child1 || !Child2 || !Root) return -1000;
    if (Child1 == Child2) return -1000;
    CobwebTreeNodeRemoveChild(Child1);
    CobwebTreeNodeRemoveChild(Child2);
    MergedMaster = CobwebClusterMergeNodes(This, Child1, Child2);
    CobwebTreeNodeAddChild(Root, MergedMaster);
    
    CobwebClusterComputeMeans(This, MergedMaster);
    /*CobwebClusterComputerMeans(This, Root); Not necessary. */
    CU11 = ClusteringCobwebCategoryUtility(This, MergedMaster);
    CU21 = CobwebTrySplit(This, MergedMaster, Child1, 0);
    CU12 = CobwebTrySplit(This, MergedMaster, Child2, 0);
    
    if ((CU21 > CU12) && (CU21 > CU11)) {
    	CobwebTrySplit(This, MergedMaster, Child1, 1);
        CU22 = CobwebTrySplit(This, MergedMaster, Child2, 0);
        if (CU22 > CU21) {
        	CobwebTrySplit(This, MergedMaster, Child2, 1);
            return CU22;
        }
        return CU21;
    }
    if ((CU12 > CU21) && (CU12 > CU11)) {
    	CobwebTrySplit(This, MergedMaster, Child2, 1);
        CU22 = CobwebTrySplit(This, MergedMaster, Child1, 0);
        if (CU22 > CU21) {
        	CobwebTrySplit(This, MergedMaster, Child1, 1);
            return CU22;
        }
        return CU12;
    }
    /* Here I assume for CU22 (both clusters split) to be greater than CU11 (none of clusters split), either CU12 or CU21 would be greater than CU11.
      I also ignore the merging after split, because it will be taken care of if new frame is added into this branch. */
    return CU11;  
}

CLUSTER_API void CobwebTreeReorganize(Clustering* This, CobwebTreeNode* Root, int DesiredClusterCount)
{
	int ChildCount;
    CobwebTreeNode *MergeNode, *MergeHost, *MergeNode1, *MergeNode2;
    CobwebTreeNode *NextBestNode, *BestCUNode;
    CobwebTreeNode *Child;
    float MergeCU, BestCU, NextBestCU, CU, RootCU;
    CTNList *MergeList, *ChildList, *Item1, *Item2, *Prev;
    
    if (DesiredClusterCount == 1) {
    	return;
    }
    if (DesiredClusterCount > Root->TotalFrameCount) {
		fprintf(stdout, "Warning: We're being asked for %d clusters, but we only have %d frames.   Ratcheting down to %d clusters!\n", DesiredClusterCount,Root->TotalFrameCount, Root->TotalFrameCount);
        DesiredClusterCount = Root->TotalFrameCount;    
    }
    for (ChildCount = 0, Child = Root->FirstChild; Child; Child = Child->Next, ChildCount++);
    if (DesiredClusterCount <= 0) {
    	fprintf(stdout, "CobwebPreCoalesce.txt contains %d children, so form %d clusters.\n", ChildCount, ChildCount);
        DesiredClusterCount = ChildCount;
    }
#ifdef DEBUG    
    if (prnlev > 5) fprintf(stdout, "ChildCount is %d, DesiredClusterCount is %d\n", 
			    ChildCount, DesiredClusterCount);   
    if (ChildCount == DesiredClusterCount) {
      if (prnlev > 5) fprintf(stdout, "ChildCount is equal to DesiredClusterCount. Return \n");
        return;
    } 
#endif
	/* Splitting */
	if (ChildCount < DesiredClusterCount) {
    	BestCU = -1000;
        ChildList = CobwebTreeNodeListChildren(Root);
        for (Item1 = ChildList; Item1; Item1 = Item1->Next) {
        	CU = CobwebTrySplit(This, Root, Item1->Node, 0);
#ifdef DEBUG    
		if (prnlev > 5) fprintf(stdout, "CU is %f, BestCU is %f\n", CU, BestCU);   
#endif
            if (CU > BestCU) {
            	BestCU = CU;
                BestCUNode = Item1->Node;
            }
        }
        CobwebTrySplit(This, Root, BestCUNode, 1);
        CTNListFree(ChildList);
        CobwebTreeReorganize(This, Root, DesiredClusterCount);
    }
    /* Merging */
    while (ChildCount > DesiredClusterCount) {
    	MergeList = NULL;
        BestCUNode = Root->FirstChild;
        NextBestNode = NULL;
        NextBestCU = BestCU = -1000;
        RootCU = ClusteringCobwebCategoryUtility(This, Root);
        /* MergeList is used to colloct children nodes that is either a leaf or higher CU. 
           We will try to merge each pairs in the list, and realize the merge for highest CU. 
           In the worst case, the list may contain all the children if the children's CU is in ascending order.*/
        ChildList = CobwebTreeNodeListChildren(Root);
        for (Item1 = ChildList; Item1; Item1 = Item1->Next) {
        	Child = Item1->Node;
            if (Child->TotalFrameCount == 1) { 
            	MergeList = CTNListAppend(MergeList, Child);    
            } else {
            	CU = ClusteringCobwebCategoryUtility(This, Child);
                if (CU > BestCU) {
                	NextBestCU = BestCU;
                    BestCU = CU;
                    BestCUNode = Child;
                    MergeList = CTNListAppend(MergeList, Child);
                } else if (CU > NextBestCU) {
                	NextBestCU = CU;
                    NextBestNode = Child;
                }
            }
        }
        /* If DesiredClusterCount < 0, keep merging until CU cannot be improved by merging. */
        if (DesiredClusterCount < 0 && RootCU >= BestCU) break;
        
        BestCU = -1000;
        for (Item1 = MergeList; Item1; Item1 = Item1->Previous) {
        	for (Item2 = ChildList; Item2; Item2 = Item2->Next) {
            	if (Item1->Node == Item2->Node) continue;
                CU = CobwebTryMerge(This, Root, Item1->Node, Item2->Node, 0);
                if (CU > BestCU) {
                	BestCU = CU;
                    MergeNode1 = Item1->Node;
                    MergeNode2 = Item2->Node;
                }
            }
        }
        Prev = NULL;
        for (Item1 = MergeList; Item1; Item1 = Item1->Previous) {
        	safe_free(Prev);
            Prev = Item1;
        }
        safe_free(Prev);
        CTNListFree(ChildList);
        CobwebTryMerge(This, Root, MergeNode1, MergeNode2, 1);
        ChildCount--;
    }
}

/* 
Now that we've produced a big tree of clusters, we want to "flatten" the tree into clusters that partition the frames.
Here is the recursive plan:
- If you want as many clusters as you have children, make each child into a cluster; you're done!
- If you want fewer clusters than you have children: MERGE your smallest (in total-frame-count) 
cluster into the best (in terms of impact on CU) other cluster.  Keep merging until the number 
of children is correct.
- If you want more clusters than you have chilren: Each child will get one cluster.  After than, 
divvy up the clusters based upon the CU of each child (a high CU gives you fewer clusters, 
since you're already pretty good).

Returns the new CTNList tail.
*/
CTNList* CobwebTreeCoalesce(Clustering* This, CobwebTreeNode* Root, CTNList* Tail, int DesiredClusterCount)
{
    int ChildCount = 0;
    int BestSize = 0;
    CobwebTreeNode* MergeNode;
    CobwebTreeNode* MergeHost;
    float OriginalCU;
    float MergedCU;
    float BestMergedCU;
    CTNList* TempList;
    CTNList* TempListItem;
    
    CTNList* NewCTNNode;
    CobwebTreeNode* BestHost;
    float* OldMeans;
    float* OldStddevs;
    CoalescentNode* CHead;
    CoalescentNode* CTail;
    CoalescentNode* CNode;
    CoalescentNode* CPrev;
    CobwebTreeNode* Child;
    CTNList* ChildNode;
    float BestCUTotal;
    CoalescentNode* ParcelNode;
    int ClustersToAssign;
    CobwebCluster* SwapClusterMeans = NULL;
    /**/
    if (DesiredClusterCount==1)
    {
        return CTNListAppend(Tail,Root);
    }
    if (DesiredClusterCount > Root->TotalFrameCount)
    {
        fprintf(stdout,"Warning: We're being asked for %d clusters, but we only have %d frames.  Ratcheting down to %d clusters!",
            DesiredClusterCount,Root->TotalFrameCount,Root->TotalFrameCount);
        DesiredClusterCount = Root->TotalFrameCount;
        /*Root->TotalFrameCount = DesiredClusterCount;*/
    }
    if (Root->FirstChild)
    {
        for (Child = Root->FirstChild; Child; Child = Child->Next)
        {
            ChildCount++;
        }
    }

    /********************************************/
    /* MERGING!  */
    while (ChildCount > DesiredClusterCount)
    {
        OriginalCU = ClusteringCobwebCategoryUtility(This,Root);
        BestMergedCU = 0;
        BestHost = NULL;
        /* Pick the smallest cluster; it will be merged into another: */
        BestSize = 0;
        for (Child = Root->FirstChild; Child; Child = Child->Next)
        {
            if (BestSize==0 || Child->TotalFrameCount < BestSize)
            {
                BestSize = Child->TotalFrameCount;
                MergeNode = Child;
            }
        }
        /* Try every other cluster as a host: */
        CobwebTreeNodeRemoveChild(MergeNode);
        TempList = CobwebTreeNodeListChildren(Root); /* So that we can iterate IN ORDER */
        for (TempListItem = TempList; TempListItem; TempListItem=TempListItem->Next)
        {
            Child = TempListItem->Node;
            if (Child==MergeNode)
            {
                continue; /* Don't try to merge with yourself */
            }
            CobwebTreeNodeRemoveChild(Child);
            MergeHost = CobwebClusterMergeNodes(This,Child,MergeNode);
            /* Temporarily insert this node, in place of its sources */
            CobwebTreeNodeAddChild(Root,MergeHost);
                SwapClusterMeans = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
               CobwebClusterCopyMeans(This,Root->Cluster,SwapClusterMeans);
            CobwebClusterComputeMeans(This,MergeHost);
            CobwebClusterComputeMeans(This,Root);
            MergedCU = ClusteringCobwebCategoryUtility(This,Root);
            if (MergedCU > BestMergedCU || BestHost == NULL)  /* CU could be negative if StddevChild > StddevParent. */
            {
                BestHost = Child;
                BestMergedCU = MergedCU;
            }
            /* Undo the merge */
            CobwebTreeNodeRemoveChild(MergeHost);
            CobwebTreeNodeAddChild(Root,Child);
            ClusteringCobwebClusterFree(This,MergeHost->Cluster);
            MergeHost->Previous=NULL;
            MergeHost->Next=NULL;
            MergeHost->Parent=NULL;
            CobwebTreeNodeFreeNR(MergeHost);
            CobwebClusterCopyMeans(This,SwapClusterMeans,Root->Cluster);          
            ClusteringCobwebClusterFree(This,SwapClusterMeans);
        }
        CTNListFree(TempList);
        /* Finalize the merge */
        CobwebTreeNodeRemoveChild(BestHost);
        MergeHost = CobwebClusterMergeNodes(This,BestHost,MergeNode);
        CobwebTreeNodeAddChild(Root,MergeHost);
        CobwebClusterComputeMeans(This,MergeHost);
        CobwebClusterComputeMeans(This,Root);
        ChildCount--;
        
	if (prnlev > 5) {    
        fprintf(stdout,"\nmerge cobweb \n");
	DebugPrintCobwebTree(This,Root,0,3);
	}
            
    }

    /* If we just did some merges - or if we had just the right number of children - then 
       we're basically done: */
    if (ChildCount == DesiredClusterCount)
    {
        for (Child = Root->FirstChild; Child; Child = Child->Next)
        {
            Tail = CobwebTreeCoalesce(This,Child,Tail,1);
        }
        return Tail;
    }

    /********************************************/
    /* SPLITTING! During each cycle, the "worst" cluster (meaning the one with 
       lowest CU) gets split */
    if (ChildCount < DesiredClusterCount)
    {
        CHead = NULL;
        CTail = NULL;
        ClustersToAssign = DesiredClusterCount; /* - ChildCount;*/
        /* Build a list CoalescentNodes to keep track of who's who: */
        for (Child = Root->FirstChild; Child; Child = Child->Next)
        {
            CNode = CoalescentNodeNew();
            CNode->Node = Child;
            CNode->PlannedClusters = 1;
            ClustersToAssign--;
            CNode->CU = ClusteringCobwebCategoryUtility(This,Child);
            CNode->CUTotal = CNode->CU;
            CNode->Previous = CTail;
            if (CTail)
            {
                CTail->Next = CNode;
            }
            else
            {
                CHead = CNode;
            }
            CTail = CNode;
        }
        /* Assign clusters to our children */
        while (ClustersToAssign > 0)
        {
            BestCUTotal = 0;
            for (CNode = CHead; CNode; CNode = CNode->Next)
            {
                /* Don't assign more than we can hold: */
                if (CNode->Node->TotalFrameCount <= CNode->PlannedClusters)
                {
                    continue;
                }
                /* Find the one with the lowest running-CU-total */
                if (BestCUTotal==0 || CNode->CUTotal < BestCUTotal)
                {
                    BestCUTotal = CNode->CUTotal;
                    ParcelNode = CNode;
                }
            }
            ParcelNode->PlannedClusters++;
            ParcelNode->CUTotal += ParcelNode->CU;
            ClustersToAssign--;
        }
        /* And now, ask our children for some nice clusters: */
        CPrev = NULL;
        for (CNode = CHead; CNode; CNode = CNode->Next)
        {
	  if (prnlev > 5) {
	    fprintf(stdout,"\nsplit cobweb into %d subnodes\n",CNode->PlannedClusters);
	    DebugPrintCobwebTree(This,CNode->Node,0,3);
	  }
	  Tail = CobwebTreeCoalesce(This,CNode->Node,Tail,CNode->PlannedClusters);

	  if (prnlev > 5) DebugPrintCobwebTree(This,CNode->Node,0,3);

            if (CPrev)
            {
                safe_free(CPrev);
            }
            CPrev = CNode;
        }
        return Tail;
    }
    return NULL; /*unreachable*/
}

void CobwebFlattenTree(Clustering* This, CobwebTreeNode* Root, int DesiredClusterCount)
{
  CTNList* CoalescentTail;
  int ClusterIndex;
  CTNList* CNode;
  int Backtracking;
  CobwebTreeNode* TreeNode;
  Cluster* TempCluster;
  /**/
 
  /*
  CoalescentTail = CobwebTreeCoalesce(This,Root,NULL,DesiredClusterCount);
  */
  CobwebTreeReorganize(This, Root, DesiredClusterCount);
  CoalescentTail = CobwebTreeNodeListChildren(Root);
  for (CNode = CoalescentTail; CNode; CNode = CNode->Next) 
  {
  	  CoalescentTail = CNode;
  }
  
  
  DebugWriteCobwebTreeToFile(This, "CobwebCoalesce.txt", Root);
  ClusterIndex = 0;
  for (CNode=CoalescentTail; CNode; CNode = CNode->Previous)
  {
      TempCluster = ClusteringGetNewCluster(This);
      ClusteringAddCluster(This,TempCluster);
      TreeNode = CNode->Node;
      Backtracking=0;        
      while (1)
      {
            /* After hitting a node with no child or next, we backtrack up the 
               tree to the next time we can move to a next-node. */
          if (Backtracking == 1)
          {
              if (TreeNode==CNode->Node)
              {
                  break;
              }
              TreeNode = TreeNode->Parent;
                if (!TreeNode || TreeNode==CNode->Node)
                {
                    break;
                }
              if (TreeNode->Next)
              {
                  TreeNode = TreeNode->Next;
                  Backtracking = 0; 
              }
              else
              {
                  continue;
              }
          }
          if (TreeNode->Cluster->PointIndex != COBWEB_NONLEAF)
          {
              ClusterAddMember(TempCluster,TreeNode->Cluster->PointIndex);
          }
          if (TreeNode->FirstChild)
          {
              TreeNode = TreeNode->FirstChild;
              continue;
          }
          if (TreeNode->Next && TreeNode!=CNode->Node)
          {
              TreeNode = TreeNode->Next;
              continue;
          }
          Backtracking=1;
      }
        ClusterIndex++;
  }

}

void BASEClusteringClusterCobweb(Clustering* This, int DesiredClusterCount)
{
    Cluster* TempCluster;
    int PointIndex;
    int PointCount = This->PointCount;
    CobwebTreeNode* Root;
    CobwebTreeNode* Frame;
    CobwebTreeNode* TreeNode;
    CobwebCluster* Cluster;
    CobwebCluster* FrameCluster;
    CTNList* CoalescentTail;
    int ClusterIndex;
    CTNList* CNode;
    int Backtracking;
    int *FinishedPoints;
    int ProcessIndex;
    /**/
    fprintf(stdout,"\n\nCOBWEB start!\n");
    Root = ReadCobwebTree(This, "CobwebPreCoalesce.txt");
    if (Root) 
    {
      if (prnlev > 5) DebugPrintCobwebTree(This,Root,0, -1);
    	CobwebFlattenTree(This, Root, DesiredClusterCount);
	
        if (prnlev > 5) {
	  fprintf(stdout,"\n\nFinal Cobweb tree:\n");
	  DebugPrintCobwebTree(This,Root,0, -1);
	}
    	return;
    }
    Root = CobwebTreeNodeNew();
    Root->Cluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
    
    FinishedPoints = (int*)malloc(sizeof(int)*PointCount);
    memset(FinishedPoints, 0, sizeof(int)*PointCount);
    
    for (ProcessIndex=0; ProcessIndex<PointCount; ProcessIndex++)
    {
    	PointIndex = ProcessIndex;
    	/*PointIndex = ChooseNextPoint(FinishedPoints, PointCount,PointCount-ProcessIndex);*/
    /*
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
    */
        if (prnlev > 5) fprintf(stdout,"\n\n\n------------\nTry point %d:\n",PointIndex);
        /* Add this frame, somehow, to the root node.  This call may recurse. */
        Frame = CobwebTreeNodeNewLeaf(This,PointIndex);
        /*Frame->Cluster = NewCobwebClusterLeafNode(This,PointIndex);*/
        ClusteringCobwebAddNewFrame(This,Root,Frame,1);
        if (prnlev > 5) {
	  fprintf(stdout,"Added frame %d:\n",PointIndex);
	  DebugPrintCobwebTree(This,Root,0, -1);
	}
    }
    DebugWriteCobwebTreeToFile(This, "CobwebPreCoalesce.txt", Root);

    /* We have our tree; we'll coalesce it into managable clusters */
    CobwebFlattenTree(This, Root, DesiredClusterCount);
    if (prnlev > 5) {
      fprintf(stdout,"\n\nFinal Cobweb tree:\n");
      DebugPrintCobwebTree(This,Root,0, -1);
    }
    
    safe_free(FinishedPoints);
}

float BASEClusteringCobwebCategoryUtility(Clustering* This, CobwebTreeNode* Parent)
{
  /* Abstract method: Do nothing! */
  return 1.0;
}

void BASECobwebClusterCopyMeans(Clustering* This, CobwebCluster* Source, CobwebCluster* Destination)
{
  /* Abstract method - does nothing! */
  
}

BayesianCluster* BayesianClusterNew(Clustering* TheClustering)
{
  BayesianCluster* This;
  This = (BayesianCluster*)malloc(sizeof(BayesianCluster));
  memset(This,0,sizeof(BayesianCluster));
  This->Members = (double*)malloc(sizeof(double)*TheClustering->PointCount);
  This->Means = (double*)malloc(sizeof(double)*ClusteringGetAttributeCount(TheClustering));
  This->Stddevs = (double*)malloc(sizeof(double)*ClusteringGetAttributeCount(TheClustering));
  return This;
}

void BayesianClusterFree(BayesianCluster* This)
{
  if (!This)
  {
    return;
  }
  safe_free(This->Means);
  safe_free(This->Stddevs);
  safe_free(This->Members);
  safe_free(This);
}

void BayesianClusterListFree(BayesianCluster* Head)
{
  BayesianCluster* Prev = NULL;
  BayesianCluster* Node;
  for (Node = Head; Node; Node=Node->Next)
  {
    if (Prev)
    { 
      BayesianClusterFree(Prev);
    }
    Prev = Node;
  }
  if (Prev)
  { 
    BayesianClusterFree(Prev);
  }
}

#define ONE_OVER_SQRT_2PI 0.398942280401432

double NormalDist(double Value,double Mean,double Stddev)
{
  float Temp;
  if (Stddev==0)
  {
    return 0;
  }
  Temp = Value-Mean;
  Temp = -(Temp*Temp)/(2 * Stddev* Stddev);
  Temp = exp(Temp);
  /*fprintf(stdout,"NormalDist: %.2f,%.2f, %.2f: %.5f\n",Value,Mean,Stddev, (ONE_OVER_SQRT_2PI * Temp / Stddev));*/
  return ONE_OVER_SQRT_2PI * Temp / Stddev;
}
/* Convert Bayesian clusters, with probabilistic membership, to a collection of ordinary
   clusters that partition the space.  Each point is assigned to the cluster where it has the
   highest membership probability.*/
void ClusteringFinalizeBayesianClusters(Clustering* This, BayesianCluster* Head)
{
  BayesianCluster* BCluster;
  BayesianCluster* BestHost;
  Cluster* NewCluster;
  int PointIndex;
  float BestProb;
  int ClusterIndex;
  
  /**/
  fprintf(stdout, "  Finalizing Bayesian clusters...\n");
  for (BCluster = Head; BCluster; BCluster = BCluster->Next)
  {
    NewCluster = ClusteringGetNewCluster(This);
    ClusteringAddCluster(This,NewCluster);
    BCluster->TrueCluster = NewCluster;
  }
  for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
  {
    BestProb = 0;
    ClusterIndex=0;
    for (BCluster = Head; BCluster; BCluster = BCluster->Next)
    {
      if (BCluster->Members[PointIndex] > BestProb)
      {
        BestHost = BCluster;
        BestProb = BCluster->Members[PointIndex];
      }
    }
    ClusterAddMember(BestHost->TrueCluster,PointIndex);
  }

  fprintf(stdout, "Bayesian clusters finalized.\n");

  /* verbose output */
  /*ClusterIndex=0;
  for (BCluster = Head; BCluster; BCluster = BCluster->Next)
  {
    fprintf(stdout,"Cluster %d centered at (%f,%f) with sigma (%f,%f)\n",ClusterIndex,BCluster->Means[0],BCluster->Stddevs[0],BCluster->Means[1],BCluster->Stddevs[1]);
    fprintf(stdout,"%d,%d,%d\n",BCluster->TrueCluster->Mask[0],BCluster->TrueCluster->Mask[1],BCluster->TrueCluster->Mask[2]);
    for (PointIndex=0; PointIndex < This->PointCount; PointIndex++)
    {
      if (ClusterIsMember(BCluster->TrueCluster,PointIndex))
      {
        fprintf(stdout,"  Has point %d at (%f,%f)\n",PointIndex,ClusteringGetAttributeValue(This,PointIndex,0),ClusteringGetAttributeValue(This,PointIndex,0));
      }
    }
    ClusterIndex++;
  }*/
}

/* Re-seed the random number generator, using the system clock: */
void Randomize()
{
  srand(time(NULL) + clock());
}

BayesianCluster* GenerateBayesianSeedClusters(Clustering* This, int ClusterCount)
{
  BayesianCluster* Head;
  BayesianCluster* Tail;
  BayesianCluster* BCluster;
  double MinAttributeValue;
  double MaxAttributeValue;
  double RandValue;
  int ClusterIndex;
  int AttributeIndex;
  double Value;
  int AttributeCount;
  AttributeCount = ClusteringGetAttributeCount(This);
  /* Build our list of bayesian clusters*/
  Head = BayesianClusterNew(This);
  Tail = Head;
  for (ClusterIndex=1; ClusterIndex<ClusterCount; ClusterIndex++)
  {
    BCluster = BayesianClusterNew(This);
    Tail->Next = BCluster;
    BCluster->Prev = Tail;
    Tail = BCluster;
  }
  /* Generate (random) parameters for our bayesian clusters */
  fprintf(stdout,"  Cycle attributes.\n");
  for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
  {
    /*fprintf(stdout,"Cycle attributes:%d\n",AttributeIndex);*/
    MinAttributeValue = (double)ClusteringGetAttributeMin(This,AttributeIndex);
    MaxAttributeValue = (double)ClusteringGetAttributeMax(This,AttributeIndex);
    /*fprintf(stdout,"Min%f,Max%f\n",MinAttributeValue,MaxAttributeValue);*/
    ClusterIndex=0;
    for (BCluster = Head; BCluster; BCluster = BCluster->Next)
    {
      RandValue = rand()/(double)RAND_MAX;
      Value = MaxAttributeValue - MinAttributeValue;
      BCluster->Means[AttributeIndex] = MinAttributeValue + RandValue*Value;
      RandValue = rand()/(double)RAND_MAX;
      /* Start with a large stddev, so that probabilities aren't vanishingly small: */
      BCluster->Stddevs[AttributeIndex] = Value*0.5 + RandValue*Value*0.5;
      BCluster->Probability = 1 / (double)ClusterCount;
      /*fprintf(stdout,"Cluster %d mean %f Stddev %f\n",ClusterIndex,BCluster->Means[AttributeIndex],BCluster->Stddevs[AttributeIndex]);*/
      ClusterIndex++;
    }
  }
  return Head;
}

void ComputeBayesianMembershipProbabilities(Clustering* This, BayesianCluster* Head, int ClusterCount)
{
  BayesianCluster* BCluster;
  BayesianCluster* ToKill;
  int ClusterIndex;
  int AttributeCount;  
  int AttributeIndex;
  int PointIndex;
  double RunningProbability;
  double Probability;
  double Value;
  double ClusteringQuality;
  double ClusteringQualityOverall;
  double PreviousQuality = 0;
  int CycleCount = 0;
  Cluster* NewCluster;
  double ScalingFactor = 1.0;
  float MaxRunningProbability;
  AttributeCount = ClusteringGetAttributeCount(This);
  /**/
  /* For each point, compute the probability that it's a member of each cluster.  */
  for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
  {
    ClusterIndex = 0;
    MaxRunningProbability = -25000.0;
    for (BCluster = Head; BCluster; BCluster = BCluster->Next)
    {     
        RunningProbability = 0.0;  /*ScalingFactor;*/
        for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
        {
          Value = (double)ClusteringGetAttributeValue(This,PointIndex,AttributeIndex);
	  /*fprintf(stdout, "Frame %d has attribute %d of %.4f\n", PointIndex, AttributeIndex, Value);*/
          Probability = NormalDist(Value,BCluster->Means[AttributeIndex],
                                   BCluster->Stddevs[AttributeIndex]);
          if (Probability == 0.0)
          {
              /*fprintf(stdout,"Resetting probability from -INF!");*/
  	      RunningProbability -= 50;
	      
          }
          else
          {
              RunningProbability += log(Probability);
          }
        }  
        MaxRunningProbability = max(RunningProbability,MaxRunningProbability);
        BCluster->Members[PointIndex] = RunningProbability;
        ClusterIndex++;
    }
    /*fprintf(stdout,"Max-running-probability for all clusters was: %.4f\n", MaxRunningProbability);*/
    /* Run through again, and convert the sums-of-logs to non-normalized probs */
    for (BCluster = Head; BCluster; BCluster = BCluster->Next)
    {  
      BCluster->Members[PointIndex] -= MaxRunningProbability;
      /*fprintf(stdout,"Prob %d before expo:%.5f\n",PointIndex,BCluster->Members[PointIndex]);*/
      if (BCluster->Members[PointIndex] < -50) /* Avoid overflow! */
      {
          BCluster->Members[PointIndex] = 0;
      }
      else
      {
          BCluster->Members[PointIndex] = exp(BCluster->Members[PointIndex]);
      }
      /*fprintf(stdout,"Prob %d after expo:%.5f\n",PointIndex,BCluster->Members[PointIndex]);*/
      
    }
    /*fprintf(stdout,"RunningProbability:%.5f\n",RunningProbability);*/
  }
  /* Normalize the membership probabilities: Replace each probability p by the quotient p / TotalProb, where
     TotalProb = sum(p for each cluster) */
  for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
  {
    Value = 0;
    for (BCluster = Head; BCluster; BCluster = BCluster->Next)
    {
      Value += BCluster->Members[PointIndex];
    }
    ClusterIndex=0;
    for (BCluster = Head; BCluster; BCluster = BCluster->Next)
    {
      BCluster->Members[PointIndex] /= Value;
      ClusterIndex++;
    }
  }

  /* Recompute the probabilities of each cluster: */
  ClusterIndex = 0;
  for (BCluster = Head; BCluster; BCluster = BCluster->Next)
  {
    BCluster->Probability = 0;
    for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
    {
        BCluster->Probability += BCluster->Members[PointIndex];
    }
    BCluster->Probability /= This->PointCount;
    fprintf(stdout,"    Total probability of cluster %d is %f\n",ClusterIndex,BCluster->Probability);
    ClusterIndex++;
  }
  /* INTEGRITY CHECK: */
  Probability = 0;
  for (BCluster = Head; BCluster; BCluster = BCluster->Next)
  {
    Probability += BCluster->Probability;
  }
  fprintf(stdout,"  Probability should be near 1.0: %f\n",Probability);

  /* REMOVE any clusters with 0 chance of containing any point: */
  ToKill = NULL;
  ClusterIndex = 0;
  for (BCluster = Head; BCluster; BCluster = BCluster->Next)
  {
    if (ToKill)
    {
      BayesianClusterFree(ToKill);
      ToKill = NULL;
    }
    if (BCluster->Probability == 0.0)
    {
      fprintf(stdout,"KILLING cluster %d, because it has no chance of holding anything!\n",ClusterIndex);
      fprintf(stdout,"Mean %f stddev %f, Mean %f stddev %f\n",BCluster->Means[0],BCluster->Stddevs[0],BCluster->Means[1],BCluster->Stddevs[1]);
      if (BCluster->Prev)
      {
        BCluster->Prev->Next = BCluster->Next;
      }
      if (BCluster->Next)
      {
        BCluster->Next->Prev = BCluster->Prev;
      }
      ToKill = BCluster;
    }
    ClusterIndex++;
  }
  if (ToKill)
  {
    BayesianClusterFree(ToKill);
    ToKill = NULL;
  }

}

CLUSTER_API void ClusteringClusterBayesian(Clustering* This, int ClusterCount)
{
  BayesianCluster* Head;
  BayesianCluster* BCluster;
  int ClusterIndex;
  int AttributeCount;  
  int AttributeIndex;
  int PointIndex;
  double RunningProbability;
  double Probability;
  double Value;
  double ClusteringQuality;
  double ClusteringQualityOverall;
  double PreviousQuality = 0;
  int CycleCount = 0;
  Cluster* NewCluster;
  double ScalingFactor = 1.0;
  float MaxRunningProbability;
  /**/
  Randomize();
  AttributeCount = ClusteringGetAttributeCount(This);
  Head = GenerateBayesianSeedClusters(This, ClusterCount);
  while (1)
  {
    
    /* Get the E-value (expected probabilities) for the clusters. */
    /* Note: The normal-distribution values get very close to zero, and that makes our log-probabilities get very close
    to -INF.  We don't want to get strange overflow behavior, so we cap our log-probabilities at -50. */
    ComputeBayesianMembershipProbabilities(This, Head, ClusterCount);

    /* Maximize the utility of our clusters. */
    ClusterIndex=0;
    for (BCluster = Head; BCluster; BCluster = BCluster->Next)
    {
      /*fprintf(stdout,"Maximizing utility of cluster %d\n",ClusterIndex);*/
      for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
      {
        /*fprintf(stdout,"Cluster %d attribute %d: Old mean %f stddev %f\n",ClusterIndex,AttributeIndex,BCluster->Means[AttributeIndex],BCluster->Stddevs[AttributeIndex]);*/
          BCluster->Means[AttributeIndex]=0;
        Probability = 0;
          for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)      
          {
            BCluster->Means[AttributeIndex] += (double)ClusteringGetAttributeValue(This,PointIndex,AttributeIndex) * BCluster->Members[PointIndex];
            /*fprintf(stdout,"Point %d with odds %f has value %f (running total %f)\n",PointIndex,BCluster->Members[PointIndex],(double)ClusteringGetAttributeValue(This,PointIndex,AttributeIndex),BCluster->Means[AttributeIndex]);*/
          Probability += BCluster->Members[PointIndex];
          }
        /*fprintf(stdout,"MeanTotal:%f\n",BCluster->Means[AttributeIndex]);*/
          BCluster->Means[AttributeIndex] /= Probability;
        /*fprintf(stdout,"Divide by %f to get %f\n",Probability,BCluster->Means[AttributeIndex]);*/
          BCluster->Stddevs[AttributeIndex]=0;
          Value = 0;
        Probability = 0;
          for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)      
          {
            Value = (double)ClusteringGetAttributeValue(This,PointIndex,AttributeIndex);
            Value = (Value - BCluster->Means[AttributeIndex]);    
            Value = Value*Value; /* C has no power operator for floats */
            BCluster->Stddevs[AttributeIndex] += BCluster->Members[PointIndex] * Value;
          Probability += BCluster->Members[PointIndex];
          }
          BCluster->Stddevs[AttributeIndex] /= Probability;
          BCluster->Stddevs[AttributeIndex] = sqrt(BCluster->Stddevs[AttributeIndex]);
        /*fprintf(stdout,"Cluster %d attribute %d: New mean %f stddev %f\n",ClusterIndex,AttributeIndex,BCluster->Means[AttributeIndex],BCluster->Stddevs[AttributeIndex]);*/
      }
      ClusterIndex++;
    }
    CycleCount++;
    fprintf(stdout,"  On cycle %d...\n",CycleCount);
    if (CycleCount>=10)
    {
      break;
    }
  }
  ClusteringFinalizeBayesianClusters(This,Head);
  fprintf(stdout,"  Free bayesian scaffolding:\n");
  BayesianClusterListFree(Head);
  fprintf(stdout,"  Bayesian clustering complete.\n");

}

CLUSTER_API int BASEClusteringGetAttributeCount(Clustering* This)
{
  return 1; /* ABSTRACT METHOD! */
}
CLUSTER_API float BASEClusteringGetAttributeMin(Clustering* This, int AttributeIndex)
{
  fprintf(stdout,"-**Abstract method!**-\n");
  return 0.0; /* ABSTRACT METHOD! */
}
CLUSTER_API float BASEClusteringGetAttributeMax(Clustering* This, int AttributeIndex)
{
  return 1.0; /* ABSTRACT METHOD! */
}
CLUSTER_API float BASEClusteringGetAttributeValue(Clustering* This, int PointIndex, int AttributeIndex)
{
  return 0.0; /* ABSTRACT METHOD! */
}

float BASEClusteringSOMGetDistance(Clustering* This, SOMNode* SOM,int PointIndex)
{
  return 1.0; /* ABSTRACT METHOD! */
}
void BASEClusteringSOMLearn(Clustering* This, SOMNode* SOM, float* LearningRate,int PointIndex)
{
  return; /* ABSTRACT METHOD! */
}

/* Random clustering */
CLUSTER_API void ClusteringClusterByChance(Clustering* This, int DesiredClusterCount)
{
    int PointIndex;
    int ClusterIndex;
	float RandValue;
    int RandInt;
    Cluster* NewCluster;
    ClusterNode* Node;
    
	Randomize();
    for (ClusterIndex=0; ClusterIndex<DesiredClusterCount; ClusterIndex++)
    {
        NewCluster = ClusteringGetNewCluster(This);
        ClusteringAddCluster(This,NewCluster);
    }
    for (PointIndex=0; PointIndex<This->PointCount; PointIndex++)
    {
        RandValue = rand()/(double)RAND_MAX;
        RandInt = RandValue * DesiredClusterCount;
        for (Node = This->Head,ClusterIndex = 0; ClusterIndex < RandInt; ClusterIndex++) {
        	Node = Node->pNext; 
        }
        ClusterAddMember(Node->Cluster, PointIndex);
        ClusterFindCentroid(Node->Cluster);
      	printf("point %d goes to cluster %d\n", PointIndex, ClusterIndex);        
    }

}

float ClusteringAverageDistancetoCentroid(Clustering* This,Cluster* TheCluster)
{
  float TotalDistance = 0;
  int PointCount = 0;
  int PointIndex;
  float Distance;
  
  for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
  {
    if (ClusterIsMember(TheCluster,PointIndex))
    {
      PointCount += 1;
      Distance = ClusteringDistanceToCentroid(This,PointIndex,TheCluster);
      TotalDistance += Distance;
    }
  }
  /*fprintf(stdout, "\n");*/
  return TotalDistance / PointCount;
}


/* ANOVA

*/
void ClusteringComputeANOVA(Clustering* This)
{
  float TotalDistanceToCentroid;
  ClusterNode* Node, *NodeB;
  int ClusterIndex, ClusterIndexB;
  float TotalClusterDistanceToCentroid;
  Cluster* TempCluster;
  int PointIndex;
  float Distance;
  float CDistance;
  float Num;
  float Den;
  /**/
  ClusteringFindAllCentroids(This); /* Make sure all centroids are up-to-date! */
  /* Form a cluster with ALL points, to get a centroid: */
  TempCluster = ClusteringGetNewCluster(This);
  ClusterAddEverything(TempCluster);
  ClusterFindCentroid(TempCluster);
  /* Iterate over all points.  Accumulate distance from the center, and from the cluster-center: */
  TotalDistanceToCentroid = 0;
  TotalClusterDistanceToCentroid = 0;
  for (PointIndex=0; PointIndex < This->PointCount; PointIndex++)
  {
    Distance = ClusteringDistanceToCentroid(This, PointIndex, TempCluster);
    Distance = Distance*Distance;
    /*fprintf(stdout, "Point %d distance to center %.4f\n", PointIndex, Distance);*/
    TotalDistanceToCentroid += Distance;
    CDistance = -1.0;
    for (Node = This->Head; Node; Node = Node->pNext)
    {
      if (ClusterIsMember(Node->Cluster, PointIndex))
      {
		CDistance = ClusteringDistanceToCentroid(This, PointIndex, Node->Cluster);
		CDistance = CDistance*CDistance;
		TotalClusterDistanceToCentroid += CDistance;
        Node->Cluster->SSEWithin += CDistance;
        /*fprintf(stdout, "Point %d distance to C-center %.4f\n", PointIndex, CDistance);*/
		break;
      }
    }
    if (CDistance < 0)
    {
      fprintf(stdout, "*** WARNING!  No cluster centroid for point %d\n", PointIndex);
    }
    /*fprintf(stdout, "Running: %.4f global, %.4f local\n", TotalDistanceToCentroid, TotalClusterDistanceToCentroid);*/
  }
  fprintf(stdout, "Pseudo-f: Total distance to centroid is %.4f\n",TotalDistanceToCentroid);
  fprintf(stdout, "Pseudo-f: Cluster distance to centroid is %.4f\n",TotalClusterDistanceToCentroid);
  This->SSE = TotalClusterDistanceToCentroid;
  This->SST = TotalDistanceToCentroid;
  Num = (TotalDistanceToCentroid - TotalClusterDistanceToCentroid) / (This->ClusterCount - 1);    
  Den = (TotalClusterDistanceToCentroid) / (This->PointCount - This->ClusterCount);
  fprintf(stdout, "Num %.4f over Den %.4f gives %.4f\n", Num, Den, (Num/Den));
  
  for (Node = This->Head,ClusterIndex=0; Node->pNext; Node = Node->pNext, ClusterIndex++)
  {
  
    for (NodeB = Node->pNext,ClusterIndexB=ClusterIndex+1; NodeB; NodeB = NodeB->pNext, ClusterIndexB++)
    {
      printf("Cluster %d, Cluster %d  \n", ClusterIndex, ClusterIndexB);
      ClusteringMeasureTScore(This, Node, NodeB);
    }
  }
}


CLUSTER_API void ClusteringBuildTotalCluster(Clustering* This) 
{
  Cluster* TempCluster;
  float TotalDistanceToCentroid;
  float Distance;
  int PointIndex;
  
  TempCluster = ClusteringGetNewCluster(This);
  ClusterAddEverything(TempCluster);
  ClusterFindCentroid(TempCluster);
  This->TotalCluster = TempCluster;

  TotalDistanceToCentroid = 0;
  for (PointIndex=0; PointIndex < This->PointCount; PointIndex++)
  {
    Distance = ClusteringDistanceToCentroid(This, PointIndex, TempCluster);
    Distance = Distance*Distance;
    /*fprintf(stdout, "Point %d distance to center %.4f\n", PointIndex, Distance);*/
    TotalDistanceToCentroid += Distance;
  }
  This->SST = TotalDistanceToCentroid;
}


float c2_a2b2(float *a, float *b, float *c, int n) 
{
  int i;
  float total = 0;
  for (i = 0; i < n; i++) {
    if (c[i] > 0)
      total += c[i]*c[i]/(a[i]*a[i]+b[i]*b[i]);
  }
  total /= n;
  return(total);
}

float c_ab(float *a, float *b, float *c, int n) 
{
  int i;
  float total = 0;
  
  for (i = 0; i < n; i++) {
    if (c[i] == 0)
      total += 0;
    else
      total += c[i]/(a[i]+b[i]);
  }
  total /= n;
  return(total);
}

float sum_c2_a2b2(float *a, float *b, float *c, int n) 
{
  int i;
  float c2,a2,b2;
  
  c2 = b2 = a2 = 0;
  for (i = 0; i < n; i++) {
    a2 += a[i]*a[i];
    b2 += b[i]*b[i];
    c2 += c[i]*c[i];
  }
  return(c2/(a2+b2));
}

float sum_c_ab(float *a, float *b, float *c, int n) 
{
  int i;
  float c2,a2,b2;
  
  c2 = b2 = a2 = 0;
  for (i = 0; i < n; i++) {
    a2 += a[i];
    b2 += b[i];
    c2 += c[i];
  }
  return(c2/(a2+b2));
}

float a2b2_c2_ab(float *a, float *b, float *c, int n) 
{
  int i;
  float total = 0;
  float temp;
  for (i = 0; i < n; i++) {
    temp = 0;
    if (a[i]*a[i]+b[i]*b[i]-c[i]*c[i] > 0)
      temp = (a[i]*a[i]+b[i]*b[i]-c[i]*c[i])/(2*a[i]*b[i]);
    total += temp;
    fprintf(stdout, "%f\t%f\t%f\t%f\n", a[i], b[i], c[i], temp);
  }
  total /= n;
  return(total);
}

float sum_a2b2_c2_ab(float *a, float *b, float *c, int n) 
{
  int i;
  float c2,a2,b2,ab,num;
  
  c2 = b2 = a2 = ab = 0;
  for (i = 0; i < n; i++) {
    a2 = a[i]*a[i];
    b2 = b[i]*b[i];
    c2 = c[i]*c[i];
    ab += a[i]*b[i];
    num += a2 + b2 -c2;
  }
  return(num/ab/2);
}
/*

*/
#define MAX_COMPARE_FILE 10
CLUSTER_API void ClusteringCompare(Clustering** Clusterings, char* FilePath) 
{
  int PointCount[MAX_COMPARE_FILE+1];
  ClusterNode* Node1;
  ClusterNode* Node2;
  ClusterNode* Nodes[MAX_COMPARE_FILE+1];
  int ClusterIndex1;
  int ClusterIndex2;
  int PointIndex;
  float Distance0, Distance1, Distance2, Distance;
  float **data;
  FILE* CompareFile;
  int FileIndex, FileCount;
  int i, j;
  Clustering *This;
  Clustering *That;
  float Distances[MAX_COMPARE_FILE+1];
  
  
  for (FileIndex = 0; Clusterings[FileIndex]; FileIndex++) {
    PointCount[FileIndex] = Clusterings[FileIndex]->PointCount;
    if (PointCount[FileIndex] != PointCount[0]) {
      fprintf(stderr, "ClusteringCompare():  Cannot compare two clusterings on different number of points!\n");
      return;
    }
  }
  FileCount = FileIndex;
  data = (float **)safe_malloc(sizeof(float*) * (1 + FileCount*(FileCount +1)/2));
  for (i=0; i<1 + FileCount*(FileCount +1)/2; i++) {
    data[i] = (float *)safe_malloc(sizeof(float) * PointCount[0]);
  }
  CompareFile = fopen(FilePath, "a");
  fprintf(CompareFile, "%10s\t%20s\t", "Point", "Distance2TotalCentroid");
  for (i = 0; i < FileCount; i++) {
    fprintf(CompareFile, "     Dist2Cluster%-3d\t", i);
  }
  for (i = 0; i < FileCount; i++) {
    for (j = i+1; j < FileCount; j++) {
      fprintf(CompareFile, "Dist2Cluster%3d-%-4d\t", i,j);
    }
  }
  fprintf(CompareFile, "\n");
  for (PointIndex=0; PointIndex < PointCount[0]; PointIndex++)
  {
    This = Clusterings[0];
    Distance = ClusteringDistanceToCentroid(This, PointIndex, This->TotalCluster);
    data[0][PointIndex] = Distance; /* Distance of a point to total centroid */
    fprintf(stdout, "PointIndex: %10d\t%20.5f\t", PointIndex, Distance);
    fprintf(CompareFile, "%10d\t%20.5f\t", PointIndex, Distance);
    for (i = 0; i < FileCount; i++) {
      This = Clusterings[i];
      for (Node1 = This->Head; Node1; Node1 = Node1->pNext)
        if (ClusterIsMember(Node1->Cluster, PointIndex))
          break;
 /*     Nodes[i] = Node1;*/
      Distances[i] = ClusteringDistanceToCentroid(This, PointIndex, Node1->Cluster);
      data[i+1][PointIndex] = Distance; /* Distance of a point to its cluster centroid */
      fprintf(stdout, "%20.5f\t", Distances[i]);
      fprintf(CompareFile, "%20.5f\t", Distances[i]);
    }
    for (i = 0; i < FileCount; i++) {
      This = Clusterings[i];
      for (Node1 = This->Head,ClusterIndex1=0; Node1; Node1 = Node1->pNext,ClusterIndex1++)
        if (ClusterIsMember(Node1->Cluster, PointIndex))
          break;
      for (j = i+1; j < FileCount; j++) { 
        That = Clusterings[j];
        for (Node2 = That->Head,ClusterIndex2=0; Node2; Node2 = Node2->pNext,ClusterIndex2++)
          if (ClusterIsMember(Node2->Cluster, PointIndex))
            break;
        /*
        Node1 = Nodes[i];
        Node2 = Nodes[j];
        */
        Distance = ClusteringCentroidToCentroidDistance(This, Node1, Node2);
        data[FileCount*(i+1)-(i+1)*i/2+j-i][PointIndex] = Distance; /* Distance of a point to its cluster centroid */
        fprintf(stdout, "%20.5f\t", Distance);
        fprintf(CompareFile, "%20.5f\t", Distance);
      }
    }
    fprintf(stdout, "\n");
    fprintf(CompareFile, "\n");
  }
  
  float *fn();
  static FloatVirtualFunction function[] = { c2_a2b2, c_ab, sum_c2_a2b2, sum_c_ab, a2b2_c2_ab, sum_a2b2_c2_ab};
  static char *functionname[] = { "avg(c2/(a2+b2))", "avg(c/(a+b))", "sum(c2)/sum(a2+b2)" , "sum(c)/sum(a+b)", "avg((a2+b2-c2)/ab)", "sum(a2+b2-c2)/sum(ab)"};
  int k;
  for (k = 0; k < 6; k++) {
    fprintf(CompareFile, "%20s", functionname[k]);
    for (i = 0; i< FileCount; i++) {
      fprintf(CompareFile, "      Comparefile%-3d ", i+1);
    }
    fprintf(CompareFile, "\n");
    for (i = 0; i < FileCount; i++) {
      fprintf(CompareFile, "      Comparefile%-3d ", i+1);
      for (j = 0; j < FileCount; j++) {
        char temp[100];
        if (j > i)
          sprintf(temp, "%20.5f", function[k](data[i+1],data[j+1],data[FileCount*(i+1)-(i+1)*i/2+j-i], PointCount[i]));
        else
          sprintf(temp, "%20s", " ");
        fprintf(CompareFile, "%s ", temp);
      }
      fprintf(CompareFile, "\n");
    }
  }  
  fclose(CompareFile);
  
}

/* The pseudo-F statistic is another measure of clustering goodness.  HIGH values are GOOD.
Generally, one selects a cluster-count that gives a peak in the pseudo-f statistic (or pSF, for short)
Formula: 
A/B, where A = (T - P)/(G-1), and B = P / (n-G)
Here n is the number of points, G is the number of clusters,
T is the total distance from the all-data centroid, and P is the sum (for all clusters) of the distances from the cluster centroid.
*/
#define DISTANCE2CLUSTER "Distance2Cluster"
CLUSTER_API float ClusteringComputePseudoF(Clustering* This)
{
  float TotalDistanceToCentroid;
  ClusterNode* Node;
  int ClusterIndex;
  float TotalClusterDistanceToCentroid;
  Cluster* TempCluster;
  int PointIndex;
  float Distance;
  float CDistance;
  float Num;
  float Den;
  float SSR = 0;
  ClusterNode* TempNode;
  FILE* DistanceFile = NULL;
  
  /**/
  ClusteringFindAllCentroids(This); /* Make sure all centroids are up-to-date! */
  /* Form a cluster with ALL points, to get a centroid: */
  
  if(This->SST == 0) {
    TempCluster = ClusteringGetNewCluster(This);
    ClusterAddEverything(TempCluster);
    ClusterFindCentroid(TempCluster);
    TotalDistanceToCentroid = 0;
  } else {
    /*TotalDistanceToCentroid = This->SST;*/
    TotalDistanceToCentroid = 0;
    TempCluster = This->TotalCluster;
  }
  /*
  TempNode = (ClusterNode*)malloc(sizeof(ClusterNode));
  memset(TempNode,0,sizeof(ClusterNode));
  TempNode->Cluster = TempCluster;
 
  DistanceFile = fopen(DISTANCE2CLUSTER,"w");
  fprintf(DistanceFile, "%10s\t%20s\t%20s\n", "Point", "Dist2TotalCentroid", "Dist2ClusterCentroid");
  */
  
  /* Iterate over all points.  Accumulate distance from the center, and from the cluster-center: */
  TotalClusterDistanceToCentroid = 0;
  for (Node = This->Head; Node; Node = Node->pNext)
  {
    Node->Cluster->SSEWithin = 0;
    /* SSR is the sum of square of regression. SSR + SSE = SST is true for regression. It is also true in this case.
    Distance = ClusteringCentroidToCentroidDistance(This, Node, TempNode);
    SSR +=  Distance * Distance * ClusterCountPointsInCluster(Node->Cluster);
    */
  }
  for (PointIndex=0; PointIndex < This->PointCount; PointIndex++)
  {
    if(This->SST >= 0) {
      Distance = ClusteringDistanceToCentroid(This, PointIndex, TempCluster);
      /*
      fprintf(DistanceFile, "%10d\t%20f\t", PointIndex, Distance);
      fprintf(stdout, "%10d\t%20f\t", PointIndex, Distance);
      */
      Distance = Distance*Distance;
      /*fprintf(stdout, "Point %d distance to center %.4f\n", PointIndex, Distance);*/
      TotalDistanceToCentroid += Distance;
    }
    CDistance = -1.0;
    for (Node = This->Head; Node; Node = Node->pNext)
    {
      if (ClusterIsMember(Node->Cluster, PointIndex))
      {
		CDistance = ClusteringDistanceToCentroid(This, PointIndex, Node->Cluster);
        /*
        fprintf(DistanceFile, "%20f\n", CDistance);
        fprintf(stdout, "%20f\n", CDistance);
        */
		CDistance = CDistance*CDistance;
        Node->Cluster->SSEWithin += CDistance;
		TotalClusterDistanceToCentroid += CDistance;
        /*fprintf(stdout, "Point %d distance to C-center %.4f\n", PointIndex, CDistance);*/
		break;
      }
    }
    if (CDistance < 0)
    {
      fprintf(stdout, "*** WARNING!  No cluster centroid for point %d\n", PointIndex);
    }
    /*fprintf(stdout, "Running: %.4f global, %.4f local\n", TotalDistanceToCentroid, TotalClusterDistanceToCentroid);*/
  }
  /*fprintf(stdout, "SST is %f, TotalDistanceToCentroid is %f\n", This->SST, TotalDistanceToCentroid);*/
  This->SST = TotalDistanceToCentroid;
  This->SSE = TotalClusterDistanceToCentroid;
  Num = (TotalDistanceToCentroid - TotalClusterDistanceToCentroid) / (This->ClusterCount - 1);    
  Den = (TotalClusterDistanceToCentroid) / (This->PointCount - This->ClusterCount);
  if(This->SST == 0) {
    This->SST = TotalDistanceToCentroid;
    fprintf(stdout, "Pseudo-f: Total distance to centroid is %.4f\n",TotalDistanceToCentroid);
    fprintf(stdout, "Pseudo-f: Cluster distance to centroid is %.4f\n",TotalClusterDistanceToCentroid);
    fprintf(stdout, "Num %.4f over Den %.4f gives %.4f\n", Num, Den, (Num/Den));
  }
  /*
  fprintf(stdout, "%f\t%f\t", TotalDistanceToCentroid, TotalClusterDistanceToCentroid);
  fprintf(stdout, "%f\t", TotalDistanceToCentroid - TotalClusterDistanceToCentroid);
  fprintf(stdout, "%f\n", SSR);
  free(TempNode);
  */
  /*fclose(DistanceFile);*/
  return (Num / Den);

}

/* The Davies-Bouldin Index (DBI) is a measure of clustering merit; the smaller the DBI,
the better.  The DBI is defined as the average, for all clusters X, of fred, where 
fred(X) = max, across other clusters Y, of (Cx + Cy)/dXY
...here Cx is the average distance from points in X to the centroid, similarly Cy, and dXY is 
the distance between cluster centroids
*/
CLUSTER_API float ClusteringComputeDBI(Clustering* This, FILE* OutputFile)
{
  ClusterNode* NodeX;
  int NodeIndexX;
  ClusterNode* NodeY;
  int NodeIndexY;
  float DBITotal = 0;
  float MaxFred = 0;
  float Fred = 0;
  float CX = 0;
  float CY = 0;
  float* AverageDistance;
  /**/
  ClusteringFindAllCentroids(This); /* Make sure all centroids are up-to-date! */
  AverageDistance = (float*)malloc(sizeof(float) * This->ClusterCount);
  memset(AverageDistance,0,sizeof(float) * This->ClusterCount);
  for (NodeX = This->Head,NodeIndexX=0; NodeX; NodeX = NodeX->pNext,NodeIndexX++)
    {
      AverageDistance[NodeIndexX] = ClusteringAverageDistancetoCentroid(This,NodeX->Cluster);
      if (OutputFile)
      {
        fprintf(OutputFile,"#Cluster %d has average-distance-to-centroid %f\n",NodeIndexX,AverageDistance[NodeIndexX]);
      }
    }
  for (NodeX = This->Head,NodeIndexX=0; NodeX; NodeX = NodeX->pNext,NodeIndexX++)
  {
    MaxFred = 0;
    for (NodeY = This->Head,NodeIndexY=0; NodeY; NodeY = NodeY->pNext,NodeIndexY++)
    {
      if (NodeX == NodeY)
      {
        continue;
      }
      Fred = AverageDistance[NodeIndexX] + AverageDistance[NodeIndexY];
      Fred /= ClusteringCentroidToCentroidDistance(This, NodeX, NodeY);
      /*fprintf(stdout,"From %d to %d, c2cdist is %f\n",NodeIndexX,NodeIndexY,ClusteringCentroidToCentroidDistance(This, NodeX, NodeY));*/
      if (Fred>MaxFred)
      {
          MaxFred = Fred;
      }
    }
    DBITotal += MaxFred;
  }
  free(AverageDistance);
  return DBITotal / This->ClusterCount;
}

#define SOM_ISOLATE		0
#define SOM_LOOP		1
#define SOM_STRING		2
#define SOM_BARBLOOP	3
#define SOM_BARBSTRING	4

void ClusteringSetSOMMap(Clustering* This,SOMNode** SOMNodes,int ClusterCount, int map)
{
  int SOMNodeIndex, Index;
  int prevIndex, nextIndex;

  if (map == SOM_LOOP) {
  	for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount; SOMNodeIndex++) {
    	SOMNodes[SOMNodeIndex]->NeighborCount = 2;
    	SOMNodes[SOMNodeIndex]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[SOMNodeIndex]->NeighborCount);
    	prevIndex = SOMNodeIndex - 1; 
        if (prevIndex == -1) prevIndex += ClusterCount;
        SOMNodes[SOMNodeIndex]->Neighbors[0] = SOMNodes[prevIndex];
        nextIndex = (SOMNodeIndex + 1)%ClusterCount;
        SOMNodes[SOMNodeIndex]->Neighbors[1] = SOMNodes[nextIndex];
    }
  } else
  if (map == SOM_STRING) {
  	for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount; SOMNodeIndex++) {
    	if (SOMNodeIndex == 0) {
	        SOMNodes[SOMNodeIndex]->NeighborCount = 1;
    		SOMNodes[SOMNodeIndex]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[SOMNodeIndex]->NeighborCount);
	        SOMNodes[SOMNodeIndex]->Neighbors[0] = SOMNodes[1];
        } else
    	if (SOMNodeIndex == ClusterCount - 1) {
	        SOMNodes[SOMNodeIndex]->NeighborCount = 1;
    		SOMNodes[SOMNodeIndex]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[SOMNodeIndex]->NeighborCount);
	        SOMNodes[SOMNodeIndex]->Neighbors[0] = SOMNodes[SOMNodeIndex - 1];
        } else {
	        SOMNodes[SOMNodeIndex]->NeighborCount = 2;
    		SOMNodes[SOMNodeIndex]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[SOMNodeIndex]->NeighborCount);
	    	prevIndex = SOMNodeIndex - 1; 
        	SOMNodes[SOMNodeIndex]->Neighbors[0] = SOMNodes[prevIndex];
	        nextIndex = SOMNodeIndex + 1;
    	    SOMNodes[SOMNodeIndex]->Neighbors[1] = SOMNodes[nextIndex];
        }
    }
  } else
  if (map == SOM_BARBLOOP) {
   	if (ClusterCount % 2 != 0 || ClusterCount < 4) {
       	fprintf(stdout, "ClusterCount has to be at least 4 and not an odd number to use barbloop. Now using ISOLATED clusters.\n");
        map = SOM_ISOLATE;
    } else {
	  	for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount / 2; SOMNodeIndex++) {
	        Index = SOMNodeIndex * 2;
            
            SOMNodes[Index]->NeighborCount = 3;
   			SOMNodes[Index]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[Index]->NeighborCount);
	    	prevIndex = Index - 2; 
    	    if (prevIndex < 0) prevIndex += ClusterCount;
        	SOMNodes[Index]->Neighbors[0] = SOMNodes[prevIndex];
        	SOMNodes[Index]->Neighbors[1] = SOMNodes[Index + 1];
    	    nextIndex = (Index + 2)%ClusterCount;
	        SOMNodes[Index]->Neighbors[2] = SOMNodes[nextIndex];
	        
            Index = Index + 1;
            SOMNodes[Index]->NeighborCount = 1;
   			SOMNodes[Index]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[Index]->NeighborCount);
    	    SOMNodes[Index]->Neighbors[0] = SOMNodes[Index - 1];
        }	
    }
  } else
  if (map == SOM_BARBSTRING) {
   	if (ClusterCount % 2 != 0 || ClusterCount < 4) {
       	fprintf(stdout, "ClusterCount has to be at least 4 and not an odd number to use barbstring. Now using ISOLATED clusters.\n");
        map = SOM_ISOLATE;
    } else {
	  	for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount / 2; SOMNodeIndex++) {
	        Index = SOMNodeIndex * 2;
            
    		if (Index == 0) {
		        SOMNodes[Index]->NeighborCount = 1;
    			SOMNodes[Index]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[Index]->NeighborCount);
	    	    SOMNodes[Index]->Neighbors[0] = SOMNodes[1];
            } else 
    		if (Index == ClusterCount - 2) {
		        SOMNodes[Index]->NeighborCount = 1;
    			SOMNodes[Index]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[Index]->NeighborCount);
	    	    SOMNodes[Index]->Neighbors[0] = SOMNodes[Index - 2];
            } else {
		        SOMNodes[SOMNodeIndex]->NeighborCount = 3;
    			SOMNodes[SOMNodeIndex]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[Index]->NeighborCount);
		    	nextIndex = Index - 2; 
		    	nextIndex = Index + 2; 
	        	SOMNodes[Index]->Neighbors[0] = SOMNodes[prevIndex];
    	    	SOMNodes[Index]->Neighbors[1] = SOMNodes[Index + 1];
    		    SOMNodes[Index]->Neighbors[2] = SOMNodes[nextIndex];
            }
            
            Index = Index + 1;
            SOMNodes[Index]->NeighborCount = 1;
   			SOMNodes[Index]->Neighbors = (SOMNode**)malloc(sizeof(SOMNode*) * SOMNodes[Index]->NeighborCount);
    	    SOMNodes[Index]->Neighbors[0] = SOMNodes[Index - 1];
        }	
    }
    	
  }

  if (map == SOM_ISOLATE) {
  	for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount; SOMNodeIndex++) {
    	SOMNodes[SOMNodeIndex]->NeighborCount = 0;
    	SOMNodes[SOMNodeIndex]->Neighbors = NULL;
    }
  } 
    

}

void ClusteringFreeSOMNodes(Clustering* This,SOMNode** SOMNodes,int ClusterCount)
{
  int SOMNodeIndex;
  for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount; SOMNodeIndex++)
    {
      safe_free(SOMNodes[SOMNodeIndex]->Attributes);
      safe_free(SOMNodes[SOMNodeIndex]->Neighbors);
      safe_free(SOMNodes[SOMNodeIndex]->Means);
      safe_free(SOMNodes[SOMNodeIndex]->Stddevs);
      safe_free(SOMNodes[SOMNodeIndex]);
    }
  safe_free(SOMNodes);
}

SOMNode** ClusteringInitSOMNodes(Clustering* This, int ClusterCount, int map)
{
  SOMNode** SOMNodes;
  int SOMNodeIndex;
  int AttributeCount;
  int AttributeIndex;
  float AttributeMin;  
  float AttributeMax;
  float RandValue;
  float Value;
  /* Allocate our SOMNodes, and initialize them: */
  SOMNodes = (SOMNode**)malloc(sizeof(SOMNode*) * ClusterCount);
  memset(SOMNodes,0,sizeof(SOMNode*) * ClusterCount);
  AttributeCount = ClusteringGetAttributeCount(This);
  for (SOMNodeIndex=0; SOMNodeIndex < ClusterCount; SOMNodeIndex++)
  {
    SOMNodes[SOMNodeIndex] = (SOMNode*)malloc(sizeof(SOMNode));
    memset(SOMNodes[SOMNodeIndex],0,sizeof(SOMNode));
    SOMNodes[SOMNodeIndex]->Attributes = (float*)malloc(sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Attributes,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->Means = (float*)malloc(sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Means,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->Stddevs = (float*)malloc(sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Stddevs,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->PointCount = 0;
  }
  for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
  {
    AttributeMin = ClusteringGetAttributeMin(This,AttributeIndex);
    AttributeMax = ClusteringGetAttributeMax(This,AttributeIndex);
    for (SOMNodeIndex=0; SOMNodeIndex < ClusterCount; SOMNodeIndex++)
    {
      RandValue = rand()/(double)RAND_MAX;
      Value = AttributeMax - AttributeMin;
      SOMNodes[SOMNodeIndex]->Attributes[AttributeIndex] = AttributeMin + Value*RandValue;
    }
  }
  ClusteringSetSOMMap(This, SOMNodes, ClusterCount, map); /* Setup map */
  return SOMNodes;
}

/* Helper function for SOM clustering: Given a point and some SOMNodes, find the 'winnner' 
   (that is, the SOMNode which is closest to the point). */
int ClusteringSOMFindWinner(Clustering* This, SOMNode** SOMNodes, int ClusterCount, int PointIndex)
{
  float WinnerDistance;
  int SOMNodeIndex;
  int WinnerIndex;
  float Distance;
  WinnerDistance = -1;
  for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount; SOMNodeIndex++)
  {
     Distance = ClusteringSOMGetDistance(This,SOMNodes[SOMNodeIndex],PointIndex);
     if (WinnerDistance < 0 || Distance < WinnerDistance)
       {
         WinnerDistance = Distance;
         WinnerIndex = SOMNodeIndex;
       }
  }
  return WinnerIndex;
}


/*
For processing points in random order: Choose a random point.  And, if it has already been
processed, scan forward (wrapping at the end of the list) until you find an unprocessed one.
(Not totally random, since point n+1 is more likely to be processed after point n has been
hit, but this greatly reduces our dependence on training-set-order)
 */
/*
int ChooseNextPoint(int* PointProcessed,int PointCount,int Dummy)
{
  int PointIndex;
  float RandValue;
  RandValue = rand()/(double)RAND_MAX;
  PointIndex = floor(RandValue * PointCount);
  while (1)
    {
      if (!PointProcessed[PointIndex])
        {
          PointProcessed[PointIndex] = 1;
          return PointIndex;
        }
      PointIndex = (PointIndex+1)%PointCount;
    }
}
*/
/*
For processing points in random order: Choose a random point.  And, if it has already been
processed, scan forward (wrapping at the end of the list) until you find an unprocessed one.
(Now is truly random, but the caller should garrantee the RemainingPointCount is correct.)
 */
int ChooseNextPoint(int* PointProcessed,int PointCount,int RemainingPointCount)
{
  int PointIndex;
  int PointChosen;
  float RandValue;
  
  RandValue = rand()/(double)RAND_MAX;
  PointChosen = floor(RandValue * RemainingPointCount);
  PointIndex = 0;
  while (PointIndex < PointCount)
    {
      if (!PointProcessed[PointIndex])
        {
	      if (PointChosen == 0)
    	    {
        	  PointProcessed[PointIndex] = 1;
	          return PointIndex;
    	    }
          PointChosen--;
         }
     PointIndex = PointIndex+1;
    }
  return -1; 
}

/*
Self-organizing maps are a learning methodology (a type of artificial neural network) 
pioneered by Kohonen et al for speech processing.  The "map" consists of a collection of
Nodes, each of which has some Neighbors.  (In our case, we use a list of Nodes, where 
each node has a previous and a next neighbor; SOMs can also be laid out as a grid or 
a higher-dimensional collection of arbitrarily zany topology)  Each SOMNode has a 
value similar to the points that we are learning about - for instance, if we're 
clustering points in the plane, each SOMNode will be a point in the plane.  The 
network 'learns' through several exposures to the training set.  For each training point 
we present, we find the point in the map which is most similar (the "winner"), and then 
make that point and its neighbors even MORE similar to the training point.  Over time, 
the system comes to "recognize" the training points.  We convert this to cluster output 
in the natural way: Each SOMNode corresponds to a cluster, and each cluster gets, as
members, every point for which its SOMNode was the winner.  (Another option would be
to use more SOMNodes that our desired cluster-count, and then cluster the SOMNodes after
training, but it's not obvious how that problem is simpler than clustering the original
points, or how the results would be better)
*/
CLUSTER_API void ClusteringClusterSOM(Clustering* This, int ClusterCount, int map)
{
  SOMNode** SOMNodes;
  float WinnerLearningRate;
  float NeighborLearningRate;
  int IterationIndex;
  int AttributeIndex;
  int AttributeCount;
  float AttributeMin;
  float AttributeMax;
  float RandValue;
  int WinnerIndex;
  float WinnerDistance;
  float Distance;
  int NeighborIndex;
  int PointIndex;
  int* PointProcessed;
  int ProcessIndex;
  /**/
  Randomize();
  SOMNodes = ClusteringInitSOMNodes(This,ClusterCount,map);
  PointProcessed = (int*)malloc(sizeof(int) * This->PointCount);
  /* Iterate! */ 
  WinnerLearningRate = 0.1;
  NeighborLearningRate = 0.05;
  for (IterationIndex = 0; IterationIndex < 1000; IterationIndex++)
  {
    /* Process the training points in random order, so that we're not sensitive to the 
       data ordering */
    memset(PointProcessed,0,sizeof(int) * This->PointCount);
    for (ProcessIndex = 0; ProcessIndex < This->PointCount; ProcessIndex++)
      {
        PointIndex = ChooseNextPoint(PointProcessed,This->PointCount,This->PointCount-ProcessIndex);
        if (PointIndex < 0) continue;
        /* For each point, find the "winner" (the SOMNode closest to it).  Adjust the
           attributes of the winner and its neighbor, to be closer to the winner */
        WinnerIndex = ClusteringSOMFindWinner(This,SOMNodes, ClusterCount, PointIndex);
        /*fprintf(stdout,"For point %d, the closest SOM is %d; learning rate %.5f:%.5f\n",PointIndex,WinnerIndex,WinnerLearningRate,NeighborLearningRate);*/
        /* We have a winnah!  Now move it, and its neighbors, closer to the point: */
        ClusteringSOMLearn(This,SOMNodes[WinnerIndex],&WinnerLearningRate,PointIndex);
        for (NeighborIndex = 0; NeighborIndex < SOMNodes[WinnerIndex]->NeighborCount; NeighborIndex++)
          {
            ClusteringSOMLearn(This,SOMNodes[WinnerIndex]->Neighbors[NeighborIndex],&NeighborLearningRate,PointIndex);
          }
        /*
        NeighborIndex = (WinnerIndex+1)%ClusterCount;
        ClusteringSOMLearn(This,SOMNodes[NeighborIndex],&NeighborLearningRate,PointIndex);
        NeighborIndex = (WinnerIndex-1);
        if (NeighborIndex < 0)
          {
            NeighborIndex = ClusterCount-1;
          }
        ClusteringSOMLearn(This,SOMNodes[NeighborIndex],&NeighborLearningRate,PointIndex);
        */
      }
    /* Decelerate your motion: */
    if (IterationIndex%10 == 0)
      {
        WinnerLearningRate *= 0.95;
        NeighborLearningRate *= 0.95;
      }
  }
  /* Use our SOMNode data to generate clusters */
  ClusteringConvertSOMNodesToClusters(This,SOMNodes,ClusterCount);
  ClusteringFreeSOMNodes(This,SOMNodes,ClusterCount);
}

/* The default output-status method is to do nothing.  */
void BASEClusteringOutputStatus(Clustering* This)
{
  
}

  
/*****************************************************************
Clustering VTable stuff 
*****************************************************************/
void ClusteringFree(Clustering* This)
{
  This->VirtualFunctions[F_CLUSTERING_FREE](This);
}
void ClusteringMergeClusters(Clustering* This,ClusterNode* NodeA,ClusterNode* NodeB)
{
  This->VirtualFunctions[F_CLUSTERING_MERGE_CLUSTERS](This,NodeA,NodeB);
}
float ClusteringDistanceToCentroid(Clustering* This,int PointIndex, Cluster* Cluster)
{
  return (float)(This->FloatVirtualFunctions[F_CLUSTERING_DISTANCE_TO_CENTROID](This,PointIndex,Cluster));
}
void ClusteringSplitCluster(Clustering* This,Cluster* ParentCluster, int PointA, int PointB)
{
  This->VirtualFunctions[F_CLUSTERING_SPLIT_CLUSTER](This,ParentCluster,PointA,PointB);
}
float ClusteringPointToPointDistance(Clustering* This, int A, int B)
{
  return (float)(This->FloatVirtualFunctions[F_CLUSTERING_POINT_TO_POINT](This, A, B));
}
float ClusteringCentroidToCentroidDistance(Clustering* This, ClusterNode* NodeA,ClusterNode* NodeB)
{
  return (float)(This->FloatVirtualFunctions[F_CLUSTERING_CENTROID_TO_CENTROID](This,NodeA,NodeB));
}
Cluster* ClusteringGetNewCluster(Clustering* This)
{
  return (Cluster*)This->VirtualFunctions[F_CLUSTERING_GET_NEW_CLUSTER](This);
}
int* ClusteringFindKmeansSeeds(Clustering* This, int Seeds)
{
  return (int*)This->VirtualFunctions[F_CLUSTERING_FIND_KMEAN_SEEDS](This,Seeds);
}
void ClusteringClusterLinkage(Clustering* This, int DesiredClusterCount,float Epsilon)
{
  This->VirtualFunctions[F_CLUSTERING_CLUSTER_LINKAGE](This,DesiredClusterCount,Epsilon);
}
void ClusteringClusterMeans(Clustering* This,int* SeedPoints, int DesiredClusterCount, int interation, int mode)
{
  This->VirtualFunctions[F_CLUSTERING_CLUSTER_MEANS](This,SeedPoints,DesiredClusterCount, interation, mode);
}

void ClusteringClusterCobweb(Clustering* This, int DesiredClusterCount)
{
  This->VirtualFunctions[F_CLUSTERING_CLUSTER_COBWEB](This,DesiredClusterCount);
}

CobwebCluster* ClusteringNewCobwebCluster(Clustering* This, int PointIndex)
{
  return This->VirtualFunctions[F_CLUSTERING_NEW_COBWEB_CLUSTER](This,PointIndex);
}

void ClusteringCobwebAddNewFrame(Clustering* This,CobwebTreeNode* Root,CobwebTreeNode* Frame,int level)
{
  This->VirtualFunctions[F_CLUSTERING_COBWEB_ADD_FRAME](This,Root,Frame,level);
}

void CobwebClusterComputeMeans(Clustering* This, CobwebTreeNode* Node)
{
  This->VirtualFunctions[F_CLUSTERING_COBWEB_COMPUTE_MEANS](This,Node);
}

float ClusteringCobwebCategoryUtility(Clustering * This, CobwebTreeNode * Root)
{
  return This->FloatVirtualFunctions[F_CLUSTERING_COBWEB_CATEGORY_UTILITY](This, Root);
}

void ClusteringCobwebClusterFree(Clustering* This, CobwebCluster* DeadGuy)
{
  This->VirtualFunctions[F_CLUSTERING_COBWEB_CLUSTER_FREE](This,DeadGuy);
}
void CobwebClusterCopyMeans(Clustering* This, CobwebCluster* Source, CobwebCluster* Destination)
{
  This->VirtualFunctions[F_CLUSTERING_COBWEB_COPY_MEANS](This,Source,Destination);
}

CLUSTER_API int ClusteringGetAttributeCount(Clustering* This)
{
  return This->IntVirtualFunctions[F_CLUSTERING_GET_ATTRIBUTE_COUNT](This);
}
CLUSTER_API float ClusteringGetAttributeMin(Clustering* This, int AttributeIndex)
{
  return This->FloatVirtualFunctions[F_CLUSTERING_GET_ATTRIBUTE_MIN](This,AttributeIndex);
}
CLUSTER_API float ClusteringGetAttributeMax(Clustering* This, int AttributeIndex)
{
  return This->FloatVirtualFunctions[F_CLUSTERING_GET_ATTRIBUTE_MAX](This,AttributeIndex);
}
float ClusteringGetAttributeValue(Clustering* This, int PointIndex, int AttributeIndex)
{
  return This->FloatVirtualFunctions[F_CLUSTERING_GET_ATTRIBUTE_VALUE](This,PointIndex,AttributeIndex);
}
float ClusteringSOMGetDistance(Clustering* This, SOMNode* SOM,int PointIndex)
{
  return This->FloatVirtualFunctions[F_CLUSTERING_SOM_GET_DISTANCE](This,SOM, PointIndex);
}
void ClusteringSOMLearn(Clustering* This, SOMNode* SOM,float* LearningRate, int PointIndex)
{
  This->VirtualFunctions[F_CLUSTERING_SOM_LEARN](This,SOM, LearningRate,PointIndex);
}
void ClusteringOutputStatus(Clustering* This)
{
  This->VirtualFunctions[F_CLUSTERING_OUTPUT_STATUS](This);
}

/*****************************************************************
General-purpose functions live here - these are not "methods" of 
a Clustering or Cluster object
*****************************************************************/
void MeasureDistances(int argc, char* argv[], int UseChiSquared)
     /* int ClusteringCount,int Points, char* FileNames[])*/
{
    Clustering* TheClustering;
    Clustering** Clusterings;
    int IndexA;
    int IndexB;
    SymmetricMatrix* PairwiseDistances;
    SymmetricMatrix* C2CDistances;
    int ClusteringCount = argc-2; /* skip prog-name and "dist" command */
    int Points;
    FILE* SourceFile;
    float Distance;
    /**/
    SourceFile = fopen(argv[2],"r");
    Points = GetLineLength(SourceFile);
    fclose(SourceFile);
    fprintf(stdout,"Point count %d\n",Points);
    Clusterings = (Clustering**)malloc(sizeof(Clustering*) * ClusteringCount);
    C2CDistances = AllocateSymmetricMatrix(ClusteringCount);
    PairwiseDistances = AllocateSymmetricMatrix(Points);
    for (IndexA=0; IndexA<ClusteringCount; IndexA++)
    {
        TheClustering= ClusteringNew(PairwiseDistances);
        ClusteringReadFromDisk(TheClustering,argv[IndexA+2]);
        Clusterings[IndexA]=TheClustering;
    }
    for (IndexA=0; IndexA<ClusteringCount; IndexA++)
    {
        for (IndexB=IndexA+1; IndexB<ClusteringCount; IndexB++)
        {
            fprintf(stdout,"\n%s vs %s:\n",argv[IndexA+2],argv[IndexB+2]);
	    if (UseChiSquared)
            {
	      /* The distance is the p-value that the two clusterings are unrelated */
	      Distance = ClusteringMeasureChiSquared(Clusterings[IndexA],Clusterings[IndexB]);
	    }
	    else
	    {
	      Distance = ClusteringMeasureDistance(Clusterings[IndexA],Clusterings[IndexB]);
	    }
            SetSymElement(C2CDistances,IndexA,IndexB,Distance);
        }
    } 
    WriteSymMatrix("C2CDistances.txt",C2CDistances, 1);
    FreeSymmetricMatrix(PairwiseDistances);
    FreeSymmetricMatrix(C2CDistances);   
}

/* return the covariance of two arrays, ith row and jth row of Matrix. 
   If Cookies is provided and n is the length of the Cookies, calculate covariance for those elements appeared in the Cookies. 
   If Cookies == NULL, n should be the rank of the matrix.
*/
float MeasureCovariance(SymmetricMatrix* Matrix, int i, int j, int* Cookies, int n) 
{
	float sumi = 0;
    float sumj = 0;
    float sumij = 0;
    
    int k, ck;
    for (k = 0; k < n; k++) {
    	if (Cookies) ck = Cookies[k];
        else ck = k;
        sumi += GetSymElement(Matrix, i, ck);
        sumj += GetSymElement(Matrix, j, ck);
        sumij += GetSymElement(Matrix, i, ck) * GetSymElement(Matrix, j, ck);
    }
    return (sumij - sumi*sumj/n)/n;
}


/*
 *  cluster matrix utilities
 */


CLUSTER_API Matrix* AllocateMatrix(int Size)
{
    Matrix* pMatrix;
    
    if (Size >= 23171) { /* sqrt(2^31/sizeof(float))*/
    	fprintf(stderr, "Error in AllocateSymmetricMatrix, Size (%d) of SymmetricMatrix is too big. Please limit the size to be less than 23171.\n", Size);
        exit(1);
    }
    pMatrix = (Matrix*)SafeMalloc(__FILE__, __LINE__, sizeof(Matrix));
    pMatrix->Size = Size;
    pMatrix->Data = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * Size * Size);
    memset(pMatrix->Data,0,sizeof(float)*Size*Size);
    return pMatrix;
}

CLUSTER_API void FreeMatrix(Matrix* pMatrix)
{
    safe_free(pMatrix->Data);
    safe_free(pMatrix);
    return;
}

CLUSTER_API SymmetricMatrix* AllocateSymmetricMatrix(int Size)
{
    SymmetricMatrix* pMatrix;
    /* */
    if (Size >= 32768) {
    	fprintf(stderr, "Error in AllocateSymmetricMatrix, Size (%d) of SymmetricMatrix is too big. Please limit the size to be less than 32768.\n", Size);
        exit(1);
    }
    unsigned int DataBlockSizeInBytes;
    pMatrix = (SymmetricMatrix*)SafeMalloc(__FILE__, __LINE__, sizeof(Matrix));
    pMatrix->Size = Size;
    DataBlockSizeInBytes = (sizeof(float) * Size * (Size+1)) / 2;
    pMatrix->Data = (float*)SafeMalloc(__FILE__, __LINE__, DataBlockSizeInBytes);
    memset(pMatrix->Data,0,DataBlockSizeInBytes);    
    return pMatrix;
}

CLUSTER_API DoubleSymmetricMatrix* AllocateDoubleSymmetricMatrix(int Size)
{
    DoubleSymmetricMatrix* pMatrix;
    /* */
    if (Size >= 23170) {
    	fprintf(stderr, "Error in AllocateSymmetricMatrix, Size (%d) of SymmetricMatrix is too big. Please limit the size to be less than 23170.\n", Size);
        exit(1);
    }
    unsigned int DataBlockSizeInBytes;
    pMatrix = (DoubleSymmetricMatrix*)SafeMalloc(__FILE__, __LINE__, sizeof(Matrix));
    pMatrix->Size = Size;
    DataBlockSizeInBytes = (sizeof(double) * Size * (Size+1)) / 2;
    pMatrix->Data = (double*)SafeMalloc(__FILE__, __LINE__, DataBlockSizeInBytes);
    memset(pMatrix->Data,0,DataBlockSizeInBytes);    
    return pMatrix;
}

CLUSTER_API void FreeSymmetricMatrix(SymmetricMatrix* pMatrix)
{
    safe_free(pMatrix->Data);
    safe_free(pMatrix);
    return;
}

CLUSTER_API void FreeDoubleSymmetricMatrix(DoubleSymmetricMatrix* pMatrix)
{
    safe_free(pMatrix->Data);
    safe_free(pMatrix);
    return;
}

CLUSTER_API float GetSymElement(SymmetricMatrix* pMatrix, int Row, int Column)
{
    int Temp;
    if (Row > Column)
    {
        Temp = Row;
        Row = Column;
        Column = Temp;
    }
    /*return pMatrix->Data[Row * pMatrix->Size + (Row - Row*Row)/2 + Column];*/
    return SymmetricMatrixElement(pMatrix,Row,Column);
}

CLUSTER_API double GetDoubleSymElement(DoubleSymmetricMatrix* pMatrix, int Row, int Column)
{
    int Temp;
    if (Row > Column)
    {
        Temp = Row;
        Row = Column;
        Column = Temp;
    }
    /*return pMatrix->Data[Row * pMatrix->Size + (Row - Row*Row)/2 + Column];*/
    return SymmetricMatrixElement(pMatrix,Row,Column);
}

CLUSTER_API void SetSymElement(SymmetricMatrix* pMatrix, int Row, int Column, float Value)
{
    int Temp;
    if (Row > Column)
    {
        Temp = Row;
        Row = Column;
        Column = Temp;
    }
    /*return pMatrix->Data[Row * pMatrix->Size + (Row - Row*Row)/2 + Column];*/
    SymmetricMatrixElement(pMatrix,Row,Column) = Value;
}

CLUSTER_API void SetDoubleSymElement(DoubleSymmetricMatrix* pMatrix, int Row, int Column, double Value)
{
    int Temp;
    if (Row > Column)
    {
        Temp = Row;
        Row = Column;
        Column = Temp;
    }
    /*return pMatrix->Data[Row * pMatrix->Size + (Row - Row*Row)/2 + Column];*/
    SymmetricMatrixElement(pMatrix,Row,Column) = Value;
}

/*
PrintMatrix prints out a matrix to stdout, for debugging purposes
*/
void PrintMatrix(Matrix* pMatrix)
{
	int Row;
	int Column;
	for (Row=0; Row<pMatrix->Size; Row++)
	{
        printf("[");
		for (Column=0; Column<pMatrix->Size; Column++)
		{
			printf("%f ",MatrixElement(pMatrix,Row,Column));
		}
		printf("]\n");
	}
	printf("\n");
}


void PrintMatrixToFile(Matrix* pMatrix, char* FileName)
{
	int Row;
	int Column;
    FILE* MatrixFile = fopen(FileName,"w");
    /* */
	for (Row=0; Row<pMatrix->Size; Row++)
	{
        fprintf(MatrixFile,"[");
		for (Column=0; Column<pMatrix->Size; Column++)
		{
            if (Column)
            {
                fprintf(MatrixFile,", ");
            }
            /* Break lines, for a more readable file: */
            if (Column%10 == 0)
            {
                fprintf(MatrixFile,"\n  ");
            }
			fprintf(MatrixFile,"%f",MatrixElement(pMatrix,Row,Column));
		}
		fprintf(MatrixFile,"\n]\n");
	}
	fprintf(MatrixFile,"\n");
    fclose(MatrixFile);
}


float DeterminantRecursive(Matrix* pMatrix,int* Rows,int* Columns, int Depth)
{
    int RowIndex;
    int ColumnIndex;
    int Sign = 1;
    float Sum = 0;

    for (RowIndex=0; RowIndex<pMatrix->Size; RowIndex++)
    {
        if (Rows[RowIndex])
        {
            continue;
        }
        Rows[RowIndex]=1;
        for (ColumnIndex=0; ColumnIndex < pMatrix->Size; ColumnIndex++)
        {            
            if (!Columns[ColumnIndex])
            {
                if (Depth==pMatrix->Size)
                {
                    Rows[RowIndex]=0;
                    return MatrixElement(pMatrix,RowIndex,ColumnIndex);
                }
                Columns[ColumnIndex]=1;
                Sum += Sign * MatrixElement(pMatrix,RowIndex,ColumnIndex) * DeterminantRecursive(pMatrix,Rows,Columns,Depth+1);
                Columns[ColumnIndex]=0;
                Sign *= -1;
            }
        }
        Rows[RowIndex]=0;
        break;
    }
    return Sum;
}

/**
Algorithm for calculating the determinant:
Proceed recursively, as follows:
Compute a sum by iterating over the elements of the topmost row.  For the element
in column ColumnIndex, add to your sum the matrix element, multiplied by Sign 
(-1 to the nth power), multipled by the determinant of the submatrix formed by removing
that row and column.  Use the Row and Column arrays to keep track of which rows and 
columns have been "removed" - a value of 1 indicates that the section has already
been removed in the current submatrix.
*/
float Determinant(Matrix* pMatrix)
{
    int* Rows;
    int* Columns;
    float Result;
    Rows = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * pMatrix->Size);
    memset(Rows,0,sizeof(int) * pMatrix->Size);
    Columns = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * pMatrix->Size);
    memset(Columns,0,sizeof(int) * pMatrix->Size);
    Result = DeterminantRecursive(pMatrix,Rows,Columns,1);
    safe_free(Rows);
    safe_free(Columns);
    return Result;
}

/**
TriangularizeMatrix takes a square matrix with a given size and changes it to lower triangular form.
It returns the number of row-switches performed on the matrix.
*/
int TriangularizeMatrix(Matrix* pMatrix)
{
    int Column;
	int Row;
    int ColumnTemp;
	float Swap;
	float Scaling;
    int RowSwitches = 0;
    int ProgressCount = 0;
    /* */
    ProgressCount = pMatrix->Size / 50;
	/*PrintMatrix(pMatrix);	*/
	for (Column=0; Column<pMatrix->Size-1; Column++)
	{		
        if (Column%ProgressCount == 0)
        {
            printf(".");
        }
		if (MatrixElement(pMatrix,Column,Column)==0)
		{
			/* Interchange this row with the first lower that has a nonzero element in this column*/
			for (Row=Column+1; Row<pMatrix->Size; Row++)
			{
				if (MatrixElement(pMatrix,Row,Column)!=0)
				{
					for (ColumnTemp=0; ColumnTemp<pMatrix->Size; ColumnTemp++)
					{
						Swap=MatrixElement(pMatrix,Row,ColumnTemp);
						MatrixElement(pMatrix,Row,ColumnTemp)=MatrixElement(pMatrix,Column,ColumnTemp);
						MatrixElement(pMatrix,Column,ColumnTemp)=Swap;
					}
					break;
				}
			}
		}
		if (MatrixElement(pMatrix,Column,Column)==0)
		{
			continue;
		}
		/* Null out each column in lower rows: */
		for (Row=Column+1; Row<pMatrix->Size; Row++)
		{			
			Scaling = MatrixElement(pMatrix,Row,Column) / MatrixElement(pMatrix,Column,Column);
			for (ColumnTemp=0; ColumnTemp<pMatrix->Size; ColumnTemp++)
			{
				MatrixElement(pMatrix,Row,ColumnTemp) -= MatrixElement(pMatrix,Column,ColumnTemp) * Scaling;
			}
		}
	}
    return RowSwitches;
}    


void PrintSymMatrix(SymmetricMatrix* Matrix)
{
    FILE* TheFile;
    int Row;
    int Column;
    float Value;
    /**/
    for (Row=0; Row<Matrix->Size; Row++)
    {
        for (Column=0; Column<Matrix->Size; Column++)
        {
		  Value = GetSymElement(Matrix,Row,Column);
	      fprintf(stdout,"%6.2f ",Value);
	    }
	    fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
}

void WriteSymMatrix(char* FilePath,SymmetricMatrix* Matrix, int TextMode)
{
    FILE* TheFile;
    int Row;
    int Column;
    float Value;
    /**/
    TheFile = fopen(FilePath,"wb");
    if (!TheFile)
    {
    	fprintf(stdout, "Can not write SymMatrix file %s\n", FilePath);
        return;
    }
    for (Row=0; Row<Matrix->Size; Row++)
    {
        for (Column=0; Column<Matrix->Size; Column++)
        {
		  Value = GetSymElement(Matrix,Row,Column);
		  if (TextMode)
		  {
		      fprintf(TheFile,"%.8f\t",Value);
		  }
		  else
		  {
		    fwrite(&Value,sizeof(float),1,TheFile);
		    /*fwrite(&SymmetricMatrixElement(Matrix,Row,Column),sizeof(float),1,TheFile);*/
		  }
	    }
	  if (TextMode)
	  {
	    fprintf(TheFile,"\n");
	  }
    }
    fclose(TheFile);
}

int ReadSymMatrix(char* FilePath,SymmetricMatrix* Matrix)
{
    FILE* TheFile;
    int Row;
    int Column;
    float Point;
    /**/
    TheFile = fopen(FilePath,"rb");
    for (Row=0; Row<Matrix->Size; Row++)
    {
        for (Column=0; Column<Matrix->Size; Column++)
        {
            if (!fread(&Point,sizeof(float),1,TheFile)) 
            {
            	warning("ReadSymMatrix()", "Existing file %s is not correct. Need more data.\n", FilePath);
                return(0);
            }
            SetSymElement(Matrix, Row, Column, Point);
            /*SymmetricMatrixElement(Matrix,Row,Column) = Point; Wrong!*/
        }
    }
    if (fread(&Point,sizeof(float),1,TheFile)) 
    {
     	warning("ReadSymMatrix()", "Existing file %s is not correct. It has more data.\n", FilePath);
        return(0);
    }
    fclose(TheFile);
    return(1);
}


/*
Test of matrix code.  No prerequisites.
*/
int UnitTestMatrices()
{
    Matrix* pMatrix;
    SymmetricMatrix* pSymmetric;
    float TestValue;
    int Index;
    /**/
    pMatrix = AllocateMatrix(3);

    /* Construct a simple matrix: */
    MatrixElement(pMatrix,0,0) = 1;
    MatrixElement(pMatrix,0,1) = 2;
    MatrixElement(pMatrix,0,2) = 3;
    MatrixElement(pMatrix,1,0) = 4;
    MatrixElement(pMatrix,1,1) = 5;
    MatrixElement(pMatrix,1,2) = 6;
    MatrixElement(pMatrix,2,0) = 7;
    MatrixElement(pMatrix,2,1) = 8;
    MatrixElement(pMatrix,2,2) = 10;
    /* Verify that we can calculate the determinant: */
    TestValue = Determinant(pMatrix);
    if (TestValue + 3.0 > 0.0001)
    {
        printf("Failed to get determinant!");
        return 0;
    }
    /* Triangularize and take determinant: */
    TriangularizeMatrix(pMatrix);
    TestValue = 1;
    for (Index=0; Index<pMatrix->Size; Index++)
    {
        TestValue = TestValue * MatrixElement(pMatrix,Index,Index);
    }
    if (TestValue + 3.0 > 0.0001)
    {
        printf("Bad triangularization!");
        return 0;
    }
    /* Construct and test a symmetric matrix */
    pSymmetric = AllocateSymmetricMatrix(4);
    SetSymElement(pSymmetric, 0, 0, 2.0);
    SetSymElement(pSymmetric, 0, 1, 3.0);
    SetSymElement(pSymmetric, 3, 3, 8.0);
    TestValue = GetSymElement(pSymmetric,1,0);
    if (TestValue - 3.0 > 0.0001)
    {
        printf("Broken symmetric matrix");
        return 0;
    }
    TestValue = GetSymElement(pSymmetric,3,3);
    if (TestValue - 8.0 > 0.0001)
    {
        printf("Broken symmetric matrix");
        return 0;
    }

    return 1;

}


CLUSTER_API SymmetricMatrix* SymmetricMatrixCopy(SymmetricMatrix* Original)
{
    SymmetricMatrix* NewMatrix;
    size_t DataBlockSizeInBytes = (sizeof(float) * Original->Size * (Original->Size+1)) / 2;
    /**/
    NewMatrix = AllocateSymmetricMatrix(Original->Size);
    memcpy(NewMatrix->Data,Original->Data,DataBlockSizeInBytes);
    return NewMatrix;
}
