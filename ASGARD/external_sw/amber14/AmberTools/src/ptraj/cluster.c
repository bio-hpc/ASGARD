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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/cluster.c,v 10.2 2009/08/17 21:40:06 sbrozell Exp $
 *
 *  Revision: $Revision: 10.2 $
 *  Date: $Date: 2009/08/17 21:40:06 $
 *  Last checked in by $Author: sbrozell $
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

#define CLUSTER_MODULE
#include "ptraj.h"




  /*
   *  Forward declaration 
   */
   void
OutputNewickTree(PtrajClustering* This, FILE* OutputFile);

   void
PtrajClusterTest(PtrajCluster* This)
{
  fprintf(stdout,"Cluster test - Ptraj override");
}

  /*
   *  PtrajClusteringGetNewCluster is reloaded for ClusteringGetNewCluster,
   *  which is defined in cluster.c and requires a return type Cluster*.
   */
   Cluster*
PtrajClusteringGetNewCluster(PtrajClustering* This)
{
  return (Cluster *) PtrajClusterNew(This);
}

  /*
   *  MDS is not finished yet. 
   */
   float
MDSDistance(PtrajClustering* This, int PointIndex, int OtherPointIndex)
{
  int AttributeIndex;
  int AttributeCount;
  int IsTor;
  float Distance = 0;
  float AttDistance;
  
  if (PointIndex < 0 || OtherPointIndex < 0 ) {
  	fprintf(stderr, "Error in MDSDistance, (Other)PointIndex is not valid.\n");
    return 0.0;
  }
  AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
  for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
  {
    IsTor = PtrajClusteringIsAttributeTorsion(This, AttributeIndex);
    if (!IsTor) {
      AttDistance = ClusteringGetAttributeValue( (Clustering *) This,PointIndex,AttributeIndex) 
	- ClusteringGetAttributeValue( (Clustering *) This,OtherPointIndex,AttributeIndex);
      Distance += (AttDistance*AttDistance);
    } else {
      AttDistance = abs(abs(ClusteringGetAttributeValue( (Clustering *) This,PointIndex,AttributeIndex)-180) 
			- ClusteringGetAttributeValue( (Clustering *) This,OtherPointIndex,AttributeIndex));
      Distance += (AttDistance*AttDistance);
    }
  }
  return sqrt(Distance/AttributeCount);
  
}

float DMEDistance(float* X, float* Y, float* Z, float* OtherX, float* OtherY, 
		  float* OtherZ, int AtomCount)
{
  int AtomIndex;
  int OtherAtomIndex;
  float AtomDistance;
  float OtherAtomDistance;
  float dX, dY, dZ;
  float dXX, dYY, dZZ;
  float PairDistance;
  float Distance = 0;
  
  for (AtomIndex=0; AtomIndex<AtomCount; AtomIndex++)
  {
    for (OtherAtomIndex=AtomIndex+1; OtherAtomIndex<AtomCount; OtherAtomIndex++)
    {
      dX = X[AtomIndex] - X[OtherAtomIndex];
      dY = Y[AtomIndex] - Y[OtherAtomIndex];
      dZ = Z[AtomIndex] - Z[OtherAtomIndex];
      dXX = OtherX[AtomIndex] - OtherX[OtherAtomIndex];
      dYY = OtherY[AtomIndex] - OtherY[OtherAtomIndex];
      dZZ = OtherZ[AtomIndex] - OtherZ[OtherAtomIndex];
      
      /*
      AtomDistance = sqrt(dX*dX + dY*dY + dZ*dZ);
      OtherAtomDistance = sqrt(dXX*dXX + dYY*dYY + dZZ*dZZ);
      PairDistance = (AtomDistance-OtherAtomDistance);
      Distance += PairDistance*PairDistance;
      
      printf("Old    AtomDistance is %f, OtherAtomDistance is %f, PairDistance is %f\n", AtomDistance, OtherAtomDistance, PairDistance);
      */
      
      /* The above section can be optimized as following, however, the AtomDistance, OtherAtomDistance, and PairDistance now has been changed to squared. -- Jianyin */
      
      AtomDistance = dX*dX + dY*dY + dZ*dZ;
      OtherAtomDistance = dXX*dXX + dYY*dYY + dZZ*dZZ;
      PairDistance = AtomDistance+OtherAtomDistance-2*sqrt(AtomDistance*OtherAtomDistance);
      Distance += PairDistance;
      
    }
  }
  Distance = Distance * 2 / (AtomCount*(AtomCount-1));
  return sqrt(Distance);	  
}


float PtrajGetDistance(PtrajClustering* This, int Frame, int Other, float* X, float* Y, float* Z, float* OtherX, 
		       float* OtherY, float* OtherZ)
{
    trajectoryInfo* trajInfo = This->trajInfo;
    int mass = This->action->iarg6;
    float rmsR[3][3];
    float rmsT[3];

    if (This->DistanceMetric == DISTANCE_METRIC_RMSD)
    {
    	if (mass)
        {
        	return (float)rmsf(trajInfo->atoms,1,trajInfo->state->masses,NULL,X,Y,Z,
                  OtherX,OtherY,OtherZ,rmsR,rmsT,0);   
         } else
         {
        	return (float)rmsf(trajInfo->atoms,1,NULL,NULL,X,Y,Z,
                  OtherX,OtherY,OtherZ,rmsR,rmsT,0);   
         }
    }
    else if (This->DistanceMetric == DISTANCE_METRIC_DME)
    {
    	return (float)DMEDistance(X,Y,Z,OtherX, OtherY, OtherZ, trajInfo->atoms);
    }
    else if (This->DistanceMetric == DISTANCE_METRIC_MDS)
    {
    	return (float)MDSDistance(This, Frame, Other);
    }

    return 0.0;
}


/*  PairwiseDistance is a matrix whose (i,j)-th element is the DISTANCE between frame i and frame j.
    The distance between two frames is the root-mean-squared distance between them.
*/
SymmetricMatrix* ComputePairwiseDistance(PtrajClustering* This, trajectoryInfo* trajInfo,  actionInformation* action, char *CacheDistanceFile)
{
    int FrameIndex;
    int OtherFrameIndex;
    int AtomIndex;
    int OtherAtomIndex;
    float* FrameX;
    float* FrameY;
    float* FrameZ;
    float* OtherFrameX;
    float* OtherFrameY;
    float* OtherFrameZ;
    float Distance;
    float AtomDistance;
    float OtherAtomDistance;
    SymmetricMatrix* PairwiseDistance;
    int FrameCount;
    FILE* DistanceFile = NULL;
    int SymMatrixReady = 0;

    FrameCount = action->iarg4;

    PairwiseDistance =  AllocateSymmetricMatrix(FrameCount);

    if (!SymMatrixReady)
    {
      fprintf(stdout, "  Calculating the PairwiseDistances\n");
      for (FrameIndex=0; FrameIndex<FrameCount; FrameIndex++)
   	  {
	  	FrameX = trajInfo->x + trajInfo->atoms*FrameIndex;
	  	FrameY = trajInfo->y + trajInfo->atoms*FrameIndex;
	  	FrameZ = trajInfo->z + trajInfo->atoms*FrameIndex;
	  	for (OtherFrameIndex=FrameIndex+1; OtherFrameIndex<FrameCount; OtherFrameIndex++)
	    {
	      OtherFrameX = trajInfo->x + trajInfo->atoms*OtherFrameIndex;
	      OtherFrameY = trajInfo->y + trajInfo->atoms*OtherFrameIndex;
	      OtherFrameZ = trajInfo->z + trajInfo->atoms*OtherFrameIndex;
	      
	      Distance=PtrajGetDistance(This, FrameIndex, OtherFrameIndex, FrameX, FrameY, FrameZ, OtherFrameX, OtherFrameY,
					OtherFrameZ);
	      SymmetricMatrixElement(PairwiseDistance,FrameIndex,OtherFrameIndex)=Distance;                
	    }
	  }  
      fprintf(stdout, "  Writing Pairwise Distances to file \"PairwiseDistances\".\n");
      WriteSymMatrix(CacheDistanceFile, PairwiseDistance, 0);
    }
    This->PairwiseDistances = PairwiseDistance;
    This->PointCount = PairwiseDistance->Size;
    return PairwiseDistance;
}

void SymmetricMatrixDistribution(SymmetricMatrix* matrix, int partition, FILE *File) {
	int i, j, l;
    int size = matrix->Size;
    float max, min, cur;
    float length;
    int* lot;
    int pairs;
    
    max = SymmetricMatrixElement(matrix,0,1);
    min = max;
    for (i=0; i<size-1; i++) {
     	for (j = i+1; j<size; j++) {
        	cur = SymmetricMatrixElement(matrix,i,j);
            if (cur > max)
            	max = cur;
            else if (cur < min)
            	min = cur;
        }
    }
    
    lot = (int *)SafeMalloc(__FILE__, __LINE__, partition * sizeof(int));
    for (i=0; i<partition; i++) 
    	lot[i] = 0;
    length = (max - min) / partition;
    for (i=0; i<size-1; i++) {
     	for (j = i+1; j<size; j++) {
        	cur = SymmetricMatrixElement(matrix,i,j);
        	l = (int)((cur - min) / length);
            if (l == partition) l--;
            lot[l]++;
        }
    }
    
    pairs = size*(size-1)/2;
    fprintf(File, "###Distribution of Distances###\n");
    fprintf(stdout, "    Distribution of Distances\n");
    for (i=0; i<partition-1; i++)  {
    	fprintf(stdout, "  [%8.3f,%8.3f] -- %6.2f%% (%6d out of %6d)\n", min+i*length, min+(i+1)*length, (float)lot[i]/pairs*100, lot[i], pairs);
    	fprintf(File, "#  [%8.3f,%8.3f] -- %6.2f%% (%6d out of %6d)\n", min+i*length, min+(i+1)*length, (float)lot[i]/pairs*100, lot[i], pairs);
        
    }
    fprintf(stdout, "  [%8.3f,%8.3f] -- %6.2f%% (%6d out of %6d)\n", min+i*length, min+(i+1)*length, (float)lot[i]/pairs*100, lot[i], pairs);
    fprintf(File, "#  [%8.3f,%8.3f] -- %6.2f%% (%6d out of %6d)\n", min+i*length, min+(i+1)*length, (float)lot[i]/pairs*100, lot[i], pairs);
}

#define BLOCK_SIZE 20480
/* Read a clustering from a text-file.  Re-reading the clustering won't 
   re-read the point definitions; we only get the membership info.  
   Useful for comparing the results of clustering algorithms. */
void PtrajClusteringReadFromDisk(Clustering* This,char* SourceFilePath)  
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
      fprintf(stdout,"Warning!  Line length %d doesn't match point count %d for file %s\n",LineSize,PointCount,SourceFilePath);
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
              PtrajClusterAddMember( (PtrajCluster *) NewCluster,PointIndex);
            }
            if (*Temp != 'X' && *Temp != '.')
            {
              fprintf(stdout,"Warning!  Bogus character %c encountered in cluster output file %s\n",*Temp,SourceFilePath);
            }
            PointIndex++;
        }
        if (PointIndex!=1 && PointIndex!=PointCount)
        {
          fprintf(stdout,"Warning!  Line length %d instead of %d found in cluster output file %s\n",PointIndex,PointCount,SourceFilePath);
        }
	/*fprintf(stdout,"Adding one cluster from disk.\n");*/
        ClusteringAddCluster(This,NewCluster);
    }
    fclose(SourceFile);
    free(FileLine);
}

void PtrajClusteringReadFromMergeFile(PtrajClustering* This,char* SourceFilePath)  
{
    char* FileLine;
    FILE* SourceFile;
    char* Result;
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    SymmetricMatrix* ClusterDistances;
    int PointCount = This->PointCount;
    int PointIndex;
    PtrajCluster* NewCluster;
    int NewNodeIndex, NodeIndex, NodeIndexA, NodeIndexB;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    ClusterNode* Node;
    float ClosestLength;
    int DesiredClusterCount = PtrajClusteringGetDesiredClusterCount(This);
    float Epsilon = This->action->darg1;
    int finish = 0;
    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    /*ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);*/

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;
    }
    
    SourceFile = fopen(SourceFilePath,"r");
    FileLine = (char*)SafeMalloc(__FILE__, __LINE__, sizeof(char)*(BLOCK_SIZE));
    while (!finish)
    {
        Result = fgets(FileLine, BLOCK_SIZE ,SourceFile);
        NodeA = NodeB = NULL;
        if (This->ClusterCount<=DesiredClusterCount+1)
        {
            finish = 1;
        }
        if (!DesiredClusterCount && ClosestLength > Epsilon)
        {
            finish = 1;
        }
        if (!Result)
        {
            printf("Error in ClusterMerging.txt, insufficient lines.\n");
            finish = 1;
        }
        sscanf(Result, "%d:\t%d\t%d\t%f", &NewNodeIndex, &NodeIndexA, &NodeIndexB, &ClosestLength);
        printf("Read from file        %i:\t%i\t%i\t%f\n", NewNodeIndex, NodeIndexA, NodeIndexB, ClosestLength);
        for (Node = This->Head,NodeIndex=0; Node && ((!NodeB) || (!NodeA)); Node=Node->pNext,NodeIndex++)
        {
        	if (!NodeB) { 
            	if (((PtrajCluster *)Node->Cluster)->IntName == NodeIndexB) {
            		NodeB = Node;
            	}
            }
            if (!NodeA) {
            	if (((PtrajCluster *)Node->Cluster)->IntName == NodeIndexA) {
            		NodeA = Node;
            	}
            }
		}
        if (!NodeA) {
        	printf("Error: Cannot find NodeA for NodeIndexA is %d\n", NodeIndexA); 
        }
        if (!NodeB) {
        	printf("Error: Cannot find NodeA for NodeIndexB is %d\n", NodeIndexB); 
        }
        
	ClusteringMergeClusters( (Clustering *) This, (ClusterNode *) NodeA, (ClusterNode *) NodeB);
        ((PtrajCluster *)NodeB->Cluster)->IntName = NewNodeIndex;        
	}
}

int PtrajClusteringGetAtomCount(PtrajClustering* This)
{
    return This->trajInfo->atoms;
}

int PtrajClusteringGetFrameCount(PtrajClustering* This)
{
    return This->PointCount; 
}

int PtrajClusteringGetDesiredClusterCount(PtrajClustering* This)
{
    return This->action->iarg2;
}
float PtrajClusteringGetDesiredEpsilon(PtrajClustering* This)
{
    return (float)This->action->darg1;
}

/*
Computes the centroid of each atom in This cluster.  The trajectory info
supplies the atom coordinates for each frame.
*/
void PtrajClusterFindCentroid(PtrajCluster* ThisCluster)
{
    int AtomIndex;
    int FrameIndex;
    int ClusteredFrameCount = 0;
    int FrameCount = ((PtrajClustering*)ThisCluster->Owner)->PointCount;
    int AtomCount  = PtrajClusteringGetAtomCount(ThisCluster->Owner);
    int ClusterIndex;
    ClusterNode* Temp;
    FILE* DebugFile;
    char FileName[1024];
    trajectoryInfo* trajInfo = ThisCluster->Owner->trajInfo;
    /* fprintf(stdout,"Current is %d\n",trajInfo->current); */
    for (Temp = ThisCluster->Owner->Head,ClusterIndex=0; Temp; Temp=Temp->pNext,ClusterIndex++)
    {
	  if (Temp->Cluster == (Cluster *) ThisCluster)
	  {
	    break;
	  }
    }
    /*fprintf(stdout,"Computing a centroid for cluster %d\n",ClusterIndex);*/
    /* Initialize centroid coordinates to 0 */
    memset(ThisCluster->CentroidX,0,sizeof(float)*AtomCount);
    memset(ThisCluster->CentroidY,0,sizeof(float)*AtomCount);
    memset(ThisCluster->CentroidZ,0,sizeof(float)*AtomCount);
    /* For each atom: Add up the x,y, and z coordinates of all frames in the cluster */
    for (FrameIndex=0; FrameIndex<FrameCount; FrameIndex++)
    {
        if (ClusterIsMember( (Cluster *) ThisCluster,FrameIndex))
        {
            ClusteredFrameCount++;
            for (AtomIndex=0; AtomIndex<AtomCount; AtomIndex++)
            {
                ThisCluster->CentroidX[AtomIndex]+=trajInfo->x[FrameIndex*trajInfo->atoms + AtomIndex];
                ThisCluster->CentroidY[AtomIndex]+=trajInfo->y[FrameIndex*trajInfo->atoms + AtomIndex];
                ThisCluster->CentroidZ[AtomIndex]+=trajInfo->z[FrameIndex*trajInfo->atoms + AtomIndex];
            }
        }
    }

    /* For each atom: Divide by the number of points in the cluster, to get the average coordinate */
    for (AtomIndex=0; AtomIndex<AtomCount; AtomIndex++)
    {
        ThisCluster->CentroidX[AtomIndex] = ThisCluster->CentroidX[AtomIndex] / ClusteredFrameCount;
        ThisCluster->CentroidY[AtomIndex] = ThisCluster->CentroidY[AtomIndex] / ClusteredFrameCount;
        ThisCluster->CentroidZ[AtomIndex] = ThisCluster->CentroidZ[AtomIndex] / ClusteredFrameCount;
    }
    /* Debugging: Dump the centroid coordinates! */
    /*sprintf(FileName,"Centroid%d",ClusterIndex);
    fprintf(stdout,"File name is %s\n",FileName);
    DebugFile = fopen(FileName,"w");
    for (AtomIndex=0; AtomIndex<AtomCount; AtomIndex++)
    {
      fprintf(DebugFile,"%.4f\t%.4f\t%.4f\n",ThisCluster->CentroidX[AtomIndex],ThisCluster->CentroidY[AtomIndex],ThisCluster->CentroidZ[AtomIndex]);
    }
    fclose(DebugFile);*/

}

PtrajCluster* PtrajClusterNew(PtrajClustering* Owner)
{
  PtrajCluster* This;
  int AtomCount = PtrajClusteringGetAtomCount(Owner);
  This = (PtrajCluster*)SafeMalloc(__FILE__, __LINE__, sizeof(PtrajCluster));
  memset(This,0,sizeof(PtrajCluster));
  This->PointCount = ClusteringGetPointCount( (Clustering *) Owner);
  This->SSEWithin = 0;
  This->Mask = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * This->PointCount);
  This->Name = NULL;
  This->IntName = -1;
  memset(This->Mask,0,sizeof(int) * This->PointCount);
  This->Owner = Owner;
  This->CentroidX = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AtomCount);
  memset(This->CentroidX, 0, sizeof(float) * AtomCount);
  This->CentroidY = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AtomCount);
  memset(This->CentroidY, 0, sizeof(float) * AtomCount);
  This->CentroidZ = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AtomCount);
  memset(This->CentroidZ, 0, sizeof(float) * AtomCount);
  This->Head = NULL;
  This->Tail = NULL;
  This->VirtualFunctions = PtrajClusterVirtualFunctions;
  This->P2C = (float *)SafeMalloc(__FILE__, __LINE__, sizeof(float) * This->PointCount);
  memset(This->P2C,0,sizeof(float) * This->PointCount);
  This->BestRep = -1;
  This->OldBestRep = -1;
  return This;
}


float ClusteringAlignAToB(PtrajClustering* This, int A, int B)
{
    float* XA;
    float* YA;
    float* ZA;
    float* Xref;
    float* Yref;
    float* Zref;
    float rmsRotation[3][3], rmsTranslation[3];
    trajectoryInfo* trajInfo = This->trajInfo;
	int* mask;
    int i;
    float xtemp, ytemp, ztemp;
    int mass = This->action->iarg6;
	int verbose = This->action->darg2;
    float rms;
        
    if (B <0) 
    {
    	warning("ClusteringAlignAToB()", "Point B is not defined");
        return (-1);
    }
    if (A <0) 
    {
    	warning("ClusteringAlignAToB()", "Point A is not defined");
        return (-1);
    }
    
    XA = &trajInfo->x[A * trajInfo->atoms]; 
    YA = &trajInfo->y[A * trajInfo->atoms]; 
    ZA = &trajInfo->z[A * trajInfo->atoms]; 
	
    Xref = &trajInfo->x[B * trajInfo->atoms]; 
    Yref = &trajInfo->y[B * trajInfo->atoms]; 
    Zref = &trajInfo->z[B * trajInfo->atoms]; 
    if (mass)
    {
    	rms = rmsf(trajInfo->atoms,1,trajInfo->state->masses,NULL,XA,YA,ZA,Xref,Yref,Zref,rmsRotation,rmsTranslation,2);
    } else
    {
    	rms = rmsf(trajInfo->atoms,1,NULL,NULL,XA,YA,ZA,Xref,Yref,Zref,rmsRotation,rmsTranslation,2);
    }
    if (verbose >= VERBOSE_ALIGN_FRAME) fprintf(stdout, "Align frame %i to frame %i,\n", A, B);
    
    return(rms);
} 

void PtrajClusterSetBestRep(PtrajCluster* This, int BestRep) 
{
	This->BestRep = BestRep;
}

int PtrajClusterGetBestRep(PtrajCluster* This) 
{
	return(This->BestRep);
}

void ClusterAddMemberP2C(PtrajCluster* This, int Point)
{
	int PointIndex;
    float Distance;
    
    if ( ClusterIsMember( (Cluster *) This, Point ) )
    {
    	warning("ClusterAddMemeberP2C()", "Point %d is already in the Cluster %d", Point, PtrajOutputIntermediateClusterStatus(This, stdout, 1));
        return;   /* Please make sure to call this function before call ClusterAddMember. */
    }
    for (PointIndex=0; PointIndex < ClusteringGetPointCount( (Clustering *) This->Owner); PointIndex++)
    {
        if (ClusterIsMember( (Cluster *) This,PointIndex) && 
	    (PointIndex != Point)) /* Skip Point, although the Distance should be 0. */
        {
            Distance = PtrajClusteringPointToPointDistance(This->Owner, Point, PointIndex);
            Distance = Distance * Distance;
            This->P2C[Point] += Distance;
            This->P2C[PointIndex] += Distance;
        }
    }
}

void ClusterSubtractMemberP2C(PtrajCluster* This, int Point)
{
	int PointIndex;
    float Distance;
    
    if (!ClusterIsMember( (Cluster *) This,Point ))
    {
    	warning("ClusterSubstractMemeberP2C()", "Point %d is not in the Cluster %d", Point, PtrajOutputIntermediateClusterStatus(This, stdout, 1));
    	return;   /* Please make sure to call this function before call ClusterRemoveMember. */
    }
    for (PointIndex=0; PointIndex < ClusteringGetPointCount( (Clustering *) This->Owner); PointIndex++)
    {
        if (ClusterIsMember( (Cluster *) This,PointIndex))
        {
            Distance = PtrajClusteringPointToPointDistance(This->Owner, Point, PointIndex);
            Distance = Distance * Distance;
            This->P2C[Point] -= Distance;
            This->P2C[PointIndex] -= Distance;
        }
    }
    This->P2C[Point] = 0;
}

void PtrajClusterAddMember(PtrajCluster* This, int Point)
{
	int OldBestRep = PtrajClusterGetBestRep(This);
    int NewBestRep;
    int C2CCentroid = This->Owner->action->darg4;
    
    ClusterAddMemberP2C(This, Point);
    ClusterAddMember( (Cluster *) This, Point);
    if (OldBestRep < 0) /* if BestRep is -1, the first member in the cluster. */
    {
    	PtrajClusterSetBestRep(This, Point);
        return;
    }
    if (C2CCentroid)
    {
	    ClusteringAlignAToB(This->Owner, Point, This->BestRep);
    }
}

void PtrajClusterRemoveMember(PtrajCluster* This, int Point)
{
  ClusterSubtractMemberP2C(This, Point);
  ClusterRemoveMember( (Cluster *) This, Point);
}


int ClusterFindBestP2C(PtrajCluster* This)
{
	int PointIndex;
	int BestRep = -1;
    float BestP2C = -100;
    
    for (PointIndex=0; PointIndex < ClusteringGetPointCount( (Clustering *) This->Owner); PointIndex++)
    {
        if (!ClusterIsMember( (Cluster *) This, PointIndex ))
        {
            continue;
        }
        if (This->P2C[PointIndex] < BestP2C || BestP2C < -99) {
        	BestP2C = This->P2C[PointIndex];
            BestRep = PointIndex;
        }
    }
    if (BestRep >= 0) {
    	This->BestRep = BestRep;
    }
    return(This->BestRep);
}


void ClusterAlignToBestRep(PtrajCluster* This)
{
	int PointIndex;
    int BestRep;
    float BestP2C = -100;
    float* XA;
    float* YA;
    float* ZA;
    float* XB;
    float* YB;
    float* ZB;
    trajectoryInfo* trajInfo = This->Owner->trajInfo;
    
    if (This->BestRep < 0) 
    {
    	warning("ClusterAlignToBestRep()", "BestRep is not defined in this cluster");
        return;
    }
    for (PointIndex=0; PointIndex < ClusteringGetPointCount( (Clustering *) This->Owner); PointIndex++)
    {
        if (!ClusterIsMember( (Cluster *) This,PointIndex) || PointIndex == This->BestRep)
        {
            continue;
        }
        if (This->P2C[PointIndex] < BestP2C || BestP2C < -99)
        {
        	ClusteringAlignAToB(This->Owner, PointIndex, This->BestRep);
        }
    }
}


void ClusteringMergeNames(PtrajClustering* This, ClusterNode* MergeNodeA, ClusterNode* MergeNodeB);
/* Merge two cluster - also merges the representatives */
void PtrajClusteringMergeClusters(PtrajClustering* This,ClusterNode* NodeA,ClusterNode* NodeB)
{
    int PointIndex;
    int OldBestRepA, OldBestRepB, NewBestRepB;
    int C2Ccentroid = This->action->darg4;
    PtrajCluster* ClusterA;
    PtrajCluster* ClusterB;
    /**/
    ClusteringMergeNames(This, NodeA, NodeB);
    ClusterA = (PtrajCluster*)NodeA->Cluster;
    ClusterB = (PtrajCluster*)NodeB->Cluster;
    for (PointIndex=0; PointIndex < ClusteringGetPointCount( (Clustering *) This); PointIndex++)
    {
        if (ClusterIsMember( (Cluster *) ClusterA,PointIndex))
        {
            PtrajClusterAddMember(ClusterB,PointIndex);
        }
    }
    
    OldBestRepA = PtrajClusterGetBestRep(ClusterA); 
    OldBestRepB = PtrajClusterGetBestRep(ClusterB); 
    NewBestRepB = ClusterFindBestP2C(ClusterB);
    /* If Distance to Cluster is measured to its BestRep, no need to rms-fit right now. 
       Rms-fit will be done at the end of the total clustering. */
    if (C2Ccentroid > 0)
    {
    	if (OldBestRepB != NewBestRepB)
	    {
			ClusterAlignToBestRep(ClusterB);
    	}
        ClusterFindCentroid( (Cluster *) ClusterB); 
    }
    /* Connecting the Rep list used in centripetal algorithm. */
    if (ClusterB->Tail && ClusterA->Head)
    {
        ClusterB->Tail->pNext = ClusterA->Head;
        ClusterA->Head->pPrev = ClusterB->Tail;
        ClusterB->Tail = ClusterA->Tail;
    }
    ClusteringRemoveCluster( (Clustering *) This, (Cluster *) ClusterA);
}


void PtrajClusterFree(PtrajCluster* This)
{
    safe_free(This->CentroidX);
    safe_free(This->CentroidY);
    safe_free(This->CentroidZ);
    safe_free(This->AverageX);
    safe_free(This->AverageY);
    safe_free(This->AverageZ);
    safe_free(This->RepX);
    safe_free(This->RepY);
    safe_free(This->RepZ);
    safe_free(This->entry);
    if (This->Name)
    {
      safe_free(This->Name);
    }
    safe_free(This->Mask);
    safe_free(This->P2C);
    safe_free(This);
}

void PtrajClusteringFree(PtrajClustering* This)
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
    safe_free(This->SieveFirstPass);
    safe_free(This);
}




PtrajClustering* PtrajClusteringNew(SymmetricMatrix* PairwiseDistances, trajectoryInfo* trajInfo, 
    actionInformation* action)
{
  PtrajClustering* This;
  This = (PtrajClustering*)SafeMalloc(__FILE__, __LINE__, sizeof(PtrajClustering));
  memset(This,0,sizeof(PtrajClustering));
  This->Head=NULL;
  This->Tail=NULL;
  This->PairwiseDistances = PairwiseDistances;
  This->VirtualFunctions = PtrajClusteringVirtualFunctions;
  This->FloatVirtualFunctions = PtrajClusteringFloatVirtualFunctions;
  This->IntVirtualFunctions = PtrajClusteringIntVirtualFunctions;
  This->trajInfo = trajInfo;
  This->action = action;
  This->PointCount = action->iarg4;
  This->SSE = 0;
  This->SST = 0;
  This->Acuity = 0.1; /* Used in incremental (cobweb) clustering */
  This->DistanceMetric = action->iarg5;
  This->StartTime = time(NULL);
  This->attributeArray = NULL;
  This->attributeArrayTorsion = NULL;
  return This;
}

float PtrajClusteringDistanceToCentroid(PtrajClustering* This,int PointIndex, PtrajCluster* Cluster)
{
    float* X;
    float* Y;
    float* Z;
    /**/
    X = This->trajInfo->x + PointIndex*This->trajInfo->atoms;
    Y = This->trajInfo->y + PointIndex*This->trajInfo->atoms;
    Z = This->trajInfo->z + PointIndex*This->trajInfo->atoms;

    return PtrajGetDistance(This, PointIndex, -1, X, Y, Z, Cluster->CentroidX, Cluster->CentroidY, Cluster->CentroidZ);
}

float PtrajClusteringDistanceToCluster(PtrajClustering* This,int PointIndex, PtrajCluster* Cluster)
{
    float* X;
    float* Y;
    float* Z;
    float* XB;
    float* YB;
    float* ZB;
  	int C2Ccentroid =  This->action->darg4; 
    /* Distance to Cluster is measured by its best representative. */
    if (Cluster->BestRep < 0) 
    {
       	ClusterFindBestP2C(Cluster);
    }
    if (C2Ccentroid == 0) 
    {
        if (This->action->performSecondPass) 
        {
		    X = This->trajInfo->x + PointIndex*This->trajInfo->atoms;
		    Y = This->trajInfo->y + PointIndex*This->trajInfo->atoms;
		    Z = This->trajInfo->z + PointIndex*This->trajInfo->atoms;
            XB = This->trajInfo->x + Cluster->BestRep*This->trajInfo->atoms;
		    YB = This->trajInfo->y + Cluster->BestRep*This->trajInfo->atoms;
		    ZB = This->trajInfo->z + Cluster->BestRep*This->trajInfo->atoms;
		    return PtrajGetDistance(This, PointIndex, -1, X, Y, Z, XB, YB,ZB);
        }
        else 
        {
        	return 	PtrajClusteringPointToPointDistance(This,PointIndex,Cluster->BestRep);
        }
    }
    else /* Distance to Cluster is measured by its centroid. */
    {
	    return PtrajClusteringDistanceToCentroid(This, PointIndex, Cluster);
    }
}

int PtrajOutputIntermediateClusterStatus(PtrajCluster* This, FILE* OutFile, int IndexOnly) 
{
    int ClusterIndex = 0;
    int PointIndex;
    PtrajClustering* Clustering = This->Owner;
    ClusterNode* Node;
    
    for (Node=Clustering->Head; Node; Node = Node->pNext, ClusterIndex++)
    {
        if (This != (PtrajCluster *) Node->Cluster) 
        {
            continue;
        }
        if (IndexOnly)
        { 
        	return(ClusterIndex);
        }
        fprintf(OutFile,"Cluster %4i: ", ClusterIndex);
        for (PointIndex=0; PointIndex<This->PointCount; PointIndex++)
        {
            if (ClusterIsMember( (Cluster *) This, PointIndex ))
            {
                fprintf(OutFile,"X");
            }
            else
            {
                fprintf(OutFile,".");
            }
        }
        fprintf(OutFile,"      With BestRep is %d\n", This->BestRep);
        return(ClusterIndex);
    }
    return (-1);
}

/*
Split the cluster, using multiple iterations to get a good new centroid
*/

#define ClusterSplittingFile "ClusterSplitting.txt"
/* Perform HIERARCHICAL clustering */
void PtrajClusteringClusterHierarchical(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    float EccentricEpsilon;
    PtrajCluster* EccentricCluster;
    int EccentricPointA;
    int EccentricPointB;
    PtrajCluster* NewCluster;
    PtrajCluster* ClusterA;
    PtrajCluster* ClusterB;
    ClusterNode* NewClusterNode;
    int CycleIndex;
    int FrameIndex, Frames, PointIndex;
    float DBI, pSF;
    int finish = 0;
    
    /**/
    NewCluster = (PtrajCluster *) ClusteringGetNewCluster( (Clustering *) This);
    ClusterAddEverything( (Cluster *) NewCluster);
    ClusteringAddCluster( (Clustering *)This, (Cluster *) NewCluster);
    CycleIndex = 0;

    FILE* SplitFile;
    SplitFile = fopen(ClusterSplittingFile,"w");
    int SplitIntName = -2;  /* -1 is uninitialized. */
    int OldIntName;
    NewCluster->IntName = SplitIntName;
    while (!finish)
    {
        CycleIndex++;
        /*PAPERGenerateClusterSnapshot(This, "Hier", CycleIndex);*/
		/*ClusteringOutputStatus(This);*/

        ClusteringFindEccentricCluster( (Clustering *) This, &EccentricEpsilon, 
					(Cluster **) &EccentricCluster,
					&EccentricPointA,  &EccentricPointB);
        
        if (!DesiredClusterCount && EccentricEpsilon<Epsilon)
        {
            finish = 1;
            continue;
        }
        
        OldIntName = EccentricCluster->IntName;
        NewClusterNode = PtrajClusteringSplitCluster( (PtrajClustering *) This,
						      (Cluster *) EccentricCluster,
						      EccentricPointA,EccentricPointB);
        ClusterA = (PtrajCluster *) NewClusterNode->Cluster;
        Frames = ClusterCountPointsInCluster( (Cluster *) ClusterA);
        if (Frames > 1)
        {
        	ClusterA->IntName = --SplitIntName;
        } else if (Frames == 1) {
        	ClusterA->IntName = FrameIndex;
        } else {
        	fprintf(stderr, "No frames in the cluster.\n");
        }
        ClusterB = (PtrajCluster *) NewClusterNode->pNext->Cluster;
        Frames = ClusterCountPointsInCluster( (Cluster *) ClusterB);
        if (Frames > 1)
        {
        	ClusterB->IntName = --SplitIntName;
        } else if (Frames == 1) {
        	ClusterB->IntName = FrameIndex;
        } else {
        	fprintf(stderr, "No frames in the cluster.\n");
        }
        fprintf(SplitFile, "%d:\t%d\t%d\t%8.3f", OldIntName, ((PtrajCluster *)ClusterA)->IntName, ((PtrajCluster *)ClusterB)->IntName, EccentricEpsilon);

        /*printf("%d:\t%d\t%d\t%8.3f", OldIntName, 
	  ((PtrajCluster *)ClusterA)->IntName, 
	  ((PtrajCluster *)ClusterB)->IntName, EccentricEpsilon);
	*/
        if (This->ClusterCount <= 50) {
        	DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
    	    fprintf(SplitFile, "\t%f\t%f\n", DBI, pSF);
	    /*printf("\t%f\t%f\n", DBI, pSF);*/
  		} else {
        	fprintf(SplitFile, "\n");
    	    printf("\n");
        }
        if (DesiredClusterCount>0 && This->ClusterCount>=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        if (DesiredClusterCount<0 && This->ClusterCount >= -DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
    }
    fclose(SplitFile);
}

#define SPLITTING_CYCLES 3
ClusterNode* PtrajClusteringSplitCluster(PtrajClustering* This, Cluster* ParentCluster, int PointA, int PointB)
{
   PtrajCluster* NewClusterA;
   PtrajCluster* NewClusterB;
   int PointCount;
   int PointIndex;
   float DistA;
   float DistB;
   int SplittingCycle;
   int C2Ccentroid = This->action->darg4;
   int OldBestRepA, OldBestRepB;
   int NewBestRepA, NewBestRepB;
   int verbose = This->action->darg2;
   ClusterNode * ReturnNode;
   /**/

   if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
   {
       fprintf(stdout, "Splitting Cluster %d with eccentric points %d, %d\n",
	       PtrajOutputIntermediateClusterStatus( (PtrajCluster *) ParentCluster, stdout, 1), PointA, PointB);
       PtrajOutputIntermediateClusterStatus( (PtrajCluster *) ParentCluster, stdout, 0); 
   }
   
   
   PointCount = ClusteringGetPointCount( (Clustering *) This);
   NewClusterA = PtrajClusterNew(This);
   NewClusterB = PtrajClusterNew(This);
   PtrajClusterSetBestRep(NewClusterA, PointA);
   PtrajClusterSetBestRep(NewClusterB, PointB);
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
	 /*ParentCluster->RemoveMember(PointIndex);*/
           PtrajClusterAddMember(NewClusterB,PointIndex);
       }
       else
       {
           PtrajClusterAddMember(NewClusterA,PointIndex);
       }
   }
   ReturnNode = ClusteringAddCluster( (Clustering *) This, (Cluster *) NewClusterA);
   ClusteringAddCluster( (Clustering *) This, (Cluster *) NewClusterB);
   /* Fix up the centroids, and re-assign points: */
   OldBestRepA = PtrajClusterGetBestRep(NewClusterA);
   OldBestRepB = PtrajClusterGetBestRep(NewClusterB);
   NewBestRepA = ClusterFindBestP2C(NewClusterA);
   NewBestRepB = ClusterFindBestP2C(NewClusterB);
   for (SplittingCycle = 0; SplittingCycle < SPLITTING_CYCLES; SplittingCycle++)
   {
	   if (C2Ccentroid > 0) {
	       if (OldBestRepA != NewBestRepA)
           {
	           ClusterAlignToBestRep(NewClusterA);
           }
   		   ClusterFindCentroid( (Cluster *) NewClusterA);
	       if (OldBestRepB != NewBestRepB)
           {
	           ClusterAlignToBestRep(NewClusterB);
           }
   		   ClusterFindCentroid( (Cluster *) NewClusterB);
 	   }
       if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
	   {
       	   fprintf(stdout, "Round %d\nNewClusterA is: \n", SplittingCycle);
	   	   PtrajOutputIntermediateClusterStatus(NewClusterA, stdout, 0); 
	   	   fprintf(stdout, "NewClusterB is:  \n");
	   	   PtrajOutputIntermediateClusterStatus(NewClusterB, stdout, 0); 
       }
       for (PointIndex=0; PointIndex<PointCount; PointIndex++)
       {
           if (!ClusterIsMember(ParentCluster,PointIndex))
           {
               continue;
           }
	   DistA = ClusteringDistanceToCentroid( (Clustering *) This, PointIndex, (Cluster *) NewClusterA);
           DistB = ClusteringDistanceToCentroid( (Clustering *) This, PointIndex, (Cluster *) NewClusterB);
           /*DistA = ClusteringDistanceToCluster(This,PointIndex,NewClusterA);
           DistB = ClusteringDistanceToCluster(This,PointIndex,NewClusterB);*/
           if (DistB < DistA && ClusterIsMember( (Cluster *) NewClusterA, PointIndex))
           {
               PtrajClusterRemoveMember(NewClusterA,PointIndex);
               PtrajClusterAddMember(NewClusterB,PointIndex);
           }
           else if (DistB >= DistA && ClusterIsMember( (Cluster *) NewClusterB, PointIndex))
           {
               PtrajClusterRemoveMember(NewClusterB,PointIndex);
               PtrajClusterAddMember(NewClusterA,PointIndex);
           }
       }
	   OldBestRepA = PtrajClusterGetBestRep(NewClusterA);
	   OldBestRepB = PtrajClusterGetBestRep(NewClusterB);
	   NewBestRepA = ClusterFindBestP2C( NewClusterA);
	   NewBestRepB = ClusterFindBestP2C( NewClusterB);
   }
   
   if (C2Ccentroid > 0) {
       if (OldBestRepA != NewBestRepA)
       {
	       ClusterAlignToBestRep(NewClusterA);
       }
	   ClusterFindCentroid( (Cluster *) NewClusterA);
       if (OldBestRepB != NewBestRepB)
       {
           ClusterAlignToBestRep(NewClusterB);
       }
	   ClusterFindCentroid( (Cluster *) NewClusterB);
   }
   ClusteringRemoveCluster( (Clustering *) This, (Cluster *) ParentCluster);
   ClusterFree( (Cluster *) ParentCluster); /* Node has been destroyed by ClusteringRemoveCluster. */
   if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
   {
	   fprintf(stdout, "into\n");
	   PtrajOutputIntermediateClusterStatus(NewClusterA, stdout, 0); 
	   PtrajOutputIntermediateClusterStatus(NewClusterB, stdout, 0); 
   	   ClusteringOutputClusterListToFile( (Clustering *) This, stdout,0); 
   }
   return(ReturnNode);
}

typedef struct StringList
{
  char* Str;
  struct StringList* next;
  struct StringList* prev;
} StringList;

/*void OutputClusteringStats(PtrajClustering* This, FILE* File);*/
/* Linkage clustering, while keeping track of centroid-to-centroid distances */
/* Special case: If DesiredClusterCount is NEGATIVE (-n), then keep joining and joining until there is only one cluster left.  
But when you are down to n clusters, start spitting out cluster.txt output at each cycle, and start building up a Newick Tree
string describing the hierarchical structure of our clusterings at each stage of the game! */
void OutputIntermediateClusters(PtrajClustering* This)
{
    char FilePath[1024];
    FILE* ClusteringFile;
    /**/
    sprintf(FilePath, "%s.%d.txt", This->action->carg1, This->ClusterCount);	  
    ClusteringFile = fopen(FilePath,"w");
    fprintf(ClusteringFile, "# Intermediate clustering tree\n");
    OutputClusteringStats(This, ClusteringFile);
    fclose(ClusteringFile);
    ClusteringAppendClusterList( (Clustering *) This, FilePath);
}


/* Produce a list of strings.  Each string is a name for the corresponding cluster.
   For now, the name is simply the number of points in the cluster. */
/*StringList* xxxNameClusters(PtrajClustering* This)
{
  StringList* Head = NULL;
  StringList* Tail = NULL;
  StringList* Node = NULL;
  ClusterNode* CNode;

  for (CNode = This->Head; CNode; CNode = CNode->pNext)
  {
    Node = (StringList*)calloc(sizeof(ClusterNode),1);
    Node->String = (char*)calloc(1024);
    sprintf(Node->String, "%d", ClusterCountPointsInCluster(CNode->Cluster));
    if (Head)
    {
      Tail->next = Node;
      Node->prev = Tail;
      Tail = Node;
    }
    else
    {
      Head = Node;
      Tail = Node;
    }
  }
  return Head;
}*/

void NameClusters(PtrajClustering* This)
{
  ClusterNode* CNode;
  PtrajCluster* pCluster;
  /**/
  for (CNode = This->Head; CNode; CNode = CNode->pNext)
  { 
    pCluster = (PtrajCluster*)CNode->Cluster;
    if (!pCluster->Name)
    {
      pCluster->Name = (char*)calloc(sizeof(char), 1024);
      sprintf(pCluster->Name, "%d", ClusterCountPointsInCluster( (Cluster *) pCluster));
    }
  }
}

void ClusteringMergeNames(PtrajClustering* This, ClusterNode* MergeNodeA, ClusterNode* MergeNodeB)
{
  char Temp[1024];
  PtrajCluster* ClusterA;
  PtrajCluster* ClusterB;
  ClusterA = (PtrajCluster*)MergeNodeA->Cluster;
  ClusterB = (PtrajCluster*)MergeNodeB->Cluster;
  if (!ClusterA->Name)
  {
    return;
  }
  sprintf(Temp, ClusterB->Name);
  sprintf(ClusterB->Name, "(%s,%s)", ClusterA->Name, Temp);

}

void PtrajClusteringOutputStatus(PtrajClustering* This)
{
    int DesiredClusterCount;
    actionInformation* action;
    action = This->action;
    DesiredClusterCount = action->iarg2;
    if (DesiredClusterCount < 0 && This->ClusterCount <= -DesiredClusterCount)
    {	  
      NameClusters(This);
      OutputIntermediateClusters(This);
    }
} 

/* Centroid-linkage clustering */
void PtrajClusteringClusterLinkage(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex;
    int PointCount = This->PointCount;
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    float Distance;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    StringList* NewickHead = NULL;
    int verbose = This->action->darg2;
    float DBI, pSF;
    int finish = 0;

    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining until you're finished */
    while (!finish)
    {
        PtrajClusteringOutputStatus(This);
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
		/*fprintf(stdout,"Linking.  Closest length %.4f vs epsilon %.4f", ClosestLength, Epsilon);*/
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
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "merge cluster %d into cluster %d\n",
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	  fprintf(stdout, "Mergee is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  fprintf(stdout, "Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
        	DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
    	    fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
  		} else {
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
        	fprintf(MergeFile, "\t%8s\t%f\n", " ", pSF);
        }
        
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
       	    fprintf(stdout, "After Merger is: \n");
	   	    PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
	
        /*ClusterFindCentroid(MergeNodeB->Cluster);  Done in ClusteringMergeClusters(). */
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update cookies */
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
            Distance = PtrajClusteringClusterToClusterDistance(This,NodeA, MergeNodeB);
            SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],Distance);
        }
        
    }
    fclose(MergeFile);
    safe_free(ClusterDistances);
}


/* Centroid-linkage-finger-print clustering */
void PtrajClusteringClusterLinkageFingerPrint(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex;
    int PointCount = This->PointCount;
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    double Distance, DistA, DistB;
    double ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    double C2CDistance;
    StringList* NewickHead = NULL;
    int verbose = This->action->darg2;
    float DBI, pSF;
    int finish = 0;
    int ClosestNodeACookie;
    DoubleSymmetricMatrix* DistanceCovariances;
    double *TotalDistances;
    int i, j;
    double xy, x, y, ii, ij, ji, jj;

    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);
    DistanceCovariances = AllocateDoubleSymmetricMatrix(This->PairwiseDistances->Size);
    TotalDistances = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double)*PointCount);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;
        
        TotalDistances[PointIndex] = 0;
        for (i = 0; i < PointCount; i++) {
        	Distance = GetSymElement(ClusterDistances, PointIndex, i);
            TotalDistances[PointIndex] += Distance;
        }
    }

    for (i = 0; i < PointCount; i++) {
    	for (j = i; j < PointCount; j++) {
        	xy = 0;
            for (PointIndex = 0; PointIndex < PointCount; PointIndex++) {
                xy += GetSymElement(ClusterDistances, PointIndex, i) * GetSymElement(ClusterDistances, PointIndex, j);
            }
  	        x = TotalDistances[i];
            y = TotalDistances[j];
	        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
        }
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining until you're finished */
    while (!finish)
    {
        PtrajClusteringOutputStatus(This);
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
        ClosestLength = GetDoubleSymElement(DistanceCovariances,ClusterCookies[ClosestNodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
        {
            for (NodeA = NodeB->pNext,NodeIndexA=NodeIndexB+1; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
            {
                C2CDistance = GetDoubleSymElement(DistanceCovariances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance > ClosestLength)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
		/*fprintf(stdout,"Linking.  Closest length %.4f vs epsilon %.4f", ClosestLength, Epsilon);*/
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
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "merge cluster %d into cluster %d\n",
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	  fprintf(stdout, "Mergee is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  fprintf(stdout, "Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
        	DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
    	    fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
  		} else {
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
        	fprintf(MergeFile, "\t%8s\t%f\n", " ", pSF);
        }
        fflush(MergeFile);
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "After Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
	
        /*ClusterFindCentroid(MergeNodeB->Cluster);  Done in ClusteringMergeClusters(). */
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update DistanceCovariance */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA]; /* Remember the old cookie for Node A. */
        for (i = 0; i < This->ClusterCount+1; i++) {
            DistA = GetSymElement(ClusterDistances,ClusterCookies[i],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[i],ClusterCookies[ClosestNodeIndexB]);
            /* Now the DistanceCovariances matrix temporarily becomes sum(xy) matrix, without including values involving merged cluster */
            for (j = i; j < This->ClusterCount+1; j++) {	
                xy = GetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j]) * (This->ClusterCount+1);
                xy += TotalDistances[ClusterCookies[i]] * TotalDistances[ClusterCookies[j]] / (This->ClusterCount+1);
                xy -= GetSymElement(ClusterDistances,ClusterCookies[j],ClosestNodeACookie) * DistA;
                xy -= GetSymElement(ClusterDistances,ClusterCookies[j],ClusterCookies[ClosestNodeIndexB]) * DistB;
                SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], xy);
            }
            TotalDistances[ClusterCookies[i]] -= (DistB + DistA); 
        }
        /* Update cookies */
        for (NodeIndexA=ClosestNodeIndexA; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
        }
        /*for (NodeIndexA=0; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            
            if (NodeIndexA>=ClosestNodeIndexA)
            {
                ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
            }
        }
        */
        if (ClosestNodeIndexB>ClosestNodeIndexA)
        {
            ClosestNodeIndexB--;
        }
        /* Update cluster-to-cluster distances */
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            Distance = PtrajClusteringClusterToClusterDistance(This,NodeA, MergeNodeB);
            SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],Distance);
            TotalDistances[ClusterCookies[NodeIndexA]] +=  Distance;
        }
        /* Fix the TotalDistances for the merged cluster */
        TotalDistances[ClusterCookies[ClosestNodeIndexB]] = 0;
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	TotalDistances[ClusterCookies[ClosestNodeIndexB]] += GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        }
        /* Add the xy value involving merged cluster, then convert it back to DistanceCovariance matrix */
        for (i = 0; i < This->ClusterCount; i++) {
            for (j = i; j < This->ClusterCount; j++) {	
                xy = GetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j]);
		        xy += GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[i]) * GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[j]);
    	        x = TotalDistances[ClusterCookies[i]];
	            y = TotalDistances[ClusterCookies[j]];
		        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
	        }
	    }
        /* Fix the covariances with the newly merged cluster. */
        for (i = 0; i < This->ClusterCount; i++) {
            xy = 0;
            for (j = 0; j < This->ClusterCount; j++) {	
        		xy += GetSymElement(ClusterDistances, ClusterCookies[i], ClusterCookies[j]) * GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[j]);
	        }
    	    x = TotalDistances[ClusterCookies[i]];
	        y = TotalDistances[ClusterCookies[ClosestNodeIndexB]];
	        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[ClosestNodeIndexB], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
	    }
        
    }
    fclose(MergeFile);
    FreeSymmetricMatrix(ClusterDistances);
    safe_free(ClusterCookies);
    FreeDoubleSymmetricMatrix(DistanceCovariances);
}

void PrintCookieSymMatrix(SymmetricMatrix* Matrix, int Count, int *Cookies)
{
    FILE* TheFile;
    int Row;
    int Column;
    float Value;
    /**/
    for (Row=0; Row<Count; Row++)
    {
        for (Column=0; Column<Count; Column++)
        {
		  Value = GetSymElement(Matrix,Cookies[Row],Cookies[Column]);
	      fprintf(stdout,"%6.2f ",Value);
	    }
	    fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
}

void PrintDoubleCookieSymMatrix(DoubleSymmetricMatrix* Matrix, int Count, int *Cookies)
{
    FILE* TheFile;
    int Row;
    int Column;
    double Value;
    /**/
    for (Row=0; Row<Count; Row++)
    {
        for (Column=0; Column<Count; Column++)
        {
		  Value = GetDoubleSymElement(Matrix,Cookies[Row],Cookies[Column]);
	      fprintf(stdout,"%6.2f ",Value);
	    }
	    fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
}


/* Average-linkage clustering */
void PtrajClusteringAverageLinkage(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int* ClusterPointCounts; /* Save for the point count for each cluster as in the cluster cookies. */
    int PointIndex, PointIndexA, PointIndexB;
    int PointCount = This->PointCount;
    PtrajCluster* NewCluster;
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
    StringList* NewickHead = NULL;
    int verbose = This->action->darg2;
    float DBI, pSF;
    int finish = 0;

    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterPointCounts = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex] = PointIndex;
        ClusterPointCounts[PointIndex] = 1;
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining until you're finished */
    while (!finish)
    {
        PtrajClusteringOutputStatus(This);
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
		/*fprintf(stdout,"Linking.  Closest length %.4f vs epsilon %.4f", ClosestLength, Epsilon);*/
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
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "merge cluster %d into cluster %d\n",
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	  fprintf(stdout, "Mergee is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  fprintf(stdout, "Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
        	DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
    	    fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
  		} else {
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
        	fprintf(MergeFile, "\t%8s\t%f\n", " ", pSF);
        }
        
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
       	    fprintf(stdout, "After Merger is: \n");
	   	    PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
	
        /*ClusterFindCentroid(MergeNodeB->Cluster);  Done in ClusteringMergeClusters(). */
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
    fclose(MergeFile);
    safe_free(ClusterDistances);
    safe_free(ClusterCookies);
    safe_free(ClusterPointCounts);
}

/* Average-linkage-finger-print clustering */
void PtrajClusteringAverageLinkageFingerPrint(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int* ClusterPointCounts; /* Save for the point count for each cluster as in the cluster cookies. */
    int PointIndex, PointIndexA, PointIndexB;
    int PointCount = This->PointCount;
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    int ClosestNodeACookie;
    int CountA, CountB;
    double DistA, DistB;
    int PointCountA, PointCountB;
    double Distance, NewDist;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    StringList* NewickHead = NULL;
    int verbose = This->action->darg2;
    float DBI, pSF;
    int finish = 0;
    DoubleSymmetricMatrix* DistanceCovariances;
    double *TotalDistances;
    int i, j;
    double xy, x, y, ii, ij, ji, jj;

    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterPointCounts = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);
    DistanceCovariances = AllocateDoubleSymmetricMatrix(This->PairwiseDistances->Size);
    TotalDistances = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double)*PointCount);


    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusterCookies[PointIndex] = PointIndex;
        ClusterPointCounts[PointIndex] = 1;
        
        TotalDistances[PointIndex] = 0;
        for (i = 0; i < PointCount; i++) {
        	Distance = GetSymElement(ClusterDistances, PointIndex, i);
            TotalDistances[PointIndex] += Distance;
        }
    }
    
    for (i = 0; i < PointCount; i++) {
    	for (j = i; j < PointCount; j++) {
        	xy = 0;
            for (PointIndex = 0; PointIndex < PointCount; PointIndex++) {
                xy += GetSymElement(ClusterDistances, PointIndex, i) * GetSymElement(ClusterDistances, PointIndex, j);
            }
            /*
            ii = GetSymElement(ClusterDistances, i, i);
            ij = GetSymElement(ClusterDistances, i, j);
            ji = GetSymElement(ClusterDistances, j, i);
            jj = GetSymElement(ClusterDistances, j, j);
            if (i != j) {
	            xy -= (ii * ij + ji * jj);
	            x = TotalDistances[i] - ii - ij;
	            y = TotalDistances[j] - ji - jj;
		        SetSymElement(DistanceCovariances, i, j, (xy - x * y / (PointCount-2)) / (PointCount-2));
            } else {
	            xy -= (ii * ij);
	            x = TotalDistances[i] - ii;
	            y = TotalDistances[j] - jj;
		        SetSymElement(DistanceCovariances, i, j, (xy - x * y / (PointCount-1)) / (PointCount-1));
            }*/
  	        x = TotalDistances[i];
            y = TotalDistances[j];
	        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
        }
    }
    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining until you're finished */
    while (!finish)
    {
        PtrajClusteringOutputStatus(This);
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
        ClosestLength = GetDoubleSymElement(DistanceCovariances,ClusterCookies[ClosestNodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
        {
            for (NodeA = NodeB->pNext,NodeIndexA=NodeIndexB+1; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
            {
                C2CDistance = GetDoubleSymElement(DistanceCovariances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance > ClosestLength)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
		/*fprintf(stdout,"Linking.  Closest length %.4f vs epsilon %.4f", ClosestLength, Epsilon);*/
        if (!DesiredClusterCount && ClosestLength < Epsilon)
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
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "merge cluster %d into cluster %d\n",
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	  fprintf(stdout, "Mergee is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  fprintf(stdout, "Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
        	DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
    	    fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
  		} else {
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
        	fprintf(MergeFile, "\t%8s\t%f\n", " ", pSF);
        }
        fflush(MergeFile);
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "After Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
	
        /*ClusterFindCentroid(MergeNodeB->Cluster);  Done in ClusteringMergeClusters(). */
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update DistanceCovariance */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA]; /* Remember the old cookie for Node A. */
        for (i = 0; i < This->ClusterCount+1; i++) {
            /*if (i == ClosestNodeIndexA) continue;*/
            DistA = GetSymElement(ClusterDistances,ClusterCookies[i],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[i],ClusterCookies[ClosestNodeIndexB]);
            /* Now the DistanceCovariances matrix temporarily becomes sum(xy) matrix, without including values involving merged cluster */
            for (j = i; j < This->ClusterCount+1; j++) {	
                xy = GetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j]) * (This->ClusterCount+1);
                xy += TotalDistances[ClusterCookies[i]] * TotalDistances[ClusterCookies[j]] / (This->ClusterCount+1);
                xy -= GetSymElement(ClusterDistances,ClusterCookies[j],ClosestNodeACookie) * DistA;
                xy -= GetSymElement(ClusterDistances,ClusterCookies[j],ClusterCookies[ClosestNodeIndexB]) * DistB;
                SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], xy);
            }
            TotalDistances[ClusterCookies[i]] -= (DistB + DistA); 
        }
        /* Update Cookies*/
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
            if (NodeIndexA == ClosestNodeIndexB) NewDist = 0;
	        SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],(float)NewDist);
            TotalDistances[ClusterCookies[NodeIndexA]] +=  NewDist;
        }
        ClusterPointCounts[ClusterCookies[ClosestNodeIndexB]] = CountA + CountB;
        /* Fix the TotalDistances for the merged cluster */
        TotalDistances[ClusterCookies[ClosestNodeIndexB]] = 0;
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	TotalDistances[ClusterCookies[ClosestNodeIndexB]] += GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        }
        /* Add the xy value involving merged cluster, then convert it back to DistanceCovariance matrix */
        for (i = 0; i < This->ClusterCount; i++) {
            for (j = i; j < This->ClusterCount; j++) {	
                xy = GetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j]);
		        xy += GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[i]) * GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[j]);
                /*
	            for (PointIndex = 0; PointIndex < This->ClusterCount; PointIndex++) {
		        	xy += GetSymElement(ClusterDistances, ClusterCookies[PointIndex], ClusterCookies[i]) * GetSymElement(ClusterDistances, ClusterCookies[PointIndex], ClusterCookies[j]);
                }
                */
	            /*ii = GetSymElement(ClusterDistances, ClusterCookies[i], ClusterCookies[i]);
    	        ij = GetSymElement(ClusterDistances, ClusterCookies[i], ClusterCookies[j]);
        	    ji = GetSymElement(ClusterDistances, ClusterCookies[j], ClusterCookies[i]);
            	jj = GetSymElement(ClusterDistances, ClusterCookies[j], ClusterCookies[j]);
	            if (i != j) {
            	    xy -= (ii * ij + ji * jj);
	    	        x = TotalDistances[ClusterCookies[i]] - ii - ij;
		            y = TotalDistances[ClusterCookies[j]] - ji - jj;
			        SetSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount-2)) / (This->ClusterCount-2));
                } else {
            	    xy -= ii * ij;
	    	        x = TotalDistances[ClusterCookies[i]] - ii;
		            y = TotalDistances[ClusterCookies[j]] - jj;
			        SetSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount-1)) / (This->ClusterCount-1));
                }*/
    	        x = TotalDistances[ClusterCookies[i]];
	            y = TotalDistances[ClusterCookies[j]];
		        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
	        }
	    }
        /* Fix the covariances with the newly merged cluster. */
        for (i = 0; i < This->ClusterCount; i++) {
            xy = 0;
            for (j = 0; j < This->ClusterCount; j++) {	
        		xy += GetSymElement(ClusterDistances, ClusterCookies[i], ClusterCookies[j]) * GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[j]);
	        }
    	    x = TotalDistances[ClusterCookies[i]];
	        y = TotalDistances[ClusterCookies[ClosestNodeIndexB]];
	        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[ClosestNodeIndexB], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
	    }
    }
    fclose(MergeFile);
    FreeSymmetricMatrix(ClusterDistances);
    safe_free(ClusterCookies);
    safe_free(ClusterPointCounts);
    FreeDoubleSymmetricMatrix(DistanceCovariances);
}

/* Edge-linkage clustering */
void PtrajClusteringEdgeLinkage(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex;
    int PointCount = This->PointCount;
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    float Distance;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    float DistA, DistB;
    int ClosestNodeACookie;
    StringList* NewickHead = NULL;
    int verbose = This->action->darg2;
    float DBI, pSF;
	int finish = 0;    

    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining until you're finished */
    while (!finish)
    {
        PtrajClusteringOutputStatus(This);
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
		/*fprintf(stdout,"Linking.  Closest length %.4f vs epsilon %.4f", ClosestLength, Epsilon);*/
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
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "merge cluster %d into cluster %d\n",
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	  fprintf(stdout, "Mergee is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  fprintf(stdout, "Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
        	DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
    	    fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
  		} else {
        	fprintf(MergeFile, "\n");
        }
        
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "After Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
	
        /*ClusterFindCentroid(MergeNodeB->Cluster);  Done in ClusteringMergeClusters(). */
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update cookies */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA];
        for (NodeIndexA=ClosestNodeIndexA; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
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
    fclose(MergeFile);
    safe_free(ClusterDistances);
}

/* Edge-linkage-finger-print clustering */
void PtrajClusteringEdgeLinkageFingerPrint(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex;
    int PointCount = This->PointCount;
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    double Distance, NewDist;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    double C2CDistance;
    double DistA, DistB;
    int ClosestNodeACookie;
    StringList* NewickHead = NULL;
    int verbose = This->action->darg2;
    float DBI, pSF;
	int finish = 0;    
    DoubleSymmetricMatrix* DistanceCovariances;
    double *TotalDistances;
    int i, j;
    double xy, x, y, ii, ij, ji, jj;
	
    
    
    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);
    DistanceCovariances = AllocateDoubleSymmetricMatrix(This->PairwiseDistances->Size);
    TotalDistances = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double)*PointCount);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;
        
        TotalDistances[PointIndex] = 0;
        for (i = 0; i < PointCount; i++) {
        	Distance = GetSymElement(ClusterDistances, PointIndex, i);
            TotalDistances[PointIndex] += Distance;
        }
    }

    for (i = 0; i < PointCount; i++) {
    	for (j = i; j < PointCount; j++) {
        	xy = 0;
            for (PointIndex = 0; PointIndex < PointCount; PointIndex++) {
                xy += GetSymElement(ClusterDistances, PointIndex, i) * GetSymElement(ClusterDistances, PointIndex, j);
            }
  	        x = TotalDistances[i];
            y = TotalDistances[j];
	        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
        }
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining until you're finished */
    while (!finish)
    {
        PtrajClusteringOutputStatus(This);
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
        ClosestLength = GetDoubleSymElement(DistanceCovariances,ClusterCookies[ClosestNodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
        {
            for (NodeA = NodeB->pNext,NodeIndexA=NodeIndexB+1; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
            {
                C2CDistance = GetDoubleSymElement(DistanceCovariances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance > ClosestLength)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
		/*fprintf(stdout,"Linking.  Closest length %.4f vs epsilon %.4f", ClosestLength, Epsilon);*/
        if (!DesiredClusterCount && ClosestLength < Epsilon)
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
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "merge cluster %d into cluster %d\n",
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	  fprintf(stdout, "Mergee is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  fprintf(stdout, "Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
	  DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
	  pSF = ClusteringComputePseudoF( (Clustering *) This);
	  fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
	} else {
	  fprintf(MergeFile, "\n");
        }
        fflush(MergeFile);
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "After Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
	
        /*ClusterFindCentroid(MergeNodeB->Cluster);  Done in ClusteringMergeClusters(). */
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update DistanceCovariance */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA]; /* Remember the old cookie for Node A. */
        for (i = 0; i < This->ClusterCount+1; i++) {
            DistA = GetSymElement(ClusterDistances,ClusterCookies[i],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[i],ClusterCookies[ClosestNodeIndexB]);
            /* Now the DistanceCovariances matrix temporarily becomes sum(xy) matrix, without including values involving merged cluster */
            for (j = i; j < This->ClusterCount+1; j++) {	
                xy = GetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j]) * (This->ClusterCount+1);
                xy += TotalDistances[ClusterCookies[i]] * TotalDistances[ClusterCookies[j]] / (This->ClusterCount+1);
                xy -= GetSymElement(ClusterDistances,ClusterCookies[j],ClosestNodeACookie) * DistA;
                xy -= GetSymElement(ClusterDistances,ClusterCookies[j],ClusterCookies[ClosestNodeIndexB]) * DistB;
                SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], xy);
            }
            TotalDistances[ClusterCookies[i]] -= (DistB + DistA); 
        }
        /* Update cookies */
        for (NodeIndexA=ClosestNodeIndexA; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
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
            NewDist = (DistA < DistB)? DistA : DistB;
            if (NodeIndexA == ClosestNodeIndexB) NewDist = 0;
	        SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],(float)NewDist);
            TotalDistances[ClusterCookies[NodeIndexA]] +=  NewDist;
        }
        /* Fix the TotalDistances for the merged cluster */
        TotalDistances[ClusterCookies[ClosestNodeIndexB]] = 0;
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	TotalDistances[ClusterCookies[ClosestNodeIndexB]] += GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        }
        /* Add the xy value involving merged cluster, then convert it back to DistanceCovariance matrix */
        for (i = 0; i < This->ClusterCount; i++) {
            for (j = i; j < This->ClusterCount; j++) {	
                xy = GetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j]);
		        xy += GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[i]) * GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[j]);
    	        x = TotalDistances[ClusterCookies[i]];
	            y = TotalDistances[ClusterCookies[j]];
		        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
	        }
	    }
        /* Fix the covariances with the newly merged cluster. */
        for (i = 0; i < This->ClusterCount; i++) {
            xy = 0;
            for (j = 0; j < This->ClusterCount; j++) {	
        		xy += GetSymElement(ClusterDistances, ClusterCookies[i], ClusterCookies[j]) * GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[j]);
	        }
    	    x = TotalDistances[ClusterCookies[i]];
	        y = TotalDistances[ClusterCookies[ClosestNodeIndexB]];
	        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[ClosestNodeIndexB], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
	    }
        
    }
    fclose(MergeFile);
    FreeSymmetricMatrix(ClusterDistances);
    safe_free(ClusterCookies);
    FreeDoubleSymmetricMatrix(DistanceCovariances);
}

/* Complete-linkage clustering */
void PtrajClusteringCompleteLinkage(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex;
    int PointCount = This->PointCount;
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    float Distance;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    float DistA, DistB;
    int ClosestNodeACookie;
    StringList* NewickHead = NULL;
    int verbose = This->action->darg2;
    float DBI, pSF;
    int finish = 0;

    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining until you're finished */
    while (!finish)
    {
        PtrajClusteringOutputStatus(This);
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
		/*fprintf(stdout,"Linking.  Closest length %.4f vs epsilon %.4f", ClosestLength, Epsilon);*/
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
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "merge cluster %d into cluster %d\n",
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	  fprintf(stdout, "Mergee is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  fprintf(stdout, "Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
        	DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
    	    fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
  		} else {
        	fprintf(MergeFile, "\n");
        }
        
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "After Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
	
        /*ClusterFindCentroid(MergeNodeB->Cluster);  Done in ClusteringMergeClusters(). */
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update cookies */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA];
        for (NodeIndexA=ClosestNodeIndexA; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
        }
        if (ClosestNodeIndexB>ClosestNodeIndexA)
        {
            ClosestNodeIndexB--;
        }
        /* Update cluster-to-cluster distances */
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            if (NodeIndexA == ClosestNodeIndexB) /*   !!!!!        Correct.
            if (NodeIndexA == ClusterCookies[ClosestNodeIndexB])   Incorrect. */
            	continue;
            DistA = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
            if (DistA > DistB)
	            SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],DistA);/*SetSymElement(ClusterDistances,NodeIndexA,ClosestNodeIndexB,DistA);*/
        }
        
    }
    fclose(MergeFile);
    safe_free(ClusterDistances);
}

/* Complete-linkage-finger-print clustering */
void PtrajClusteringCompleteLinkageFingerPrint(PtrajClustering* This,int DesiredClusterCount,float Epsilon)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex;
    int PointCount = This->PointCount;
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    double Distance, NewDist;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    double C2CDistance;
    double DistA, DistB;
    int ClosestNodeACookie;
    StringList* NewickHead = NULL;
    int verbose = This->action->darg2;
    float DBI, pSF;
    int finish = 0;
    DoubleSymmetricMatrix* DistanceCovariances;
    double *TotalDistances;
    int i, j;
    double xy, x, y, ii, ij, ji, jj;

    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);
    DistanceCovariances = AllocateDoubleSymmetricMatrix(This->PairwiseDistances->Size);
    TotalDistances = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double)*PointCount);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;
        
        TotalDistances[PointIndex] = 0;
        for (i = 0; i < PointCount; i++) {
        	Distance = GetSymElement(ClusterDistances, PointIndex, i);
            TotalDistances[PointIndex] += Distance;
        }
    }

    for (i = 0; i < PointCount; i++) {
    	for (j = i; j < PointCount; j++) {
        	xy = 0;
            for (PointIndex = 0; PointIndex < PointCount; PointIndex++) {
                xy += GetSymElement(ClusterDistances, PointIndex, i) * GetSymElement(ClusterDistances, PointIndex, j);
            }
  	        x = TotalDistances[i];
            y = TotalDistances[j];
	        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
        }
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining until you're finished */
    while (!finish)
    {
        PtrajClusteringOutputStatus(This);
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
        ClosestLength = GetDoubleSymElement(DistanceCovariances,ClusterCookies[ClosestNodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
        {
            for (NodeA = NodeB->pNext,NodeIndexA=NodeIndexB+1; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
            {
                C2CDistance = GetDoubleSymElement(DistanceCovariances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance > ClosestLength)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
		/*fprintf(stdout,"Linking.  Closest length %.4f vs epsilon %.4f", ClosestLength, Epsilon);*/
        if (!DesiredClusterCount && ClosestLength < Epsilon)
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
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "merge cluster %d into cluster %d\n",
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	  fprintf(stdout, "Mergee is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  fprintf(stdout, "Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
	  DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
	  pSF = ClusteringComputePseudoF( (Clustering *) This);
	  fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
	} else {
	  fprintf(MergeFile, "\n");
        }
        fflush(MergeFile);
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	  fprintf(stdout, "After Merger is: \n");
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
	
        /*ClusterFindCentroid(MergeNodeB->Cluster);  Done in ClusteringMergeClusters(). */
        if (This->ClusterCount<=DesiredClusterCount)
        {
            finish = 1;
            continue;
        }
        /* Update DistanceCovariance */
        ClosestNodeACookie = ClusterCookies[ClosestNodeIndexA]; /* Remember the old cookie for Node A. */
        for (i = 0; i < This->ClusterCount+1; i++) {
            DistA = GetSymElement(ClusterDistances,ClusterCookies[i],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[i],ClusterCookies[ClosestNodeIndexB]);
            /* Now the DistanceCovariances matrix temporarily becomes sum(xy) matrix, without including values involving merged cluster */
            for (j = i; j < This->ClusterCount+1; j++) {	
                xy = GetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j]) * (This->ClusterCount+1);
                xy += TotalDistances[ClusterCookies[i]] * TotalDistances[ClusterCookies[j]] / (This->ClusterCount+1);
                xy -= GetSymElement(ClusterDistances,ClusterCookies[j],ClosestNodeACookie) * DistA;
                xy -= GetSymElement(ClusterDistances,ClusterCookies[j],ClusterCookies[ClosestNodeIndexB]) * DistB;
                SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], xy);
            }
            TotalDistances[ClusterCookies[i]] -= (DistB + DistA); 
        }
        /* Update cookies */
        for (NodeIndexA=ClosestNodeIndexA; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
        }
        if (ClosestNodeIndexB>ClosestNodeIndexA)
        {
            ClosestNodeIndexB--;
        }
        /* Update cluster-to-cluster distances */
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            if (NodeIndexA == ClosestNodeIndexB) /*   !!!!!        Correct.
            if (NodeIndexA == ClusterCookies[ClosestNodeIndexB])   Incorrect. */
            	continue;
            DistA = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClosestNodeACookie);
            DistB = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
            NewDist = (DistA > DistB)? DistA : DistB;
            if (NodeIndexA == ClosestNodeIndexB) NewDist = 0;
	        SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],(float)NewDist);
            TotalDistances[ClusterCookies[NodeIndexA]] +=  NewDist;
        }
        /* Fix the TotalDistances for the merged cluster */
        TotalDistances[ClusterCookies[ClosestNodeIndexB]] = 0;
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	TotalDistances[ClusterCookies[ClosestNodeIndexB]] += GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB]);
        }
        /* Add the xy value involving merged cluster, then convert it back to DistanceCovariance matrix */
        for (i = 0; i < This->ClusterCount; i++) {
            for (j = i; j < This->ClusterCount; j++) {	
                xy = GetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j]);
		        xy += GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[i]) * GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[j]);
    	        x = TotalDistances[ClusterCookies[i]];
	            y = TotalDistances[ClusterCookies[j]];
		        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[j], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
	        }
	    }
        /* Fix the covariances with the newly merged cluster. */
        for (i = 0; i < This->ClusterCount; i++) {
            xy = 0;
            for (j = 0; j < This->ClusterCount; j++) {	
        		xy += GetSymElement(ClusterDistances, ClusterCookies[i], ClusterCookies[j]) * GetSymElement(ClusterDistances, ClusterCookies[ClosestNodeIndexB], ClusterCookies[j]);
	        }
    	    x = TotalDistances[ClusterCookies[i]];
	        y = TotalDistances[ClusterCookies[ClosestNodeIndexB]];
	        SetDoubleSymElement(DistanceCovariances, ClusterCookies[i], ClusterCookies[ClosestNodeIndexB], (xy - x*y / (This->ClusterCount)) / (This->ClusterCount));
	    }
        
    }
    fclose(MergeFile);
    FreeSymmetricMatrix(ClusterDistances);
    safe_free(ClusterCookies);
    FreeDoubleSymmetricMatrix(DistanceCovariances);
}

void PtrajClusteringBuildTotalCluster(Clustering* This) 
{
  PtrajCluster* TempCluster;
  float TotalDistanceToCentroid;
  float Distance;
  int PointIndex;
  
  TempCluster = (PtrajCluster *) ClusteringGetNewCluster( (Clustering *) This);
  ClusterAddEverything( (Cluster *) TempCluster);
  ClusterFindCentroid( (Cluster *)  TempCluster);
  This->TotalCluster = (Cluster *) TempCluster;

  TotalDistanceToCentroid = 0;
  for (PointIndex=0; PointIndex < This->PointCount; PointIndex++)
  {
    Distance = ClusteringDistanceToCentroid( (Clustering *) This, PointIndex, (Cluster *) TempCluster);
    Distance = Distance*Distance;
    /* fprintf(stdout, "Point %d distance to center %.4f\n", PointIndex, Distance); */
    TotalDistanceToCentroid += Distance;
  }
  This->SST = TotalDistanceToCentroid;
}

float PtrajClusteringPointToPointDistance(PtrajClustering* This, int A, int B)
{	
	if (This->PairwiseDistances)
	    return GetSymElement(This->PairwiseDistances,A,B);
    else {
    	float *XA, *YA, *ZA, *XB, *YB, *ZB;
        trajectoryInfo* trajInfo;
        trajInfo = This->action->carg2;
        XA = trajInfo->x + trajInfo->atoms*A;
        YA = trajInfo->y + trajInfo->atoms*A;
        ZA = trajInfo->z + trajInfo->atoms*A;
        XB = trajInfo->x + trajInfo->atoms*B;
        YB = trajInfo->y + trajInfo->atoms*B;
        ZB = trajInfo->z + trajInfo->atoms*B;
        return  PtrajGetDistance(This, A, B, XA, YA, ZA, XB, YB, ZB);
     }
}

float PtrajClusteringCentroidToCentroidDistance(PtrajClustering* This,ClusterNode* NodeA,ClusterNode* NodeB)
{
    PtrajCluster* ClusterA = (PtrajCluster*)NodeA->Cluster;
    PtrajCluster* ClusterB = (PtrajCluster*)NodeB->Cluster;

    return PtrajGetDistance(This, -1, -1, ClusterA->CentroidX, ClusterA->CentroidY, ClusterA->CentroidZ,
			    ClusterB->CentroidX, ClusterB->CentroidY, ClusterB->CentroidZ);
      /*return (float)rmsf(This->trajInfo->atoms,1,This->trajInfo->state->masses,NULL,
        ClusterA->CentroidX,ClusterA->CentroidY,ClusterA->CentroidZ,
        ClusterB->CentroidX,ClusterB->CentroidY,ClusterB->CentroidZ,
        NULL,NULL,0);*/

}

float PtrajClusteringClusterToClusterDistance(PtrajClustering* This,ClusterNode* NodeA,ClusterNode* NodeB)
{
    PtrajCluster* ClusterA = (PtrajCluster*)NodeA->Cluster;
    PtrajCluster* ClusterB = (PtrajCluster*)NodeB->Cluster;
	int C2Ccentroid = This->action->darg4;
    float* X;
    float* Y;
    float* Z;
    float* XB;
    float* YB;
    float* ZB;
    
    if (C2Ccentroid == 0) 
    {
    	if (This->action->performSecondPass) 
        {
		    /*
            X = This->trajInfo->x + ClusterA->BestRep*This->trajInfo->atoms;
		    Y = This->trajInfo->y + ClusterA->BestRep*This->trajInfo->atoms;
		    Z = This->trajInfo->z + ClusterA->BestRep*This->trajInfo->atoms;
            XB = This->trajInfo->x + ClusterB->BestRep*This->trajInfo->atoms;
		    YB = This->trajInfo->y + ClusterB->BestRep*This->trajInfo->atoms;
		    ZB = This->trajInfo->z + ClusterB->BestRep*This->trajInfo->atoms;
		    return PtrajGetDistance(This, X, Y, Z, XB, YB,ZB);
            */
        	return 	GetSymElement(This->PairwiseDistances,ClusterA->OldBestRep,ClusterB->OldBestRep);
        }
        else
        {
        	return 	GetSymElement(This->PairwiseDistances,ClusterA->BestRep,ClusterB->BestRep);
        }
    }
    return PtrajGetDistance(This, -1, -1, ClusterA->CentroidX, ClusterA->CentroidY, ClusterA->CentroidZ,
			    ClusterB->CentroidX, ClusterB->CentroidY, ClusterB->CentroidZ);
      /*return (float)rmsf(This->trajInfo->atoms,1,This->trajInfo->state->masses,NULL,
        ClusterA->CentroidX,ClusterA->CentroidY,ClusterA->CentroidZ,
        ClusterB->CentroidX,ClusterB->CentroidY,ClusterB->CentroidZ,
        NULL,NULL,0);*/

}

void PtrajClusterExpandClusterSize(PtrajCluster* This,int Multiplier, int NewFrameCount)
{
    int* NewMask;
    int FrameIndex;
    int OldMaskIndex;
    /**/
    NewMask = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*NewFrameCount);
    memset(NewMask,0,sizeof(int)*NewFrameCount);
    OldMaskIndex = 0;
    for (FrameIndex=0; FrameIndex<NewFrameCount; FrameIndex++)
    {
      if (This->Owner->SieveFirstPass[FrameIndex])
      {
        NewMask[FrameIndex] = This->Mask[OldMaskIndex];
		OldMaskIndex++;
      }
    }
    safe_free(This->Mask);
    This->Mask = NewMask;
    This->PointCount = NewFrameCount;
}

int PtrajClusterGetNewIndex(PtrajCluster* This, int OldIndex, int NewFrameCount) 
{
	int FrameIndex;
    int OldMaskIndex;
    
    OldMaskIndex = 0;
    for (FrameIndex=0; FrameIndex<NewFrameCount; FrameIndex++)
    {
      if (This->Owner->SieveFirstPass[FrameIndex]) {
		OldMaskIndex++;
  	    if (OldMaskIndex == OldIndex + 1)
          return(FrameIndex);
      }
    }
    return (OldIndex); /* Did not find new index, probably -1. */
    
}
void PtrajClusteringExpandClusterListSizes(PtrajClustering* This,int Multiplier, int NewFrameCount)
{
    ClusterNode* OldNode;
    int AtomCount = PtrajClusteringGetAtomCount(This);
    PtrajCluster* TheCluster;
    int OldBestRep, NewBestRep;
    /**/
    This->PointCount = NewFrameCount;
    for (OldNode=This->Head; OldNode; OldNode = OldNode->pNext)
    {
      TheCluster = ((PtrajCluster*)(OldNode->Cluster));
      OldBestRep = PtrajClusterGetBestRep(TheCluster);
      /* When OldBestRep is less than 0 (-1), the Cluster contains only one item and is initalized by ClusterAddMember,
         so no BestRep has been updated. The following code try to fix this problem. */
      if (OldBestRep < 0) 
      {
      	int i;
        for (i = 0; i < TheCluster->PointCount; i++)
        {
        	if (TheCluster->Mask[i])
            {
            	if (OldBestRep < 0) 
                	OldBestRep = i;
                else
                	warning("PtrajClusteringExpandClusterListSizes()", "Cluster has more than 1 memeber but without proper BestRep.");
            }
        }
      
      }
      PtrajClusterExpandClusterSize(TheCluster,Multiplier,NewFrameCount);
      PtrajClusterSetBestRep(TheCluster, PtrajClusterGetNewIndex(TheCluster, OldBestRep, NewFrameCount));
      TheCluster->OldBestRep = OldBestRep;
    }

}

void PtrajClusteringAddFrameToBestCluster(PtrajClustering* This,int FrameIndex, float* x, float* y, float* z)
{
    ClusterNode* Node;
    PtrajCluster* TheCluster;
    PtrajCluster* BestCluster;
    int ClusterIndex;    
    float BestDistance = -1;
    float Distance;
    int AtomIndex;
    char FileName[1024];
    FILE* DebugFile;
  	int C2Ccentroid =  This->action->darg4; 
    float* X;
    float* Y;
    float* Z;
    int BestRep;
    trajectoryInfo* trajInfo = This->trajInfo;
    int mass = This->action->iarg6;
	int verbose = This->action->darg2;
    float rmsT[3],rmsR[3][3];
    /**/
    /*fprintf(stdout, "Add %d to best cluster:\n", FrameIndex);*/
    for (Node = This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
    {
        TheCluster = (PtrajCluster*)Node->Cluster;
        if (TheCluster->BestRep < 0)
        {
        	error("PtrajClusteringAddFrameToBestCluster", "BestRep is less than 0, no member is included in cluster %d\n", ClusterIndex);
        }
        else 
        {
        	BestRep = TheCluster->OldBestRep;
        }
		if (C2Ccentroid == 0) 
        {
            X = This->trajInfo->x + BestRep*This->trajInfo->atoms;
		    Y = This->trajInfo->y + BestRep*This->trajInfo->atoms;
		    Z = This->trajInfo->z + BestRep*This->trajInfo->atoms;
        }
        else
        {
        	X = TheCluster->CentroidX;
        	Y = TheCluster->CentroidY;
        	Z = TheCluster->CentroidZ;
        }
        Distance = PtrajGetDistance(This, BestRep, -1, X,Y,Z, x,y,z);
	/*fprintf(stdout,"Assigning %d: Cluster %d has distance %.4f\n",FrameIndex,ClusterIndex,Distance);*/
	/*fprintf(stdout," %.4f ",Distance);*/
        if (BestDistance<0 || Distance<BestDistance)
        {
	  /*fprintf(stdout,"Cluster %d is the best so far\n",ClusterIndex);*/
            BestDistance = Distance;
            BestCluster = TheCluster;
        }
    }

    /*fprintf(stdout, "Add point %d to cluster %d\n", FrameIndex, PtrajOutputIntermediateClusterStatus(BestCluster, stdout, 1));*/
   	if (verbose >= 6)
    {
       PtrajOutputIntermediateClusterStatus(BestCluster, stdout, 0); 
    }
    ClusterAddMember( (Cluster *) BestCluster,FrameIndex);

    
    /* Caution: Also the BestRep in each cluster has been updated due to second pass expanding, 
       the trajectories in trajInfo has not been expanded, so the OldBestRep hold their positions in the trajInfo. 
       So we can compare new frames to the bestrep(old numbering) of each cluster.      --Jianyin */
    X = &trajInfo->x[BestCluster->OldBestRep * trajInfo->atoms]; 
    Y = &trajInfo->y[BestCluster->OldBestRep * trajInfo->atoms]; 
    Z = &trajInfo->z[BestCluster->OldBestRep * trajInfo->atoms]; 
    if (mass) {
    	rmsf(trajInfo->atoms,1,trajInfo->state->masses,NULL,x,y,x,X,Y,Z,rmsR,rmsT,2);
    } else {
    	rmsf(trajInfo->atoms,1,NULL,NULL,x,y,x,X,Y,Z,rmsR,rmsT,2);
    }
    if (verbose >= VERBOSE_ALIGN_FRAME) fprintf(stdout, "Align frame %i to frame %i,\n", FrameIndex, BestCluster->BestRep);
    
    /* Debugging: Dump the frame coordinates! */
    /*sprintf(FileName,"Frame%d",FrameIndex);
    fprintf(stdout,"File name is %s\n",FileName);
    DebugFile = fopen(FileName,"w");
    for (AtomIndex=0; AtomIndex<This->trajInfo->atoms; AtomIndex++)
    {
      fprintf(DebugFile,"%.4f\t%.4f\t%.4f\n",x[AtomIndex],y[AtomIndex],z[AtomIndex]);
    }
    fclose(DebugFile);*/
}

Representative* RepresentativeNew(int PointIndex,int AtomCount)
{
  Representative* This;
  This = (Representative*)SafeMalloc(__FILE__, __LINE__, sizeof(Representative));
  memset(This,0,sizeof(Representative));
    This->PointIndex = PointIndex;
    This->MovedFlag = 0;
    This->X = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount);
    memset(This->X,0,sizeof(float)*AtomCount);
    This->Y = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount);
    memset(This->Y,0,sizeof(float)*AtomCount);
    This->Z = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount);
    memset(This->Z,0,sizeof(float)*AtomCount);
    return This;
}

void RepresentativeFree(Representative* This)
{
    safe_free(This->X);
    safe_free(This->Y);
    safe_free(This->Z);
    safe_free(This);
}

#define MOVE_TO_CENTROID_FACTOR 0.25

/* Take the representative points and move them toward the centroid */
void PtrajClusterRepositionRepresentatives(PtrajCluster* This)
{
    RepresentativeNode* Node;
    Representative* Rep;
    int AtomIndex;
    int AtomCount = PtrajClusteringGetAtomCount(This->Owner);
    /*int PointCount = This->Owner->FrameCount;*/
    trajectoryInfo* trajInfo = This->Owner->trajInfo;
    float CentroidValue;
    float PointValue;
    /* */
    for (Node = This->Head; Node; Node = Node->pNext)
    {
        Rep = Node->Point;
        Rep->MovedFlag = 1;
        for (AtomIndex=0; AtomIndex<AtomCount; AtomIndex++)
        {
            CentroidValue = This->CentroidX[AtomIndex];
            PointValue = trajInfo->x[Rep->PointIndex*trajInfo->atoms + AtomIndex];
            Rep->X[AtomIndex] = (float)(PointValue + (CentroidValue-PointValue)*MOVE_TO_CENTROID_FACTOR);
            CentroidValue = This->CentroidY[AtomIndex];
            PointValue = trajInfo->y[Rep->PointIndex*trajInfo->atoms + AtomIndex];
            Rep->Y[AtomIndex] = (float)(PointValue + (CentroidValue-PointValue)*MOVE_TO_CENTROID_FACTOR);
            CentroidValue = This->CentroidZ[AtomIndex];
            PointValue = trajInfo->z[Rep->PointIndex*trajInfo->atoms + AtomIndex];
            Rep->Z[AtomIndex] = (float)(PointValue + (CentroidValue-PointValue)*MOVE_TO_CENTROID_FACTOR);
        }
    }
}

int PtrajClusterGetRepresentativeCount(PtrajCluster* This)
{
    int RepresentativeCount = 0;
    RepresentativeNode* Node;
    for (Node = This->Head; Node; Node = Node->pNext)
    {
        RepresentativeCount++;
    }
    return RepresentativeCount;
}

void PtrajClusterPruneRepresentatives(PtrajCluster* This,int MaxRepresentatives)
{
    int RepresentativeCount;
    int RepIndex;
    int OtherRepIndex;
    RepresentativeNode* Node;
    RepresentativeNode* OtherNode;
    int* RepCookies;
    float Distance;
    float* RepSpreading;
    SymmetricMatrix* RepDistance = NULL;
    float MinDistance;
    int MinDistanceRep;
    /* */
    /* Check to see how many we need to prune, and return if we're okay already */
    RepresentativeCount = PtrajClusterGetRepresentativeCount(This);
    if (RepresentativeCount<=MaxRepresentatives)
    {
        return;
    }

    RepDistance = AllocateSymmetricMatrix(RepresentativeCount);

    /* Find our representatives' distance from each other */
    RepCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*RepresentativeCount);
    RepSpreading = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*RepresentativeCount);
    memset(RepCookies,0,sizeof(int)*RepresentativeCount);
    /*fprintf(stdout,"RepCookies reset%d\n",This->Owner->FrameCount);*/
    for (RepIndex=0; RepIndex<RepresentativeCount; RepIndex++)
    {
        RepCookies[RepIndex]=RepIndex;
        RepSpreading[RepIndex]=0;
    }
    for (Node=This->Head,RepIndex=0; Node; Node=Node->pNext,RepIndex++)
    {
      /*fprintf(stdout,"Get spread for %d:%d",RepIndex,RepCookies[RepIndex]);*/
        OtherNode=Node->pNext;
        OtherRepIndex=RepIndex+1; 
		while (OtherNode) 
        {
            Distance = PtrajClusteringRepToRepDistance(This->Owner,Node->Point,OtherNode->Point);
            SymmetricMatrixElement(RepDistance,RepIndex,OtherRepIndex)  = Distance;
            RepSpreading[RepIndex] += Distance;
            RepSpreading[OtherRepIndex] += Distance;
	    	OtherNode=OtherNode->pNext;
		    OtherRepIndex++;
        }
    }

    /* Now trim off the first few reps */
    while (RepresentativeCount>MaxRepresentatives)
    {
        /* Find the rep to kill: */
        MinDistance = -1;
        for (RepIndex=0; RepIndex<RepresentativeCount; RepIndex++)
        {
            if (MinDistance<0 || RepSpreading[RepCookies[RepIndex]] < MinDistance)
            {
                MinDistance = RepSpreading[RepCookies[RepIndex]];
                MinDistanceRep = RepIndex;
            }
        }
        /* Update the total-distances of other reps: */
        for (OtherRepIndex=0; OtherRepIndex<RepresentativeCount; OtherRepIndex++)
        {
            Distance = GetSymElement(RepDistance,RepCookies[MinDistanceRep],RepCookies[OtherRepIndex]);
            RepSpreading[RepCookies[OtherRepIndex]] -= Distance;
        }
        /* Yank the rep: */
        for (Node=This->Head,RepIndex=0; Node; Node = Node->pNext,RepIndex++)
        {
            if (RepIndex==MinDistanceRep)
            {
                if (Node->pNext)
                {
                    Node->pNext->pPrev = Node->pPrev;
                }
                if (Node->pPrev)
                {
                    Node->pPrev->pNext = Node->pNext;
                }
				if (This->Head==Node)
                {
				    This->Head = Node->pNext;
				}
				if (This->Tail==Node)
                {
		    		This->Tail = Node->pPrev;
				}
                RepresentativeCount--;
                RepresentativeFree(Node->Point);
                safe_free(Node);
                break;
            }
        }
        /* Update the cookies: */
        for (RepIndex=0; RepIndex<RepresentativeCount; RepIndex++)
        {
	        /*fprintf(stdout,"Update rep %d from cookie old %d to cookie new %d",RepIndex,RepCookies[RepIndex],RepCookies[RepIndex+1]);*/
            if (RepIndex >= MinDistanceRep)
            {
                RepCookies[RepIndex] = RepCookies[RepIndex+1];
            }
        }
    }
cleanup:
    safe_free(RepCookies);
    safe_free(RepSpreading);
    FreeSymmetricMatrix(RepDistance);
}

void PointToCoordinates(Representative* Rep, trajectoryInfo* trajInfo, float** X,float** Y, float** Z)
{
    if (Rep->MovedFlag)
    {
        *X = Rep->X;
        *Y = Rep->Y;
        *Z = Rep->Z;
    }
    else
    {
        *X = &trajInfo->x[Rep->PointIndex*trajInfo->atoms];
        *Y = &trajInfo->y[Rep->PointIndex*trajInfo->atoms];
        *Z = &trajInfo->z[Rep->PointIndex*trajInfo->atoms];
    }

}

float PtrajClusteringRepToRepDistance(PtrajClustering* This,Representative* RepA, Representative* RepB)
{
    float* XA;
    float* YA;
    float* ZA;
    float* XB;
    float* YB;
    float* ZB;
    float Distance;
    trajectoryInfo* trajInfo = This->trajInfo;
    /**/
    PointToCoordinates(RepA,This->trajInfo,&XA,&YA,&ZA);
    PointToCoordinates(RepB,This->trajInfo,&XB,&YB,&ZB);

    Distance=PtrajGetDistance(This, -1, -1, XA,YA,ZA, XB,YB,ZB);
    /* Distance=(float)rmsf(trajInfo->atoms,1,trajInfo->state->masses,NULL,XA,YA,ZA,
       XB,YB,ZB,NULL,NULL,0);*/
    return Distance;
}

void CopyDoubleToFloat(float* Target, double* Source, int Count)
{
  int Index;
  for (Index = 0; Index < Count; Index++)
    {
      Target[Index] = (float)Source[Index];
    }
}

/* Select a new collection of representatives for a cluster of points.  Choose points that are maximally distant
from each other: The first point is taken from the old representative list.  For each other representative, choose
the cluster point whose total distance from the other representatives is maximal.  (The first point may not
be as "edgy" as we would like, but choosing n good points is much harder than choosing n-1 good points and 1 decent 
point, so if you want a better spread, just increase the number of representatives) */

void PtrajClusterSelectRepresentatives(PtrajClustering* This, PtrajCluster* MergedCluster, int Representatives)
{
  int RepresentativeIndex;
  float BestTotalDistance;
  int BestDistanceIndex;
  int* UsedAlready;
  float Distance = 0;
  RepresentativeNode* Node;
  RepresentativeNode* OtherNode;
  RepresentativeNode* BestNode = NULL;
  SymmetricMatrix* PointDistances = This->PairwiseDistances;
  Representative* NewRep;
  RepresentativeNode* NewRepNode;
  int PointIndex;
  float TempDistance;
  trajectoryInfo* trajInfo = This->trajInfo;
  /**/

  UsedAlready = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * This->PointCount);
  memset(UsedAlready, 0, sizeof(int) * This->PointCount);

  /*****************************************/
  /* Reset the representative list: Keep the maximally-distant one*/
  BestTotalDistance = 0;
  for (Node = MergedCluster->Head; Node; Node = Node->pNext)
  {
      Distance = 0;
      for (OtherNode = MergedCluster->Head; OtherNode; OtherNode = OtherNode->pNext)
	  {
	    if (Node!=OtherNode)
	    {
	      Distance += GetSymElement(PointDistances, Node->Point->PointIndex, OtherNode->Point->PointIndex);
	    }
	  }      
      if (Distance >= BestTotalDistance)
      {
	  	BestTotalDistance = Distance;
	  	BestNode = Node;
	  }
  }
  for (Node = MergedCluster->Head; Node;)
  {
      if (Node != BestNode)
      {
          if (Node->pNext) 
          {
              Node = Node->pNext;
              RepresentativeFree(Node->pPrev->Point);
              safe_free(Node->pPrev);
          }
          else 
          {
              RepresentativeFree(Node->Point);
              safe_free(Node);
              Node = NULL;
          }	 	  
      }
      else
      {
          Node = Node->pNext;
      }
  } 
  MergedCluster->Head = BestNode;
  MergedCluster->Tail = BestNode;
  BestNode->pNext = NULL;
  BestNode->pPrev = NULL;
  /*CopyDoubleToFloat(BestNode->Point->X, trajInfo->x + (BestNode->Point->PointIndex * trajInfo->atoms), trajInfo->atoms);*/
  memcpy(BestNode->Point->X, trajInfo->x + (BestNode->Point->PointIndex * trajInfo->atoms), sizeof(float) * trajInfo->atoms);
  memcpy(BestNode->Point->Y, trajInfo->y + (BestNode->Point->PointIndex * trajInfo->atoms), sizeof(float) * trajInfo->atoms);
  memcpy(BestNode->Point->Z, trajInfo->z + (BestNode->Point->PointIndex * trajInfo->atoms), sizeof(float) * trajInfo->atoms);

  UsedAlready[BestNode->Point->PointIndex] = 1;
  /*************************************************************/
  /* Choose new points! */
  for (RepresentativeIndex = 1; RepresentativeIndex < Representatives; RepresentativeIndex++)
  {
    BestTotalDistance = 0;
    for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
    {
		if ( ! ClusterIsMember( (Cluster *) MergedCluster, PointIndex ) )
		{
	    	continue;
		}
		if (UsedAlready[PointIndex])
		{
		    continue;
		}
		Distance = 0;
		for (Node = MergedCluster->Head; Node; Node = Node->pNext)
        {
		  /*Distance += GetSymElement(PointDistances, PointIndex, Node->Point->PointIndex);	    */
		  TempDistance = GetSymElement(PointDistances, PointIndex, Node->Point->PointIndex);
		  if (Distance == 0 || (TempDistance > 0.0 && TempDistance < Distance))
		  {
	       	Distance = TempDistance;
	      }
		}	
		if (Distance > BestTotalDistance)
	    {
	    	BestTotalDistance = Distance;
	    	BestDistanceIndex = PointIndex;
	    }
    }
    /* If we've used all the points, stop now: */
    if (BestTotalDistance == 0.0)
    {
		break;
    }
    /* Add representative point: */
    NewRep = RepresentativeNew(BestDistanceIndex,PtrajClusteringGetAtomCount(This));
    memcpy(NewRep->X, trajInfo->x + (BestDistanceIndex * trajInfo->atoms), sizeof(float) * trajInfo->atoms);
    memcpy(NewRep->Y, trajInfo->y + (BestDistanceIndex * trajInfo->atoms), sizeof(float) * trajInfo->atoms);
    memcpy(NewRep->Z, trajInfo->z + (BestDistanceIndex * trajInfo->atoms), sizeof(float) * trajInfo->atoms);
    
    NewRepNode = (RepresentativeNode*)SafeMalloc(__FILE__, __LINE__, sizeof(RepresentativeNode));
    NewRepNode->pNext = NULL;
    NewRepNode->pPrev = MergedCluster->Tail;
    MergedCluster->Tail->pNext = NewRepNode;
    MergedCluster->Tail = NewRepNode;
    NewRepNode->Point = NewRep;
    UsedAlready[BestDistanceIndex] = 1;
  }
  safe_free(UsedAlready);
}


/* A CURE-like clustering algorithm.  Start by putting each point in its own teeny cluster, and making it the
representative of the cluster.  Then cycle through, merging a pair of clusters at each step of the cycle,
until the number of clusters is as low as requested.

Deciding who to merge: Find the pair of representative points (X,Y) such that X and Y are as close as possible
where X and Y are in different clusters.

How to merge: All points from the two clusters are put into the new cluster.  Recompute the centroid.  Then choose
up to <Representatives> new representative points for the new cluster, and shrink them by <Alpha> toward the 
centroid.  
*/
void PtrajClusteringClusterCentripetal(PtrajClustering* This,int DesiredClusterCount, float Epsilon, int Representatives)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex;
    int PointCount = This->PointCount;
    int AtomCount = PtrajClusteringGetAtomCount(This);
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    float Distance;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    RepresentativeNode* NewNode;
    Representative* NewRep;
    PtrajCluster* MergedCluster;
    int OldBestRep, NewBestRep;
    int C2Ccentroid = This->action->darg4;
    int verbose = This->action->darg2;
    float DBI, pSF;
    int finish = 0;
    
    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);*/
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;

        /* Add representative point: */
        NewRep = RepresentativeNew(PointIndex,AtomCount);
        NewNode = (RepresentativeNode*)SafeMalloc(__FILE__, __LINE__, sizeof(RepresentativeNode));
        NewNode->pNext = NULL;
        NewNode->pPrev = NULL;
        NewNode->Point = NewRep;
        NewCluster->Head = NewNode;
        NewCluster->Tail = NewNode;        
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining clusters until you're finished */
    while (!finish)
    {
        /*fprintf(stdout,"%d,",This->ClusterCount);*/
        ClusteringOutputStatus( (Clustering *) This);
		if (This->ClusterCount < 2)
		{
	  	  finish = 1;
          continue;
		}
        /* Look at each pair of clusters, to find the closet pair */
        ClosestLength = 0;
        for (NodeA = This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            for (NodeB = NodeA->pNext,NodeIndexB=NodeIndexA+1; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
            {
                C2CDistance = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance < ClosestLength || ClosestLength == 0)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
        /* If they are "too far", stop agglomerating now */
        if (!DesiredClusterCount && C2CDistance>Epsilon)
        {
	  	  finish = 1;
          continue;
        }
	    if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
      	{
	  fprintf(stdout, "Merge cluster %d to cluster %d.\n", 
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1), 
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1)); 
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        OldBestRep = PtrajClusterGetBestRep( (PtrajCluster *) MergeNodeB->Cluster);
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
	  DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
	  pSF = ClusteringComputePseudoF( (Clustering *) This);
	  fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
	} else {
	  fprintf(MergeFile, "\n");
        }
        
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        MergedCluster = (PtrajCluster*)MergeNodeB->Cluster;
        ClusterFindCentroid( (Cluster *) MergedCluster);
        NewBestRep = ClusterFindBestP2C(MergedCluster);
	    if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
      	{
        	PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
		/* Done in ClusteringMergeClusters().
        if (C2Ccentroid > 0) 
    	{
	       if (OldBestRep != NewBestRep)
       	   {
	    	  ClusterAlignToBestRep(MergedCluster);
       	   }
	   	   ClusterFindCentroid(MergedCluster);
   		}*/
        if (This->ClusterCount<=DesiredClusterCount)
        {
	  	  finish = 1;
          continue;
        }

		PtrajClusterSelectRepresentatives(This, MergedCluster, Representatives);

        /* Move the representatives to their new locations */
        PtrajClusterRepositionRepresentatives(MergedCluster);
		
        /* Trim the list of representatives to the best */
        /*PtrajClusterPruneRepresentatives(MergedCluster,Representatives);       */

        /* Update cookies */
        for (NodeIndexA=ClosestNodeIndexA; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            
	  /*if (NodeIndexA>=ClosestNodeIndexA)*/
            ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
        }
        if (ClosestNodeIndexB>ClosestNodeIndexA)
        {
            ClosestNodeIndexB--;
        }

        /* Update cluster-to-cluster distances */
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            Distance = PtrajClusteringMinimumRepresentativeDistance(This,(PtrajCluster*)(NodeA->Cluster), (PtrajCluster*)(MergeNodeB->Cluster));
            SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],Distance);
        }
        
    }
    fclose(MergeFile);
    safe_free(ClusterDistances);
}

float PtrajClusteringMinimumRepresentativeDistance(PtrajClustering* This,PtrajCluster* ClusterA,PtrajCluster* ClusterB)
{
    float MinDistance = -1;
    float Distance;
    RepresentativeNode* NodeA;
    RepresentativeNode* NodeB;
    /* */
    for (NodeA = ClusterA->Head; NodeA; NodeA = NodeA->pNext)
    {
        for (NodeB = ClusterB->Head; NodeB; NodeB = NodeB->pNext)
        {
            Distance = PtrajClusteringRepToRepDistance(This,NodeA->Point,NodeB->Point);
            if (Distance<MinDistance || MinDistance<0)
            {
                MinDistance = Distance;
            }
        }
    }
    return MinDistance;
}


/* A CURE-like clustering algorithm.  Start by putting each point in its own teeny cluster, and making it the
representative of the cluster.  Then cycle through, merging a pair of clusters at each step of the cycle,
until the number of clusters is as low as requested.

Deciding who to merge: Find the pair of cluster (A, B) such that A and B are as close as possible
where X and Y are in different clusters. The distance between cluster A and B is defined as that in the
Complete algorithm, which is the largest distance among the pairs of points from A and B. This is not same 
as in the Centripetal algorithm SWT define before, which is more like the Edge algorithm.

How to merge: All points from the two clusters are put into the new cluster.  Recompute the centroid.  Then choose
up to <Representatives> new representative points for the new cluster, and shrink them by <Alpha> toward the 
centroid.  
*/
void PtrajClusteringClusterCentripetalComplete(PtrajClustering* This,int DesiredClusterCount, float Epsilon, int Representatives)
{
    int* ClusterCookies; /* Maps cluster LIST index to cluster distance ARRAY index */
    int PointIndex;
    int PointCount = This->PointCount;
    int AtomCount = PtrajClusteringGetAtomCount(This);
    PtrajCluster* NewCluster;
    SymmetricMatrix* ClusterDistances;
    int NodeIndexA;
    int NodeIndexB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    float Distance;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    float C2CDistance;
    RepresentativeNode* NewNode;
    Representative* NewRep;
    PtrajCluster* MergedCluster;
    int OldBestRep, NewBestRep;
    int C2Ccentroid = This->action->darg4;
    int verbose = This->action->darg2;
    float DBI, pSF;
    int finish = 0;
    /* */
    ClusterCookies = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    ClusterDistances = SymmetricMatrixCopy(This->PairwiseDistances);

    /* Build initial clusters */
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);*/
        PtrajClusterAddMember(NewCluster,PointIndex);
        NewCluster->IntName = PointIndex;
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        ClusterCookies[PointIndex]=PointIndex;

        /* Add representative point: */
        NewRep = RepresentativeNew(PointIndex,AtomCount);
        NewNode = (RepresentativeNode*)SafeMalloc(__FILE__, __LINE__, sizeof(RepresentativeNode));
        NewNode->pNext = NULL;
        NewNode->pPrev = NULL;
        NewNode->Point = NewRep;
        NewCluster->Head = NewNode;
        NewCluster->Tail = NewNode;        
    }

    FILE* MergeFile;
    MergeFile = fopen(ClusterMergingFile,"w");
    int MergeIntName = -2;  /* -1 is uninitialized. */
    /* Keep joining clusters until you're finished */
    while (!finish)
    {
        /*fprintf(stdout,"%d,",This->ClusterCount);*/
        ClusteringOutputStatus( (Clustering *) This);
		if (This->ClusterCount < 2)
		{
	  	  finish = 1;
          continue;
		}
        /* Look at each pair of clusters, to find the closet pair */
        ClosestLength = 0;
        for (NodeA = This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            for (NodeB = NodeA->pNext,NodeIndexB=NodeIndexA+1; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
            {
                C2CDistance = GetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[NodeIndexB]);
                if (C2CDistance < ClosestLength || ClosestLength == 0)
                {
                    MergeNodeA=NodeA;
                    MergeNodeB=NodeB;
                    ClosestLength = C2CDistance;
                    ClosestNodeIndexA=NodeIndexA;
                    ClosestNodeIndexB=NodeIndexB;
                }
            }
        }
        /* If they are "too far", stop agglomerating now */
        if (!DesiredClusterCount && C2CDistance>Epsilon)
        {
	  	  finish = 1;
          continue;
        }
	    if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
      	{
	  fprintf(stdout, "Merge cluster %d to cluster %d.\n", 
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1), 
		  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1)); 
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	  PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
        OldBestRep = PtrajClusterGetBestRep( (PtrajCluster *) MergeNodeB->Cluster);
        fprintf(MergeFile, "%d:\t%d\t%d\t%f", MergeIntName, 
		((PtrajCluster *)MergeNodeA->Cluster)->IntName, 
		((PtrajCluster *)MergeNodeB->Cluster)->IntName, ClosestLength);
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
        if (This->ClusterCount <= 50) {
        	DBI = ClusteringComputeDBI( (Clustering *) This, NULL);
        	pSF = ClusteringComputePseudoF( (Clustering *) This);
    	    fprintf(MergeFile, "\t%f\t%f\n", DBI, pSF);
  		} else {
        	fprintf(MergeFile, "\n");
        }
        
        ((PtrajCluster *)MergeNodeB->Cluster)->IntName = MergeIntName;
        MergeIntName--;
        MergedCluster = (PtrajCluster*)MergeNodeB->Cluster;
        ClusterFindCentroid( (Cluster *) MergedCluster);
        NewBestRep = ClusterFindBestP2C(MergedCluster);
	    if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
      	{
        	PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
        }
		/* Done in ClusteringMergeClusters().
        if (C2Ccentroid > 0) 
    	{
	       if (OldBestRep != NewBestRep)
       	   {
	    	  ClusterAlignToBestRep(MergedCluster);
       	   }
	   	   ClusterFindCentroid(MergedCluster);
   		}*/
        if (This->ClusterCount<=DesiredClusterCount)
        {
	  	  finish = 1;
          continue;
        }

		PtrajClusterSelectRepresentatives(This, MergedCluster, Representatives);

        /* Move the representatives to their new locations */
        PtrajClusterRepositionRepresentatives(MergedCluster);
		
        /* Trim the list of representatives to the best */
        /*PtrajClusterPruneRepresentatives(MergedCluster,Representatives);       */

        /* Update cookies */
        for (NodeIndexA=ClosestNodeIndexA; NodeIndexA<This->ClusterCount; NodeIndexA++)
        {
            
	  /*if (NodeIndexA>=ClosestNodeIndexA)*/
            ClusterCookies[NodeIndexA] = ClusterCookies[NodeIndexA+1];
        }
        if (ClosestNodeIndexB>ClosestNodeIndexA)
        {
            ClosestNodeIndexB--;
        }

        /* Update cluster-to-cluster distances */
        for (NodeA=This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
            Distance = PtrajClusteringMaximumRepresentativeDistance(This,(PtrajCluster*)(NodeA->Cluster), (PtrajCluster*)(MergeNodeB->Cluster));
            SetSymElement(ClusterDistances,ClusterCookies[NodeIndexA],ClusterCookies[ClosestNodeIndexB],Distance);
        }
        
    }
    fclose(MergeFile);
    safe_free(ClusterDistances);
}

float PtrajClusteringMaximumRepresentativeDistance(PtrajClustering* This,PtrajCluster* ClusterA,PtrajCluster* ClusterB)
{
    float MaxDistance = -1;
    float Distance;
    RepresentativeNode* NodeA;
    RepresentativeNode* NodeB;
    /* */
    for (NodeA = ClusterA->Head; NodeA; NodeA = NodeA->pNext)
    {
        for (NodeB = ClusterB->Head; NodeB; NodeB = NodeB->pNext)
        {
            Distance = PtrajClusteringRepToRepDistance(This,NodeA->Point,NodeB->Point);
            if (Distance>MaxDistance || MaxDistance<0)
            {
                MaxDistance = Distance;
            }
        }
    }
    return MaxDistance;
}


trajectoryInfo* PtrajClusteringReadDecoy(PtrajClustering* This, char* filename, int ClusterCount)
{
    coordinateInfo *DecoyCoordInfo;
    actionInformation *DecoyAction;
    ptrajState* DecoyState;
    int readCoordinates, processCoordinates;
    int i; 
    double box[6], boxnew[6], boxfixed[6];
    double *X, *Y, *Z;
    trajectoryInfo* DecoyTrajInfo;
    
    DecoyAction = ptrajCopyAction(&(This->action));
    DecoyState = This->action->state; /* See TRANSFORM_CLUSTER for DecoyAction's assignments of parameters. */
    DecoyCoordInfo = checkCoordinates(filename, DecoyState->atoms);
    if (ptrajPreprocessInputCoordinates(DecoyCoordInfo)|| DecoyCoordInfo == NULL) {
    	fprintf(stdout, "Cannot read Decoy coordinates!\n");
        return 0;
    }
    if (DecoyCoordInfo->stop != ClusterCount) {
    	fprintf(stdout, "Warning: Cluster count (%d) is not consistent with Decoy coordinates(%d frames)!\n", ClusterCount, DecoyCoordInfo->stop);
        DecoyAction->iarg2 = DecoyCoordInfo->stop;
        This->action->iarg2 = DecoyCoordInfo->stop;
    	fprintf(stdout, "         Cluster count is set to %d !\n", DecoyCoordInfo->stop);
    }
    DecoyAction->iarg3 = 0; /* No sieve for decoy trajectory. */
    DecoyAction->iarg4 = DecoyAction->iarg2; 
    /* The user should make sure that the trajin files and the decoy files use same prmtop file, 
       especially in this case, have same box info. */
    if (DecoyState->IFBOX) {
    	for (i=0; i<6; i++) {
      		box[i] = DecoyState->box[i];
      		boxfixed[i] = DecoyState->boxfixed[i];
    	}
    } else {
	    box[0] = 0; box[1] = 0; box[2] = 0;
	    box[3] = 90; box[4] = 90; box[5] = 90;
	    boxfixed[0] = 0; boxfixed[1] = 0; boxfixed[2] = 0;
	    boxfixed[3] = 0; boxfixed[4] = 0; boxfixed[5] = 0;
    }
	readCoordinates = 1;
    processCoordinates=0;
    X = (double *) SafeMalloc(__FILE__, __LINE__, sizeof(double) * DecoyState->atoms);
    Y = (double *) SafeMalloc(__FILE__, __LINE__, sizeof(double) * DecoyState->atoms);
    Z = (double *) SafeMalloc(__FILE__, __LINE__, sizeof(double) * DecoyState->atoms);
    
    
    DecoyTrajInfo = SafeMalloc(__FILE__, __LINE__, sizeof(trajectoryInfo));
    INITIALIZE_trajectoryInfo(DecoyTrajInfo);
    DecoyTrajInfo->state = NULL;
    DecoyTrajInfo->state = ((trajectoryInfo*)DecoyAction->carg2)->state;
    DecoyTrajInfo->atoms = DecoyTrajInfo->state->atoms;
    
    int set = 1;
    while (readCoordinates) {
        for (i=0; i<6; i++) boxnew[i] = box[i];
    	ptrajProcessInputCoordinates(DecoyCoordInfo, DecoyState, X, Y, Z, boxnew, set++, &readCoordinates, &processCoordinates);
	   	if (readCoordinates) fprintf(stdout, "  Read in decoy structure %d.\n", set-1);
        for (i=0; i<6; i++)	if (boxfixed[i] == 0)  box[i] = boxnew[i];
        if (readCoordinates) AccumulateCoordinates(DecoyTrajInfo, DecoyAction, X, Y, Z, NULL);
    }
    return DecoyTrajInfo;
}

void PtrajClusteringAssignDecoyToCentroid(PtrajClustering *This, int DesiredClusterCount, trajectoryInfo *DecoyTrajInfo)
{
    int ClusterIndex, i;
    PtrajCluster* NewCluster;
    ClusterNode* Node;
    
    for (ClusterIndex=0; ClusterIndex<DesiredClusterCount; ClusterIndex++)
    {
        NewCluster = (PtrajCluster *) ClusteringGetNewCluster( (Clustering *) This);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        for (i = 0; i < DecoyTrajInfo->atoms; i++) {
        	NewCluster->CentroidX[i] = DecoyTrajInfo->x[i + ClusterIndex * DecoyTrajInfo->atoms];
        	NewCluster->CentroidY[i] = DecoyTrajInfo->y[i + ClusterIndex * DecoyTrajInfo->atoms];
        	NewCluster->CentroidZ[i] = DecoyTrajInfo->z[i + ClusterIndex * DecoyTrajInfo->atoms];
        }
    }
	
}

coordinateInfo* CloneCoordinateInfo(coordinateInfo* Source, char* ExtraName)
{
    coordinateInfo* NewCoordinateInfo;
    /*int FileNameLength;*/
    /*char* Buffer;*/
    /* */
    NewCoordinateInfo = (coordinateInfo *) SafeMalloc(__FILE__, __LINE__, sizeof(coordinateInfo));
    INITIALIZE_coordinateInfo(NewCoordinateInfo);
    NewCoordinateInfo->isBox = Source->isBox;
    NewCoordinateInfo->type = Source->type;
    /*FileNameLength = strlen(Source->filename);
    FileNameLength += strlen(ExtraName)+1;
    NewCoordinateInfo->filename = SafeMalloc(__FILE__, __LINE__, sizeof(char)*FileNameLength);
    strcpy(NewCoordinateInfo->filename,Source->filename);
    strcat(NewCoordinateInfo->filename,".");
    strcat(NewCoordinateInfo->filename,ExtraName);*/
    return NewCoordinateInfo;
}

/* Allocate once, to save memory juking */
static coordinateInfo* SecondPassCoordinateInfo = NULL;

void PtrajClusteringSetupSecondPassInfo(PtrajClustering* This, actionInformation* action)
{
    ClusterNode* Node;
    PtrajCluster* Cluster;
    int ClusterIndex;
    int AtomCount = PtrajClusteringGetAtomCount(This);
    trajectoryInfo* OutputTrajInfo;    
    int FullAtomCount; /* = action->state->atoms;*/
    /**/
    if (action->carg3)
    {
		OutputTrajInfo = (trajectoryInfo*)action->carg3;
    }
    else
    {
		OutputTrajInfo = (trajectoryInfo*)action->carg2;
    }
    FullAtomCount = OutputTrajInfo->atoms;
    for (Node = This->Head,ClusterIndex = 0; Node!=NULL; Node = Node->pNext,ClusterIndex++)
    {
        Cluster = (PtrajCluster*)Node->Cluster;
        /* Allocate space to store our average */
        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE] != CLUSTER_OUTPUT_NONE)
        {
            Cluster->AverageX = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * FullAtomCount);
            memset(Cluster->AverageX, 0, sizeof(float) * FullAtomCount);
            Cluster->AverageY = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * FullAtomCount);
            memset(Cluster->AverageY, 0, sizeof(float) * FullAtomCount);
            Cluster->AverageZ = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * FullAtomCount);
            memset(Cluster->AverageZ, 0, sizeof(float) * FullAtomCount);
        }

        /* Allocate space to store our representative */
        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_REPRESENTATIVE] != CLUSTER_OUTPUT_NONE)
        {
            Cluster->RepX = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * FullAtomCount);
            memset(Cluster->RepX, 0, sizeof(float) * FullAtomCount);
            Cluster->RepY = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * FullAtomCount);
            memset(Cluster->RepY, 0, sizeof(float) * FullAtomCount);
            Cluster->RepZ = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * FullAtomCount);
            memset(Cluster->RepZ, 0, sizeof(float) * FullAtomCount);
            Cluster->MinRepFrameDistance = -1.0;
        }
    }
}


/* Cope with zany: PtrajOutputCoordinates takes a "first" flag which determines whether it opens a 
new file, or appends to its already-open file.  We need to set the 'first' flag to true each time
we start a new cluster.  And, we need to save the file object references! */
static int* FirstOutputFlag = NULL;
FILE** OutputFiles = NULL;
static int* FileNumber = NULL;

#define FILE_SIZE_LIMIT 524288000
/*#define FILE_SIZE_LIMIT 1024*/

/*
There are some steps (like outputting average frames) that we do to a frame *after* it's clustered.
When we're sieving, and making a second pass through the frames, we do that processing *during* the 
second pass, while we still have the frame data available. 
Note: The arrays, dX,dY,dZ include all coordinates.  The arrays X,Y,Z include only coordinates
for the selected atoms.
*/
void PtrajClusteringSecondPassOutputInfo(PtrajClustering* This, actionInformation* action, 
    int FrameIndex, double* dX, double* dY, double* dZ, float* X, float* Y, float* Z)
{
    coordinateInfo* CoordinateInfo = NULL;
    coordinateInfo* CoordinateInfoSwap = NULL;
    coordinateInfo* outInfo = NULL;
    ClusterNode* Node;
    PtrajCluster* pCluster;
    int FirstFlag;
    int ClusterIndex;
    char Buffer[BUFFER_SIZE];    
    int AtomCount = This->trajInfo->atoms;
    trajectoryInfo* OutputTrajInfo;
    int AtomIndex;
    float Distance;
    /**/

    if (action->carg3)
    {
		OutputTrajInfo = (trajectoryInfo*)action->carg3;
    }
    else
    {
		OutputTrajInfo = (trajectoryInfo*)action->carg2;
    }

    if (FirstOutputFlag==NULL)
    {
      FirstOutputFlag = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * This->ClusterCount);
      memset(FirstOutputFlag, 1, sizeof(int) * This->ClusterCount);
      OutputFiles = (FILE**)SafeMalloc(__FILE__, __LINE__, sizeof(FILE*) * This->ClusterCount);
      memset(OutputFiles, 0, sizeof(FILE*) * This->ClusterCount);
      FileNumber = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * This->ClusterCount);
      memset(FileNumber, 0, sizeof(int) * This->ClusterCount);
    }

    /* Temporarily switch the value of the global variable outInfo */
    CoordinateInfoSwap = outInfo;
    if (!SecondPassCoordinateInfo)
    {
        SecondPassCoordinateInfo = (coordinateInfo *) SafeMalloc(__FILE__, __LINE__, sizeof(coordinateInfo));
        INITIALIZE_coordinateInfo(SecondPassCoordinateInfo);
        SecondPassCoordinateInfo->isBox = 0; 
        SecondPassCoordinateInfo->type = ((coordType*)action->carg4)[CLUSTER_FILEFORMAT_ALL];
    }
    outInfo = SecondPassCoordinateInfo;
        
    for (Node = This->Head,ClusterIndex = 0; Node!=NULL; Node = Node->pNext,ClusterIndex++)
    {
        pCluster = (PtrajCluster*)Node->Cluster;
        if (!ClusterIsMember( (Cluster *) pCluster,FrameIndex))
        {
            continue;
        }

          /*
	   *  Output this frame in this cluster
	   */
        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_ALL] != CLUSTER_OUTPUT_NONE)
        {
	  FirstFlag = FirstOutputFlag[ClusterIndex];
	  if (!FirstFlag)
	    {
	      outInfo->file = OutputFiles[ClusterIndex];
	           /*
		    *  fprintf(stdout, "We're at file position %d\n",ftell(outInfo->file));
		    */
	      if (ftell(outInfo->file) > FILE_SIZE_LIMIT)
	        {
		    /*
		     *  This is a *LARGE* trajectory.  Split the cluster into separate files, 
		     *  to avoid crashing out if we hit the OS-level file-size limit: 
		     */
		  FileNumber[ClusterIndex] = FileNumber[ClusterIndex] + 1;
		  safe_fclose(OutputFiles[ClusterIndex]);
		  FirstOutputFlag[ClusterIndex] = 1;
		  FirstFlag = 1;      
	        }
	    }

	  if (FileNumber[ClusterIndex])
	    sprintf(Buffer, "%s.c%i.%i", action->carg1, ClusterIndex, FileNumber[ClusterIndex]);
	  else
	    sprintf(Buffer, "%s.c%i", action->carg1, ClusterIndex);

	  SecondPassCoordinateInfo->filename = Buffer;
	    

	    /*
	     *  fprintf(stdout,"Output frame %d to cluster file %s firstflag %d\n",FrameIndex, Buffer, FirstFlag);
	     */

	  ptrajOutputCoordinates(SecondPassCoordinateInfo, 
				 This->trajInfo->state, 
				 This->trajInfo->current, 
				 SecondPassCoordinateInfo->append, FirstFlag, 
				 SecondPassCoordinateInfo->append, 
				 OutputTrajInfo->atoms, 
				 dX,dY,dZ, 
				 This->trajInfo->state->box);
	  
	  OutputFiles[ClusterIndex] = outInfo->file;
	  FirstOutputFlag[ClusterIndex] = 0;
	}
 
          /*
	   *  Accumulate averages, if requested 
	   */
        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE] != CLUSTER_OUTPUT_NONE)
        {
	  for (AtomIndex=0; AtomIndex < OutputTrajInfo->atoms; AtomIndex++)
            {
	      pCluster->AverageX[AtomIndex] += dX[AtomIndex];
	      pCluster->AverageY[AtomIndex] += dY[AtomIndex];
	      pCluster->AverageZ[AtomIndex] += dZ[AtomIndex];
            }
        }


          /*
	   *  Find 'representative', if requested 
	   */
        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_REPRESENTATIVE] != CLUSTER_OUTPUT_NONE)
        {
	  Distance = PtrajGetDistance(This, -1, -1, 
				      pCluster->CentroidX,pCluster->CentroidY,pCluster->CentroidZ,
				      X,Y,Z);
	  if (pCluster->MinRepFrameDistance < 0 || Distance < pCluster->MinRepFrameDistance)
            {
	      pCluster->MinRepFrameDistance = Distance;
	      pCluster->MinRepFrameIndex = FrameIndex;
	      for (AtomIndex=0; AtomIndex < OutputTrajInfo->atoms; AtomIndex++)
		{
		  pCluster->RepX[AtomIndex] = (float)dX[AtomIndex];
		  pCluster->RepY[AtomIndex] = (float)dY[AtomIndex];
		  pCluster->RepZ[AtomIndex] = (float)dZ[AtomIndex];
		}
            }
        }
    }

    outInfo = CoordinateInfoSwap;

}



void PtrajClusteringChooseSievePoints(PtrajClustering* This, int FullPointCount)
{
  int PointsToChoose;
  float RandomValue;
  int PointIndex;
  /**/
  This->SieveFirstPass = calloc(FullPointCount, sizeof(int));
  PointIndex = (int)This->action->darg3;


  if ( PointIndex >= 0) 
  {
  	for (;PointIndex < FullPointCount; PointIndex += This->action->iarg3)
    {
      This->SieveFirstPass[PointIndex] = 1;	
      /* printf("PointIndex is %i\n", PointIndex);*/
    }
	This->action->iarg4 = (PointIndex - This->action->darg3) /  This->action->iarg3;
    return;
  }
  
  /*  randomly picking */
  srand(time(NULL) + clock());

  PointsToChoose = FullPointCount / This->action->iarg3;
  This->PointCount = PointsToChoose;
  This->action->iarg4 = PointsToChoose;
  PointsToChoose = max(PointsToChoose, 1);
  while (PointsToChoose)
  {
    RandomValue = rand() / (double)RAND_MAX;
    PointIndex = (int)(RandomValue * FullPointCount);
    PointIndex = min(FullPointCount-1, PointIndex);
    while (This->SieveFirstPass[PointIndex])
    {
      PointIndex++;
      if (PointIndex >= FullPointCount)
      {
		PointIndex = 0;
      }
    }
    This->SieveFirstPass[PointIndex] = 1;
    PointsToChoose--;
  }
  This->SievedIndex = 0;
}


void PtrajClusteringSecondPass(actionInformation* action,double* x, double* y, double* z)
{
    PtrajClustering* This;
    trajectoryInfo* trajInfo;
    float* fX = NULL;
    float* fY = NULL;
    float* fZ = NULL;
    int AtomIndex;
    int ValidIndex = 0;
    /* */
    This = (PtrajClustering*)action->carg5;
    trajInfo = (trajectoryInfo *) action->carg2;
    /* The x,y,z double arrays are the coordinates for ALL the atoms in this frame.  For most of the 
       clustering code, we want to have floats, and we want floats for just the selected atoms. 
       So, construct (possibly shorter-length) fX,fY,fZ arrays. */
    fX = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*This->trajInfo->atoms);
    fY = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*This->trajInfo->atoms);
    fZ = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*This->trajInfo->atoms);

    for (AtomIndex=0; AtomIndex<action->state->atoms; AtomIndex++)
      {
		if (action->mask[AtomIndex])
		  {
	        fX[ValidIndex] = (float)x[AtomIndex];
    	    fY[ValidIndex] = (float)y[AtomIndex];
	   	    fZ[ValidIndex] = (float)z[AtomIndex];
			ValidIndex++;
		  }
	  }
    /* If we're clustering, then all the frames except for every n-th frame will need to be added
    to a cluster now: */
    /*PtrajClusteringAddFrameToBestCluster(This,trajInfo->current,fX,fY,fZ);*/
    /*if (trajInfo->current%action->iarg3)*/
    if (!This->SieveFirstPass[trajInfo->current])
    {
        PtrajClusteringAddFrameToBestCluster(This,trajInfo->current,fX,fY,fZ);
    }
    PtrajClusteringSecondPassOutputInfo(This, action, trajInfo->current, x, y, z, fX, fY, fZ);
    safe_free(fX);
    safe_free(fY);
    safe_free(fZ);

    trajInfo->current += 1;
}

void OutputNewickTree(PtrajClustering* This, FILE* OutputFile)
{
  ClusterNode* Node;
  PtrajCluster* Cluster;
  Cluster = (PtrajCluster*)This->Head->Cluster;
  if (!Cluster->Name)
  {
    return;
  }
  fprintf(OutputFile,"#(");
  for (Node = This->Head; Node; Node = Node->pNext)
  {
    if (Node != This->Head)
    {
      fprintf(OutputFile,",");
    }
    Cluster = (PtrajCluster*)Node->Cluster;
    fprintf(OutputFile,"%s",Cluster->Name);
  }
  fprintf(OutputFile,")\n");
}


void OutputClusteringStats(PtrajClustering* This, FILE* File)
{
    actionInformation* action;
    float DBI = 0.0;
    float pSF = 0.0;
    action = This->action;
    int SavedC2C;
    float SSRSST;
    int PointCount = This->PointCount;
    int ClusterCount = This->ClusterCount;
    
    /* ONLY compute DBI if we aren't sieving.  (DBI requires another pass through the data) */
    if (!action->iarg3)
    {
        SavedC2C = action->darg4;
        action->darg4 = 1;         /* When calculate the DBI or pSF, always use centroid, not bestrep. */
        DBI = ClusteringComputeDBI( (Clustering *) This, File);
	pSF = ClusteringComputePseudoF( (Clustering *) This);
        action->darg4 = SavedC2C;
    }
    fprintf(File,"#DBI: %.5f\n",DBI);
    fprintf(File,"#pSF: %.5f\n",pSF);
    SSRSST = pSF*(ClusterCount-1)/(PointCount-ClusterCount+pSF*(ClusterCount-1));
    fprintf(File,"#SSR/SST: %.5f\n",SSRSST);
    /* If there are names, output our tree: */
    OutputNewickTree(This, File);
}

/*
Write header txt file for clustering, such as algorithm, clusters, time,
*/
void PtrajClusteringOutputHeaderToFile(PtrajClustering* This, FILE* File)
{
    actionInformation* action;
	time_t t;
    
    action = This->action;
    fprintf(File,"#Desired Clusters: %d\n",action->iarg2);
    fprintf(File,"#Clustering Algorithm: %s\n",CLUSTER_ALGORITHM_NAMES[action->iarg1]);
    fprintf(File,"#Distance Metric: %s\n",DISTANCE_METRIC_NAMES[action->iarg5]);
    fprintf(File,"#Matrix Calculation Time: %d seconds\n",This->MatrixTime-This->StartTime);
    fprintf(stdout,"  Clustering elapsed time: %d seconds\n",time(NULL) - This->StartTime);
    if (action->iarg3 > 0)
	    fprintf(File,"#First Pass Time: %d seconds\n",This->FirstPassTime - This->MatrixTime);
    t = time(NULL);
    fprintf(File,"#Total Calculation Time: %d seconds\n",t - This->StartTime);

}

/* Align the FullTrajInfo to its trajInfo position, which may be changed during clustering */
void AlignFullTraj(PtrajClustering* This)
{
    ClusterNode* Node;
    PtrajCluster* TheCluster;
    PtrajCluster* MaxCluster;
    int ClusterIndex;
    int Points;    
    int MaxPoints = 0;  
    int PointIndex;

    float rmsRotation[3][3], rmsTranslation[3];
    trajectoryInfo* trajInfo = This->trajInfo;
    trajectoryInfo* OutputTrajInfo;
    int i, j;
    float xtemp, ytemp, ztemp;
    float* XA;
    float* YA;
    float* ZA;
    float* Xref;
    float* Yref;
    float* Zref;
    int mass = This->action->iarg6;
    
	actionInformation* action = This->action;
        
    if (action->carg3)
    {
        OutputTrajInfo = (trajectoryInfo*)action->carg3;
    }
    else
    {
        OutputTrajInfo = (trajectoryInfo*)action->carg2;
        return; /* TrajInfo is already full-atom-selection */
    }
    trajInfo = (trajectoryInfo*)action->carg2;

    XA = (float *)SafeMalloc(__FILE__, __LINE__, trajInfo->atoms * sizeof(float));
    YA = (float *)SafeMalloc(__FILE__, __LINE__, trajInfo->atoms * sizeof(float));
    ZA = (float *)SafeMalloc(__FILE__, __LINE__, trajInfo->atoms * sizeof(float));
    for (PointIndex=0; PointIndex<ClusteringGetPointCount( (Clustering *) This); PointIndex++)
    {
		for (i=0,j=0; i<OutputTrajInfo->atoms;i++)
		{
	        if (action->mask[i])
       		{
				XA[j] = OutputTrajInfo->x[PointIndex * OutputTrajInfo->atoms + i];
				YA[j] = OutputTrajInfo->y[PointIndex * OutputTrajInfo->atoms + i];
				ZA[j] = OutputTrajInfo->z[PointIndex * OutputTrajInfo->atoms + i];
                j++;
            }
		}
	    Xref = &trajInfo->x[PointIndex * trajInfo->atoms]; 
	    Yref = &trajInfo->y[PointIndex * trajInfo->atoms]; 
   		Zref = &trajInfo->z[PointIndex * trajInfo->atoms]; 
        float fitX, fitY, fitZ;
        fitX = XA[0];
   	    fitY = YA[0];
       	fitZ = ZA[0];
	    if (mass)
        {
	    	rmsf(trajInfo->atoms,1,trajInfo->state->masses,NULL,XA,YA,ZA,Xref,Yref,Zref,rmsRotation,rmsTranslation,1);
	    }
        else
        {
           	rmsf(trajInfo->atoms,1,NULL,NULL,XA,YA,ZA,Xref,Yref,Zref,rmsRotation,rmsTranslation,1);
   		}
        fitX -= XA[0];
        fitY -= YA[0];
        fitZ -= ZA[0];
	    for (i=0; i<OutputTrajInfo->atoms;i++)
	    {
   			OutputTrajInfo->x[PointIndex * OutputTrajInfo->atoms + i] -= fitX;
   			OutputTrajInfo->y[PointIndex * OutputTrajInfo->atoms + i] -= fitY;
   			OutputTrajInfo->z[PointIndex * OutputTrajInfo->atoms + i] -= fitZ;
            VOP_3x3_TIMES_COORDS(rmsRotation, OutputTrajInfo->x[PointIndex * OutputTrajInfo->atoms + i],OutputTrajInfo->y[PointIndex * OutputTrajInfo->atoms + i],OutputTrajInfo->z[PointIndex * OutputTrajInfo->atoms + i], xtemp, ytemp, ztemp);    
	    	OutputTrajInfo->x[PointIndex * OutputTrajInfo->atoms + i] += rmsTranslation[0];
	    	OutputTrajInfo->y[PointIndex * OutputTrajInfo->atoms + i] += rmsTranslation[1];
	    	OutputTrajInfo->z[PointIndex * OutputTrajInfo->atoms + i] += rmsTranslation[2];
	    }
    }
    ClusteringFindAllCentroids( (Clustering *) This);
    safe_free(XA);
    safe_free(YA);
    safe_free(ZA);
}


/* Align all the BestRep of each cluster to the most populated cluster's BestRep for purpose of good centroid of all frames. 
   Since this is the final steps in a clustering algorithm, centroid of each cluster will be updated. */
void AlignBestReps(PtrajClustering* This)
{
    ClusterNode* Node;
    PtrajCluster* TheCluster;
    PtrajCluster* MaxCluster;
    int ClusterIndex;
    int Points;    
    int MaxPoints = 0;  
    int i;  
    int PointIndex;
    float rms;
    
    
    for (Node = This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
    {
        TheCluster = (PtrajCluster*)Node->Cluster;
	Points = ClusterCountPointsInCluster( (Cluster *) TheCluster);
        if (Points > MaxPoints)
        {
        	MaxPoints = Points;
            MaxCluster = TheCluster;
            i = ClusterIndex;
        }
    }

    fprintf(stdout, "  Aligning the most representative frames from each clusters to cluster %d.\n", i);

    for (Node = This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
    {
        TheCluster = (PtrajCluster*)Node->Cluster;
        if (TheCluster == MaxCluster) continue;
        ClusteringAlignAToB(This, TheCluster->BestRep, MaxCluster->BestRep);
        for (PointIndex=0; PointIndex<ClusteringGetPointCount( (Clustering *) This); PointIndex++)
        {
	        if (ClusterIsMember( (Cluster *) TheCluster,PointIndex))
    	    {
        		rms = ClusteringAlignAToB(This, PointIndex, TheCluster->BestRep);
            }
        }
        ClusterFindCentroid( (Cluster *) TheCluster);
    }
    
	AlignFullTraj(This);    
}


/*
Perform end-of-pass operations on our clusters.
*/
void PtrajClusteringOutputClusterInfo(PtrajClustering* This,actionInformation* action)
{
    ClusterNode* Node;
    PtrajCluster* pCluster;
    FILE* OutputFile;
    int FrameIndex;
    int ClusterIndex = 0;
    double* dX;
    double* dY;
    double* dZ;
    float* X;
    float* Y;
    float* Z;
        
    double* AverageX;
    double* AverageY;
    double* AverageZ;
	float *avgX, *avgY, *avgZ;
    
    char Buffer[BUFFER_SIZE];        
    coordinateInfo* CoordinateInfo = NULL;
    coordinateInfo* CoordinateInfoAverage = NULL;
    coordinateInfo* CoordinateInfoRepresentative = NULL;
    coordinateInfo* CoordinateInfoSwap = NULL;    
    coordinateInfo* outInfo = NULL;    
    int SetIndex;
    int FirstFlag;
    int AtomIndex;
    float ClosestDistance;
    float Distance;
    int ClosestFrameIndex;
    char FilePath[1024];
    FILE* File;
    float DBI = 0.0;
    float pSF = 0.0;
    trajectoryInfo* OutputTrajInfo; 
    trajectoryInfo* TrajInfo; 
    /* */  
    
    if (action->carg3)
        OutputTrajInfo = (trajectoryInfo*)action->carg3;
    else
        OutputTrajInfo = (trajectoryInfo*)action->carg2;
        
    TrajInfo = (trajectoryInfo*)action->carg2;
   
    strcpy(FilePath,(char*)action->carg1);
    strcat(FilePath,".txt");
    File = fopen(FilePath,"w");
    /*
    fprintf(stdout,"Sorting Clusters ..\n");
    ClusteringSortCluster(This);
    */
    /*AlignBestReps(This);*/
    PtrajClusteringOutputHeaderToFile(This, File);
    OutputClusteringStats(This, File); /* get the DBI and pSF stuff. */
    printf("  Printing distribution of distances to file %s\n", FilePath);
    if (This->PairwiseDistances) SymmetricMatrixDistribution(This->PairwiseDistances,20,File);
    fclose(File);
    ClusteringAppendClusterList( (Clustering *) This,FilePath);
    
    /* Set up the CoordinateInfo structure for use later */
    CoordinateInfo = (coordinateInfo *) SafeMalloc(__FILE__, __LINE__, sizeof(coordinateInfo));
    INITIALIZE_coordinateInfo(CoordinateInfo);
    CoordinateInfo->isBox = 0; 
    CoordinateInfo->filename = (char*)action->carg1;
    CoordinateInfo->type = ((coordType*)action->carg4)[CLUSTER_FILEFORMAT_ALL];

    if (((int*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE] != CLUSTER_OUTPUT_NONE)
    {
        AverageX = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);
        AverageY = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);
        AverageZ = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);
        CoordinateInfoAverage = CloneCoordinateInfo(CoordinateInfo,"avg");
        CoordinateInfoAverage->type = ((coordType*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE];
    }
    if (((int*)action->carg4)[CLUSTER_FILEFORMAT_REPRESENTATIVE] != CLUSTER_OUTPUT_NONE)
    {
        CoordinateInfoRepresentative = CloneCoordinateInfo(CoordinateInfo,"rep");
        CoordinateInfoRepresentative->type = ((coordType*)action->carg4)[CLUSTER_FILEFORMAT_REPRESENTATIVE];
    }

    /* If we performed a second pass, then we've already done most of the processing: */
    if (action->iarg3)
    {
        /* Output representatives and averages */
        dX = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);
        dY = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);
        dZ = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);
        CoordinateInfoSwap = outInfo;
        for (Node = This->Head,ClusterIndex=0; Node!=NULL; Node = Node->pNext,ClusterIndex++)
        {
            pCluster = (PtrajCluster*)Node->Cluster;
            if (((int*)action->carg4)[CLUSTER_FILEFORMAT_REPRESENTATIVE] != CLUSTER_OUTPUT_NONE)
            {
                for (AtomIndex=0; AtomIndex<action->state->atoms; AtomIndex++)
                {
                    dX[AtomIndex] = pCluster->RepX[AtomIndex];
                    dY[AtomIndex] = pCluster->RepY[AtomIndex];
                    dZ[AtomIndex] = pCluster->RepZ[AtomIndex];
                }
                sprintf(Buffer, "%s.rep.c%i", action->carg1, ClusterIndex);
                CoordinateInfoRepresentative->filename = Buffer;
                outInfo = CoordinateInfoRepresentative;
                ptrajOutputCoordinates(outInfo,action->state, -1, 0, 1, 0, action->state->atoms, dX,dY,dZ, action->state->box);
            }
            if (((int*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE] != CLUSTER_OUTPUT_NONE)
            {
                int count = ClusterCountPointsInCluster( (Cluster *) pCluster);
                for (AtomIndex=0; AtomIndex<action->state->atoms; AtomIndex++)
                {
                    AverageX[AtomIndex] = pCluster->AverageX[AtomIndex] / count;
                    AverageY[AtomIndex] = pCluster->AverageY[AtomIndex] / count;
                    AverageZ[AtomIndex] = pCluster->AverageZ[AtomIndex] / count;
                }
                sprintf(Buffer, "%s.avg.c%i", action->carg1, ClusterIndex);
                CoordinateInfoAverage->filename = Buffer;
                outInfo = CoordinateInfoAverage;
                ptrajOutputCoordinates(outInfo, action->state, -1, 0, 1, 0, action->state->atoms, AverageX,AverageY,AverageZ, action->state->box);
            }
        }
        safe_free(dX);
        safe_free(dY);
        safe_free(dZ);
        outInfo = CoordinateInfoSwap;
        /* Done!  The rest of our processing occurred while we were cycling through the frames */
        if (CoordinateInfoAverage)
        {
            safe_free(CoordinateInfoAverage);
        }
        if (CoordinateInfoRepresentative)
        {
            safe_free(CoordinateInfoRepresentative);
        }
        if (CoordinateInfo)
        {
            safe_free(CoordinateInfo);
        }

        return;
    }
    
    /*
     * If second pass is not required, print out average structure or representative structure or cluster trajectory.
     */
    CoordinateInfoSwap = outInfo;

    dX = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);
    dY = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);
    dZ = (double*)SafeMalloc(__FILE__, __LINE__, sizeof(double) * OutputTrajInfo->atoms);


    for (Node = This->Head; Node!=NULL; Node = Node->pNext)
    {
        pCluster = (PtrajCluster*)Node->Cluster;
        
	    avgX = ((PtrajCluster*)pCluster)->CentroidX;
	    avgY = ((PtrajCluster*)pCluster)->CentroidY;
	    avgZ = ((PtrajCluster*)pCluster)->CentroidZ;

        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE] != CLUSTER_OUTPUT_NONE)
        {
            memset(AverageX,0,sizeof(double) * OutputTrajInfo->atoms);
            memset(AverageY,0,sizeof(double) * OutputTrajInfo->atoms);
            memset(AverageZ,0,sizeof(double) * OutputTrajInfo->atoms);            
        }
        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_REPRESENTATIVE] != CLUSTER_OUTPUT_NONE)
        {            
            ClosestDistance = 0;
        }
        SetIndex = 0;
        for (FrameIndex=0; FrameIndex<pCluster->PointCount; FrameIndex++)
        {
            if (!ClusterIsMember( (Cluster *) pCluster,FrameIndex))
            {
                continue;
            }            
            if (action->iarg3 && This->SieveFirstPass[FrameIndex])
            {
                continue;
            }

            /* Output type 'all': Spit out each frame of the cluster individually */
            if (((int*)action->carg4)[CLUSTER_FILEFORMAT_ALL] != CLUSTER_OUTPUT_NONE)
            {
                /* Convert float-array to double-array. This is necessary because ptrajOutputCoordinates() want double array. */
                for (AtomIndex=0; AtomIndex<OutputTrajInfo->atoms; AtomIndex++)
                {
                    dX[AtomIndex] = (double)OutputTrajInfo->x[FrameIndex*OutputTrajInfo->atoms + AtomIndex];
                    dY[AtomIndex] = (double)OutputTrajInfo->y[FrameIndex*OutputTrajInfo->atoms + AtomIndex];
                    dZ[AtomIndex] = (double)OutputTrajInfo->z[FrameIndex*OutputTrajInfo->atoms + AtomIndex];
                }
                FirstFlag = (SetIndex==0);
                sprintf(Buffer, "%s.c%i", action->carg1, ClusterIndex);
                CoordinateInfo->filename = Buffer;

                outInfo = CoordinateInfo;

                ptrajOutputCoordinates(outInfo, 
                                       OutputTrajInfo->state, 
                                       SetIndex, 
                                       1-FirstFlag, /* if not first, then append */
                                       FirstFlag,   /* is first? */
                                       0,           /* is last? */    
                                       OutputTrajInfo->atoms, 
                                       dX,dY,dZ, OutputTrajInfo->state->box);
            } 
            /* Output type 'average': Accumulate coordinates over the cluster frames */
            if (((int*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE] != CLUSTER_OUTPUT_NONE)
            {
                for (AtomIndex=0; AtomIndex<OutputTrajInfo->atoms; AtomIndex++)
                {
                    AverageX[AtomIndex] += OutputTrajInfo->x[FrameIndex*OutputTrajInfo->atoms + AtomIndex];
                    AverageY[AtomIndex] += OutputTrajInfo->y[FrameIndex*OutputTrajInfo->atoms + AtomIndex];
                    AverageZ[AtomIndex] += OutputTrajInfo->z[FrameIndex*OutputTrajInfo->atoms + AtomIndex];
                }
            }
            /* Output type 'representative': Choose the frame closest to the centroid */
            if (((int*)action->carg4)[CLUSTER_FILEFORMAT_REPRESENTATIVE] != CLUSTER_OUTPUT_NONE)
            {
                X = TrajInfo->x + FrameIndex*TrajInfo->atoms;
                Y = TrajInfo->y + FrameIndex*TrajInfo->atoms;
                Z = TrajInfo->z + FrameIndex*TrajInfo->atoms;
                Distance = PtrajGetDistance(This, -1, -1, avgX,avgY,avgZ, X,Y,Z);
                if (ClosestDistance==0 || Distance < ClosestDistance)
                {
                    ClosestDistance = Distance;
                    ClosestFrameIndex = FrameIndex;
                }
            }

            SetIndex++;

        } /* loop on frames */

        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_REPRESENTATIVE] != CLUSTER_OUTPUT_NONE)
        {
            for (AtomIndex=0; AtomIndex<OutputTrajInfo->atoms; AtomIndex++)
            {
                dX[AtomIndex] = OutputTrajInfo->x[ClosestFrameIndex*OutputTrajInfo->atoms + AtomIndex];
                dY[AtomIndex] = OutputTrajInfo->y[ClosestFrameIndex*OutputTrajInfo->atoms + AtomIndex];
                dZ[AtomIndex] = OutputTrajInfo->z[ClosestFrameIndex*OutputTrajInfo->atoms + AtomIndex];
            }
            sprintf(Buffer, "%s.rep.c%i", action->carg1, ClusterIndex);
            CoordinateInfoRepresentative->filename = Buffer;

            outInfo = CoordinateInfoRepresentative;
            ptrajOutputCoordinates(outInfo, OutputTrajInfo->state, -1, 0, 1, 0, OutputTrajInfo->atoms, dX,dY,dZ, OutputTrajInfo->state->box);

        }
        if (((int*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE] != CLUSTER_OUTPUT_NONE)
        {
            for (AtomIndex=0; AtomIndex<OutputTrajInfo->atoms; AtomIndex++)
            {
                AverageX[AtomIndex] = AverageX[AtomIndex] / (SetIndex);
                AverageY[AtomIndex] = AverageY[AtomIndex] / (SetIndex);
                AverageZ[AtomIndex] = AverageZ[AtomIndex] / (SetIndex);
            }
            sprintf(Buffer, "%s.avg.c%i", action->carg1, ClusterIndex);
            CoordinateInfoAverage->filename = Buffer;
            outInfo = CoordinateInfoAverage;
            ptrajOutputCoordinates(outInfo, OutputTrajInfo->state, -1, 0, 1, 0, OutputTrajInfo->atoms, AverageX,AverageY,AverageZ, OutputTrajInfo->state->box);
        }
        ClusterIndex++;
    }

    /* Free data */
    if (((int*)action->carg4)[CLUSTER_FILEFORMAT_AVERAGE] != CLUSTER_OUTPUT_NONE)
    {
        safe_free(AverageX);
        safe_free(AverageY);
        safe_free(AverageZ);
    }    
    safe_free(dX);
    safe_free(dY);
    safe_free(dZ);

    if (CoordinateInfoAverage)
    {
        safe_free(CoordinateInfoAverage);
    }
    if (CoordinateInfoRepresentative)
    {
        safe_free(CoordinateInfoRepresentative);
    }
    if (CoordinateInfo)
    {
        safe_free(CoordinateInfo);
    }

    outInfo = CoordinateInfoSwap;
}

float PtrajClusteringGetStddev(PtrajClustering* This, PtrajCobwebCluster* Node, int StddevIndex)
{
    if (Node->PointIndex==-1 && Node->Stddevs[StddevIndex]>This->Acuity)
    {
        return Node->Stddevs[StddevIndex];
    }
    else
    {
        return This->Acuity; /* To avoid ever returning 0 */
    }
}


/* 1 / 2*sqrt(pi) */
#define UTILITY_MULTIPLIER 0.28209479177387814

void CobwebTreeNodeRemoveChild(CobwebTreeNode* Child); /* forward */
PtrajCobwebCluster* NewCobwebClusterLeafNode(PtrajClustering* Clustering, int FrameIndex); /* forward*/


void PrintCobwebTreeNodeName(CobwebTreeNode* Node)
{
    CobwebTreeNode* Child;
    char Name[256];
    char Temp[256];
    memset(Name,0,256);
    memset(Temp,0,256);
    if (Node->FirstChild)
    {
        fprintf(stdout,"N(");
        for (Child=Node->FirstChild; Child; Child=Child->Next)
        {
            PrintCobwebTreeNodeName(Child);
            fprintf(stdout,",");
        }
        fprintf(stdout,")");
    }
    else
    {
        /*printf(Temp,"L%d",Node->Cluster->FrameIndex);
        strcpy(Name,Temp);*/
        fprintf(stdout,"L%d",Node->Cluster->PointIndex);
    }
}

float PtrajClusteringCobwebCategoryUtility(PtrajClustering* This, CobwebTreeNode* Parent)
{
    int ClusterCount = 0;
    int DimensionIndex;
    int AtomIndex;
    int AtomCount;
    CobwebTreeNode* Child;
    float ProbabilityThisCluster;
    float StddevSum;
    float UtilitySum = 0;
    float StddevChild;
    float StddevParent;
    
    int AttributeIndex;
    int AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
    /**/
    AtomCount = PtrajClusteringGetAtomCount(This);
    /*fprintf(stdout,"CATUTIL: ");
    PrintCobwebTreeNodeName(Parent);
    fprintf(stdout,"\n");*/
    for (Child = Parent->FirstChild; Child; Child = Child->Next)
    {
        ClusterCount++;
        ProbabilityThisCluster = Child->TotalFrameCount / (float)Parent->TotalFrameCount;
        /* Iterate over all the atoms and all 3 dimensions, adding the difference-in-standard-deviations */
        StddevSum = 0;
        for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
        {
            StddevChild = PtrajClusteringGetStddev(This, (PtrajCobwebCluster *) Child->Cluster,AttributeIndex);
            StddevParent = PtrajClusteringGetStddev(This, (PtrajCobwebCluster *) Parent->Cluster,AttributeIndex);
            StddevSum += (1/StddevChild) - (1/StddevParent);
        /*
        for (DimensionIndex=0; DimensionIndex<3; DimensionIndex++)
        {
            for (AtomIndex=0; AtomIndex<AtomCount; AtomIndex++)
            {
                StddevChild = PtrajClusteringGetStddev(This,Child->Cluster,DimensionIndex*AtomCount + AtomIndex);
                StddevParent = PtrajClusteringGetStddev(This,Parent->Cluster,DimensionIndex*AtomCount + AtomIndex);
                StddevSum += (1/StddevChild) - (1/StddevParent);
                /*StddevSum += (1/PtrajClusteringGetStddev(This,Child,DimensionIndex*AtomCount + AtomIndex)) - \
                    (1/PtrajClusteringGetStddev(This,Parent,DimensionIndex*AtomCount + AtomIndex));*/
            
            /*
            }
            */
        }
        /*PrintCobwebTreeNodeName(Child);
        fprintf(stdout,"Clustersize %d odds %.4f Stddevsum %.4f net %.4f\n",Child->TotalFrameCount,ProbabilityThisCluster,StddevSum,ProbabilityThisCluster * StddevSum * UTILITY_MULTIPLIER);*/
        UtilitySum += ProbabilityThisCluster * StddevSum;
    }
    /*fprintf(stdout,"Total Catutil: %.4f\n",UtilitySum/ClusterCount);*/
    return (UtilitySum * UTILITY_MULTIPLIER / ClusterCount);
}

/*
void CobwebTreeNodeFree(CobwebTreeNode* Head)
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

void CobwebTreeNodeFreeNR(CobwebTreeNode* Head)
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
*/


/*
CobwebTreeNode* CobwebTreeNodeNew()
{
    CobwebTreeNode* This;
    This = (CobwebTreeNode*)SafeMalloc(__FILE__, __LINE__, sizeof(CobwebTreeNode));
    memset(This,0,sizeof(CobwebTreeNode));
    return This;
}

CobwebTreeNode* CobwebTreeNodeNewLeaf(PtrajClustering* This, int FrameIndex)
{
    CobwebTreeNode* CTN;
    CTN = (CobwebTreeNode*)SafeMalloc(__FILE__, __LINE__, sizeof(CobwebTreeNode));
    memset(CTN,0,sizeof(CobwebTreeNode));
    CTN->Cluster = NewCobwebClusterLeafNode(This,FrameIndex);
    CTN->TotalFrameCount = 1;
    return CTN;
}
*/


    
CobwebTreeNode* CobwebTreeCopy(CobwebTreeNode* Head)
{
    CobwebTreeNode* NewHead;
    CobwebTreeNode* Child;
    NewHead = CobwebTreeNodeNew();
    NewHead->Cluster = Head->Cluster;
    /*NewHead->TotalFrameCount = Head->TotalFrameCount;*/
    for (Child=Head->FirstChild; Child; Child=Child->Next)
    {
        CobwebTreeNodeAddChild(NewHead,CobwebTreeCopy(Child));
    }
    return NewHead;
}

/* Old PtrajClusteringNewCobwebCluster()
PtrajCobwebCluster* PtrajClusteringNewCobwebCluster(PtrajClustering* Clustering, int PointIndex)
{
    PtrajCobwebCluster* This;
    int AtomCount = PtrajClusteringGetAtomCount(Clustering);
    This = (PtrajCobwebCluster*)SafeMalloc(__FILE__, __LINE__, sizeof(PtrajCobwebCluster));
    memset(This,0,sizeof(PtrajCobwebCluster));
    if (PointIndex==-1)
    {
        This->Stddevs = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount*3);
        memset(This->Stddevs,0,sizeof(float)*AtomCount*3);
        This->Means = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount*3);
        memset(This->Means,0,sizeof(float)*AtomCount*3);
    }
    This->PointIndex=PointIndex;
    return This;
}
*/

PtrajCobwebCluster* PtrajClusteringNewCobwebCluster(PtrajClustering* pClustering, int PointIndex)
{
    PtrajCobwebCluster* This;
    int AttributeIndex;
    int AttributeCount = ClusteringGetAttributeCount( (Clustering *) pClustering);
    
    This = (PtrajCobwebCluster*)SafeMalloc(__FILE__, __LINE__, sizeof(PtrajCobwebCluster));
    memset(This,0,sizeof(PtrajCobwebCluster));
    if (PointIndex==-1)
    {
        This->Stddevs = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AttributeCount);
        memset(This->Stddevs,0,sizeof(float)*AttributeCount);
        This->Means = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AttributeCount);
        memset(This->Means,0,sizeof(float)*AttributeCount);
        This->SumX = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AttributeCount);
        memset(This->SumX,0,sizeof(float)*AttributeCount);
        This->SumX2 = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AttributeCount);
        memset(This->SumX2,0,sizeof(float)*AttributeCount);
        This->Cos = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AttributeCount);
        memset(This->Cos,0,sizeof(float)*AttributeCount);
        This->Cos2 = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AttributeCount);
        memset(This->Cos2,0,sizeof(float)*AttributeCount);
        This->Sin = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AttributeCount);
        memset(This->Sin,0,sizeof(float)*AttributeCount);
        This->Sin2 = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AttributeCount);
        memset(This->Sin2,0,sizeof(float)*AttributeCount);
    }
    This->PointIndex=PointIndex;
    return This;
}


/* Free JUST THIS NODE, not its children */
void PtrajClusteringCobwebClusterFree(PtrajClustering* Clustering, PtrajCobwebCluster* This) 
{
    if (This)
    {
        safe_free(This->Stddevs);
        safe_free(This->Means);
        safe_free(This->SumX);
        safe_free(This->SumX2);
        safe_free(This->Sin);
        safe_free(This->Sin2);
        safe_free(This->Cos);
        safe_free(This->Cos2);
    }
    safe_free(This);
}

/* Old  PtrajClusteringCobwebCopyMeans()
void PtrajClusteringCobwebCopyMeans(PtrajClustering* This, PtrajCobwebCluster* Source, PtrajCobwebCluster* Destination)
{
    int AtomCount = PtrajClusteringGetAtomCount(This);
    memcpy(Destination->Means,Source->Means,sizeof(float)*AtomCount*3);
    memcpy(Destination->Stddevs,Source->Stddevs,sizeof(float)*AtomCount*3);
}
*/
void PtrajClusteringCobwebCopyMeans(PtrajClustering* This, PtrajCobwebCluster* Source, PtrajCobwebCluster* Destination)
{
    int AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
    memcpy(Destination->Means,Source->Means,sizeof(float)*AttributeCount);
    memcpy(Destination->Stddevs,Source->Stddevs,sizeof(float)*AttributeCount);
    memcpy(Destination->SumX,Source->SumX,sizeof(float)*AttributeCount);
    memcpy(Destination->SumX2,Source->SumX2,sizeof(float)*AttributeCount);
    memcpy(Destination->Cos,Source->Cos,sizeof(float)*AttributeCount);
    memcpy(Destination->Cos2,Source->Cos2,sizeof(float)*AttributeCount);
    memcpy(Destination->Sin,Source->Sin,sizeof(float)*AttributeCount);
    memcpy(Destination->Sin2,Source->Sin2,sizeof(float)*AttributeCount);
}



/* Accumulate standard-deviation in the node TopParent, for the subtree headed by CurrentNode. */
void CobwebClusterStddevRecurse(CobwebTreeNode* TopParent, CobwebTreeNode* CurrentNode, int AtomCount, trajectoryInfo* trajInfo, PtrajClustering* This)
{
    int AtomIndex;
    float Distance;
    int AttributeIndex;
    int AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
    CobwebTreeNode* Child;
    PtrajCobwebCluster* Cluster = (PtrajCobwebCluster *) CurrentNode->Cluster;
    PtrajCobwebCluster* TopParentCluster = (PtrajCobwebCluster *) TopParent->Cluster;
    /* If this is a leaf, add on the standard deviations and finish */
    if (Cluster->PointIndex!=-1)
    {
        for (AttributeIndex=0; AttributeIndex<AttributeCount; AttributeIndex++)
        {
        	if (!PtrajClusteringIsAttributeTorsion(This, AttributeIndex)) continue;
            Distance = 180-abs(180-abs(ClusteringGetAttributeValue( (Clustering *) This,Cluster->PointIndex,AttributeIndex) - TopParentCluster->Means[AttributeIndex]));
            TopParentCluster->Stddevs[AttributeIndex] += Distance*Distance / (float)TopParent->TotalFrameCount;
        	
        }
        /*
        for (AtomIndex=0; AtomIndex<AtomCount; AtomIndex++)
        {
            Distance = trajInfo->x[Cluster->PointIndex*AtomCount + AtomIndex] - TopParentCluster->Means[AtomIndex];
            TopParentCluster->Stddevs[AtomIndex] += Distance*Distance / (float)TopParent->TotalFrameCount;
            Distance = trajInfo->y[Cluster->PointIndex*AtomCount + AtomIndex] - TopParentCluster->Means[AtomCount+AtomIndex];
            TopParentCluster->Stddevs[AtomCount+AtomIndex] += Distance*Distance / (float)TopParent->TotalFrameCount;
            Distance = trajInfo->z[Cluster->PointIndex*AtomCount + AtomIndex] - TopParentCluster->Means[2*AtomCount+AtomIndex];
            TopParentCluster->Stddevs[2*AtomCount+AtomIndex] += Distance*Distance / (float)TopParent->TotalFrameCount;
        }
        */
        return;
    }

    /* If this wasn't a leaf: Handle all of the children */
    for (Child = CurrentNode->FirstChild; Child; Child=Child->Next)
    {
        CobwebClusterStddevRecurse(TopParent,Child, AtomCount, trajInfo, This);
    }
}


/* 
Compute the MEANS for this PtrajCobwebCluster node, as well as the STANDARD DEVIATIONS.  Computing the means for our three
is pretty easy: A child leaf's coordinates are simply added to the running total.  A child node with mean X and 50 nodes 
can be treated as a lump (we add 50*X to the running total).  Computing the standard deviations is a bit more difficult,
because we must traverse the tree to do so.
*/
void PtrajCobwebClusterComputeMeans(PtrajClustering* This, PtrajCobwebTreeNode* Node)
{
    CobwebTreeNode* Child;
    PtrajCobwebCluster* Cluster = Node->Cluster;
    int AtomIndex;
    int DummyIndex;
    int AtomCount = PtrajClusteringGetAtomCount(This);
    trajectoryInfo* trajInfo = This->trajInfo;
    float Distance;
    int AttributeIndex;
    int AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
    float value;
    float Sin, Cos;
    
    /**/
    if (Node->Cluster->PointIndex != -1)
    {
        return;
    }
    /* Compute the MEANS! */
    /*memset(Cluster->Means,0,sizeof(float)*AtomCount*3);*/
    memset(Cluster->Means,0,sizeof(float)*AttributeCount);
    memset(Cluster->Stddevs,0,sizeof(float)*AttributeCount);
    memset(Cluster->SumX,0,sizeof(float)*AttributeCount);
    memset(Cluster->SumX2,0,sizeof(float)*AttributeCount);
    
    memset(Cluster->Cos,0,sizeof(float)*AttributeCount);
    memset(Cluster->Cos2,0,sizeof(float)*AttributeCount);
    memset(Cluster->Sin,0,sizeof(float)*AttributeCount);
    memset(Cluster->Sin2,0,sizeof(float)*AttributeCount);
    
    for (Child = Node->FirstChild; Child; Child=Child->Next)
    {
        if (Child->Cluster->PointIndex==-1)
        {
            /*
            for (DummyIndex=0; DummyIndex<AtomCount*3; DummyIndex++)
            {
                Cluster->Means[DummyIndex] += ((PtrajCobwebCluster*)Child->Cluster)->Means[DummyIndex]*Child->TotalFrameCount/Node->TotalFrameCount;
            }*/
            for (AttributeIndex=0; AttributeIndex<AttributeCount; AttributeIndex++)
            {
                if (PtrajClusteringIsAttributeTorsion(This, AttributeIndex))
                {
	                Cluster->Sin[AttributeIndex] += ((PtrajCobwebCluster*)Child->Cluster)->Sin[AttributeIndex];
    	            Cluster->Sin2[AttributeIndex] += ((PtrajCobwebCluster*)Child->Cluster)->Sin2[AttributeIndex];
	                Cluster->Cos[AttributeIndex] += ((PtrajCobwebCluster*)Child->Cluster)->Cos[AttributeIndex];
    	            Cluster->Cos2[AttributeIndex] += ((PtrajCobwebCluster*)Child->Cluster)->Cos2[AttributeIndex];
                } else {
	                Cluster->SumX[AttributeIndex] += ((PtrajCobwebCluster*)Child->Cluster)->SumX[AttributeIndex];
    	            Cluster->SumX2[AttributeIndex] += ((PtrajCobwebCluster*)Child->Cluster)->SumX2[AttributeIndex];
                }
                /*
                Cluster->Means[AttributeIndex] += ((PtrajCobwebCluster*)Child->Cluster)->Means[AttributeIndex]*Child->TotalFrameCount/Node->TotalFrameCount;
                */
            }    
        }
        else
        {
            /*
            for (AtomIndex=0; AtomIndex<AtomCount; AtomIndex++)
            {
                Cluster->Means[AtomIndex]             += trajInfo->x[Child->Cluster->PointIndex*AtomCount + AtomIndex] / (float)Node->TotalFrameCount;
                Cluster->Means[AtomIndex+AtomCount]   += trajInfo->y[Child->Cluster->PointIndex*AtomCount + AtomIndex] / (float)Node->TotalFrameCount;
                Cluster->Means[AtomIndex+2*AtomCount] += trajInfo->z[Child->Cluster->PointIndex*AtomCount + AtomIndex] / (float)Node->TotalFrameCount;
            }*/
            for (AttributeIndex=0; AttributeIndex<AttributeCount; AttributeIndex++)
            {
                if (PtrajClusteringIsAttributeTorsion(This, AttributeIndex))
                {
	                value = PtrajClusteringGetAttributeArraySin(This,Child->Cluster->PointIndex,AttributeIndex);
	                Cluster->Sin[AttributeIndex] += value;
    	            Cluster->Sin2[AttributeIndex] += value*value;
	                value = PtrajClusteringGetAttributeArrayCos(This,Child->Cluster->PointIndex,AttributeIndex);
	                Cluster->Cos[AttributeIndex] += value;
    	            Cluster->Cos2[AttributeIndex] += value*value;
                } else {
	                value = ClusteringGetAttributeValue( (Clustering *) This,
							     Child->Cluster->PointIndex,AttributeIndex);
	                Cluster->SumX[AttributeIndex] += value;
    	            Cluster->SumX2[AttributeIndex] += value*value;
                }
                /*
                Cluster->Means[AttributeIndex] += ClusteringGetAttributeValue(This,Child->Cluster->PointIndex,AttributeIndex)/Node->TotalFrameCount;*/
            }    
        }
    }
    /* Now that we have the means, compute the STANDARD DEVIATIONS! */
    /*memset(Cluster->Stddevs,0,sizeof(float)*AtomCount*3);*/
    memset(Cluster->Stddevs,0,sizeof(float)*AttributeCount);
    for (AttributeIndex=0; AttributeIndex<AttributeCount; AttributeIndex++)
    {
        if (PtrajClusteringIsAttributeTorsion(This, AttributeIndex))
        {
    	    Sin = Cluster->Sin[AttributeIndex];
    	    Cos = Cluster->Cos[AttributeIndex];
            Cluster->Means[AttributeIndex] = atan2(Sin, Cos) * RADDEG;
    	    /*
            Cluster->Stddevs[AttributeIndex] = (Cluster->Sin2[AttributeIndex] + Cluster->Cos2[AttributeIndex])/Node->TotalFrameCount - (Sin*Sin + Cos*Cos)/(Node->TotalFrameCount*Node->TotalFrameCount);
            Cluster->Stddevs[AttributeIndex] = acos(1-Cluster->Stddevs[AttributeIndex]/2) * RADDEG;
            */
            
            /*
            value = sqrt(Sin*Sin+Cos*Cos)/Node->TotalFrameCount;
            if (value < 0.001) 
            	Cluster->Stddevs[AttributeIndex] = - log(0.001);
            else
            	Cluster->Stddevs[AttributeIndex] = - log(value);
            */
            
        
        } else {
    	    Cluster->Means[AttributeIndex] = Cluster->SumX[AttributeIndex]/Node->TotalFrameCount;
	        Cluster->Stddevs[AttributeIndex] = sqrt(Cluster->SumX2[AttributeIndex]/Node->TotalFrameCount - Cluster->Means[AttributeIndex]*Cluster->Means[AttributeIndex]);
        }
    }
    
    for (Child = Node->FirstChild; Child; Child=Child->Next)
    {
        CobwebClusterStddevRecurse((CobwebTreeNode *) Node,Child,AtomCount,trajInfo,This);
    }
    for (AttributeIndex=0; AttributeIndex<AttributeCount; AttributeIndex++)
    {
       	if (!PtrajClusteringIsAttributeTorsion(This, AttributeIndex)) continue;
        Cluster->Stddevs[AttributeIndex] = sqrt(Cluster->Stddevs[AttributeIndex]);
    }
    
    
    /* Take the square roots (to go from variance to std-dev) */
    /*
    for (DummyIndex=0; DummyIndex<AtomCount*3; DummyIndex++)
    {
        Cluster->Stddevs[DummyIndex] = sqrt(Cluster->Stddevs[DummyIndex]);
    }
    */
    /*
    for (AttributeIndex=0; AttributeIndex<AttributeCount; AttributeIndex++)
    {
        Cluster->Stddevs[AttributeIndex] = sqrt(Cluster->Stddevs[AttributeIndex]);
    }
    */
}




void CobwebTreeNodeStripChildren(CobwebTreeNode* Parent)
{
    Parent->TotalFrameCount = 0;
    Parent->FirstChild = NULL;
    Parent->LastChild = NULL;
}

void PtrajClusteringCobwebAddNewFrame(PtrajClustering* This, CobwebTreeNode* Root, CobwebTreeNode* FrameNode, int level)
{
	float BestHostCU, CU, RootCU, NextBestHostCU;
    float BestSplitCU, BestMergeCU, BestCU;
    int ChildCount;
    CobwebTreeNode *Child, *PrevChild, *Parent;
    CobwebTreeNode *BestHost, *NextBestHost;
    CobwebTreeNode *BestSplitNode, *MergeNode1, *MergeNode2;
    CobwebTreeNode *NewNode;
    CTNList *ChildList, *MergeList, *Item1, *Item2;
    CobwebCluster *SwapRootMeans;
    int verbose = This->action->darg2;
    
    /* If there is 0 or 1 child, simply add a child. */
    if (Root->FirstChild == Root->LastChild || !Root->FirstChild) {
    	if (Root->Cluster->PointIndex == COBWEB_NONLEAF) {
        	CobwebTreeNodeAddChild(Root, FrameNode);
        	CobwebClusterComputeMeans( (Clustering *) This, Root);
        } else {
        	NewNode = CobwebTreeNodeNew();
            NewNode->Cluster = ClusteringNewCobwebCluster( (Clustering *) This, COBWEB_NONLEAF);
            Parent = Root->Parent;
            CobwebTreeNodeRemoveChild(Root);
            CobwebTreeNodeAddChild(NewNode, Root);
            CobwebTreeNodeAddChild(NewNode, FrameNode);
            CobwebTreeNodeAddChild(Parent, NewNode);
        }
        if (Root->Parent) {
        	for (Parent = Root->Parent; Parent; Parent = Parent->Parent) {
            	CobwebClusterComputeMeans( (Clustering *) This, Parent);
            }
        }
        return;
    }
    
    ChildList = CobwebTreeNodeListChildren( (CobwebTreeNode *) Root);
    MergeList = NULL;
    BestHost = NextBestHost = NULL;
    BestHostCU = NextBestHostCU = -1000;
    RootCU = TryCategoryUtility( (Clustering *) This, Root, Root, FrameNode);
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\nTry CU:\n");
    for (Item1 = ChildList; Item1; Item1 = Item1->Next) {
    	Child = Item1->Node;
    	CU = TryCategoryUtility( (Clustering *) This, Root, Child, FrameNode);
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "%10.4f    ", CU);
        if (CU > BestHostCU) {
        	NextBestHost = BestHost;
            NextBestHostCU = BestHostCU;
            BestHost = Child;
            BestHostCU = CU;
            MergeList = CTNListAppend(MergeList, Child);
        } else {
        	if (CU > NextBestHostCU) {
        		NextBestHost = Child;
            	NextBestHostCU = CU;
	        }
            if (Child->TotalFrameCount == 1) {
    	    	MergeList = CTNListAppend(MergeList, Child);
        	}
        }
    }
	if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\n");
    MergeList = CTNListAppend(MergeList, NextBestHost);
	/* The MergeList will store a list of nodes which may have high CU if NewFrame is added to. 
       This list consists of the NextBestNode and all node with only one frame, as well as all the nodes have been selected as BestNode temporarily. 
       I will use this list to find the BestMerge nodes instead of just the BestNode and NextBestNode as in classic COBWEB. */
    
    
    /* print MergeList. */
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) {
    	fprintf(stdout, "\nMergeList:");
    	int i;
    	for (Item2 = MergeList; Item2->Previous; Item2 = Item2->Previous);
    	for (Item1 = Item2; Item1; Item1 = Item1->Next) {
			Child = Item1->Node;
		    for (Item2 = ChildList, i = 0; Item2; Item2 = Item2->Next, i++) {
    	    	if (Item2->Node == Item1->Node) fprintf(stdout, "%d, ", i);
    	    }
	    }
    	fprintf(stdout, "\n");
    }
    /* Temporally add FrameNode as Root's child. So we can try Split and Merge. */
    CobwebTreeNodeAddChild(Root, FrameNode);
    CobwebClusterComputeMeans( (Clustering *) This, Root);
    MergeNode1 = BestHost;
    MergeNode2 = NextBestHost;
    BestMergeCU = CobwebTryMerge( (Clustering *) This, Root, BestHost, NextBestHost, 0);
    BestSplitNode = BestHost;
    BestSplitCU = CobwebTrySplit( (Clustering *) This, Root, BestHost, 0);
    /* In classic COBWEB method, only BestNode is subject to splitting.
       Here, we will try to split all nodes, and find the BestSplitNode.
       As a matter of fact, the BestSplitNode is not always BestNode, sometimes is NextBestNode, but sometimes is not even in the MergeList. */
    /*
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\nTry splitting:\n");
    for (Item1 = ChildList; Item1; Item1 = Item1->Next) {
    	Child = Item1->Node;
    	if (Child == FrameNode) continue;
        CU = CobwebTrySplit(This, Root, Child, 0);
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "%10.4f    ", CU);
        if(CU > BestSplitCU)  {
        	BestSplitCU = CU;
            BestSplitNode = Child;
        }
    }*/
    /* In classic COBWEB method, only the BestNode and the NextBestNode will be tried for merging.
       Here, we try those nodes in the MergeList, which have higher probability to find higher CU.
       As a matter of fact, the BestMergeNode is not necessarily the merge of the BestNode and the NextBestNode. */
    /*
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\nTry merging:\n");
    BestMergeCU =  -1000;
    for (Item1 = MergeList; Item1; Item1 = Item1->Previous) {
    	for (Item2 = Item1->Previous; Item2; Item2 = Item2->Previous) {
        	if (Item1->Node == Item2->Node) continue;
            CU = CobwebTryMerge(This, Root, Item1->Node, Item2->Node, 0);
		    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "%10.4f    ", CU);
            if (CU > BestMergeCU) {
            	BestMergeCU = CU;
                MergeNode1 = Item1->Node;
                MergeNode2 = Item2->Node;
            }
        }
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\n");
    }
    */
    CobwebTreeNodeRemoveChild(FrameNode);
    CobwebClusterComputeMeans( (Clustering *) This, Root);
	
    /* MergeList point to the last Node of the list, and CTNListFree need the Head of the list. */
    for (Item1 = MergeList; Item1->Previous; Item1 = Item1->Previous);
    CTNListFree(Item1);
    CTNListFree(ChildList);

    
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER)
	    fprintf(stdout, "RootCU is %.4f, BestHostCU is %.4f, BestMergeCU is %.4f, BestSplitCU is %.4f\n", RootCU, BestHostCU, BestMergeCU, BestSplitCU);
    BestCU = max(max(RootCU, BestMergeCU), BestSplitCU);
    if (BestCU == BestMergeCU && BestMergeCU > 0) {
    	/* Use CobwebMerge instead of CobwebTryMerge, because CobwebMerge will try to split the MergeNode if appropriate. */
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "Merging..\n");
    	CobwebMerge( (Clustering *) This, Root, MergeNode1, MergeNode2);
        PtrajClusteringCobwebAddNewFrame(This, Root, FrameNode, level);
    } else if (BestCU == BestSplitCU && BestSplitCU > 0) {
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "Splitting..\n");
    	CobwebTrySplit( (Clustering *) This, Root, BestSplitNode, 1);
        PtrajClusteringCobwebAddNewFrame(This, Root, FrameNode, level);
    } else if (RootCU > BestHostCU && RootCU > 0) {
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "At root ..\n");
    	CobwebTreeNodeAddChild(Root, FrameNode);
    } else {
#define MAX_LEVEL 50
	    fprintf(stdout, "Goto host ..\n");
    	if (level < MAX_LEVEL) 
    		PtrajClusteringCobwebAddNewFrame(This, BestHost, FrameNode, level+1);
        else
	    	CobwebTreeNodeAddChild(Root, FrameNode);
    }
    CobwebClusterComputeMeans( (Clustering *) This, Root);
}

void PtrajClusteringCobwebAddNewFrame_MergeList(PtrajClustering* This, CobwebTreeNode* Root, CobwebTreeNode* FrameNode, int level)
{
	float BestHostCU, CU, RootCU, NextBestHostCU;
    float BestSplitCU, BestMergeCU, BestCU;
    int ChildCount;
    CobwebTreeNode *Child, *PrevChild, *Parent;
    CobwebTreeNode *BestHost, *NextBestHost;
    CobwebTreeNode *BestSplitNode, *MergeNode1, *MergeNode2;
    CobwebTreeNode *NewNode;
    CTNList *ChildList, *MergeList, *Item1, *Item2;
    CobwebCluster *SwapRootMeans;
    int verbose = This->action->darg2;
    
    /* If there is 0 or 1 child, simply add a child. */
    if (Root->FirstChild == Root->LastChild || !Root->FirstChild) {
    	if (Root->Cluster->PointIndex == COBWEB_NONLEAF) {
        	CobwebTreeNodeAddChild(Root, FrameNode);
        	CobwebClusterComputeMeans( (Clustering *) This, Root);
        } else {
        	NewNode = CobwebTreeNodeNew();
            NewNode->Cluster = ClusteringNewCobwebCluster( (Clustering *) This, COBWEB_NONLEAF);
            Parent = Root->Parent;
            CobwebTreeNodeRemoveChild(Root);
            CobwebTreeNodeAddChild(NewNode, Root);
            CobwebTreeNodeAddChild(NewNode, FrameNode);
            CobwebTreeNodeAddChild(Parent, NewNode);
        }
        if (Root->Parent) {
        	for (Parent = Root->Parent; Parent; Parent = Parent->Parent) {
            	CobwebClusterComputeMeans( (Clustering *) This, Parent);
            }
        }
        return;
    }
    
    ChildList = CobwebTreeNodeListChildren(Root);
    MergeList = NULL;
    BestHost = NextBestHost = NULL;
    BestHostCU = NextBestHostCU = -1000;
    RootCU = TryCategoryUtility( (Clustering *) This, Root, Root, FrameNode);
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\nTry CU:\n");
    for (Item1 = ChildList; Item1; Item1 = Item1->Next) {
    	Child = Item1->Node;
    	CU = TryCategoryUtility( (Clustering *) This, Root, Child, FrameNode);
	    fprintf(stdout, "%10.4f    ", CU);
        if (CU > BestHostCU) {
        	NextBestHost = BestHost;
            NextBestHostCU = BestHostCU;
            BestHost = Child;
            BestHostCU = CU;
            MergeList = CTNListAppend(MergeList, Child);
        } else {
        	if (CU > NextBestHostCU) {
        		NextBestHost = Child;
            	NextBestHostCU = CU;
	        }
            if (Child->TotalFrameCount == 1) {
    	    	MergeList = CTNListAppend(MergeList, Child);
        	}
        }
    }
    MergeList = CTNListAppend(MergeList, NextBestHost);
	/* The MergeList will store a list of nodes which may have high CU if NewFrame is added to. 
       This list consists of the NextBestNode and all node with only one frame, as well as all the nodes have been selected as BestNode temporarily. 
       I will use this list to find the BestMerge nodes instead of just the BestNode and NextBestNode as in classic COBWEB. */
    
    
    /* print MergeList. */
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) {
    	fprintf(stdout, "\nMergeList:");
    	int i;
    	for (Item2 = MergeList; Item2->Previous; Item2 = Item2->Previous);
    	for (Item1 = Item2; Item1; Item1 = Item1->Next) {
			Child = Item1->Node;
		    for (Item2 = ChildList, i = 0; Item2; Item2 = Item2->Next, i++) {
    	    	if (Item2->Node == Item1->Node) fprintf(stdout, "%d, ", i);
    	    }
	    }
    	fprintf(stdout, "\n");
    }
    /* Temporally add FrameNode as Root's child. So we can try Split and Merge. */
    CobwebTreeNodeAddChild(Root, FrameNode);
    CobwebClusterComputeMeans( (Clustering *) This, Root);
    BestSplitCU = -1000;
    BestSplitNode = NULL;
    
    /* In classic COBWEB method, only BestNode is subject to splitting.
       Here, we will try to split all nodes, and find the BestSplitNode.
       As a matter of fact, the BestSplitNode is not always BestNode, sometimes is NextBestNode, but sometimes is not even in the MergeList. */
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\nTry splitting:\n");
    for (Item1 = ChildList; Item1; Item1 = Item1->Next) {
    	Child = Item1->Node;
    	if (Child == FrameNode) continue;
        CU = CobwebTrySplit( (Clustering *) This, Root, Child, 0);
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "%10.4f    ", CU);
        if(CU > BestSplitCU)  {
        	BestSplitCU = CU;
            BestSplitNode = Child;
        }
    }
    /* In classic COBWEB method, only the BestNode and the NextBestNode will be tried for merging.
       Here, we try those nodes in the MergeList, which have higher probability to find higher CU.
       As a matter of fact, the BestMergeNode is not necessarily the merge of the BestNode and the NextBestNode. */
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\nTry merging:\n");
    BestMergeCU =  -1000;
    for (Item1 = MergeList; Item1; Item1 = Item1->Previous) {
    	for (Item2 = Item1->Previous; Item2; Item2 = Item2->Previous) {
        	if (Item1->Node == Item2->Node) continue;
            CU = CobwebTryMerge( (Clustering *) This, Root, Item1->Node, Item2->Node, 0);
		    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "%10.4f    ", CU);
            if (CU > BestMergeCU) {
            	BestMergeCU = CU;
                MergeNode1 = Item1->Node;
                MergeNode2 = Item2->Node;
            }
        }
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "\n");
    }
    CobwebTreeNodeRemoveChild(FrameNode);
    CobwebClusterComputeMeans( (Clustering *) This, Root);
	
    /* MergeList point to the last Node of the list, and CTNListFree need the Head of the list. */
    for (Item1 = MergeList; Item1->Previous; Item1 = Item1->Previous);
    CTNListFree(Item1);
    CTNListFree(ChildList);

    
    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER)
	    fprintf(stdout, "RootCU is %.4f, BestHostCU is %.4f, BestMergeCU is %.4f, BestSplitCU is %.4f\n", RootCU, BestHostCU, BestMergeCU, BestSplitCU);
    BestCU = max(max(RootCU, BestMergeCU), BestSplitCU);
    if (BestCU == BestMergeCU && BestMergeCU > 0) {
    	/* Use CobwebMerge instead of CobwebTryMerge, because CobwebMerge will try to split the MergeNode if appropriate. */
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "Merging..\n");
    	CobwebMerge( (Clustering *) This, Root, MergeNode1, MergeNode2);
        PtrajClusteringCobwebAddNewFrame(This, Root, FrameNode, level);
    } else if (BestCU == BestSplitCU && BestSplitCU > 0) {
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "Splitting..\n");
    	CobwebTrySplit( (Clustering *) This, Root, BestSplitNode, 1);
        PtrajClusteringCobwebAddNewFrame(This, Root, FrameNode, level);
    } else if (RootCU > BestHostCU && RootCU > 0) {
	    if (verbose > VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "At root ..\n");
    	CobwebTreeNodeAddChild(Root, FrameNode);
    } else {
#define MAX_LEVEL 50
	    fprintf(stdout, "Goto host ..\n");
    	if (level < MAX_LEVEL) 
    		PtrajClusteringCobwebAddNewFrame(This, BestHost, FrameNode, level+1);
        else
	    	CobwebTreeNodeAddChild(Root, FrameNode);
    }
    CobwebClusterComputeMeans( (Clustering *) This, Root);
}


/*
This is the workhorse function of COBWEB clustering.  It handles simple insertions, merges, splits, and 
the selection of hosts.  This function RECURSES often: Assume the root has a leaf-child A and child B with its
own sub-tree.  If B is the best parent, we call PtrajClusteringCobwebAddNewFrame again to insert our frame into B.
*/
void PtrajClusteringCobwebAddNewFrame_old(PtrajClustering* This,CobwebTreeNode* Root,CobwebTreeNode* FrameNode)
{
    int AtomCount = PtrajClusteringGetAtomCount(This);
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
    PtrajCobwebCluster* SwapMeans = PtrajClusteringNewCobwebCluster(This,-1);
    /**/
    /* If there are 0 or 1 children, simply add a child! */
    if (!Root->FirstChild || Root->FirstChild==Root->LastChild)
    {
        CobwebTreeNodeAddChild(Root,FrameNode);
        CobwebClusterComputeMeans( (Clustering *) This, Root);
        if (Root->Parent)
        {
            for (Parent = Root->Parent; Parent; Parent=Parent->Parent)
            {
                CobwebClusterComputeMeans( (Clustering *) This,Parent);
            }
        }
        PtrajClusteringCobwebClusterFree(This,SwapMeans);
        return;
    }
    /* Consider adding a child to the root, as that may be best: */
    BestHost = Root;
    BestHostCU = TryCategoryUtility( (Clustering *) This,Root, Root,FrameNode);
    /* 
    Ok: We have at least 2 children.  Consider the CUs of adding the new node to each of our children!  Keep track
    of the best host, and the runner-up. 
    */
    for (Child=Root->FirstChild; Child; Child=Child->Next)
    {
        CU = TryCategoryUtility( (Clustering *) This, Root, Child, FrameNode); /* Manufactures a new host which MAY later replace our child */
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
        StatusQuoCU = PtrajClusteringCobwebCategoryUtility(This,Root);       
        if (NextBestHost!=Root && RootChildCount>2)
        {
            /***************************/
            /* Consider merging the two best hosts, if that would improve CU */
            /*TempListA = CopyChildren(BestHost);
            TempListB = CopyChildren(NextBestHost);*/
            CobwebTreeNodeRemoveChild(BestHost);
            CobwebTreeNodeRemoveChild(NextBestHost);
            MergedMaster = CobwebClusterMergeNodes( (Clustering *) This,BestHost, NextBestHost);
            /* Temporarily insert this node, in place of its sources */
            CobwebTreeNodeAddChild(Root,MergedMaster);
            CobwebClusterComputeMeans( (Clustering *) This,MergedMaster);
            PtrajClusteringCobwebCopyMeans( This, (PtrajCobwebCluster *) Root->Cluster,SwapMeans);
            /*OldMeans = Root->Cluster->Means;
            OldStddevs = Root->Cluster->Stddevs;
            Root->Cluster->Means = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount*3);
            Root->Cluster->Stddevs = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount*3);*/
            CobwebClusterComputeMeans( (Clustering *) This,Root);
            CU = PtrajClusteringCobwebCategoryUtility(This,Root);
            if (CU <= StatusQuoCU)
            {
                /* No - merging will not improve our CU.  Undo the mischief we've done to our tree: */
                CobwebTreeNodeRemoveChild(MergedMaster);
                PtrajClusteringCobwebClusterFree(This, (PtrajCobwebCluster *) MergedMaster->Cluster);
                CobwebTreeNodeFreeNR(MergedMaster);
                PtrajClusteringCobwebCopyMeans(This,SwapMeans, (PtrajCobwebCluster *) Root->Cluster);
                CobwebTreeNodeAddChild(Root,BestHost);
                CobwebTreeNodeAddChild(Root,NextBestHost);
            }
            else
            {
                /* The merge was a GOOD idea!  Finalize it, and add our frame as a child of the merged node */
                /*fprintf(stdout,"---We did a merge.  Here's where we stand:\n");
		if (prnlev > 5) DebugPrintCobwebTree(This,Root,0);*/
                PtrajClusteringCobwebAddNewFrame_old(This,MergedMaster,FrameNode);
                PtrajClusteringCobwebClusterFree(This,SwapMeans);
                return;
            }
        }
        /*************************************************/
        /* Consider splitting the best host!  It may be a good idea to replace it in the tree with all-its-children */
        if (BestHost->Cluster->PointIndex==-1)
        {
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

            /*OldMeans = Root->Cluster->Means;
            OldStddevs = Root->Cluster->Stddevs;
            Root->Cluster->Means = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount*3);
            Root->Cluster->Stddevs = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float)*AtomCount*3);*/
            PtrajClusteringCobwebCopyMeans(This, (PtrajCobwebCluster *) Root->Cluster,SwapMeans);
            CobwebClusterComputeMeans( (Clustering *) This,Root);
            CU = PtrajClusteringCobwebCategoryUtility(This,Root);
            if (CU <= StatusQuoCU)
            {
                /* No - splitting will not improve our CU.  Undo the mischief we've done to our tree: */
                /*CobwebTreeNodeStripChildren(BestHost);*/
                for (TempListItem=TempList; TempListItem; TempListItem=TempListItem->Next)
                {
                    CobwebTreeNodeRemoveChild(TempListItem->Node);
                    CobwebTreeNodeAddChild(BestHost,TempListItem->Node);
                }
                CTNListFree(TempList);
                /*safe_free(Root->Cluster->Means);
                Root->Cluster->Means = OldMeans;
                safe_free(Root->Cluster->Stddevs);
                Root->Cluster->Stddevs = OldStddevs;*/
                PtrajClusteringCobwebCopyMeans(This,SwapMeans, (PtrajCobwebCluster *) Root->Cluster);
                CobwebTreeNodeAddChild(Root,BestHost);
                /*fprintf(stdout,"---We considered splitting, but didn't:\n");
                if (prnlev > 5) DebugPrintCobwebTree(Root,0);*/
            }
            else
            {
                /* Looks like that zany splitting scheme of ours was a good plan, after all!  Finalize it, and begin again: */
                /*safe_free(OldMeans);
                safe_free(OldStddevs);*/
                PtrajClusteringCobwebClusterFree(This, (PtrajCobwebCluster *) BestHost->Cluster);
                /*CobwebClusterFree(BestHost->Cluster);*/
                CobwebTreeNodeFreeNR(BestHost);
                /*fprintf(stdout,"---We split!  Here's where we stand:\n");
                if (prnlev > 5) DebugPrintCobwebTree(This,Root,0);*/
                PtrajClusteringCobwebAddNewFrame_old(This,Root,FrameNode);
                PtrajClusteringCobwebClusterFree(This,SwapMeans);
                return;
            }
        }
    }
    /**************************************************/
    /* Add the frame node as a child of the best host */
    if (BestHost==Root)
    {
        CobwebTreeNodeAddChild(Root,FrameNode);
        CobwebClusterComputeMeans( (Clustering *) This,Root);
        if (Root->Parent)
        {
            for (Parent = Root->Parent; Parent; Parent=Parent->Parent)
            {
                /*Parent->TotalFrameCount += 1;*/
                CobwebClusterComputeMeans( (Clustering *) This,Parent);
            }
        }

    }
    else
    {
        CobwebClusterHostNode(This,BestHost,FrameNode);
    }
    /*fprintf(stdout,"---At the end of Incremental, here's what we've got:\n");
    if (prnlev > 5) DebugPrintCobwebTree(Root,0);*/

}

void PtrajCobwebFlattenTree(Clustering* This, CobwebTreeNode* Root, int DesiredClusterCount)
{
  CTNList* CoalescentTail;
  int ClusterIndex;
  CTNList* CNode;
  int Backtracking;
  CobwebTreeNode* TreeNode;
  Cluster* TempCluster;
  /**/
 
  CobwebTreeReorganize(This, Root, DesiredClusterCount);
  CoalescentTail = CobwebTreeNodeListChildren(Root);
  for (CNode = CoalescentTail; CNode; CNode = CNode->Next) 
  	  CoalescentTail = CNode;
  
  DebugWriteCobwebTreeToFile(This, "CobwebCoalesce.txt", Root);
  ClusterIndex = 0;
  for (CNode=CoalescentTail; CNode; CNode = CNode->Previous)
  {
      TempCluster = ClusteringGetNewCluster(This);
      ClusteringAddCluster(This,TempCluster);
      TreeNode = CNode->Node;
      Backtracking=0;        
      while (1) {
            /* After hitting a node with no child or next, we backtrack up the 
               tree to the next time we can move to a next-node. */
          if (Backtracking == 1) {
              if (TreeNode==CNode->Node) break;
              TreeNode = TreeNode->Parent;
              if (!TreeNode || TreeNode==CNode->Node) break;
              if (TreeNode->Next) {
                  TreeNode = TreeNode->Next;
                  Backtracking = 0; 
              } else {
                  continue;
              }
          }
          if (TreeNode->Cluster->PointIndex != COBWEB_NONLEAF) {
              PtrajClusterAddMember( (PtrajCluster *) TempCluster,TreeNode->Cluster->PointIndex);
	          ClusterFindCentroid(TempCluster);
          }
          if (TreeNode->FirstChild) {
              TreeNode = TreeNode->FirstChild;
              continue;
          }
          if (TreeNode->Next && TreeNode!=CNode->Node) {
              TreeNode = TreeNode->Next;
              continue;
          }
          Backtracking=1;
      }
      ClusterIndex++;
  }
}

/* mode == 1, sequentially modify each points.
   mode == 2, randomly pick point for modification.
*/
void PtrajClusteringClusterCobweb(Clustering* This, int DesiredClusterCount, int mode)
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
    int* FinishedPoints = NULL;
    int ProcessIndex;

    /**/
    Root = ReadCobwebTree(This, "CobwebPreCoalesce.txt");
    if (Root) {
        if (prnlev > 5) DebugPrintCobwebTree(This,Root,0, -1);
    	PtrajCobwebFlattenTree(This, Root, DesiredClusterCount);
	
        if (prnlev > 5) fprintf(stdout,"\n\nFinal Cobweb tree:\n");
        if (prnlev > 5) DebugPrintCobwebTree(This,Root,0, -1);
    	return;
    } else {
	    Root = CobwebTreeNodeNew();
    	Root->Cluster = ClusteringNewCobwebCluster(This,COBWEB_NONLEAF);
	    
	    FinishedPoints = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
	    memset(FinishedPoints, 0, sizeof(int)*PointCount);
	    
	    for (ProcessIndex=0; ProcessIndex<PointCount; ProcessIndex++) {
    		if (mode == 1) PointIndex = ProcessIndex;
    		else PointIndex = ChooseNextPoint(FinishedPoints, PointCount,PointCount-ProcessIndex);
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
	}
    /* We have our tree; we'll coalesce it into managable clusters */
    PtrajCobwebFlattenTree(This, Root, DesiredClusterCount);
    if (prnlev > 5) {
      fprintf(stdout,"\n\nFinal Cobweb tree:\n");
      DebugPrintCobwebTree(This,Root,0, -1);
    }

    safe_free(FinishedPoints);
}





typedef struct CoalescentNode
{
    struct CobwebTreeNode* Node;
    struct CoalescentNode* Previous;
    struct CoalescentNode* Next;
    int PlannedClusters;
    float CU;
    float CUTotal;
} CoalescentNode;


int PtrajClusteringGetAttributeCount(PtrajClustering* This)
{
  if (This->attributeArray || This->attributeArrayTorsion) {
    return PtrajClusteringGetAttributeArrayCount(This);
  } else
  return 3 * PtrajClusteringGetAtomCount(This);  
}

float PtrajClusteringGetAttributeMin(PtrajClustering* This, int AttributeIndex)
{
  int FrameIndex;
  int FrameCount = This->PointCount;
  float MinValue = 200.0;
  float* Attribute;
  int AtomCount = PtrajClusteringGetAtomCount(This); 

  if (This->attributeArray || This->attributeArrayTorsion) {
    return PtrajClusteringGetAttributeArrayMin(This, AttributeIndex);
  }

  trajectoryInfo* trajInfo = This->trajInfo;
  for (FrameIndex = 0; FrameIndex<FrameCount; FrameIndex++)
  {
    if (AttributeIndex < AtomCount)
    {
      Attribute = trajInfo->x + AtomCount*FrameIndex + AttributeIndex;
    }
    else if (AttributeIndex < AtomCount*2)
    {
      Attribute = trajInfo->y + AtomCount*FrameIndex + AttributeIndex - AtomCount;
    }
    else 
    {
      Attribute = trajInfo->z + AtomCount*FrameIndex + AttributeIndex - (AtomCount*2);
    }
    MinValue = min(MinValue, *Attribute);
  }
  return MinValue;
}

float PtrajClusteringGetAttributeMax(PtrajClustering* This, int AttributeIndex)
{
  int FrameIndex;
  int FrameCount = This->PointCount;
  float MaxValue = -200.0;
  float* Attribute;
  int AtomCount = PtrajClusteringGetAtomCount(This); 
  
  if (This->attributeArray || This->attributeArrayTorsion) {
    return PtrajClusteringGetAttributeArrayMax(This, AttributeIndex);
  }
  
  trajectoryInfo* trajInfo = This->trajInfo;
  for (FrameIndex = 0; FrameIndex<FrameCount; FrameIndex++)
  {
    if (AttributeIndex < AtomCount)
    {
      Attribute = trajInfo->x + AtomCount*FrameIndex + AttributeIndex;
    }
    else if (AttributeIndex < AtomCount*2)
    {
      Attribute = trajInfo->y + AtomCount*FrameIndex + AttributeIndex - AtomCount;
    }
    else 
    {
      Attribute = trajInfo->z + AtomCount*FrameIndex + AttributeIndex - (AtomCount*2);
    }
    MaxValue = max(MaxValue, *Attribute);
  }
  return MaxValue;
}
float PtrajClusteringGetAttributeValue(PtrajClustering* This, int PointIndex, int AttributeIndex)
{
  trajectoryInfo* trajInfo = This->trajInfo;
  int AtomCount = PtrajClusteringGetAtomCount(This); 
  float* Attribute;
  /**/
  if (This->attributeArray || This->attributeArrayTorsion) {
    return PtrajClusteringGetAttributeArrayValue(This, PointIndex, AttributeIndex);
  }

  if (AttributeIndex < AtomCount)
  {
    Attribute = trajInfo->x + AtomCount*PointIndex + AttributeIndex;
  }
  else if (AttributeIndex < AtomCount*2)
  {
    Attribute = trajInfo->y + AtomCount*PointIndex + AttributeIndex - AtomCount;
  }
  else 
  {
    Attribute = trajInfo->z + AtomCount*PointIndex + AttributeIndex - (AtomCount*2);
  }
  return *Attribute;
}

int PtrajClusteringGetAttributeArrayCount(PtrajClustering* This)
{
  /*return This->AttributeMatrixCount;  */
  int i = 0;
  int j = 0;

  if (This->attributeArray) 
  	i = This->attributeArray->length;
  if (This->attributeArrayTorsion) 
  	j = This->attributeArrayTorsion->length;
  
  if (i+j) 
    return (i+j);

  return -1;
}

/* Old
float PtrajClusteringGetAttributeArrayValue(PtrajClustering* This, int PointIndex, int AttributeIndex)
{
  double ** Attribute;
  int i = 0;
  int j = 0;

  if (This->attributeArray) 
  	i = This->attributeArray->length;
  if (This->attributeArrayTorsion) 
  	j = This->attributeArrayTorsion->length;

  if (AttributeIndex < i) {
    Attribute = (double **)This->attributeArray->entry;
    return (float)Attribute[AttributeIndex][PointIndex];
  } else {
    AttributeIndex -= i;
    Attribute = (double **)This->attributeArrayTorsion->entry;
    return (float)Attribute[AttributeIndex][PointIndex];
  }
}

float PtrajClusteringGetAttributeArrayMax(PtrajClustering* This, int AttributeIndex)
{
  int FrameIndex;
  int FrameCount = This->PointCount;
  float MaxValue;
  double** Attribute;
  int i = 0;
  int j = 0;

  if (This->attributeArray) 
  	i = This->attributeArray->length;
  if (This->attributeArrayTorsion) 
  	j = This->attributeArrayTorsion->length;

  if (AttributeIndex < i) {
    Attribute = (double **)This->attributeArray->entry;
  } else {
    AttributeIndex -= i;
    Attribute = (double **)This->attributeArrayTorsion->entry;
  }
  MaxValue = Attribute[AttributeIndex][0];
  for (FrameIndex = 0; FrameIndex<FrameCount; FrameIndex++)
  {
    MaxValue = max(MaxValue, Attribute[AttributeIndex][FrameIndex]);
  }
  return MaxValue;
}

float PtrajClusteringGetAttributeArrayMin(PtrajClustering* This, int AttributeIndex)
{
  int FrameIndex;
  int FrameCount = This->PointCount;
  float MinValue;
  double** Attribute;
  int i = 0;
  int j = 0;

  if (This->attributeArray) 
  	i = This->attributeArray->length;
  if (This->attributeArrayTorsion) 
  	j = This->attributeArrayTorsion->length;

  if (AttributeIndex < i) {
    Attribute = (double **)This->attributeArray->entry;
  } else {
    AttributeIndex -= i;
    Attribute = (double **)This->attributeArrayTorsion->entry;
  }
  MinValue = Attribute[AttributeIndex][0];
  for (FrameIndex = 0; FrameIndex<FrameCount; FrameIndex++)
  {
    MinValue = min(MinValue, Attribute[AttributeIndex][FrameIndex]);
  }
  return MinValue;
}
*/

float PtrajClusteringGetAttributeArrayGetMean(PtrajClustering* This, int AttributeIndex)
{
  scalarInfo ** Attribute;
  int i = 0;
  int j = 0;

  if (This->attributeArray) {
    Attribute = (scalarInfo **)This->attributeArray->entry;
    return (float)(Attribute[AttributeIndex]->mean); 
  }
  else /* using the old XYZ attributes */
    return 0;
}

float PtrajClusteringGetAttributeArrayGetStddev(PtrajClustering* This, int AttributeIndex)
{
  scalarInfo ** Attribute;
  int i = 0;
  int j = 0;

  if (This->attributeArray) {
    Attribute = (scalarInfo **)This->attributeArray->entry;
    return (float)(Attribute[AttributeIndex]->stddev); 
  }
  else /* using the old XYZ attributes */
    return 0;
}

char * PtrajClusteringGetAttributeArrayGetName(PtrajClustering* This, int AttributeIndex)
{
  scalarInfo ** Attribute;
  int i = 0;
  int j = 0;

  if (This->attributeArray) {
    Attribute = (scalarInfo **)This->attributeArray->entry;
    return (char *)(Attribute[AttributeIndex]->name); 
  }
  else /* using the old XYZ attributes */
    return 0;
}

float PtrajClusteringGetAttributeArrayValue(PtrajClustering* This, int PointIndex, int AttributeIndex)
{
  scalarInfo ** Attribute;
  int i = 0;
  int j = 0;

  if (This->attributeArray) {
    Attribute = (scalarInfo **)This->attributeArray->entry;
    return (float)(Attribute[AttributeIndex]->value)[PointIndex]; 
  }
  else /* using the old XYZ attributes */
    return PtrajClusteringGetAttributeValue(This, PointIndex, AttributeIndex);
}

int PtrajClusteringIsAttributeTorsion(PtrajClustering* This, int AttributeIndex)
{
  scalarInfo ** Attribute;
  if (This->attributeArray) {
    Attribute = (scalarInfo **)This->attributeArray->entry;
    return (Attribute[AttributeIndex]->cos?1:0);
  } else {
    return 0;
  }
}
float PtrajClusteringGetAttributeArrayCos(PtrajClustering* This, int PointIndex, int AttributeIndex)
{
  scalarInfo ** Attribute;
  int i = 0;
  int j = 0;

  Attribute = (scalarInfo **)This->attributeArray->entry;
  return (float)(Attribute[AttributeIndex]->cos)[PointIndex]; 
}

float PtrajClusteringGetAttributeArraySin(PtrajClustering* This, int PointIndex, int AttributeIndex)
{
  scalarInfo ** Attribute;
  int i = 0;
  int j = 0;

  Attribute = (scalarInfo **)This->attributeArray->entry;
  return (float)(Attribute[AttributeIndex]->sin)[PointIndex]; 
}

float PtrajClusteringGetAttributeArrayMax(PtrajClustering* This, int AttributeIndex)
{
  int FrameIndex;
  int FrameCount = This->PointCount;
  float MaxValue;

  MaxValue = PtrajClusteringGetAttributeArrayValue(This, 0, AttributeIndex);
  for (FrameIndex = 0; FrameIndex<FrameCount; FrameIndex++)
  {
    MaxValue = max(MaxValue, PtrajClusteringGetAttributeArrayValue(This, FrameIndex, AttributeIndex));
  }
  return MaxValue;
}

float PtrajClusteringGetAttributeArrayMin(PtrajClustering* This, int AttributeIndex)
{
  int FrameIndex;
  int FrameCount = This->PointCount;
  float MinValue;

  MinValue = PtrajClusteringGetAttributeArrayValue(This, 0, AttributeIndex);
  for (FrameIndex = 0; FrameIndex<FrameCount; FrameIndex++)
  {
    MinValue = min(MinValue, PtrajClusteringGetAttributeArrayValue(This, FrameIndex, AttributeIndex));
  }
  return MinValue;
}

PtrajSOMNode** PtrajClusteringInitSOMNodes(PtrajClustering* This, int ClusterCount, int map)
{
  PtrajSOMNode** SOMNodes;
  int SOMNodeIndex;
  int AttributeCount;
  int AttributeIndex;
  float AttributeMin;  
  float AttributeMax;
  float RandValue;
  float Value;
  int IsTor;
  
  /* Allocate our SOMNodes, and initialize them: */
  SOMNodes = (PtrajSOMNode**)SafeMalloc(__FILE__, __LINE__, sizeof(PtrajSOMNode*) * ClusterCount);
  memset(SOMNodes,0,sizeof(PtrajSOMNode*) * ClusterCount);
  AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
  for (SOMNodeIndex=0; SOMNodeIndex < ClusterCount; SOMNodeIndex++)
  {
    SOMNodes[SOMNodeIndex] = (PtrajSOMNode*)SafeMalloc(__FILE__, __LINE__, sizeof(PtrajSOMNode));
    memset(SOMNodes[SOMNodeIndex],0,sizeof(PtrajSOMNode));
    SOMNodes[SOMNodeIndex]->Attributes = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Attributes,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->Means = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Means,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->Stddevs = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Stddevs,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->Cos = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Cos,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->Sin = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Sin,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->Cos2 = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Cos2,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->Sin2 = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->Sin2,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->PointCount = 0;
    SOMNodes[SOMNodeIndex]->MeanCos = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->MeanCos,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->MeanSin = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    memset(SOMNodes[SOMNodeIndex]->MeanSin,0,sizeof(float) * AttributeCount);
    SOMNodes[SOMNodeIndex]->PointCount = 0;
  }
  for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
  {
    AttributeMin = ClusteringGetAttributeMin( (Clustering *) This,AttributeIndex);
    AttributeMax = ClusteringGetAttributeMax( (Clustering *) This,AttributeIndex);
    IsTor = PtrajClusteringIsAttributeTorsion(This, AttributeIndex);
    for (SOMNodeIndex=0; SOMNodeIndex < ClusterCount; SOMNodeIndex++)
    {
      RandValue = rand()/(double)RAND_MAX;
      Value = AttributeMax - AttributeMin;
      SOMNodes[SOMNodeIndex]->Attributes[AttributeIndex] = AttributeMin + Value*RandValue;
      if (IsTor) {
        /*SOMNodes[SOMNodeIndex]->Attributes[AttributeIndex] = 360*RandValue;*/
        SOMNodes[SOMNodeIndex]->Cos[AttributeIndex] = cos(SOMNodes[SOMNodeIndex]->Attributes[AttributeIndex]/RADDEG);
        SOMNodes[SOMNodeIndex]->Sin[AttributeIndex] = sin(SOMNodes[SOMNodeIndex]->Attributes[AttributeIndex]/RADDEG);
      }
    }
  }
  ClusteringSetSOMMap(This, (SOMNode**)SOMNodes, ClusterCount, map); /* Setup map */
  return SOMNodes;
}

void PtrajClusteringFreeSOMNodes(Clustering* This,PtrajSOMNode** SOMNodes,int ClusterCount)
{
  int SOMNodeIndex;
  for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount; SOMNodeIndex++)
    {
      safe_free(SOMNodes[SOMNodeIndex]->Attributes);
      safe_free(SOMNodes[SOMNodeIndex]->Neighbors);
      safe_free(SOMNodes[SOMNodeIndex]->Means);
      safe_free(SOMNodes[SOMNodeIndex]->Stddevs);
      safe_free(SOMNodes[SOMNodeIndex]->Cos);
      safe_free(SOMNodes[SOMNodeIndex]->Sin);
      safe_free(SOMNodes[SOMNodeIndex]->Cos2);
      safe_free(SOMNodes[SOMNodeIndex]->Sin2);
      safe_free(SOMNodes[SOMNodeIndex]->MeanCos);
      safe_free(SOMNodes[SOMNodeIndex]->MeanSin);
      safe_free(SOMNodes[SOMNodeIndex]);
    }
  safe_free(SOMNodes);
}

void PtrajClusteringPrintTransformMap(PtrajClustering* This)
{
    int* TransformMatrix;
    int i, j, PointIndex;
    int PointCount = This->PointCount;
    int ClusterCount =  This->ClusterCount;
    ClusterNode* Node;

	/* One Clustering should have been done. So calculate the transform matrix. */
    TransformMatrix = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int*) * ClusterCount*ClusterCount);
    for (i = 0; i < ClusterCount * ClusterCount; i++)
    	TransformMatrix[i] = 0;
    
    for (Node = This->Head, i = 0; Node; Node = Node->pNext,i++) {
    	if (Node->Cluster->Mask[0] == 1) break;
    }
    if (i >= ClusterCount) {
    	fprintf(stderr, "GetSOMMap(): Point 0 does not belong to any Cluster!\n");
        exit(1);
    }
    for (PointIndex = 1; PointIndex < PointCount; PointIndex++){
    	for (Node = This->Head, j = 0; Node; Node = Node->pNext,j++) {
	    	if (Node->Cluster->Mask[PointIndex] == 1) break;
	    }
        if (j >= ClusterCount) {
    		fprintf(stderr, "GetSOMMap(): Point %d does not belong to any Cluster!\n", PointIndex);
        	exit(1);
    	}
    	TransformMatrix[i*ClusterCount+j]++;
        i = j;
    }
    
/*    if (This->action->darg2 > VERBOSE_INTERMEDIATE_CLUSTER) {*/
    if (This->action->darg2) {

      if (prnlev > 5) {
	fprintf(stdout, "Cluster transformation matrix:\nCluster                     ");
        for (i = 0; i < ClusterCount; i++) {
        	fprintf(stdout, "%8d", i);
        }
        fprintf(stdout, "\n");

        for (i = 0, Node = This->Head; i < ClusterCount; i++, Node= Node->pNext) {
        	fprintf(stdout, "Cluster %4d (%5d points) ", i, ClusterCountPointsInCluster(Node->Cluster));
            for (j = 0; j < ClusterCount; j++) {
            	fprintf(stdout, "%8d", TransformMatrix[i*ClusterCount+j]);
            }
            fprintf(stdout, "\n");
        }
      }        
    }

    safe_free(TransformMatrix); 

}


PtrajSOMNode** PtrajClusteringGetSOMMap(PtrajClustering* This)
{
	int SOMNodeIndex, PointIndex;
    int prevIndec, nextIndex, Index;
    int* TransformMatrix;
	int ClusterCount = This->ClusterCount;
    int PointCount = This->PointCount;
    int i, j;
    ClusterNode* Node;
    PtrajSOMNode** SOMNodes;
    int AttributeCount, AttributeIndex;
    int IsTor;
    
	/* One Clustering should have been done. So calculate the transform matrix. */
    TransformMatrix = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int*) * ClusterCount*ClusterCount);
    for (i = 0; i < ClusterCount * ClusterCount; i++)
    	TransformMatrix[i] = 0;
    
    for (Node = This->Head, i = 0; Node; Node = Node->pNext,i++) {
    	if (Node->Cluster->Mask[0] == 1) break;
    }
    if (i >= ClusterCount) {
    	fprintf(stderr, "GetSOMMap(): Point 0 does not belong to any Cluster!\n");
        exit(1);
    }
    for (PointIndex = 1; PointIndex < PointCount; PointIndex++){
    	for (Node = This->Head, j = 0; Node; Node = Node->pNext,j++) {
	    	if (Node->Cluster->Mask[PointIndex] == 1) break;
	    }
        if (j >= ClusterCount) {
    		fprintf(stderr, "GetSOMMap(): Point %d does not belong to any Cluster!\n", PointIndex);
        	exit(1);
    	}
    	TransformMatrix[i*ClusterCount+j]++;
        i = j;
    }
    
/*    if (This->action->darg2 > VERBOSE_INTERMEDIATE_CLUSTER) {*/
    if (This->action->darg2) {
    	fprintf(stdout, "Cluster transformation matrix:\nCluster                     ");
        for (i = 0; i < ClusterCount; i++) {
        	fprintf(stdout, "%8d", i);
        }
        fprintf(stdout, "\n");
        for (i = 0, Node = This->Head; i < ClusterCount; i++, Node= Node->pNext) {
        	fprintf(stdout, "Cluster %4d (%5d points) ", i, ClusterCountPointsInCluster(Node->Cluster));
            for (j = 0; j < ClusterCount; j++) {
            	fprintf(stdout, "%8d", TransformMatrix[i*ClusterCount+j]);
            }
            fprintf(stdout, "\n");
        }
        
    }
    
    /* Allocate our SOMNodes, and initialize them: */
    SOMNodes = (PtrajSOMNode**)SafeMalloc(__FILE__, __LINE__, sizeof(PtrajSOMNode*) * ClusterCount);
    memset(SOMNodes,0,sizeof(PtrajSOMNode*) * ClusterCount);
    AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
    for (SOMNodeIndex=0; SOMNodeIndex < ClusterCount; SOMNodeIndex++)
    {
	    SOMNodes[SOMNodeIndex] = (PtrajSOMNode*)SafeMalloc(__FILE__, __LINE__, sizeof(PtrajSOMNode));
    	memset(SOMNodes[SOMNodeIndex],0,sizeof(PtrajSOMNode));
	    SOMNodes[SOMNodeIndex]->Attributes = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    	memset(SOMNodes[SOMNodeIndex]->Attributes,0,sizeof(float) * AttributeCount);
	    SOMNodes[SOMNodeIndex]->Means = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    	memset(SOMNodes[SOMNodeIndex]->Means,0,sizeof(float) * AttributeCount);
	    SOMNodes[SOMNodeIndex]->Stddevs = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    	memset(SOMNodes[SOMNodeIndex]->Stddevs,0,sizeof(float) * AttributeCount);
	    SOMNodes[SOMNodeIndex]->Cos = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    	memset(SOMNodes[SOMNodeIndex]->Cos,0,sizeof(float) * AttributeCount);
	    SOMNodes[SOMNodeIndex]->Sin = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    	memset(SOMNodes[SOMNodeIndex]->Sin,0,sizeof(float) * AttributeCount);
	    SOMNodes[SOMNodeIndex]->Cos2 = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    	memset(SOMNodes[SOMNodeIndex]->Cos2,0,sizeof(float) * AttributeCount);
	    SOMNodes[SOMNodeIndex]->Sin2 = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
	    memset(SOMNodes[SOMNodeIndex]->Sin2,0,sizeof(float) * AttributeCount);
    	SOMNodes[SOMNodeIndex]->PointCount = 0;
	    SOMNodes[SOMNodeIndex]->MeanCos = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    	memset(SOMNodes[SOMNodeIndex]->MeanCos,0,sizeof(float) * AttributeCount);
	    SOMNodes[SOMNodeIndex]->MeanSin = (float*)SafeMalloc(__FILE__, __LINE__, sizeof(float) * AttributeCount);
    	memset(SOMNodes[SOMNodeIndex]->MeanSin,0,sizeof(float) * AttributeCount);
	    SOMNodes[SOMNodeIndex]->PointCount = 0;
    }
    
   	for (Node = This->Head, j = 0; Node; Node = Node->pNext,j++) {
    	for (PointIndex = 1, i = 0; PointIndex < PointCount; PointIndex++){ 
           	if (Node->Cluster->Mask[PointIndex] == 0) continue; 
		    i++;
            for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
   			{
   			 	IsTor = PtrajClusteringIsAttributeTorsion(This, AttributeIndex);
                if (!IsTor) {
	               	SOMNodes[j]->Attributes[AttributeIndex] += 
			  ClusteringGetAttributeValue( (Clustering *) This,PointIndex,AttributeIndex);
                } else {
	               	SOMNodes[j]->Cos[AttributeIndex] += 
			  PtrajClusteringGetAttributeArrayCos(This,PointIndex,AttributeIndex);
	               	SOMNodes[j]->Sin[AttributeIndex] += 
			  PtrajClusteringGetAttributeArraySin(This,PointIndex,AttributeIndex);
                }
            }
	    }
        for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
		{
		 	IsTor = PtrajClusteringIsAttributeTorsion(This, AttributeIndex);
            if (!IsTor) {
	           	SOMNodes[j]->Attributes[AttributeIndex] /= i;
            } else {
	           	SOMNodes[j]->Cos[AttributeIndex] /= i;
	           	SOMNodes[j]->Sin[AttributeIndex] /= i;
		        SOMNodes[j]->Attributes[AttributeIndex] = atan2(SOMNodes[j]->Sin[AttributeIndex], SOMNodes[j]->Cos[AttributeIndex]) * RADDEG;
		        if (SOMNodes[j]->Attributes[AttributeIndex]<0)
		        	SOMNodes[j]->Attributes[AttributeIndex] += 360;
            }
		}        
    }
    
  	for (SOMNodeIndex = 0; SOMNodeIndex < ClusterCount; SOMNodeIndex++) {
    	safe_free(SOMNodes[SOMNodeIndex]->Neighbors);
        SOMNodes[SOMNodeIndex]->NeighborCount = 0;
        for (j = 0; j < ClusterCount; j++) {
        	if (j == SOMNodeIndex) continue;
        	if (TransformMatrix[SOMNodeIndex*ClusterCount+j] > 0) SOMNodes[SOMNodeIndex]->NeighborCount += 1;
        }
        SOMNodes[SOMNodeIndex]->Neighbors = (SOMNode **)
	  malloc(sizeof(SOMNode) * SOMNodes[SOMNodeIndex]->NeighborCount);
        for (j = 0, i = 0; j < ClusterCount; j++) {
        	if (j == SOMNodeIndex) continue;
        	if (TransformMatrix[SOMNodeIndex*ClusterCount+j] > 0) 
		  SOMNodes[SOMNodeIndex]->Neighbors[i++] = (SOMNode *) SOMNodes[j];
        }
    }
    return(SOMNodes);        
}


float PtrajClusteringSOMGetDistance(PtrajClustering* This, PtrajSOMNode* SOM, int PointIndex)
{
  float Distance = 0;
  float AttDistance;
  int AttributeIndex;
  int AttributeCount;
  int IsTor;
  float sin, cos;
  
  AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
  for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
  {
    IsTor = PtrajClusteringIsAttributeTorsion(This, AttributeIndex);
    if (!IsTor) {
      AttDistance = ClusteringGetAttributeValue( (Clustering *) This,PointIndex,AttributeIndex) 
	- SOM->Attributes[AttributeIndex];
      Distance += (AttDistance*AttDistance);
    } else {
      sin = PtrajClusteringGetAttributeArraySin(This,PointIndex,AttributeIndex);
      cos = PtrajClusteringGetAttributeArrayCos(This,PointIndex,AttributeIndex);
      AttDistance = atan2(sin, cos);
      AttDistance = AttDistance * RADDEG;
      if (AttDistance < 0) 
        /*AttDistance += 360 / RADDEG; /* RADDEG = 180/pi */
      /*AttDistance = AttDistance - SOM->Attributes[AttributeIndex] / RADDEG;*/
        AttDistance += 360;
      AttDistance = 180 - abs(180- abs(AttDistance - SOM->Attributes[AttributeIndex]));
      Distance += (AttDistance*AttDistance); /* Distance is the sum of the difference in radian, not in degree. */
    }
  }
  return sqrt(Distance/AttributeCount);

}

void PtrajClusteringSOMLearn(PtrajClustering* This, PtrajSOMNode* SOM, float* TheLearningRate, int PointIndex)
{
  float Distance;
  float AttDistance;    
  float AttValue;
  float NewValue; 
  float LearningRate = *TheLearningRate;
  int AttributeIndex;
  int AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
  int IsTor;
  float sin, cos;
  
  for (AttributeIndex = 0; AttributeIndex < AttributeCount; AttributeIndex++)
  {
    IsTor = PtrajClusteringIsAttributeTorsion(This, AttributeIndex);
    if (!IsTor) {
      AttValue = ClusteringGetAttributeValue( (Clustering *) This,PointIndex,AttributeIndex);
      AttDistance = AttValue - SOM->Attributes[AttributeIndex];
      SOM->Attributes[AttributeIndex] += AttDistance * LearningRate;
    } else {
      sin = PtrajClusteringGetAttributeArraySin(This,PointIndex,AttributeIndex);
      cos = PtrajClusteringGetAttributeArrayCos(This,PointIndex,AttributeIndex);
      SOM->Cos[AttributeIndex] += cos * LearningRate;
      SOM->Sin[AttributeIndex] += sin * LearningRate;
      SOM->Attributes[AttributeIndex] = atan2(SOM->Sin[AttributeIndex], SOM->Cos[AttributeIndex]) * RADDEG;
      if (SOM->Attributes[AttributeIndex]<0)
        SOM->Attributes[AttributeIndex] += 360;
    
    }
  }
} 

void PtrajClusteringSOMExtractDifference(PtrajClustering* This, PtrajSOMNode** SOMNodes) 
{
  ClusterNode* Node;
  Cluster* NewCluster;
  int PointInCluster;
  int ClusterCount = This->ClusterCount;
  int AttributeCount = PtrajClusteringGetAttributeCount(This);
  int i, j, k, diff;
  int* AttributeDifference;
  
  AttributeDifference = SafeMalloc(__FILE__, __LINE__, sizeof(int*) * AttributeCount);
  memset(AttributeDifference, 0, sizeof(int*) * AttributeCount);
  
  for (i = 0; i < ClusterCount; i++) {
    PointInCluster = SOMNodes[i]->PointCount;
    if (PointInCluster < 10) {
      fprintf(stdout, "ClusterCount in Cluster %d is only %d, less than 10.\n", i, ClusterCount);
      continue;
    }
    for (j = 0; j < i; j++) {
      PointInCluster = SOMNodes[j]->PointCount;
      if (PointInCluster < 10) {
        fprintf(stdout, "ClusterCount in Cluster %d is only %d, less than 10.\n", i, ClusterCount);
        continue;
      }
      diff = 0;
      for (k = 0; k <AttributeCount; k++) {
        if (((SOMNodes[i]->Means[k] - SOMNodes[j]->Means[k]) / (SOMNodes[i]->Stddevs[k] + SOMNodes[j]->Stddevs[k])) > 2) { /* Difference is greater than 2 times the sum of 2 standard deviation. */
          if (diff == 0) fprintf(stdout, "Cluster %d (%d points) and Cluster %d (%d points) differs significantly in the following attribute(s).\n", i, SOMNodes[i]->PointCount, j, SOMNodes[j]->PointCount);
          diff = 1;
          fprintf(stdout, "  attribute %i(%s) differs: %.3f+/-%.3f, %.3f+/-%.3f.   ", k, PtrajClusteringGetAttributeArrayGetName(This, k), SOMNodes[i]->Means[k], SOMNodes[i]->Stddevs[k], SOMNodes[j]->Means[k], SOMNodes[j]->Stddevs[k]);
          fprintf(stdout, "for all points: %.3f+/-%.3f.\n", PtrajClusteringGetAttributeArrayGetMean(This, k), PtrajClusteringGetAttributeArrayGetStddev(This, k));
          AttributeDifference[k]++;
        }
      } 
    }
  }
  
  for (k = 0; k <AttributeCount; k++) {
    if (AttributeDifference[k])
      fprintf(stdout, "attribute %d (%s)--%d\n", k, PtrajClusteringGetAttributeArrayGetName(This, k), AttributeDifference[k]);
  }
    
    
}


/* Helper function for SOM clustering: Once we've generated and trained the map, it's time 
   to convert the map into clusters */
void ClusteringConvertSOMNodesToClusters(Clustering* This, SOMNode** SOMNodes, int ClusterCount)
{
  ClusterNode* Node;
  Cluster* NewCluster;
  int ClusterIndex;
  int PointIndex;
  int WinnerIndex;
  float X,Y;
  int i;
  int PointCount;
  
  /**/

    /*
     *  First, build the empty clusters: 
     */
  for (ClusterIndex = 0; ClusterIndex < ClusterCount; ClusterIndex++)
    {
      NewCluster = ClusteringGetNewCluster(This);
      ClusteringAddCluster(This,NewCluster);
      if (prnlev > 5) {
	fprintf(stdout,"SOMNode %d is at (%.2f", ClusterIndex,SOMNodes[ClusterIndex]->Attributes[0]);
	for (i = 1; i < PtrajClusteringGetAttributeArrayCount( (PtrajClustering *) This) - 1; i++) 
	  {
	    if (i < 6)
              fprintf(stdout, ", %.2f",SOMNodes[ClusterIndex]->Attributes[i]);
	    else 
	      continue;
	  }
	if (i < 6) 
          fprintf(stdout, ")\n"); 
	if (i == 6) 
          fprintf(stdout, ", %.2f)\n",SOMNodes[ClusterIndex]->Attributes[i]); 
	if (i > 6) 
          fprintf(stdout, " ... %.2f)\n",SOMNodes[ClusterIndex]->Attributes[i]); 
      }
    }

    /*
     *  Now iterate over all points, find the WINNER for each point, and add the point as a member
     *  of the corresponding Cluster 
     */
  for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
    {
      WinnerIndex = ClusteringSOMFindWinner(This,SOMNodes,ClusterCount,PointIndex);
      if (prnlev > 5) 
	fprintf(stdout,"Point %d (%.2f", PointIndex, ClusteringGetAttributeValue(This,PointIndex,0)); 
      SOMNodes[WinnerIndex]->PointCount += 1;
      for (i = 0; i < PtrajClusteringGetAttributeArrayCount( (PtrajClustering *) This); i++) 
        {
          SOMNodes[WinnerIndex]->Means[i] += ClusteringGetAttributeValue(This,PointIndex,i);
      	  SOMNodes[WinnerIndex]->Stddevs[i] += ClusteringGetAttributeValue(This,PointIndex,i) *  
	    ClusteringGetAttributeValue(This,PointIndex,i);       
        }

      if (prnlev > 5) {
	for (i = 1; i < PtrajClusteringGetAttributeArrayCount( (PtrajClustering *) This) - 1; i++) 
	  {
	    if (i < 6)
              fprintf(stdout, ", %.2f",ClusteringGetAttributeValue(This,PointIndex,i));
	    else 
	      continue;
	  }
	if (i <= 6) 
          fprintf(stdout, ", %.2f) closest to SOMNode %d",
		  ClusteringGetAttributeValue(This,PointIndex,i), WinnerIndex); 
	if (i > 6) 
          fprintf(stdout, " ... %.2f) closest to SOMNode %d",
		  ClusteringGetAttributeValue(This,PointIndex,i), WinnerIndex); 
	fprintf(stdout, " with the distance of %f.\n",
		ClusteringSOMGetDistance(This,SOMNodes[WinnerIndex],PointIndex));
      }

      for (Node = This->Head,ClusterIndex=0; ClusterIndex<ClusterCount; Node = Node->pNext,ClusterIndex++)
        {
          if (ClusterIndex == WinnerIndex)
              /*ClusterAddMember(Node->Cluster,PointIndex);*/
              PtrajClusterAddMember( (PtrajCluster *) Node->Cluster,PointIndex);
        }
    }
  for (ClusterIndex = 0; ClusterIndex < ClusterCount; ClusterIndex++)
    {
      PointCount = SOMNodes[ClusterIndex]->PointCount;
      if (prnlev > 5) fprintf(stdout, "SOMNode %d (%d points): \n", ClusterIndex, PointCount);
      for (i = 0; i < PtrajClusteringGetAttributeArrayCount( (PtrajClustering *) This); i++) 
        {
          if (PointCount == 0) continue;
          if (PtrajClusteringIsAttributeTorsion( (PtrajClustering *) This, i))
            SOMNodes[ClusterIndex]->Means[i] = SOMNodes[ClusterIndex]->Means[i] / PointCount; 
            SOMNodes[ClusterIndex]->Stddevs[i] = SOMNodes[ClusterIndex]->Stddevs[i] / PointCount - (SOMNodes[ClusterIndex]->Means[i]) * (SOMNodes[ClusterIndex]->Means[i]);
          SOMNodes[ClusterIndex]->Stddevs[i] = sqrt(SOMNodes[ClusterIndex]->Stddevs[i]);
          if (prnlev > 5) fprintf(stdout, "  mean is %.2f, stddev is %.2f. \n",
				  SOMNodes[ClusterIndex]->Means[i], SOMNodes[ClusterIndex]->Stddevs[i]);
        }
    }
    
  /* Remove cluster which did not contain any points */  
  for (Node = This->Head,ClusterIndex=0; ClusterIndex < ClusterCount; ClusterIndex++)
    {
      if (ClusterCountPointsInCluster(Node->Cluster) == 0 )
        {
          if (Node->pNext) {
            Node = Node->pNext;
            ClusteringRemoveCluster(This, Node->pPrev->Cluster);
          } else {
            ClusteringRemoveCluster(This, Node->Cluster);
          }
          fprintf(stdout, "  SOM: cluster %d does not contain any points. Deleting it. \n", ClusterIndex);
          continue;
        }
      ClusterFindBestP2C( (PtrajCluster *) Node->Cluster);
      Node = Node->pNext;  
    }
}

/* Helper function for SOM clustering: Once we've generated and trained the map, it's time 
   to convert the map into clusters */
void PtrajClusteringConvertSOMNodesToClusters(Clustering* This, PtrajSOMNode** SOMNodes, int ClusterCount)
{
  ClusterNode* Node;
  Cluster* NewCluster;
  int ClusterIndex;
  int PointIndex;
  int WinnerIndex;
  float X,Y;
  int i;
  int PointCount;
  
  /**/
  /* First, build the empty clusters: */
  fprintf(stdout,"\n\n--------------------------------------------------------------------\n");
  for (ClusterIndex = 0; ClusterIndex < ClusterCount; ClusterIndex++)
    {
      NewCluster = ClusteringGetNewCluster(This);
      ClusteringAddCluster(This,NewCluster);

      if (prnlev > 5) {
      fprintf(stdout,"SOMNode %d is at (%.5f", ClusterIndex,SOMNodes[ClusterIndex]->Attributes[0]);
      for (i = 1; i < PtrajClusteringGetAttributeCount( (PtrajClustering *) This) - 1; i++) 
        {
          if (i < 6)
              fprintf(stdout, ", %.2f",SOMNodes[ClusterIndex]->Attributes[i]);
          else 
             continue;
        }
      if (i <= 6) 
          fprintf(stdout, ", %.2f)\n",SOMNodes[ClusterIndex]->Attributes[i]); 
      if (i > 6) 
          fprintf(stdout, " ... %.2f)\n",SOMNodes[ClusterIndex]->Attributes[i]); 
      }
    }
  /* Now iterate over all points, find the WINNER for each point, and add the point as a member
     of the corresponding Cluster */
  for (PointIndex = 0; PointIndex < This->PointCount; PointIndex++)
    {
      WinnerIndex = ClusteringSOMFindWinner(This,SOMNodes,ClusterCount,PointIndex);
      if (prnlev > 5) fprintf(stdout,"Point %d (%.2f", PointIndex, ClusteringGetAttributeValue(This,PointIndex,0)); 
      SOMNodes[WinnerIndex]->PointCount += 1;
      for (i = 0; i < PtrajClusteringGetAttributeCount( (PtrajClustering *) This); i++) 
        {
          if (PtrajClusteringIsAttributeTorsion( (PtrajClustering *) This, i))
            {
              SOMNodes[WinnerIndex]->MeanSin[i] += 
		PtrajClusteringGetAttributeArraySin( (PtrajClustering *) This,PointIndex,i);
              SOMNodes[WinnerIndex]->MeanCos[i] += 
		PtrajClusteringGetAttributeArrayCos( (PtrajClustering *) This,PointIndex,i);
              SOMNodes[WinnerIndex]->Sin2[i] += 
		PtrajClusteringGetAttributeArraySin( (PtrajClustering *) This,PointIndex,i) * 
		PtrajClusteringGetAttributeArraySin( (PtrajClustering *) This,PointIndex,i);
              SOMNodes[WinnerIndex]->Cos2[i] += 
		PtrajClusteringGetAttributeArrayCos( (PtrajClustering *) This,PointIndex,i) * 
		PtrajClusteringGetAttributeArrayCos( (PtrajClustering *) This,PointIndex,i);
              
            } else
            {
              SOMNodes[WinnerIndex]->Means[i] += ClusteringGetAttributeValue(This,PointIndex,i);
	      SOMNodes[WinnerIndex]->Stddevs[i] += 
		ClusteringGetAttributeValue(This,PointIndex,i) * 
		ClusteringGetAttributeValue(This,PointIndex,i);
            } 
        }

      if (prnlev > 5) {
      for (i = 1; i < PtrajClusteringGetAttributeCount( (PtrajClustering *) This) - 1; i++) 
        {
          if (i < 6)
              fprintf(stdout, ", %.2f",ClusteringGetAttributeValue(This,PointIndex,i));
          else 
             continue;
        }
      if (i <= 6) 
          fprintf(stdout, ", %.2f) closest to SOMNode %d",
		  ClusteringGetAttributeValue(This,PointIndex,i), WinnerIndex); 
      if (i > 6) 
          fprintf(stdout, " ... %.2f) closest to SOMNode %d",
		  ClusteringGetAttributeValue(This,PointIndex,i), WinnerIndex); 
      fprintf(stdout, " with the distance of %f.\n", 
	      ClusteringSOMGetDistance( (Clustering *) This, 
					(SOMNode *) SOMNodes[WinnerIndex], PointIndex));
      }

      for (Node = This->Head,ClusterIndex=0; ClusterIndex<ClusterCount; Node = Node->pNext,ClusterIndex++)
        {
          if (ClusterIndex == WinnerIndex)
              /*ClusterAddMember(Node->Cluster,PointIndex);*/
              PtrajClusterAddMember( (PtrajCluster *) Node->Cluster,PointIndex);
        }
    }
  for (ClusterIndex = 0; ClusterIndex < ClusterCount; ClusterIndex++)
    {
      PointCount = SOMNodes[ClusterIndex]->PointCount;
      if (prnlev > 5) fprintf(stdout, "SOMNode %d (%d points): \n", ClusterIndex, PointCount);
      for (i = 0; i < PtrajClusteringGetAttributeArrayCount( (PtrajClustering *) This); i++) 
        {
          if (PointCount == 0) continue;
          if (PtrajClusteringIsAttributeTorsion( (PtrajClustering *) This, i))
            {
              SOMNodes[ClusterIndex]->Means[i] = atan2(SOMNodes[ClusterIndex]->MeanSin[i], SOMNodes[ClusterIndex]->MeanCos[i]) *RADDEG; 
              if (SOMNodes[ClusterIndex]->Means[i] < 0) SOMNodes[ClusterIndex]->Means[i] += 360; 
              SOMNodes[ClusterIndex]->Stddevs[i] = (SOMNodes[ClusterIndex]->Sin2[i] + SOMNodes[ClusterIndex]->Cos2[i])/PointCount - (SOMNodes[ClusterIndex]->MeanSin[i] * SOMNodes[ClusterIndex]->MeanSin[i] + SOMNodes[ClusterIndex]->MeanCos[i]* SOMNodes[ClusterIndex]->MeanCos[i]) /( PointCount * PointCount) ; 
              SOMNodes[ClusterIndex]->Stddevs[i] = acos(1 - SOMNodes[ClusterIndex]->Stddevs[i]/2)*RADDEG;               
            } else {
              SOMNodes[ClusterIndex]->Means[i] = SOMNodes[ClusterIndex]->Means[i] / PointCount; 
          	  SOMNodes[ClusterIndex]->Stddevs[i] = SOMNodes[ClusterIndex]->Stddevs[i] / PointCount - (SOMNodes[ClusterIndex]->Means[i]) * (SOMNodes[ClusterIndex]->Means[i]);
              SOMNodes[ClusterIndex]->Stddevs[i] = sqrt(SOMNodes[ClusterIndex]->Stddevs[i]);
          }
          char *AttrName = PtrajClusteringGetAttributeArrayGetName( (PtrajClustering *) This, i);
          if (AttrName && prnlev > 5) fprintf(stdout, "  attribute %4i: %-20s--", i, AttrName);
          if (prnlev > 5) fprintf(stdout, "  mean is %.2f, stddev is %.2f. \n", 
				  SOMNodes[ClusterIndex]->Means[i], SOMNodes[ClusterIndex]->Stddevs[i]);
        }
    }
    
  /* Remove cluster which did not contain any points */  
  for (Node = This->Head,ClusterIndex=0; ClusterIndex < ClusterCount; ClusterIndex++)
    {
      if (ClusterCountPointsInCluster(Node->Cluster) == 0 )
        {
          if (Node->pNext) {
            Node = Node->pNext;
            ClusteringRemoveCluster(This, Node->pPrev->Cluster);
          } else {
            ClusteringRemoveCluster(This, Node->Cluster);
          }
          fprintf(stdout, "Cluster %d does not contain any points. Delete it. \n", ClusterIndex);
          continue;
        }
      ClusterFindBestP2C( (PtrajCluster *) Node->Cluster);
      Node = Node->pNext;  
    }
    
    if (prnlev >6) PtrajClusteringSOMExtractDifference( (PtrajClustering *) This, SOMNodes); /* print out significant differences in attributes between clusters. */
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
void PtrajClusteringClusterSOM(PtrajClustering* This, int ClusterCount, int map)
{
  PtrajSOMNode** SOMNodes;
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
  SOMNodes = PtrajClusteringInitSOMNodes(This,ClusterCount,map);
  PointProcessed = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * This->PointCount);
  /* Iterate! */ 
  WinnerLearningRate = 0.2;
  NeighborLearningRate = 0.1;
  for (IterationIndex = 0; IterationIndex < 500; IterationIndex++)
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
        PtrajClusteringSOMLearn(This,SOMNodes[WinnerIndex],&WinnerLearningRate,PointIndex);
        for (NeighborIndex = 0; NeighborIndex < SOMNodes[WinnerIndex]->NeighborCount; NeighborIndex++)
          {
            PtrajClusteringSOMLearn(This, (PtrajSOMNode *) 
				    SOMNodes[WinnerIndex]->Neighbors[NeighborIndex],
				    &NeighborLearningRate,PointIndex);
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
  PtrajClusteringConvertSOMNodesToClusters( (Clustering *) This, SOMNodes,ClusterCount);
  PtrajClusteringFreeSOMNodes( (Clustering *) This,SOMNodes,ClusterCount);
  safe_free(PointProcessed);
}

/* This method (SOM2) will do a SOM-isolate first, then determine the SOM map, then run SOM algorithm based on the SOM map. */
void PtrajClusteringClusterSOM2(PtrajClustering* This, int ClusterCount)
{
  PtrajSOMNode** SOMNodes;
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
  int PointCount = This->PointCount;
  ClusterNode* Node;
  /**/
  Randomize();
  /*SOMNodes = PtrajClusteringInitSOMNodes(This,ClusterCount,0); /* SOM_ISOLATE 0 */
  PointProcessed = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * This->PointCount);
  WinnerLearningRate = 0.2;
  NeighborLearningRate = 0.1;
  /* Get the SOM Map from the existing clustering. */
  SOMNodes = PtrajClusteringGetSOMMap(This);
  /* Reset the clusters. And remove those clusters. */
  ClusterCount = This->ClusterCount; /* Save ClusterCount if old clusters removed. */
  for (PointIndex = 0; PointIndex < PointCount; PointIndex++){
  	for (Node = This->Head; Node; Node = Node->pNext) {
	  	Node->Cluster->Mask[PointIndex] = 0;
	}
  }
  fprintf(stdout, "Delete old clusters.\n");
  ClusteringRemoveEmptyCluster( (Clustering *) This);
  PtrajClusteringConvertSOMNodesToClusters( (Clustering *) This,SOMNodes,ClusterCount); 
  for (PointIndex = 0; PointIndex < PointCount; PointIndex++){
  	for (Node = This->Head; Node; Node = Node->pNext) {
	  	Node->Cluster->Mask[PointIndex] = 0;
	}
  }
  ClusteringRemoveEmptyCluster( (Clustering *) This);
  /* Iterate! */ 
  for (IterationIndex = 0; IterationIndex < 500; IterationIndex++)
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
        PtrajClusteringSOMLearn(This,SOMNodes[WinnerIndex],&WinnerLearningRate,PointIndex);
        for (NeighborIndex = 0; NeighborIndex < SOMNodes[WinnerIndex]->NeighborCount; NeighborIndex++)
          {
            PtrajClusteringSOMLearn(This, 
				    (PtrajSOMNode *) SOMNodes[WinnerIndex]->Neighbors[NeighborIndex],
				    &NeighborLearningRate,PointIndex);
          }
      }
    /* Decelerate your motion: */
    if (IterationIndex%10 == 0)
      {
        WinnerLearningRate *= 0.95;
        NeighborLearningRate *= 0.95;
      }
  }
  /* Use our SOMNode data to generate clusters */
  PtrajClusteringConvertSOMNodesToClusters( (Clustering *) This,SOMNodes,ClusterCount);
  PtrajClusteringFreeSOMNodes( (Clustering *) This,SOMNodes,ClusterCount);
  safe_free(PointProcessed);
}

void PtrajClusteringFinalizeBayesianClusters(Clustering* This, BayesianCluster* Head)
{
  BayesianCluster* BCluster;
  BayesianCluster* BestHost;
  Cluster* NewCluster;
  int PointIndex;
  float BestProb;
  int ClusterIndex;
  
  /**/
  fprintf(stdout, "  Finalizing bayesian clusters...\n");
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
    PtrajClusterAddMember( (PtrajCluster *) BestHost->TrueCluster,PointIndex);
  }

  fprintf(stdout, "  Bayesian clusters finalized.\n");

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


void PtrajClusteringClusterBayesian(PtrajClustering* This, int ClusterCount)
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
  AttributeCount = ClusteringGetAttributeCount( (Clustering *) This);
  Head = GenerateBayesianSeedClusters( (Clustering *) This, ClusterCount);
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
            BCluster->Means[AttributeIndex] += 
	      (double) ClusteringGetAttributeValue( (Clustering *) This,PointIndex,AttributeIndex) * 
	      BCluster->Members[PointIndex];
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
            Value = (double) ClusteringGetAttributeValue( (Clustering *) This,PointIndex,AttributeIndex);
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
  PtrajClusteringFinalizeBayesianClusters( (Clustering *) This,Head);
  fprintf(stdout,"  Free bayesian scaffolding:\n");
  BayesianClusterListFree(Head);
  fprintf(stdout,"  Bayesian clustering complete.\n");

}

int* PtrajClusteringFindKmeansSeeds(PtrajClustering* This)
{
    int* SeedIndices;
    int SeedIndex;
    int Seeds;
    float BestDistance;
    int CandidateFrameIndex, FrameIndex;
    int FrameCount = This->PointCount;
    int SkipFlag;
    int CheckIndex;
    float NearestDistance;
    float Distance;
    int BestDistanceIndex;
    /**/
    Seeds = PtrajClusteringGetDesiredClusterCount(This);
    SeedIndices = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int) * Seeds);
    SeedIndices[0]=1; /* Arbitrary first-choice */
    
    BestDistance = 0;
    for (FrameIndex = 0; FrameIndex < FrameCount; FrameIndex++) 
    {
        for (CandidateFrameIndex = FrameIndex; CandidateFrameIndex < FrameCount; CandidateFrameIndex++)
        {
            Distance = PtrajClusteringPointToPointDistance(This,FrameIndex,CandidateFrameIndex);
            if (Distance>BestDistance)
            {
                BestDistance = Distance;
                SeedIndices[0] = FrameIndex;
                SeedIndices[1] = CandidateFrameIndex;
            }
        }
    }
	
    for (SeedIndex=2; SeedIndex<Seeds; SeedIndex++)
    /*
    for (SeedIndex=1; SeedIndex<Seeds; SeedIndex++)
    */
    {
        BestDistance = 0;
        for (CandidateFrameIndex=0; CandidateFrameIndex<FrameCount; CandidateFrameIndex++)
        {
            /* Make sure This candidate isn't already a seed */
            SkipFlag=0;
            for (CheckIndex=0; CheckIndex<SeedIndex; CheckIndex++)
            {
                if (SeedIndices[CheckIndex]==CandidateFrameIndex)
                {
                    SkipFlag=1;
                    break; /* added by jianyin */
                }
            }
            if (SkipFlag)
            {
                continue;
            }
            /* Get the closest distance from This candidate to a current seed */
            NearestDistance = 0;
            for (CheckIndex=0; CheckIndex<SeedIndex; CheckIndex++)
            {
                Distance = PtrajClusteringPointToPointDistance(This,SeedIndices[CheckIndex],CandidateFrameIndex);
                if (Distance<NearestDistance || NearestDistance==0)
                {
                    NearestDistance = Distance;
                }
            }
            /* Is This the best so far? */
            if (NearestDistance > BestDistance)
            {
                BestDistance = NearestDistance;
                BestDistanceIndex = CandidateFrameIndex;
            }
        }
        SeedIndices[SeedIndex] = BestDistanceIndex;
    }
    return SeedIndices;
}

/* In this function, points are assigned to closest cluster based on the distance to each cluster's centroid.
   However, the cluster centroid is not updated by adding members. The updating will be done after all the points are assigned. 
   If FinishedPoints is NULL, all points will added as in a decoy 
   OldClusterMap will record the cluster each point belongs to. 
*/
void PtrajClusteringClusterDecoy(PtrajClustering* This, int* SeedPoints, int DesiredClusterCount, int iteration) 
{
    int PointIndex, PointCount;
    int ClusterIndex, IterationIndex;
    int verbose = This->action->darg2;
    int C2Ccentroid = This->action->darg4;
    ClusterNode* Node;
    PtrajCluster* ClosestCluster;
    PtrajCluster* NewCluster;
    float Distance;
    float Closest;
    int *OldClusterMap, *FinishedPoints; 
    int ClosestClusterIndex, i;
	int Changed = 0;
        
    PointCount = This->PointCount;
   	FinishedPoints = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
   	OldClusterMap = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
   	for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
       	OldClusterMap[PointIndex] = -1;
       	FinishedPoints[PointIndex] = 0;
    }
    i = 0;
    for (PointIndex=0; PointIndex<DesiredClusterCount; PointIndex++)
    {
    	if (!SeedPoints) break; /* Use Decoy structure. */
        NewCluster = (PtrajCluster *) ClusteringGetNewCluster( (Clustering *) This);
        PtrajClusterAddMember(NewCluster,SeedPoints[PointIndex]);
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        FinishedPoints[SeedPoints[PointIndex]]=1; 
        if (verbose)
        	fprintf(stdout, "Put frame %i in cluster%i.\n", SeedPoints[PointIndex], i++);
    }
   
    for (IterationIndex = 0; IterationIndex < iteration; IterationIndex++)
    {
        if (verbose)
        {
        	 fprintf(stdout, "Round %d\n", IterationIndex);
        }
        Changed = 0;
       	if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
   	    {
            fprintf(stdout, "OldClusterMap: ");
	    	for (PointIndex=0; PointIndex<PointCount; PointIndex++)
		    {
            	fprintf(stdout, "%3d", OldClusterMap[PointIndex]);
            }
            fprintf(stdout, "\n");
        }	
    	for (PointIndex=0; PointIndex<PointCount; PointIndex++)
	    {	
    	    if (IterationIndex==0 && FinishedPoints[PointIndex])
	        {
            	continue;
        	}
            if (IterationIndex > 0) 
            {
            	for (Node = This->Head, ClusterIndex = 0; Node; Node = Node->pNext, ClusterIndex++)
            	{
        	        if (ClusterIsMember(Node->Cluster,PointIndex))
    	            {
		                PtrajClusterRemoveMember( (PtrajCluster *) Node->Cluster,PointIndex);
				        FinishedPoints[PointIndex] = 0;
                    	OldClusterMap[PointIndex] = ClusterIndex;
                    }
                	
            	}
            }
    	    Closest = -1;
	        for (Node = This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
        	{
    	    	Distance = ClusteringDistanceToCentroid( (Clustering *) This,PointIndex,Node->Cluster);
	            /*printf(" point %d, cluster %d, distance %f; \n", PointIndex, ClusterIndex, Distance);*/
	            if (Closest < 0 || Distance<Closest)
            	{
        	        Closest = Distance;
    	            ClosestCluster = (PtrajCluster *) Node->Cluster;
                    ClosestClusterIndex = ClusterIndex;
	            }
        	}
            /*printf("\n");*/
    	    PtrajClusterAddMember(ClosestCluster,PointIndex);
	        FinishedPoints[PointIndex] = 1;
            if (ClosestClusterIndex != OldClusterMap[PointIndex]) 
            {
            	Changed++;
                if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) fprintf(stdout, "Remove Frame %d from cluster %d, but add to cluster %d.\n", PointIndex, OldClusterMap[PointIndex], ClosestClusterIndex);
            }
        	if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
    	    {
				fprintf(stdout, "Put Frame %d in cluster %d.\n", PointIndex, PtrajOutputIntermediateClusterStatus(ClosestCluster, stdout, 1)); 
            	PtrajOutputIntermediateClusterStatus(ClosestCluster, stdout, 0); 
        	}
	    }
        
    	for (Node = This->Head,ClusterIndex=0; ClusterIndex<DesiredClusterCount; Node = Node->pNext, ClusterIndex++)
	    {	
			if (C2Ccentroid > 0) 
	        {
        		if (ClusterCountPointsInCluster(Node->Cluster) != 0 )
        		{
			ClusterFindBestP2C( (PtrajCluster *) Node->Cluster);
    	    		ClusterAlignToBestRep( (PtrajCluster *) Node->Cluster);
	            	ClusterFindCentroid(Node->Cluster);
        		}
    	    }
	    }
        if (Changed == 0)
        {
        	fprintf(stdout, "  Means round %2d: No change.  Skip the rest of iterations.\n\n", IterationIndex);
            break;
        } else {
        	fprintf(stdout, "  Means round %2d: %5d points changing cluster assignment...\n", 
			IterationIndex, Changed);
        }
    }
    ClusteringRemoveEmptyCluster( (Clustering *) This);
cleanup:
    safe_free(OldClusterMap);
    safe_free(FinishedPoints);
}

/* mode == 1, sequentially modify each points.
   mode == 2, randomly pick point for modification.
*/
void PtrajClusteringClusterMeans(PtrajClustering* This, int* SeedPoints, int DesiredClusterCount, int interation, int mode)
{
    int PointIndex;
    int PointCount, UnprocessedPointCount;
    int* FinishedPoints;
    PtrajCluster* NewCluster;
    float Distance;
    float Closest;
    PtrajCluster* ClosestCluster;
    ClusterNode* Node;
    PtrajCluster* OldCluster;
    int Changed;
    int IterationIndex;
    int ClusterIndex, ProcessIndex;
    int C2Ccentroid = This->action->darg4;
    int verbose = This->action->darg2;
    int OldBestRep;
    int NewBestRep;
    int OldClusterIndex, ClosestClusterIndex;
    int Yanked;
    /**/
    
    PointCount = This->PointCount;
    FinishedPoints = (int*)SafeMalloc(__FILE__, __LINE__, sizeof(int)*PointCount);
    memset(FinishedPoints,0,sizeof(int) * PointCount);
    /*  randomly picking */
    srand(time(NULL) + clock());

    /* Add the "seed-clusters" */
    int i = 0;
    for (PointIndex=0; PointIndex<DesiredClusterCount; PointIndex++)
    {
    	if (!SeedPoints) break; /* Use Decoy structure. */
        NewCluster = (PtrajCluster *) ClusteringGetNewCluster( (Clustering *) This);
        /*PtrajClusterSetBestRep(NewCluster, SeedPoints[PointIndex]);*/
        PtrajClusterAddMember(NewCluster,SeedPoints[PointIndex]);
        ClusterFindCentroid( (Cluster *) NewCluster);
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        FinishedPoints[SeedPoints[PointIndex]]=1; 
        if (verbose)
        	fprintf(stdout, "Put frame %i in cluster%i.\n", SeedPoints[PointIndex], i++);
    }
    UnprocessedPointCount = PointCount - DesiredClusterCount;
	OldClusterIndex = -1;
    /* Assign points in iteration passes.  If a point looked like it belonged to cluster A
       at first, but then we added many other points and altered our cluster shapes, it's
       possible that we'll want to reassign it to cluster B. */
    for (IterationIndex = 0; IterationIndex < interation; IterationIndex++)
    {
      /* Add each point to an existing cluster, and re-compute centroid*/
      
      if (verbose)
      {
      	 fprintf(stdout, "Round %d\n", IterationIndex);
      }
      /* Process the training points in random order, so that we're not sensitive to the 
         data ordering */
      if (IterationIndex != 0) {
          memset(FinishedPoints,0,sizeof(int) * PointCount);
          UnprocessedPointCount = PointCount;
      }
      Changed = 0;
      for (ProcessIndex = 0; ProcessIndex < PointCount; ProcessIndex++)
      /*for (PointIndex=0; PointIndex<PointCount; PointIndex++)*/
      {
          if (mode == 1)
            PointIndex = ProcessIndex;
          else if (mode == 2)
            PointIndex = ChooseNextPoint(FinishedPoints,This->PointCount,UnprocessedPointCount-ProcessIndex);
          
          if (PointIndex < 0) continue;
          if (IterationIndex==0 && mode == 1 && FinishedPoints[PointIndex])
          {
              continue;
          }
          Yanked = 1;
          if (IterationIndex > 0)
          {
            /* Yank this point out of its cluster, recompute the centroid. Yanked is 1 unless the cluster has only one point. */
            for (Node = This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
            {
              if (ClusterIsMember(Node->Cluster,PointIndex))
              {
  	        	/* If this point is alone in its cluster, it's in the right place already! */
	        	if (ClusterCountPointsInCluster(Node->Cluster)==1)
				{
				  Yanked = 0;
                  continue;
				}
	        OldBestRep = PtrajClusterGetBestRep( (PtrajCluster *) Node->Cluster);
                PtrajClusterRemoveMember( (PtrajCluster *) Node->Cluster,PointIndex);
                OldClusterIndex = ClusterIndex;
		NewBestRep = ClusterFindBestP2C( (PtrajCluster *) Node->Cluster);
                ClusterFindCentroid(Node->Cluster);
	            if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
      		    {
          		  	/*fprintf(stdout, "Remove Frame %d from cluster %d.\n", PointIndex, PtrajOutputIntermediateClusterStatus(Node->Cluster, stdout, 1));*/ 
          		  	fprintf(stdout, "Remove Frame %d from cluster %d.\n", PointIndex, ClusterIndex); 
              		PtrajOutputIntermediateClusterStatus( (PtrajCluster *) Node->Cluster, stdout, 0); 
          		}
			    if (C2Ccentroid > 0) 
    	        {
	        	   if (OldBestRep != NewBestRep)
       			   {
	       			  ClusterAlignToBestRep( (PtrajCluster *) Node->Cluster);
       		       }
	   		  	   ClusterFindCentroid(Node->Cluster);
   		        }
              }
            }
          }
          if (!Yanked) continue;
          Closest = -1;
          for (Node = This->Head,ClusterIndex=0; Node; Node = Node->pNext,ClusterIndex++)
          {
              Distance = ClusteringDistanceToCentroid( (Clustering *) This,PointIndex,Node->Cluster);
              if (Closest < 0 || Distance<Closest)
              {
                  Closest = Distance;
                  ClosestCluster = (PtrajCluster *) Node->Cluster;
                  ClosestClusterIndex = ClusterIndex;
              }
          }
          OldBestRep = PtrajClusterGetBestRep(ClosestCluster);
          PtrajClusterAddMember(ClosestCluster,PointIndex);
		  NewBestRep = ClusterFindBestP2C(ClosestCluster);
          ClusterFindCentroid( (Cluster *) ClosestCluster);
          if (ClosestClusterIndex != OldClusterIndex) { 
              Changed++;
	          if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
      		  {
              	fprintf(stdout, "Remove Frame %d from cluster %d, but add to cluster %d.\n", PointIndex, OldClusterIndex, ClosestClusterIndex);
          	  }	
          }    
		  if (C2Ccentroid > 0) 
          {
	          if (OldBestRep != NewBestRep)
       		  {
	       		ClusterAlignToBestRep(ClosestCluster);
       		  }
                        ClusterFindCentroid( (Cluster *) ClosestCluster);
   		  }
          
          if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
          {
          	  fprintf(stdout, "Put Frame %d in cluster %d.\n", PointIndex, PtrajOutputIntermediateClusterStatus(ClosestCluster, stdout, 1)); 
              PtrajOutputIntermediateClusterStatus(ClosestCluster, stdout, 0); 
          }
       }
       if (Changed == 0)
       {
            fprintf(stdout, "  Means round %2d: No change.  Skipping the rest of the iterations!\n\n", IterationIndex);
            break;
       } else {
            fprintf(stdout, "  Means round %2d: %5d points changed cluster assignment...\n", 
		    IterationIndex, Changed);
       }
    }
    
    /* Remove cluster which did not contain any points */  
    ClusteringRemoveEmptyCluster( (Clustering *) This);
    for (Node = This->Head; Node; Node = Node->pNext)
    {
	    if (C2Ccentroid == 0) /* Align the cluster if using bestrep*/
    	{
	  ClusterFindBestP2C( (PtrajCluster *) Node->Cluster);
	  ClusterAlignToBestRep( (PtrajCluster *) Node->Cluster);
        }
    }

cleanup:
    safe_free(FinishedPoints);
}

/* Random clustering */
void PtrajClusteringClusterByChance(PtrajClustering* This, int DesiredClusterCount)
{
    int PointIndex;
    int ClusterIndex;
	float RandValue;
    int RandInt;
    Cluster* NewCluster;
    ClusterNode* Node;
    int C2Ccentroid = This->action->darg4;
    int verbose = This->action->darg2;
    
	Randomize();
    for (ClusterIndex=0; ClusterIndex<DesiredClusterCount; ClusterIndex++)
    {
        NewCluster = ClusteringGetNewCluster( (Clustering *) This);
        ClusteringAddCluster( (Clustering *) This,NewCluster);
    }
    for (PointIndex=0; PointIndex<This->PointCount; PointIndex++)
    {
        RandValue = rand()/(double)RAND_MAX;
        RandInt = RandValue * DesiredClusterCount;
        for (Node = This->Head,ClusterIndex = 0; ClusterIndex < RandInt; ClusterIndex++) {
        	Node = Node->pNext; 
        }
        PtrajClusterAddMember( (PtrajCluster *) Node->Cluster, PointIndex);
		if (C2Ccentroid > 0) 
        {
	    	ClusterAlignToBestRep( (PtrajCluster *) Node->Cluster);
        }
        ClusterFindCentroid(Node->Cluster);
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER)
        {
      		printf("point %d goes to cluster %d\n", PointIndex, ClusterIndex);        
        }
    }

}


void PtrajClusteringInitEdge(PtrajClustering* This)
{
	int PointCount, PointIndex;
    int C2Ccentroid = This->action->darg4;
    PtrajCluster* NewCluster;

    PointCount = This->PointCount;
    for (PointIndex=0; PointIndex<PointCount; PointIndex++)
    {
        NewCluster = PtrajClusterNew(This);
        /*PtrajClusterSetBestRep(NewCluster, PointIndex);   It will be done in PtrajClusterAddMember() for new cluster. */
        PtrajClusterAddMember(NewCluster,PointIndex);
        if (C2Ccentroid)
        {	
        	ClusterFindCentroid( (Cluster *) NewCluster);
        }
        ClusteringAddCluster( (Clustering *) This, (Cluster *) NewCluster);
        NewCluster->entry = (ClosestCluster *)SafeMalloc(__FILE__, __LINE__, sizeof(ClosestCluster));
        ((ClosestCluster *)NewCluster->entry)->Distance = -1;
        ((ClosestCluster *)NewCluster->entry)->ClosestNode = NULL;
    }
}

void PtrajClusteringOldEdgeLinkage(PtrajClustering* This, int DesiredClusterCount, float Epsilon)
{
    int PointIndex;
    int PointCount = This->PointCount;
    int NodeIndexA;
    int NodeIndexB;
    int PointA;
    int PointB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    float Distance;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    ClusterNode* ClosestNode = NULL;
    float C2CDistance;
    int verbose = This->action->darg2;
    int C2Ccentroid = This->action->darg4;
	ClosestCluster *temp;	
    /* */
    PtrajClusteringInitEdge(This);
	/* Initialize the ClosestNode and closest Distance for each Cluster. */
    for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
    {
        ClosestLength = -1;
        for (NodeA = This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	if (NodeA == NodeB)
       	    	continue;
           	C2CDistance = GetSymElement(This->PairwiseDistances,NodeIndexA,NodeIndexB);
            if (ClosestLength < 0 || C2CDistance < ClosestLength)
            {
            	ClosestLength = C2CDistance;
                ClosestNode = NodeA; 
            }
        }
        temp = ((PtrajCluster *)NodeB->Cluster)->entry;
        temp->ClosestNode = ClosestNode;
        temp->Distance = ClosestLength;
        /*NodeB->Cluster->entry->ClosestNode = ClosestNode;
        NodeB->Cluster->entry->Distance = ClosestLength;*/
    }
    
    /* Keep joining until you're finished */
    while (1)
    {
        PtrajClusteringOutputStatus(This);
        /* If there's 1 cluster left, then we're definitely done: */
		if (This->ClusterCount < 2)
	  		return;
        
        /* Find two closest Node to merge */
        MergeNodeA = NULL;
        MergeNodeB = NULL;
        ClosestLength = -1;
        for (NodeA = This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	temp = ((PtrajCluster *)NodeA->Cluster)->entry;
            Distance = temp->Distance;
            if (ClosestLength < 0 || Distance < ClosestLength)
            {
            	ClosestLength = Distance;
                MergeNodeA = NodeA;
                MergeNodeB = temp->ClosestNode;
            }
        }
        
        if (ClosestLength<0) /* We found nothing to merge, in this case */
            return;
        /* If the pair are too far apart, stop! */
        if (!DesiredClusterCount && ClosestLength>Epsilon)
            return;
        /* Fix up distances for all the nodes */
		ClosestLength = -1;
        for (NodeA = This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	temp = ((PtrajCluster *)NodeA->Cluster)->entry;
            if (NodeA == MergeNodeA) 
            	continue;
            else if (NodeA == MergeNodeB)
            	continue;
            else if (temp->ClosestNode == MergeNodeA || temp->ClosestNode == MergeNodeB)
            {
                if (ClosestLength < 0 || ClosestLength > temp->Distance)
                {
                	ClosestLength = temp->Distance;
                    ClosestNode = NodeA;
                }
            	temp->ClosestNode = MergeNodeB;
            } 
            else /* Those who are not involved with merge, we still need to check their distances to the newly merged Node, i.e. MergeNodeB*/
            {
            	Distance = -1;
                for (PointA = 0; PointA < PointCount; PointA++)
                {
			    	if (!ClusterIsMember(MergeNodeA->Cluster, PointA) && 
                    	!ClusterIsMember(MergeNodeB->Cluster, PointA)) /* For every point in MergeNodeA or MergeNodeB, check mininum distance to NodeA */
			        {
			        	continue;
        			}
                	Distance = ClusteringMinDistanceToCluster( (Clustering *) This,PointA,NodeA->Cluster);
	                if ((ClosestLength < 0) || (Distance < ClosestLength))
	                {
	                	ClosestLength = Distance;
	                    ClosestNode = NodeA;
	                }
                }
            } 
        }
        temp = ((PtrajCluster*)MergeNodeB->Cluster)->entry;
        temp->Distance = ClosestLength;
        temp->ClosestNode = ClosestNode;  
        /*((ClosestCluster *)((PtrajCluster*)MergeNodeB->Cluster)->entry)->Distance = ClosestLength;
        ((ClosestCluster *)((PtrajCluster*)MergeNodeB->Cluster)->entry)->ClosestNode = ClosestNode;  */
        
        /* Merge!  Remove cluster A, and corresponding linkage node A */
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	        fprintf(stdout, "merge cluster %d into cluster %d\n",
			PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
			PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
	       	fprintf(stdout, "Mergee is: \n");
		   	PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	       	fprintf(stdout, "Merger is: \n");
		   	PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0);   
        }
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
       	if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
        	fprintf(stdout, "After Merger is: \n");
	   		PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
	        fprintf(stdout, "Cluster count is %i\n", This->ClusterCount);
        } 
        if (C2Ccentroid) 
        {	
        	ClusterFindCentroid(MergeNodeB->Cluster);
        }
        if (This->ClusterCount<=DesiredClusterCount)
        {
            return; /* We've merged enough, we're done */
        }
    }
}
void PtrajClusteringOldCompleteLinkage(PtrajClustering* This, int DesiredClusterCount, float Epsilon)
{
    int PointIndex;
    int PointCount = This->PointCount;
    int NodeIndexA;
    int NodeIndexB;
    int PointA;
    int PointB;
    ClusterNode* MergeNodeA;
    ClusterNode* MergeNodeB;
    int ClosestNodeIndexA;
    int ClosestNodeIndexB;
    float Distance;
    float ClosestLength;
    ClusterNode* NodeA;
    ClusterNode* NodeB;
    ClusterNode* ClosestNode = NULL;
    float C2CDistance, MergeNodeDistance;
    int verbose = This->action->darg2;
    int C2Ccentroid = This->action->darg4;
	ClosestCluster *temp;	
    /* */
    PtrajClusteringInitEdge(This);
	/* Initialize the ClosestNode and closest Distance for each Cluster. */
    for (NodeB = This->Head,NodeIndexB=0; NodeB; NodeB=NodeB->pNext,NodeIndexB++)
    {
        ClosestLength = -1;
        for (NodeA = This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	if (NodeA == NodeB)
   	        {
       	    	continue;
           	}
           	C2CDistance = GetSymElement(This->PairwiseDistances,NodeIndexA,NodeIndexB);
            if (ClosestLength < 0 || C2CDistance < ClosestLength)
            {
            	ClosestLength = C2CDistance;
                ClosestNode = NodeA; 
            }
        }
        temp = ((PtrajCluster *)NodeB->Cluster)->entry;
        temp->ClosestNode = ClosestNode;
        temp->Distance = ClosestLength;
        /*NodeB->Cluster->entry->ClosestNode = ClosestNode;
        NodeB->Cluster->entry->Distance = ClosestLength;*/
    }
    /* Keep joining until you're finished */
    while (1)
    {
        PtrajClusteringOutputStatus(This);
        /* If there's 1 cluster left, then we're definitely done: */
		if (This->ClusterCount < 2)
        {
	  		return;
		}
        
        /* Find two closest Node to merge */
        MergeNodeA = NULL;
        MergeNodeB = NULL;
        ClosestLength = -1;
        for (NodeA = This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	temp = ((PtrajCluster *)NodeA->Cluster)->entry;
            Distance = temp->Distance;
            /*printf("ClosestLength is %f, Epsilon is %f\n", ClosestLength, Epsilon);*/
            if (ClosestLength < 0 || Distance < ClosestLength)
            {
            	ClosestLength = Distance;
                MergeNodeA = NodeA;
                MergeNodeB = temp->ClosestNode;
            }
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
        /* Fix up distances for all the nodes */
		ClosestLength = -1;
        for (NodeA = This->Head,NodeIndexA=0; NodeA; NodeA=NodeA->pNext,NodeIndexA++)
        {
        	temp = ((PtrajCluster *)NodeA->Cluster)->entry;
            if (NodeA == MergeNodeA) 
            {
            	continue;
            }
            else if (NodeA == MergeNodeB)
            {
            	continue;
            }
            /* For any clusters' closest node is MergeNodeA, we need to find the closest node for them. 
               C2CDistance holds the largest distance of NodeA to each NodeB,
               MergeNodeDistance holds the larger of the distance of NodeA to MergeNodeA or MergeNodeB, i.e., the distance of NodeA to newly merged node.
               ClosestLength holds the closest length to newly merged node MergeNodeB */
            else if (temp->ClosestNode == MergeNodeA) 
            {
                MergeNodeDistance = temp->Distance; 
                temp->Distance = -1;
                for (NodeB=This->Head,NodeIndexB=0;NodeB;NodeB=NodeB->pNext,NodeIndexB++)
                {
                	if (NodeB == NodeA || NodeB == MergeNodeA) continue;
                    C2CDistance = -1;
                    for (PointA = 0; PointA < PointCount; PointA++)
                    {
                    	if (!ClusterIsMember(NodeA->Cluster, PointA)) continue;
                        Distance = ClusteringMaxDistanceToCluster( (Clustering *) This,PointA,NodeB->Cluster);
		                if (Distance > C2CDistance)	C2CDistance = Distance;
	        	    }
                    if (C2CDistance < temp->Distance || temp->Distance < 0) 
                    {
                    	temp->Distance = C2CDistance;
                        temp->ClosestNode = NodeB;
                    }
                	if ((NodeB == MergeNodeB) && (C2CDistance > MergeNodeDistance))
            	    {
        	        	MergeNodeDistance = C2CDistance;
	                }
                }
               	if (ClosestLength > MergeNodeDistance || ClosestLength < 0)
           	    {
       	        	ClosestLength = MergeNodeDistance;
   	                ClosestNode = NodeA;
                }
            } 
            /* For any clusters' closest node is MergeNodeB, we need to find the closest node for them. 
               C2CDistance holds the largest distance of NodeA to each NodeB,
               MergeNodeDistance holds the larger of the distance of NodeA to MergeNodeA or MergeNodeB, i.e., the distance of NodeA to newly merged node.
               ClosestLength holds the closest length to newly merged node MergeNodeB */
            else if (temp->ClosestNode == MergeNodeB) 
            {
                MergeNodeDistance = temp->Distance; 
                temp->Distance = -1;
                for (NodeB=This->Head,NodeIndexB=0;NodeB;NodeB=NodeB->pNext,NodeIndexB++)
                {
                	if (NodeB == NodeA || NodeB == MergeNodeB) continue;
                    C2CDistance = -1;
                    for (PointA = 0; PointA < PointCount; PointA++)
                    {
                    	if (!ClusterIsMember(NodeA->Cluster, PointA)) continue;
                        Distance = ClusteringMaxDistanceToCluster( (Clustering *) This,PointA,NodeB->Cluster);
		                if (Distance > C2CDistance)	C2CDistance = Distance;
	        	    }
                    if (C2CDistance < temp->Distance || temp->Distance < 0) 
                    {
                    	temp->Distance = C2CDistance;
                        temp->ClosestNode = NodeB;
                        if (NodeB == MergeNodeA) temp->ClosestNode = MergeNodeB;
                    }
                	if ((NodeB == MergeNodeA) && (C2CDistance > MergeNodeDistance))
            	    {
        	        	MergeNodeDistance = C2CDistance;
	                }
                }
               	if (ClosestLength > MergeNodeDistance || ClosestLength < 0)
           	    {
       	        	ClosestLength = MergeNodeDistance;
   	                ClosestNode = NodeA;
                }
            } 
            /* Those who are not involved with merge, we still need to check their distances to the newly merged Node, i.e. MergeNodeB*/
            else 
            {
            	C2CDistance = -1;
                for (PointA = 0; PointA < PointCount; PointA++)
                {
			    	if (!ClusterIsMember(MergeNodeA->Cluster, PointA) && 
                    	!ClusterIsMember(MergeNodeB->Cluster, PointA)) /* For every point in MergeNodeA or MergeNodeB, check maxium distance to NodeA */
			        {
			        	continue;
        			}
                	Distance = ClusteringMaxDistanceToCluster( (Clustering *) This,PointA,NodeA->Cluster);
	                if ((Distance > C2CDistance) || (C2CDistance < 0) )
	                {
	                	C2CDistance = Distance;
	                }
                }
                if ((C2CDistance < ClosestLength) || (ClosestLength < 0))
                {
                	ClosestLength = C2CDistance;
                    ClosestNode = NodeA;
                }
            } 
        }
        temp = ((PtrajCluster*)MergeNodeB->Cluster)->entry;
        temp->Distance = ClosestLength;
        temp->ClosestNode = ClosestNode;  
        /*((ClosestCluster *)((PtrajCluster*)MergeNodeB->Cluster)->entry)->Distance = ClosestLength;
        ((ClosestCluster *)((PtrajCluster*)MergeNodeB->Cluster)->entry)->ClosestNode = ClosestNode;  */
        
        /* Merge!  Remove cluster A, and corresponding linkage node A */
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	        fprintf(stdout, "merge cluster %d into cluster %d\n",
			PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 1),
			PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 1));
        }
        if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
	       	fprintf(stdout, "Mergee is: \n");
		   	PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeA->Cluster, stdout, 0); 
	       	fprintf(stdout, "Merger is: \n");
		   	PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0);   
        }
        ClusteringMergeClusters( (Clustering *) This,MergeNodeA,MergeNodeB);
       	if (verbose >= VERBOSE_INTERMEDIATE_CLUSTER) 
        {
        	fprintf(stdout, "After Merger is: \n");
	   		PtrajOutputIntermediateClusterStatus( (PtrajCluster *) MergeNodeB->Cluster, stdout, 0); 
	        fprintf(stdout, "Cluster count is %i\n", This->ClusterCount);
        } 
        if (C2Ccentroid) 
        {	
        	ClusterFindCentroid(MergeNodeB->Cluster);
        }
        if (This->ClusterCount<=DesiredClusterCount)
        {
            return; /* We've merged enough, we're done */
        }
    }
}


/* print Series --- For debugging */

ClusterNode* printClusterNode(PtrajClustering *This, int Cluster)
{
	int i;
    ClusterNode* temp;
    
    for (i=0,temp=This->Head; i<This->ClusterCount;i++,temp=temp->pNext)
    {
    	if (i == Cluster) 
        	return temp;
    }

    return NULL;
}

PtrajCluster* printCluster(PtrajClustering *This, int Cluster)
{
	int i;
    ClusterNode* temp;
    
    for (i=0,temp=This->Head; i<This->PointCount;i++,temp=temp->pNext)
    {
    	if (i == Cluster) 
        	return (PtrajCluster*)temp->Cluster;
    }

    return NULL;
}

ClusterNode* printClosestNode(PtrajClustering *This, int Cluster)
{
	int i;
    ClusterNode* temp;
    
    for (i=0,temp=This->Head; i<This->PointCount;i++,temp=temp->pNext)
    {
    	if (i == Cluster) 
        	return (ClusterNode*)((ClosestCluster*)((PtrajCluster*)temp->Cluster)->entry)->ClosestNode;
    }

    return NULL;
}

void printClusterNodes(PtrajClustering *This, char mode)
{
	int i;
    ClusterNode *temp;
    
    for (temp=This->Head,i=0;temp;temp=temp->pNext,i++)
    {
    	if(mode == 'n')
        	printf("cluster%i--%x\n", i, temp);
    	else if(mode == 'c')
        	printf("cluster%i--%x\n", i, (PtrajCluster*)temp->Cluster);
    	else if(mode == 'd')
        	printf("cluster%i--%f\n", i, ((ClosestCluster*)((PtrajCluster*)temp->Cluster)->entry)->Distance);
    	else if(mode == 'l')
        	printf("cluster%i--%x\n", i, (ClusterNode*)((ClosestCluster*)((PtrajCluster*)temp->Cluster)->entry)->ClosestNode);
    }	
}

/* Used in debugging, print out the Cluster/Node number and so on. */
void printClusters(PtrajClustering *This)
{
	int i,j;   
    ClusterNode *temp;
    
    printf("clusters--Node   \tCluste \tClosest    \tDistance\tBestRep   \tPoints\tMask\n");
    for (temp=This->Head,i=0;temp;temp=temp->pNext,i++)
    {
       	printf("cluster%i--%x\t", i, temp);
       	printf("%x\t", (PtrajCluster*)temp->Cluster);
       	if (((PtrajCluster*)temp->Cluster)->entry)
        {
        	printf("%x(%i)\t", (ClusterNode*) 
		       ((ClosestCluster*)((PtrajCluster*)temp->Cluster)->entry)->ClosestNode, 
		       PtrajOutputIntermediateClusterStatus( (PtrajCluster *) ((ClosestCluster*)((PtrajCluster*)temp->Cluster)->entry)->ClosestNode->Cluster, stdout, 1));
        printf("%f\t", ((ClosestCluster*)((PtrajCluster*)temp->Cluster)->entry)->Distance);
        } else {
	       	printf("n.a.          \t");
	       	printf("n.a.     \t");
       	}
       	printf("%i\t\t", ((PtrajCluster*)temp->Cluster)->BestRep);
        printf("%d\t", ClusterCountPointsInCluster((Cluster *) (PtrajCluster*)temp->Cluster));
        printf("("); 
        for (j=0; j<This->PointCount;)
        {
        	printf("%i", temp->Cluster->Mask[j++]);
            if (j % 5 == 0) printf(" ");  
        }
        printf(")\n");
        
    }	
}

/* Used for debugging, print out frame x, y, z of atom a, b, c*/
void printFrameXYZ(trajectoryInfo *trajInfo, int frame, int atom, int number)
{
	int i;
    for (i = atom; i < atom + number; i++)
    {
    	printf("%.3f \t%.3f \t%.3f\n", trajInfo->x[frame*trajInfo->atoms+i],trajInfo->y[frame*trajInfo->atoms+i],trajInfo->z[frame*trajInfo->atoms+i]);
    }
}


/* end of print Series */
