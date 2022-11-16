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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/cluster.h,v 10.0 2008/04/15 23:24:11 case Exp $
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

#define CLUSTER_OUTPUT_NONE -1
#define CLUSTER_OUTPUT_AVERAGE 1
#define CLUSTER_OUTPUT_REPRESENTATIVE 2
#define CLUSTER_OUTPUT_ALL 4

#define CLUSTER_FILEFORMAT_AVERAGE 0
#define CLUSTER_FILEFORMAT_REPRESENTATIVE 1
#define CLUSTER_FILEFORMAT_ALL 2
#define CLUSTER_FILEFORMAT_LEN 3

#define ClusterMergingFile "ClusterMerging.txt"
typedef struct Representative
{
    int MovedFlag;
    int PointIndex; /* The point in the original dataset that we correspond to */
    float* X;
    float* Y;
    float* Z;
} Representative;

typedef struct RepresentativeNode
{
    struct RepresentativeNode* pNext;
    struct RepresentativeNode* pPrev;
    struct Representative* Point;
} RepresentativeNode;

typedef struct PtrajCluster
{
  int* Mask; /* Mask[n] is true iff point n is in this cluster*/
  int PointCount; /* Size of the Mask array */
  int CentroidIndex; /* Index of our centroid point */
  int IntName;  /* Use an integer as its name. */
  float SSEWithin;     /* Sum of squared error */
  struct PtrajClustering* Owner;
  /* Function table: */
  VirtualFunction* VirtualFunctions;
  FloatVirtualFunction* FloatVirtualFunctions;
  float* CentroidX;
  float* CentroidY;
  float* CentroidZ;
  RepresentativeNode* Head;
  RepresentativeNode* Tail;
  /* These members are used when processing sieved clusters, for sanity: */
  float* AverageX;
  float* AverageY;
  float* AverageZ;
  float* RepX;
  float* RepY;
  float* RepZ;
  float MinRepFrameDistance;
  int MinRepFrameIndex;
  char* Name;
  /* Indicate the number of frame of the best representative, whose sum of distance to other members in the cluster is smallest.
     The cluster should be aligned with the best representative. 
     In the case of two-member cluster, the first point will be picked as the BestRep. */
  int BestRep;
  /* If the second pass is used for sieving, OldBestRep saves the BestRep numbering for the first pass */
  int OldBestRep;
     /* The Point2Cluster array stores the sum of the distances of the point to the other members in the cluster.
        If only one member is in the cluster, the value should be 0. 
        The values of the points in different cluster may not be comparable since the numbers of points in the cluster are different.   */
  float* P2C;
  /* entry can store any data structure.
     Currently it is store a structure of ClosestNode,
  */
  void* entry;
} PtrajCluster;


/* Here are all the available distance metrics for clustering */
#define DISTANCE_METRIC_RMSD 0
#define DISTANCE_METRIC_DME 1
#define DISTANCE_METRIC_MDS 2

typedef struct PtrajClustering
{
    SymmetricMatrix* PairwiseDistances; /* Element [X,Y] is the distance from point X to point Y */
    ClusterNode* Head; /*First cluster*/
    ClusterNode* Tail; /*Last cluster*/
    int PointCount;
    int ClusterCount; /* Length of the cluster list */
    float SSE;  /* sum of squared of error */
    float SST;  /* sum of squared of total */
    float Acuity;
    Cluster* TotalCluster;
    arrayType *attributeArray;
    arrayType *attributeArrayTorsion;
/* Function table: */
    VirtualFunction* VirtualFunctions;
    FloatVirtualFunction* FloatVirtualFunctions;
    IntVirtualFunction* IntVirtualFunctions;
    trajectoryInfo* trajInfo;
    actionInformation* action;
    time_t StartTime;
    time_t MatrixTime;
    time_t FirstPassTime;
    int DistanceMetric;
    /* For sieving purposes, this array tracks the frames which will be processed
     during the first pass.  We randomly choose 1 in every n points (where n is the 
     sieve size).  We don't choose every nth point, because we might be confused 
     by harmonic motion of our structure (imagine sampling the peak of each wave in 
     a sine curve, then clustering the points; they'd all be the same, and wounldn't
     represent the full range of motion well at all) */
     int* SieveFirstPass;
     int SievedIndex;
     
} PtrajClustering;

typedef struct PtrajCobwebCluster
{
  int PointIndex;
  float CU;
  float* Means;
  float* Stddevs;
  float* SumX;
  float* SumX2;
  int numbers;
  float* Cos;
  float* Sin;
  float* Cos2;
  float* Sin2;
} PtrajCobwebCluster;

typedef struct PtrajCobwebTreeNode
{
    PtrajCobwebCluster* Cluster;
    int TotalFrameCount; /* How many points are in this cluster and all its children */
    struct CobwebTreeNode* Next;
    struct CobwebTreeNode* Previous;
    struct CobwebTreeNode* Parent;
    struct CobwebTreeNode* FirstChild;
    struct CobwebTreeNode* LastChild;
} PtrajCobwebTreeNode;

typedef struct PtrajBayesianCluster
{
  int AttributeCount;
  double* Means; 
  double* Stddevs;
  double Probability;
  double* Members;
  struct BayesianCluster* Next;
  struct BayesianCluster* Prev;
  Cluster* TrueCluster; /* Convert to this */
  float* SumX;
  float* SumX2;
  int numbers;
  float* Cos;
  float* Sin;
  float* Cos2;
  float* Sin2;
} PtrajBayesianCluster;

struct PtrajSOMNodeStruct
{
  float* Attributes;
  int NeighborCount;
  struct SOMNodeStruct** Neighbors;
 /* The following variables are used for calculating the actual count, means, stddevs for the cluster. They will get the actual values when we assign each points to the closest SOMNodes. */
  float* Means;
  float* Stddevs;
  int PointCount;
  float* Cos;
  float* Sin;
  float* Cos2;
  float* Sin2;
  /* The MeanCos and MeanSin will be used for calculate the actual mean angles. They will get actual value after learning process. */
  float* MeanCos;
  float* MeanSin;
};
typedef struct PtrajSOMNodeStruct PtrajSOMNode;

/* Representative functions */
Representative* RepresentativeNew(int PointIndex,int AtomCount);
void RepresentativeFree(Representative* This);

/* Cluster functions */
PtrajCluster* PtrajClusterNew(PtrajClustering* Owner);
void PtrajClusterAddMember(PtrajCluster* This, int Point);
void PtrajClusterRemoveMember(PtrajCluster* This, int Point);
void PtrajClusterFree(PtrajCluster* This);
void ExpandClusterSize(PtrajCluster* This, int Multiplier, int NewFrameCount);
void PtrajClusterFindCentroid(PtrajCluster* This);
void ClusterAddMemberP2C(PtrajCluster* This, int Point);
void ClusterSubtractMemberP2C(PtrajCluster* This, int Point);
int ClusterFindBestP2C(PtrajCluster* This);
int PtrajOutputIntermediateClusterStatus(PtrajCluster* This, FILE* OutFile, int IndexOnly); 

/* Clustering functions */
PtrajClustering* PtrajClusteringNew(SymmetricMatrix* PairwiseDistances, trajectoryInfo* trajInfo, 
				    actionInformation* action);
void PtrajClusteringFree(PtrajClustering* This);
void PtrajClusteringClusterLinkage(PtrajClustering* This,int DesiredClusterCount,float Epsilon);
Cluster* PtrajClusteringGetNewCluster(PtrajClustering* This);
void PtrajClusteringMergeClusters(PtrajClustering* This,ClusterNode* NodeA,ClusterNode* NodeB);
int PtrajClusteringGetFrameCount(PtrajClustering* This);
int PtrajClusteringGetAtomCount(PtrajClustering* This);
float PtrajClusteringDistanceToCentroid(PtrajClustering* This,int PointIndex, PtrajCluster* Cluster);
float PtrajClusteringDistanceToCluster(PtrajClustering* This,int PointIndex, PtrajCluster* Cluster);
ClusterNode* PtrajClusteringSplitCluster(PtrajClustering* This,Cluster* ParentCluster, int PointA, int PointB);
float PtrajClusteringPointToPointDistance(PtrajClustering* This,int A,int B);
float PtrajClusteringCentroidToCentroidDistance(PtrajClustering* This,ClusterNode* NodeA,ClusterNode* NodeB);
float PtrajClusteringClusterToClusterDistance(PtrajClustering* This,ClusterNode* NodeA,ClusterNode* NodeB);
int PtrajClusteringGetDesiredClusterCount(PtrajClustering* This);
float PtrajClusteringGetDesiredEpsilon(PtrajClustering* This);
int* PtrajClusteringFindKmeansSeeds(PtrajClustering* This);
void PtrajClusteringAddFrameToBestCluster(PtrajClustering* This,int FrameIndex, float* x, float* y, float* z);

void PtrajClusteringClusterCentripetal(PtrajClustering* This,int DesiredClusterCount, float Epsilon, int Representatives);
void PtrajClusteringClusterCentripetalComplete(PtrajClustering* This,int DesiredClusterCount, float Epsilon, int Representatives);
void PtrajClusterRepositionRepresentatives(PtrajCluster* This);
int PtrajClusterGetRepresentativeCount(PtrajCluster* This);
void PtrajClusterPruneRepresentatives(PtrajCluster* This,int MaxRepresentatives);
float PtrajClusteringRepToRepDistance(PtrajClustering* This,Representative* RepA, Representative* RepB);
float PtrajClusteringMinimumRepresentativeDistance(PtrajClustering* This,PtrajCluster* ClusterA,PtrajCluster* ClusterB);
float PtrajClusteringMaximumRepresentativeDistance(PtrajClustering* This,PtrajCluster* ClusterA,PtrajCluster* ClusterB);

void PtrajClusteringOutputClusterInfo(PtrajClustering* This,actionInformation* action);
/*void PtrajClusteringClusterBayesian(PtrajClustering* This, int DesiredClusterCount, float Epsilon); // like AutoClass */
PtrajCobwebCluster* PtrajClusteringNewCobwebCluster(PtrajClustering* Clustering, int PointIndex);
void PtrajClusteringCobwebAddNewFrame(PtrajClustering* This,CobwebTreeNode* Root,CobwebTreeNode* FrameNode, int level); 
void PtrajCobwebClusterComputeMeans(PtrajClustering* This, PtrajCobwebTreeNode* Node);
void PtrajClusteringCobwebClusterFree(PtrajClustering* Clustering, PtrajCobwebCluster* This);
void PtrajClusteringCobwebCopyMeans(PtrajClustering* This, PtrajCobwebCluster* Source, PtrajCobwebCluster* Destination);
float PtrajClusteringCobwebCategoryUtility(PtrajClustering* This, CobwebTreeNode* Parent);
void PtrajClusteringClusterCobweb(Clustering* This, int DesiredClusterCount, int mode);
void PtrajCobwebFlattenTree(Clustering* This, CobwebTreeNode* Root, int DesiredClusterCount);

int PtrajClusteringGetAttributeCount(PtrajClustering* This);
float PtrajClusteringGetAttributeMin(PtrajClustering* This, int AttributeIndex);
float PtrajClusteringGetAttributeMax(PtrajClustering* This, int AttributeIndex);
float PtrajClusteringGetAttributeValue(PtrajClustering* This, int PointIndex, int AttributeIndex);
int PtrajClusteringGetAttributeArrayCount(PtrajClustering* This);
float PtrajClusteringGetAttributeArrayMin(PtrajClustering* This, int AttributeIndex);
float PtrajClusteringGetAttributeArrayMax(PtrajClustering* This, int AttributeIndex);
float PtrajClusteringGetAttributeArrayValue(PtrajClustering* This, int PointIndex, int AttributeIndex);
float PtrajClusteringGetAttributeArraySin(PtrajClustering* This, int PointIndex, int AttributeIndex);
float PtrajClusteringGetAttributeArrayCos(PtrajClustering* This, int PointIndex, int AttributeIndex);

float PtrajClusteringSOMGetDistance(PtrajClustering* This, PtrajSOMNode* SOM, int PointIndex);
void PtrajClusteringSOMLearn(PtrajClustering* This, PtrajSOMNode* SOM, float* LearningRate, int PointIndex);
void PtrajClusteringClusterSOM(PtrajClustering* This, int ClusterCount, int map);
extern PtrajSOMNode** PtrajClusteringInitSOMNodes(PtrajClustering* This, int ClusterCount, int map);
void PtrajClusteringFreeSOMNodes(Clustering* This,PtrajSOMNode** SOMNodes,int ClusterCount);
float PtrajClusteringSOMGetDistance(PtrajClustering* This, PtrajSOMNode* SOM, int PointIndex);
PtrajSOMNode** PtrajClusteringGetSOMMap(PtrajClustering* This);

void PtrajClusteringClusterMeans(PtrajClustering* This,int* SeedPoints, int DesiredClusterCount, int iteration, int mode);
extern void PtrajClusteringClusterHierarchical(PtrajClustering* This,int DesiredClusterCount,float Epsilon);

void PtrajClusteringEdgeLinkage(PtrajClustering* This, int DesiredClusterCount, float Epsilon);
void PtrajClusteringCompleteLinkage(PtrajClustering* This, int DesiredClusterCount, float Epsilon);
void PtrajClusteringAverageLinkage(PtrajClustering* This, int DesiredClusterCount, float Epsilon);

void PtrajClusteringClusterBayesian(PtrajClustering* This, int DesiredClusterCount); /* like AutoClass */
void PtrajClusteringFinalizeBayesianClusters(Clustering* This, BayesianCluster* Head);


trajectoryInfo* PtrajClusteringReadDecoy(PtrajClustering* This, char* filename, int ClusterCount);
void PtrajClusteringClusterDecoy(PtrajClustering* This,int* SeedPoints, int DesiredClusterCount, int iteration);

/* For Sieve */
void PtrajClusteringChooseSievePoints(PtrajClustering* This, int FullPointCount);
extern void ExpandClusterListSizes(PtrajClustering* This,int Multiplier, int NewFrameCount);


/* For Output */
void PtrajClusteringOutputHeaderToFile(PtrajClustering* This, FILE* File);
void PtrajClusteringOutputStatus(PtrajClustering* This);
extern void OutputClusteringStats(PtrajClustering* This, FILE* File);
void PtrajClusteringPrintTransformMap(PtrajClustering* This);

/* print Series --- For debugging */
ClusterNode* printClusterNode(PtrajClustering *This, int Cluster);
PtrajCluster* printCluster(PtrajClustering *This, int Cluster);
PtrajCluster* printCluster(PtrajClustering *This, int Cluster);
void printClusterNodes(PtrajClustering *This, char mode);
void printClusters(PtrajClustering *This);
void printFrameXYZ(trajectoryInfo *trajInfo, int frame, int atom, int number);

  /*
   *  Vtables 
   */
static VirtualFunction PtrajClusterVirtualFunctions[FCOUNT_CLUSTER] = 
  { (void *) PtrajClusterFree,
    (void *) PtrajClusterFindCentroid
  };

  /*
   *  Remember: Don't point to an un-prefixed method like ClusteringClusterLinkage,
   *  because it'll just call itself again, leading to stack overflow, and you will cry. 
   */
static VirtualFunction PtrajClusteringVirtualFunctions[FCOUNT_CLUSTERING] = 
  { (void *) PtrajClusteringFree, 
    (void *) PtrajClusteringMergeClusters,
    (void *) PtrajClusteringSplitCluster,
    (void *) PtrajClusteringGetNewCluster,
    (void *) PtrajClusteringFindKmeansSeeds,
    (void *) PtrajClusteringClusterLinkage,
    (void *) PtrajClusteringClusterMeans,
    (void *) BASEClusteringClusterCobweb,
    (void *) PtrajClusteringNewCobwebCluster, 
    (void *) PtrajClusteringCobwebAddNewFrame,
    (void *) PtrajCobwebClusterComputeMeans,
    (void *) PtrajClusteringCobwebClusterFree, 
    (void *) PtrajClusteringCobwebCopyMeans,
    (void *) PtrajClusteringSOMLearn, 
    (void *) PtrajClusteringOutputStatus
  }; 

static FloatVirtualFunction PtrajClusteringFloatVirtualFunctions[FFCOUNT_CLUSTERING] =
  { PtrajClusteringPointToPointDistance, 
    PtrajClusteringCentroidToCentroidDistance, 
    PtrajClusteringDistanceToCluster,
    PtrajClusteringCobwebCategoryUtility, 
    PtrajClusteringGetAttributeValue,
    PtrajClusteringGetAttributeMin,
    PtrajClusteringGetAttributeMax, 
    PtrajClusteringSOMGetDistance
  };
  
static IntVirtualFunction PtrajClusteringIntVirtualFunctions[FICOUNT_CLUSTERING] =
  { PtrajClusteringGetAttributeCount
  };

/* Other functions */
void OutputClusterList(SymmetricMatrix* PairwiseDistance,trajectoryInfo* trajInfo,actionInformation* action);
extern SymmetricMatrix* ComputePairwiseDistance(PtrajClustering* This, trajectoryInfo* trajInfo,  actionInformation* action, char *CacheDistanceFile);

extern float PtrajGetDistance(PtrajClustering*, int, int, float*, float*, float*, 
			      float*, float*, float*);


static char* CLUSTER_ALGORITHM_NAMES[] = 
  { "Hierarchical",
    "Linkage",
    "Means(1)",
    "Centripetal",
    "Cobweb(1)",
    "Bayesian",
    "SOM",
    "DBI",
    "EdgeLink",
    "CompleteLink",
    "ANOVA",
    "ByChance",
    "Centripetal_Complete",
    "ReadFromClusterMerging", 
    "ClusterComparing", 
    "AverageLink",
    "Decoy",
    "Means2",
    "Means3",
    "ReadTxt",
    "Cobweb2",
    "AverageLinkFingerPrint",  
    "EdgeLinkFingerPrint",  
    "CompleteLinkFingerPrint",  
    "CentripetalFingerPrint",  
    "CentripetalCompleteFingerPrint",  
    "LinkFingerPrint",
    "SOM2"
  };

static char* DISTANCE_METRIC_NAMES[] = 
  { "RMS",
    "DME",
    "MDS"
  };
 
#define VERBOSE_INTERMEDIATE_CLUSTER 6
#define VERBOSE_ALIGN_FRAME 7

