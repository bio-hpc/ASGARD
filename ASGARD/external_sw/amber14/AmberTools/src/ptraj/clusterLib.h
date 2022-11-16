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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/clusterLib.h,v 10.1 2010/02/22 21:08:34 sbrozell Exp $
 *
 *  Revision: $Revision: 10.1 $
 *  Date: $Date: 2010/02/22 21:08:34 $
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


#define abs(X) ((X)>0 ? (X) : -(X))
#define max(X,Y) (((X)>(Y)) ? (X) : (Y))
#define min(X,Y) (((X)>(Y)) ? (Y) : (X))


  /*
   *  Clustering algorithms group like points.  The Cluster and Clustering objects
   *  in this library are used in general-purpose clustering algorithms.  A Cluster is
   *  a set of points.  A Clustering is essentially just a list of clusters.  The Clusters
   *  in a Clustering partition the set of all points.  Most points have many attributes,
   *  which are generally real-valued.  (For instance, points in the plane have two
   *  attributes, their x and y coordinates)
   *
   *  Before you can perform clustering on any kind of data, you must create extension 
   *  Clustering and Cluster objects, and alter any functions that operate differently on
   *  your particular data.  For an example of how to do this, see the ClusterTest 
   *  executable; it performs clustering on the canonical example, points in the plane.
   *
   *  In order to cluster a group of points, you must generally first create a 
   *  SymmetricMatrix (see matrix.h)  storing the pairwise distances between any 
   *  two points, then use the pairwise distance matrix to create an instance of 
   *  the Clustering class.  (Not all algorithms require a pairwise distance matrix 
   *  - see the algorithm comments for the few exceptions)
   *
   *  Brief overview of clustering algorithms:
   *  Linkage - Start with singleton clusters, and merge close-by clusters into larger ones.
   *    The clusters with shortest centroid-to-centroid distance are merged.
   *  Hierarchical - Start with one big cluster, and split it into smaller ones
   *  K-Means - Start with n seeds, and grow clusters around each seed
   *  Centripetal - (Based on CURE) Merge smaller clusters to larger ones, using centralized
   *    representatives to merge more sensibly than in pure Linkage
   *  Cobweb - Build a tree of clusters and sub-clusters, minimizing standard deviation
   *    within a cluster
   *  Bayesian - Expectation Maximization (EM) clustering.  Start with randomly-defined 
   *    normal distributions for attributes of each cluster, and iteratively improve 
   *    the clusters.
   *  SOM - Self-organizing maps.  Start with a collection random map points.  In each iteration,
   *    run through all the data points.  The "winning" (most similar) map point draws closer to 
   *    the data point.  This represents "training".  Over the course of many iterations, 
   *    each map point will gravitate toward the centroid of its neighbors.
   *  EdgeLink - Just link linkage, but join clusters when edge-points are close
   **/

#ifndef CLUSTER_H
#define CLUSTER_H

  /* 
   *  Use CLUSTER_API to mark exports from a .dll in Windows.
   *  (On unix, all functions are visible from the library)
   */

#ifdef WIN32
#ifdef CLUSTER_EXPORTS
#define CLUSTER_API __declspec(dllexport)
#else
#define CLUSTER_API __declspec(dllimport)
#endif
#else
#define CLUSTER_API
#endif

#define CLUSTER_HIERARCHICAL 		 0
#define CLUSTER_LINKAGE 		 1
#define CLUSTER_MEANS 			 2
#define CLUSTER_CENTRIPETAL 		 3
#define CLUSTER_COBWEB 			 4
#define CLUSTER_BAYESIAN 		 5
#define CLUSTER_SOM 		         6
/* Special case: Run DBI on an already-finished sieved clustering */
#define CLUSTER_DBI 		         7 
#define CLUSTER_EDGELINK 		 8 
#define CLUSTER_COMPLETELINK		 9 
#define CLUSTER_ANOVA 			10 
#define CLUSTER_BYCHANCE	 	11 
#define CLUSTER_CENTRIPETAL_COMPLETE 	12
#define CLUSTER_READMERGE 	       	13
#define CLUSTER_COMPARE        		14
#define CLUSTER_AVERAGELINK 		15
#define CLUSTER_DECOY 			16
#define CLUSTER_MEANS2 			17
#define CLUSTER_MEANS3 		        18
#define CLUSTER_READTXT		       	19
#define CLUSTER_COBWEB2			20
#define CLUSTER_AVERAGELINK_FP		21
#define CLUSTER_EDGELINK_FP		22
#define CLUSTER_COMPLETELINK_FP		23
#define CLUSTER_CENTRIPETAL_FP		24
#define CLUSTER_CC_FP			25
#define CLUSTER_LINKAGE_FP		26
#define CLUSTER_SOM2    		27

  /*
   *  MEANS1 is same as MEANS, SEQUENTIALLY add or adjust points and update centroid after each point adjustment. 
   *  MEANS2 RANDOMLY add or adjust points and update centroid after each point adjustment. 
   *  MEANS3 SEQUENTIALLY add or remove points but update the centroid ONLY after each round. 
   */


  /*
   *  Symmetric matrices are stingy on space: We only store the upper triangle.
   *  Indexing into the matrix
   *  data block is: Row*Size - Row**2/2 + Row + Column
   */

typedef struct CLUSTER_API SymmetricMatrix {
  int Size;
  float* Data;
} SymmetricMatrix;

  /*
   * 2D Matrix with float elements 
   */
typedef struct CLUSTER_API Matrix {
    int Size;
    float* Data;
} Matrix;

CLUSTER_API Matrix* AllocateMatrix(int Size);
CLUSTER_API void FreeMatrix(Matrix* pMatrix);

CLUSTER_API int TriangularizeMatrix(Matrix* pMatrix);
CLUSTER_API float Determinant(Matrix* pMatrix);

/*
PrintMatrix prints out a matrix to stdout, for debugging purposes
*/
CLUSTER_API void PrintMatrix(Matrix* pMatrix);

/*
PrintMatrixToFile prints out a matrix to a file, for debugging purposes
*/
CLUSTER_API void PrintMatrixToFile(Matrix* pMatrix, char* FileName);

/*int UnitTestMatrices();*/
/*int LapackSolveMatrix(Matrix* pMatrix);*/

#define MatrixElement(matrix,row,column) matrix->Data[((row)*(matrix->Size)) + (column)]

CLUSTER_API SymmetricMatrix* AllocateSymmetricMatrix(int Size);
CLUSTER_API SymmetricMatrix* SymmetricMatrixCopy(SymmetricMatrix* Original);
CLUSTER_API void FreeSymmetricMatrix(SymmetricMatrix* pMatrix);
CLUSTER_API float GetSymElement(SymmetricMatrix* pMatrix, int Row, int Column);
CLUSTER_API void SetSymElement(SymmetricMatrix* pMatrix, int Row, int Column, float Value);

typedef struct CLUSTER_API DoubleSymmetricMatrix {
    int Size;
    double* Data;
} DoubleSymmetricMatrix;

CLUSTER_API DoubleSymmetricMatrix* AllocateDoubleSymmetricMatrix(int Size);
CLUSTER_API DoubleSymmetricMatrix* SymmetricDoubleMatrixCopy(SymmetricMatrix* Original);
CLUSTER_API void FreeDoubleSymmetricMatrix(DoubleSymmetricMatrix* pMatrix);
CLUSTER_API double GetDoubleSymElement(DoubleSymmetricMatrix* pMatrix, int Row, int Column);
CLUSTER_API void SetDoubleSymElement(DoubleSymmetricMatrix* pMatrix, int Row, int Column, double Value);

/*float* ClapackEigenvalues(Matrix* pMatrix);*/
/*float ClapackDeterminant(Matrix* pMatrix);*/

/*
This macro only works if row <= column!
*/
#define SymmetricMatrixElement(matrix,row,column) matrix->Data[(row)*(matrix)->Size + ((row) - (row)*(row))/2 + ((column)-(row))]

CLUSTER_API void PrintSymMatrix(SymmetricMatrix* Matrix);
CLUSTER_API void WriteSymMatrix(char* FilePath,SymmetricMatrix* Matrix, int TextMode);
CLUSTER_API int ReadSymMatrix(char* FilePath,SymmetricMatrix* Matrix);



  /* 
   *  A brief explanation of what's going on with the virtual functions:
   *  We want to be able to use the clustering code as a foundation upon which
   *  to build various custom clusterings.  (We can cluster atoms, cluster frames,
   *  cluster points in the plane, cluster anything).  The natural way to do this
   *  is with a class and subclasses.  The C equivalent is to create a table of
   *  virtual functions, with one entry for each overridable method of the 
   *  clustering object.  "Subclasses" (or rather, structs that contain the same
   *  members, possibly with some appended) should have their own virtual-function 
   *  table, and point the VirtualFunctions member of each instance at the
   *  table. 
   * 
   *  An UGLY HACK is that we cannot return floats in a VirtualFunction.
   *  Same problem for ints. So, we have separate tables for the separate return
   *  types.
   */

typedef void* (*VirtualFunction)();
typedef float (*FloatVirtualFunction)();
typedef int (*IntVirtualFunction)();

  /*
   *  These constants index into the VirtualFunction tables of a Cluster 
   */
#define F_CLUSTER_FREE 0
#define F_CLUSTER_FIND_CENTROID 1
#define FCOUNT_CLUSTER 2

  /*
   *  Indices into the VirtualFunctions of Clustering 
   */
#define F_CLUSTERING_FREE 0
#define F_CLUSTERING_MERGE_CLUSTERS 1
#define F_CLUSTERING_SPLIT_CLUSTER 2
#define F_CLUSTERING_GET_NEW_CLUSTER 3
#define F_CLUSTERING_FIND_KMEAN_SEEDS 4
#define F_CLUSTERING_CLUSTER_LINKAGE 5
#define F_CLUSTERING_CLUSTER_MEANS 6
#define F_CLUSTERING_CLUSTER_COBWEB 7
#define F_CLUSTERING_NEW_COBWEB_CLUSTER 8
#define F_CLUSTERING_COBWEB_ADD_FRAME 9
#define F_CLUSTERING_COBWEB_COMPUTE_MEANS 10
#define F_CLUSTERING_COBWEB_CLUSTER_FREE 11
#define F_CLUSTERING_COBWEB_COPY_MEANS 12
#define F_CLUSTERING_SOM_LEARN 13
#define F_CLUSTERING_OUTPUT_STATUS 14
#define FCOUNT_CLUSTERING 15

  /*
   *  Indices into the FloatVirtualFunctions of Clustering 
   */
#define F_CLUSTERING_POINT_TO_POINT 0
#define F_CLUSTERING_CENTROID_TO_CENTROID 1
#define F_CLUSTERING_DISTANCE_TO_CENTROID 2
#define F_CLUSTERING_COBWEB_CATEGORY_UTILITY 3
#define F_CLUSTERING_GET_ATTRIBUTE_VALUE 4
#define F_CLUSTERING_GET_ATTRIBUTE_MIN 5
#define F_CLUSTERING_GET_ATTRIBUTE_MAX 6
#define F_CLUSTERING_SOM_GET_DISTANCE 7
#define FFCOUNT_CLUSTERING 8

  /*
   *  Indices into the IntVirtualFunctions of Clustering 
   */
#define F_CLUSTERING_GET_ATTRIBUTE_COUNT 0
#define FICOUNT_CLUSTERING 1

typedef struct CLUSTER_API Cluster
{
  int* Mask; /* Mask[n] is true iff point n is in this cluster*/
  int PointCount; /* Size of the Mask array */
  /* Index of our centroid point.  (Most subclasses ignore CentroidIndex; they
     construct a centroid point instead of choosing an existing point) */
  int CentroidIndex; 
  int IntName;  /* Use an integer as its name. */
  /* Sum of squared error within cluster, sum( dist( x(i), avg(x(i)) )^2 ) */
  float SSEWithin;
  /* Pointer to our clustering */
  struct Clustering* Owner;
  /* Function tables: */
  VirtualFunction* VirtualFunctions;
  FloatVirtualFunction* FloatVirtualFunctions;
} Cluster;

  /*
   *  A ClusterNode is a member of a doubly-linked list; it holds one cluster. 
   */
typedef struct CLUSTER_API ClusterNode
{
    struct Cluster* Cluster;
    struct ClusterNode *pNext;
    struct ClusterNode *pPrev;
} ClusterNode;


typedef struct CLUSTER_API Clustering
{
    /* Element [X,Y] in PairwiseDistances is the distance from point X to point Y */
    SymmetricMatrix* PairwiseDistances; 
    ClusterNode* Head; /*First cluster*/
    ClusterNode* Tail; /*Last cluster*/
    int PointCount; /* Number of points being clustered */
    int ClusterCount; /* Length of the cluster list */
    float SSE;
    float SST;
    float Acuity;
    Cluster* TotalCluster;
    arrayType *attributeArray;
    arrayType *attributeArrayTorsion;
    /* Function tables: */
    VirtualFunction* VirtualFunctions;
    FloatVirtualFunction* FloatVirtualFunctions;
    IntVirtualFunction* IntVirtualFunctions;
} Clustering;

  /*
   *  A LinkageNode is used internally for Linkage clustering 
   */
typedef struct LinkageNode
{
    struct LinkageNode* pNext;
    struct LinkageNode* pPrev;
    ClusterNode* pCluster;
    float ClosestDistance;
} LinkageNode;

typedef struct ClosestCluster
{
    ClusterNode* ClosestNode;	
	float Distance;
} ClosestCluster;


#define COBWEB_NONLEAF -1

  /*
   *  The CobwebCluster and CobwebTreeNode structs are used only in COBWEB
   *  clustering; a CobwebCluster may be involved in several trees at any
   *  one time. 
   */
typedef struct CobwebCluster
{
    int PointIndex; /* COBWEB_NON_LEAF for non-leaf nodes */
    float CU;
    float* Means;
    float* Stddevs;
    float* SumX;
    float* SumX2;
    int numbers;
} CobwebCluster;

typedef struct CobwebTreeNode
{
    CobwebCluster* Cluster;
    int TotalFrameCount; /* How many points are in this cluster and all its children */
    struct CobwebTreeNode* Next;
    struct CobwebTreeNode* Previous;
    struct CobwebTreeNode* Parent;
    struct CobwebTreeNode* FirstChild;
    struct CobwebTreeNode* LastChild;
} CobwebTreeNode;

  /*
   *  A simple linked-list of CobwebTreeNodes.  This list doesn't
   *  relate to the tree structure that the CobwebTreeNodes themselves
   *  track.  The CTNList struct lets us make a list of nodes from all over the 
   *  tree, without disrupting the tree structure.  
   */
typedef struct CTNList
{
    CobwebTreeNode* Node;
    struct CTNList* Next;
    struct CTNList* Previous;
} CTNList;

  /* 
   *  A bayesian cluster has a distribution (mean and standard deviation) for
   *  each attribute in our space.  A "typical" point for the cluster is one
   *  whose first attribute is near Means[0], et cetera.  The probabilities 
   *  for a collection of BayesianCluster instances sum to 1.  The Members
   *  array has one entry for each point, and each member gives the probability 
   *  that a particular point is a member of THIS bayesian cluster.
   */
typedef struct BayesianCluster
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
} BayesianCluster;

/*
typedef struct SOMNode
{
  float* Attributes;
} SOMNode;
*/
struct SOMNodeStruct
{
  float* Attributes;
  int NeighborCount;
  struct SOMNodeStruct** Neighbors;
 /* The following variables are used for calculating the actual count, means, stddevs for the cluster. They will get the actual values when we assign each points to the closest SOMNodes. */
  float* Means;
  float* Stddevs;
  int PointCount;
};
typedef struct SOMNodeStruct SOMNode;

/*****************************************************************
 *  Functions pertaining to the cluster struct follow - these are
 *  the "members" of the Cluster "class". 
 *****************************************************************/

  /*
   *  Returns the number of points available for clustering (NOT the 
   *  number of members in this cluster) 
   */
CLUSTER_API int ClusterGetPointCount(Cluster* This);

  /*
   *  Returns true if the specified point is in this cluster 
   */
CLUSTER_API int ClusterIsMember(Cluster* This, int PointIndex);

  /*
   *  Adds the specified point to this cluster 
   */
CLUSTER_API void ClusterAddMember(Cluster* This, int PointIndex);

  /*
   *  Add all points to this cluster 
   */
CLUSTER_API void ClusterAddEverything(Cluster* This);

  /*
   *  Get the index of the point being used as this cluster's centroid.
   *  (Most subclasses build a new centroid point, and so don't use this method) 
   */
CLUSTER_API int ClusterGetCentroidIndex(Cluster* This);

  /*
   *  Remove a point from this cluster 
   */
CLUSTER_API void ClusterRemoveMember(Cluster* This, int PointIndex);

  /*
   *  Allocate a new Cluster object; free with ClusterFree 
   */
CLUSTER_API Cluster* ClusterNew(Clustering* Owner);

  /*
   *  Free a Cluster object 
   */
CLUSTER_API void ClusterFree(Cluster* This);

  /*
   *  Find the centroid of this cluster.  After altering the cluster's 
   *  membership (adding and removing points), you must call ClusterFindCentroid
   *  to keep the centroid up-to-date 
   */
CLUSTER_API void ClusterFindCentroid(Cluster* This);

  /*
   *  Return the number of points in this cluster 
   */
CLUSTER_API int ClusterCountPointsInCluster(Cluster* This);

  /*
   *  Base-class implementations of overridable functions: 
   */
CLUSTER_API void BASEClusterFree(Cluster* This);
CLUSTER_API void BASEClusterFindCentroid(Cluster* This);


/*****************************************************************
 *  Functions pertaining to the clustering struct follow - these are
 *  the "members" of the Clustering "class". 
 *****************************************************************/

  /*
   *  Allocate a new cluster object for this clustering.  (Afterwards, you still
   *  must call ClusteringAddClustering if you want the Cluster to our children) 
   */
CLUSTER_API Cluster* ClusteringGetNewCluster(Clustering* This);

  /*
   *  Returns the PairwiseDistances (set at initialization time) for a Clustering 
   */
CLUSTER_API SymmetricMatrix* ClusteringGetPairwiseDistances(Clustering* This);

  /*
   *  Returns the number of points being clustered 
   */
CLUSTER_API int ClusteringGetPointCount(Clustering* This);

  /*
   *  Allocates a new Clustering object.  Free with ClusteringFree. 
   */
CLUSTER_API Clustering* ClusteringNew(SymmetricMatrix* PairwiseDistances);

  /*
   *  Frees a Clustering object.  Also frees its clusters! 
   */
CLUSTER_API void ClusteringFree(Clustering* This);

  /*
   *  Adds a Cluster (probably created by ClusteringGetNewCluster) to our children 
   */
CLUSTER_API ClusterNode* ClusteringAddCluster(Clustering* This, Cluster* NewCluster);

  /*
   *  Removes a cluster from our list of children.  (e.g. if merging two 
   *  clusters into one, get rid of an old one).  Does *not* free the 
   *  old cluster, just removes its node. 
   */
CLUSTER_API void ClusteringRemoveCluster(Clustering* This, Cluster* OldCluster);

  /*
   *  Remove empty clusters. 
   */
CLUSTER_API void ClusteringRemoveEmptyCluster(Clustering* This);

  /*
   *  Swap the order of two cluster next to each, used in the bubble sorting 
   */
CLUSTER_API void ClusteringSwapCluster(Clustering* This, ClusterNode* Node);

  /*
   *  Bubble sorting of the Cluster, largest to smallest 
   */
CLUSTER_API void ClusteringSortCluster(Clustering* This);

  /*
   *  Iterate over each cluster, and call FindCentroid on each one.  Useful if you
   *  have just preformed major restructuring on them. 
   */
CLUSTER_API void ClusteringFindAllCentroids(Clustering* This);

  /*
   *  Read a clustering from a text-file.  Re-reading the clustering won't re-read
   *  the point definitions; we only get the membership info.  Useful for comparing
   *  the results of clustering algorithms. 
   */
CLUSTER_API void ClusteringReadFromDisk(Clustering* This,char* SourceFilePath);

  /*
   *  Return an array whose nth entry is the cluster index that the nth point has
   *  been assigned to. Used for comparing results.
   */
CLUSTER_API int* ClusteringGetClusterNumberArray(Clustering* This);

  /*
   *  Helper function for comparing results: Assuming random distribution (that is,
   *  using only our cluster sizes), what are the odds that two points fall in
   *  the same cluster? 
   */
CLUSTER_API float ClusteringGetOddsSameCluster(Clustering* This);

  /*
   *  Measure the distance between this cluster and another.  Our distance metric
   *  is defined as follows: Say two clusterings "agree" on a pair of points if they 
   *  both assign them to the same cluster, or if they both assign them to different
   *  clusters.  We iterate over all pairs of points, and compute the number of times
   *  our two clusterings disagree.  Then, we take the ratio of this sum with the 
   *  "expected" number of disagreements.  Two random clusterings will have a distance
   *  near 1.0; a clustering has distance 0.0 from itself; it is possible (though not
   *  likely) to have distances >1 
   */
CLUSTER_API float ClusteringMeasureDistance(Clustering* This,Clustering* Other);

  /*
   *  Print out a clustering to stdout, for debugging purposes.  Prints lines of the 
   *  form "...X..X", where an X in the xth position on the yth row means that point 
   *  x is a member of cluster y 
   */
CLUSTER_API void ClusteringDebugPrintClusterList(Clustering* This);

  /*
   *  Creates a file at the specified path, and calls ClusteringOutputClusterListToFile 
   */
CLUSTER_API void ClusteringOutputClusterList(Clustering* This,char* FilePath);

  /*
   *  Calls ClusteringOutputClusterListToFile, *appending* to the specified file 
   */
CLUSTER_API void ClusteringAppendClusterList(Clustering* This,char* FilePath);

  /*
   *  Prints a cluster list to a file.  If IncludeSeparators is true, include extra
   *  lines (which we can't parse later) for better readability 
   */
CLUSTER_API void ClusteringOutputClusterListToFile(Clustering* This, FILE* OutFile,int IncludeSeparators);
CLUSTER_API void ClusteringOutputClusterCondenseToFile(Clustering* This, FILE* OutFile);
CLUSTER_API void ClusteringOutputClusterConsensusToFile(Clustering* This, FILE* OutFile,int IncludeSeparators);

  /*
   *  Perform HIERARCHICAL clustering 
   */
CLUSTER_API void ClusteringClusterHierarchical(Clustering* This,int DesiredClusterCount,float Epsilon);

  /*
   *  Perform K-means clustering 
   */
CLUSTER_API void ClusteringClusterMeans(Clustering* This,int* SeedPoints, int DesiredClusterCount, int iteration, int mode);

  /*
   *  Perform clustering by chance 
   */
CLUSTER_API void ClusteringClusterByChance(Clustering* This, int DesiredClusterCount);

  /*
   *  Helper for Linkage clustering: Find the closest cluster to NodeA (where "closest"
   *  means minimum point-to-point distance between any two cluster members) 
   */
CLUSTER_API void ClusteringLinkageSeekClosest(Clustering* This,ClusterNode* NodeA,LinkageNode* Node);

  /*
   *  Single-linkage clustering - start with each point in its own singleton-cluster,
   *  and merge the closest clusters until the number of clusters is low enough (or
   *  cluster sizes are at the max allowed limit).  To make the whole process efficient,
   *  we do NOT compute all the cluster-to-cluster distances at every step.  Rather,
   *  we track the current nearest-neighbor for each cluster.  When you merge clusters,
   *  every other cluster just has to decide whether the new merged cluster will 
   *  displace its old nearest-neighbor.
   */
CLUSTER_API void ClusteringClusterLinkage(Clustering* This, int DesiredClusterCount,float Epsilon);

  /*
   *  Variant of Linkage clustering: Instead of merging the cluster-pair with shortest
   *  centroid-to-centroid distance, merge the pair of clusters with minimal point-to-point
   *  distance.
   */
CLUSTER_API void ClusteringClusterEdgeLink(Clustering* This, int DesiredClusterCount,float Epsilon);

  /*
   *  The next three linkage clustering methods are optimized. 
   *  A copy of PairwiseDistances is used to keep track the distance between two clusters. 
   *  An array of cookies is used to register the position of a cluster in the copy of PairwiseDistances.
   *  AverageLinkage, the distance between two clusters is the average of all distances between 
   *  individial points of the two clusters.  
   *  CompleteLinkage, the distance between two clusters is the longest inter-cluster point-to-point distance.
   *  EdgeLinkage, the distance between two clusters is the shortest inter-cluster point-to-point distance.
   */
CLUSTER_API void ClusteringClusterAverageLinkage(Clustering* This, int DesiredClusterCount, float Epsilon);
CLUSTER_API void ClusteringClusterCompleteLinkage(Clustering* This, int DesiredClusterCount, float Epsilon);
CLUSTER_API void ClusteringClusterEdgeLinkage(Clustering* This, int DesiredClusterCount, float Epsilon);

  /*
   *  Merge two clusters into one large cluster containing all the points from the 
   *  original.  (We remove NodeA from our list, and put the new cluster into NodeB) 
   */
CLUSTER_API void ClusteringMergeClusters(Clustering* This,ClusterNode* NodeA,ClusterNode* NodeB);

  /*
   *  Split a cluster.  Each of the two specified points becomes a "seed" for one of
   *  the two children. Each point in ParentCluster is assigned to the child-cluster
   *  whose seed it is nearest to.  
   */
CLUSTER_API void ClusteringSplitCluster(Clustering* This,Cluster* ParentCluster, int PointA, int PointB);

  /*
   *  Find the most eccentric cluster.  The eccentricity of a cluster is defined as
   *  the maximal point-to-point distance for points in the cluster.  This function 
   *  also tracks the distance (Epsilon), and the indices of the two eccentric 
   *  points. 
   */
CLUSTER_API void ClusteringFindEccentricCluster(Clustering* This, float* Epsilon, Cluster** EccentricCluster, int* PointA, int* PointB);

  /*
   *  Return the head of our cluster-list 
   */
CLUSTER_API ClusterNode* ClusteringGetHead(Clustering* This);

  /*
   *  Find some seed-points for K-means clustering.  Take the first point as an 
   *  arbitrary first choice.  Then, at each iteration, add the point whose total 
   *  distance from our set of seeds is as large as possible.  
   */
CLUSTER_API int* ClusteringFindKmeansSeeds(Clustering* This, int Seeds);

  /*
   *  Return the distance between two clusters' points 
   */
CLUSTER_API float ClusteringPointToPointDistance(Clustering* This, int A, int B);

  /*
   *  Return the distance between two clusters' centroids 
   */
CLUSTER_API float ClusteringCentroidToCentroidDistance(Clustering* This, ClusterNode* NodeA,ClusterNode* NodeB);

  /*
   *  Measure the distances between two or more clusterings.  Invokable via the 
   *  ClusterTest program!  Assumes that each command-line argument is a filename,
   *  of a file created by ClusteringOutputClusterListToFile 
   */
void MeasureDistances(int argc, char* argv[], int UseChiSquared);

  /*
   *  Helper function for COBWEB.
   *  Try to split a Node. Promote all its children to replace the parent node.
   *  Return the CU.
   *  Will restore the original cobweb tree.
   *  If execute == 1, do not restore. 
   */
CLUSTER_API float CobwebTrySplit(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* Node, int execute);

  /*
   *  Helper function. 
   *  Try merge two children of Root, report the new CU for Root. 
   *  Return the CU.
   *  Will free the newly merged node and restore the original state of the cobweb tree. 
   *  If execute == 1, merge.
   */
CLUSTER_API float CobwebTryMerge(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* Child1, CobwebTreeNode* Child2, int execute);

  /*
   *  Helper function. 
   *  Merge two children of Root, report the new CU for Root. 
   *  Will try to split one of the two nodes. If CU increases, try to split both.
   */
CLUSTER_API float CobwebMerge(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* Child1, CobwebTreeNode* Child2);

  /*
   *  Helper function in CobwebCoalscece
   *  Reorganize the CobwebTree by splitting or merging, so that the number of the children 
   *  of Root is equal to DesiredClusterCount.
   */
CLUSTER_API void CobwebTreeReorganize(Clustering* This, CobwebTreeNode* Root, int DesiredClusterCount);

  /* 
   *  Perform COBWEB clustering. 
   *  A tree of clusters is built up incrementally.  At each stage, we try to maximize the 
   *  Category Utility (see below) of the tree.  Each node in the cobweb tree is either a 
   *  leaf node (corresponding to ONE point), or a parent node (corresponding to a cluster 
   *  of all the points in its descendant leaf nodes).  We consider MERGING and SPLITTING
   *  nodes, as we go, in order to re-model the tree.  After every frame has been added, 
   *  we COALESCE the tree into ordinary "flat" clusters, to get something we can compare
   *  with other algorithm output.
   **/
CLUSTER_API void ClusteringClusterCobweb(Clustering* This, int DesiredClusterCount);

  /*
   *  Create a new Cobweb cluster object.  You should override this,
   *  to allocate the right number of means and stddevs.  
   */
CLUSTER_API CobwebCluster* ClusteringNewCobwebCluster(Clustering* This, int FrameIndex);

  /*
   *  Helper function for COBWEB clustering: Add the specified point to the
   *  cluster tree.  Called recursively with a lower-in-the-tree root. 
   */
CLUSTER_API void ClusteringCobwebAddNewFrame(Clustering* This,CobwebTreeNode* Root,CobwebTreeNode* Frame,int level);

  /*
   *  Helper function for COBWEB clustering: Compute means and standard-deviations
   *  for the specified tree node.  This function should be called on a non-leaf 
   *  node whenever there is a change in its chldren.  For permanent changes to the
   *  tree, you should *also* call this method for all of the node's parents. 
   */
CLUSTER_API void CobwebClusterComputeMeans(Clustering* This, CobwebTreeNode* Node);

  /*
   *  Compute category-utility for a cobweb tree node.  The category-utility of a
   *  node is a measure of its "usefulness".  The CU is computed by iterating over
   *  all attributes, and over all child nodes.  For each child node, take 
   *  take (1/a - 1/b), where a is the standard deviation of the child, and b is 
   *  the standard deviation of the parent.  The ACUITY of cobweb clustering (defined
   *  in the clustering sub-class) is a floor on the standard deviations we use in
   *  this computation.  (A child with only one point in it has a standard deviation
   *  of zero, but we use the acuity - say, 0.1 - to avoid dividing by zero).  
   */
CLUSTER_API float ClusteringCobwebCategoryUtility(Clustering* This, CobwebTreeNode* Parent);

  /*
   *  Free memory from a cobweb cluster.  Frees JUST this node, not its children. 
   */
CLUSTER_API void ClusteringCobwebClusterFree(Clustering* This, CobwebCluster* DeadGuy);

  /*
   *  Copy the means and standard deviations from one cobweb cluster to another.
   *  Useful for saving old means, when doing temporary rearrangements of the tree 
   */
CLUSTER_API void CobwebClusterCopyMeans(Clustering* This, CobwebCluster* Source, CobwebCluster* Destination);

CLUSTER_API CTNList* CobwebTreeNodeListChildren(CobwebTreeNode*);

CLUSTER_API CTNList* CTNListAppend(CTNList*, CobwebTreeNode*);

  /*
   *  Returns the number of attributes of each point 
   */
CLUSTER_API int ClusteringGetAttributeCount(Clustering* This);

  /*
   *  Returns the minimum value of the specified attribute, across *all* points.
   *  Used in bayesian clustering, when we are assigning each cluster a mean that is
   *  somewhere between the min and max attribute values 
   */
CLUSTER_API float ClusteringGetAttributeMin(Clustering* This, int AttributeIndex);

  /*
   *  Returns the maximum value of the specified attribute, across *all* points. 
   */
CLUSTER_API float ClusteringGetAttributeMax(Clustering* This, int AttributeIndex);

  /*
   *  Returns the value of a specified attribute for a specified point 
   */
CLUSTER_API float ClusteringGetAttributeValue(Clustering* This, int PointIndex, int AttributeIndex);

  /*
   *  Performs bayesian clustering 
   */
CLUSTER_API void ClusteringClusterBayesian(Clustering* This, int ClusterCount);

  /*
   *  The Davies-Bouldin Index (DBI) is a measure of clustering merit; the smaller the DBI,
   *  the better.  The DBI is defined as the average, for all clusters X, of fred, where 
   *  fred(X) = max, across other clusters Y, of (Cx + Cy)/dXY
   *  ...here Cx is the average distance from points in X to the centroid, similarly Cy, and dXY is 
   *  the distance between cluster centroids
   */
CLUSTER_API float ClusteringComputeDBI(Clustering* This, FILE* OutputFile);

  /*
   *  ClusteringBuildTotalCluster is the auxilary function for computing psF. 
   *  It will build a cluster consisting of all points, and calculate the 
   *  TotalDistanceToCentroid stored in This->SST. 
   */
CLUSTER_API void ClusteringBuildTotalCluster(Clustering* This);
CLUSTER_API float ClusteringComputePseudoF(Clustering* This);

  /*
   *  Compare two Clusterings, write out file 
   */
CLUSTER_API void ClusteringCompare(Clustering** Clusterings, char* FilePath); 

CLUSTER_API float ClusteringSOMGetDistance(Clustering* This, SOMNode* SOM, int PointIndex);
CLUSTER_API void ClusteringSOMLearn(Clustering* This, SOMNode* SOM, float* LearningRate, int PointIndex);
CLUSTER_API void ClusteringClusterSOM(Clustering* This, int ClusterCount, int map);
CLUSTER_API void ClusteringOutputStatus(Clustering* This);

/*****************************************************************/
/* The following Clustering methods are overridable, using the vtables: */
/* Compute the distance from the specified point to the centroid of the specified
   cluster.  Must override this, if you construct a centroid! */
CLUSTER_API float ClusteringDistanceToCentroid(Clustering* This,int PointIndex, Cluster* Cluster);
CLUSTER_API float ClusteringMinDistanceToCluster(Clustering* This,int PointIndex, Cluster* Cluster);
CLUSTER_API float ClusteringMaxDistanceToCluster(Clustering* This,int PointIndex, Cluster* Cluster);



  /*
   *  Base-class implementations: 
   */
CLUSTER_API void BASEClusteringFree(Clustering* This);
CLUSTER_API void BASEClusteringMergeClusters(Clustering* This,ClusterNode* NodeA,ClusterNode* NodeB);
CLUSTER_API float BASEClusteringDistanceToCentroid(Clustering* This,int PointIndex, Cluster* Cluster);
CLUSTER_API void BASEClusteringSplitCluster(Clustering* This,Cluster* ParentCluster, int PointA, int PointB);
CLUSTER_API float BASEClusteringPointToPointDistance(Clustering* This, int A,int B);
CLUSTER_API float BASEClusteringCentroidToCentroidDistance(Clustering* This, ClusterNode* NodeA,ClusterNode* NodeB);
CLUSTER_API Cluster* BASEClusteringGetNewCluster(Clustering* This);
CLUSTER_API int* BASEClusteringFindKmeansSeeds(Clustering* This, int Seeds);
CLUSTER_API void BASEClusteringClusterLinkage(Clustering* This, int DesiredClusterCount,float Epsilon);
CLUSTER_API void BASEClusteringClusterMeans(Clustering* This,int* SeedPoints, int DesiredClusterCount, int iteration, int mode);
CLUSTER_API void BASEClusteringClusterCobweb(Clustering* This, int DesiredClusterCount);
CLUSTER_API CobwebCluster* BASEClusteringNewCobwebCluster(Clustering* This, int FrameIndex);
CLUSTER_API void BASECobwebClusterComputeMeans(Clustering* This, CobwebTreeNode* Node);
CLUSTER_API void BASEClusteringCobwebAddNewFrame(Clustering* This,CobwebTreeNode* Root,CobwebTreeNode* Frame, int level);
CLUSTER_API void BASECobwebClusterComputeMeans(Clustering* This, CobwebTreeNode* Node);
CLUSTER_API float BASEClusteringCobwebCategoryUtility(Clustering* This, CobwebTreeNode* Parent);
CLUSTER_API void BASEClusteringCobwebClusterFree(Clustering* This, CobwebCluster* DeadGuy);
CLUSTER_API void BASECobwebClusterCopyMeans(Clustering* This, CobwebCluster* Source, CobwebCluster* Destination);
CLUSTER_API CobwebTreeNode* CobwebClusterMergeNodes(Clustering*, CobwebTreeNode*,CobwebTreeNode*);
CLUSTER_API CobwebTreeNode* ReadCobwebTree(Clustering *, char *);

CLUSTER_API BayesianCluster* GenerateBayesianSeedClusters(Clustering*, int);

CLUSTER_API int BASEClusteringGetAttributeCount(Clustering* This);
CLUSTER_API float BASEClusteringGetAttributeMin(Clustering* This, int AttributeIndex);
CLUSTER_API float BASEClusteringGetAttributeMax(Clustering* This, int AttributeIndex);
CLUSTER_API float BASEClusteringGetAttributeValue(Clustering* This, int PointIndex, int AttributeIndex);
CLUSTER_API float BASEClusteringSOMGetDistance(Clustering* This, SOMNode* SOM, int PointIndex);
CLUSTER_API void BASEClusteringSOMLearn(Clustering* This, SOMNode* SOM, float* LearningRate, int PointIndex);

  /*
   *  Called at the start of each cycle through the clustering algorithm.
   *  Used for doing things like building up a cluster-tree. 
   */
CLUSTER_API void BASEClusteringOutputStatus(Clustering* This);

  /*
   *  VTables 
   */
static VirtualFunction ClusterVirtualFunctions[FCOUNT_CLUSTER] = 
  { (void *) BASEClusterFree,
    (void *) BASEClusterFindCentroid
  };

static VirtualFunction ClusteringVirtualFunctions[FCOUNT_CLUSTERING] = 
  { (void *) BASEClusteringFree, 
    (void *) BASEClusteringMergeClusters,
    (void *) BASEClusteringSplitCluster,
    (void *) BASEClusteringGetNewCluster,
    (void *) BASEClusteringFindKmeansSeeds,
    (void *) BASEClusteringClusterLinkage,
    (void *) BASEClusteringClusterMeans,
    (void *) BASEClusteringClusterCobweb,
    (void *) BASEClusteringNewCobwebCluster, 
    (void *) BASEClusteringCobwebAddNewFrame,
    (void *) BASECobwebClusterComputeMeans,
    (void *) BASEClusteringCobwebClusterFree, 
    (void *) BASECobwebClusterCopyMeans,
    (void *) BASEClusteringSOMLearn, 
    (void *) BASEClusteringOutputStatus,
  };

static FloatVirtualFunction ClusteringFloatVirtualFunctions[FFCOUNT_CLUSTERING] = \
   {BASEClusteringPointToPointDistance,
    BASEClusteringCentroidToCentroidDistance,
    BASEClusteringDistanceToCentroid,
	BASEClusteringCobwebCategoryUtility, 
    BASEClusteringGetAttributeValue,
	BASEClusteringGetAttributeMin,
    BASEClusteringGetAttributeMax,
    BASEClusteringSOMGetDistance};

static IntVirtualFunction ClusteringIntVirtualFunctions[FICOUNT_CLUSTERING] = \
   {BASEClusteringGetAttributeCount};

/******************************************************************/
  /*
   *  Other functions (not clustering or cluster methods) 
   */

  /*
   *  Allocate a new (non-leaf) CobwebTreeNode 
   */
CLUSTER_API CobwebTreeNode* CobwebTreeNodeNew();

  /*
   *  Allocate a leaf CobwebTreeNode corresponding to the given point
   */
CLUSTER_API CobwebTreeNode* CobwebTreeNodeNewLeaf(Clustering* This, int PointIndex);

  /*
   *  Free memory from a CobwebTreeNode and its children 
   */
CLUSTER_API void CobwebTreeNodeFree(CobwebTreeNode* Head);

  /*
   *  Free memory from a CobwebTreeNode, but DO NOT eliminate the children 
   */
CLUSTER_API void CobwebTreeNodeFreeNR(CobwebTreeNode* Head);

  /*
   *  Add a child to a CobwebTreeNode. 
   */
CLUSTER_API void CobwebTreeNodeAddChild(CobwebTreeNode* Parent, CobwebTreeNode* Child);
CLUSTER_API float TryCategoryUtility(Clustering* This, CobwebTreeNode* Root, CobwebTreeNode* Host, CobwebTreeNode* Frame); 

#endif /* CLUSTER_H */








