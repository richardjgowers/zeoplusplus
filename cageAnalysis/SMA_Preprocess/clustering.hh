/*
 * clustering.hh
 *
 *  Created on: 10/3/2016
 *      Author: Ismael
 */

#ifndef CLUSTERING_HH_
#define CLUSTERING_HH_

#include <cmath>
#include <list>

#include <iostream>

#include <boost/utility.hpp>

//#include "../SMA_Core/complex.hh"

//#include "../SMA_Interface/complex_io.hh"

#include "../SMA_Core/molecule_info.hh"


using namespace std;

using namespace boost;

#define DISTANCE_THRES	0.85

typedef struct Cluster_
{
	list<Point *> 	PointList;
	int				ClusterId;

} Cluster;

#define FREE_EDGE_WEIGHT	1000.


/* ClusterizeComplex
 * 	Given a complex, it applies density-based clustering to merge those points that are close enough to each other to be considered the same.
 * 	Edges are recomputed considering if they connected points within clusters (i.e., there will be an edge between two clusters if there was
 * 	an edge between any point of the first cluster and any point of the second).
 *
 * Inputs:
 * 	TargetComplex - Complex to be clusterized (typically, the Voronoi graph of the molecule).
 *
 * Outputs:
 * 	Clusterized complex.
 *
 */
Complex *ClusterizeComplex (	Complex *TargetComplex,
								Complex *Chemical,
								double DistanceThreshold);



#endif /* CLUSTERING_HH_ */
