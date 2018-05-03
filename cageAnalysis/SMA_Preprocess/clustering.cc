/*
 * clustering.cc
 *
 *  Created on: 10/3/2016
 *      Author: Ismael
 */


#include "clustering.hh"

bool PointsWithinDistance (	double Point[DIM_SMA],
							double CPoint[DIM_SMA],
							double DistThreshold)
{

	double V12[DIM_SMA] = {CPoint[0] - Point[0], CPoint[1] - Point[1], CPoint[2] - Point[2]};
	double Dist = sqrt( pow(V12[0], 2.) + pow(V12[1], 2.) + pow(V12[2], 2.) );
	//double Dist = 0.0;

	if (Dist <= DistThreshold)
		return true;

	return false;
}

/*
 * Check if the given points are related (meaning they belong to the same cluster).
 * Two criterions:
 * 	Distance - Points have to be close enough to be in the same cluster.
 * 	Point type - Points must have same type to be in the same cluster.
 * Warning: Two points may not be related and belong to the same cluster (through other points).
 */
bool PointsRelated(	Point *LPoint,
			Point *CPoint,
			double DistThreshold)
{

	double CPCoord[DIM_SMA], LPCoord[DIM_SMA];

	LPoint->PointCoordinates(LPCoord);
	CPoint->PointCoordinates(CPCoord);

	bool C1 = PointsWithinDistance(LPCoord, CPCoord, DistThreshold);
	bool C2 = (LPoint->GetType() == CPoint->GetType());

	return (C1 && C2);

}


bool ClusterPointsDensityBased (	list<Point *> &Points,
									list <Cluster *> &Clustering,
									double DistanceThreshold)
{

	int ClusterCounter = 0;

	while (Points.size() > 0)
	{
		bool ClusterFlag = true;

		// Initialize cluster and iterators
		Cluster *ClusterT = new Cluster;
		list<Point *>::iterator CPointsIt, CPointsStop, PointsFirst;

		// Insert the first unclassified point into cluster
		// Initialize iterators for first round
		PointsFirst = Points.begin();
		Point *NewCPoint; //= new Point();
		NewCPoint =	*PointsFirst;

		ClusterT->PointList.push_back(NewCPoint);
		Points.pop_front();

		// Iterator for rounds
		list<Point *>::iterator CPointsItSwap;

		CPointsItSwap = ClusterT->PointList.begin();
		CPointsStop = ClusterT->PointList.end();

		// Iterate until cluster is complete (all points have been added)
		while (ClusterFlag)
		{
			ClusterFlag = false;

			// Flag for point search round control
			bool RoundFlag = true;

			// Assign the iterators for next round
			CPointsIt = CPointsItSwap;
			CPointsStop = ClusterT->PointList.end();

			// Iterate over points on current cluster
			// Iteration starts in the first point added in last round,
			//	ending at the last element before starting that round.
			while (	CPointsIt != CPointsStop)
			{
				Point *CPoint; // = new Point();
				CPoint = *CPointsIt;

				//double CPCoord[DIM_SMA];
				//CPoint->PointCoordinates(CPCoord);

				// Iterator for the point list
				list<Point *>::iterator PointsIt;

				// Iterate over points in list
				for (PointsIt = Points.begin(); PointsIt != Points.end(); PointsIt++)
				{
					Point *LPoint; // = new Point();
					LPoint = *PointsIt;

					//double LPCoord[DIM_SMA];

					//LPoint->PointCoordinates(LPCoord);


					if (PointsRelated(CPoint, LPoint, DistanceThreshold) )
					{
						// Insert the new point into the cluster
						ClusterT->PointList.push_back(LPoint);

						// Erase the point from point list
						PointsIt = Points.erase(PointsIt);
						PointsIt--;

						// New point has been added: keep iterating over same cluster
						ClusterFlag = true;

						// First time an element is inserted in current round, set the next round's first
						if (RoundFlag)
						{
							// Point the swap pointer to the last element of the cluster
							//	This is only happening at every first insertion of each round
							CPointsItSwap = ClusterT->PointList.end();
							CPointsItSwap--;
							RoundFlag = false;
						}

					}

				}	// END_FOR (All unclassified points have been checked for current cluster)

				CPointsIt++;


			} // END_WHILE (All cluster points recently added have been checked)

		}	// END_WHILE (ClusterFlag is false)

		ClusterT->ClusterId = ClusterCounter;
		ClusterCounter++;

		Clustering.push_back(ClusterT);

	}	// END_WHILE (Points is empty)

	//cout << Clustering.size() << " clusters identified" << endl;

	return true;

}



/*	
 *	Computes the centroid of the cluster and returns it as a point. 
 */
Point *GetClusterCentroid(Cluster *ClusterT)
{

	list<Point *>::iterator CPointsIt;

	double AvgCoords[DIM_SMA] = {0., 0., 0.};

	int PointId = ClusterT->ClusterId;
	double Radius = 0.0;
	PointType Type = undefined_point;

	// Run through the list of points, get the coordinates, and accumulate them
	for (CPointsIt = ClusterT->PointList.begin(); CPointsIt != ClusterT->PointList.end(); CPointsIt++)
	{
		Point *PAux = *CPointsIt;

		// Store Id and Type from first point
		if (CPointsIt == ClusterT->PointList.begin())
		{
			//PointId = PAux->GetPointId();
			Type = PAux->GetType();
		}

		double LocalCoords[DIM_SMA];
		PAux->PointCoordinates(LocalCoords);

		AvgCoords[X_COORD] += LocalCoords[X_COORD];
		AvgCoords[Y_COORD] += LocalCoords[Y_COORD];
		AvgCoords[Z_COORD] += LocalCoords[Z_COORD];

		// TODO: Fix the way to compute the radius!!!!
		//Radius += PAux->GetRadius();

	}

	int ClusterNPoints = ClusterT->PointList.size();

	AvgCoords[X_COORD] /= ClusterNPoints;
	AvgCoords[Y_COORD] /= ClusterNPoints;
	AvgCoords[Z_COORD] /= ClusterNPoints;

	// TODO: Fix this!!!!
	//Radius /= ClusterNPoints;
	Radius = 0.0;

	Point *ClusterCentroid = new Point(PointId, AvgCoords, Radius, Type);

	return ClusterCentroid;

}

bool ConnectedClusters(	Cluster *Cluster1,
			Cluster *Cluster2,
			Complex *TargetComplex)
{
	list<Point *>::iterator CIt1, CIt2;

	for (CIt1 = Cluster1->PointList.begin(); CIt1 != Cluster1->PointList.end(); CIt1++)
	{
		Point *Cluster1P = *CIt1;

		for (CIt2 = Cluster2->PointList.begin(); CIt2 != Cluster2->PointList.end(); CIt2++)
		{
			Point *Cluster2P = *CIt2;

			// If there is an edge connecting both points, the clusters are connected
			//if (TargetComplex->GetEdgeByPoints(Cluster1P, Cluster2P) != NULL)
			if (TargetComplex->GetEdgeByIds( Cluster1P->GetPointId(), Cluster2P->GetPointId() ) != NULL)
			{
				return true;
			}
		}
	}

	return false;

}

/*
 * 	Computes the radius of a cluster center based on the distance to the atoms (the minimum distance among the distances to the atoms is the chosen one).
 */

void CentroidRadius(	Point *Centroid,
						Complex *Chemical)
{
	double Radius = 0.;

	for (int i = 0; i < Chemical->GetNPoints(); i++)
	{
		Point *ChemPoint = Chemical->GetPointList()[i];

		///		Compute distance to atoms taking into account their radius
		double DtP = ChemPoint->DistanceToPoint(Centroid) - ChemPoint->GetRadius();

		if (DtP < 0.)
		{
			DtP = 0.;
		}

		if (i == 0)
		{
			Radius = DtP;
		}
		else if (Radius > DtP)
		{
			Radius = DtP;
		}

	}

	//Centroid->SetRadius(Radius - AtomRadius);
	Centroid->SetRadius(Radius);
}

/*
 * 	Computes the radius of a cluster edge based on the distance to the atoms (the most restrictive one is chosen).
 * 	TODO - TAKE INTO ACCOUNT ATOM RADII
 */

void ClusterEdgeRadius(	Edge *ClusterEdge,
						Complex *Chemical)
{
	double EdgeRadius = FREE_EDGE_WEIGHT;

	double EdgeLength = ClusterEdge->GetPSource()->DistanceToPoint(ClusterEdge->GetPTarget());

	double EdgeVector[DIM_SMA] = {	ClusterEdge->GetPTarget()->GetCoordinate(X_COORD) -  ClusterEdge->GetPSource()->GetCoordinate(X_COORD),
									ClusterEdge->GetPTarget()->GetCoordinate(Y_COORD) -  ClusterEdge->GetPSource()->GetCoordinate(Y_COORD),
									ClusterEdge->GetPTarget()->GetCoordinate(Z_COORD) -  ClusterEdge->GetPSource()->GetCoordinate(Z_COORD)	};

	for (int i = 0; i < Chemical->GetNPoints(); i++)
	{
		Point *ChemPoint = Chemical->GetPointList()[i];

		double PDist = ClusterEdge->GetPSource()->DistanceToPoint(ChemPoint);

		double VectorToPoint[DIM_SMA] = {	ChemPoint->GetCoordinate(X_COORD) -  ClusterEdge->GetPSource()->GetCoordinate(X_COORD),
											ChemPoint->GetCoordinate(Y_COORD) -  ClusterEdge->GetPSource()->GetCoordinate(Y_COORD),
											ChemPoint->GetCoordinate(Z_COORD) -  ClusterEdge->GetPSource()->GetCoordinate(Z_COORD)	};

		double PAngle = Angle(EdgeVector, VectorToPoint);

		double ProjectionDist = PDist*cos(PAngle);

		if (ProjectionDist <= EdgeLength && PAngle <= PI/2.)
		{
			double PHeight = PDist*sin(PAngle);

			if (EdgeRadius > PHeight)
			{
				EdgeRadius = PHeight;
			}

		}



	}

	/*if (EdgeRadius == FREE_EDGE_WEIGHT)
	{
		ClusterEdge->GetPSource()->SetType(hybrid_point);
		ClusterEdge->GetPTarget()->SetType(hybrid_point);
	}*/

	ClusterEdge->SetChannelRadius(EdgeRadius);

}

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
								double DistanceThreshold)
{

	// List that will contain the clusters after applying the clustering
	list <Cluster *> Clustering;

	// Store the points of the complex into a list for processing
	// Pointers will point towards original objects in the complex

	list <Point *> Points;
	for (int i = 0; i < TargetComplex->GetNPoints(); i++)
	{
		Points.push_back(TargetComplex->GetPointList()[i]);
	}

	// Compute clusters
	ClusterPointsDensityBased(Points, Clustering, DistanceThreshold);

	//cout << "Clustering has " << Clustering.size() << " clusters" << endl;

	/////	CREATE CLUSTERIZED COMPLEX
	///		Create the complex
	//		Insert the points corresponding to the centroids of every cluster
	//		Compute the edges between clusters and insert them

	Complex *ClusteredComplex = new Complex();


	// Inserting points

	list<Cluster *>::iterator ClusterIt;

	for (ClusterIt = Clustering.begin(); ClusterIt != Clustering.end(); ClusterIt++)
	{
		Cluster *TCluster = *ClusterIt;

		Point *Centroid = GetClusterCentroid(TCluster);

		CentroidRadius(Centroid, Chemical);

		ClusteredComplex->InsertPointByValue(Centroid);
	}

	// Inserting edges
	list<Cluster *>::iterator CIt1, CIt2;

	for (CIt1 = Clustering.begin(); CIt1 != Clustering.end(); CIt1++)
	{
		Cluster *Cluster1 = *CIt1;

		for (CIt2 = next(CIt1, 1); CIt2 != Clustering.end(); CIt2++)
		{
			Cluster *Cluster2 = *CIt2;

			if ( ConnectedClusters(Cluster1, Cluster2, TargetComplex) )
			{
				Edge *ClusterEdge = new Edge(ClusteredComplex->GetPointById(Cluster1->ClusterId), ClusteredComplex->GetPointById(Cluster2->ClusterId), 0.0);

				ClusterEdgeRadius(ClusterEdge, Chemical);

				ClusteredComplex->InsertEdgeByValue(ClusterEdge);
			}

		}
	}


	list <Cluster *>::iterator CIt;

	for (CIt = Clustering.begin(); CIt != Clustering.end(); CIt++)
	{
		Cluster *Cluster0 = *CIt;
		delete Cluster0;
	}

	//Clustering.erase(Clustering.begin(), Clustering.end() );
	//Clustering.clear();
	Points.clear();

	// TODO: Remove this (only for debugging)
	/*ClusteredComplex->PrintComplex();

	char Filename[FILENAME_MAX];

	getwd(Filename);
	strncat(Filename, "/Molecules/Clustered.mol", FILENAME_MAX);
	PrintToMolFile (ClusteredComplex, Filename);*/
	// Till here



	return ClusteredComplex;


}





