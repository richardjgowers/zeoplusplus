//	Software for single molecule analysis
//
//	windows_aux.cc
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : February 5th 2016

#include "windows_aux.hh"


/*************************************************************/
/***** FUNCTIONS FOR CHECKING IF THE COMPLEX IS A WINDOW *****/
/*************************************************************/

/* GetUnvisitedPoint
 *
 * Auxiliary function for the IsWindow function computation. 
 * Given a point, and the list of visited points, returns true if the point is contained in the list.
 *
 */
/*bool PointInList (	int PointId,
			int *VisitedList,
			int NumElements)
{

	for (int i = 0; i < NumElements; i++)
		if (VisitedList[i] == PointId)
			return true;

	return false;

}*/

/* GetUnvisitedPoint
 *
 * Auxiliary function for the IsWindow function computation. 
 * Given a point, its neighboring edges (2), and the previous point visited when travelling through the complex, it returns the one which was unvisited.
 * Example: Given A as local point, A-B, D-A as edges, and B as previous point, it returns D. 
 *
 */
/*Point *GetUnvisitedPoint(Edge **NeighborEdges, Point *LocalPoint, Point *PreviousPoint)
{

	// Check the first edge
	Point *PAux = NeighborEdges[0]->GetOpposite(LocalPoint);

	// If the point is the previous point, return the opposite of the second edge
	if (PAux->GetPointId() == PreviousPoint->GetPointId())
		return NeighborEdges[1]->GetOpposite(LocalPoint);

	// If the opposite point on the first edge was not the previous point, return the second (its assumed that it will be correct)
	return NeighborEdges[0]->GetOpposite(LocalPoint);
	
	


}*/

/* IsWindow
 *
 * Checks if the edges and vertices of the complex form a simple cycle with no trivial edges. This is a necessary condition for it to be a window of the molecule if its formed by bonding edges.
 * The function starts in the first point of the point list, and travels around the list until is not able to continue or returns to a point already visited. Then, if the point is the initial one, returns true (false in case any condition fails).
 *
 * Input:
 *	- Window: A complex containing the information of the hypothetical window edges and vertices.
 * Output:
 *	- True if the complex is actually a cycle.
 */

/*bool IsWindow(Complex *Window)
{

	// Error check
	if (Window->GetNPoints() == 0)
		return false;

	// List of points visited	
	int *PointsVisited = new int[Window->GetNPoints()];
	
	// Counters for the number of points and edges already considered
	int PCount = 0;

	// Auxiliary edge list
	Edge **NeighborEdges = new Edge*[2];

	int NumEdgesRead = 0;

	// Start at the initial point
	
	Point *AuxP = Window->GetPointList()[0], *PrevP = NULL, *SwapP = NULL;
	
	NeighborEdges = Window->GetEdgesByPoint(AuxP, NumEdgesRead);

	if (NumEdgesRead != 2)
		return false;

	PointsVisited[0] = AuxP->GetPointId();

	// Take one of the directions by reading the opposite point of the first edge
	PrevP = AuxP;
	AuxP = NeighborEdges[0]->GetOpposite(AuxP);

	PCount++;

	// Iterate over the whole list
	while (NumEdgesRead == 2 && !PointInList(AuxP->GetPointId(), PointsVisited, PCount) )
	{

		// Get the neighboring edges
		NeighborEdges = Window->GetEdgesByPoint(AuxP, NumEdgesRead);
		
		// Check if the number is correct
		if (NumEdgesRead != 2)
			return false;

		// Compute the point of interest
		SwapP = AuxP;
		AuxP = GetUnvisitedPoint(NeighborEdges, AuxP, PrevP);
		PrevP = SwapP;

		// Insert the new point visited into the visited list
		PointsVisited[PCount] = PrevP->GetPointId();
		
		PCount++;

		// Security check
		if (PCount > Window->GetNPoints() + 1)
			return false;

	}

	// End of loop: a cycle has been found
	
	// The number of points in the cycle should be the same as the number of vertices in the complex
	if (PCount != Window->GetNPoints())
		return false;

	return true;
	

}*/

/*********************************************/
/********** WINDOW AREA COMPUTATION **********/
/*********************************************/




/*********************************************/
/************* SPHERIC EXPANSION *************/
/*********************************************/


/*
 *
 */
Complex *CompressToSphere (	Complex *TargetComplex)
{

	Complex *CompressedComplex = new Complex();

	double Centre[DIM_SMA] = {0., 0., 0.};

	/////	COMPUTE THE CENTRE OF THE COMPLEX
	///	Get the mean position of all the points
	//	(1) Sum all the coordinates
	//	(2) Average among the number of points

	for (int i = 0; i < TargetComplex->GetNPoints(); i++)
	{
		for (int j = 0; j < DIM_SMA; j++)
		{		
			Centre[j] += TargetComplex->GetPointList()[i]->GetCoordinate(j);
		}

	}

	for (int i = 0; i < DIM_SMA; i++) {
		Centre[i] /= TargetComplex->GetNPoints();
	}


	/////	COMPUTE RADIUS OF THE PROJECTED SPHERE 
	///	Use the minimum distance among the distances between the centre and any point
	//	Certain goodness in distribution is assumed.

	/*double Coordinates[DIM_SMA];
	TargetComplex->GetPointList()[0]->PointCoordinates(Coordinates);
		
	double Radius = Distance( Coordinates, Centre );*/

	double Radius = 0.;

	for (int i = 0; i < TargetComplex->GetNPoints(); i++)
	{
		double Coordinates[DIM_SMA];
		
		TargetComplex->GetPointList()[i]->PointCoordinates(Coordinates);

		double PDist = Distance(Coordinates, Centre);

		if ( PDist > Radius) {
		//if (PDist < Radius) {
			Radius = PDist;
		}

	}

	/////	PROJECT THE POINTS
	///	Projection in the direction marked by the C-P vector, until the distance is the computed radius.
	//	(1) Get the vector in the right direction (Centre-Point)
	//	(2) Normalize the vector and multiply by the desired radius
	// 	(3) Sum the direction vector to the centre to get the position of the compressed point
	//	(4) Insert the recomputed point into the compressed complex
	//	(5) Copy the edges from the old complex into the new one (with adequate references to points)

	for (int i = 0; i < TargetComplex->GetNPoints(); i++)
	{
		// (1)
		double CPVector[DIM_SMA];

		Point *PAux = TargetComplex->GetPointList()[i];

		for (int j = 0; j < DIM_SMA; j++) {
			CPVector[j] = PAux->GetCoordinate(j) - Centre[j];
		}

		// (2)
		double Factor = Radius/Module(CPVector);
	
		for (int j = 0; j < DIM_SMA; j++) {
			CPVector[j] *= Factor;
		}

		// (3)
		double NPCoordinates[DIM_SMA];

		for (int j = 0; j < DIM_SMA; j++) {
			NPCoordinates[j] = Centre[j] + CPVector[j];
		}

		// (4)
		CompressedComplex->InsertPointByValue(PAux->GetPointId(), NPCoordinates, PAux->GetRadius(), PAux->GetType() );

	}

	// (5)
	for (int j = 0; j < TargetComplex->GetNEdges(); j++)
	{
		Edge *EAux = TargetComplex->GetEdgeList()[j];

		// Get the pointers to the points with same Id's in the Compressed complex
		Point *PSource, *PTarget;

		PSource = CompressedComplex->GetPointById( EAux->GetPSource()->GetPointId() );

		PTarget = CompressedComplex->GetPointById( EAux->GetPTarget()->GetPointId() );

		CompressedComplex->InsertEdgeByValue( PSource, PTarget, EAux->GetChannelRadius() );

	}


	// Insert centre
	//CompressedComplex->InsertPointByValue(0, Centre, 0.0, atom_point );

	cout << "COMPRESSING TO SPHERE: " << endl << "CENTRE -> (" << Centre[0] << ", " << Centre[1] << ", " << Centre[2] << ")" << endl << "RADIUS -> " << Radius << endl;
	/*for (int i = 0; i < CompressedComplex->GetNPoints(); i++)
	{
		cout << "Point : " << CompressedComplex->GetPointList()[i]->GetPointId() << " - Distance to centre: " << Distance(CompressedComplex->GetPointList()[i]->GetCoordinates(), Centre) << endl;
	}*/

	return CompressedComplex;

}













// WINDOW CC
