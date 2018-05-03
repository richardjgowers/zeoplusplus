//	Software for single molecule analysis
//
//	windows.cc
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : January 25th 2016

#include "window_comp.hh"


/*************************************************************/
/***** FUNCTIONS FOR CHECKING IF THE COMPLEX IS A WINDOW *****/
/*************************************************************/

/* GetUnvisitedPoint
 *
 * Auxiliary function for the IsWindow function computation.
 * Given a point, and the list of visited points, returns true if the point is contained in the list.
 *
 */
bool PointInList (	int PointId,
			int *VisitedList,
			int NumElements)
{

	for (int i = 0; i < NumElements; i++)
		if (VisitedList[i] == PointId)
			return true;

	return false;

}

/* GetUnvisitedPoint
 *
 * Auxiliary function for the IsWindow function computation.
 * Given a point, its neighboring edges (2), and the previous point visited when travelling through the complex, it returns the one which was unvisited.
 * Example: Given A as local point, A-B, D-A as edges, and B as previous point, it returns D.
 *
 */
Point *GetUnvisitedPoint(Edge **NeighborEdges, Point *LocalPoint, Point *PreviousPoint)
{

	// Check the first edge
	Point *PAux = NeighborEdges[0]->GetOpposite(LocalPoint);

	// If the point is the previous point, return the opposite of the second edge
	if (PAux->GetPointId() == PreviousPoint->GetPointId())
		return NeighborEdges[1]->GetOpposite(LocalPoint);

	// If the opposite point on the first edge was not the previous point, return the second (its assumed that it will be correct)
	return NeighborEdges[0]->GetOpposite(LocalPoint);




}

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

bool IsWindow(Complex *Window)
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


}

/*********************************************/
/********** WINDOW AREA COMPUTATION **********/
/*********************************************/

/**double FacetArea (	double Point1[DIM_SMA],
					double Point2[DIM_SMA],
					double Point3[DIM_SMA])
{

	// Vectors for distance calculation
	double U[DIM_SMA], E[DIM_SMA], ProyE[DIM_SMA], V[DIM_SMA];

	// Modules for vectors
	double Base, Height;

	// Lengths
	double Lambda1;

	/////	COMPUTE THE BASE AND HEIGHT OF THE TRIANGLE
	///	Take base as the distance between Point 1 and Point 2
	//	The height will be the length of the orthogonal proyection of Point 3 into the P1-P2 vector.


	//	VECTOR SUBSPACE
	//	Straight line connecting Point 1 and Point 2
	// 	U vector: unitary vector base of the vector subspace: S = <U>
	for (int i = 0; i < DIM_SMA; i++)
		U[i]=Point2[i]-Point1[i];

	Base = sqrt(U[X_COORD]*U[X_COORD] + U[Y_COORD]*U[Y_COORD] + U[Z_COORD]*U[Z_COORD]);

	for (int i = 0; i < DIM_SMA; i++)
		U[i]=U[i]/Base;

	//	ORTHOGONAL PROJECTION
	//	Vector formed by Point1-Centre is projected into <U>
	//	E = Centre-Point1
	//	U := p(E)

	// E vector: vector connecting Point1 with center (not unitary)
	for (int i = 0; i < DIM_SMA; i++)
		E[i]=Point3[i]-Point1[i];

	// Compute the orthogonal projection of the centre point over the Point1-Point2 line
	//Lambda1 = (U[X_COORD]*E[X_COORD] + U[Y_COORD]*E[Y_COORD] + U[Z_COORD]*E[Z_COORD])/(U[X_COORD]*U[X_COORD] + U[Y_COORD]*U[Y_COORD] + U[Z_COORD]*U[Z_COORD]);
	Lambda1 = U[X_COORD]*E[X_COORD] + U[Y_COORD]*E[Y_COORD] + U[Z_COORD]*E[Z_COORD];

	// U vector redefined: orthogonal projection of E into the subspace generated by U
	for (int i = 0; i < DIM_SMA; i++)
		ProyE[i] = U[i]*Lambda1;

	// 	DISTANCE
	//	Compute V = E-p(E)
	//	Compute ||V||

	for (int i = 0; i < DIM_SMA; i++)
		V[i]=E[i]-ProyE[i];

	Height = sqrt( V[X_COORD]*V[X_COORD] + V[Y_COORD]*V[Y_COORD] + V[Z_COORD]*V[Z_COORD] );


	/////	AREA OF THE TRIANGLE
	///	Base*Height/2

	return Base*Height/2;


}*/


double MaximumNeighDist( 	Complex *Chemical,
							Point *TargetPoint)
{

	int NBonds = -1;

	Edge **PointBonds = Chemical->GetEdgesByPoint(TargetPoint, NBonds);

	double MaxLength = 0.0;

	double TPCoordinates[DIM_SMA] = {TargetPoint->GetCoordinate(X_COORD), TargetPoint->GetCoordinate(Y_COORD), TargetPoint->GetCoordinate(Z_COORD)};

	for (int i = 0; i < NBonds; i++)
	{
		Point *Opposite = PointBonds[i]->GetOpposite(TargetPoint);

		double OpCoordinates[DIM_SMA] = {Opposite->GetCoordinate(X_COORD), Opposite->GetCoordinate(Y_COORD), Opposite->GetCoordinate(Z_COORD)};

		double DistToOp = Distance(TPCoordinates, OpCoordinates);

		if (MaxLength < DistToOp)
		{
			MaxLength = DistToOp;
		}

	}

	return MaxLength;

}


double MinimumNeighDist( 	Complex *Chemical,
							Point *TargetPoint)
{

	int NBonds = -1;

	Edge **PointBonds = Chemical->GetEdgesByPoint(TargetPoint, NBonds);

	// TODO: Use a maximum const defined somewhere
	double MinDist = 1000.0;

	double TPCoordinates[DIM_SMA] = {TargetPoint->GetCoordinate(X_COORD), TargetPoint->GetCoordinate(Y_COORD), TargetPoint->GetCoordinate(Z_COORD)};

	for (int i = 0; i < NBonds; i++)
	{
		Point *Opposite = PointBonds[i]->GetOpposite(TargetPoint);

		double OpCoordinates[DIM_SMA] = {Opposite->GetCoordinate(X_COORD), Opposite->GetCoordinate(Y_COORD), Opposite->GetCoordinate(Z_COORD)};

		double DistToOp = Distance(TPCoordinates, OpCoordinates);

		if (MinDist > DistToOp)
		{
			MinDist = DistToOp;
		}

	}

	return MinDist;

}

/*	NewWindowFromFacet
 *		Auxiliary function for EntryWindows. Creates a complex for containing the information about the window.
 *
 */

void NewWindowFromFacet( 	MoleculeInfo *Molecule,
							Facet OriginFacet)
{

	// Create the point list
	//Point **PointList = new Point*[3];

	int OppositeVertId = OriginFacet.second;

	//Complex *WindowT = new Complex();
	Window *WindowT = new Window();

	///// 	INSERT POINTS IN WINDOW
	///	Run around the cell and identify the point not remaining to the facet.
	//	Insert point by coordinates (radius and type are not important here).

	for (int i = 0; i < CELL_N_V; i++)
	{
	
		if (i != OppositeVertId)
		{
			// Create coordinates and insert points
			double Coordinates[DIM_SMA] = {OriginFacet.first->vertex(i)->point().x(), OriginFacet.first->vertex(i)->point().y(), OriginFacet.first->vertex(i)->point().z() };

			//int PointId = Molecule->GetChemical()->GetPointByCoordinates(Coordinates)->GetPointId();
			//WindowT->InsertPointByValue(PointId, Coordinates, 0.0, undefined_point);

			if (Molecule->GetChemical()->GetPointByCoordinates(Coordinates) != NULL)
			{
				int PointId = Molecule->GetChemical()->GetPointByCoordinates(Coordinates)->GetPointId();
				WindowT->InsertPointByValue(PointId, Coordinates, 0.0, undefined_point);
			}

			
		}
				
	}

	///// 	INSERT EDGES IN WINDOW
	///	Run around the cell and identify the edges not remaining to the facet.
	//	Insert edges only if they are bonding edges

	for (int i = 0; i < CELL_N_V; i++) 
	{
		for (int j = i+1; j < CELL_N_V; j++)
		{
			if ( (i != OppositeVertId) && (j != OppositeVertId) )
			{
				double Coordinates1[DIM_SMA] = {OriginFacet.first->vertex(i)->point().x(), OriginFacet.first->vertex(i)->point().y(), OriginFacet.first->vertex(i)->point().z()};
				double Coordinates2[DIM_SMA] = {OriginFacet.first->vertex(j)->point().x(), OriginFacet.first->vertex(j)->point().y(), OriginFacet.first->vertex(j)->point().z()};

				// Check if the edge is a bonding edge, then insert it into the window.
				if (Molecule->GetChemical()->GetEdgeByCoordinates(Coordinates1, Coordinates2) != NULL)
				{
					WindowT->InsertEdgeByValue(WindowT->GetPointByCoordinates(Coordinates1), WindowT->GetPointByCoordinates(Coordinates2));
				}
			}
		}
				
	}

	Molecule->InsertWindow(WindowT);

}

/*	NewWindowFromFacet
 *		ALTERNATIVE FORM - Auxiliary function for EntryWindows. Creates a new window containing vertices and edges (only bonding ones). Does NOT insert the window into the MoleculeInfo structure.
 *
 */

void NewWindowFromFacet( 	MoleculeInfo *Molecule,
							Facet OriginFacet,
							Window *WindowT)
{

	int OppositeVertId = OriginFacet.second;

	//WindowT = new Window();

	///// 	INSERT POINTS IN WINDOW
	///	Run around the cell and identify the point not remaining to the facet.
	//	Insert point by coordinates (radius and type are not important here).

	for (int i = 0; i < CELL_N_V; i++)
	{

		if (i != OppositeVertId)
		{
			// Create coordinates and insert points
			double Coordinates[DIM_SMA] = {OriginFacet.first->vertex(i)->point().x(), OriginFacet.first->vertex(i)->point().y(), OriginFacet.first->vertex(i)->point().z() };

			//int PointId = Molecule->GetChemical()->GetPointByCoordinates(Coordinates)->GetPointId();
			//WindowT->InsertPointByValue(PointId, Coordinates, 0.0, undefined_point);

			if (Molecule->GetChemical()->GetPointByCoordinates(Coordinates) != NULL)
			{
				int PointId = Molecule->GetChemical()->GetPointByCoordinates(Coordinates)->GetPointId();
				WindowT->InsertPointByValue(PointId, Coordinates, 0.0, undefined_point);
			}

		}

	}

	///// 	INSERT EDGES IN WINDOW
	///	Run around the cell and identify the edges not remaining to the facet.
	//	Insert edges only if they are bonding edges

	for (int i = 0; i < CELL_N_V; i++)
	{
		for (int j = i+1; j < CELL_N_V; j++)
		{
			if ( (i != OppositeVertId) && (j != OppositeVertId) )
			{
				double Coordinates1[DIM_SMA] = {OriginFacet.first->vertex(i)->point().x(), OriginFacet.first->vertex(i)->point().y(), OriginFacet.first->vertex(i)->point().z()};
				double Coordinates2[DIM_SMA] = {OriginFacet.first->vertex(j)->point().x(), OriginFacet.first->vertex(j)->point().y(), OriginFacet.first->vertex(j)->point().z()};

				// Check if the edge is a bonding edge, then, in such case, insert it into the window
				if (Molecule->GetChemical()->GetEdgeByCoordinates(Coordinates1, Coordinates2) != NULL)
				{
					WindowT->InsertEdgeByValue(WindowT->GetPointByCoordinates(Coordinates1), WindowT->GetPointByCoordinates(Coordinates2));
				}
			}
		}

	}

	//Molecule->InsertWindow(WindowT);

}

/*	NeighboringFacetP
 *		Checks if facets 1 and 2 are neighbours, and, in that case, if they are separated by a bonding edge. If the first happens and the second doesn't, returns the id of the point opposite to the shared edge (referring to facet 2).
 */

int NeighboringFacetP(	Facet Facet1,
						Facet Facet2,
						Complex *Chemical)
{

	int IdPoints[CELL_N_V]; 	// Identifiers of the points
	int Count = 0;			// Counter of the number of coinciding points found

	int OppositeVert1 = Facet1.second, OppositeVert2 = Facet2.second;

	/////	CHECK IF FACETS ARE NEIGHBOURS
	///	They are neighbours if they share exactly one edge (i.e., exactly two points)
	//	If they share the 4 points, they are the same edge (should never happen)

	for (int i = 0; i < CELL_N_V; i++)
	{
		if (i != OppositeVert1)
		{
			for (int j = 0; j < CELL_N_V; j++)
			{
				if (j != OppositeVert2)
				{
					// Check if point i of facet1 and point j of face2 are the same
					if (	(Facet1.first->vertex(i)->point().x() == Facet2.first->vertex(j)->point().x()) &&
							(Facet1.first->vertex(i)->point().y() == Facet2.first->vertex(j)->point().y()) &&
							(Facet1.first->vertex(i)->point().z() == Facet2.first->vertex(j)->point().z()) )
					{
						IdPoints[Count] = j;	// Store points shared by facet2 with facet1
						Count++;
					}
				}
			}
		}
	}


	//////	CHECK IF SHARED BONDS
	///	Look at the molecule for the edge (if there's any)
	// 	Since the facets are triangles, they share an edge if Count >= 2 (N edges shared == Count -1), and are the same if Count == 4 (Count == 3 is not possible).s

	if (Count == 2)
	{
		double 	Coordinates1[DIM_SMA] = {Facet2.first->vertex(IdPoints[0])->point().x(), Facet2.first->vertex(IdPoints[0])->point().y(), Facet2.first->vertex(IdPoints[0])->point().z()}, 
				Coordinates2[DIM_SMA] = {Facet2.first->vertex(IdPoints[1])->point().x(), Facet2.first->vertex(IdPoints[1])->point().y(), Facet2.first->vertex(IdPoints[1])->point().z()};

		// Check if the edge is a bonding edge (and return -1 if that's the case)
		if ( Chemical->GetEdgeByPoints(Chemical->GetPointByCoordinates(Coordinates1), Chemical->GetPointByCoordinates(Coordinates2)) != NULL )
			return -1;
		else
		{
			// Return the index of the point of facet2 that was not in the edge
			for (int i = 0; i < CELL_N_V; i++)
				if (i != IdPoints[0] && i!= IdPoints[1] && i != OppositeVert2)
					return i;
		}
	}

	return -2;

}

/* 	AddPointToWindow
 *		Adds the point (belonging to the facet) to the window. Then it adds all the edges that connects that point to the other points of the window if they are bonding edges.
 *
 *	Inputs:	
 *
 *	Return:
 */
void AddPointToWindow( 	int PointIdF,
						Facet FacetT,
						Window *WindowT,
						Complex *Chemical)
{

	/////	INSERT THE POINT INTO THE WINDOW
	///		Compute the coordinates through the facet.
	//		Insert only the coordinates.

	double Coordinates[DIM_SMA] = {	FacetT.first->vertex(PointIdF)->point().x(), FacetT.first->vertex(PointIdF)->point().y(), FacetT.first->vertex(PointIdF)->point().z()};

	//int PointId = Chemical->GetPointByCoordinates(Coordinates)->GetPointId();
	//WindowT->InsertPointByValue(PointId, Coordinates, 0.0, undefined_point);

	if (Chemical->GetPointByCoordinates(Coordinates) != NULL)
	{
		int PointId = Chemical->GetPointByCoordinates(Coordinates)->GetPointId();
		WindowT->InsertPointByValue(PointId, Coordinates, 0.0, undefined_point);
	}
	else
	{
		return;
	}

	/////	INSERT THE BONDING EDGES LEFT
	///		Get the bonding edges connected to the point.
	//		Check if the points computed from those edges are into the window, then add the corresponding edges.

	int NConnBonds = -1;

	Point *MLocalPoint = Chemical->GetPointByCoordinates(Coordinates);

	//Edge **ConnectingBonds = Chemical->GetEdgesByPoint(Chemical->GetPointByCoordinates(Coordinates), NConnBonds);
	Edge **ConnectingBonds = Chemical->GetEdgesByPoint(MLocalPoint, NConnBonds);


	for (int i = 0; i < NConnBonds; i++)
	{
		Point *MOppositePoint = ConnectingBonds[i]->GetOpposite(MLocalPoint);
		Point *WOppositePoint = WindowT->GetPointByCoordinates(MOppositePoint->GetCoordinates());
		
		// Check if the opposite point of every bond is already in the window
		if (WOppositePoint != NULL)
		{
			// Insert a bond into the window connecting both the local and the opposite point
			WindowT->InsertEdgeByValue(WOppositePoint, WindowT->GetPointByCoordinates(Coordinates));
		}
	}


}


void AddFacetToWindow(	Facet FacetT,
						Window *WindowT)
{

	Simplex *WinFacet = new Simplex();

	//WinFacet->SetDim(2);
	WinFacet->SetDim(3);

	int OppositeVertId = FacetT.second;

	for (int i = 0; i < CELL_N_V; i++)
	{

		Point *SPoint;

		if (i != OppositeVertId)
		{
			SPoint = new Point();

			SPoint->SetCoordinate( X_COORD, FacetT.first->vertex(i)->point().x() );
			SPoint->SetCoordinate( Y_COORD, FacetT.first->vertex(i)->point().y() );
			SPoint->SetCoordinate( Z_COORD, FacetT.first->vertex(i)->point().z() );

			WinFacet->InsertPoint(SPoint);
		}

	}

	WindowT->InsertSimplex(WinFacet);



}


/* 	EntryWindows
 *		ALTERNATIVE FORM - Computes the entry windows to the molecule from the alpha shape, using the information about the bonds and the viccinity of the facets.
 *
 *	Inputs:
 *		Molecule 	- Pointer to molecule info object containing relevant chemical information (and eventually others).
 *		AlphaShape 	- Alpha shape of the given molecule (general alpha-shape: all simplices with their alpha interval).
 *		Alpha		- Alpha value chosen for the alpha shape.
 *	Return:
 *		True if the computation of the windows was performed succesfully.
 *		The windows computed are stored into the MoleculeInfo object passed as an argument.
 */


bool EntryWindows ( 	MoleculeInfo *Molecule,
						Alpha_shape_3 &AlphaShape,
						double Alpha,
						bool Relaxed)
{

	bool OnlyWindowsDetected = true;

	// List of facets for containing the total of facets of the alpha shape
	list<Facet> FacetList, SingFacetList;
	//list<Facet>::const_iterator FacetIt1, FacetIt2;

	// FOR LINUX
	list<Facet>::iterator FacetIt1, FacetIt2;

	// Get the list of facets of the alpha shape (Singular and Regular); append both lists.
	AlphaShape.get_alpha_shape_facets(back_inserter(FacetList), Alpha_shape_3::REGULAR, Alpha);
	AlphaShape.get_alpha_shape_facets(back_inserter(SingFacetList), Alpha_shape_3::SINGULAR, Alpha);

	FacetList.splice(FacetList.end(), SingFacetList);

	//cout << "ALPHA SHAPE TOTAL FACETS: " << FacetList.size() << endl;


	/////	COMPUTE WINDOWS
	///		Iterate over facets. Two facets belong to the same window if:
	//			(1) They share a non-bonding edge.
	//			(2) They are both of the same type (i.e., internal, regular, etc) -> This is already provided.
	//		Construct windows. Each window is a Complex to be inserted in the complex list.
	//			(3) When identifying a facet as being part of a window, add it to the same complex.

	//int FacetCount = 0;

	while (!FacetList.empty())
	{
		//int WFacetCount = 0;

		FacetIt1 = FacetList.begin();

		Facet WFacet0 = *FacetIt1;

		///		Create new window from original facet.
		//		Insert points, edges (only bonding ones) and facet (to facet list)

		//Window *WindowT = NULL;
		Window *WindowT = new Window();
		NewWindowFromFacet(Molecule, WFacet0, WindowT);
		//TODO - make this part of previous function
		AddFacetToWindow(WFacet0, WindowT);
		//WFacetCount++;

		///		Create a list of window facets to explore the rest of the window.
		//		(1) Insert the original facet into it.
		//		(2) Erase that facet from previous list.

		list<Facet> WinFacetList;
		//list<Facet>::const_iterator CurrentWinFacet;

		// FOR LINUX
		list<Facet>::iterator CurrentWinFacet;

		WinFacetList.push_back(WFacet0);

		FacetList.erase(FacetIt1);

		// FOR LINUX
		//FacetList.erase(FacetIt1, FacetIt1);

		///		Add window facets
		//		(1)	For each facet in FacetList, compare with the current window facet
		//		(2)	If they are neighbors (in the sense defined above), add the facet to the window.
		//		(3)	When all facets have been compared, step to the next window facet. If it was the last, stop.

		bool CurrentEqualsLast= false;

		CurrentWinFacet = WinFacetList.begin();

		while (!CurrentEqualsLast)
		{
			//CurrentEqualsLast = true;

			Facet WinFacet = *CurrentWinFacet;

			// (1)
			for (FacetIt2 = FacetList.begin(); FacetIt2 != FacetList.end(); FacetIt2++)
			{
				Facet Facet2 = *FacetIt2;

				// (2)
				int PointToAdd = NeighboringFacetP(WinFacet, Facet2, Molecule->GetChemical());

				if (PointToAdd >= 0)
				{
					// Add the point to the window (and the edges that are necessary to be add)
					AddPointToWindow(PointToAdd, Facet2, WindowT, Molecule->GetChemical());
					AddFacetToWindow(Facet2, WindowT);
					//WFacetCount++;

					// Erase the facet from the facet list moving iterator back
					FacetList.erase(FacetIt2);
					--FacetIt2;

					// Add the facet into the window facets list
					WinFacetList.push_back(Facet2);

					// Facet has been add (keep processing)
					//CurrentEqualsLast = false;

				}
			}	// End for (FacetList It2)

//			CurrentEqualsLast = next(CurrentWinFacet) == WinFacetList.end();

			// FOR LINUX
			CurrentEqualsLast = next(CurrentWinFacet, 1) == WinFacetList.end();

			// (3)
			if (!CurrentEqualsLast)
			{
				CurrentWinFacet++;
			}

		}	// End while (CurrentEqualsLast)

		//cout << "WINDOW HAS " << WFacetCount << " FACETS" << endl;
		//FacetCount += WFacetCount;


		///		Insert the window once completely constructed
		//		Step now into the next iteration (new window will be generated).

		if (IsWindow(WindowT))
		{
			Molecule->InsertWindow(WindowT);
		}
		else if (Relaxed)
		{
			//cout << "Relaxed window accepted" << endl;
			if (WindowT->Pruned())
			{
				Molecule->InsertWindow(WindowT, true);
			}
		}
		else
		{
			//cout << "Pruning window" << endl;
			/*Window *WindowPruned = WindowT->PruneWindow();

			if (IsWindow(WindowPruned))
			{
				Molecule->InsertWindow(WindowPruned);
			}
			else
			{
				OnlyWindowsDetected = false;
			}*/

			OnlyWindowsDetected = false;

			delete WindowT;

		}
		//Molecule->InsertWindow(WindowT);



	}	// End while (FacetList is empty)

	//cout << "TOTAL FACETS PROCESSED: " << FacetCount << endl;



	return OnlyWindowsDetected;

}

/*
 * 	Given a Voronoi graph, computes all the points in between all pairs of internal points. A copy of the Voronoi graph with the new points inside is returned.
 */
Complex *GetExpandedVoronoi (	Complex *VoroGraph)
{
	int IdCount = VoroGraph->GetNPoints();

	Complex *ExpandedVoro = VoroGraph->ComplexCopy();

	for (int i = 0; i < VoroGraph->GetNPoints(); i++)
	{
		Point *P1 = VoroGraph->GetPointList()[i];

		for (int j = i+1; j < VoroGraph->GetNPoints(); j++)
		{
			Point *P2 = VoroGraph->GetPointList()[j];

			Point *MidPoint;

			if (P1->GetType() == internal_point && P2->GetType() == internal_point)
			{
				double MidCoords[DIM_SMA] = { 	(P1->GetCoordinate(X_COORD) + P2->GetCoordinate(X_COORD))/2. ,
												(P1->GetCoordinate(Y_COORD) + P2->GetCoordinate(Y_COORD))/2. ,
												(P1->GetCoordinate(Z_COORD) + P2->GetCoordinate(Z_COORD))/2. };

				double MidRad = P1->DistanceToPoint(P2)/2.;

				MidPoint = new Point(IdCount++, MidCoords, MidRad, internal_point);

				ExpandedVoro->InsertPointByValue(MidPoint);
			}
		}

	}

	return ExpandedVoro;


}


/*	AllWindowsDetected
 *		Check if all entrying Voronoi edges do cross a facet.
 *
 */

bool AllWindowsDetected(	MoleculeInfo *Molecule)
{

	Complex *VoroGraph = Molecule->GetVoronoiGraph();


	for (int i = 0; i < VoroGraph->GetNEdges(); i++)
	{
		Edge *VoroEdge = VoroGraph->GetEdgeList()[i];

		bool EdgeCrossesWindow = false;

		if (VoroEdge->IsEntryEdge())
		{

			/// Get edge points coordinates
			double SourceCoords[DIM_SMA] = {	VoroEdge->GetPSource()->GetCoordinate(X_COORD),
												VoroEdge->GetPSource()->GetCoordinate(Y_COORD),
												VoroEdge->GetPSource()->GetCoordinate(Z_COORD)};

			double TargetCoords[DIM_SMA] = {	VoroEdge->GetPTarget()->GetCoordinate(X_COORD),
												VoroEdge->GetPTarget()->GetCoordinate(Y_COORD),
												VoroEdge->GetPTarget()->GetCoordinate(Z_COORD)};


			for (int j = 0; j < Molecule->GetNWindows(); j++)
			{
				Window *MolWindow = Molecule->GetWindows()[j];

				list <Simplex *> *Simplices = MolWindow->GetSimplices();
				list <Simplex *>::iterator SIt;

				for (SIt = Simplices->begin(); SIt != Simplices->end(); SIt++)
				{

					Simplex *WinSimplex = *(SIt);

					Point 	*P1 = WinSimplex->ReadPoint(0),
							*P2 = WinSimplex->ReadPoint(1),
							*P3 = WinSimplex->ReadPoint(2);

					double 	P1_Coords[DIM_SMA] = {P1->GetCoordinate(X_COORD), P1->GetCoordinate(Y_COORD), P1->GetCoordinate(Z_COORD)},
							P2_Coords[DIM_SMA] = {P2->GetCoordinate(X_COORD), P2->GetCoordinate(Y_COORD), P2->GetCoordinate(Z_COORD)},
							P3_Coords[DIM_SMA] = {P3->GetCoordinate(X_COORD), P3->GetCoordinate(Y_COORD), P3->GetCoordinate(Z_COORD)};

					if (RayCrossesTriangle ( SourceCoords, TargetCoords, P1_Coords, P2_Coords, P3_Coords) )
					{
						EdgeCrossesWindow = true;
						break;
					}


				}	// End FOR (Window simplices)

				if (EdgeCrossesWindow)
				{
					break;
				}

			}	// End FOR (Windows)

			///		Return false if the entry edge does not cross any facet from any window.
			if (!EdgeCrossesWindow)
			{
				//cout << "Failing edge: " << endl;
				//VoroEdge->GetPSource()->SetType(hybrid_point);
				//VoroEdge->GetPTarget()->SetType(hybrid_point);
				return false;
			}


		}		// END IF (Edge is entry edge)

	}

	return true;
}




/*	ComputeWindows_UW_Opt
 * 		Compute windows using a unweighted alpha-shape of alpha=optimal*1.1.
 *
 *
 */
bool ComputeWindows_UW_Opt(MoleculeInfo *Molecule, bool Relaxed)
{

	///		List of weighted points.
	list<Weighted_point> WeightedPList;


	///		Insert atoms into weighted point list.
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		Point *PAux = Molecule->GetChemical()->GetPointList()[i];


		///		Error control: the point is an atom
		if (PAux->GetType() == atom_point)
		{
			///		Compute the weight: maximum distance to neighbors
			//double Weight = MaximumNeighDist(Molecule->GetChemical(), PAux);

			//Weight = pow(Weight, 2.);
			double Weight = 0.0;

			///		Insert point into the weighted point list (weighted alpha-shapes)
			WeightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );
		}

	}

	///		Compute the alpha shape in general mode
	Alpha_shape_3 MoleculeAlphaShape(WeightedPList.begin(), WeightedPList.end(), 0, Alpha_shape_3::GENERAL);
	Alpha_iterator 	opt = MoleculeAlphaShape.find_optimal_alpha(1);
	double Alpha = *opt;

	Alpha += 0.1*Alpha;


	//////		WEIGHTED ALPHA SHAPE FOR WINDOW COMPUTING
	///			Alpha = 0.0
	//			Weights = square maximum distance among neighbors

	EntryWindows(Molecule, MoleculeAlphaShape, Alpha, Relaxed);

	///		Check if the windows have been detected.
	if(AllWindowsDetected(Molecule))
	{
		return true;
	}

	return false;
}



/*	ComputeWindows_Weighted0
 * 		Compute windows using a weighted alpha-shape of alpha=0 and weights equal to the maximum distance from any atom to it's (chemical) neighbours.
 * 		Useful method to detect small windows, especially if they're partially hidden.
 *
 */
bool ComputeWindows_Weighted0(MoleculeInfo *Molecule, bool Relaxed)
{

	///		List of weighted points.
	list<Weighted_point> WeightedPList;


	///		Insert atoms into weighted point list.
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		Point *PAux = Molecule->GetChemical()->GetPointList()[i];


		///		Error control: the point is an atom
		if (PAux->GetType() == atom_point)
		{
			///		Compute the weight: maximum distance to neighbors
			double Weight = MaximumNeighDist(Molecule->GetChemical(), PAux);

			Weight = pow(Weight, 2.);

			///		Insert point into the weighted point list (weighted alpha-shapes)
			WeightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );
		}

	}

	///		Compute the alpha shape in general mode
	Alpha_shape_3 MoleculeAlphaShape(WeightedPList.begin(), WeightedPList.end(), 0, Alpha_shape_3::GENERAL);


	//////		WEIGHTED ALPHA SHAPE FOR WINDOW COMPUTING
	///			Alpha = 0.0
	//			Weights = square maximum distance among neighbors

	EntryWindows(Molecule, MoleculeAlphaShape, 0.0, Relaxed);

	///		Check if the windows have been detected.
	if(AllWindowsDetected(Molecule))
	{
		return true;
	}

	return false;
}

/*	ComputeWindows_ExploratoryAlphas
 * 		Compute windows using an unweighted alpha-shape with variable alphas.
 * 		Useful method to detect average regular windows, but of low performance.
 *
 */
bool ComputeWindows_ExploratoryAlphas(MoleculeInfo *Molecule, bool Relaxed)
{

	///		List of weighted points.
	list<Weighted_point> PList;


	///		Insert atoms into weighted point list.
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		Point *PAux = Molecule->GetChemical()->GetPointList()[i];


		///		Error control: the point is an atom
		if (PAux->GetType() == atom_point)
		{
			///		Insert point into the weighted point list (weighted alpha-shapes)
			PList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), 0.0  ) );
		}

	}

	//////		ALPHA SHAPE FOR WINDOW COMPUTING
	///			Initial alpha = alpha_lower_bound
	//			Weights = square maximum distance among neighbors

	///		Compute the alpha shape in general mode
	Alpha_shape_3 MoleculeAlphaShape(PList.begin(), PList.end(), 0, Alpha_shape_3::GENERAL);
	NT alpha_NT = MoleculeAlphaShape.find_alpha_solid();
	Alpha_iterator AlphaIt;

	int Counter = 0;

	///		Iterate over all alphas
	for (AlphaIt = MoleculeAlphaShape.alpha_lower_bound(alpha_NT); AlphaIt != MoleculeAlphaShape.alpha_end(); AlphaIt++)
	{
		double Alpha = *AlphaIt;

		///		Compute the windows.
		EntryWindows(Molecule, MoleculeAlphaShape, Alpha, Relaxed);

		///		After certain number of iterations, check if allt the windows have been detected.
		if (Counter%MAX_ITER_CHECK == 0)
		{
			if(AllWindowsDetected(Molecule)) {	return true; }
		}
		Counter++;
	}

	///		Check if all the windows have been detected.
	return AllWindowsDetected(Molecule);
}

/*	ComputeWindows_ExploratoryAlphasW
 * 		Compute windows using weighted alpha-shapes with variable alphas. The weights are fixed to the maximum distance among the bonding neighbors.
 * 		Useful method to detect average regular windows, but of low performance.
 *
 */
bool ComputeWindows_ExploratoryAlphasW(MoleculeInfo *Molecule, bool Relaxed)
{

	///		List of weighted points.
	list<Weighted_point> PList;


	///		Insert atoms into weighted point list.
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		Point *PAux = Molecule->GetChemical()->GetPointList()[i];


		///		Error control: the point is an atom
		if (PAux->GetType() == atom_point)
		{
			///		Compute the weight: maximum distance to neighbors
			double Weight = MaximumNeighDist(Molecule->GetChemical(), PAux);

			///		Insert point into the weighted point list (weighted alpha-shapes)
			PList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );
		}

	}

	//////		ALPHA SHAPE FOR WINDOW COMPUTING
	///			Initial alpha = alpha_lower_bound
	//			Weights = square maximum distance among neighbors

	///		Compute the alpha shape in general mode
	Alpha_shape_3 MoleculeAlphaShape(PList.begin(), PList.end(), 0, Alpha_shape_3::GENERAL);
	NT alpha_NT = MoleculeAlphaShape.find_alpha_solid();
	Alpha_iterator AlphaIt;

	int Counter = 0;

	///		Iterate over all alphas
	for (AlphaIt = MoleculeAlphaShape.alpha_lower_bound(alpha_NT); AlphaIt != MoleculeAlphaShape.alpha_end(); AlphaIt++)
	{
		double Alpha = *AlphaIt;

		///		Compute the windows.
		EntryWindows(Molecule, MoleculeAlphaShape, Alpha, Relaxed);

		///		After certain number of iterations, check if allt the windows have been detected.
		if (Counter%MAX_ITER_CHECK == 0)
		{
			if(AllWindowsDetected(Molecule)) {	return true; }
		}
		Counter++;
	}

	///		Check if all the windows have been detected.
	return AllWindowsDetected(Molecule);
}


/*	ComputeWindows_GenericMesh
 * 		Compute windows using a weighted alpha-shape (alpha = alpha_optimal_1) where atoms and a generic mesh are used as reference points.
 * 		Useful method to detect huge irregular windows.
 *
 */
bool ComputeWindows_GenericMesh(MoleculeInfo *Molecule, bool Relaxed)
{

	///		List of weighted points.
	list<Weighted_point> PList;

	///		Insert atoms into weighted point list.
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		Point *PAux = Molecule->GetChemical()->GetPointList()[i];


		///		Error control: the point is an atom
		if (PAux->GetType() == atom_point)
		{
			///		Compute the weight: maximum distance to neighbors
			double Weight = MaximumNeighDist(Molecule->GetChemical(), PAux);

			///		Insert point into the weighted point list (weighted alpha-shapes)
			PList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );
		}

	}

	///	Compute grid and insert internal points into the Delaunay triangulation

	///		List for grid points. Each point has a weight that connects it with its neighbours.
	list <Point *> *PointList = new list<Point *>();
	list <Point *>::iterator PointIt;

	///		Compute grid
	InternalPointsGrid(Molecule, PointList);

	///		Inser grid points into point list
	for (PointIt = PointList->begin(); PointIt != PointList->end(); PointIt++)
	{
		Point *PAux = *(PointIt);

		// Check that the point is of desired type
		if (PAux->GetType() == internal_point)
		{
			// Insert point into the Delaunay triangulation (unweighted alpha-shapes)
			double Weight = PAux->GetRadius();

			Weight = pow(Weight, 2.);

			// Insert point into the weighted point list (weighted alpha-shapes)
			PList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );
		}

	}



	//////		ALPHA SHAPE FOR WINDOW COMPUTING
	///			Initial alpha = alpha_lower_bound
	//			Weights = square maximum distance among neighbors

	///		Compute the alpha shape in general mode
	Alpha_shape_3 MoleculeAlphaShape(PList.begin(), PList.end(), 0, Alpha_shape_3::GENERAL);
	//Alpha_iterator AlphaIt = MoleculeAlphaShape.find_optimal_alpha(1);

	//double Alpha = *(AlphaIt);

	///		Compute entry windows
	EntryWindows(Molecule, MoleculeAlphaShape, 0.0, Relaxed);

	///		Check if all the windows have been detected.
	return AllWindowsDetected(Molecule);
}

/*	ComputeWindows_VoroGrid
 * 		Compute windows using a weighted alpha-shape (alpha = alpha_optimal_1) where atoms and internal Voronoi nodes are used as reference points.
 * 		Useful method to detect average regular windows, medium performance.
 *
 */
bool ComputeWindows_VoroGrid(MoleculeInfo *Molecule, bool Relaxed)
{

	///		List of weighted points.
	list<Weighted_point> PList;

	///		Insert atoms into weighted point list.
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		Point *PAux = Molecule->GetChemical()->GetPointList()[i];


		///		Error control: the point is an atom
		if (PAux->GetType() == atom_point)
		{
			///		Compute the weight: maximum distance to neighbors
			double Weight = MaximumNeighDist(Molecule->GetChemical(), PAux);

			///		Insert point into the weighted point list (weighted alpha-shapes)
			PList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );
		}

	}

	///		Insert atoms into weighted point list.
	Complex *ExpandedVoro =  Molecule->GetVoronoiGraph();

	///	Insert Voronoi internal points into the Delaunay triangulation
	for (int i = 0; i < ExpandedVoro->GetNPoints(); i++)
	{
		Point *PAux = ExpandedVoro->GetPointList()[i];

		// Check that the point is of desired type
		if (PAux->GetType() == internal_point)
		{
			// Insert point into the Delaunay triangulation (unweighted alpha-shapes)
			double Weight = PAux->GetRadius();

			Weight = pow(Weight, 2.);

			// Insert point into the weighted point list (weighted alpha-shapes)
			PList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );
		}
	}

	//////		ALPHA SHAPE FOR WINDOW COMPUTING
	///			Initial alpha = alpha_lower_bound
	//			Weights = square maximum distance among neighbors

	///		Compute the alpha shape in general mode
	Alpha_shape_3 MoleculeAlphaShape(PList.begin(), PList.end(), 0, Alpha_shape_3::GENERAL);
	Alpha_iterator AlphaIt = MoleculeAlphaShape.find_optimal_alpha(1);

	double Alpha = *(AlphaIt);

	///		Compute entry windows
	EntryWindows(Molecule, MoleculeAlphaShape, Alpha, Relaxed);

	///		Check if all the windows have been detected.
	return AllWindowsDetected(Molecule);
}


/* 	ComputeWindows
 *		Calculates the windows of the molecule through its alpha shape, storing the information into the MoleculeInfo object pointed by the only argument of the function.
 *
 *	Input:
 *		Molecule - Pointer to complex containing atom and bonding information.
 *	Return:
 *		(In Molecule object) - Complex with added information: holes (nodes, internal and external) and channels between holes.
 */
bool ComputeWindows( MoleculeInfo *Molecule, int Options)
{

	bool CWFlag = false;

	/////	COMPUTE THE WINDOWS (STRICT MODE)
	///			Window computing in strict mode (only cyclic windows are considered).
	//			Different methods (according to top layer requirements) are applied.

	///		Compute windows: uw alpha shape with alpha optimal
	if (Options%UW_ALPHA_OPT == 0)
	{
		CWFlag = ComputeWindows_UW_Opt(Molecule, false);
		if (CWFlag) { return true; }
	}

	///		Compute windows: weighted alpha 0 criterion.
	if (Options%F_ALPHA_W0 == 0)
	{
		CWFlag = ComputeWindows_Weighted0(Molecule, false);
		if (CWFlag) { return true; }
	}

	///		Compute windows: Voronoi grid criterion.
	if (Options%VORO_GRID == 0)
	{
		CWFlag = ComputeWindows_VoroGrid(Molecule, false);
		if (CWFlag) { return true; }
	}

	///		Compute windows: Moving weighted alpha's criterion.
	if (Options%MV_ALPHA_W == 0)
	{
		CWFlag = ComputeWindows_ExploratoryAlphasW(Molecule, false);
		if (CWFlag) { return true; }
	}

	///		Compute windows: Moving unweighted alpha's criterion.
	if (Options%MV_ALPHA_UW == 0)
	{
		CWFlag = ComputeWindows_ExploratoryAlphas(Molecule, false);
		if (CWFlag) { return true; }
	}

	///		Compute windows: Generic mesh alpha shape criterion.
	if (Options%GEN_MESH == 0)
	{
		CWFlag = ComputeWindows_GenericMesh(Molecule, false);
		if (CWFlag) { return true; }
	}

	/////	COMPUTE THE WINDOWS (RELAXED MODE)
	///			Window computing in relaxed mode (non cyclic windows may be accepted).
	//			Different methods (according to top layer requirements) are applied.

//	cout << "Window computing entered in relaxed mode" << endl;

	///		Compute windows: weighted alpha 0 criterion.
	if (Options%F_ALPHA_W0 == 0)
	{
		CWFlag = ComputeWindows_Weighted0(Molecule, true);
		if (CWFlag) { return true; }
	}

	///		Compute windows: Voronoi grid criterion.
	if (Options%VORO_GRID == 0)
	{
		CWFlag = ComputeWindows_VoroGrid(Molecule, true);
		if (CWFlag) { return true; }
	}

	///		Compute windows: Moving weighted alpha's criterion.
	if (Options%MV_ALPHA_W == 0)
	{
		CWFlag = ComputeWindows_ExploratoryAlphasW(Molecule, true);
		if (CWFlag) { return true; }
	}

	///		Compute windows: Moving unweighted alpha's criterion.
	if (Options%MV_ALPHA_UW == 0)
	{
		CWFlag = ComputeWindows_ExploratoryAlphas(Molecule, true);
		if (CWFlag) { return true; }
	}

	///		Compute windows: Generic mesh alpha shape criterion.
	if (Options%GEN_MESH == 0)
	{
		CWFlag = ComputeWindows_GenericMesh(Molecule, true);
		if (CWFlag) { return true; }
	}

	return false;

}



/********** FUNCTIONS FOR CALCULATION OF VORONOI ENTRY RAYS **********/

int IsEntryEdge(PointType Type1, PointType Type2)
{
	if (Type1 == internal_point && Type2 == boundary_point)
		return 0;
	if (Type2 == internal_point && Type1 == boundary_point)
		return 0;
	if (Type1 == internal_point && Type2 == external_point)
		return 0;
	if (Type2 == internal_point && Type1 == external_point)
		return 0;

	return -1;

}

int VoroEntryPaths( Complex *MoleculeC)
{

	int i;
	int NPaths = 0;

	for (i = 0; i < MoleculeC->GetNEdges(); i++)
	{
		PointType AuxType1, AuxType2;

		AuxType1 = MoleculeC->GetEdgeList()[i]->GetPSource()->GetType();
		AuxType2 = MoleculeC->GetEdgeList()[i]->GetPTarget()->GetType();

		if (IsEntryEdge(AuxType1, AuxType2) == 0)
			NPaths++;
	}

	return NPaths;

}





