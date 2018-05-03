//	Software for single molecule analysis
//
//	internal_cells.cc
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : February 15th 2016

#include "internal_cells.hh"

/*
 * 	Checks if the given edge (that must be an entry edge, i.e., have one internal node and one external/boundary node) fits the restriction criterion to be an actual indicator of a window.
 */
double RestrictionCriterion (	Edge *TEdge)
{


	///		Compute the restrictions and apply the criterion

	Point 	*PSource = TEdge->GetPSource(),
			*PTarget = TEdge->GetPTarget();

	double Inner, Outer, Mid;

	///		Store the outer as the less restrictive
	if (PSource->GetRadius() > PTarget->GetRadius())
	{
		Outer = PSource->GetRadius();
		Inner = PTarget->GetRadius();
	}
	else
	{
		Inner = PSource->GetRadius();
		Outer = PTarget->GetRadius();
	}

	Mid = TEdge->GetChannelRadius();

	//cout << Outer << " " << Mid << " " << Inner << endl;

	///		Compute the restriction

	if (Inner >= Mid)
	{
		//cout << "Restriction found: " << Mid << endl;
		return Mid;
	}

	return -1.0;

}

/*
 * 	Computes alpha for the alpha shape, as the highest restricting value for the molecule
 */
double ComputeAlphaRestrictionCriterion(Complex *Voronoi)
{
	double Alpha = 0;

	for (int i = 0; i < Voronoi->GetNEdges(); i++)
	{
		Edge *VoroEdge = Voronoi->GetEdgeList()[i];

		double Restriction = 0.;

		Restriction = RestrictionCriterion(VoroEdge);

		if (Alpha < Restriction)
		{
			Alpha = Restriction;
		}

	}

	return 2.*Alpha;

}

/*
 *
 */
bool RayCrossFacet (	double RayEnd1[DIM_SMA],
						double RayEnd2[DIM_SMA],
						Simplex *TFacet)
{

	double Point1[DIM_SMA], Point2[DIM_SMA], Point3[DIM_SMA];

	TFacet->ReadPoint(0)->PointCoordinates(Point1);
	TFacet->ReadPoint(1)->PointCoordinates(Point2);
	TFacet->ReadPoint(2)->PointCoordinates(Point3);

	return RayCrossesTriangle(RayEnd1, RayEnd2, Point1, Point2, Point3);

}


bool EdgeCrossWindow(	Edge *TEdge,
						Window *TWindow)
{

	/////	CHECK IF EDGE CROSSES WINDOW
	///		For each facet in window:
	//		Check if the edge crosses it, then return true.
	//		If none of the facets are crossed, return false.

	list <Simplex *>::iterator FacetIt;

	for (FacetIt = TWindow->GetSimplices()->begin(); FacetIt != TWindow->GetSimplices()->end(); FacetIt++)
	{
		double RayEnd1[DIM_SMA], RayEnd2[DIM_SMA];

		TEdge->GetPSource()->PointCoordinates(RayEnd1);
		TEdge->GetPTarget()->PointCoordinates(RayEnd2);

		Simplex *TFacet = *(FacetIt);

		if (RayCrossFacet(RayEnd1, RayEnd2, TFacet))
		{
			return true;
		}
	}

	return false;


}

/*
 * Checks if the edge is an entry (using point types as criterion)
 */
bool EntryEdgeTypeCriterion (Edge *TEdge)
{

	PointType SType = TEdge->GetPSource()->GetType();
	PointType TType = TEdge->GetPTarget()->GetType();

	// TODO: Add error control for the case where points are undefined (shouldn't be any)

	// Both types are the same -> Is not an entry
	if (SType == TType) {
		return false;
	}

	// The first end is internal -> Edge will be an entry iff the other end is either external or boundary point
	if (SType == internal_point) {
		if (TType == external_point || TType == boundary_point) {
			return true;
		}
		else {
			return false;
		}
	}

	// The first end is external/boundary -> Edge will be an entry iff the other end is internal
	if (SType == external_point || SType == boundary_point) {
		if (TType == internal_point) {
			return true;
		}
		else {
			return false;
		}
	}

	return false;

}


void InternalPointReclassify (Edge *TEdge)
{

	if (TEdge->GetPSource()->GetType() == internal_point)
	{
		TEdge->GetPSource()->SetType(hybrid_point);
	}

	if (TEdge->GetPTarget()->GetType() == internal_point)
	{
		TEdge->GetPTarget()->SetType(hybrid_point);
	}

}


/* WindowCriterionCheck
 * 	Checks if the internal/external classification was correctly applied using the window criterion: any external node should be connected to an internal one only by means of an edge that crosses a window.
 *
 * Input:
 * 	Molecule	- Molecule information, including Voronoi complex.
 * Output:
 * 	True if the Voronoi nodes fulfill the "window criterion" explained above.
 */
bool WindowCriterionCheck (MoleculeInfo *Molecule)
{

	/////	APPLY WINDOW CROSSING CRITERION
	///		For each edge in the Voronoi graph:
	//		(1) Check if the edge crosses any window
	//		(2) - If the 2 points are classified as different (internal/external or internal/boundary) and the edge does not cross a window, return false.
	//			- If the 2 points are classified as equal (boundary/external or same types) and the edge does cross a window, return false.
	//		If the loop ends without finding any inconsistency, return true.

	Complex *VoronoiGraph = Molecule->GetVoronoiGraph();

	bool Reclassified = false;

	for (int i = 0; i < VoronoiGraph->GetNEdges(); i++)
	{
		Edge *EAux = VoronoiGraph->GetEdgeList()[i];

		bool FlagW = false;

		// (1)
		if (EntryEdgeTypeCriterion(EAux))
		{
			for (int j = 0; j < Molecule->GetNWindows(); j++)
			{
				Window *WAux = Molecule->GetWindows()[j];

				if (EdgeCrossWindow(EAux, WAux))
				{
					FlagW = true;
					break;
				}
				else
				{
					InternalPointReclassify(EAux);
					Reclassified = true;
				}
			}
		}

		// (2)
		//bool FlagT = EntryEdgeTypeCriterion(EAux);

		/*if (FlagW != FlagT)
		{
			cout << "Window criterion violated: " << endl;
			EAux->GetPSource()->PrintPoint();
			EAux->GetPTarget()->PrintPoint();
			if (FlagT) cout << "Edge connects internal and external." << endl;
			if (FlagW) cout << "Window crossed." << endl;

			InternalPointReclassify(EAux);
			Reclassified = true;
			//return false;
		}*/

	}

	return Reclassified;

}

/* InternalCellsAlphaCriterion
 * 	Computes which are the internal Voronoi nodes of the molecule, based on the following criterion: a node is internal if its contained inside an internal cell of the alpha shape of the molecule.
 *
 * 
 */

bool ComputeInternalCellsAlpha(		Complex	*MolVoronoi,
									Alpha_shape_3 &AlphaShape,
									double Alpha)
{


	/////	COMPUTE VORONOI NODES' TYPE
	///		Check, for every node, if its contained in any cell of the alpha shape.
	//		The ones contained are internal nodes, the ones not contained, external.

	list <Cell_handle> CellList;
	list <Cell_handle>::const_iterator CellIt;

	// Retrieve all the internal cells of the alpha shape	
	AlphaShape.get_alpha_shape_cells(back_inserter(CellList), Alpha_shape_3::INTERIOR, Alpha);

	// For every voronoi node and every cell in the alpha shape, check if the voronoi node is inside that cell
	for (int i = 0; i < MolVoronoi->GetNPoints(); i++)
	{
		if (MolVoronoi->GetPointList()[i]->GetType() == undefined_point)
		{
			
			// Set the point as external point (if its internal, it will be modified later)
			MolVoronoi->GetPointList()[i]->SetType(external_point);

			// Get the coordinates of the point
			double VoroNode[DIM_SMA];
			MolVoronoi->GetPointList()[i]->PointCoordinates(VoroNode);
			
			// Iterate over all the cells, checking if the node is inside any of them
			for (CellIt = CellList.begin(); CellIt != CellList.end(); CellIt++)
			{
				if ( PointInsideCell(VoroNode, *CellIt) )
				{
					MolVoronoi->GetPointList()[i]->SetType(internal_point);
					break;
				}
			} // END FOR (Alpha-shape cell list)
		}

	} // END FOR (Voronoi nodes list)
	

	return true;

}




/* MoleculeInternalCells
 * 	Computes internal cells of the molecule using the alpha-shape cells criterion, then checks if the computation is right using the window-edge criterion.
 *
 *
 */

bool MoleculeInternalCells (	MoleculeInfo *Molecule)
{

	/////	COMPUTE THE ALPHA-SHAPE OF THE MOLECULE
	///		First compute the Delaunay Triangulation (cgal-based).
	//		Alpha-shape computed for every possible alpha.
	//		Alpha calculated in such a way that the molecule is connex.


	list<Weighted_point> UnweightedPList;

	Point *PAux;

	//////////

	// Insert points into delaunay triangulaton object
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		PAux = Molecule->GetChemical()->GetPointList()[i];


		// Check that the point is of desired type
		if (PAux->GetType() == atom_point)
		{
			// Insert point into the unweighted point list (weighted alpha-shapes)
			UnweightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), 0.0	  ) );

		}


	}

	//Alpha_shape_3 MoleculeAlphaShape(DT);
	Alpha_shape_3 MoleculeAlphaShape(UnweightedPList.begin(), UnweightedPList.end(), 0, Alpha_shape_3::GENERAL);

	double Alpha;

	Alpha_iterator 	opt = MoleculeAlphaShape.find_optimal_alpha(1),
					end = MoleculeAlphaShape.alpha_end();

	// If the alpha iterator is placed at the end, no optimal alpha was found --> The method failed to identify internal cells (there are none)
	if (opt == MoleculeAlphaShape.alpha_end())
	{
		//cout << "Ups" << endl;
		return false;
	}

	///		Fix alpha to optimal

	Alpha = *opt;
	Alpha += 0.1*Alpha;

	//Alpha = ComputeAlphaRestrictionCriterion(Molecule->GetVoronoiGraph());

	//cout << "Alpha-pre: " << Alpha << endl;

	double NewAlpha = AlphaSelectionAllBonds(Molecule->GetChemical(), MoleculeAlphaShape);

	if (NewAlpha > Alpha)
	{
		//cout << "Warning: Alpha selected did not include all bonds" << endl;
		Alpha = NewAlpha;
	}

	//Alpha = NewAlpha;

	//cout << "Alpha: " << Alpha << endl;


	//Alpha *= 5;

	/////	CLASSIFY POINTS OF VORONOI TESSELLATION
	///		Compute internal points using alpha-cells criterion, then use window-cross criterion to check if the points were correctly classified.
	//		(1) Apply alpha-cells criterion: a point is inside if it's inside an alpha-cell (includes border).
	//		(2) Clusterize complex: after classification, extremely close nodes are clusterized using distance and node-type criteria.
	//		(3)	Window-cross criterion: any edge connecting inside with outside should cross and edge.

	Complex	*MolVoronoi = Molecule->GetVoronoiGraph();

	// (1)

	end--;
	Alpha = *end;

	ComputeInternalCellsAlpha(MolVoronoi, MoleculeAlphaShape, Alpha);

	// (2)
	Complex *Clusterized = ClusterizeComplex(MolVoronoi, Molecule->GetChemical(), ATOMIC_CLUSTER_THRESHOLD);

	delete Molecule->GetVoronoiGraph();

	Molecule->SetVoronoiGraph(Clusterized);

	//Molecule->SetVoronoiGraph(MolVoronoi);

	// (3)
	//if (WindowCriterionCheck(Molecule))
	//{
		// TODO Send the error message only after several iterations
	//	cout << "Error: computing internal cells failed!!" << endl;
	//	return false;

	//}

	/////	COMPUTE CONVEX HULL INTERNAL CELLS
	///		Same procedure, but alpha is set to the maximum possible value.

	end--;
	Alpha = *end;

	//cout << "Alha end: " << Alpha << endl;

	Complex	*VoroCH = Molecule->GetVoroConvexHull();

	// (1)
	ComputeInternalCellsAlpha(VoroCH, MoleculeAlphaShape, Alpha);

	// (2)
	Complex *VoroCHCluster = ClusterizeComplex(VoroCH, Molecule->GetChemical(), ATOMIC_CLUSTER_THRESHOLD);

	delete Molecule->GetVoroConvexHull();

	Molecule->SetVoroConvexHull(VoroCHCluster);


	return true;
}

/*
 *
 *
 */

bool MoleculeConvexHullCells (MoleculeInfo *Molecule)
{

	/////	COMPUTE CONVEX HULL INTERNAL CELLS
	///		Copy Voronoi graph after internal cells were computed
	//		Compute internal cells for convex hull.
	//		Compare both graphs and mark as hybrid cells that are inside the convex hull but not inside the molecule.


	return true;
}


/////////////////////////////////////////////////////////////
/////	FUNCTIONS FOR GENERAL POINT-INSIDE-CELL STUDY	/////
/////////////////////////////////////////////////////////////

/*
 *	Check if the given point is inside an alpha-shape
 */
bool IsInternalPoint(	Point *TargetPoint,
						Alpha_shape_3 &AlphaShape,
						double Alpha)
{

	list <Cell_handle> CellList;
	list <Cell_handle>::const_iterator CellIt;

	// Retrieve all the internal cells of the alpha shape
	AlphaShape.get_alpha_shape_cells(back_inserter(CellList), Alpha_shape_3::INTERIOR, Alpha);

	double PointCoords[DIM_SMA];
	TargetPoint->GetCoordinates(PointCoords);

	// Iterate over all the cells, checking if the node is inside any of them
	for (CellIt = CellList.begin(); CellIt != CellList.end(); CellIt++)
	{
		if ( PointInsideCell(PointCoords, *CellIt) )
		{
			return true;
		}
	} // END FOR (Alpha-shape cell list)

	return false;

}

/*
 *	Checks if the points in given list are inside an alpha-shape.
 *
 *	Points are reclassified as internal/external in depends on this.
 *	If FastAnswer is true, no reclassification happens. Instead, true is returned if at least one point is inside the chemical.
 *
 */
bool AreInternalPoints(	list <Point *> &PointList,
						Complex *Chemical,
						bool FastAnswer)
{

	/////	COMPUTE THE ALPHA-SHAPE OF THE MOLECULE
	///		First compute the Delaunay Triangulation (cgal-based).
	//		Alpha-shape computed for every possible alpha.
	//		Alpha calculated in such a way that the molecule is connex.


	list<Weighted_point> UnweightedPList;

	Point *PAux;

	//////////

	// Insert points into delaunay triangulaton object
	for (int i = 0; i < Chemical->GetNPoints(); i++)
	{
		PAux = Chemical->GetPointList()[i];

		// Check that the point is of desired type
		if (PAux->GetType() == atom_point)
		{
			// Insert point into the unweighted point list (weighted alpha-shapes)
			UnweightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), 0.0	  ) );

		}


	}

	//Alpha_shape_3 MoleculeAlphaShape(DT);
	Alpha_shape_3 MoleculeAlphaShape(UnweightedPList.begin(), UnweightedPList.end(), 0, Alpha_shape_3::GENERAL);

	double Alpha;

	Alpha_iterator 	opt = MoleculeAlphaShape.find_optimal_alpha(1);
					//end = MoleculeAlphaShape.alpha_end();

	// If the alpha iterator is placed at the end, no optimal alpha was found --> The method failed to identify internal cells (there are none)
	if (opt == MoleculeAlphaShape.alpha_end())
	{
		//cout << "Ups" << endl;
		return false;
	}

	///		Fix alpha to optimal

	Alpha = *opt;
	Alpha += 0.1*Alpha;

	/////	DECIDE WHETHER OR NOT POINTS IN LIST ARE INTERNAL
	///		Check point by point.
	//		Internal points are reclassified as internal_point, non-internal are reclassified as external_point

	list <Point *>::iterator PIt;

	bool Answer = false;

	for (PIt = PointList.begin(); PIt != PointList.end(); PIt++)
	{
		Point *TargetPoint = *PIt;

		if (IsInternalPoint(TargetPoint, MoleculeAlphaShape, Alpha) )
		{
			if (FastAnswer)
			{
				return true;
			}
			else
			{
				Answer = true;
				TargetPoint->SetType(internal_point);
			}

		}
		else
		{
			if (!FastAnswer)
			{
				TargetPoint->SetType(external_point);
			}

		}

	}

	return Answer;

}


