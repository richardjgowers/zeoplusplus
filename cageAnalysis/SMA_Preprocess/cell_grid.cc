/*
 * cell_grid.cc
 *
 *  Created on: Jul 1, 2016
 *      Author: ismael.gomez
 */

#include "cell_grid.hh"

/*
 * 	Auxiliary function that computes the minimum distance among the distances to the atoms.
 */
double MinimumAtomDistance (	Point *GridPoint,
								Complex *Chemical)
{

	double MinDist = Chemical->GetPointList()[0]->DistanceToPoint(GridPoint);

	for (int i = 1; i < Chemical->GetNPoints(); i++)
	{
		MinDist = min(MinDist, Chemical->GetPointList()[i]->DistanceToPoint(GridPoint));
	}

	return MinDist;
}

/*
 * Generates a rectangular grid where all the points are separated from each other a distance proportional to the length of the box used to contain the molecule divided by the
 * 	number of points required by the user.
 */
void RectangularGrid (	double xlim[2],
						double ylim[2],
						double zlim[2],
						int NSteps,
						Complex *Chemical,
						list <Point *> *Grid)
{

	double 	dx = (xlim[1]-xlim[0])/NSteps,
			dy = (ylim[1]-ylim[0])/NSteps,
			dz = (zlim[1]-zlim[0])/NSteps;

	double BaseWeight;

	// Weight = max(dx, dy, dz)
	BaseWeight = max(dx, dy);
	BaseWeight = max(BaseWeight, dz);

	// Distance to the "farest neighbor" --> 3D Pitagoras.
	//Weight = sqrt(pow(dx, 2.) + pow(dy, 2.) + pow(dz, 2.));

	for (int i = 0; i < NSteps; i++)
	{
		for (int j = 0; j < NSteps; j++)
		{
			for (int k = 0; k < NSteps; k++)
			{
				double Coordinates[DIM_SMA] = {xlim[0]+i*dx, ylim[0]+j*dy, zlim[0]+k*dz};

				Point *GridPoint = new Point(0, Coordinates, 0.0, undefined_point);

				double Weight = min(BaseWeight, MinimumAtomDistance(GridPoint, Chemical) );

				GridPoint->SetRadius(Weight);

				Grid->push_back(GridPoint);
			}
		}
	}
}

/*
 * 	Given a molecule's chemical description and a list of points, computes which of those are internal. Information is stored at point's Type field.
 */
void MoleculeInternalPoints (	Complex *Chemical,
								list <Point *> *PointList)
{

	/////	GENERATE THE ALPHA SHAPE
	///		Use all points in the chemical description of the molecule.
	//		Use the optimal value for alpha: minimum alpha that generates 1 connected component.

	list<Weighted_point> UnweightedPList;

	Point *PAux;

	// Insert points into delaunay triangulaton object
	for (int i = 0; i < Chemical->GetNPoints(); i++)
	{
		PAux = Chemical->GetPointList()[i];

		// Insert point into the unweighted point list (weighted alpha-shapes)
		UnweightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), 0.0	  ) );

	}

	//Alpha_shape_3 MoleculeAlphaShape(DT);
	Alpha_shape_3 MoleculeAlphaShape(UnweightedPList.begin(), UnweightedPList.end(), 0, Alpha_shape_3::GENERAL);
	Alpha_iterator AlphaIt = MoleculeAlphaShape.find_optimal_alpha(1);
	//Alpha_iterator AlphaIt = MoleculeAlphaShape.alpha_end();

	//AlphaIt--;

	double Alpha = *(AlphaIt);

	Alpha = 1.1*Alpha;

	/////	CLASSIFY THE POINTS
	///		Check if they are inside any cell of the alpha shape.
	//			- Get all the alpha cells
	//			- Check point by point for all the cells
	//		In that case, change the Type field to internal_point.

	/// Retrieve all the internal cells of the alpha shape
	list <Cell_handle> CellList;
	list <Cell_handle>::iterator CellIt;

	MoleculeAlphaShape.get_alpha_shape_cells(back_inserter(CellList), Alpha_shape_3::INTERIOR, Alpha);

	///	Run for all the points in list
	list <Point *>::iterator PointIt;

	for (PointIt = PointList->begin(); PointIt != PointList->end(); PointIt++)
	{
		Point *P0 = *(PointIt);

		for (CellIt = CellList.begin(); CellIt != CellList.end(); CellIt++)
		{
			Cell_handle Cell = *(CellIt);

			double Coordinates[DIM_SMA] = {P0->GetCoordinate(X_COORD), P0->GetCoordinate(Y_COORD), P0->GetCoordinate(Z_COORD)};

			if ( PointInsideCell(Coordinates, Cell) )
			{
				P0->SetType(internal_point);
				break;
			}
		}


	}

}

/*
 * Compute boundaries.
 *
 * TODO - ERASE THIS!!!!!! -> Duplicated in voro_interface.cc (CalculateBoundaries)
 */
void CalculateGridBoundaries (	double *BoundRef,
								int Coordinate,
								Complex *Chemical)
{

	double Max, Min;

	Max = Min = Chemical->GetPointList()[0]->GetCoordinate(Coordinate);

	for (int i = 1; i < Chemical->GetNPoints(); i++)
	{
		if (Chemical->GetPointList()[i]->GetCoordinate(Coordinate) > Max) {
			Max = Chemical->GetPointList()[i]->GetCoordinate(Coordinate);
		}
		if (Chemical->GetPointList()[i]->GetCoordinate(Coordinate) < Min) {
			Min = Chemical->GetPointList()[i]->GetCoordinate(Coordinate);
		}
	}

	double Diff = 0.1*abs(Max-Min);

	// Store the boundaries: Maximum distances + 10%
	BoundRef[0] = Min-Diff;
	BoundRef[1] = Max+Diff;

}


/*
 * 	Generates a grid with points separated an equal distance, then computes which of them are internal.
 */
void InternalPointsGrid (	MoleculeInfo *Molecule,
							list <Point *> *PointList)
{

	double xlim[2], ylim[2], zlim[2];

	CalculateGridBoundaries(xlim, X_COORD, Molecule->GetChemical());
	CalculateGridBoundaries(ylim, Y_COORD, Molecule->GetChemical());
	CalculateGridBoundaries(zlim, Z_COORD, Molecule->GetChemical());

//	cout << "xlim = (" << xlim[0] << ", " << xlim[1] << ")" << endl;
//	cout << "ylim = (" << ylim[0] << ", " << ylim[1] << ")" << endl;
//	cout << "zlim = (" << zlim[0] << ", " << zlim[1] << ")" << endl;


//	int NPoints = Molecule->GetChemical()->GetNPoints();

	//RectangularGrid(xlim, ylim, zlim, NPoints, PointList);
	RectangularGrid(xlim, ylim, zlim, 10, Molecule->GetChemical(), PointList);

	MoleculeInternalPoints(Molecule->GetChemical(), PointList);


}


/*
 *
 */

void CellCenterGrid (	Complex *Chemical,
						list <Point *> *PointList)
{

	/////	GENERATE THE ALPHA SHAPE
	///		Use all points in the chemical description of the molecule.
	//		Use the optimal value for alpha: minimum alpha that generates 1 connected component.

	list<Weighted_point> UnweightedPList;

	Point *PAux;

	// Insert points into delaunay triangulaton object
	for (int i = 0; i < Chemical->GetNPoints(); i++)
	{
		PAux = Chemical->GetPointList()[i];

		// Insert point into the unweighted point list (weighted alpha-shapes)
		UnweightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), 0.0	  ) );

	}

	//Alpha_shape_3 MoleculeAlphaShape(DT);
	Alpha_shape_3 MoleculeAlphaShape(UnweightedPList.begin(), UnweightedPList.end(), 0, Alpha_shape_3::GENERAL);
	Alpha_iterator AlphaIt = MoleculeAlphaShape.find_optimal_alpha(1);

	double Alpha = *(AlphaIt);

	Alpha = 1.5*Alpha;

	/////	CLASSIFY THE POINTS
	///		Check if they are inside any cell of the alpha shape.
	//			- Get all the alpha cells
	//			- Check point by point for all the cells
	//		In that case, change the Type field to internal_point.

	/// Retrieve all the internal cells of the alpha shape
	list <Cell_handle> CellList;
	list <Cell_handle>::iterator CellIt;

	MoleculeAlphaShape.get_alpha_shape_cells(back_inserter(CellList), Alpha_shape_3::INTERIOR, Alpha);

	for (CellIt = CellList.begin(); CellIt != CellList.end(); CellIt++)
	{
		Cell_handle Cell = *(CellIt);

		double Pseudocenter[DIM_SMA];

		double Weight = CellPseudocenter(Cell, Pseudocenter);

		Point *GridPoint = new Point(0, Pseudocenter, Weight, internal_point);

		PointList->push_back(GridPoint);

	}


}


/*	InternalNodeSphereCheck
 * 		Checks if a Voronoi node is internal by building a set of points around it (in spherical shape) at a distance equal to 0.95 times the radius of that Voronoi node.
 * 		Points around should be inside the alpha-shape.
 *
 */
bool InternalNodeSphereCheck ( 	MoleculeInfo *Molecule,
								Point *VoroInternalPoint)
{

	//double VIPCoords[DIM_SMA] = {VoroInternalPoint->GetCoordinate(X_COORD), VoroInternalPoint->GetCoordinate(Y_COORD), VoroInternalPoint->GetCoordinate(Z_COORD)  };
	double VIPCoords[DIM_SMA];
	VoroInternalPoint->GetCoordinates(VIPCoords);

	///		Build the sphere

	double Radius = 0.95*VoroInternalPoint->GetRadius();

	list <Point *> *SpherePoints = new list <Point *>();

	///		Points in the axis
	for (int i = 0; i < 2; i++)
	{
		Point *SpherePointX = new Point();

		SpherePointX->SetCoordinate(X_COORD, VIPCoords[X_COORD] + pow(-1, i)*Radius);
		SpherePointX->SetCoordinate(Y_COORD, VIPCoords[Y_COORD]);
		SpherePointX->SetCoordinate(Z_COORD, VIPCoords[Z_COORD]);

		SpherePoints->push_back(SpherePointX);

		Point *SpherePointY = new Point();

		SpherePointY->SetCoordinate(X_COORD, VIPCoords[X_COORD]);
		SpherePointY->SetCoordinate(Y_COORD, VIPCoords[Y_COORD] + pow(-1, i)*Radius);
		SpherePointY->SetCoordinate(Z_COORD, VIPCoords[Z_COORD]);

		SpherePoints->push_back(SpherePointY);

		Point *SpherePointZ = new Point();

		SpherePointZ->SetCoordinate(X_COORD, VIPCoords[X_COORD]);
		SpherePointZ->SetCoordinate(Y_COORD, VIPCoords[Y_COORD]);
		SpherePointZ->SetCoordinate(Z_COORD, VIPCoords[Z_COORD] + pow(-1, i)*Radius);

		SpherePoints->push_back(SpherePointZ);
	}


	///		Check if the generated points are internal

	MoleculeInternalPoints ( Molecule->GetChemical(), SpherePoints);

	list <Point *>::iterator SPIt;

	for (SPIt = SpherePoints->begin(); SPIt != SpherePoints->end(); SPIt++)
	{
		Point *SP = *SPIt;

		///		If any point of the sphere is not contained in the molecule, return false
		if (SP->GetType() != internal_point)
		{
			return false;
		}

	}

	///		All the sphere points are inside the molecule (according to the alpha-shape criterion)
	return true;

}

/*	SphericalGridInternalCorrection
 * 		Checks biggest internal cell and adjusts its type using the spherical grid criterion. If the cell is reclassified as external, the process is repeated.
 */
void SphericalGridInternalCorrection(MoleculeInfo *Molecule)
{

	int 	Count = 0,
			MaxCount = Molecule->GetVoronoiGraph()->GetNPoints();

	///		Until a legal internal point is found, apply the correction
	bool MaxIntSphereLegal = false;

	while (!MaxIntSphereLegal && Count < MaxCount)
	{
		///		Get the point with maximum internal volume.
		Point *MVPoint = Molecule->GetGreatestVoroNode(internal_point);

		if (MVPoint == NULL)
		{
			return;
		}

		///		Check if is legal internal point.
		MaxIntSphereLegal = InternalNodeSphereCheck(Molecule, MVPoint);

		if (!MaxIntSphereLegal)
		{

			cout << "Internal cell remodeled "<< endl;

			MVPoint->SetType(external_point);

			//MVPoint = Molecule->GetGreatestVoroNode(internal_point);
			Count++;


		}
	}

}
