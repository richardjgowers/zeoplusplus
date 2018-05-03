/*
 * alpha_shape.cc
 *
 *  Created on: Apr 4, 2016
 *      Author: ismael.gomez
 */

#include "alpha_shape.hh"

////// FUNCTIONS FOR STORING ALPHA SHAPE INFO INTO COMPLEX

void InsertAlphaEdges (	Alpha_shape_3 &MoleculeAlphaShape,
						Classification_type CType,
						double Alpha,
						Complex *CAlphaShape)
{


	list<CG_edge> EdgeList;
	list<CG_edge>::const_iterator EIt;

	// Get the list of edges of the alpha shape
	MoleculeAlphaShape.get_alpha_shape_edges(back_inserter(EdgeList), CType, Alpha);

	cout << "Alpha shape edges: " << EdgeList.size() << endl;

	PointType PT;

	switch(CType) {

		case Alpha_shape_3::INTERIOR:
			PT = internal_point;
			break;
		case Alpha_shape_3::EXTERIOR:
			PT = external_point;
			break;
		case Alpha_shape_3::REGULAR:
			PT = boundary_point;
			break;
		case Alpha_shape_3::SINGULAR:
			PT = hybrid_point;
			break;
		default:
			break;
	}

	// Go though the list of edges, printing vertices belonging to each cell

	for (EIt = EdgeList.begin(); EIt != EdgeList.end(); EIt++)
	{

		CG_edge TargetEdge = *EIt;

		int Vert1 = TargetEdge.second, Vert2 = TargetEdge.third;

		double 	CoordinatesS[DIM_SMA] = {TargetEdge.first->vertex(Vert1)->point().x(), TargetEdge.first->vertex(Vert1)->point().y(), TargetEdge.first->vertex(Vert1)->point().z()},
				CoordinatesT[DIM_SMA] = {TargetEdge.first->vertex(Vert2)->point().x(), TargetEdge.first->vertex(Vert2)->point().y(), TargetEdge.first->vertex(Vert2)->point().z()};

		Point 	*PSource = CAlphaShape->GetPointByCoordinates(CoordinatesS),
				*PTarget = CAlphaShape->GetPointByCoordinates(CoordinatesT);

		if (PSource != NULL && PTarget != NULL)
		{
			PSource->SetType(PT);
			PTarget->SetType(PT);
			CAlphaShape->InsertEdgeByValue(PSource, PTarget);
		}
	}

}


void InsertAlphaVertices (	Alpha_shape_3 &MoleculeAlphaShape,
							Classification_type CType,
							double Alpha,
							Complex *CAlphaShape)
{


	list<Vertex_handle> VertexList;
	list<Vertex_handle>::const_iterator VertexListIt;

	// Get the list of vertices of the alpha shape
	MoleculeAlphaShape.get_alpha_shape_vertices(back_inserter(VertexList), CType, Alpha);

	// Go though the list of facets, printing vertices belonging to each cell

	for (VertexListIt = VertexList.begin(); VertexListIt != VertexList.end(); VertexListIt++)
	{

		Vertex_handle TargetF = *VertexListIt;

		double Coordinates[DIM_SMA] = {TargetF->point().x(), TargetF->point().y(), TargetF->point().z()};

		//CAlphaShape->InsertPointByValue(Coordinates, 0., PT);
		CAlphaShape->InsertPointByValue(Coordinates);

	 }

}

void InsertFromAlphaCells(	Alpha_shape_3 &MoleculeAlphaShape,
							Classification_type CType,
							double Alpha,
							Complex *CAlphaShape)
{

	list<Cell_handle> CellList;
	list<Cell_handle>::const_iterator CellIt;

	MoleculeAlphaShape.get_alpha_shape_cells(back_inserter(CellList), CType, Alpha);

	for (CellIt = CellList.begin(); CellIt != CellList.end(); CellIt++)
	{
		Cell_handle TCell = *CellIt;

		///		INSERT VERTICES

		for (int i = 0; i < CELL_N_V; i++)
		{
			double Coordinates[DIM_SMA] = {TCell->vertex(i)->point().x(), TCell->vertex(i)->point().y(), TCell->vertex(i)->point().z()};

			if (CAlphaShape->GetPointByCoordinates(Coordinates) == NULL)
			{
				CAlphaShape->InsertPointByValue(Coordinates);
			}

		}	// End for (cell nodes)


		///		INSERT EDGES

		for (int i = 0; i < CELL_N_V-1; i++)
		{
			double Coordinates[DIM_SMA] = {TCell->vertex(i)->point().x(), TCell->vertex(i)->point().y(), TCell->vertex(i)->point().z()};

			for (int j = i+1; j < CELL_N_V; j++)
			{
				double Coordinates2[DIM_SMA] = {TCell->vertex(j)->point().x(), TCell->vertex(j)->point().y(), TCell->vertex(j)->point().z()};

				//cout << "Printing " << i << " " << j << endl;

				if (CAlphaShape->GetEdgeByCoordinates(Coordinates, Coordinates2) == NULL)
				{
					CAlphaShape->InsertEdgeByValue(CAlphaShape->GetPointByCoordinates(Coordinates), CAlphaShape->GetPointByCoordinates(Coordinates2));

					//if (CAlphaShape->GetPointByCoordinates(Coordinates) != NULL && CAlphaShape->GetPointByCoordinates(Coordinates2) != NULL)
					//{}
				}
			}	// End for (cell nodes [internal])
		}	// End for (cell nodes [external])


	}	// End for (cells)

}


void InsertFromAlphaFacets(	Alpha_shape_3 &MoleculeAlphaShape,
							Classification_type CType,
							double Alpha,
							Complex *CAlphaShape)
{

	// List of facets for containing the total of facets of the alpha shape
	list<Facet> FacetList;
	list<Facet>::const_iterator FacetIt;

	// Get the list of facets of the alpha shape
	MoleculeAlphaShape.get_alpha_shape_facets(back_inserter(FacetList), CType, Alpha);

	cout << "ALPHA SHAPE TOTAL FACETS: " << FacetList.size() << endl;

	for (FacetIt = FacetList.begin(); FacetIt != FacetList.end(); FacetIt++)
	{
		Facet FacetT = *FacetIt;

		int OppositeVertId = FacetT.second;

		int Indices[FACET_N_V];
		int j = 0;

		////	Get the indices of the facet points
		for (int i = 0; i < CELL_N_V; i++)
		{
			if (i != OppositeVertId)
			{
				Indices[j] = i;
				j++;
			}

		}

		////	Insert points and edges
		double 	Point1[DIM_SMA] = {FacetT.first->vertex(Indices[0])->point().x(), FacetT.first->vertex(Indices[0])->point().y(), FacetT.first->vertex(Indices[0])->point().z()},
				Point2[DIM_SMA] = {FacetT.first->vertex(Indices[1])->point().x(), FacetT.first->vertex(Indices[1])->point().y(), FacetT.first->vertex(Indices[1])->point().z()},
				Point3[DIM_SMA] = {FacetT.first->vertex(Indices[2])->point().x(), FacetT.first->vertex(Indices[2])->point().y(), FacetT.first->vertex(Indices[2])->point().z()};

		CAlphaShape->InsertPointByValue(Point1);
		CAlphaShape->InsertPointByValue(Point2);
		CAlphaShape->InsertPointByValue(Point3);

		CAlphaShape->InsertEdgeByValue(CAlphaShape->GetPointByCoordinates(Point1), CAlphaShape->GetPointByCoordinates(Point2));
		CAlphaShape->InsertEdgeByValue(CAlphaShape->GetPointByCoordinates(Point1), CAlphaShape->GetPointByCoordinates(Point3));
		CAlphaShape->InsertEdgeByValue(CAlphaShape->GetPointByCoordinates(Point2), CAlphaShape->GetPointByCoordinates(Point3));

	}

}


double ComputeWeight( 	Complex *Chemical,
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

double ComputeWeightMax ( 	Complex *Chemical,
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


/*	AlphaSelectionAllBonds
 * 		Given the alpha-shape and the molecule's chemical description, returns the minimum value of alpha that guarantees all bonds to be part of the alpha shape.
 */

/*
 *
 */

bool EdgePartOfCell (	Edge *TargetEdge,
						Cell_handle &Cell)
{


	//double 	CoordsPSource[DIM_SMA],
	//		CoordsPTarget[DIM_SMA];

	//TargetEdge->GetPSource()->GetCoordinates(CoordsPSource);
	//TargetEdge->GetPTarget()->GetCoordinates(CoordsPTarget);

	bool 	PSFlag = false,
			PTFlag = false;

	for (int i = 0; i < CELL_N_V; i++)
	{
		double CellCoords[DIM_SMA] = {Cell->vertex(i)->point().x(), Cell->vertex(i)->point().y(), Cell->vertex(i)->point().z()};

		if (!PSFlag) {
			if ( TargetEdge->GetPSource()->SameCoordinates(CellCoords)) {
				PSFlag = true;
			}
		}

		if (!PTFlag) {
			if ( TargetEdge->GetPTarget()->SameCoordinates(CellCoords)) {
				PTFlag = true;
			}
		}

	}

	return PSFlag && PTFlag;

}

/*	AlphaSelectionAllBonds
 * 		Given the alpha-shape and the molecule's chemical description, returns the minimum value of alpha that guarantees all bonds to be part of the alpha shape.
 */
double AlphaSelectionAllBonds( 	Complex *Chemical,
								Alpha_shape_3 &AlphaShape)
{


	/////	BONDS IN ALPHA-SHAPE CHECK
	///		Set alpha to maximum possible value and create alpha-shape.
	//		For the given alpha-shape, check that every bond is in it:
	//			Retrieve non-external cell.
	//			Check that there is a cell that has the bond as edge.
	//		When a bond is not found as being part of the alpha-shape, return to the prior alpha.

	bool AllBondsIn = true;

	Alpha_iterator AlphaIt = AlphaShape.alpha_end();

	AlphaIt--;

	double 	Alpha = *AlphaIt,
			AlphaPrev = Alpha;

	//cout << "Alpha end:" << Alpha << endl;

	while (AllBondsIn && Alpha > 0.0)
	{

		list <Cell_handle> CellList;
		list <Cell_handle>::iterator CellIt;

		///	 Retrieve all the internal cells of the alpha shape
		AlphaShape.get_alpha_shape_cells(back_inserter(CellList), Alpha_shape_3::INTERIOR, Alpha);

		AllBondsIn = true;

		///		Run over molecule's bonds
		for (int i = 0; i < Chemical->GetNEdges(); i++)
		{

			///		Retrieve the edge
			Edge *Bond = Chemical->GetEdgeList()[i];

			bool BondIn = false;

			///		Check cells
			for (CellIt = CellList.begin(); CellIt != CellList.end(); CellIt++)
			{
				Cell_handle Cell = *CellIt;

				if (EdgePartOfCell(Bond, Cell))
				{
					BondIn = true;
					break;
				}
			}	// End for (Cells)

			if (!BondIn)
			{
				AllBondsIn = false;
				return AlphaPrev;
				//return Alpha;
				//break;
			}

		}	// End for (bonds)

		///		If the bond was found to be in, continue, otherwise finalize
		if (AllBondsIn)
		{
			AlphaPrev = Alpha;
			--AlphaIt;
			Alpha = *AlphaIt;

			//cout << "alpha prev:" << AlphaPrev << " - alpha " << Alpha << endl;
			//AllBondsIn = false;
		}
		else
		{
			return AlphaPrev;
		}


	} // End while



	cout << "Ups" << endl;

	return AlphaPrev;

}

/*
 *
 *
 *
 */
Complex *AlphaShapeToComplex (	MoleculeInfo *Molecule)
{

	//Delaunay_triangulation DT;

	list<Weighted_point> WeightedPList;

	Point *PAux;

	// Insert points into Delaunay triangulation object
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		PAux = Molecule->GetChemical()->GetPointList()[i];


		// Check that the point is of desired type
		if (PAux->GetType() == atom_point)
		{
			CG_point AuxCGP(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2));
			//DT.insert(AuxCGP);

			//double Weight = ComputeWeight(Molecule->GetChemical(), PAux)/2.;
			//double Weight = pow(ComputeWeight(Molecule->GetChemical(), PAux), 2.)/2.;
			//double Weight = ComputeWeightMax(Molecule->GetChemical(), PAux);

			//cout << "Weight chosen for point " << Weight << endl;

			//Weight = pow(Weight, 2.);
			double Weight = 0.;

			// Insert point into the weighted point list (weighted alpha-shapes)
			WeightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );

		}
	}

	//Alpha_shape_3 MoleculeAlphaShape(DT);
	Alpha_shape_3 MoleculeAlphaShape(WeightedPList.begin(), WeightedPList.end(), 0, Alpha_shape_3::GENERAL);
	//Alpha_shape_3 MoleculeAlphaShape(WeightedPList.begin(), WeightedPList.end(), Alpha_shape_3::REGULARIZED);
	//Alpha_shape_3 MoleculeAlphaShape(WeightedPList.begin(), WeightedPList.end(), Alpha_shape_3::GENERAL);

	/////	INSERT ALPHA SHAPE POINTS AND EDGES INTO A COMPLEX
	///		Insert the points.
	//		Insert the edges with adequate references.

	Complex *CAlphaShape = new Complex();

	double Alpha;
	//Alpha_iterator opt = MoleculeAlphaShape.find_optimal_alpha(1);
	//NT alpha_NT = MoleculeAlphaShape.find_alpha_solid();
	//Alpha_iterator first = MoleculeAlphaShape.alpha_lower_bound(alpha_NT);

	//Alpha = *opt;

	Alpha = AlphaSelectionAllBonds(Molecule->GetChemical(), MoleculeAlphaShape);

	//Alpha = Alpha + 0.1*Alpha;


	//cout << "ALPHA OPTIMAL " << Alpha << endl;

	//Alpha = 0.;
	//Alpha *= 10.;

	//InsertAlphaVertices(MoleculeAlphaShape, Alpha_shape_3::INTERIOR, Alpha, CAlphaShape);
	//InsertAlphaVertices(MoleculeAlphaShape, Alpha_shape_3::EXTERIOR, Alpha, CAlphaShape);
	//InsertAlphaVertices(MoleculeAlphaShape, Alpha_shape_3::REGULAR, Alpha, CAlphaShape);
	//InsertAlphaVertices(MoleculeAlphaShape, Alpha_shape_3::SINGULAR, Alpha, CAlphaShape);


	//InsertAlphaEdges(MoleculeAlphaShape, Alpha_shape_3::REGULAR, Alpha, CAlphaShape);
	//InsertAlphaEdges(MoleculeAlphaShape, Alpha_shape_3::INTERIOR, Alpha, CAlphaShape);
	//InsertAlphaEdges(MoleculeAlphaShape, Alpha_shape_3::EXTERIOR, Alpha, CAlphaShape);
	//InsertAlphaEdges(MoleculeAlphaShape, Alpha_shape_3::SINGULAR, Alpha, CAlphaShape);

	InsertFromAlphaCells(MoleculeAlphaShape, Alpha_shape_3::INTERIOR, Alpha, CAlphaShape);
	//InsertFromAlphaCells(MoleculeAlphaShape, Alpha_shape_3::REGULAR, Alpha, CAlphaShape);

	//InsertFromAlphaFacets(MoleculeAlphaShape, Alpha_shape_3::REGULAR, Alpha, CAlphaShape);
	//InsertFromAlphaFacets(MoleculeAlphaShape, Alpha_shape_3::SINGULAR, Alpha, CAlphaShape);
	//InsertFromAlphaFacets(MoleculeAlphaShape, Alpha_shape_3::EXTERIOR, Alpha, CAlphaShape);
	//InsertFromAlphaFacets(MoleculeAlphaShape, Alpha_shape_3::INTERIOR, Alpha, CAlphaShape);

	CAlphaShape->Normalize();


	return CAlphaShape;

}


/*
 *
 *
 *
 */
Complex *AlphaShapeToComplex_ExtraPoints (	MoleculeInfo *Molecule, list <Point *> &ExtraPoints)
{

	//Delaunay_triangulation DT;

	list<Weighted_point> WeightedPList;

	Point *PAux;

	// Insert points into Delaunay triangulation object
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		PAux = Molecule->GetChemical()->GetPointList()[i];


		// Check that the point is of desired type
		if (PAux->GetType() == atom_point)
		{
			CG_point AuxCGP(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2));
			//DT.insert(AuxCGP);

			double Weight = 0.;

			// Insert point into the weighted point list (weighted alpha-shapes)
			WeightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );

		}
	}

	list <Point *>::iterator PIt;

	for (PIt = ExtraPoints.begin(); PIt != ExtraPoints.end(); PIt++)
	{
		PAux = *PIt;

		// Check that the point is of desired type
		if (PAux->GetType() == atom_point)
		{
			CG_point AuxCGP(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2));
			//DT.insert(AuxCGP);

			double Weight = 0.;

			// Insert point into the weighted point list (weighted alpha-shapes)
			WeightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );

		}
	}

	Alpha_shape_3 MoleculeAlphaShape(WeightedPList.begin(), WeightedPList.end(), 0, Alpha_shape_3::GENERAL);

	/////	INSERT ALPHA SHAPE POINTS AND EDGES INTO A COMPLEX
	///		Insert the points.
	//		Insert the edges with adequate references.

	Complex *CAlphaShape = new Complex();

	double Alpha;

	Alpha_iterator opt = MoleculeAlphaShape.find_optimal_alpha(1);

	Alpha = *opt;

	//Alpha = AlphaSelectionAllBonds(Molecule->GetChemical(), MoleculeAlphaShape);

	Alpha = Alpha + 0.1*Alpha;


	InsertFromAlphaCells(MoleculeAlphaShape, Alpha_shape_3::INTERIOR, Alpha, CAlphaShape);

	CAlphaShape->Normalize();


	return CAlphaShape;

}

// TODO add functions to test this (insert weighted points, edges)
/*Complex *WeightedAlphaShapeToComplex (	MoleculeInfo *Molecule)
{

	list<Weighted_point> WeightedPList;

	Point *PAux;

	// Insert points into Delaunay triangulation object
	for (int i = 0; i < Molecule->GetChemical()->GetNPoints(); i++)
	{
		PAux = Molecule->GetChemical()->GetPointList()[i];


		// Check that the point is of desired type
		if (PAux->GetType() == atom_point)
		{
			//double Weight = ComputeWeight(Molecule->GetChemical(), PAux)/2.;
			//double Weight = pow(ComputeWeight(Molecule->GetChemical(), PAux), 2.)/2.;
			double Weight = ComputeWeightMax(Molecule->GetChemical(), PAux);

			//cout << "Weight chosen for point " << Weight << endl;

			// Insert point into the weighted point list (weighted alpha-shapes)
			WeightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), Weight  ) );

		}
	}

	Weighted_alpha_shape_3 MoleculeAlphaShape(WeightedPList.begin(), WeightedPList.end(), 0, Weighted_alpha_shape_3::GENERAL);

	/////	INSERT ALPHA SHAPE POINTS AND EDGES INTO A COMPLEX
	///		Insert the points.
	//		Insert the edges with adequate references.

	Complex *CAlphaShape = new Complex();

	double Alpha = 0.0;
	//W_alpha_iterator opt = MoleculeAlphaShape.find_optimal_alpha(1);
	//Alpha = *opt;

	Alpha -= .75*Alpha;

	//Alpha = Alpha + 0.1*Alpha;

	//Alpha = 0.0*Alpha;
	cout << "ALPHA OPTIMAL " << Alpha << endl;

	InsertAlphaVertices(MoleculeAlphaShape, Weighted_alpha_shape_3::INTERIOR, Alpha, CAlphaShape);
	InsertAlphaVertices(MoleculeAlphaShape, Weighted_alpha_shape_3::EXTERIOR, Alpha, CAlphaShape);
	InsertAlphaVertices(MoleculeAlphaShape, Weighted_alpha_shape_3::REGULAR, Alpha, CAlphaShape);
	InsertAlphaVertices(MoleculeAlphaShape, Weighted_alpha_shape_3::SINGULAR, Alpha, CAlphaShape);


	//InsertAlphaEdges(MoleculeAlphaShape, Weighted_alpha_shape_3::REGULAR, Alpha, CAlphaShape);
	//InsertAlphaEdges(MoleculeAlphaShape, Alpha_shape_3::INTERIOR, Alpha, CAlphaShape);
	//InsertAlphaEdges(MoleculeAlphaShape, Alpha_shape_3::EXTERIOR, Alpha, CAlphaShape);
	//InsertAlphaEdges(MoleculeAlphaShape, Alpha_shape_3::SINGULAR, Alpha, CAlphaShape);

	//InsertFromAlphaCells(MoleculeAlphaShape, Alpha_shape_3::INTERIOR, Alpha, CAlphaShape);
	//InsertFromAlphaCells(MoleculeAlphaShape, Alpha_shape_3::REGULAR, Alpha, CAlphaShape);


	return CAlphaShape;

}*/
