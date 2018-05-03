/*
 * node_characteristic.cc
 *
 *  Created on: Nov 7, 2016
 *      Author: ismael.gomez
 */

#include "node_characteristic.hh"


/*	DivideTriangles
 *
 * 	Function that, given a triangle from the initial icosahedron, divides it into the object equivalent to the self division applied NDivisions times.
 * 	This is constructed iteratively (it could be done recursively) in order to easily build the edges and triangles without duplicates.
 */
void DivideTriangles (	double V1[DIM_SMA],
						double V2[DIM_SMA],
						double V3[DIM_SMA],
						int NDivisions,
						int NTriangle,
						Window *Grid)
{


	///		Number of layers of points
	const int NLayers = pow(2, NDivisions) + 1;
	const int NSteps = NLayers - 1;

	int PointsPerLayer = NLayers;

	/////		INSERT THE POINTS		/////

	///		Number of points
	int NPoints = NLayers*(NLayers+1)/2;

	///		Initial id
	int Id0 = NTriangle*NPoints;

	//cout << "Id0 " << Id0 << endl;

	int Id = Id0;

	//cout << "Layers" << NLayers << endl;

	///		Progression directions
	double	Dir1[DIM_SMA] = { (V2[X_COORD] - V1[X_COORD])/NSteps, (V2[Y_COORD] - V1[Y_COORD])/NSteps, (V2[Z_COORD] - V1[Z_COORD])/NSteps },
			Dir2[DIM_SMA] = { (V3[X_COORD] - V1[X_COORD])/NSteps, (V3[Y_COORD] - V1[Y_COORD])/NSteps, (V3[Z_COORD] - V1[Z_COORD])/NSteps },
			Dir3[DIM_SMA] = { (V3[X_COORD] - V2[X_COORD])/NSteps, (V3[Y_COORD] - V2[Y_COORD])/NSteps, (V3[Z_COORD] - V2[Z_COORD])/NSteps };

	///		Initial point

	double CurrentCoordinates[DIM_SMA];

	//cout << "New triangle" << endl;

	for (int i = 0; i < NLayers; i++)
	{

		CurrentCoordinates[X_COORD] = V1[X_COORD] + i*Dir2[X_COORD];
		CurrentCoordinates[Y_COORD] = V1[Y_COORD] + i*Dir2[Y_COORD];
		CurrentCoordinates[Z_COORD] = V1[Z_COORD] + i*Dir2[Z_COORD];


		for (int j = 0; j < PointsPerLayer; j++)
		{

			if (j != 0 )
			{
				CurrentCoordinates[X_COORD] += Dir1[X_COORD];
				CurrentCoordinates[Y_COORD] += Dir1[Y_COORD];
				CurrentCoordinates[Z_COORD] += Dir1[Z_COORD];
			}

			Point *P = new Point();
			P->SetPointId(Id++);
			P->SetCoordinates(CurrentCoordinates);

			//cout << "Point Id " << P->GetPointId() << endl;

			Grid->InsertPointByValue(P, true);

		}

		PointsPerLayer--;

	}

	/////		INSERT THE EDGES		/////

	///		Horizontal edges

	int IdPrev = Id0,
		IdPos;

	PointsPerLayer = NLayers;

	int CumulativePoints = 0;

	for (int i = 0; i < NLayers; i++)
	{
		IdPrev = Id0 + CumulativePoints;

		for (int j = 1; j < PointsPerLayer; j++)
		{
			IdPos = IdPrev+1;

			Point *PSource = Grid->GetPointById(IdPrev);
			Point *PTarget = Grid->GetPointById(IdPos);

			Grid->InsertEdgeByValue(PSource, PTarget);

			IdPrev++;
		}

		CumulativePoints += PointsPerLayer;
		PointsPerLayer--;

	}

	///		Right oblicuous edges

	PointsPerLayer = NLayers;

	for (int i = 0; i < NLayers; i++)
	{
		IdPrev = Id0 + i;

		double NextLayerJump = NLayers;

		for (int j = 1; j < PointsPerLayer; j++)
		{
			IdPos = IdPrev + NextLayerJump;

			Point *PSource = Grid->GetPointById(IdPrev);
			Point *PTarget = Grid->GetPointById(IdPos);

			Grid->InsertEdgeByValue(PSource, PTarget);

			//cout << "Edge: " << IdPrev << " " << IdPos << endl;

			IdPrev += NextLayerJump;
			NextLayerJump--;
		}

		PointsPerLayer--;
		//cout << "New layer - PPL" << PointsPerLayer << endl;

	}

	//cout << "---" << endl;

	///		Left oblicuous edges

	PointsPerLayer = NLayers;

	for (int i = NLayers-1; i >= 0; i--)
	{
		IdPrev = Id0 + i;

		double NextLayerJump = NLayers-1;

		for (int j = 1; j < PointsPerLayer; j++)
		{
			IdPos = IdPrev + NextLayerJump;

			Point *PSource = Grid->GetPointById(IdPrev);
			Point *PTarget = Grid->GetPointById(IdPos);

			Grid->InsertEdgeByValue(PSource, PTarget);

			//cout << "Edge: " << IdPrev << " " << IdPos << endl;

			IdPrev += NextLayerJump;
			NextLayerJump--;
		}

		PointsPerLayer--;
		//cout << "New layer - PPL" << PointsPerLayer << endl;


	}


	/////		INSERT THE TRIANGLES 		/////

	///		Non-inverted triangles
	int IdTop;

	//IdPrev = Id0;

	//cout << "Id0 end " << Id0 << endl;

	PointsPerLayer = NLayers;

	CumulativePoints = 0;

	for (int i = 0; i < NLayers; i++)
	{
		IdPrev = Id0 + CumulativePoints;

		for (int j = 1; j < PointsPerLayer; j++)
		{
			IdPos = IdPrev+1;
			IdTop = IdPrev+PointsPerLayer;

			Point *P1 = Grid->GetPointById(IdPrev);
			Point *P2 = Grid->GetPointById(IdPos);
			Point *P3 = Grid->GetPointById(IdTop);

			//cout << "Estim ids " << IdPrev << " " << IdPos << " " << IdTop << endl;


			///		Insert the triangle
			Simplex *GridTriangle = new Simplex();
			GridTriangle->InsertPoint(P1);
			GridTriangle->InsertPoint(P2);
			GridTriangle->InsertPoint(P3);

			//cout << "Triangle ids" << P1->GetPointId() << " " << P2->GetPointId() << " " << P3->GetPointId() << endl;

			Grid->InsertSimplex(GridTriangle);

			IdPrev++;
		}

		CumulativePoints += PointsPerLayer;
		PointsPerLayer--;

	}


	///		Inverted triangles

	IdPrev = Id0;

	PointsPerLayer = NLayers-1;
	//nt PointsPrevLayer = NLayers;

	CumulativePoints = NLayers;

	for (int i = 1; i < NLayers-1; i++)
	{
		IdPrev = Id0 + CumulativePoints;

		for (int j = 1; j < PointsPerLayer; j++)
		{
			IdPos = IdPrev+1;
			IdTop = IdPrev-PointsPerLayer;

			Point *P1 = Grid->GetPointById(IdPrev);
			Point *P2 = Grid->GetPointById(IdPos);
			Point *P3 = Grid->GetPointById(IdTop);

			//cout << "T: " << IdPrev << " " << IdPos << " " << IdTop << endl;


			///		Insert the triangle
			Simplex *GridTriangle = new Simplex();
			GridTriangle->InsertPoint(P1);
			GridTriangle->InsertPoint(P2);
			GridTriangle->InsertPoint(P3);

			//cout << "Triangle ids" << P1->GetPointId() << " " << P2->GetPointId() << " " << P3->GetPointId() << endl;

			Grid->InsertSimplex(GridTriangle);

			IdPrev++;
		}

		//cout << "New layer - PPL" << PointsPerLayer << endl;

		CumulativePoints += PointsPerLayer;
		PointsPerLayer--;
		//PointsPrevLayer--;

	}

	PrintToVTK(Grid, "Molecules_Output/Triangle_edges");


}

/*		IcosahedronSphericGrid
 *
 * 		Creates a icosahedron-based spherical grid. The number of recursive iterations (providing n^2 triangles per icosahedron facet) is given in Iterates.
 *
 */
Window *IcosahedronSphericGrid (	double Center[DIM_SMA],
									double Radius,
									int NDivisions)
{


	Window *Grid = new Window();

	//		TODO - Brute-force icosahedron coordinates. Then translate (by Center) and resize (by Radius)

	/////	BASE ICOSAHEDRON

	double Phi = (1. + sqrt(5.))/2.;

	//cout << Phi << endl;

	///		Icosahedron coordinate points

	double 	P0[DIM_SMA]  = {-1.,  Phi, 0},
			P1[DIM_SMA]  = { 1.,  Phi, 0},
			P2[DIM_SMA]  = {-1., -Phi, 0},
			P3[DIM_SMA]  = { 1., -Phi, 0},

			P4[DIM_SMA]  = { 0, -1.,  Phi},
			P5[DIM_SMA]  = { 0,  1.,  Phi},
			P6[DIM_SMA]  = { 0, -1., -Phi},
			P7[DIM_SMA]  = { 0,  1., -Phi},

			P8[DIM_SMA]  = {  Phi, 0, -1},
			P9[DIM_SMA]  = {  Phi, 0,  1},
			P10[DIM_SMA] = { -Phi, 0, -1},
			P11[DIM_SMA] = { -Phi, 0,  1};

	///		Triangles and subdivisions

	//int NDivisions = 2;
	int TriangleId = 0;

	DivideTriangles (P0, P11, P5, NDivisions, TriangleId++, Grid);
	DivideTriangles (P0, P1,  P5, NDivisions, TriangleId++, Grid);
	DivideTriangles (P0, P1,  P7, NDivisions, TriangleId++, Grid);
	DivideTriangles (P0, P10, P7, NDivisions, TriangleId++, Grid);
	DivideTriangles (P0, P10, P11, NDivisions, TriangleId++, Grid);

	//cout << "Hey!" << endl;

	DivideTriangles (P1,  P5,  P9, NDivisions, TriangleId++, Grid);
	DivideTriangles (P11, P5,  P4, NDivisions, TriangleId++, Grid);
	DivideTriangles (P11, P10, P2, NDivisions, TriangleId++, Grid);
	DivideTriangles (P10, P7,  P6, NDivisions, TriangleId++, Grid);
	DivideTriangles (P1,  P7,  P8, NDivisions, TriangleId++, Grid);

	DivideTriangles (P3, P4, P9, NDivisions, TriangleId++, Grid);
	DivideTriangles (P3, P4, P2, NDivisions, TriangleId++, Grid);
	DivideTriangles (P3, P2, P6, NDivisions, TriangleId++, Grid);
	DivideTriangles (P3, P8, P6, NDivisions, TriangleId++, Grid);
	DivideTriangles (P3, P8, P9, NDivisions, TriangleId++, Grid);

	DivideTriangles (P4, P5, P9, NDivisions, TriangleId++, Grid);
	DivideTriangles (P2, P4, P11, NDivisions, TriangleId++, Grid);
	DivideTriangles (P2, P6, P10, NDivisions, TriangleId++, Grid);
	DivideTriangles (P8, P6, P7, NDivisions, TriangleId++, Grid);
	DivideTriangles (P8, P9, P1, NDivisions, TriangleId++, Grid);

	Grid->Normalize();


	/////		PROJECTION TO SPHERE

	for (int i = 0; i < Grid->GetNPoints(); i++)
	{

		Point *GridP = Grid->GetPointList()[i];

		double Cartesian[DIM_SMA], Spherical[DIM_SMA];

		//	Get point adjusted coordinates
		GridP->GetCoordinates(Cartesian);

		//	Get point module
		double VMod = Module(Cartesian);

		//	Adjust the point to lie in the unit sphere
		Cartesian[X_COORD] = Radius*Cartesian[X_COORD]/VMod + Center[X_COORD];
		Cartesian[Y_COORD] = Radius*Cartesian[Y_COORD]/VMod + Center[Y_COORD];
		Cartesian[Z_COORD] = Radius*Cartesian[Z_COORD]/VMod + Center[Z_COORD];

		//	Cartesian to spheric
		/*CartesianToSpheric(Cartesian, Spherical);

		//	Project to sphere surface
		//GridP->SetCoordinate(R_COORD, Radius);
		Spherical[R_COORD] = Radius;

		//	New point to cartesian
		SphericToCartesian(Cartesian, Spherical);*/

		//	Set new coordinates (projected to sphere and displaced)
		GridP->SetCoordinates(Cartesian);

		//GridP->TranslatePoint(Center);

	}


	//PrintWindowToVTK(Grid, "Molecules_Output/Grid_Triangles");
	//PrintToVTK(Grid, "Molecules_Output/Grid_Complex");

	///		Update surface computing after adjusting all triangles
	Grid->SurfaceAreaUpdate();

	return Grid;

}


/*	SphericalGrid
 * 	Generate a spheric grid around a given point with given radius. The grid is constructed in stereographic projection. Number of meridians and cartesian is predefined.
 *
 * 	TODO:	Move to grid file
 */

Window *SphericGrid(double Center[DIM_SMA], double Radius)
{

	int N = 13, M = 14;

	int PCount = 0;

	Window *Grid = new Window();

	//cout << "Init surface: " << Grid->GetSurface() << endl;

	///		Generate the points in spheric coordinates

	//	Poles
	double 	CoordsSouthPole[DIM_SMA] = {0., 0., Radius},
			CoordsNorthPole[DIM_SMA] = {PI, 0., Radius};
			//CoordsNorthPole[DIM_SMA] = {acos(-1), 0., Radius};

	Point *SouthPole = new Point(PCount++, CoordsSouthPole, 0., boundary_point);
	Point *NorthPole = new Point(PCount++, CoordsNorthPole, 0., boundary_point);

	Grid->InsertPointByValue(SouthPole);
	Grid->InsertPointByValue(NorthPole);

	///// TODO - USE INTEGER INDICES
	//	Rest of the points
	for (double Phi = PI/N; Phi < PI; Phi += PI/N)
	//for (int n = 1; n < N; n++)
	{

		//double PhiStep = (double) n/N;
		//double Phi = acos(2*PhiStep - 1);

		for (double Theta = 0; Theta < 2.*PI; Theta += 2.*PI/M)
		{
			//cout << Phi << " " << Theta << endl;

			double Coords[DIM_SMA];

			Coords[PHI_COORD] = Phi;
			Coords[THETA_COORD] = Theta;
			Coords[R_COORD] = Radius;

			Point *GridPoint = new Point(PCount++, Coords, 0., boundary_point);

			Grid->InsertPointByValue(GridPoint);

		}
	}

	///		Connect the points forming edges and windows


	//	Poles' edges and triangles

	for (int i = 0; i < 2; i++)
	{
		double Phi = 0.;
		double PhiStep = 0.;
		Point *Pole;

		//	North pole (only edges from pole to next meridian and triangles formed by groups of three points
		if (i == 0) {
			Phi = PI-(PI/N);
			//PhiStep = (double) 1/N;
			//Phi = acos(2*PhiStep-1);
			Pole = NorthPole;
		}
		//	South pole (same)
		else {
			Phi = PI/N;
			//PhiStep = (double) (N-1)/N;
			//Phi = acos(2*PhiStep-1);
			Pole = SouthPole;
		}

		for (double Theta = 0.; Theta < 2.*PI; Theta += 2.*PI/M)
		{

			double Theta2 = Theta + 2.*PI/M;
			if (Theta2 >= 2.*PI)
			{
				Theta2 = 0.;
			}

			double 	Coords1[DIM_SMA] = {Phi, Theta, Radius},
					Coords2[DIM_SMA] = {Phi, Theta2, Radius};

			Point *P1 = Grid->GetPointByCoordinates(Coords1);
			Point *P2 = Grid->GetPointByCoordinates(Coords2);

			Grid->InsertEdgeByValue(Pole, P1);
			Grid->InsertEdgeByValue(Pole, P2);

			Simplex *GridTriangle = new Simplex();
			GridTriangle->InsertPoint(P1);
			GridTriangle->InsertPoint(P2);
			GridTriangle->InsertPoint(Pole);

			Grid->InsertSimplex(GridTriangle);

		}

	}


	//cout << "---------------------" << endl;

	//	Rest of edges and triangles
	for (double Phi = PI/N; Phi < PI; Phi += PI/N)
	//for (int n = 1; n < N; n++)
	{
		// Next layer
		double Phi2 = Phi + PI/N;
		double Phi0 = Phi - PI/N;

		/*double PhiStep = (double) n/N;
		double Phi = acos(2*PhiStep - 1);

		double Phi0Step = (double) (n+1)/N;
		double Phi0 = acos(2*Phi0Step - 1);

		double Phi2Step = (double) (n-1)/N;
		double Phi2 = acos(2*Phi2Step - 1);*/



		for (double Theta = 0.; Theta < 2.*PI; Theta += 2.*PI/M)
		{

			//cout << "Hey" <<  Phi << " " << Theta << endl;

			double Theta2 = Theta + 2.*PI/M;
			if (Theta2 >= 2.*PI)
			{
				Theta2 = 0.;
			}

			double 	Coords1[DIM_SMA] = {Phi, Theta, Radius},
					Coords2[DIM_SMA] = {Phi, Theta2, Radius},
					Coords3[DIM_SMA] = {Phi2, Theta2, Radius},
					Coords4[DIM_SMA] = {Phi0, Theta, Radius};

			Point *P1 = Grid->GetPointByCoordinates(Coords1);
			Point *P2 = Grid->GetPointByCoordinates(Coords2);

			Grid->InsertEdgeByValue(P1, P2);


			if (Phi2 < PI)
			//if (Phi2 < acos(-1))
			{
				Point *P3 = Grid->GetPointByCoordinates(Coords3);
				Grid->InsertEdgeByValue(P1, P3);

				///		Insert also the triangle
				Simplex *GridTriangle = new Simplex();
				GridTriangle->InsertPoint(P1);
				GridTriangle->InsertPoint(P2);
				GridTriangle->InsertPoint(P3);

				Grid->InsertSimplex(GridTriangle);

			}

			if (Phi0 > 0.)
			{
				Point *P4 = Grid->GetPointByCoordinates(Coords4);
				Grid->InsertEdgeByValue(P1, P4);

				///		Insert the inferior triangle
				Simplex *GridTriangle = new Simplex();
				GridTriangle->InsertPoint(P1);
				GridTriangle->InsertPoint(P2);
				GridTriangle->InsertPoint(P4);

				Grid->InsertSimplex(GridTriangle);
			}

		}
	}

	///		Transform into cartesian coordinates

	for (int i = 0; i < Grid->GetNPoints(); i++)
	{
		double CartesianCoords[DIM_SMA];

		Point *GridPoint = Grid->GetPointList()[i];

		double SphericCoords[DIM_SMA];

		GridPoint->GetCoordinates(SphericCoords);

		SphericToCartesian (CartesianCoords, SphericCoords);

		GridPoint->SetCoordinate(X_COORD, CartesianCoords[X_COORD] + Center[X_COORD]);
		GridPoint->SetCoordinate(Y_COORD, CartesianCoords[Y_COORD] + Center[Y_COORD]);
		GridPoint->SetCoordinate(Z_COORD, CartesianCoords[Z_COORD] + Center[Z_COORD]);



	}

	///		Update surface area after turning coordinates back into cartesian

	Grid->SurfaceAreaUpdate();

	return Grid;

}


/*
 * 	Erases all triangles from list if some edge is missing.
 */
void EraseTrianglesNoEdge (Window *Characteristic)
{

	list <Simplex *>::iterator TrIt;

	//int Count = 0;

	for (TrIt = Characteristic->GetSimplices()->begin(); TrIt != Characteristic->GetSimplices()->end(); TrIt++)
	{

		Simplex *Triangle = *TrIt;

		if (Triangle != NULL)
		{
			Point 	*P0 = Triangle->ReadPoint(0),
					*P1 = Triangle->ReadPoint(1),
					*P2 = Triangle->ReadPoint(2);

			///		If any edge is missing, erase the triangle
			if (Characteristic->GetEdgeByPoints(P0, P1) == NULL)
			{
				//Characteristic->GetSimplices()->erase(TrIt);


				//////
				//vector <Point *> *SimplexPoints = Triangle->GetSimplexPoints();

				//Point 	*P0 = SimplexPoints->at(0),
				//		*P1 = SimplexPoints->at(1),
				//		*P2 = SimplexPoints->at(2);

				//cout << "--- " << endl;
				//P0->PrintPoint();
				//P1->PrintPoint();
				//P2->PrintPoint();

				/////

				TrIt = Characteristic->EraseSimplex(TrIt);
			}
			else if (Characteristic->GetEdgeByPoints(P0, P2) == NULL)
			{
				//Characteristic->GetSimplices()->erase(TrIt);
				TrIt = Characteristic->EraseSimplex(TrIt);
			}
			else if (Characteristic->GetEdgeByPoints(P1, P2) == NULL)
			{
				//Characteristic->GetSimplices()->erase(TrIt);
				TrIt = Characteristic->EraseSimplex(TrIt);
			}
		}




	}



}


/*
 *	Checks whether the bond is a candidate to cross the Boundary-Point triangle.
 */

bool DistanceCrossingEdgeTest(	Edge *Bond,
								Edge *Boundary,
								Point *Node)
{

	///		Get the critical distances

	double 	PointToBound1 = Node->DistanceToPoint(Boundary->GetPSource()),
			PointToBound2 = Node->DistanceToPoint(Boundary->GetPTarget());

	double 	D1A = Bond->GetPSource()->DistanceToPoint(Boundary->GetPSource()),
			D1B = Bond->GetPSource()->DistanceToPoint(Boundary->GetPTarget()),
			D2A = Bond->GetPTarget()->DistanceToPoint(Boundary->GetPSource()),
			D2B = Bond->GetPTarget()->DistanceToPoint(Boundary->GetPTarget());

	double D1 = 0.,
			D2 = 0.;

	if (D1A < D1B)
	{
		D1 = D1A;
		D2 = D2B;
	}
	if (D2A > D2B)
	{
		D2 = D2A;
		D1 = D1B;

	}

	///		Check if the bond is candidate

	double LPDist = 0.;
	if (PointToBound1 < PointToBound2)
	{
		LPDist = PointToBound1;
	}
	else
	{
		LPDist = PointToBound2;
	}

	return (LPDist > D1 && LPDist > D2);
}

/*
 * Given the complete molecule's information, and information of one node (preferably an internal Voronoi node), provides a graph (Window) with the map of all boundary nodes, connected among themselves, in a manner that reflects where the molecular structure
 * divides the space. This gives an intuition of big windows, their size and shape. Unexpected behavior may come out from multi-cage molecules.
 */
Window *NodeCharacteristic(	Point *Node,
							MoleculeInfo *Molecule,
							Window *Boundary,
							double &GridSurface)
{

	Window *Characteristic = new Window();

	//Window *Boundary = SphericGrid(Center, Radius);

	GridSurface = Boundary->GetSurface();

	///		Extract node coordinates
	double NodeCoords[DIM_SMA];
	Node->GetCoordinates(NodeCoords);

	for (int i = 0; i < Boundary->GetNEdges(); i++)
	{

		//	Point 1
		Point *BNode1 = Boundary->GetEdgeList()[i]->GetPSource();

		double 	BNode1Coords[DIM_SMA];
		BNode1->GetCoordinates(BNode1Coords);

		//	Point 2
		Point *BNode2 = Boundary->GetEdgeList()[i]->GetPTarget();

		double BNode2Coords[DIM_SMA];
		BNode2->GetCoordinates(BNode2Coords);

		bool RCTflag = false;

		///		Check all bonds and test whether they cross the triangle or not
		for (int k = 0; k < Molecule->GetChemical()->GetNEdges(); k++)
		{
			Edge *Bond = Molecule->GetChemical()->GetEdgeList()[k];

//time_t t;
//t = clock();

			bool TestCross = DistanceCrossingEdgeTest(Bond, Boundary->GetEdgeList()[i], Node );

//cout << "DECT: " << ((float) (clock() - t) )/CLOCKS_PER_SEC << endl;

			//bool TestCross = true;

			if (TestCross)
			{

				double 	PSourceCoords[DIM_SMA],
						PTargetCoords[DIM_SMA];

				Bond->GetPSource()->GetCoordinates(PSourceCoords);
				Bond->GetPTarget()->GetCoordinates(PTargetCoords);


//t = clock();

				RCTflag = RayCrossesTriangle(PSourceCoords, PTargetCoords, BNode1Coords, BNode2Coords, NodeCoords);

//cout << "RCT: " << ((float) (clock() - t) )/CLOCKS_PER_SEC << endl;

				if (RCTflag)
				{
					break;
				}
			}



		}		///	end for (edges)

		if (!RCTflag)
		{
			///		Add the points to the complex
			Point 	*BNode1New = Characteristic->GetPointByCoordinates(BNode1Coords),
					*BNode2New = Characteristic->GetPointByCoordinates(BNode2Coords);

			if (BNode1New == NULL)
			{
				int Id = Characteristic->GetNPoints();
				BNode1New = new Point(Id, BNode1Coords, 0., boundary_point);
				Characteristic->InsertPointByValue(BNode1New);
			}

			if (BNode2New == NULL)
			{
				int Id = Characteristic->GetNPoints();
				BNode2New = new Point(Id, BNode2Coords, 0., boundary_point);
				Characteristic->InsertPointByValue(BNode2New);
			}

			///		Add the edge
			Edge *BEdge = new Edge(BNode1New, BNode2New, 0.0);

			Characteristic->InsertEdgeByValue(BEdge);

		}

		//}		/// end for (internal loop)

	}		///	end for (external loop)


	///		Copy all triangles from grid to characteristic (if points are present in characteristic)

	list <Simplex *>::iterator TriangleIt;

	for (TriangleIt = Boundary->GetSimplices()->begin(); TriangleIt != Boundary->GetSimplices()->end(); TriangleIt++)
	{
		//	Simplex and points from grid
		Simplex *Triangle = *TriangleIt;

		Point 	*P0 = Triangle->ReadPoint(0),
				*P1 = Triangle->ReadPoint(1),
				*P2 = Triangle->ReadPoint(2);

		double 	CoordsP0[DIM_SMA],
				CoordsP1[DIM_SMA],
				CoordsP2[DIM_SMA];

		P0->GetCoordinates(CoordsP0);
		P1->GetCoordinates(CoordsP1);
		P2->GetCoordinates(CoordsP2);

		//	Use id's to identify points in the new grid
		Point 	*ChP0 = Characteristic->GetPointByCoordinates(CoordsP0),
				*ChP1 = Characteristic->GetPointByCoordinates(CoordsP1),
				*ChP2 = Characteristic->GetPointByCoordinates(CoordsP2);

		if (ChP0 != NULL && ChP1 != NULL && ChP2 != NULL)
		{
			Simplex *CharTriangle = new Simplex();

			CharTriangle->InsertPoint(ChP0);
			CharTriangle->InsertPoint(ChP1);
			CharTriangle->InsertPoint(ChP2);

			Characteristic->InsertSimplex(CharTriangle);
		}


	}


	EraseTrianglesNoEdge(Characteristic);

	Point *NewNode = Node->PointCopy();

	NewNode->SetPointId(Characteristic->GetNPoints());

	Characteristic->InsertPointByValue(NewNode);

	return Characteristic;

}

/*
 *
 */
bool PointInList(Point *TargetPoint, list <Point *> &Points)
{

	list <Point *>::iterator PIt;

	for (PIt = Points.begin(); PIt != Points.end(); PIt++)
	{
		Point *LP = *PIt;

		if (TargetPoint->GetPointId() == LP->GetPointId())
		{
			return true;
		}
	}


	return false;
}

/*
 *
 */

void EraseFromList(Point *TargetPoint, list <Point *> &Points)
{

	list <Point *>::iterator PIt;

	for (PIt = Points.begin(); PIt != Points.end(); PIt++)
	{
		Point *LP = *PIt;

		if (TargetPoint->GetPointId() == LP->GetPointId())
		{
			PIt = Points.erase(PIt);
			delete LP;
		}
	}
}


/*
 *
 */
void InsertTrianglesInComponent ( 	Window *Characteristic,
									Window *Component)
{

	//	If there are less than three points in component, return (no triangle is possibly inside)
	if (Component->GetNPoints() < 3)
	{
		return;
	}

	///		Insert triangles
	list <Simplex *>::iterator CharTrianglesIt;

	for (CharTrianglesIt = Characteristic->GetSimplices()->begin(); CharTrianglesIt != Characteristic->GetSimplices()->end(); CharTrianglesIt++)
	{

		///		Characteristic triangle and points
		Simplex *CharT = *CharTrianglesIt;

		Point 	*P0 = CharT->ReadPoint(0),
				*P1 = CharT->ReadPoint(1),
				*P2 = CharT->ReadPoint(2);

		double 	P0Coords[DIM_SMA],
				P1Coords[DIM_SMA],
				P2Coords[DIM_SMA];

		P0->GetCoordinates(P0Coords);
		P1->GetCoordinates(P1Coords);
		P2->GetCoordinates(P2Coords);

		///		Component points

		Point 	*CompP0 = Component->GetPointByCoordinates(P0Coords),
				*CompP1 = Component->GetPointByCoordinates(P1Coords),
				*CompP2 = Component->GetPointByCoordinates(P2Coords);

		if ( CompP0 != NULL && CompP1 != NULL && CompP2 != NULL)
		{
			Simplex *CompT = new Simplex();

			CompT->InsertPoint(CompP0);
			CompT->InsertPoint(CompP1);
			CompT->InsertPoint(CompP2);

			Component->InsertSimplex(CompT);

		}

	}


}

/*
 *
 */
void InsertEdgesInComponent ( 	Window *Characteristic,
								Window *Component)
{

	//	If there are less than three points in component, return (no triangle is possibly inside)
	if (Component->GetNPoints() < 3)
	{
		return;
	}

	///		Insert edges


	for (int i = 0; i < Characteristic->GetNEdges(); i++)
	{

		///		Characteristic triangle and points
		Edge *CharE = Characteristic->GetEdgeList()[i];

		Point 	*P0 = CharE->GetPSource(),
				*P1 = CharE->GetPTarget();

		double 	P0Coords[DIM_SMA],
				P1Coords[DIM_SMA];

		P0->GetCoordinates(P0Coords);
		P1->GetCoordinates(P1Coords);

		///		Component points

		Point 	*CompP0 = Component->GetPointByCoordinates(P0Coords),
				*CompP1 = Component->GetPointByCoordinates(P1Coords);

		if ( CompP0 != NULL && CompP1 != NULL)
		{
			Component->InsertEdgeByValue(CompP0, CompP1);
		}

	}

}

/*
 * Given a graph (Window), the set of its connected components is computed, each returned in the same form, in a vector
 */
int ConnectedComponents (	Window *Characteristic,
							vector <Window *> &Components)
{

	/////	COMPUTE CONNECTED COMPONENTS
	///		For each point, check connected points and add them to the same component.
	//		Stop when all points from characteristic has been processed.
	//		Store each connected component as a Window * in the provided vector.


	///		Insert all graph's points into a list (to be erased later)

	list <Point *> 	GraphPoints;

	for (int i = 0; i < Characteristic->GetNPoints(); i++)
	{
		Point *LP = Characteristic->GetPointList()[i]->PointCopy();
		GraphPoints.push_back(LP);
	}

	///		Create the different connected components from the edges and points

	list <Point *>::iterator CCPIt, GraphPIt;

	///		Process all points in graph
	while (GraphPoints.size() > 0)
	{
		///		Lists of points for processing
		list <Point *> 	ProcPoints;

		///		Insert first point of general list into unprocessed list

		GraphPIt = GraphPoints.begin();

		Point *GraphPoint = (*GraphPIt)->PointCopy();

		ProcPoints.push_back(GraphPoint);

		Point *PAux = *GraphPIt;
		GraphPIt = GraphPoints.erase(GraphPIt);
		delete PAux;
		//GraphPIt--;

		///		Process connected component
		CCPIt = ProcPoints.begin();

		while (CCPIt != ProcPoints.end())
		{

			Point *CCPoint = *CCPIt;

			//	Get point neighbors
			vector <Point *> Neighbors;
			Characteristic->GetPointNeighbors(CCPoint, Neighbors);

			for (int i = 0; i < Neighbors.size(); i++)
			{
				Point *Neighbor = Neighbors.at(i);

				if (PointInList(Neighbor, GraphPoints))
				{
					//	Add point to processed points list
					ProcPoints.push_back(Neighbor->PointCopy());
					//	Erase the point with equal Id from graph points list
					EraseFromList(Neighbor, GraphPoints);
				}

			}

			//	Step to next element in the list of processed points, end if is the last (none was added)
			CCPIt++;

		}

		///		Store all the points of the new connected component into a new complex (no edge information) and store complex into complex's list
		//Complex *NewConnectedComponent = new Complex();
		Window *NewConnectedComponent = new Window();

		list <Point *>::iterator CCIt;

		for (CCIt = ProcPoints.begin(); CCIt != ProcPoints.end(); CCIt++)
		{
			Point *CCP = *CCIt;
			NewConnectedComponent->InsertPointByValue(CCP);
		}

		Components.push_back(NewConnectedComponent);
	}


	/////	INSERT TRIANGLES IN CONNECTED COMPONENTS
	///		Check every triangle in characteristic
	//		If the three points of the triangle are in the component, add that triangle

	for (int i = 0; i < Components.size(); i++)
	{
		Window *Component = Components.at(i);

		InsertTrianglesInComponent(Characteristic, Component);

		InsertEdgesInComponent(Characteristic, Component);

		//cout << "Surface: " << Component->GetSurface() << endl;

	}


	return Components.size();

}


/*
 * 	Check if the biggest component has relative surface area under certain threshold (with respect to the total surface area of the spheric grid). If beyond this threshold, the criterion fails (false returned).
 * 	Criterion failure potentially results in pore rejection.
 */
double RelativeComponentsSurfaceCriterion (	double GridSurface,
											vector <Window *> &Components)
{

	double MaxSurface = 0.;

	Window *MaxComponent = NULL;

	///		Retrieve component with maximum surface

	for (int i = 0; i < Components.size(); i++)
	{
		Window *Component = Components.at(i);

		if (Component->GetSurface() > MaxSurface)
		{
			MaxSurface = Component->GetSurface();
			MaxComponent = Component;
		}

	}

	///		Check if criterion is fullfilled

	//cout << "Rate: " << MaxSurface/GridSurface << endl;

	return MaxSurface/GridSurface;

}

/*
 * 	Process the shadow of all internal Voronoi nodes.
 * 	Store the rate of the biggest internal Voronoi node as molecule's shadow ratio.
 */

int ProcessShadow (	MoleculeInfo *Molecule)
{

	int Count = 0;

	///		Compute the grid around molecule's center
	/*double Center[DIM_SMA];
	Molecule->GetChemical()->GetComplexCenter(Center);
	double Radius = Molecule->GetChemical()->GetComplexRadius();

	Radius *= 1.5;

	Window *Boundary = SphericGrid(Center, Radius);*/

	///

	double MaxNodeRadius = 0.0;
	double MoleculeRatio = 0.0;

	const int NDivisions = 2;

	//cout << "Radius: " << Radius << endl;



	for (int i = 0; i < Molecule->GetVoronoiGraph()->GetNPoints(); i++)
	{
		Point *VPoint = Molecule->GetVoronoiGraph()->GetPointList()[i];

		//if (VPoint->GetType() == internal_point)
		if (VPoint->GetType() != boundary_point)
		{

			/////
			double Center[DIM_SMA];
			VPoint->GetCoordinates(Center);

			double Radius = Molecule->GetChemical()->GetComplexRadius(VPoint);

			//Window *Boundary = SphericGrid(Center, Radius);
			Window *Boundary = IcosahedronSphericGrid(Center, Radius, NDivisions);

			/////


			double GridSurface;

			Window *Characteristic = NodeCharacteristic(VPoint, Molecule, Boundary, GridSurface);

			vector <Window *> Components;

			ConnectedComponents(Characteristic, Components);



			double Rate = RelativeComponentsSurfaceCriterion(GridSurface, Components);

			//cout << "Rate: " << Rate << endl;

			if (Rate > RELATIVE_SURFACE_LIMIT)
			{
				//cout << "Setting point to external" << endl;
				VPoint->SetType(external_point);
			}
			else
			{
				VPoint->SetType(internal_point);

				Count++;

				//VPoint->PrintPoint();
				//cout << "SR: " << Rate << endl;


				if (VPoint->GetRadius() > MaxNodeRadius)
				{
					MaxNodeRadius = VPoint->GetRadius();
					MoleculeRatio = Rate;
				}

			}

			///		Release characteristic and components
			delete Characteristic;

			for (int i = 0; i < Components.size(); i++)
			{
				Window *Component = Components.at(i);

				delete Component;
			}


			/////
			delete Boundary;
			/////

		}

	}

	Molecule->SetShadowCharact(MoleculeRatio);

	//delete Boundary;

	return Count;

}


/*
 *		Computes shadow of the molecule from its center and applies the test: it succeeds (TRUE) if the relative area (compared with the area of the spherical grid constructed) of the biggest connected component is lesser than certain threshold.
 *		Threshold defined in RELATIVE_SURFACE_LIMIT (node_characteristic.hh).
 */
bool CenterShadowTest(MoleculeInfo *Molecule)
{


	/////		CENTER SHADOW TEST
	///			Compute center point of the box (approximation to the center of the molecule)
	//			Compute spheric grid and point characteristic (connected regions of triangles on the surface of the spheric grid).
	//			Apply test: maximum region is not bigger than certain threshold.


	///		Compute center

	double Center[DIM_SMA];
	Molecule->GetChemical()->GetComplexCenter(Center);
	double Radius = Molecule->GetChemical()->GetComplexRadius();

	//cout << "Center: " << Center[X_COORD] << " " << Center[Y_COORD] << " " << Center[Z_COORD] << endl;

	Radius *= 1.5;

	const int NDivisions = 3;

	//cout << "Radius: " << Radius << endl;

	//Window *Boundary = SphericGrid(Center, Radius);
	Window *Boundary = IcosahedronSphericGrid(Center, Radius, NDivisions);

PrintWindowToVTK(Boundary, "Molecules_Output/SphericGrid");

	///		Point characteristic and connected components

	Point *CenterPoint = new Point();

	CenterPoint->SetCoordinate(X_COORD, Center[X_COORD]);
	CenterPoint->SetCoordinate(Y_COORD, Center[Y_COORD]);
	CenterPoint->SetCoordinate(Z_COORD, Center[Z_COORD]);

	double GridSurface;

	Window *Characteristic = NodeCharacteristic(CenterPoint, Molecule, Boundary, GridSurface);

PrintWindowToVTK(Characteristic, "Molecules_Output/MolCharact");

	vector <Window *> Components;

	ConnectedComponents(Characteristic, Components);

	///		Printing to file
/*	for (int i = 0; i < Components.size(); i++)
	{
		char CCName[FILENAME_MAX];

		strncpy(CCName, "Molecules_Output/Grid_CC_", FILENAME_MAX);

		///		int to string
		char String[FILENAME_MAX];

		int LastDigit = i%10;

		int FirstDigit = (i-LastDigit)/10;

		String[0] = (char) (FirstDigit+48);
		String[1] = (char) (LastDigit+48);
		String[2] = '\0';

		///	copy number

		strncat(CCName, String, FILENAME_MAX);

		PrintWindowToVTK(Components.at(i), CCName);
	}
*/
	///

	///		Apply test

	double Rate = RelativeComponentsSurfaceCriterion(GridSurface, Components);

	//cout << "Center SR: " << Rate << endl;

	//		Store the characteristic
	Molecule->SetShadowCharact(Rate);

	bool Test = Rate < RELATIVE_SURFACE_LIMIT;

	///		Release taken resources

	delete CenterPoint;

	delete Characteristic;

	for (int i = 0; i < Components.size(); i++)
	{
		Window *Component = Components.at(i);

		delete Component;
	}

	delete Boundary;

	return Test;

}

/*
 *
 */
bool MaxIntSphereShadowTest (MoleculeInfo *Molecule)
{
	///		Compute the maximum internal sphere

	Point *VoroMax = NULL;

	for (int i = 0; i < Molecule->GetVoronoiGraph()->GetNPoints(); i++)
	{
		Point *VoroPoint = Molecule->GetVoronoiGraph()->GetPointList()[i];

		if (VoroPoint->GetType() == internal_point)
		{
			if (VoroMax == NULL)
			{
				VoroMax = VoroPoint;
			}
			else if (VoroPoint->GetRadius() > VoroMax->GetRadius())
			{
				VoroMax = VoroPoint;
			}
		}
	}

	///		Compute the shadow for the given point

	double GridSurface;

	double Center[DIM_SMA];
	Molecule->GetChemical()->GetComplexCenter(Center);
	double Radius = Molecule->GetChemical()->GetComplexRadius();

	Radius *= 1.5;

	Window *Boundary = SphericGrid(Center, Radius);

	Window *Characteristic = NodeCharacteristic(VoroMax, Molecule, Boundary, GridSurface);

	vector <Window *> Components;

	ConnectedComponents(Characteristic, Components);

	///		Apply test

	double Rate = RelativeComponentsSurfaceCriterion(GridSurface, Components);

	//		Store the characteristic (if improves the existing one, when it has been computed)
	if (Rate < Molecule->GetShadowCharact() && Molecule->GetShadowCharact() > 0.)
	{
		Molecule->SetShadowCharact(Rate);
	}


	bool Test = Rate < RELATIVE_SURFACE_LIMIT;

	///		Release taken resources

	delete Characteristic;

	for (int i = 0; i < Components.size(); i++)
	{
		Window *Component = Components.at(i);

		delete Component;
	}

	delete Boundary;

	return Test;


}


/*
 *		Test if the points in the list are inside the molecule. Reclassify them as internal/external (if Reclassify == TRUE) and return the number of points that were inside.
 */

int PointListShadowTest (	MoleculeInfo *Molecule,
							list <Point *> &PointList,
							bool Reclassify)
{
	int Count = 0;

	///		Compute the grid around molecule's center
	double Center[DIM_SMA];
	Molecule->GetChemical()->GetComplexCenter(Center);
	double Radius = Molecule->GetChemical()->GetComplexRadius();

	Radius *= 1.5;

	//cout << "Radius: " << Radius << endl;

	Window *Boundary = SphericGrid(Center, Radius);

	list<Point *>::iterator PLIt;

	//for (int i = 0; i < Molecule->GetVoronoiGraph()->GetNPoints(); i++)
	for(PLIt = PointList.begin(); PLIt != PointList.end(); PLIt++)
	{
		Point *VPoint = *PLIt;

		double GridSurface;


//clock_t t;

//cout << "Node char" << endl;

//t = clock();

		Window *Characteristic = NodeCharacteristic(VPoint, Molecule, Boundary, GridSurface);

//cout << ((float) (clock() - t) )/CLOCKS_PER_SEC << endl;

		vector <Window *> Components;

//cout << "CC" << endl;

//t = clock();

		ConnectedComponents(Characteristic, Components);

//cout << ((float) (clock() - t) )/CLOCKS_PER_SEC << endl;
//cout << "Rate " << endl;

//t = clock();

		double Rate = RelativeComponentsSurfaceCriterion(GridSurface, Components);

//cout << ((float) (clock() - t) )/CLOCKS_PER_SEC << endl;


		if (Rate > RELATIVE_SURFACE_LIMIT)
		{
			if (Reclassify)
			{
				VPoint->SetType(external_point);
			}

		}
		else
		{
			if (Reclassify)
			{
				VPoint->SetType(internal_point);
			}

			Count++;

			//VPoint->PrintPoint();

		}

		///		Release characteristic and components
		delete Characteristic;

		for (int i = 0; i < Components.size(); i++)
		{
			Window *Component = Components.at(i);

			delete Component;
		}


	}


	delete Boundary;

	return Count;



}




