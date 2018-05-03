// Voro++, a 3D cell-based Voronoi library
//
// geometry.cc
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : February 22nd 2015

#include "geometry.hh"

double Distance (	double Point1[DIM_SMA],
					double Point2[DIM_SMA])
{

	double Vector[DIM_SMA] = {Point2[X_COORD] - Point1[X_COORD], Point2[Y_COORD] - Point1[Y_COORD], Point2[Z_COORD] - Point1[Z_COORD]};

	return sqrt(Vector[X_COORD]*Vector[X_COORD] + Vector[Y_COORD]*Vector[Y_COORD] + Vector[Z_COORD]*Vector[Z_COORD]);
}

double Module (	double Vector[DIM_SMA])
{

	return sqrt(Vector[X_COORD]*Vector[X_COORD] + Vector[Y_COORD]*Vector[Y_COORD] + Vector[Z_COORD]*Vector[Z_COORD]);

}

double DotProduct(	double Vector1[DIM_SMA],
					double Vector2[DIM_SMA])
{
	return Vector1[X_COORD]*Vector2[X_COORD] + Vector1[Y_COORD]*Vector2[Y_COORD] + Vector1[Z_COORD]*Vector2[Z_COORD];
}

void EdgeHalf( 	double End1[DIM_SMA],
				double End2[DIM_SMA],
				double *HalfPoint)
{

	HalfPoint[X_COORD] = (End2[X_COORD] - End1[X_COORD])/2.;

	HalfPoint[Y_COORD] = (End2[Y_COORD] - End1[Y_COORD])/2.;

	HalfPoint[Z_COORD] = (End2[Z_COORD] - End1[Z_COORD])/2.;
}

double Angle (	double Vector1[DIM_SMA],
				double Vector2[DIM_SMA])
{

	double Module1 = sqrt(pow(Vector1[X_COORD], 2) + pow(Vector1[Y_COORD], 2) + pow(Vector1[Z_COORD], 2) );
	double Module2 = sqrt(pow(Vector2[X_COORD], 2) + pow(Vector2[Y_COORD], 2) + pow(Vector2[Z_COORD], 2) );

	double DotProduct = Vector1[X_COORD]*Vector2[X_COORD] + Vector1[Y_COORD]*Vector2[Y_COORD] + Vector1[Z_COORD]*Vector2[Z_COORD];

	return acos(DotProduct/(Module1*Module2));

}


/*double Angle2PI (	double Vector1[DIM_SMA],
					double Vector2[DIM_SMA])
{

	//	Get the standard angle among the vectors (-Pi, Pi).
	double AngleStd = Angle(Vector1, Vector2);

	if (AngleStd < 0)
	{
		return PI - AngleStd;
	}

	return AngleStd;

}*/


bool RayCrossesPlane (	double RayEnd1[DIM_SMA],
						double RayEnd2[DIM_SMA],
						double DVect1[DIM_SMA],
						double DVect2[DIM_SMA],
						double RefPoint[DIM_SMA])
{

/*	CG_vector_3 w_vect(DVect1[X_COORD], DVect1[Y_COORD], DVect1[Z_COORD]), v_vect(DVect2[X_COORD], DVect2[Y_COORD], DVect2[Z_COORD]);

	// Vector normal to the facet
	CG_vector_3 n_vect = CGAL::cross_product(w_vect,v_vect);

	// Vectors from reference point (Point 3) to both ends of the ray
	CG_vector_3 	u_vect (RayEnd1[X_COORD] - RefPoint[X_COORD], RayEnd1[Y_COORD] - RefPoint[Y_COORD], RayEnd1[Z_COORD] - RefPoint[Z_COORD]),
					up_vect(RayEnd2[X_COORD] - RefPoint[X_COORD], RayEnd2[Y_COORD] - RefPoint[Y_COORD], RayEnd2[Z_COORD] - RefPoint[Z_COORD]);

	double 	Scalar1 = u_vect*n_vect,
			Scalar2 = up_vect*n_vect;

	// Collinear case: either one of the scalar products is 0 (the ray end is placed exactly in the facet)
	if (Scalar1 == 0 || Scalar2 == 0)
	{
		return false;
	}

	// If both ends of the ray are at different sides of the facet, their scalar products with n will have different signs
	if ( (Scalar1 > 0 && Scalar2 < 0) || (Scalar1 < 0 && Scalar2 > 0) )
	{
		return true;
	}
*/
	return false;

}

/*
 * Compute the point of intersection between the given ray and plane
 */
bool RayPlaneIntersection ( double RayEnd1[DIM_SMA],
							double RayEnd2[DIM_SMA],
							double DVect1[DIM_SMA],
							double DVect2[DIM_SMA],
							double RefPoint[DIM_SMA],
							double *IntersectingPoint )
{


	/////	COMPUTE INTERSECTION BETWEEN PLANE AND RAY
	///		Point of intersection satisfies the following conditions:
	//		(a)	P = Lambda1*U + Lambda2*V + RefP 	[U,V director vectors of plane]
	//		(b) P = RayEnd1 + Alpha*E				[E director vector of ray]
	//		Then, it follows:
	//		(c) Lambda1*U + Lambda2*V - Alpha*E = RayEnd1 - RefP
	//		This gives the values for the intersection between the infinite line and the plane.
	//		But the ray may not intersect the plane.
	//		(d)	RayEnd2 = Beta*E + RayEnd1
	//		If abs(Beta) < abs(Alpha) the ray does not cross the plane.
	//
	//		Algorithm:
	//		(1) Compute RayEnd1 - RefP.
	//		(2)	Solve the system to obtain Lambda1, Lambda2, Alpha.
	//		(3) Use the parameters to compute the intersection point.
	//		(4)	Compute Beta and compare it with Alpha.


	// (1)
/*	Vector3d b( RayEnd1[X_COORD] - RefPoint[X_COORD], RayEnd1[Y_COORD] - RefPoint[Y_COORD], RayEnd1[Z_COORD] - RefPoint[Z_COORD]);

	Matrix3d A;
	Vector3d x;

	// Fill the matrix with the coordinates of the three segments that will allow computing the point
	for (int i = 0; i < DIM_SMA; i++)
	{
		A(i,0) = DVect1[i];
	}
	for (int i = 0; i < DIM_SMA; i++)
	{
		A(i,1) = DVect2[i];
	}
	for (int i = 0; i < DIM_SMA; i++)
	{
		A(i,2) = RayEnd2[i] - RayEnd1[i];
	}

	// (2)

	x = A.colPivHouseholderQr().solve(b);

	bool a_solution_exists = (A*x).isApprox(b, LINEAR_ERR);

	if (!a_solution_exists)
	{
		return false;
	}

	double 	Lambda1 = x(0),
			Lambda2 = x(1),
			Alpha = x(2);

	//cout << "Lambda 1: " << Lambda1 << " - Lambda 2: " << Lambda2 << " - Alpha: " << Alpha << endl;

	// (3)

	for (int i = 0; i < DIM_SMA; i++)
	{
		IntersectingPoint[i] = RefPoint[i] + Lambda1*DVect1[i] + Lambda2*DVect2[i];
	}

	// (4)
	//if (Alpha > 1 || Alpha < 0)
	if (abs(Alpha) > 1 || Alpha == 0)
	{
		return false;
	}


	//////////


	return true;*/

	/////		ALTERNATIVE METHOD
	///			TODO: Test and erase previous (if this is better)

	/////		COMPUTE INTERSECTION BETWEEN PLANE AND RAY
	///			(1) Compute the normal to the plane (N = UxV)
	//			(2) Points in the ray are given by the equation P = RayEnd1 + r(RayEnd2 - RayEnd1)
	//			(3) Intersection point happens for r = N*(RefPoint - RayEnd1)/(N*(RayEnd2 - RayEnd1))
	//			Method based on the ideas and developments from: http://geomalgorithms.com/a06-_intersect-2.html

	///		(1)
	double Normal[DIM_SMA] = {0., 0., 0.};

	Normal[X_COORD] = DVect1[Y_COORD]*DVect2[Z_COORD] - DVect1[Z_COORD]*DVect2[Y_COORD];

	Normal[Y_COORD] = -DVect1[X_COORD]*DVect2[Z_COORD] + DVect1[Z_COORD]*DVect2[X_COORD];

	Normal[Z_COORD] = DVect1[X_COORD]*DVect2[Y_COORD] - DVect1[Y_COORD]*DVect2[X_COORD];



	///		(3)

	double 	V1[DIM_SMA] = { RefPoint[X_COORD] - RayEnd1[X_COORD],
							RefPoint[Y_COORD] - RayEnd1[Y_COORD],
							RefPoint[Z_COORD] - RayEnd1[Z_COORD]};

	double	VP[DIM_SMA] = {	RayEnd2[X_COORD] - RayEnd1[X_COORD],
							RayEnd2[Y_COORD] - RayEnd1[Y_COORD],
							RayEnd2[Z_COORD] - RayEnd1[Z_COORD]};

	double 	N_V1 = DotProduct(Normal, V1),
			N_VP = DotProduct(Normal, VP);

	if (N_VP == 0)
	{
		return false;
	}

	double r = N_V1/N_VP;

	if (r < 0 || r > 1)
	{
		return false;
	}

	///		(2)

	IntersectingPoint[X_COORD] = RayEnd1[X_COORD] + r*VP[X_COORD];
	IntersectingPoint[Y_COORD] = RayEnd1[Y_COORD] + r*VP[Y_COORD];
	IntersectingPoint[Z_COORD] = RayEnd1[Z_COORD] + r*VP[Z_COORD];

	return true;



}


/*	RetrievePointFromFacet
 * 		Auxiliary function for RayCrossFacet. Retrieves coordinates of a point from a CGAL-type facet.
 *
 */
void RetrievePointFromFacet (	Facet TargetFacet,
								int FPIndex,
								double *Point)
{

	int Opposite = TargetFacet.second;

	int Count = 0;

	for (int i = 0; i < CELL_N_V; i++)
	{
		if (i != Opposite) {
			if(Count == FPIndex)
			{
				Point[X_COORD] = TargetFacet.first->vertex(i)->point().x();
				Point[Y_COORD] = TargetFacet.first->vertex(i)->point().y();
				Point[Z_COORD] = TargetFacet.first->vertex(i)->point().z();
				return;
			}
			else {
				Count++;
			}
		}
	}

}

bool RayCrossFacetPlane (	double RayEnd1[DIM_SMA],
							double RayEnd2[DIM_SMA],
							Facet TFacet)
{

	// Get the points from the facets
	double Point1[DIM_SMA], Point2[DIM_SMA], Point3[DIM_SMA];

	RetrievePointFromFacet(TFacet, 0, Point1);
	RetrievePointFromFacet(TFacet, 1, Point2);
	RetrievePointFromFacet(TFacet, 2, Point3);

	// Director vectors for the plane formed by the facet
	double 	DVect1[DIM_SMA] = {Point3[X_COORD] - Point1[X_COORD], Point3[Y_COORD] - Point1[Y_COORD], Point3[Z_COORD] - Point1[Z_COORD]},
			DVect2[DIM_SMA] = {Point3[X_COORD] - Point2[X_COORD], Point3[Y_COORD] - Point2[Y_COORD], Point3[Z_COORD] - Point2[Z_COORD]};

	CG_vector_3 w_vect(DVect1[X_COORD], DVect1[Y_COORD], DVect1[Z_COORD]), v_vect(DVect2[X_COORD], DVect2[Y_COORD], DVect2[Z_COORD]);

	// Vector normal to the facet
	CG_vector_3 n_vect = CGAL::cross_product(w_vect,v_vect);

	// Vectors from reference point (Point 3) to both ends of the ray
	CG_vector_3 	u_vect(RayEnd1[X_COORD] - Point3[X_COORD], RayEnd1[Y_COORD] - Point3[Y_COORD], RayEnd1[Z_COORD] - Point3[Z_COORD]),
					up_vect(RayEnd2[X_COORD] - Point3[X_COORD], RayEnd2[Y_COORD] - Point3[Y_COORD], RayEnd2[Z_COORD] - Point3[Z_COORD]);

	double Scalar1 = u_vect*n_vect, Scalar2 = up_vect*n_vect;

	// Collinear case: either one of the scalar products is 0 (the ray end is placed exactly in the facet)
	if (Scalar1 == 0 || Scalar2 == 0)
		return true;

	// If both ends of the ray are at different sides of the facet, their scalar products with n will have different signs
	if ( (Scalar1 > 0 && Scalar2 < 0) || (Scalar1 < 0 && Scalar2 > 0) )
		return true;

	return false;

	//return RayCrossesTriangle(RayEnd1, RayEnd2, Point1, Point2, Point3);





}


/*
 * Checks if a point is inside a tetrahedral cell
 */
bool PointInsideCell(	double Point[DIM_SMA],
						Cell_handle Cell)
{

	for (int i = 0; i < CELL_N_V; i++)
	{
		// Get the facet
		Facet TFacet(Cell, i);

		// Get the point with which the ray is formed (one node of the cell)
		double CellPoint[DIM_SMA] = { Cell->vertex(i)->point().x(), Cell->vertex(i)->point().y(), Cell->vertex(i)->point().z() };

		// Check if the ray formed by the cell point and the node of interest crosses the facet (if it does it, the point is outside the cell
		if (RayCrossFacetPlane(Point, CellPoint, TFacet))
		{
			return false;
		}
	}

	// If none of the points crossed the facet, the voronoi node is inside the given cell
	return true;

}

/*
 * Computes cell centroid or center of mass of a cell (tetrahedron). Returns the distance from the centroid to the most distant vertex of the cell.
 */
double CellPseudocenter(	Cell_handle Cell,
							double *Pseudocenter)
{
	Pseudocenter[X_COORD] = 0.;
	Pseudocenter[Y_COORD] = 0.;
	Pseudocenter[Z_COORD] = 0.;

	for (int i = 0; i < CELL_N_V; i++)
	{
		// Get the point with which the ray is formed (one node of the cell)
		double CellPoint[DIM_SMA] = { Cell->vertex(i)->point().x(), Cell->vertex(i)->point().y(), Cell->vertex(i)->point().z() };

		Pseudocenter[X_COORD] += CellPoint[X_COORD];
		Pseudocenter[Y_COORD] += CellPoint[Y_COORD];
		Pseudocenter[Z_COORD] += CellPoint[Z_COORD];

	}

	Pseudocenter[X_COORD] /= 4.;
	Pseudocenter[Y_COORD] /= 4.;
	Pseudocenter[Z_COORD] /= 4.;

	double MaxDistance = 0.;

	for (int i = 0; i < CELL_N_V; i++)
	{
		// Get the point with which the ray is formed (one node of the cell)
		double CellPoint[DIM_SMA] = { Cell->vertex(i)->point().x(), Cell->vertex(i)->point().y(), Cell->vertex(i)->point().z() };

		double DistanceToCP = Distance(CellPoint, Pseudocenter);

		if (MaxDistance < DistanceToCP)
		{
			MaxDistance = DistanceToCP;
		}

	}

	return MaxDistance;

}

/*
 * Checks if a point is inside a triangle PROVIDED THAT IS IN THE SAME PLANE.
 */
bool PointInsideTriangle(	double TriangleP1[DIM_SMA],
							double TriangleP2[DIM_SMA],
							double TriangleP3[DIM_SMA],
							double TargetPoint[DIM_SMA])
{

	// Compute distances between the point and the triangle vertices
	/*double 	D1 = Distance(TriangleP1, TargetPoint),
			D2 = Distance(TriangleP2, TargetPoint),
			D3 = Distance(TriangleP3, TargetPoint);

	// Compute triangle segment lengths
	double 	L12 = Distance(TriangleP1, TriangleP2),
			L13 = Distance(TriangleP1, TriangleP3),
			L23 = Distance(TriangleP2, TriangleP3);

	/////	DISTANCE CRITERION
	///		Point is inside the triangle iff it's distance to any vertex is lesser than the distance between that vertex and the other two.
	//		NOTE: The point is assumed to be in the same plane than the triangle.
	//		(1) Compare each distance with the corresponding lengths, returning false if any of them is lesser than the edge length.
	//		(2) If criterion does not fail, the point is inside the triangle (return true).

	// (1)
	if (D1 > L12 || D1 > L13)
	{
		return false;
	}

	if (D2 > L12 || D2 > L23)
	{
		return false;
	}

	if (D3 > L13 || D3 > L23)
	{
		return false;
	}*/

	/////	ANGLE CRITERION
	///		Check the angles [0, 2*pi] between the sides and the vectors that connect triangle points with the point of interest.
	//		The angle formed with the point must be less than the angles between the sides of the triangle.


	///		Get the side vectors
	/*double 	SideA[DIM_SMA] = {TriangleP2[X_COORD] - TriangleP1[X_COORD], TriangleP2[Y_COORD] - TriangleP1[Y_COORD], TriangleP2[Z_COORD] - TriangleP1[Z_COORD]},
			SideB[DIM_SMA] = {TriangleP3[X_COORD] - TriangleP2[X_COORD], TriangleP3[Y_COORD] - TriangleP2[Y_COORD], TriangleP3[Z_COORD] - TriangleP2[Z_COORD]},
			SideC[DIM_SMA] = {TriangleP3[X_COORD] - TriangleP1[X_COORD], TriangleP3[Y_COORD] - TriangleP1[Y_COORD], TriangleP3[Z_COORD] - TriangleP1[Z_COORD]};

	///		Get the side angles
	double 	AngleAC = Angle(SideA, SideC),
			AngleBC = Angle(SideB, SideC),
			AngleAB = Angle(SideA, SideB);*/

	///		Get the Point-to-Vertex vectors
	double 	TP1toPoint[DIM_SMA] = {TriangleP1[X_COORD] - TargetPoint[X_COORD], TriangleP1[Y_COORD] - TargetPoint[Y_COORD], TriangleP1[Z_COORD] - TargetPoint[Z_COORD]},
			TP2toPoint[DIM_SMA] = {TriangleP2[X_COORD] - TargetPoint[X_COORD], TriangleP2[Y_COORD] - TargetPoint[Y_COORD], TriangleP2[Z_COORD] - TargetPoint[Z_COORD]},
			TP3toPoint[DIM_SMA] = {TriangleP3[X_COORD] - TargetPoint[X_COORD], TriangleP3[Y_COORD] - TargetPoint[Y_COORD], TriangleP3[Z_COORD] - TargetPoint[Z_COORD]};

	///		Get the angles among sides and ptv vectors
	double	Angle1 = Angle(TP1toPoint, TP2toPoint),
			Angle2 = Angle(TP2toPoint, TP3toPoint),
			Angle3 = Angle(TP3toPoint, TP1toPoint);

	///		If the three angles sum 2PI, the point is inside the triangle. Otherwise it's not.

	//cout << "Angles: " << Angle1 << " " << Angle2 << " " << Angle3 << " - Sum: " << Angle1 + Angle2 + Angle3 << endl;

	if ( abs(Angle1 + Angle2 + Angle3 - 2.*PI) < PI_ERR)
	//if ( Angle1 + Angle2 + Angle3 == 2.*PI)
	{
		return true;
	}

	return false;


}


bool RayCrossesTriangle (	double RayEnd1[DIM_SMA],
							double RayEnd2[DIM_SMA],
							double TriangleP1[DIM_SMA],
							double TriangleP2[DIM_SMA],
							double TriangleP3[DIM_SMA])
{

	/////	RAY CROSSES TRIANGLE
	///		Check if the ray (defined by its points) crosses the surface trapped by the triangle (given by the coordinates of it's vertices)
	//		(1) Decide whether the ray crosses the plane in which the triangle is inscribed.
	//		(2) In case it crosses the plane, compute the point P that intersects with that plane.
	//		(3) Once computed the point, check whether it's contained into the triangle (provided 1 and 2).
	//		(4) Return true if all the previous criteria didn't fail.

	// (1)

	// Compute director vectors for the plane

	double 	DVect1[DIM_SMA] = { TriangleP2[X_COORD] - TriangleP1[X_COORD], TriangleP2[Y_COORD] - TriangleP1[Y_COORD], TriangleP2[Z_COORD] - TriangleP1[Z_COORD]},
			DVect2[DIM_SMA] = { TriangleP3[X_COORD] - TriangleP1[X_COORD], TriangleP3[Y_COORD] - TriangleP1[Y_COORD], TriangleP3[Z_COORD] - TriangleP1[Z_COORD]};


//	if( !RayCrossesPlane ( RayEnd1, RayEnd2, DVect1, DVect2, TriangleP1) )
//	{
//		return false;
//	}


	// (2)

	double Intersect[DIM_SMA];

	if ( !RayPlaneIntersection ( RayEnd1, RayEnd2, DVect1, DVect2, TriangleP1, Intersect ) )
	{
		//cout << "Warning: No intersecting point was found between plane and ray (RayCrossTriangle)!!!" << endl;
		return false;
	}

	//cout << "Point crossed: " << Intersect[0] << " " << Intersect[1] << " " << Intersect[2] << endl;

	// (3)

	if ( !PointInsideTriangle(TriangleP1, TriangleP2, TriangleP3, Intersect) )
	{
		//cout << Intersect[X_COORD] << " " << Intersect[Y_COORD] << " " << Intersect[Z_COORD] << " " << endl;

		//cout << "No PIT" << endl;

		return false;
	}

	// (4)

	return true;


	// Director vectors for the plane formed by the facet
/*	double 	DVect1[DIM_SMA] = {TriangleP3[X_COORD] - TriangleP1[X_COORD], TriangleP3[Y_COORD] - TriangleP1[Y_COORD], TriangleP3[Z_COORD] - TriangleP1[Z_COORD]},
			DVect2[DIM_SMA] = {TriangleP3[X_COORD] - TriangleP2[X_COORD], TriangleP3[Y_COORD] - TriangleP2[Y_COORD], TriangleP3[Z_COORD] - TriangleP2[Z_COORD]};

	CG_vector_3 	w_vect(DVect1[X_COORD], DVect1[Y_COORD], DVect1[Z_COORD]),
					v_vect(DVect2[X_COORD], DVect2[Y_COORD], DVect2[Z_COORD]);

	// Vector normal to the facet
	CG_vector_3 n_vect = CGAL::cross_product(w_vect,v_vect);

	// Vectors from reference point (Point 3) to both ends of the ray
	CG_vector_3 	u_vect(RayEnd1[X_COORD] - TriangleP3[X_COORD], RayEnd1[Y_COORD] - TriangleP3[Y_COORD], RayEnd1[Z_COORD] - TriangleP3[Z_COORD]),
					up_vect(RayEnd2[X_COORD] - TriangleP3[X_COORD], RayEnd2[Y_COORD] - TriangleP3[Y_COORD], RayEnd2[Z_COORD] - TriangleP3[Z_COORD]);

	double 	Scalar1 = u_vect*n_vect,
			Scalar2 = up_vect*n_vect;

	// Collinear case: either one of the scalar products is 0 (the ray end is placed exactly in the facet)
	if (Scalar1 == 0 || Scalar2 == 0)
		return true;

	// If both ends of the ray are at different sides of the facet, their scalar products with n will have different signs
	if ( (Scalar1 > 0 && Scalar2 < 0) || (Scalar1 < 0 && Scalar2 > 0) )
		return true;

	return false;
*/

}

/*
 *
 */
double TriangleArea (	double Point1[DIM_SMA],
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


}


/*
 *
 */

void CartesianToSpheric (double *CartCoords, double *SphericCoords)
{

	// Calculate the phi coordinate (interval [0, 2*PI] )
	SphericCoords[PHI_COORD] = atan2(CartCoords[Y_COORD], CartCoords[X_COORD]);
	if (SphericCoords[PHI_COORD] < 0)
	{
		SphericCoords[PHI_COORD] = 2*PI + SphericCoords[PHI_COORD];
	}

	// Calcuate the theta coordinate (interval [0, PI] )
	double Aux = sqrt( CartCoords[X_COORD]*CartCoords[X_COORD] + CartCoords[Y_COORD]*CartCoords[Y_COORD] );

	SphericCoords[THETA_COORD] = atan2(Aux, CartCoords[Z_COORD]);

	if (SphericCoords[THETA_COORD] < 0)
	{
		SphericCoords[THETA_COORD] = PI + SphericCoords[THETA_COORD];
	}

	// Calculate the radius
	SphericCoords[R_COORD] = sqrt( pow(CartCoords[X_COORD], 2.) + pow(CartCoords[Y_COORD], 2.) + pow(CartCoords[Z_COORD], 2.) );

}


void SphericToCartesian (double *CartCoords, double *SphericCoords)
{

	double R, Phi, Theta;

	R = SphericCoords[R_COORD];
	Phi = SphericCoords[PHI_COORD];
	Theta = SphericCoords[THETA_COORD];

	CartCoords[X_COORD] = R*sin(Phi)*cos(Theta);
	CartCoords[Y_COORD] = R*sin(Phi)*sin(Theta);
	CartCoords[Z_COORD] = R*cos(Phi);

}





