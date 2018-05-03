// Voro++, a 3D cell-based Voronoi library
//
// geometry.hh
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_


#include <cmath>

#include <iostream>

#include "cgal_general.hh"


#include "../eigen/Eigen/Dense"

using namespace Eigen;

using namespace std;


#define FLOAT_ERR	-111111111.

#define X_COORD	0
#define Y_COORD	1
#define Z_COORD 2

#define PHI_COORD 		0
#define THETA_COORD 	1
#define R_COORD			2

#define EQUAL_DIST 0.0000000000001

#define LINEAR_ERR	0.0001

#define PI 3.14159265
#define PI_ERR 0.000001

#define	DIM_SMA	3

// TODO Move this to a Chemistry module
#define BOND_THRESHOLD 		0.4
#define HIDROGEN_BOND_LEN	0.5
#define GENERIC_BOND_LEN	0.8

double Distance 		(	double Point1[DIM_SMA],
							double Point2[DIM_SMA]);

double Module 			(	double Vector[DIM_SMA]);

double DotProduct(	double Vector1[DIM_SMA],
					double Vector2[DIM_SMA]);

void EdgeHalf( 	double End1[DIM_SMA],
				double End2[DIM_SMA],
				double *HalfPoint);

bool RayCrossesPlane 	(	double RayEnd1[DIM_SMA],
							double RayEnd2[DIM_SMA],
							double DVect1[DIM_SMA],
							double DVect2[DIM_SMA],
							double RefPoint[DIM_SMA]);

bool RayPlaneIntersection ( double RayEnd1[DIM_SMA],
							double RayEnd2[DIM_SMA],
							double DVect1[DIM_SMA],
							double DVect2[DIM_SMA],
							double RefPoint[DIM_SMA],
							double *IntersectingPoint );

bool RayCrossFacet (	double RayEnd1[DIM_SMA],
						double RayEnd2[DIM_SMA],
						Facet TFacet);

bool PointInsideCell(	double Point[DIM_SMA],
						Cell_handle Cell);

bool PointInsideTriangle(	double TriangleP1[DIM_SMA],
							double TriangleP2[DIM_SMA],
							double TriangleP3[DIM_SMA],
							double TargetPoint[DIM_SMA]);

bool RayCrossesTriangle (	double RayEnd1[DIM_SMA],
							double RayEnd2[DIM_SMA],
							double TriangleP1[DIM_SMA],
							double TriangleP2[DIM_SMA],
							double TriangleP3[DIM_SMA]);

double Angle (	double Vector1[DIM_SMA],
				double Vector2[DIM_SMA]);

double CellPseudocenter(	Cell_handle Cell,
							double *Pseudocenter);

double TriangleArea (	double Point1[DIM_SMA],
						double Point2[DIM_SMA],
						double Point3[DIM_SMA]);


void CartesianToSpheric (	double *CartCoords,
							double *SphericCoords);

void SphericToCartesian (	double *CartCoords,
							double *SphericCoords);

#endif
