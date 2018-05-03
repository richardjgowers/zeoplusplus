// Voro++, a 3D cell-based Voronoi library
//
// point.hh
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015

#ifndef _POINT_H_
#define _POINT_H_

#include "geometry.hh"

#include <cmath>

#include <iostream>
using namespace std;

typedef enum ptype_
{
	internal_point,		// Point is inside the molecule
	external_point,		// Point is outside the molecule but not in the boundary of the cage
	boundary_point,		// Point is part of the boundary of the cage
	hybrid_point,		// Point is unclassifiable
	//molecule_point,		// Point is part of the molecule
	atom_point,
	undefined_point		// Point type is unknown (typically internal or external)

} PointType;


class Point 
{

public:

	// Constructor and destructor
	Point(int PointId, int Degree, double Coordinates[DIM_SMA], double Radius, PointType Type);
	Point(int PointId, double Coordinates[DIM_SMA], double Radius, PointType Type);
	Point();
	~Point();

	// Getters and setters
	int GetPointId();
	double *GetCoordinates();
	void GetCoordinates(double *Coordinates);
	double GetCoordinate(int Position);
	double GetRadius();
	PointType GetType();
	int GetDegree();
	void PointCoordinates(double *Coordinates);

	void SetPointId(int PointId);
	void SetCoordinate(int Position, double Value);
	void SetCoordinates(double *Coordinates);
	void SetRadius(double Radius);
	void SetType(PointType Type);
	void SetDegree(int Degree);

	// Point manipulation
	void IncreaseDegree(int Value);

	// I/O
	void PrintPoint();

	// Comparison
	bool SameCoordinates(double Coordinates[DIM_SMA]);

	// Object copy
	Point *PointCopy();

	// Distance
	double DistanceToPoint(Point *TargetP);

	void TranslatePoint(double *Coordinates);



private:

	/*****	GRAPHICAL INFORMATION */
	/***	Information about the graph structure associated with the complex:
	 *	- Vertex identifier
	 *	- Degree of the vertex (number of edges touching it)
	 */
	
	// Point identifier - Must be unique (in complex)
	// Not mutable (constant)
	int PointId;

	int Degree;

	/***** Geometrical information *****/

	// Spatial coordinates of the point (its center)
	double Coordinates[DIM_SMA];

	// For the case of different radii (otherwise Radius = 1)
	double Radius;

	// Type of point (for processing)
	PointType Type;

	
};

#endif
