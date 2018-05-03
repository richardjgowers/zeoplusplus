// Voro++, a 3D cell-based Voronoi library
//
// edge.hh
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015

#ifndef _EDGE_H_
#define _EDGE_H_

//#include "geometry.hh"
#include "point.hh"

#include <iostream>
using namespace std;

class Edge 
{

public:

	// Constructor and destructor
	Edge(Point *PSource, Point *PTarget, double ChannelRadius);
	~Edge();

	// Getters and setters
	Point* GetPSource();
	Point* GetPTarget();
	double GetChannelRadius();

	void SetPSource(Point *PSource);
	void SetPTarget(Point *PTarget);
	void SetChannelRadius(double ChannelRadius);

	// Special getters and setters
	Point *GetOpposite(Point *TargetPoint);

	void GetMiddlePoint(double *Coordinates);

	// Other methods
	void PrintEdge();
	bool IsEntryEdge();

private:

	Point *PSource;
	Point *PTarget;

	double ChannelRadius;
	
};

#endif
