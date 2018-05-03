// Voro++, a 3D cell-based Voronoi library
//
// complex.hh
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015

#ifndef _COMPLEX_H_
#define _COMPLEX_H_

//#include "point.hh"
#include "edge.hh"

#include <stdlib.h>
#include <iostream>

#include <vector>

using namespace std;

// Maximum number of points and edges (rough solution to avoid use of realloc)
#define COMPLEX_MAX_POINTS 10000
#define COMPLEX_MAX_EDGES 10000

#define COMPLEX_NAME_SZ 100

class Complex
{

public:

	Complex(Point **PointList, Edge **EdgeList, int NPoints, int NEdges, char *CName);
	Complex(Point **PointList, Edge **EdgeList, int NPoints, int NEdges);
	Complex();
	~Complex();

	// Getters and setters
	/***** *****/
	Point** GetPointList();
	Edge** GetEdgeList();
	int GetNPoints();
	int GetNEdges();
	char *GetCName();

	void SetPointList(Point **PointList);
	void SetEdgeList(Edge **EdgeList);
	void SetNPoints(int NPoints);
	void SetNEdges(int NEdges);
	void SetCName(char *CName);

	// Useful methods for complex processing
	int CalculateNeighborIdList(Point *Point, int *NeighborIdList);

	// Point getters
	Point *GetPointById(int PointId);
	Point *GetPointByCoordinates(double Coordinates[DIM_SMA]);

	void GetPointNeighbors(	Point *TargetPoint, vector <Point *> &Neighbors);

	// Edge getters
	Edge *GetEdgeByIds(int IdSource, int IdTarget);
	Edge *GetEdgeByPoints(Point *Point1, Point *Point2);
	Edge *GetEdgeByCoordinates(double Coordinates1[DIM_SMA], double Coordinates2[DIM_SMA]);

	Edge **GetEdgesByPoint(Point *TargetPoint, int &Size);
	Edge **GetEdgesByPointOppositeRestricted(Point *TargetPoint, PointType OppositeType, int &Size);

	//	Other getters

	void GetComplexBoundary(int Coordinate, double &Min, double &Max);

	void GetComplexCenter(double *Coordinates);

	double GetComplexRadius();
	double GetComplexRadius(Point *P);

	// Complex visualization
	void PrintComplex();
	void PrintComplexPoints(PointType PT);

	void PrintComplexEdges();

	// Complex modificaton
	void InsertPointByValue(int PointId, double Coordinates[DIM_SMA], double Radius, PointType Type);
	void InsertPointByValue(double Coordinates[DIM_SMA], double Radius, PointType Type);
	void InsertPointByValue(double Coordinates[DIM_SMA]);
	void InsertPointByValue(Point *NewPoint);
	bool InsertPointByValue(Point *NewPoint, bool Update);

	void InsertEdgeByValue(Point *SourceP, Point *TargetP, double ChannelRadius);
	void InsertEdgeByValue(Point *SourceP, Point *TargetP);
	void InsertEdgeByValue(Edge *NewEdge);

	void DestroyComplex();

	// Complex copy
	Complex *ComplexCopy();
	Complex *ComplexCopy(PointType PT);

	// Complex prunning
	Complex *PruneComplexOnce();
	bool Pruned();
	Complex *PruneComplex();

	//	Complex normalizing
	void Normalize();

	void Externalize();

	// Special complex checks
	bool IsWindow();


private:

	// Set of points included in the complex
	int NPoints;
	Point **PointList;

	// Set of edges included in the comples
	int NEdges;	
	Edge **EdgeList;

	//	Complex identifier
	char *CName;

};




#endif
