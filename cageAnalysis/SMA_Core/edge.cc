// Voro++, a 3D cell-based Voronoi library
//
// edge.cc
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015


#include "edge.hh"

/***** CONSTRUCTOR AND DESTRUCTOR *****/

Edge::Edge (Point *PSource, Point *PTarget, double ChannelRadius)
{
	this->PSource = PSource;
	this->PTarget = PTarget;
	this->ChannelRadius = ChannelRadius;
}


Edge::~Edge()
{
	//cout << "Edge destructor called" << endl;
}

/***** GETTERS AND SETTERS *****/

Point* Edge::GetPSource()
{
	return this->PSource;
}


Point* Edge::GetPTarget()
{
	return this->PTarget;
}

double Edge::GetChannelRadius()
{
	return this->ChannelRadius;
}

void Edge::SetPSource(Point *PSource)
{
	this->PSource=PSource;
}

void Edge::SetPTarget(Point *PTarget)
{
	this->PTarget=PTarget;
}


void Edge::SetChannelRadius(double ChannelRadius)
{
	this->ChannelRadius=ChannelRadius;
}

/***** OTHER METHODS *****/


Point* Edge::GetOpposite(Point *TargetPoint)
{

	if (this->PSource->GetPointId() == TargetPoint->GetPointId())
		return this->PTarget;

	if (this->PTarget->GetPointId() == TargetPoint->GetPointId())
		return this->PSource;

	return NULL;


}

void Edge::GetMiddlePoint(double *Coordinates)
{

	double 	SCoords[DIM_SMA],
			TCoords[DIM_SMA];

	this->PSource->GetCoordinates(SCoords);
	this->PTarget->GetCoordinates(TCoords);

	Coordinates[X_COORD] = (SCoords[X_COORD] + TCoords[X_COORD])/2.;
	Coordinates[Y_COORD] = (SCoords[Y_COORD] + TCoords[Y_COORD])/2.;
	Coordinates[Z_COORD] = (SCoords[Z_COORD] + TCoords[Z_COORD])/2.;

}


void Edge::PrintEdge()
{

	cout << "Edge: ( " << this->PSource->GetPointId() << " , " << this->PTarget->GetPointId() << " )" << " - Channel Radius " << this->ChannelRadius << "\n"; 


}

/*Edge* Edge::EdgeCopy()
{
	Edge *NewEdge = new Edge(	this->PSource->PointCopy(), 
					this->PTarget->PointCopy, 
					this->ChannelRadius);
}*/

bool Edge::IsEntryEdge()
{

	PointType 	T1 = this->PSource->GetType(),
				T2 = this->PTarget->GetType();

	if (T1 == external_point || T1 == boundary_point)
	{
		if (T2 == internal_point)
		{
			return true;
		}
	}

	if (T2 == external_point || T2 == boundary_point)
	{
		if (T1 == internal_point)
		{
			return true;
		}
	}

	return false;
}




