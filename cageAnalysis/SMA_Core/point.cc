// Voro++, a 3D cell-based Voronoi library
//
// point.cc
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015


#include "point.hh"

/***** CONSTRUCTOR AND DESTRUCTOR ****/
Point::Point(	int PointId, 
		int Degree,
		double Coordinates[DIM_SMA], 
		double Radius, 
		PointType Type)
{

	this->PointId = PointId;

	this->Degree = Degree;	

	for (int i=0; i<DIM_SMA; i++)
		this->Coordinates[i] = Coordinates[i];
	
	this->Radius = Radius;

	this->Type = Type;


}

Point::Point(	int PointId, 
				double Coordinates[DIM_SMA],
				double Radius,
				PointType Type)
{

	this->PointId = PointId;

	// When constructed without degree, this is set to 0
	this->Degree = 0;	

	for (int i=0; i<DIM_SMA; i++)
		this->Coordinates[i] = Coordinates[i];
	
	this->Radius = Radius;

	this->Type = Type;


}

Point::Point()
{

	this->PointId = -1;
	this->Degree = 0;

	for (int i=0; i<DIM_SMA; i++)
		this->Coordinates[i] = 0.;

	this->Radius = 0.;

	this->Type = undefined_point;

}


Point::~Point()
{
	//cout << "Point destructor called" << endl;
}

/***** GETTERS *****/

int Point::GetPointId()
{
	return this->PointId;

}


double* Point::GetCoordinates()
{
	return this->Coordinates;

}

void Point::GetCoordinates(double *Coordinates)
{
	Coordinates[X_COORD] = this->Coordinates[X_COORD];
	Coordinates[Y_COORD] = this->Coordinates[Y_COORD];
	Coordinates[Z_COORD] = this->Coordinates[Z_COORD];
}


double Point::GetCoordinate(int Position)
{
	if (Position >= 0 && Position < DIM_SMA)	
		return this->Coordinates[Position];

	return FLOAT_ERR;

}

double Point::GetRadius()
{
	return this->Radius;
}

PointType Point::GetType()
{
	return this->Type;
}

int Point::GetDegree()
{

	return this->Degree;

}


/***** SETTERS *****/

void Point::SetPointId(int PointId)
{
	this->PointId = PointId;
}

void Point::SetCoordinate(int Position, double Value)
{
	if (Position < DIM_SMA)	
		this->Coordinates[Position] = Value;
}

void Point::SetCoordinates(double *Coordinates)
{
	for (int i = 0; i < DIM_SMA; i++)
	{
		this->Coordinates[i] = Coordinates[i];
	}

}


void Point::SetRadius(double Radius)
{
	this->Radius = Radius;
}

void Point::SetType(PointType Type)
{
	this->Type = Type;
}

/***** OTHER METHODS *****/

void Point::PointCoordinates(double *Coordinates)
{
	for (int i = 0; i < DIM_SMA; i++)
		Coordinates[i] = this->Coordinates[i];

}


void Point::IncreaseDegree(int Value)
{

	this->Degree += Value;

}

bool Point::SameCoordinates(double Coordinates[DIM_SMA])
{

	/*if (	this->Coordinates[X_COORD] == Coordinates[X_COORD] &&
		this->Coordinates[Y_COORD] == Coordinates[Y_COORD] &&
		this->Coordinates[Z_COORD] == Coordinates[Z_COORD])*/

	if (	fabs(this->Coordinates[X_COORD] - Coordinates[X_COORD]) < EQUAL_DIST &&
		fabs(this->Coordinates[Y_COORD] - Coordinates[Y_COORD]) < EQUAL_DIST &&
		fabs(this->Coordinates[Z_COORD] - Coordinates[Z_COORD]) < EQUAL_DIST )
			return true;

	return false;
	

}


void Point::PrintPoint()
{

	cout << "Point: " << this->GetPointId() << "; Degree: " << this->Degree << "; Coordinates: (" << this->GetCoordinate(0) << ", " << this->GetCoordinate(1) << ", " << this->GetCoordinate(2) << "); Radius: " << this->GetRadius() << "; Type: ";

		switch ( this->GetType() )
		{
			case internal_point:
				cout << "internal";
				break;
			case external_point:
				cout << "external";
				break;
			case boundary_point:
				cout << "boundary";
				break;
			case atom_point:
				cout << "atom";
				break;
			case undefined_point:
				cout << "undefined";
				break;
			default:
				cout << "uninitialized";
				break;
		}
	
		cout << "\n";


}

Point* Point::PointCopy()
{
	Point *NewPoint;
	double NewCoordinates[DIM_SMA];

	for (int i=0; i < DIM_SMA; i++)
		NewCoordinates[i] = this->Coordinates[i];

	// Create the new point using the values of the older one.
	// NOTE: PointId is also copied. This may cause conflict if points are sharing Complex and should be sorted out at higher levels.
	NewPoint = new Point(this->PointId, NewCoordinates, this->Radius, this->Type);
	
	return NewPoint;

}

double Point::DistanceToPoint(Point *TargetP)
{

	double TPCoordinates[DIM_SMA] = {TargetP->GetCoordinate(X_COORD), TargetP->GetCoordinate(Y_COORD), TargetP->GetCoordinate(Z_COORD)};

	return Distance(this->Coordinates, TPCoordinates);

}

void Point::TranslatePoint(double *Coordinates)
{

	this->Coordinates[X_COORD] -= Coordinates[X_COORD];
	this->Coordinates[Y_COORD] -= Coordinates[Y_COORD];
	this->Coordinates[Z_COORD] -= Coordinates[Z_COORD];

}








