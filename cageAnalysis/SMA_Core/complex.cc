// Voro++, a 3D cell-based Voronoi library
//
// complex.cc
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015

#include "complex.hh"

Complex::Complex(Point **PointList, Edge **EdgeList, int NPoints, int NEdges, char *CName)
{

	this->PointList = PointList;
	this->EdgeList = EdgeList;

	this->NPoints = NPoints;
	this->NEdges = NEdges;

	this->CName = CName;

}

Complex::Complex(Point **PointList, Edge **EdgeList, int NPoints, int NEdges)
{

	this->PointList = PointList;
	this->EdgeList = EdgeList;

	this->NPoints = NPoints;
	this->NEdges = NEdges;

	this->CName = NULL;

}

Complex::Complex()
{

	this->PointList = NULL;
	this->EdgeList = NULL;

	this->NPoints = 0;
	this->NEdges = 0;

	this->CName = NULL;

}

Complex::~Complex() 
{

	//cout << "Complex destructor called" << endl;


	///		Erase points
	for (int i = 0; i < this->NPoints; i++)
	{
		delete this->PointList[i];
	}

	delete[] this->PointList;


	///		Erase edges
	for (int i = 0; i < this->NEdges; i++)
	{
		delete this->EdgeList[i];
	}

	delete[] this->EdgeList;


	///		Erase name
	if (this->CName != NULL)
	{
		delete[] this->CName;
	}

}

/***** *****/
Point** Complex::GetPointList()
{
	return this->PointList;
}

Edge** Complex::GetEdgeList()
{
	return this->EdgeList;
}

int Complex::GetNPoints()
{
	return this->NPoints;
}

int Complex::GetNEdges()
{
	return this->NEdges;
}

char * Complex::GetCName()
{
	return this->CName;
}


/***** *****/
void Complex::SetPointList(Point **PointList)
{
	this->PointList = PointList;
}

void Complex::SetEdgeList(Edge **EdgeList)
{
	this->EdgeList = EdgeList;
}

void Complex::SetNPoints(int NPoints)
{
	this->NPoints = NPoints;
}

void Complex::SetNEdges(int NEdges)
{
	this->NEdges = NEdges;
}

void Complex::SetCName(char *CName)
{
	this->CName = CName;
}

/***** *****/

int Complex::CalculateNeighborIdList(Point *TargetPoint, int *NeighborIdList)
{

	int i,j=0;
	Point *p1, *p2;

	//int *NeighborIdList = new int[this->NEdges];

	for (i = 0; i < this->NEdges; i++)
	{
		p1 = this->EdgeList[i]->GetPSource();
		p2 = this->EdgeList[i]->GetPTarget();

		// If any of the points has the index of TargetPoint, store the other's Id on IdList
		if (p1->GetPointId() == TargetPoint->GetPointId())
			NeighborIdList[j++] = p2->GetPointId();
		else if (p2->GetPointId() == TargetPoint->GetPointId())
			NeighborIdList[j++] = p1->GetPointId();
		
	}

	return j;

}


Point* Complex::GetPointById(int PointId)
{

	for (int i=0; i < this->NPoints; i++)
	{
		if (this->PointList[i]->GetPointId() == PointId)
		{
			return this->PointList[i];
		}
	}

	//cout << "Point not found: " << PointId << endl;

	return NULL;


}


Point* Complex::GetPointByCoordinates(double Coordinates[DIM_SMA])
{
	for (int i = 0; i < this->NPoints; i++)
		if(this->PointList[i]->SameCoordinates(Coordinates) )
			return this->PointList[i];

	return NULL;



}


/*************************/
/***** EDGE GETTTERS *****/
/*************************/

Edge* Complex::GetEdgeByPoints(Point *Point1, Point *Point2)
{

	int Id1, Id2;

	if (Point1 == NULL || Point2 == NULL)
	{
		return NULL;
	}

	for (int i = 0; i < this->NEdges; i++)
	{
		Id1 = this->EdgeList[i]->GetPSource()->GetPointId();
		Id2 = this->EdgeList[i]->GetPTarget()->GetPointId();
		if (Id1 == Point1->GetPointId() && Id2 == Point2->GetPointId())
			return this->EdgeList[i];
		if (Id1 == Point2->GetPointId() && Id2 == Point1->GetPointId())
			return this->EdgeList[i];
	}

	return NULL;

}


Edge* Complex::GetEdgeByIds(int IdSource, int IdTarget)
{

	int PSourceId, PTargetId;

	for (int i = 0; i < NEdges; i++)
	{
		PSourceId = this->EdgeList[i]->GetPSource()->GetPointId();
		PTargetId = this->EdgeList[i]->GetPTarget()->GetPointId();
		if ( 	(PSourceId == IdSource && PTargetId == IdTarget) || 
			(PSourceId == IdTarget && PTargetId == IdSource) )
			return this->EdgeList[i];
	}

	return NULL;


}

/* Complex::GetEdgeByCoordinates
 *
 * Returns a pointer to the edge connecting the two points with given coordinates. 
 * Returns NULL if there's no such edge.
 *
 */

Edge* Complex::GetEdgeByCoordinates(double Coordinates1[DIM_SMA], double Coordinates2[DIM_SMA])
{

	//double SourceCoordinates[DIM_SMA], TargetCoordinates[DIM_SMA];
	Point *SourcePoint, *TargetPoint;


	for (int i = 0; i < this->NEdges; i++)
	{
		SourcePoint = this->EdgeList[i]->GetPSource();

		TargetPoint = this->EdgeList[i]->GetPTarget();

		if (SourcePoint->SameCoordinates(Coordinates1) && TargetPoint->SameCoordinates(Coordinates2))
		{
			return this->EdgeList[i];
		}

		if (SourcePoint->SameCoordinates(Coordinates2) && TargetPoint->SameCoordinates(Coordinates1))
		{
			return this->EdgeList[i];
		}
	}


	return NULL;


}

Edge** Complex::GetEdgesByPoint(Point *TargetPoint, int &Size)
{

	int Count = 0;
	Point *PSAux, *PTAux;

	Edge **EdgeListAux = new Edge*[this->NEdges];

	for (int i = 0; i < this->NEdges; i++)
	{
		PSAux = this->EdgeList[i]->GetPSource();
		PTAux = this->EdgeList[i]->GetPTarget();
		if ((PSAux->GetPointId() == TargetPoint->GetPointId()) || (PTAux->GetPointId() == TargetPoint->GetPointId()) )
		{
			EdgeListAux[Count] = this->EdgeList[i];
			Count++;
		}	
	}

	Size = Count;

	if (Count == 0)
	{
		delete(EdgeListAux);
		return NULL;
	}

	return EdgeListAux;

}

/**** Gets the edges touching the point when the other end is of a given type ****/
Edge** Complex::GetEdgesByPointOppositeRestricted(	Point *TargetPoint, 
							PointType OppositeType,
							int &Size)
{

	int Count = 0;
	Point *PSAux, *PTAux;

	Edge **EdgeListAux = new Edge*[this->NEdges];

	for (int i = 0; i < this->NEdges; i++)
	{
		PSAux = this->EdgeList[i]->GetPSource();
		PTAux = this->EdgeList[i]->GetPTarget();

		// Check if any of the points corresponds with the point of interest
		if (PSAux->GetPointId() == TargetPoint->GetPointId())
		{
			// Check if the point is connected to a point of the type of interest
			if (PTAux->GetType() == OppositeType)
			{
				EdgeListAux[Count] = this->EdgeList[i];
				Count++;
			}
		}
		else if (PTAux->GetPointId() == TargetPoint->GetPointId())
		{
			// Check if the point is connected to a point of the type of interest
			if (PSAux->GetType() == OppositeType)
			{
				EdgeListAux[Count] = this->EdgeList[i];
				Count++;
			}
		}	
	}

	Size = Count;

	if (Count == 0)
	{
		delete(EdgeListAux);
		return NULL;
	}

	return EdgeListAux;

}

/*
 *
 */
void Complex::GetPointNeighbors(	Point *TargetPoint,
									vector <Point *> &Neighbors)
{

	for (int i = 0; i < this->NEdges; i++)
	{
		Point *PSAux, *PTAux;

		PSAux = this->EdgeList[i]->GetPSource();
		PTAux = this->EdgeList[i]->GetPTarget();

		if (PSAux->GetPointId() == TargetPoint->GetPointId())
		{
			Neighbors.push_back(PTAux);
		}

		if (PTAux->GetPointId() == TargetPoint->GetPointId())
		{
			Neighbors.push_back(PSAux);
		}

	}

}

/*
 *	Returns the minimum and maximum values for the given coordinate within the complex
 */

void Complex::GetComplexBoundary(int Coordinate, double &Min, double &Max)
{

	Max = Min = this->GetPointList()[0]->GetCoordinate(Coordinate);

	for (int i = 1; i < this->GetNPoints(); i++)
	{

		Point *TPoint = this->GetPointList()[i];

		if (TPoint->GetCoordinate(Coordinate) > Max)
		{
			Max = TPoint->GetCoordinate(Coordinate);
		}
		if (TPoint->GetCoordinate(Coordinate) < Min)
		{
			Min = TPoint->GetCoordinate(Coordinate);
		}
	}

}

void Complex::GetComplexCenter(double *Coordinates)
{

	double Min, Max;

	///		Center x coordinate
	this->GetComplexBoundary(X_COORD, Min, Max);

	Coordinates[X_COORD] = (Max + Min)/2.;

	///		Center y coordinate
	this->GetComplexBoundary(Y_COORD, Min, Max);

	Coordinates[Y_COORD] = (Max + Min)/2.;

	///		Center z coordinate
	this->GetComplexBoundary(Z_COORD, Min, Max);

	Coordinates[Z_COORD] = (Max + Min)/2.;

}

double Complex::GetComplexRadius()
{

	double Min, Max;

	///		X limits
	this->GetComplexBoundary(X_COORD, Min, Max);
	double DiffX = abs(Max-Min);

	///		Y limits
	this->GetComplexBoundary(Y_COORD, Min, Max);
	double DiffY = abs(Max-Min);

	///		Z limits
	this->GetComplexBoundary(Z_COORD, Min, Max);
	double DiffZ = abs(Max-Min);

	///		Return maximum difference divided by 2

	return max(DiffX, max(DiffY, DiffZ))/2.;


}

double Complex::GetComplexRadius(Point *P)
{

	double Radius = 0.;

	for (int i = 0; i < this->NPoints; i++)
	{
		double Dist = P->DistanceToPoint(this->PointList[i]);

		if (Dist > Radius)
		{
			Radius = Dist;
		}
	}

	return Radius;

}

/*************************/
/****** STD OUTPUT *******/
/*************************/

/* Complex::PrintComplex
 *
 * Prints complex through stdout, showing first the points and then the edges, with all the information relevant to each one.
 *
 */

void Complex::PrintComplex()
{

	Edge *EAux;

	//	PRINT POINTS
	//

	for (int i=0; i < this->NPoints; i++)
	{
		this->PointList[i]->PrintPoint();
		
	}

	//	PRINT EDGES
	//

	for (int i=0; i < this->NEdges; i++)
	{
		EAux = this->EdgeList[i];
		cout << "Edge: " << EAux->GetPSource()->GetPointId() << " - " << EAux->GetPTarget()->GetPointId() << "; Radius: " << EAux->GetChannelRadius() << "\n";
	}

	//	PRINT GENERAL DATA
	//

	cout << "Complex has " << this->NPoints << " points and " << this->NEdges << " edges." << endl;

}


void Complex::PrintComplexPoints(PointType PT)
{

	Point *PAux;

	int Count = 0;

	for (int i = 0; i < this->NPoints; i++)
	{
		PAux = this->PointList[i];
		//if (PAux->GetType() == PT || PT == undefined_point)
		if (PAux->GetType() == PT)
		{
			PAux->PrintPoint();

			Count++;
		}
	}

	cout << "There are " << Count << " points of type " << PT << "\n";

}

void Complex::PrintComplexEdges()
{
	
	for (int i = 0; i < this->NEdges; i++)
	{
		cout << "Edge: ( " << this->EdgeList[i]->GetPSource()->GetPointId() << " , " << this->EdgeList[i]->GetPTarget()->GetPointId() << " )" << " - Channel Radius " << this->EdgeList[i]->GetChannelRadius() << "\n"; 
	}

}


/* Complex::InsertPointByValue
 *
 * Brief: Inserts a point into the complex, creating the point object through the values given as arguments.
 *
 * Inputs:
 *	- PointId (opt): Identifier of the point --> WARNING! The identifier must be unique (the program won't check it)
 *	- Coordinates: doubles containing the x,y,z coordinates of the point to be inserted.
 * 	- Radius (opt): double value containing the radius of the point.
 *	- Type (opt): Enumerate containing the type of the point.
 */


void Complex::InsertPointByValue(int PointId, double Coordinates[DIM_SMA], double Radius, PointType Type)
{

	Point *NewPoint;

	if (this->PointList == NULL)
	{
		this->PointList = new Point*[COMPLEX_MAX_POINTS];
	}
	/*else
		this->PointList = (Point **) realloc ((this->NPoints+1)*sizeof(Point *));*/

	if (this->NPoints >= COMPLEX_MAX_POINTS)
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Maximum number of points exceeded \n";
		return;
	}

	//	WARNING!
	if ( this->GetPointByCoordinates(Coordinates) != NULL )
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Point already exists" << endl;
		return;
	}
	//
		
	NewPoint =  new Point(PointId, Coordinates, Radius, Type);
	this->PointList[NPoints] = NewPoint;
	this->NPoints++;

}

void Complex::InsertPointByValue(double Coordinates[DIM_SMA], double Radius, PointType Type)
{

	Point *NewPoint;

	if (this->PointList == NULL)
	{
		this->PointList = new Point*[COMPLEX_MAX_POINTS];
	}
	/*else
		this->PointList = (Point **) realloc ((this->NPoints+1)*sizeof(Point *));*/

	if (this->NPoints >= COMPLEX_MAX_POINTS)
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Maximum number of points exceeded \n";
		return;
	}

	//	WARNING!
	if ( this->GetPointByCoordinates(Coordinates) != NULL )
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Point already exists" << endl;
		return;
	}
	//
		
	NewPoint =  new Point(this->NPoints, Coordinates, Radius, Type);
	this->PointList[NPoints] = NewPoint;
	this->NPoints++;

}


void Complex::InsertPointByValue(double Coordinates[DIM_SMA])
{

	Point *NewPoint;

	if (this->PointList == NULL)
		this->PointList = new Point*[COMPLEX_MAX_POINTS];
	/*else
		this->PointList = (Point **) realloc ((this->NPoints+1)*sizeof(Point *));*/

	if (this->NPoints >= COMPLEX_MAX_POINTS)
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Maximum number of points exceeded \n";
		return;
	}

	//	WARNING!
	if ( this->GetPointByCoordinates(Coordinates) != NULL )
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Point already exists" << endl;
		return;
	}
	//
		
	NewPoint =  new Point(this->NPoints, Coordinates, 0.0, undefined_point);
	this->PointList[NPoints] = NewPoint;
	this->NPoints++;

}

void Complex::InsertPointByValue(Point *NewPoint)
{
	if (NewPoint == NULL)
	{
		//cout << "Error: trying to insert null point into complex" << endl;
		return;
	}

	if (this->NPoints >= COMPLEX_MAX_POINTS)
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Maximum number of points exceeded \n";
		return;
	}

	//	WARNING!
	if (  this->GetPointByCoordinates( NewPoint->GetCoordinates() ) != NULL  )
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Point already exists" << endl;
		delete NewPoint;
		return;
	}
	//

	if (this->PointList == NULL) {
		this->PointList = new Point*[COMPLEX_MAX_POINTS];
	}

	this->PointList[NPoints] = NewPoint;
	this->NPoints++;

}

/*
 *
 */

bool Complex::InsertPointByValue(Point *NewPoint, bool Update)
{
	if (NewPoint == NULL)
	{
		//cout << "Error: trying to insert null point into complex" << endl;
		return false;
	}

	if (this->NPoints >= COMPLEX_MAX_POINTS)
	{
		//cout << "Warning: InsertPointByValue returned without insertion - Maximum number of points exceeded \n";
		return false;
	}

	//	Update point id if required by the user
	Point *P = this->GetPointByCoordinates( NewPoint->GetCoordinates() );
	if (  P != NULL  )
	{
		//cout << "Updating point ID" << endl;
		if (Update) P->SetPointId(NewPoint->GetPointId());
		delete NewPoint;
		return false;
	}
	//

	if (this->PointList == NULL) {
		this->PointList = new Point*[COMPLEX_MAX_POINTS];
	}

	this->PointList[NPoints] = NewPoint;
	this->NPoints++;

	return true;

}

/* Complex::InsertPointByValue
 *
 * Brief: Inserts a point into the complex, creating the point object through the values given as arguments.
 *
 * Inputs:
 *	- SourceP, TargetP: Point connected by the edge
 * 	- Radius (opt): double value containing the radius of the edge.
 */

void Complex::InsertEdgeByValue(Point *SourceP, Point *TargetP, double ChannelRadius)
{


	Edge *NewEdge;

	if (this->EdgeList == NULL)
		this->EdgeList = new Edge*[COMPLEX_MAX_EDGES];

	if (SourceP == NULL || TargetP == NULL)
	{
		//cout << "Warning: InsertEdgeByValue returned without insertion - One or more points do not exist \n";
		return;
	}

	if (this->NEdges >= COMPLEX_MAX_EDGES)
	{
		//cout << "Warning: InsertEdgeByValue returned without insertion - Maximum number of edges exceeded \n";
		return;
	}

	if (this->GetEdgeByPoints(SourceP, TargetP) != NULL)
	{
		//cout << "Warning: InsertEdgeByValue returned without insertion - Edge already exists" << endl;
		return;
	}

	NewEdge = new Edge(SourceP, TargetP, ChannelRadius);
	this->EdgeList[NEdges] = NewEdge;
	this->NEdges++;

	// Increase the degree of the points
	SourceP->IncreaseDegree(1);
	TargetP->IncreaseDegree(1);

}

void Complex::InsertEdgeByValue(Point *SourceP, Point *TargetP)
{


	Edge *NewEdge;

	if (this->EdgeList == NULL)
	{
		this->EdgeList = new Edge*[COMPLEX_MAX_EDGES];
	}

	if (SourceP == NULL || TargetP == NULL)
	{
		//cout << "Warning: InsertEdgeByValue returned without insertion - One or more points do not exist \n";
		return;
	}

	if (this->NEdges >= COMPLEX_MAX_EDGES)
	{
		//cout << "Warning: InsertEdgeByValue returned without insertion - Maximum number of edges exceeded \n";
		return;
	}

	if (this->GetEdgeByPoints(SourceP, TargetP) != NULL)
	{
		//cout << "Warning: InsertEdgeByValue returned without insertion - Edge already exists" << endl;
		return;
	}

	NewEdge = new Edge(SourceP, TargetP, 0.0);
	this->EdgeList[NEdges] = NewEdge;
	this->NEdges++;

	// Increase the degree of the points
	SourceP->IncreaseDegree(1);
	TargetP->IncreaseDegree(1);

}

void Complex::InsertEdgeByValue(Edge *NewEdge)
{

	if (this->EdgeList == NULL)
	{
		this->EdgeList = new Edge*[COMPLEX_MAX_EDGES];
	}

	if (this->GetEdgeByPoints(NewEdge->GetPSource(), NewEdge->GetPTarget()) != NULL)
	{
		//cout << "Warning: InsertEdgeByValue returned without insertion - Edge already exists" << endl;
		return;
	}

	this->EdgeList[NEdges] = NewEdge;
	this->NEdges++;

}


/* Complex::CopyComplex
 *
 * Copies the complex and returns it into a new object. All the memory is copied, creating new arrays if neccessary.
 *
 */

Complex* Complex::ComplexCopy()
{

//	cout << "Copy complex " << endl;
//	if (this->CName == NULL) cout << "NULL name!" << endl;

	Complex *NewComplex;

	int NPointsCopy = this->NPoints;
	//Point **PointListCopy = new Point*[NPointsCopy];
	Point **PointListCopy = new Point*[COMPLEX_MAX_POINTS];

	int NEdgesCopy = this->NEdges;	
	//Edge **EdgeListCopy = new Edge*[NEdgesCopy];
	Edge **EdgeListCopy = new Edge*[COMPLEX_MAX_EDGES];

	int AuxPSourceId, AuxPTargetId;
	Point *AuxPSource, *AuxPTarget;

	// Copy all the points
	for (int i=0; i<NPointsCopy; i++)
	{
		PointListCopy[i] = this->PointList[i]->PointCopy();
	}

	///		Copy the name of the complex
	char *CName = NULL;

	if ((this->CName) != NULL)
	{
		CName = new char[COMPLEX_NAME_SZ];
		strncpy(CName, this->CName, COMPLEX_NAME_SZ);
	}

	// Create the new complex (edgelist is still not filled)
	NewComplex = new Complex(PointListCopy, EdgeListCopy, NPointsCopy, NEdgesCopy, CName);

	for (int i=0; i < NEdgesCopy; i++)
	{
		// Get the point ID's forming the edge from the current complex
		AuxPSourceId = this->EdgeList[i]->GetPSource()->GetPointId();
		AuxPTargetId = this->EdgeList[i]->GetPTarget()->GetPointId();

		// Get the new points addresses from the copied complex
		AuxPSource = NewComplex->GetPointById(AuxPSourceId);
		AuxPTarget = NewComplex->GetPointById(AuxPTargetId);

		// Create the new edge
		// Point it towards the adequate points
		EdgeListCopy[i] = new Edge(AuxPSource, AuxPTarget, this->EdgeList[i]->GetChannelRadius() );

	}


	// Return the copied complex
	return NewComplex;
	
}


Complex* Complex::ComplexCopy(PointType PT)
{

	Complex *NewComplex;

	int NPointsCopy = 0;
	Point **PointListCopy = new Point*[COMPLEX_MAX_POINTS];

	int NEdgesCopy = 0;
	Edge **EdgeListCopy = new Edge*[COMPLEX_MAX_EDGES];

	///		Copy the name of the complex
	char *CName = NULL;

	if ((this->CName) != NULL)
	{
		CName = new char[COMPLEX_NAME_SZ];
		strncpy(CName, this->CName, COMPLEX_NAME_SZ);
	}


	/////	FILL THE POINTS
	///		Add points that have adequate type

	for (int i=0; i<this->NPoints; i++)
	{
		if (this->PointList[i]->GetType() == PT)
		{
			PointListCopy[NPointsCopy] = this->PointList[i]->PointCopy();
			NPointsCopy++;
		}

	}

	// Create the new complex (edgelist is still not filled)
	NewComplex = new Complex(PointListCopy, EdgeListCopy, NPointsCopy, NEdgesCopy, CName);

	/////	FILL THE EDGES
	///		Add edges pointing to the points.
	//		Check if they exist

	for (int i=0; i < this->NEdges; i++)
	{
		int AuxPSourceId, AuxPTargetId;
		Point *AuxPSource, *AuxPTarget;

		// Get the point ID's forming the edge from the current complex
		AuxPSourceId = this->EdgeList[i]->GetPSource()->GetPointId();
		AuxPTargetId = this->EdgeList[i]->GetPTarget()->GetPointId();

		// Get the new points addresses from the copied complex
		AuxPSource = NewComplex->GetPointById(AuxPSourceId);
		AuxPTarget = NewComplex->GetPointById(AuxPTargetId);

		// Create the new edge
		// Point it towards the adequate points
		if (AuxPSource != NULL && AuxPTarget !=NULL) {
			EdgeListCopy[NEdgesCopy] = new Edge(AuxPSource, AuxPTarget, this->EdgeList[i]->GetChannelRadius() );
			NEdgesCopy++;
		}

	}

	NewComplex->SetNEdges(NEdgesCopy);

	return NewComplex;

}


/* Complex::PruneComplexOnce
 *
 * Copies the complex prunning all the nodes with degree 1 or less (and the edges connecting to such nodes. 
 *
 */

Complex * Complex::PruneComplexOnce()
{

//	cout << "Prune complex once" << endl;
//	if (this->CName == NULL) cout << "NULL name!" << endl;

	Complex *NewComplex;

	int NPointsCopy = 0;
	//Point **PointListCopy = new Point*[NPointsCopy];
	Point **PointListCopy = new Point*[COMPLEX_MAX_POINTS];

	int NEdgesCopy = 0;	
	//Edge **EdgeListCopy = new Edge*[NEdgesCopy];
	Edge **EdgeListCopy = new Edge*[COMPLEX_MAX_EDGES];

	int AuxPSourceId, AuxPTargetId;
	Point *AuxPSource, *AuxPTarget;

	///		Copy the name of the complex
	char *CName = NULL;

	if (this->CName != NULL)
	{
		CName = new char[COMPLEX_NAME_SZ];
		strncpy(CName, this->CName, COMPLEX_NAME_SZ);
	}

	// Copy all the points
	for (int i = 0; i < this->NPoints; i++)
	{
		if (this->PointList[i]->GetDegree() > 1)
		{
			PointListCopy[NPointsCopy] = this->PointList[i]->PointCopy();
			NPointsCopy++;
		}
	}

	// Copy all the edges

	// Create the new complex (edgelist is still not filled)
	NewComplex = new Complex(PointListCopy, EdgeListCopy, NPointsCopy, NEdgesCopy, CName);

	for (int i = 0; i < this->NEdges; i++)
	{
		// Get the point ID's forming the edge from the current complex
		AuxPSourceId = this->EdgeList[i]->GetPSource()->GetPointId();
		AuxPTargetId = this->EdgeList[i]->GetPTarget()->GetPointId();

		// Get the new points addresses from the copied complex 
		AuxPSource = NewComplex->GetPointById(AuxPSourceId);
		AuxPTarget = NewComplex->GetPointById(AuxPTargetId);

		// Check if the points have been previously inserted (otherwise the edge shouldn't be add)		
		if (AuxPSource != NULL && AuxPTarget != NULL)
		{

			// Create the new edge
			// Point it towards the adequate points
			//EdgeListCopy[NEdgesCopy] = new Edge(AuxPSource, AuxPTarget, this->EdgeList[i]->GetChannelRadius() );
			NewComplex->InsertEdgeByValue(AuxPSource, AuxPTarget);
			NEdgesCopy++;

		}
	}
	NewComplex->SetNEdges(NEdgesCopy);

	// Return the copied complex
	return NewComplex;



}



bool Complex::Pruned()
{

	if (this->NPoints == 0)
		return true;

	for (int i = 0; i < this->NPoints; i++)
		if (this->PointList[i]->GetDegree() <= 1)
			return false;

	return true;



}


Complex* Complex::PruneComplex()
{

	if (this->Pruned())
	{
		return this->ComplexCopy();
	}

	Complex *NewComplex, *SwapComplex;	
	int Count = 0;

	bool Pruned = false;

	while (!Pruned)
	{
		if (Count == 0)	
		{
			SwapComplex = this->PruneComplexOnce();
		}
		else 
		{
			SwapComplex = NewComplex->PruneComplexOnce();
			delete NewComplex;
		}		
		NewComplex = SwapComplex;

		Pruned = NewComplex->Pruned();

		//cout << "New: " << NewComplex->GetNPoints() << " - this " << this->NPoints << endl;
	
		if (Count > this->NPoints)
		{
			//cout << "Error pruning complex: maximum number of possible prune steps surpassed" << endl;
			delete NewComplex;
			return NULL;
		}
		Count++;
	}

	if (NewComplex->GetNPoints() == 0 || NewComplex->GetNEdges() == 0)
	{
		delete NewComplex;
		return NULL;
	}

	NewComplex->Normalize();

	return NewComplex;

}

/*
 * 	Redefines complex keeping the topology, with changing indices from 0 to NPoints-1
 */

void Complex::Normalize()
{

	for (int i = 0; i < this->NPoints; i++)
	{
		this->PointList[i]->SetPointId(i);
	}

}

void Complex::Externalize()
{
	for (int i = 0; i < this->NPoints; i++)
	{
		Point *P0 = this->PointList[i];

		if (P0->GetType() == internal_point)
		{
			P0->SetType(external_point);
		}
	}
}

/*
 *
 */

/*void Complex::DestroyComplex()
{

	for (int i=0; i < this->NEdges; i++)
		delete this->EdgeList[i];

	for (int i=0; i < this->NPoints; i++)
		delete this->PointList[i];

	delete this->EdgeList;
	delete this->PointList;

	delete this; 


}*/







// COMPLEX CC






