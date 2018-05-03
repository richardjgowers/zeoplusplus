/*
 * simplex.cc
 *
 *  Created on: Mar 18, 2016
 *      Author: ismael.gomez
 */

#include "simplex.hh"


Simplex::Simplex(int Dim, vector <Point *> *SimplexPoints)
{

	this->Dim = Dim;

	this->SimplexPoints = SimplexPoints;

}


Simplex::Simplex()
{

	this->Dim = 0;

	this->SimplexPoints = new vector<Point *>();

}

Simplex::~Simplex()
{
	//cout << "Simplex destructor called" << endl;

	delete this->SimplexPoints;
}


void Simplex::InsertPoint (Point *SimplexPoint)
{

	this->SimplexPoints->push_back(SimplexPoint);

}

void Simplex::SetDim(int Dim)
{

	this->Dim = Dim;

}

int Simplex::GetDim()
{

	return this->Dim;

}

vector <Point *> * Simplex::GetSimplexPoints()
{
	return this->SimplexPoints;
}

Point * Simplex::ReadPoint(int Position)
{

	/*if (Position > (this->Dim)+1 )
	{
		cout << "Error: trying to read point for simplex out of its dimension" << endl;
		return NULL;
	}*/


	return this->SimplexPoints->at(Position);

}

double Simplex::GetSurface()
{
	Point 	*P1 = this->SimplexPoints->at(0),
			*P2 = this->SimplexPoints->at(1),
			*P3 = this->SimplexPoints->at(2);

	double 	P1Coords[DIM_SMA] = {P1->GetCoordinate(X_COORD), P1->GetCoordinate(Y_COORD), P1->GetCoordinate(Z_COORD)},
			P2Coords[DIM_SMA] = {P2->GetCoordinate(X_COORD), P2->GetCoordinate(Y_COORD), P2->GetCoordinate(Z_COORD)},
			P3Coords[DIM_SMA] = {P3->GetCoordinate(X_COORD), P3->GetCoordinate(Y_COORD), P3->GetCoordinate(Z_COORD)};

	return TriangleArea(P1Coords, P2Coords, P3Coords);
}

/**********/
void Simplex::PrintSimplex()
{

	vector <Point *>::iterator PIt;

	for (PIt = this->SimplexPoints->begin(); PIt != this->SimplexPoints->end(); PIt++)
	{
		Point *PAux = *(PIt);
		PAux->PrintPoint();
	}

}
