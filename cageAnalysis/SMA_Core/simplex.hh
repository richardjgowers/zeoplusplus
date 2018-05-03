/*
 * simplex.hh
 *
 *  Created on: Mar 18, 2016
 *      Author: ismael.gomez
 */

#ifndef SIMPLEX_HH_
#define SIMPLEX_HH_


#include <vector>
#include <iostream>

#include "point.hh"

using namespace std;

class Simplex
{

public:

	///		Constructors and destructors
	Simplex(int Dim, vector <Point *> *SimplexPoints);
	Simplex();
	~Simplex();

	///		Setters
	void SetDim(int Dim);

	///		Getters
	int GetDim();
	vector <Point *> *GetSimplexPoints();
	double GetSurface();

	///		Other methods
	void InsertPoint (Point *SimplexPoint);
	Point *ReadPoint(int Position);

	void PrintSimplex();


private:

	int Dim;

	vector <Point *> *SimplexPoints;

};

#endif /* SIMPLEX_HH_ */
