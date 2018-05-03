/*
 * surface_area.cc
 *
 *  Created on: Feb 27, 2017
 *      Author: ismael
 */

#include "surface_area.hh"


/*
 *	SurfaceArea_RandomPoint
 */

Point *SurfaceArea_RandomPoint ( Point *Atom)
{

	Point *RandomPoint = NULL;

	/////	RANDOM POINT
	///		Point with fixed radius, variable phi and theta (in spheric coordinates).
	//		Transform into cartesian just after

	///		Generate random coords (except radius)
	double Radius = Atom->GetRadius();

	double V = rand() % 100;
	V /= 100.;
	double Phi = acos(2*V - 1);

	double U = rand() % 100;
	U /= 100.;
	double Theta = 2.*PI*U;

	double 	CartesianCoords[DIM_SMA],
			SphericCoords[DIM_SMA];

	SphericCoords[R_COORD] = Radius;
	SphericCoords[PHI_COORD] = Phi;
	SphericCoords[THETA_COORD] = Theta;

	//cout << "Phi: " << Phi << " Theta: " << Theta << endl;

	///		Generate cartesian coordinates and insert into point

	SphericToCartesian(CartesianCoords, SphericCoords);

	RandomPoint = new Point();

	RandomPoint->SetCoordinate(X_COORD, CartesianCoords[X_COORD] + Atom->GetCoordinate(X_COORD));
	RandomPoint->SetCoordinate(Y_COORD, CartesianCoords[Y_COORD] + Atom->GetCoordinate(Y_COORD));
	RandomPoint->SetCoordinate(Z_COORD, CartesianCoords[Z_COORD] + Atom->GetCoordinate(Z_COORD));

	return RandomPoint;

}

/*
 *	SurfaceArea_MonteCarlo
 *
 * 	Computes molecule's surface area via Monte Carlo sampling of points around atom's surface (random positions).
 *
 * 	Sampled points on atoms' surfaces are checked for being internal (fast alpha-shape based check) and partial against
 * 		total values are counted for surface area value.
 */

double SurfaceArea_MonteCarlo( MoleculeInfo *Molecule)
{
	//int SampleSize = SAMPLES_PER_ATOM*Molecule->GetChemical()->GetNPoints();

	/////	SURFACE AREA BY MONTE CARLO
	///		For each atom, randomly generate points at its surface
	//		Count internal points from previous group and store total sample size
	//		SurfaceArea = InternalSamplePoints/TotalSamplePoints

	Complex *Chemical = Molecule->GetChemical();

	list <Point *>	SamplePointsList;

	for (int i = 0; i < Chemical->GetNPoints(); i++)
	{

		Point *Atom = Chemical->GetPointList()[i];

		//Atom->PrintPoint();

		for (int j = 0; j < SAMPLES_PER_ATOM; j++)
		{
			Point *PAux = SurfaceArea_RandomPoint(Atom);
			PAux->SetRadius(Atom->GetRadius());

			SamplePointsList.push_back(PAux);

			//if (j == 0) cout << Atom->DistanceToPoint(PAux) << endl;
		}

	}

	///		Count for positive points

	int NInternalPoints = IEClas_AlphaShape_AtomsIntSRVoro(SamplePointsList, Molecule, true);

	///		Sum each point contribution to area
	//		Point contribution PC = PI*R^2/SAMPLES_PER_ATOM

	double MC_Surface = 0.0;

	if (NInternalPoints == 0)
	{
		return 0.0;
	}

	list <Point *>::iterator PLIt;

	for (PLIt = SamplePointsList.begin(); PLIt != SamplePointsList.end(); PLIt++)
	{
		Point *P = *PLIt;

		if (P->GetType() == internal_point)
		{
			double R = P->GetRadius();
			MC_Surface += PI*pow(R, 2.)/SAMPLES_PER_ATOM;
		}

		//	Free the point once its been accounted for
		delete P;
	}

//	cout << "Internal surface: " << MC_Surface << endl;

	///		Free points
	//		TODO

	return MC_Surface;

}
