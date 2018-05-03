/*
 * volume.cc
 *
 *  Created on: Mar 3, 2016
 *      Author: ismael.gomez
 */

#include "volume.hh"


/* LargestInternalSphere
 *	Searches for the largest internal node (in terms of radius), returning that value.
 *
 * Inputs:
 * 		MolVoronoi - Voronoi tessellation of the molecule of interest.
 *
 * Outputs:
 * 		Radius of the Voronoi node with maximum radius (biggest cell's index).
 */
double LargestInternalSphere (	Complex *MolVoronoi)
{

	int i;
	double Maximum = 0.0;

	for (i=0; i < MolVoronoi->GetNPoints(); i++)
	{
		if (MolVoronoi->GetPointList()[i]->GetType() == internal_point)
		{
			if (MolVoronoi->GetPointList()[i]->GetRadius() > Maximum)
			{
				Maximum = MolVoronoi->GetPointList()[i]->GetRadius();
			}
		}
	}

	//double SphereVolume = 3.*PI/4.*pow(Maximum, 3.);

	return Maximum;
}


/* LargestInternalNodeIndex
 *
 * Inputs:
 * 		MolVoronoi - Voronoi tessellation of the molecule of interest.
 *
 * Outputs:
 * 		Index of the Voronoi node with maximum radius (biggest cell's index).
 */
int LargestInternalNodeIndex(	Complex *MolVoronoi)
{
	int i;
	int IndexMaximum = 0;
	double Maximum = 0.0;

	for (i=0; i < MolVoronoi->GetNPoints(); i++)
		if (MolVoronoi->GetPointList()[i]->GetType() == internal_point)
			if (MolVoronoi->GetPointList()[i]->GetRadius() > Maximum)
			{
				IndexMaximum = i;
				Maximum = MolVoronoi->GetPointList()[i]->GetRadius();
			}

	return IndexMaximum;

}


/* LargestInternalNodeCoordinates
 * 	Looks for the largest internal node and loads its coordinates into the vector passed as an argument
 *
 * Inputs:
 * 		MolVoronoi - Voronoi tessellation of the molecule of interest.
 * 		Coordinates - Vector for storing the coordinates.
 *
 * Outputs:
 * 		(In Coordinates vector) Coordinates of the largest node.
 */
void LargestInternalNodeCoordinates(	Complex *MolVoronoi,
										double *Coordinates)
{
	int i;
	Point *LargestP, *AuxP;
	double Maximum = 0.0;

	for (i=0; i < MolVoronoi->GetNPoints(); i++)
	{
		AuxP = MolVoronoi->GetPointList()[i];
		if (AuxP->GetType() == internal_point)
		{
			if (AuxP->GetRadius() > Maximum)
			{
				LargestP = AuxP;
				Maximum = AuxP->GetRadius();
			}
		}
	}

	// Get the coordinates
	Coordinates[X_COORD] = AuxP->GetCoordinate(X_COORD);
	Coordinates[Y_COORD] = AuxP->GetCoordinate(Y_COORD);
	Coordinates[Z_COORD] = AuxP->GetCoordinate(Z_COORD);

}

/* VoroTotalVolume
 * 	Computes total molecule's volume based on Voronoi tessellation (internal cells are decided using the alpha-shape based criterion implemented at internal_cells.cc)
 *
 * Inputs:
 * 		MolVoronoi - Voronoi tessellation of the molecule of interest.
 *
 * Outputs:
 * 		Total volume trapped by the molecule.
 */

double VoroTotalVolume(	Complex *MolVoronoi)
{
	int i;
	double TotalVolume = 0.0;

	/////	COMPUTE VOLUMES AND SUM
	///		Check every point and see if its internal.
	//		Calculate the volume of the maximum sphere fitting in that cell
	//		Add that volume to the total

	//MolVoronoi->PrintComplex();

	for (i=0; i < MolVoronoi->GetNPoints(); i++)
	{
		if (MolVoronoi->GetPointList()[i]->GetType() == internal_point)
		{
			// Get the radius of the cell
			double CellRadius = MolVoronoi->GetPointList()[i]->GetRadius();

			// Add the volume of the maximum sphere fitting into the cell
			TotalVolume += (4/3)*PI*pow(CellRadius, 3.);
		}
	}

	return TotalVolume;

}







// VOLUME_CC
