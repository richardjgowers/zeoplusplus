/*
 * internal_external_class.cc
 *
 *  Created on: Feb 27, 2017
 *      Author: ismael
 */

#include "internal_external_class.hh"

/*
 *	Check if the given point is inside the given alpha shape
 */
bool IEClas_AlphaShape_ClasPoint(	Point *TargetPoint,
									Alpha_shape_3 &AlphaShape,
									double Alpha)
{

	list <Cell_handle> CellList;
	list <Cell_handle>::const_iterator CellIt;

	// Retrieve all the internal cells of the alpha shape
	AlphaShape.get_alpha_shape_cells(back_inserter(CellList), Alpha_shape_3::INTERIOR, Alpha);

	double PointCoords[DIM_SMA];
	TargetPoint->GetCoordinates(PointCoords);

	// Iterate over all the cells, checking if the node is inside any of them
	for (CellIt = CellList.begin(); CellIt != CellList.end(); CellIt++)
	{
		if ( PointInsideCell(PointCoords, *CellIt) )
		{
			return true;
		}
	} // END FOR (Alpha-shape cell list)

	return false;

}


/*
 * 	IEClas_AlphaShape_AtomsIntSRVoro
 *
 * 	Given a list of points, checks if they are inside the molecule, according to the following criterion:
 * 		Points are inside if they are inside the alpha-shape formed by the set of atoms and the internal voronoi
 * 		nodes computed with help of the shadow rate.
 *
 * 	If Reclassify TRUE, points in list are reclassified as internal points.
 *
 * 	Returns the number of points that have been considered internal.
 *
 */

int IEClas_AlphaShape_AtomsIntSRVoro(	list <Point *> &PointList,
										MoleculeInfo *Molecule,
										bool Reclassify)
{

	/////	COMPUTE THE ALPHA-SHAPE OF THE MOLECULE
	///		Add molecule's atoms and Voronoi's internal nodes
	//		Alpha = Optimal

	list<Weighted_point> UnweightedPList;

	///		ADD INTERNAL VORONOI NODES

	Complex *Voro = Molecule->GetVoronoiGraph();

	int InternalCount = 0;

	for (int i = 0; i < Voro->GetNPoints(); i++)
	{
		Point *PAux;

		PAux = Voro->GetPointList()[i];

		// Check that the point is of desired type
		if (PAux->GetType() == internal_point)
		{
			// Insert point into the unweighted point list (weighted alpha-shapes)
			UnweightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), 0.0	  ) );
			InternalCount++;

		}
	}

	if (InternalCount == 0)
	{
		return 0;
	}


	///		ADD THE ATOMS

	Complex *Chemical = Molecule->GetChemical();

	// Insert points into delaunay triangulaton object
	for (int i = 0; i < Chemical->GetNPoints(); i++)
	{
		Point *PAux;

		PAux = Chemical->GetPointList()[i];

		// Check that the point is of desired type
		if (PAux->GetType() == atom_point)
		{
			// Insert point into the unweighted point list (weighted alpha-shapes)
			UnweightedPList.push_back( Weighted_point(  Bare_point(PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2)), 0.0	  ) );

		}

	}


	Alpha_shape_3 MoleculeAlphaShape(UnweightedPList.begin(), UnweightedPList.end(), 0, Alpha_shape_3::GENERAL);

	double Alpha;

	Alpha_iterator 	opt = MoleculeAlphaShape.find_optimal_alpha(1);

	// If the alpha iterator is placed at the end, no optimal alpha was found --> The method failed to identify internal cells (there are none)
	if (opt == MoleculeAlphaShape.alpha_end())
	{
		//cout << "Ups" << endl;
		return false;
	}

	//	Fix alpha to optimal

	Alpha = *opt;
	Alpha += 0.1*Alpha;


	/////	CHECK POINTS TO BE INTERNAL OR EXTERNAL
	///		Count for total internal points
	//		Reclassify in depends on flag

	int NInternalPoints = 0;

	list <Point *>::iterator PIt;

	for (PIt = PointList.begin(); PIt != PointList.end(); PIt++)
	{
		Point *TargetPoint = *PIt;

		if (IEClas_AlphaShape_ClasPoint(TargetPoint, MoleculeAlphaShape, Alpha) )
		{
			NInternalPoints++;

			if (Reclassify)
			{
				TargetPoint->SetType(internal_point);
			}
		}

	}

	return NInternalPoints;

}


/*	IEClas_MoleculesOverlap
 *
 *	Checks if two molecules overlap by testing if any of the atoms of the second are inside the first. Based on shadow ratio to decide whether points are inside.
 */

bool IEClas_MoleculesOverlap(	Complex *Chemical1,
								Complex *Chemical2)
{

	//		Load Chemical2 points into list

	list <Point *> PointList;

	for (int i = 0; i < Chemical2->GetNPoints(); i++)
	{
		PointList.push_back(Chemical2->GetPointList()[i]);
	}

	MoleculeInfo *Chem1MI = new MoleculeInfo();

	Chem1MI->SetChemical(Chemical1);

	int NInternal = PointListShadowTest (	Chem1MI, PointList, false);

	if (NInternal == 0)
	{
		return false;
	}

	return true;

}

/*
 *
 */

int IEClas_PointsInMolecule( 	Complex *Chemical,
								list <Point *> &PointList)
{

	MoleculeInfo *ChemMI = new MoleculeInfo();

	Complex *PrunedChemical = Chemical->PruneComplex();

	ChemMI->SetChemical(PrunedChemical);

	int NInternal = PointListShadowTest( ChemMI, PointList, true);

	delete ChemMI;

	return NInternal;



}
