/*
 * window_PER.cc
 *
 *  Created on: 03/07/2017
 *      Author: Ismael
 */

#include "window_PER.hh"


/*	ClosestBond
 *
 * 	Computes and returns the closest bond to a given edge.
 *
 */

Edge *ClosestBond(	Complex *Chemical,
					Edge *BoundaryE)
{

	double MinDist = 1000.;

	double ECoords[DIM_SMA];

	BoundaryE->GetMiddlePoint(ECoords);

	Edge *BondCandidate = NULL;

	for (int i = 0; i < Chemical->GetNEdges(); i++)
	{

		Edge *Bond = Chemical->GetEdgeList()[i];

		double BCoords[DIM_SMA];

		Bond->GetMiddlePoint(BCoords);

		double Dist = Distance(BCoords, ECoords);

		if (Dist < MinDist)
		{
			BondCandidate = Bond;
			MinDist = Dist;
		}

	}

	return BondCandidate;

}

/*
 *
 */

void ComponentBoundaryEdges( Window *CC,
							 list <Edge *> &BoundEdges)
{

	///	Run over all the edges

	for (int i = 0; i < CC->GetNEdges(); i++)
	{

		Edge *TargetE = CC->GetEdgeList()[i];

		BoundEdges.push_back(TargetE);

		/*double 	SCoords[DIM_SMA],
				TCoords[DIM_SMA];

		TargetE->GetPSource()->GetCoordinates(SCoords);
		TargetE->GetPTarget()->GetCoordinates(TCoords);

		int SimplicesCounter = 0;

		///	Run over simplices list and see if it belongs to only one
		list <Simplex *>::iterator SIt;

		for (SIt = CC->GetSimplices()->begin(); SIt != CC->GetSimplices()->end(); SIt++)
		{
			int SCounter = 0;

			Simplex *BS = *SIt;

			for (int j = 0; j < BS->GetSimplexPoints()->size(); j++)
			{
				Point *SP = BS->GetSimplexPoints()->at(j);

				if (SP->SameCoordinates(SCoords) || SP->SameCoordinates(TCoords))
				{
					SCounter++;
				}

			}

			//	 If 2 or more points of simplex coincide with edge points, edge is part of the simplex
			if (SCounter == 2)
			{
				SimplicesCounter++;
			}

		}

		///		If only one simplex was found to contain the edge, add the edge to the boundary list
		if (SimplicesCounter == 1)
		{
			BoundEdges.push_back(TargetE);
		}
	*/

	}

}

/*
 *
 */

void PointToBoundaryEdges (	Point *RefP,
							list <Edge *> &BoundaryEdges,
							list <Edge *> &BTCEdges)
{

	list <Edge *>::iterator CBEIt;

	for (CBEIt = BoundaryEdges.begin(); CBEIt != BoundaryEdges.end(); CBEIt++)
	{

		Edge *BE = *CBEIt;

		Edge *BTCE = new Edge(BE->GetPSource(), RefP, 0.);

		BTCEdges.push_back(BTCE);

	}

}

/*	ComputeWindowsPER
 *
 * 	Retrieve molecule's windows using the same approach as the one used for computing the internal pores. Each component of the external grid (after erasing triangles) is assumed to be associated to one of the molecule's windows.
 * 	Closest atoms and bonds to a given component must form a chemical window entering the molecule. This criterion is used to retrieve a list of windows, being each window a connected set of bonds and atoms. Windows are directly inserted
 * 	into the MoleculeInfo object provided as an argument.
 */

bool ComputeWindowsPER( MoleculeInfo *Molecule)
{

	///		Get the biggest pore (reference for window computing)

	Point *Pore = Molecule->GetGreatestVoroNode(internal_point);

	double Maximum = 0.;

	///			Construct characteristic for the pore

	double Center[DIM_SMA];

	Pore->GetCoordinates(Center);
	double Radius = Molecule->GetChemical()->GetComplexRadius(Pore);

	Radius *= 1.5;

	const int NDivisions = 4;

	Window *Boundary = IcosahedronSphericGrid(Center, Radius, NDivisions);

	double GridSurface;

	Window *Characteristic = NodeCharacteristic(Pore, Molecule, Boundary, GridSurface);

	vector <Window *> Components;

	ConnectedComponents(Characteristic, Components);


	/////	PER-BASED WINDOW RECONSTRUCTION
	///		Window computing using tools from pore exposure. Each connected component from the shadow is considered to have 1 window associated
	//		(1) 	Identify boundary edges for CC.
	//		(2)		For each boundary edge, identify molecule's closest bond.
	//		(3)		Reconstruct the window with given resolution

	for (int i = 0; i < Components.size(); i++)
	{

		Window *CC = Components.at(i);

		///		Create the new window
		Window *MolWindow = new Window();

		/////	BOUNDARY EDGES

		list <Edge *> BoundaryEdges;

		ComponentBoundaryEdges(CC, BoundaryEdges);

		list <Edge *>::iterator CBEIt;

		for (CBEIt = BoundaryEdges.begin(); CBEIt != BoundaryEdges.end(); CBEIt++)
		{

			Edge *BE = *CBEIt;

			Edge *WE = ClosestBond(Molecule->GetChemical(), BE);

			//	Insert points

			double 	SCoords[DIM_SMA],
					TCoords[DIM_SMA];

			WE->GetPSource()->GetCoordinates(SCoords);
			WE->GetPTarget()->GetCoordinates(TCoords);

			MolWindow->InsertPointByValue(SCoords);
			MolWindow->InsertPointByValue(TCoords);

			//	Insert edge
			Point 	*Source = MolWindow->GetPointByCoordinates(SCoords),
					*Target = MolWindow->GetPointByCoordinates(TCoords);

			if (Source != NULL && Target != NULL)
			{
				MolWindow->InsertEdgeByValue(Source, Target);
			}


		}

		/////	BCT EDGES
		list <Edge *> BCTEdges;

		PointToBoundaryEdges(Pore, BoundaryEdges, BCTEdges);

		list <Edge *>::iterator BCTIt;

		for (BCTIt = BCTEdges.begin(); BCTIt != BCTEdges.end(); BCTIt++)
		{

			Edge *BE = *BCTIt;

			Edge *WE = ClosestBond(Molecule->GetChemical(), BE);

			//	Insert points

			double 	SCoords[DIM_SMA],
					TCoords[DIM_SMA];

			WE->GetPSource()->GetCoordinates(SCoords);
			WE->GetPTarget()->GetCoordinates(TCoords);

			MolWindow->InsertPointByValue(SCoords);
			MolWindow->InsertPointByValue(TCoords);

			//	Insert edge
			Point 	*Source = MolWindow->GetPointByCoordinates(SCoords),
					*Target = MolWindow->GetPointByCoordinates(TCoords);

			if (Source != NULL && Target != NULL)
			{
				MolWindow->InsertEdgeByValue(Source, Target);
			}

			// TODO - ERASE BCT EDGES (but not edges from previous loop! these are boundary edges erased elsewhere)


		}


		MolWindow->Normalize();

		if (MolWindow->GetNPoints() <= 6)
		{
			delete MolWindow;
		}
		else
		{
			Molecule->InsertWindow(MolWindow);
		}


	}


	return true;
}
