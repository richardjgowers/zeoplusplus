/*
 * window.cc
 *
 *  Created on: Mar 7, 2016
 *      Author: ismael.gomez
 */

#include "window.hh"

Window::Window(Point **PointList, Edge **EdgeList, int NPoints, int NEdges, double Surface, list <Simplex *> *Simplices) : Complex (PointList, EdgeList, NPoints, NEdges)
{

	this->Surface = Surface;

	this->Simplices = Simplices;

	//this->Simplices = new list<Simplex *>();

}

Window::Window() : Complex()
{

	this->Surface = 0.0;

	this->Simplices = new list<Simplex *>();

}

Window::~Window()
{

	//cout << "Window destructor called" << endl;

	//this->Simplices->clear();
	list <Simplex *>::iterator SIt;

	for (SIt = this->Simplices->begin(); SIt != this->Simplices->end(); SIt++)
	{
		Simplex *S = *SIt;

		delete S;

	}

	delete this->Simplices;

	//delete this;

}

void Window::SetSurface(double Surface)
{
	this->Surface = Surface;
}

double Window::GetSurface()
{
	if (this->Surface == 0.0 && this->Simplices->size() > 0)
	{
		// Compute surface
	}

	return this->Surface;
}

void Window::PrintWindow()
{

	this->PrintComplex();

	cout << "Window surface is: " << this->Surface << endl;

}

void Window::InsertSimplex(Simplex *NewSimplex)
{

	this->Surface += NewSimplex->GetSurface();

	this->Simplices->push_back(NewSimplex);

}

list <Simplex *>::iterator Window::EraseSimplex(list <Simplex *>::iterator Target)
{
	Simplex *TargetSimplex = *Target;

	list <Simplex *>::iterator SIt;

	this->Surface -= TargetSimplex->GetSurface();

	SIt = this->Simplices->erase(Target);
	SIt--;

	delete TargetSimplex;

	return SIt;

}

list <Simplex *> * Window::GetSimplices()
{

	return this->Simplices;

}


bool Window::SameWindow(Window *WindowT)
{

	if (this->GetNPoints() != WindowT->GetNPoints())
	{
		return false;
	}

	for (int i = 0; i < WindowT->GetNPoints(); i++)
	{
		Point *WPoint = WindowT->GetPointList()[i];

		double Coordinates[DIM_SMA] = { WPoint->GetCoordinate(X_COORD), WPoint->GetCoordinate(Y_COORD), WPoint->GetCoordinate(Z_COORD) };

		// Try to get the point from this. If is not contained, windows are not the same.
		if (this->GetPointByCoordinates(Coordinates) == NULL)
		{
			return false;
		}
	}

	// If all the points of target window were found in this, return true.
	// Note: Number of points was checked and repeated pointes are not allowed in complex.

	//cout << "Windows are the same!" << endl;

	return true;

}

/*
 * 	Recomputes entire surface area for window (necessary in some cases)
 */
void Window::SurfaceAreaUpdate()
{
	list <Simplex *>::iterator SIt;

	double Surface = 0.;

	for (SIt = this->Simplices->begin(); SIt != this->Simplices->end(); SIt++)
	{
		Simplex *S = *SIt;

		Surface += S->GetSurface();
	}

	this->Surface = Surface;
}


Complex * Window::JointSimplices()
{

	Complex *JointSimplices = new Complex();

	list <Simplex *>::iterator SIt;

	for (SIt = this->Simplices->begin(); SIt != this->Simplices->end(); SIt++)
	{
		Simplex *SAux = *SIt;

		/////	INSERT SIMPLEX POINTS
		///		Assume that simplices are of dimension 2 (in this case it fits)

		Point *PAux1 = SAux->GetSimplexPoints()->at(0);
		double Coordinates1[DIM_SMA] = {PAux1->GetCoordinate(X_COORD), PAux1->GetCoordinate(Y_COORD), PAux1->GetCoordinate(Z_COORD)};
		if (JointSimplices->GetPointByCoordinates(Coordinates1) == NULL)
		{
			JointSimplices->InsertPointByValue(Coordinates1);
		}

		Point *PAux2 = SAux->GetSimplexPoints()->at(1);
		double Coordinates2[DIM_SMA] = {PAux2->GetCoordinate(X_COORD), PAux2->GetCoordinate(Y_COORD), PAux2->GetCoordinate(Z_COORD)};
		if (JointSimplices->GetPointByCoordinates(Coordinates2) == NULL)
		{
			JointSimplices->InsertPointByValue(Coordinates2);
		}

		Point *PAux3 = SAux->GetSimplexPoints()->at(2);
		double Coordinates3[DIM_SMA] = {PAux3->GetCoordinate(X_COORD), PAux3->GetCoordinate(Y_COORD), PAux3->GetCoordinate(Z_COORD)};
		if (JointSimplices->GetPointByCoordinates(Coordinates3) == NULL)
		{
			JointSimplices->InsertPointByValue(Coordinates3);
		}

		//cout << "Facet points are " << PAux1->GetPointId() << " " << PAux2->GetPointId() << " " << PAux3->GetPointId() << endl;

		/////	INSERT SIMPLEX EDGES
		///

		Point 	*JP1 = JointSimplices->GetPointByCoordinates(Coordinates1),
				*JP2 = JointSimplices->GetPointByCoordinates(Coordinates2),
				*JP3 = JointSimplices->GetPointByCoordinates(Coordinates3);


		if (JointSimplices->GetEdgeByPoints(JP1, JP2) == NULL)
		{
			JointSimplices->InsertEdgeByValue(JP1, JP2);
		}

		if (JointSimplices->GetEdgeByPoints(JP1, JP3) == NULL)
		{
			JointSimplices->InsertEdgeByValue(JP1, JP3);
		}

		if (JointSimplices->GetEdgeByPoints(JP2, JP3) == NULL)
		{
			JointSimplices->InsertEdgeByValue(JP2, JP3);
		}

	}

	return JointSimplices;
}

///////// 	WARNING!!!!!
//////		NOT WORKING PROPERLY!!! - FIX BEFORE USE
/*Window * Window::PruneWindowOnce()
{

	Window *NewWindow;

	int NPointsCopy = 0;
	//Point **PointListCopy = new Point*[NPointsCopy];
	Point **PointListCopy = new Point*[COMPLEX_MAX_POINTS];

	int NEdgesCopy = 0;
	//Edge **EdgeListCopy = new Edge*[NEdgesCopy];
	Edge **EdgeListCopy = new Edge*[COMPLEX_MAX_EDGES];

	int AuxPSourceId, AuxPTargetId;
	Point *AuxPSource, *AuxPTarget;

	// Copy all the points
	for (int i = 0; i < this->GetNPoints(); i++)
	{
		if (this->GetPointList()[i]->GetDegree() > 1)
		{
			PointListCopy[NPointsCopy] = this->GetPointList()[i]->PointCopy();
			NPointsCopy++;
		}
	}

	// Copy all the edges

	double Surface = 0.;

	// Create the new complex (edgelist is still not filled)
	NewWindow = new Window(PointListCopy, EdgeListCopy, NPointsCopy, NEdgesCopy, Surface);

	for (int i = 0; i < this->GetNEdges(); i++)
	{
		// Get the point ID's forming the edge from the current complex
		AuxPSourceId = this->GetEdgeList()[i]->GetPSource()->GetPointId();
		AuxPTargetId = this->GetEdgeList()[i]->GetPTarget()->GetPointId();

		// Get the new points addresses from the copied complex
		AuxPSource = NewWindow->GetPointById(AuxPSourceId);
		AuxPTarget = NewWindow->GetPointById(AuxPTargetId);

		// Check if the points have been previously inserted (otherwise the edge shouldn't be add)
		if (AuxPSource != NULL && AuxPTarget != NULL)
		{

			// Create the new edge
			// Point it towards the adequate points
			//EdgeListCopy[NEdgesCopy] = new Edge(AuxPSource, AuxPTarget, this->EdgeList[i]->GetChannelRadius() );
			NewWindow->InsertEdgeByValue(AuxPSource, AuxPTarget);
			NEdgesCopy++;

		}
	}
	NewWindow->SetNEdges(NEdgesCopy);

	/////////// WARNING!!!! ERROR HERE!!!
	// Copy the simplices that should still be part of the window (i.e., simplices with at least one edge in the window)
//	list <Simplex *>::iterator SIt;
	//
	//	int Count = 0;

	//	for (SIt = this->Simplices->begin(); SIt != this->Simplices->end(); SIt++)
	//	{
	//		Simplex *Simplex0 = *SIt;
	//		Simplex *NewSimplex = new Simplex();

	//		NewSimplex->SetDim(2);

	//		Point 	*P0 = Simplex0->ReadPoint(0),
	//				*P1 = Simplex0->ReadPoint(1),
	//				*P2 = Simplex0->ReadPoint(2);

	//	double 	Coords0[DIM_SMA] = {P0->GetCoordinate(X_COORD), P0->GetCoordinate(Y_COORD), P0->GetCoordinate(Z_COORD)},
	//			Coords1[DIM_SMA] = {P1->GetCoordinate(X_COORD), P1->GetCoordinate(Y_COORD), P1->GetCoordinate(Z_COORD)},
	//			Coords2[DIM_SMA] = {P2->GetCoordinate(X_COORD), P2->GetCoordinate(Y_COORD), P2->GetCoordinate(Z_COORD)};
				//
	//	if (	this->GetEdgeByCoordinates(Coords0, Coords1) != NULL ||
	//			this->GetEdgeByCoordinates(Coords0, Coords2) != NULL ||
	//			this->GetEdgeByCoordinates(Coords1, Coords2) != NULL)
	//	{
	//		NewSimplex->InsertPoint(this->GetPointByCoordinates(Coords0));
	//		NewSimplex->InsertPoint(this->GetPointByCoordinates(Coords1));
	//		NewSimplex->InsertPoint(this->GetPointByCoordinates(Coords2));
			//
	//		NewWindow->InsertSimplex(NewSimplex);
	//	}

//	cout << Count << endl;
	//	Count++;


	//	}

	// Return the copied complex
	return NewWindow;

}

Window* Window::PruneWindow()
{

	Window 	*NewWindow,
			*SwapWindow;
	int Count = 0;

	if (this->Pruned())
	{
		return this;
	}
	else
	{
		NewWindow = this->PruneWindowOnce();
	}


	while (!NewWindow->Pruned())
	{

		//cout << "NewWindow simplices: " << NewWindow->GetSimplices()->size() << endl;

		SwapWindow = NewWindow->PruneWindowOnce();
		delete NewWindow;
		NewWindow = SwapWindow;

		if (Count > this->GetNPoints())
		{
			cout << "Error pruning Window: maximum number of possible prune steps surpassed" << endl;
			delete NewWindow;
			return NULL;
		}
		Count++;


	}

	return NewWindow;

}*/
