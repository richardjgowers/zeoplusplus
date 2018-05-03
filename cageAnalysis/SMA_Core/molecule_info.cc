//	Software for single molecule analysis
//
//	complex_set.cc
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : February 5th 2016

#include "molecule_info.hh"

MoleculeInfo::MoleculeInfo(Complex *Chemical,
							Window **Windows,
							Complex *VoronoiGraph,
							Complex *VoroConvexHull,
							int NWindows,
							bool WinCorrect,
							double ShadowChar)
{

	this->Chemical = Chemical;
	this->Windows = Windows;
	this->VoronoiGraph = VoronoiGraph;

	this->VoroConvexHull = VoroConvexHull;

	this->NWindows = NWindows;

	this->WinCorrect = WinCorrect;

	this->ShadowCharact = ShadowCharact;

}

MoleculeInfo::MoleculeInfo()
{

	this->Chemical = NULL;
	this->Windows = NULL;
	this->VoronoiGraph = NULL;

	this->VoroConvexHull = NULL;

	this->NWindows = 0;

	this->WinCorrect = true;

	this->ShadowCharact = -1.;

}


MoleculeInfo::~MoleculeInfo()
{

	//cout << "MoleculeInfo destructor called" << endl;

	///		DELETE CHEMICAL
	if (this->Chemical != NULL)
	{
		delete this->Chemical;
	}


	///		DELETE WINDOWS
	for (int i = 0; i < this->NWindows; i++)
	{
		Window *LW = this->Windows[i];
		delete LW;
	}

	///		DELETE VORONOI GRAPHS

	delete this->VoronoiGraph;

	delete this->VoroConvexHull;

}

/*******************/
/***** GETTERS *****/
/*******************/

Complex* MoleculeInfo::GetChemical()
{
	return this->Chemical;
}

Window ** MoleculeInfo::GetWindows()
{
	return this->Windows;
}

Complex * MoleculeInfo::GetVoronoiGraph()
{
	return this->VoronoiGraph;
}

Complex * MoleculeInfo::GetVoroConvexHull()
{
	return this->VoroConvexHull;
}

int MoleculeInfo::GetNWindows()
{
	return this->NWindows;
}

bool MoleculeInfo::GetWinCorrect()
{
	return this->WinCorrect;
}

double MoleculeInfo::GetShadowCharact()
{
	return this->ShadowCharact;
}

/*
 * Get the biggest node (Radius based) with the given type from the Voro graph.
 */
Point * MoleculeInfo::GetGreatestVoroNode(PointType Type)
{
	if (this->VoronoiGraph == NULL)
	{
		return NULL;
	}

	if (this->VoronoiGraph->GetNPoints() <= 0)
	{
		return NULL;
	}

	Point *MaxPoint = NULL; // = this->VoronoiGraph->GetPointList()[0];
	double MaxRadius = 0.0;

	for (int i = 0; i < this->VoronoiGraph->GetNPoints(); i++)
	{
		Point *VoroPoint = this->VoronoiGraph->GetPointList()[i];

		if (VoroPoint->GetType() == Type)
		{
			if (VoroPoint->GetRadius() > MaxRadius)
			{
				MaxPoint = VoroPoint;
				MaxRadius = VoroPoint->GetRadius();
			}
		}

	}

	return MaxPoint;
}



/*******************/
/***** SETTERS *****/
/*******************/

void MoleculeInfo::SetChemical(Complex *Chemical)
{
	this->Chemical = Chemical;
}

void MoleculeInfo::SetWindows(Window **Windows)
{
	this->Windows = Windows;
}

void MoleculeInfo::SetVoronoiGraph(Complex *VoronoiGraph)
{
	this->VoronoiGraph = VoronoiGraph;
}

void MoleculeInfo::SetVoroConvexHull(Complex *VoroConvexHull)
{
	this->VoroConvexHull = VoroConvexHull;
}

void MoleculeInfo::SetNWindows(int NWindows)
{
	this->NWindows = NWindows;
}

void MoleculeInfo::SetWinCorrect(bool WinCorrect)
{
	this->WinCorrect = WinCorrect;
}

void MoleculeInfo::SetShadowCharact(double ShadowCharact)
{
	this->ShadowCharact = ShadowCharact;
}


/*************************/
/***** OTHER METHODS *****/
/*************************/

/*	InsertWindow
 * 		Insert the given window into the molecule info structure if it was not already there (same set of bonds)
 */
void MoleculeInfo::InsertWindow(Window *WindowT)
{

	if (this->WindowIn(WindowT) )
	{
		//cout << "Warning: Window not inserted (already in MoleculeInfo object)" << endl;
		return;
	}

	if (this->Windows == NULL) {
		this->Windows = new Window*[MAX_WINDOWS];
	}

	this->Windows[NWindows] = WindowT;
	this->NWindows++;

}

/*	InsertWindow
 * 		Insert the given window into the molecule info structure if it was not already there (or if there were other windows that together form the same window).
 */
void MoleculeInfo::InsertWindow(	Window *WindowT,
									bool NonCyclic)
{

	///		If the flag is false, use the simpler method.

	//cout << "HOLA" << endl;

	if (!NonCyclic)
	{
		this->InsertWindow(WindowT);
	}

	/////	CHECK: WINDOW = UNION OF WINDOWS ALREADY IN MOL_INFO
	///		For every edge in the window, check if there is a window that contains it.
	//		If there is at least one edge that wasn't present in another window, the window gets inserted.

	bool IsComposite = true;

	///		Go through all the edges of the window
	for (int i = 0; i < WindowT->GetNEdges(); i++)
	{
		Edge *WEdge = WindowT->GetEdgeList()[i];

		///		Check every other window in the structure
		for (int j = 0; j < this->GetNWindows(); j++)
		{
			Window *Win0 = this->GetWindows()[j];

			Point 	*WPSource = WEdge->GetPSource(),
					*WPTarget = WEdge->GetPTarget();

			double 	PSCoords[DIM_SMA] = {WPSource->GetCoordinate(X_COORD), WPSource->GetCoordinate(Y_COORD), WPSource->GetCoordinate(Z_COORD)},
					PTCoords[DIM_SMA] = {WPTarget->GetCoordinate(X_COORD), WPTarget->GetCoordinate(Y_COORD), WPTarget->GetCoordinate(Z_COORD)};

			if (Win0->GetEdgeByCoordinates(PSCoords, PTCoords) == NULL)
			{
				IsComposite = false;
				break;
			}
		}

		///		If some edge was found that is not present in other window, insert the given window.
		if (!IsComposite)
		{
			break;
		}

	}

	///		If the window is built by other windows, don't insert.
	if (IsComposite)
	{
		return;
	}

	this->InsertWindow(WindowT);



}

void MoleculeInfo::DeleteWindows()
{

	delete[] this->Windows;

	this->Windows = NULL;

	this->NWindows = 0;

}

void MoleculeInfo::PrintWindows()
{

	for (int i = 0; i < this->NWindows; i++)
	{
		cout << "Printing window: " << i << endl;
		this->Windows[i]->PrintWindow();
		cout << endl;
	}

	cout << "The molecule has " << this->NWindows << " entry windows" << endl;

}


/*
 * True if given window is already inside the molecule info object. Checks only points (edges are assumed to be the same).
 */
bool MoleculeInfo::WindowIn(Window *WindowT)
{

	for (int i = 0; i < this->NWindows; i++)
	{

		Window *LocalWindow = this->Windows[i];

		if (WindowT->SameWindow(LocalWindow))
		{
			return true;
		}

	}

	return false;



}


// MOLECULE_INFO_CC

