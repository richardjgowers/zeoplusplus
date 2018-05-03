/*
 * SMA.cc
 *
 *  Created on: Mar 4, 2016
 *      Author: ismael.gomez
 */

#include "SMA.hh"

/////	GENERAL MOLECULE INFO
///		Object that will store relevant information to avoid unnecesary recomputations under user's queries.
MoleculeInfo *Molecule = NULL;

char CurrentFileName[MAX_FILENAME] = "\0";

void LoadMolecule(char *Filename, int ExecMode)
{
	//cout << "Loading molecule from " << Filename << endl;

	// If the filename has changed, recompute
	if (Molecule != NULL)
	{
		if ( strcmp(CurrentFileName, Filename) )
		{
			strcpy(Filename, CurrentFileName);
			delete Molecule;
			//Molecule = Preprocess(Filename);
			Molecule = MolecularCageAnalysis(Filename, ExecMode);
		}
	}
	// If molecule info is empty, fill it
	else	{
		//Molecule = Preprocess(Filename);
		Molecule = MolecularCageAnalysis(Filename, ExecMode);
	}

}


void ReleaseMolecule()
{
	if (Molecule == NULL)
	{
		//cout << "Warning: trying to release non existing molecule (ReleaseMolecule() at SMA.cc)" << endl;
		return;
	}

	delete Molecule;

}

/*void LoadMolecule(char *Filename, MoleculeInfo *NewMolecule)
{
	cout << "Loading molecule from " << Filename << endl;

	// If the filename has changed, recompute
	NewMolecule = Preprocess(Filename);

}*/

/*
 *
 */
void IntToString (int Number, char *String)
{

	int LastDigit = Number%10;

	int FirstDigit = (Number-LastDigit)/10;

	String[0] = (char) (FirstDigit+48);
	String[1] = (char) (LastDigit+48);
	String[2] = '\0';

}


char *MoleculeName()
{

	/*if (Molecule == NULL)
	{
		return NULL;
	}*/

	return Molecule->GetChemical()->GetCName();


}

char *MoleculeName(char *Name)
{
	if (Molecule == NULL)
	{
		return NULL;
	}
	strncpy(Name, Molecule->GetChemical()->GetCName(), FILENAME_MAX);

	return Molecule->GetChemical()->GetCName();
}


/******************************/
/****** INTERNAL VOLUMES ******/
/******************************/

double MaximumInternalSphere()
{

	if (Molecule == NULL)
	{
		//cout << "Error calculating maximum internal sphere size: molecule has not been correctly loaded" << endl;
		return -1.;
	}
	return LargestInternalSphere(Molecule->GetVoronoiGraph());

}

void MaximumInternalSphereCoordinates(double *Coordinates)
{
	if (Molecule == NULL)
	{
		//cout << "Error calculating maximum internal node coordinates: molecule has not been correctly loaded" << endl;
		return;
	}

	return LargestInternalNodeCoordinates(Molecule->GetVoronoiGraph(), Coordinates);

}

double TotalInternalVolume ()
{

	//PreCompute(Filename);
	if (Molecule == NULL)
	{
		//cout << "Error calculating total internal volume: molecule has not been correctly loaded" << endl;
		return -1.;
	}

	return VoroTotalVolume(Molecule->GetVoronoiGraph());

}

double TotalConvexHullVolume ()
{

	//PreCompute(Filename);
	if (Molecule == NULL)
	{
		//cout << "Error calculating total internal volume: molecule has not been correctly loaded" << endl;
		return -1.;
	}

	return VoroTotalVolume(Molecule->GetVoroConvexHull());

}

void InternalCells(vector <double *> &Cells)
{

	// Read the cell coordinates and radius from Voronoi graph


	for (int i = 0; i < Molecule->GetVoronoiGraph()->GetNPoints(); i++)
	{
		Point *PAux = Molecule->GetVoronoiGraph()->GetPointList()[i];
		double *CellInfo = new double[DIM_SMA+1];

		if (PAux->GetType() == internal_point)
		{
			CellInfo[X_COORD] = PAux->GetCoordinate(X_COORD);
			CellInfo[Y_COORD] = PAux->GetCoordinate(Y_COORD);
			CellInfo[Z_COORD] = PAux->GetCoordinate(Z_COORD);

			CellInfo[DIM_SMA] = PAux->GetRadius();

			Cells.push_back(CellInfo);
		}

	}

}

/******************************/
/********** WINDOWS ***********/
/******************************/

bool NumberOfWindows(int &NWindows)
{

	//PreCompute(Filename);

	if (Molecule == NULL)
	{
		//cout << "Error calculating number of windows: molecule has not been correctly loaded" << endl;
		return -1;
	}

	NWindows = Molecule->GetNWindows();

	return Molecule->GetWinCorrect();
}

/*
 * 		Return the number of windows (in parameter) that have a number of simplices equal or bigger than WinNTriangles
 */
bool NumberOfWindows(int &NWindows, int WinNTriangles)
{
	if (Molecule == NULL)
	{
		//cout << "Error calculating number of windows: molecule has not been correctly loaded" << endl;
		return -1;
	}

	Window **MoleculeWindows = Molecule->GetWindows();
	int NTotWin = Molecule->GetNWindows();

	NWindows = NTotWin;

/*	NWindows = 0;

	for (int i = 0; i < NTotWin; i++)
	{
		Window *Win = MoleculeWindows[i];

		if (Win->GetNPoints() >= WinNTriangles)
		{
			NWindows++;
			//cout << "Window size: " << Win->GetSurface() << " - Triangles: " << Win->GetSimplices()->size() << endl;
		}

	}*/



	return Molecule->GetWinCorrect();

}


int NumberOfEntryPaths()
{
	if (Molecule == NULL)
	{
		//cout << "Error calculating number of entry paths: molecule has not been correctly loaded" << endl;
		return -1;
	}

	return VoroEntryPaths(Molecule->GetVoronoiGraph());

}

void MoleculeNthWindow (	int Index,
							std::vector<double *> &WindowPoints,
							double &WindowSurface)
{

	if (Molecule == NULL)
	{
		//cout << "Error calculating Nth window: molecule has not been correctly loaded" << endl;
		return;
	}

	Window *WindowT = NULL;

	if (Index < Molecule->GetNWindows() )
		WindowT = Molecule->GetWindows()[Index];
	else
		return;

	for (int j = 0; j < WindowT->GetNPoints(); j++)
	{
		Point *WindowPoint = WindowT->GetPointList()[j];

		// Allocate three positions for each vector
		double *WindowInfo = new double[DIM_SMA];

		// Store the information about the coordinates of the points
		for (int k = 0; k < DIM_SMA; k++)
		{
			WindowInfo[k] = WindowPoint->GetCoordinate(k);
		}

		WindowPoints.push_back(WindowInfo);
	}

	WindowSurface = WindowT->GetSurface();


}

double TotalSurface ()
{
	double TotSurf = 0.0;

	for (int i = 0; i < Molecule->GetNWindows(); i++)
	{
		TotSurf += Molecule->GetWindows()[i]->GetSurface();
	}

	return TotSurf;

}

double MaximumAccessibleSize()
{
	if (Molecule == NULL)
	{
		//cout << "Error calculating accessible size: molecule has not been correctly loaded" << endl;
		return -1.;
	}

	return MoleculeAccessibleSize(Molecule->GetVoronoiGraph());
}

double TargetAccessibleSize(int CellId)
{
	Point *PAux = Molecule->GetVoronoiGraph()->GetPointById(CellId);

	return CellAccessibleSize(Molecule->GetVoronoiGraph(), PAux);
}

/*
 * 	Given the molecule, returns the shadow characteristic of its center
 */

double ShadowCharacteristic()
{
	return Molecule->GetShadowCharact();
}

/*
 * 	Return surface area
 */
double InternalSurfaceArea()
{
	return SurfaceArea_MonteCarlo(Molecule);
}



/**********/

void PrintWindowsToFile(char *Basename)
{

	vector <Complex *> AS_boundary;
	char AlphaBoundName[MAX_FILENAME];
	strncpy(AlphaBoundName, Basename, MAX_FILENAME);
	strncat(AlphaBoundName, "/Alpha_bound.mol", MAX_FILENAME);

	for (int i = 0; i < Molecule->GetNWindows(); i++)
	{
		char WindowName[MAX_FILENAME];
		strncpy(WindowName, Basename, MAX_FILENAME);
		strncat(WindowName, "/Window_", MAX_FILENAME);


		char a[3] = {'\0', '\0', '\0'};
		IntToString(i, a);

		strncat(WindowName, a, MAX_FILENAME);
		//strncat(WindowName, ".mol", MAX_FILENAME);
		//strncat(WindowName, ".vtk", MAX_FILENAME);

		Window *TWind = Molecule->GetWindows()[i];

		Complex *SC = TWind->JointSimplices();

		//PrintToMolFile(TWind->JointSimplices(), WindowName);
		//PrintToVTK(TWind->JointSimplices(), WindowName);
		PrintToVTK(TWind, WindowName);

		AS_boundary.push_back( SC );
	}
	PrintToMolFile(AS_boundary, AlphaBoundName);
	//PrintToVTK(AS_boundary, AlphaBoundName);

	/////	PRINT WINDOWS
	///

	for (int i = 0; i < Molecule->GetNWindows(); i++)
	{
		///		Build window name
		char WindowName[MAX_FILENAME];
		strncpy(WindowName, Basename, MAX_FILENAME);
		strncat(WindowName, "/WindowT_", MAX_FILENAME);

		char a[3] = {'\0', '\0', '\0'};
		IntToString(i, a);

		strncat(WindowName, a, MAX_FILENAME);

		///		Print window to VTK

		Window *TWind = Molecule->GetWindows()[i];


		PrintWindowToVTK(TWind, WindowName);

	}



	/////	PRINT MOLECULE
	///

	//cout << "Printing molecule description" << endl;

	char MoleculeVTKName[MAX_FILENAME];

	strncpy(MoleculeVTKName, Basename, MAX_FILENAME);
	MoleculeVTKName[strlen(MoleculeVTKName)] = '/';
	strncat(MoleculeVTKName, "MoleculeDesc", MAX_FILENAME);

	Molecule->GetChemical()->Normalize();

	PrintToVTK(Molecule->GetChemical(), MoleculeVTKName);

	//strncat(MoleculeVTKName, ".mol", MAX_FILENAME);

	//PrintToMolFile(Molecule->GetChemical(), MoleculeVTKName);


}


void PrintVoronoiGraphMolecule(char *FileBaseName)
{
	if (Molecule == NULL)
	{
		//cout << "Error printing Voronoi graph: molecule has not been correctly loaded" << endl;
		return;
	}

	char Filename[MAX_FILENAME];

	///		Print interior of voronoi network

	strncpy(Filename, FileBaseName, MAX_FILENAME);

	strncat(Filename, "/Voronoi_interior", MAX_FILENAME);

	Complex *InteriorVoro = Molecule->GetVoronoiGraph()->ComplexCopy(internal_point);

	InteriorVoro->Normalize();

	PrintToVTK(InteriorVoro, Filename);

	delete InteriorVoro;

	///		Print boundary network

	strncpy(Filename, FileBaseName, MAX_FILENAME);

	strncat(Filename, "/Voronoi_boundary", MAX_FILENAME);

	Complex *BoundVoro = Molecule->GetVoronoiGraph()->ComplexCopy(boundary_point);

	BoundVoro->Normalize();

	PrintToVTK(BoundVoro, Filename);

	delete BoundVoro;

	///		Print voronoi network complete

	strncpy(Filename, FileBaseName, MAX_FILENAME);

	strncat(Filename, "/Voronoi", MAX_FILENAME);

	PrintToVTK(Molecule->GetVoronoiGraph(), Filename);


	///		Print voronoi convex hull

	strncat(Filename, "_CH.mol", MAX_FILENAME);

	PrintToMolFile(Molecule->GetVoroConvexHull(), Filename);


	/// TODO: Set this up
/*	Complex *InteriorVoroCH = Molecule->GetVoroConvexHull()->ComplexCopy(internal_point);

	strncat(Filename, "_interior", MAX_FILENAME);

	PrintToVTK(InteriorVoroCH, Filename);*/

}


void PrintChemical (char *FileBaseName)
{
	if (Molecule == NULL)
	{
		//cout << "Error printing Voronoi graph: molecule has not been correctly loaded" << endl;
		return;
	}

	char Filename[MAX_FILENAME];

	strncpy(Filename, FileBaseName, MAX_FILENAME);

	///		Print to VTK

	strncat(Filename, "/MoleculeDescription", MAX_FILENAME);

	PrintToVTK(Molecule->GetChemical(), Filename);

	///		Print to XYZ

	char xyzFilename[MAX_FILENAME];

	strncpy(xyzFilename, FileBaseName, MAX_FILENAME);

	//strncat(xyzFilename, "/MoleculeDescription_xyz.xyz", MAX_FILENAME);

	//cout << xyzFilename << endl;

	PrintToXYZFile ( Molecule->GetChemical(), xyzFilename);


	///		Print to mol

	strncat(Filename, ".mol", MAX_FILENAME);

	PrintToMolFile(Molecule->GetChemical(), Filename);

}

void PrintChemicalSDF (char *Location)
{
	if (Molecule == NULL)
	{
		//cout << "Error printing Voronoi graph: molecule has not been correctly loaded" << endl;
		return;
	}

	if (Molecule->GetChemical()->GetCName() == NULL)
	{
		return;
	}

	char Filename[MAX_FILENAME];

	strncpy(Filename, Location, MAX_FILENAME);

	strncat(Filename, Molecule->GetChemical()->GetCName(), MAX_FILENAME);

	PrintToVTK(Molecule->GetChemical(), Filename);


}

void PrintAlphaShape (char *Basename)
{

	if (Molecule == NULL)
	{
		//cout << "Error printing Voronoi graph: molecule has not been correctly loaded" << endl;
		return;
	}

	char Filename[MAX_FILENAME];

	strncpy(Filename, Basename, MAX_FILENAME);

	///		Print to VTK

	strncat(Filename, "/AlphaShape", MAX_FILENAME);

	Complex *AShapeC = AlphaShapeToComplex(Molecule);

	PrintToVTK(AShapeC, Filename);

	delete AShapeC;


}

void PrintAlphaShape_Extended (char *Basename)
{

	if (Molecule == NULL)
	{
		//cout << "Error printing Voronoi graph: molecule has not been correctly loaded" << endl;
		return;
	}

	char Filename[MAX_FILENAME];

	strncpy(Filename, Basename, MAX_FILENAME);

	///		Print to VTK

	strncat(Filename, "/AlphaShape_Extended", MAX_FILENAME);

	list <Point *> InternalPoints;

	for (int i = 0; i < Molecule->GetVoronoiGraph()->GetNPoints(); i++)
	{
		Point *VP = Molecule->GetVoronoiGraph()->GetPointList()[i];

		if (VP->GetType() == internal_point)
		{
			InternalPoints.push_back(VP);
		}
	}

	Complex *AShapeC = AlphaShapeToComplex_ExtraPoints(Molecule, InternalPoints);

	PrintToVTK(AShapeC, Filename);

	delete AShapeC;


}

/*****************************************/
/******* INTERNAL POINTS COMPUTING *******/
/*****************************************/

/*
 * 	Given an XYZ file input with points' coordinates, produces an output with every point classified as internal/external.
 */
void ClassifyXYZInternal (	char *Filename)
{

	list <Point *> *Points = new list <Point *>();

	ReadXYZPoints(Points, Filename);

	MoleculeInternalPoints (Molecule->GetChemical(), Points);

	///	TODO - PRINT OUTPUT

	delete Points;


}


/****************************/
/******* SDF HANDLING *******/
/****************************/



//vector <Complex *> *Chemicals = new vector <Complex *>();
vector <Complex *> *Chemicals = NULL;

int LoadSDFMolecule (char *Filename, bool SDF_Loosely)
{

	Chemicals = new vector <Complex *>();

	///		Read molecules from SDF file
	if (!SDF_Loosely)
	{
		ReadSDF(Chemicals, Filename);
	}
	else
	{
		ReadSDFL(Chemicals, Filename);
	}

	return Chemicals->size();

}


/*	Set ith molecule from list as current chemical and compute all the parameters related to it. Release previous molecule information if present.
 *
 */
bool SetCurrentMolecule (int Position, int ExecMode)
{

	Complex *Chemical = Chemicals->at(Position);

	if (Molecule != NULL)
	{
		delete Molecule;
	}

	Molecule = PerformMCA (Chemical, ExecMode);

	if (Molecule == NULL)
	{
		return false;
	}

	return true;

}

/*
 *
 */
void ReleaseSDFChemicals()
{


	for (int i = 0; i < Chemicals->size(); i++)
	{

		//SetCurrentMolecule(i);

		Complex *LocalChem = Chemicals->at(i);

		delete LocalChem;

	}

	delete Chemicals;

	delete Molecule;

}


/*************************************************/
/******* POINTS INSIDE MOLECULES STUDIES   *******/
/*************************************************/

bool MoleculesOverlap()
{

	for (int i = 0; i < Chemicals->size(); i++)
	{
		Complex *Chemical1 = Chemicals->at(i);

		for (int j = i+1; j < Chemicals->size(); j++)
		{
			Complex *Chemical2 = Chemicals->at(j);

			if (IEClas_MoleculesOverlap(Chemical1, Chemical2))
			{
				return true;
			}
		}

	}

	return false;

}

/*	PointsInMolecules
 *
 * 		Given a list of points, checks if its contained in any of the molecules loaded in Chemicals, printing it in case it is.
 *
 */
void PointsInMolecules(list <Point *> &PointList)
{

	int NInternalPoints = 0;

	int j = 0;

	///		Vector to store the id of the molecule where each point is stored in (if none, -1 is kept)
	vector <int> WhichMolecule(PointList.size(), -1);

	///		Run over all the molecules read
	for (int i = 0; i < Chemicals->size(); i++)
	{
		Complex *Chemical = Chemicals->at(i);

		int MolInternalPoints = 0;

		IEClas_PointsInMolecule(Chemical, PointList);

		list <Point *>::iterator PIt;

		///		Run over points after processing for local molecule
		for (PIt = PointList.begin(); PIt != PointList.end(); PIt++)
		{
			Point *P = *PIt;

			///		Print and erase points classified as internal
			if (P->GetType() == internal_point)
			{
				/*if (Chemical->GetCName() != NULL)
				{
					cout << j++ << " " << P->GetCoordinate(X_COORD) << " " << P->GetCoordinate(Y_COORD) << " " << P->GetCoordinate(Z_COORD) << " " << i << " " << Chemical->GetCName() << endl;
				}
				else
				{
					cout << P->GetCoordinate(X_COORD) << " " << P->GetCoordinate(Y_COORD) << " " << P->GetCoordinate(Z_COORD) << " " << i << endl;
				}

				PIt = PointList.erase(PIt);
				PIt--;*/

				WhichMolecule[j] = i;
				NInternalPoints++;
				MolInternalPoints++;

			}
			j++;
		}

		cout << "Points in " << Chemical->GetCName() << ": " << MolInternalPoints << endl;

	}

	cout << "Total internal points (inside any molecule): " << NInternalPoints << endl;

	///		Print points that are not inside any of the molecules (return if none of them are)
	/*if (PointList.size() == 0)
	{
		return;
	}

	list <Point *>::iterator PIt2;

	for (PIt2 = PointList.begin(); PIt2 != PointList.end(); PIt2++)
	{
		Point *P = *PIt2;

		cout << P->GetCoordinate(X_COORD) << " " << P->GetCoordinate(Y_COORD) << " " << P->GetCoordinate(Z_COORD) << " -1" << endl;
	}*/

	///		Print the points along with the molecule they belong to

	list <Point *>::iterator PIt;

	j = 0;

	for (PIt = PointList.begin(); PIt != PointList.end(); PIt++)
	{
		Point *P = *PIt;

		//	Get chemical's index
		int Index = WhichMolecule[j];

		//	Get the chemical (if the index is positive)
		if (Index >= 0)
		{
			Complex *Chemical = Chemicals->at(Index);

			//	Print the point along with the chemical it belongs to
			if (Chemical->GetCName() != NULL)
			{
				cout << j << " " << P->GetCoordinate(X_COORD) << " " << P->GetCoordinate(Y_COORD) << " " << P->GetCoordinate(Z_COORD) << " " << Chemical->GetCName() << endl;
			}
			else
			{
				cout << P->GetCoordinate(X_COORD) << " " << P->GetCoordinate(Y_COORD) << " " << P->GetCoordinate(Z_COORD) << " " << j << endl;
			}
		}
		//	Print the point along with the -1 code (point is not inside any of the molecules given)
		else
		{
			cout << j << " " << P->GetCoordinate(X_COORD) << " " << P->GetCoordinate(Y_COORD) << " " << P->GetCoordinate(Z_COORD) << " " << -1 << endl;
		}

		//	Erase the point after processing
		PIt = PointList.erase(PIt);
		PIt--;

		j++;

	}

}

