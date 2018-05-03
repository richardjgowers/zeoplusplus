/*
 * vtk_out.cc
 *
 *  Created on: Mar 6, 2017
 *      Author: ismael.gomez
 */

#include "vtk_out.hh"

void PrintWindowsToFile(char *Basename, MoleculeInfo *Molecule)
{

	vector <Complex *> AS_boundary;
	char AlphaBoundName[MAX_FILENAME];
	strncpy(AlphaBoundName, Basename, MAX_FILENAME);
	strncat(AlphaBoundName, "/Alpha_bound.mol", MAX_FILENAME);
	//strncat(AlphaBoundName, "/Alpha_bound.vtk", MAX_FILENAME);

	cout << "Printing alpha bound" << endl;

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
		PrintToVTK(TWind->JointSimplices(), WindowName);

		AS_boundary.push_back( SC );
	}
	PrintToMolFile(AS_boundary, AlphaBoundName);
	//PrintToVTK(AS_boundary, AlphaBoundName);

	/////	PRINT WINDOWS
	///

	cout << "Printing windows (individually)" << endl;

	//cout << "WARNING: Print window to VTK - triangles - is disabled due to an error" << endl;

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

	cout << "Printing molecule description" << endl;

	char MoleculeVTKName[MAX_FILENAME];

	strncpy(MoleculeVTKName, Basename, MAX_FILENAME);
	MoleculeVTKName[strlen(MoleculeVTKName)] = '/';
	strncat(MoleculeVTKName, "MoleculeDesc", MAX_FILENAME);

	Molecule->GetChemical()->Normalize();

	PrintToVTK(Molecule->GetChemical(), MoleculeVTKName);

	//strncat(MoleculeVTKName, ".mol", MAX_FILENAME);

	//PrintToMolFile(Molecule->GetChemical(), MoleculeVTKName);


}


void PrintVoronoiGraphMolecule(char *FileBaseName, MoleculeInfo *Molecule)
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


/*
 *
 */
void PrintChemical (char *FileBaseName, MoleculeInfo *Molecule)
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


	///		Print to mol

	strncat(Filename, ".mol", MAX_FILENAME);

	PrintToMolFile(Molecule->GetChemical(), Filename);

}


/*
 *
 */
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


/*
 *
 */
void PrintAlphaShape (char *Basename, MoleculeInfo *Molecule)
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
