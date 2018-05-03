// Voro++, a 3D cell-based Voronoi library
//
// single_mol_analysis.cc
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>

/*#include "SMA_Analysis/dijkstra_path.hh"
#include "SMA_Analysis/windows.hh"
#include "SMA_Analysis/internal_cells.hh"

#include "SMA_GeoGraphComp/voro_interface.hh"


#include "SMA_Interface/preprocess.hh"*/


//#include "SMA_Preprocess/voroedge_classificator.hh"

#include "SMA.hh"

//#include "SMA_Interface/complex_io.hh"

//#include "SMA_GeoGraphComp/clustering.hh"

//#include "SMA_GeoGraphComp/alpha_shape.hh"

//#include "SMA_Core/geometry.hh"

//#include "face_compute.hh"

#include <iostream>

using namespace std;

#define INPUT_EXT_NONE		0
#define INPUT_EXT_MOL		1
#define INPUT_EXT_FRAMEID	2
#define	INPUT_EXT_SDF		3


#define OUTPUT_LIMIT		1024	//	Maximum size of text output



/*
 *
 */
void ExtractBasename(char *Filename, char *Basename, char *Extension)
{

	char Separator = '.';

	int Len = strlen(Filename);

	bool ExtFlag = false;

	int j = 0;

	for (int i = 0; i < Len; i++)
	{
		if (Filename[i] == Separator)
		{
			Basename[i] = '\0';
			ExtFlag = true;
		}
		else if (!ExtFlag)
		{
			Basename[i] = Filename[i];
		}
		else
		{
			Extension[j] = Filename[i];
			j++;
		}

	}

	Extension[j] = '\0';

}

/*
 *
 */


int ExecutionOptions = 1;

int ProcessOption (char *OptString)
{

	if (strncmp(OptString, "-w", MAX_FILENAME) == 0)
	{
		ExecutionOptions *= WINDOW_COMP;
	}

	if (strncmp(OptString, "-p", MAX_FILENAME) == 0)
	{
		ExecutionOptions *= PRINT_OUTPUT;
	}

	if (strncmp(OptString, "-e", MAX_FILENAME) == 0)
	{
		ExecutionOptions *= EXPLORATORY;
	}

	if (strncmp(OptString, "-np", MAX_FILENAME) == 0)
	{
		ExecutionOptions *= NO_PRUNE;
	}

	if (strncmp(OptString, "-olt", MAX_FILENAME) == 0)
	{
		ExecutionOptions *= MOLS_OVERLAP;
	}

	if (strncmp(OptString, "-iep", MAX_FILENAME) == 0)
	{
		ExecutionOptions *= IE_ANALYSIS;
		return IE_ANALYSIS;
	}

	if (strncmp(OptString, "-sdfl", MAX_FILENAME) == 0)
	{
		ExecutionOptions *= SDFL;
		return SDFL;
	}

	return 0;
}


/*
 *		Prints output line with all parameters required by the user. Also creates a folder with molecule name and fills it with visual information if that option was chosen.
 */
void PrintOutput (	char *Basename,
					int ExecMode)
{


	///	OUTPUT LEGEND
	cout << "Maxium Internal Sphere (MIS); Total Internal Volume (TIV); Convex Hull Volume (CHV); Pore Exposure Ratio (PER); Internal Surface Area (ISA)";

	if (ExecMode%WINDOW_COMP == 0)
	{
		cout << " ; Number of Windows (NW); Number of Entry Paths (NEP); Internal Accessible Size (IAS)";
	}

	cout << endl;

	///		PRINT MOLECULE'S NAME

	//cout << Basename;

	//cout << Molecule->GetChemical()->GetCName();

	char *Name = MoleculeName();

	//cout << Name << endl;

	cout << "& ";

	if (Name != NULL) {
		cout << Name;
	}
	else {
		cout << Basename;
	}

	/// MAXIMUM INTERNAL SPHERE AND TOTAL INTERNAL VOLUME

	double MaxIntSphere = MaximumInternalSphere();

	cout << " MIS " << MaxIntSphere;

	double TotIntVol = TotalInternalVolume();

	cout << " TIV " << TotIntVol;

	double TotCHVol = TotalConvexHullVolume();

	cout << " CHV " << TotCHVol;

	double ShadChar = ShadowCharacteristic();

	cout << " PER " << ShadChar;

	double IntSurfArea = InternalSurfaceArea();

	cout << " ISA " << IntSurfArea;

	if (ExecMode%WINDOW_COMP != 0)
	{
		cout << endl;
	}

	///		Print molecule structure and voronoi graph (if selected by the user)
	char Location[OUTPUT_LIMIT];

	if (ExecMode%PRINT_OUTPUT == 0)
	{
		char SystOut[OUTPUT_LIMIT];

		strncpy(SystOut, "mkdir ", OUTPUT_LIMIT);

		strncpy(Location, "Molecules_Output/", OUTPUT_LIMIT);
		if (Name != NULL)
		{
			strncat(Location, Name, OUTPUT_LIMIT);
		}
		else
		{
			strncat(Location, "Molecule", OUTPUT_LIMIT);
		}

		strncat(SystOut, Location, OUTPUT_LIMIT);
		system(SystOut);

		PrintVoronoiGraphMolecule(Location);

		PrintChemical(Location);

		PrintAlphaShape_Extended(Location);
	}

	///		Window information
	if (ExecMode%WINDOW_COMP == 0)
	{

		///	TOTAL PATHS AND WINDOWS

		int NWindowsTot = 0;

		bool WC = NumberOfWindows(NWindowsTot, 0);		///		Print only windows with 5 or more simplices
		cout << " NW " << NWindowsTot;

		if (!WC) {
			cout << "*";
		}


		int NPathsTot = NumberOfEntryPaths();

		cout << " NEP " << NPathsTot;

		///		Total surface 		///

		double TotSurf = TotalSurface();

		cout << " TCS " << TotSurf ;

		/// 	Maximum accessible size		///

		double IntAccessSize = MaximumAccessibleSize();
		cout << " IAS " << IntAccessSize;

		cout << endl;

		///	CREATE FOLDER FOR INFO STORAGE

		///		Print windows and alpha shape (if selected by the user)
		if (ExecMode%PRINT_OUTPUT == 0)
		{
			PrintAlphaShape(Location);

			///		Store windows		///
			PrintWindowsToFile(Location);

		}


	}



	/*if (TotIntVol > 0 || MaxIntSphere > 0)
	{
		char Location[MAX_FILENAME];

		strncpy(Location, "Molecules/", MAX_FILENAME);

		PrintChemicalSDF(Location);
	}*/

	//cout << endl;



}

int main( int argc, char *argv[])
{

	/////	PARSE INPUT
	///		Check the number of arguments.
	//		Store the filename (right now only .xyz accepted).
	//		TODO
	//			Accept .xyz or .mol molecules.

	///		Read program arguments

	if (argc < 2)
	{
		// TODO Write proper help
		cout << "Error: input format is not correct" << endl;
		cout << "./single_mol_analysis [options] filename" << endl;
		return -1;
	}
	
	char Filename[MAX_FILENAME];

	strncpy(Filename, argv[argc-1], MAX_FILENAME);

	///		Build file name

	char Extension[MAX_FILENAME];
	char Basename[MAX_FILENAME];

	char PointsFilename[MAX_FILENAME];

	ExtractBasename(Filename, Basename, Extension);

	strncpy(Filename, argv[argc-1], MAX_FILENAME);


	///		Define input extension

	int InputExt = INPUT_EXT_NONE;

	if (strncmp(Extension, "mol", MAX_FILENAME) == 0)
	{
		InputExt = INPUT_EXT_MOL;
	}
	if (strncmp(Extension, "frameid", MAX_FILENAME) == 0)
	{
		InputExt = INPUT_EXT_FRAMEID;
	}
	if (strncmp(Extension, "sdf", MAX_FILENAME) == 0)
	{
		InputExt = INPUT_EXT_SDF;
	}

	///		Process user's options
	for (int i = 1; i < (argc-1); i++)
	{
		int Option = ProcessOption(argv[i]);

		if (Option == IE_ANALYSIS)
		{
			strncpy(PointsFilename, argv[i+1], MAX_FILENAME);
			i++;
		}
	}


	if (InputExt == INPUT_EXT_FRAMEID)
	{

		cout << "Processing " << Filename << " file ---> This will create a .mol file, execute again for analysis!" << endl;

		vector <Complex *> MoleculeChemList;

		if (!ReadFrameinfoFile(MoleculeChemList, Filename) )
		{
			cout << "Error reading " << Filename << ": file not found" << endl;
			return -1;
		}

		//Complex *Chem = MoleculeChemList.at(0);

		char MolFilename[MAX_FILENAME];

		//getwd(MolFilename);

		strncpy(MolFilename, Basename, MAX_FILENAME);

		//strncat(MolFileName, "_0", MAX_FILENAME);

		strncat(MolFilename, ".mol", MAX_FILENAME);

		cout << "File created: " << MolFilename << endl;

		//PrintToMolFile_Molecule(Chem, MolFilename);

		FrameidToMol(Filename, MolFilename);

	}


	////////////////////////////////////////////////////
	/////	MAIN COMPUTATIONS WITH STANDARD OUTPUT /////
	////////////////////////////////////////////////////


	///		Single molecule computing
	if (InputExt == INPUT_EXT_MOL)
	{

		///		LOAD MOLECULE FROM FILE

		LoadMolecule(Filename, ExecutionOptions);

		PrintOutput(Basename, ExecutionOptions);

		ReleaseMolecule();

	}

	///		Multiple molecule computing
	if (InputExt == INPUT_EXT_SDF && ExecutionOptions%MOLS_OVERLAP != 0 && ExecutionOptions%IE_ANALYSIS != 0)
	{

		bool sdfl = (ExecutionOptions%SDFL) == 0;

		int NMolecules = LoadSDFMolecule(Filename, sdfl);

		//cout << "SDF_Read" << endl;

		for (int i = 0; i < NMolecules; i++)
		{
			//cout << "Analyzing new molecule " << endl;

			bool CM = SetCurrentMolecule(i, ExecutionOptions);

			if (CM) {
				PrintOutput(Basename, ExecutionOptions);
			}

		}

		ReleaseSDFChemicals();

		//cout << "Number of molecules read: " << NMolecules << endl;

	}

	if (ExecutionOptions%MOLS_OVERLAP == 0)
	{
		//	TODO - CHECK SDF

		bool sdfl = (ExecutionOptions%SDFL) == 0;

		int NMolecules = LoadSDFMolecule(Filename, sdfl);

		if (MoleculesOverlap())
		{
			cout << Basename << " " << "1" << endl;
			//cout << Filename << endl;
		}
		else
		{
			cout << Basename << " " << "0" << endl;

		}

		ReleaseSDFChemicals();

	}


	if (ExecutionOptions%IE_ANALYSIS == 0)
	{
		bool sdfl = (ExecutionOptions%SDFL) == 0;

		int NMolecules = LoadSDFMolecule(Filename, sdfl);

		list <Point *> PointList;

		cout << PointsFilename << endl;

		ReadPoints(PointsFilename, PointList);

		//cout << "Points read" << endl;

		PointsInMolecules(PointList);

		ReleaseSDFChemicals();
	}


	return 0;

}
