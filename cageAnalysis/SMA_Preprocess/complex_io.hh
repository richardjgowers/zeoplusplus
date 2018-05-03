//	Software for single molecule analysis
//
//	complex_io.hh
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : January 25th 2016

#include <iostream>
#include <fstream>
#include <vector>

#include <string.h>

#include <stdlib.h>

//#include "../SMA_Core/complex.hh"
#include "../SMA_Core/molecule_info.hh"

#include "chem_info.hh"

#define PINFO_SZ (DIM_SMA+1)

#define SDF_EOF				0
#define	SDF_OK	 			1
#define SDF_READ_ERR 		-1
#define SDF_READ_ERR_EOF 	-2


using namespace std;

bool ReadMoleculeInfo(	char *Filename,
						Complex *Molecule);

bool PrintToXYZFile (	Complex *TargetComplex,
						char *Filename);

bool PrintToMolFile (	Complex *TargetComplex,
						char *Filename);

bool PrintToMolFile (	vector <Complex *> &TargetComplexes,
						char *Filename);

bool PrintToMolFile_Molecule (	Complex *Chemical,
								char *Filename);

bool PrintToVTK (	Complex *TargetComplex,
					char *Filename);

bool PrintWindowToVTK (	Window *TWin,
						char *Basename);

bool ReadMolFile(	Complex *MoleculeChem,
					char *Filename);

bool ReadXYZFile(	Complex *MoleculeChem,
					char *Filename);

bool ReadXYZPoints (	list <Point *> *Points,
						char *Filename);

bool ReadFrameinfoFile(	vector <Complex *> &MoleculeChemList,
						char *Filename);

bool FrameidToMol (	char *FrameidFilename,
					char *MolFilename);

bool ReadPoints( 	char *Filename,
					list <Point *> &PointList);

bool ReadSDF (	vector <Complex *> *Chemicals,
				char *Filename);

bool ReadSDFL (	vector <Complex *> *Chemicals,
				char *Filename);
