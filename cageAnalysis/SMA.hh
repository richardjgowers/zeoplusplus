/*
 * SMA.hh
 *
 *  Created on: Mar 4, 2016
 *      Author: ismael.gomez
 */

#ifndef SMA_HH_
#define SMA_HH_

#include <string.h>
#include <vector>

//#include "SMA_Analysis/internal_cells.hh"
//#include "SMA_Analysis/volume.hh"
//#include "SMA_Analysis/window_comp.hh"

#include "SMA_Postprocess/volume.hh"
#include "SMA_Postprocess/dijkstra_path.hh"

#include "SMA_Postprocess/surface_area.hh"

#include "SMA_Postprocess/internal_external_class.hh"

#include "SMA_Preprocess/preprocess.hh"
#include "SMA_Preprocess/complex_io.hh"

#include "SMA_Preprocess/alpha_shape.hh"

#define MAX_FILENAME 1000

// TODO - Make it unaccessible
void IntToString (int Number, char *String);

void LoadMolecule(char *Filename, int ExecMode);

void LoadMolecule(char *Filename, MoleculeInfo *NewMolecule);

void ReleaseMolecule();

char *MoleculeName();

char *MoleculeName(char *Name);

/****************************************/
/****** INTERNAL VOLUMES AND CELLS ******/
/****************************************/

double MaximumInternalSphere();

void MaximumInternalSphereCoordinates(double *Coordinates);

double TotalInternalVolume ();

double TotalConvexHullVolume ();

void InternalCells(vector <double *> &Cells);

double ShadowCharacteristic();

double InternalSurfaceArea();

/***************************************/
/*************** WINDOWS ***************/
/***************************************/

bool NumberOfWindows(int &NWindows);

bool NumberOfWindows(int &NWindows, int WinNTriangles);

int NumberOfEntryPaths();

void MoleculeNthWindow (	int Index,
							std::vector<double *> &WindowPoints,
							double &WindowSurface);

double TotalSurface ();


/***************************************/
/******* ACCESS ROUTES AND SIZES *******/
/***************************************/

double MaximumAccessibleSize();

double TargetAccessibleSize(int CellId);


/***************************************/
/********* VISUALIZATION TOOLS *********/
/***************************************/

void PrintVoronoiGraphMolecule(char *FileBaseName);

void PrintChemical (char *FileBaseName);

void PrintChemicalSDF (char *Location);

void PrintWindowsToFile(char *Basename);

void PrintAlphaShape (char *Basename);

void PrintAlphaShape_Extended(char *Basename);


/***************************************/
/************ SDF HANDLING *************/
/***************************************/

int LoadSDFMolecule (char *Filename, bool SDF_Loosely);

void ReleaseSDFChemicals();

bool SetCurrentMolecule (int Position, int ExecMode);

/*************************************************/
/******* POINTS INSIDE MOLECULES STUDIES   *******/
/*************************************************/

bool MoleculesOverlap();

void PointsInMolecules(list <Point *> &PointList);


#endif /* SMA_HH_ */
