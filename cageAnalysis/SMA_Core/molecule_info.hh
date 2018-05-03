//	Software for single molecule analysis
//
//	complex_set.hh
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : February 5th 2016

#ifndef _MOLECULE_INFO_H_
#define _MOLECULE_INFO_H_

#include "complex.hh"
#include "window.hh"

#include <iostream>

#define MAX_WINDOWS	10000

class MoleculeInfo
{

public:

	MoleculeInfo(Complex *Chemical, Window **Windows, Complex *VoronoiGraph, Complex *VoroConvexHull, int NWindows, bool WinCorrect, double ShadowChar);
	MoleculeInfo();
	~MoleculeInfo();

	// Getters

	Complex *GetChemical();
	Complex *GetVoronoiGraph();
	Complex *GetVoroConvexHull();
	Window 	**GetWindows();
	int 	GetNWindows();
	bool 	GetWinCorrect();
	double 	GetShadowCharact();

	Point *GetGreatestVoroNode(PointType Type);

	// Setters

	void SetChemical(Complex *Chemical);
	void SetWindows(Window **Windows);
	void SetVoronoiGraph(Complex *VoronoiGraph);
	void SetVoroConvexHull(Complex *VoroConvexHull);
	void SetNWindows(int NWindows);
	void SetWinCorrect(bool WinCorrect);
	void SetShadowCharact(double ShadowCharact);

	// Additional methods
	void InsertWindow(Window *WindowT);
	void InsertWindow(Window *WindowT, bool NonCyclic);

	void DeleteWindows();

	void PrintWindows();

	bool WindowIn(Window *WindowT);

private:

	// Atoms and bonds information
	Complex *Chemical;

	// Windows of the molecule (external)
	int NWindows;
	Window **Windows;
	bool WinCorrect;

	// Voronoi graph of the molecule
	Complex *VoronoiGraph;

	Complex *VoroConvexHull;

	//	Molecule's shadow characteristic
	double ShadowCharact;

};

#endif
