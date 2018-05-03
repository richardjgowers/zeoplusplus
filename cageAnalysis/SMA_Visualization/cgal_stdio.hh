//	Software for single molecule analysis
//
//	cgal_stdio.hh
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : February 1st 2016


#ifndef _CGAL_STDIO_H_
#define _CGAL_STDIO_H_

// C++ Libraries
#include <iostream>

using namespace std;

// SMA Libraries
#include "../SMA_Core/cgal_general.hh"

// Functions

void PrintAlphaCells (Alpha_shape_3 &MoleculeAlphaShape, Classification_type CType, double Alpha);
void PrintAlphaFacets (Alpha_shape_3 &MoleculeAlphaShape, Classification_type CType, double Alpha);
void PrintAlphaEdges (Alpha_shape_3 &MoleculeAlphaShape, Classification_type CType, double Alpha);
void PrintAlphaVertices (Alpha_shape_3 &MoleculeAlphaShape, Classification_type CType, double Alpha);

#endif
