// Voro++, a 3D cell-based Voronoi library
//
// dijkstra_path.hh
//
// Routines for single molecule analysis
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015

//#include "../SMA_Core/complex.hh"
//#include "../SMA_Core/geometry.hh"


#include "dijkstra_aux.hh"

#include <stdlib.h>
#include <stdio.h>

#include <vector>

/* ClasifyPoins
 *
 * Checks the given complex (with the information of one molecule trapped in a cage).
 * Updates the information about the types of points deciding if they are internal or external.
 */
void ClassifyPoints(	Complex *MoleculeC);

void ClassifyPointsRestrictionCriterion (Complex *MoleculeC);


/********** FUNCTIONS FOR CALCULATION OF CELL ACCESIBLE SIZE **********/

/*double GetChannelSize(	Complex *MoleculeC,
						int SourcePointId,
						int TargetPointId);

void VisitNode(	int PointId,
				Complex *MoleculeC,
				PriorityQueue **PQList,
				VisitedList *VList,
				WeightList *WList);*/

/* CellAccessibleSize
 *
 * Checks a particular cell inside the given complex (with the information of one molecule trapped in a cage).
 * Returns the radius of the maximum sphere capable of escaping from that cell to the outside.
 */
double CellAccessibleSize (	Complex *MoleculeC,
							Point *TargetCell);


/*
 *
 */
double MoleculeAccessibleSize(Complex *MolVoronoi);


/* CellEscapePath
 *
 * Checks a particular cell inside the given complex (with the information of one molecule trapped in a cage).
 * Returns the escape path for the biggest molecule trapped in that cell.
 */
Point **CellEscapePath (	Complex *MoleculeC,
							Point *TargetCell);

/* CellProbeEscapePath
 *
 * Checks a particular cell inside the given complex (with the information of one molecule trapped in a cage).
 * Returns the escape path for a probe of given size trapped in that cell (Size_probe < Size_cell).
 */
Point **CellProbeEscapePath (	Complex *MoleculeC,
								Point *TargetCell);

