/*
 * node_characteristic.hh
 *
 *  Created on: Nov 7, 2016
 *      Author: ismael.gomez
 */

#ifndef NODE_CHARACTERISTIC_HH_
#define NODE_CHARACTERISTIC_HH_

#include "../SMA_Core/molecule_info.hh"

#include "complex_io.hh"

#include <time.h>

//#include <list>

///		Criterion ranges

#define RELATIVE_SURFACE_LIMIT	0.35
//#define RELATIVE_SURFACE_LIMIT	0.18
//#define RELATIVE_SURFACE_LIMIT	0.7


///		Functions

Window *SphericGrid(double Center[DIM_SMA], double Radius);

Window *IcosahedronSphericGrid (	double Center[DIM_SMA],
									double Radius,
									int NDivisions);

int ConnectedComponents (	Window *Characteristic,
							vector <Window *> &Components);

Window *NodeCharacteristic(	Point *Node,
							MoleculeInfo *Molecule,
							Window *Boundary,
							double &GridSurface);

int ProcessShadow (	MoleculeInfo *Molecule);

bool CenterShadowTest(	MoleculeInfo *Molecule);

bool MaxIntSphereShadowTest (MoleculeInfo *Molecule);

int PointListShadowTest (	MoleculeInfo *Molecule,
							list <Point *> &PointList,
							bool Reclassify);




#endif /* NODE_CHARACTERISTIC_HH_ */
