// Voro++, a 3D cell-based Voronoi library
//
// voro_interface.hh
//
// Routines for interfacing with voro++ classes
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 23rd 2015

#ifndef _VORO_INTERFACE_H_
#define _VORO_INTERFACE_H_

#include <string.h>

#include <cmath>

#include "../voro/src/voro++.hh"

//#include "/Users/ismael.gomez/Documents/Software/Voro++//voro++-0.4.6/src/v_compute.hh"

//#include "../SMA_Core/complex.hh"
#include "../SMA_Core/molecule_info.hh"


using namespace voro;

#define FILENAME_SIZE 1000
#define POV_EXT ".pov"
#define GNU_EXT ".gnu"

#define RADIUS_COORD 3


Complex *ComputeVoronoi(			Complex *Chemical);

Complex *ComputeVoronoiFromFile(	double x_boundary[2],
									double y_boundary[2],
									double z_boundary[2],
									char *Filename);

#endif
