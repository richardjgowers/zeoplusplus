/*
 * volume.hh
 *
 *  Created on: Mar 3, 2016
 *      Author: ismael.gomez
 */

#ifndef VOLUME_HH_
#define VOLUME_HH_

#include "../SMA_Core/molecule_info.hh"

/* LargestInternalSphere
 *	Searches for the largest internal node (in terms of radius), returning that value.
 *
 * Inputs:
 * 		MolVoronoi - Voronoi tessellation of the molecule of interest.
 *
 * Outputs:
 * 		Radius of the Voronoi node with maximum radius (biggest cell's index).
 */
double LargestInternalSphere (	Complex *MoleculeC);


/* LargestInternalNodeIndex
 *
 * Inputs:
 * 		MolVoronoi - Voronoi tessellation of the molecule of interest.
 *
 * Outputs:
 * 		Index of the Voronoi node with maximum radius (biggest cell's index).
 */
int LargestInternalNodeIndex(	Complex *MolVoronoi);

void LargestInternalNodeCoordinates(	Complex *MolVoronoi,
										double *Coordinates);


double VoroTotalVolume(	Complex *MolVoronoi);

#endif /* VOLUME_HH_ */
