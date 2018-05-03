/*
 * internal_external_class.hh
 *
 *  Created on: Feb 27, 2017
 *      Author: ismael
 */

#ifndef INTERNAL_EXTERNAL_CLASS_HH_
#define INTERNAL_EXTERNAL_CLASS_HH_

#include "../SMA_Preprocess/internal_cells.hh"

#include "../SMA_Preprocess/node_characteristic.hh"

//#include "../SMA_Core/molecule_info.hh"

int IEClas_AlphaShape_AtomsIntSRVoro(	list <Point *> &PointList,
										MoleculeInfo *Molecule,
										bool Reclassify);

bool IEClas_MoleculesOverlap(	Complex *Chemical1,
								Complex *Chemical2);

int IEClas_PointsInMolecule( 	Complex *Chemical,
								list <Point *> &PointList);

#endif /* INTERNAL_EXTERNAL_CLASS_HH_ */
