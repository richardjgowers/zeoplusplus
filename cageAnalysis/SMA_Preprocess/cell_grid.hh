/*
 * cell_grid.hh
 *
 *  Created on: Jul 1, 2016
 *      Author: ismael.gomez
 */

#ifndef CELL_GRID_HH_
#define CELL_GRID_HH_

#include "../SMA_Core/molecule_info.hh"

void MoleculeInternalPoints (	Complex *Chemical,
								list <Point *> *PointList);

void InternalPointsGrid (	MoleculeInfo *Molecule,
							list <Point *> *PointList);

void CellCenterGrid (	Complex *Chemical,
						list <Point *> *PointList);

void SphericalGridInternalCorrection(MoleculeInfo *Molecule);

#endif /* CELL_GRID_HH_ */
