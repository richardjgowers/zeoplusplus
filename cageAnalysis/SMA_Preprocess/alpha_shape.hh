/*
 * alpha_shape.hh
 *
 *  Created on: Apr 4, 2016
 *      Author: ismael.gomez
 */

#ifndef ALPHA_SHAPE_HH_
#define ALPHA_SHAPE_HH_

//#include "cgal_general.hh"

//#include "../SMA_Core/geometry.hh"

#include "../SMA_Core/molecule_info.hh"

Complex *AlphaShapeToComplex (	MoleculeInfo *Molecule);

Complex *AlphaShapeToComplex_ExtraPoints (	MoleculeInfo *Molecule, list <Point *> &ExtraPoints);

double AlphaSelectionAllBonds( 	Complex *Chemical,
								Alpha_shape_3 &AlphaShape);


#endif /* ALPHA_SHAPE_HH_ */
