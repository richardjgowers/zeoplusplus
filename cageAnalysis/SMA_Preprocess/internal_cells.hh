//	Software for single molecule analysis
//
//	internal_cells.hh
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : February 15th 2016

//#include "../SMA_Core/molecule_info.hh"

#include "clustering.hh"

#include "alpha_shape.hh"

#define		INTERNAL_CRIT_OPT			0
#define		INTERNAL_CRIT_CONV_HULL		1

#define		ATOMIC_CLUSTER_THRESHOLD 	0.85

/* Function: 
 *
 */
bool MoleculeInternalCells (	MoleculeInfo *Molecule);

/*
 *
 */
bool AreInternalPoints(	list <Point *> &PointList,
						Complex *Chemical,
						bool FastAnswer);


