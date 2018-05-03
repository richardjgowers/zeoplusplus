/*
 * preprocess.hh
 *
 *  Created on: Mar 3, 2016
 *      Author: ismael.gomez
 */

#ifndef PREPROCESS_HH_
#define PREPROCESS_HH_

#include "complex_io.hh"

//#include "../SMA_Core/molecule_info.hh"


#include "alpha_shape.hh"

#include "internal_cells.hh"
#include "window_comp.hh"
//#include "node_characteristic.hh"
#include "window_PER.hh"
#include "voro_interface.hh"

//#define EXEC_MODE_NONE		0		//	No execution: just show usage help
//#define EXEC_MODE_EXPLORE	1			//	Execution mode for db exploration: internal volumes only
//#define EXEC_MODE_DETAILED	2		//	Execution mode for detailed molecule description: window computing

///		Execution options (by default only internal regions are computed)
#define EXEC_MODE_NONE	1
#define WINDOW_COMP		2
#define	PRINT_OUTPUT	3
#define EXPLORATORY		5
#define NO_PRUNE		7
#define MOLS_OVERLAP	11
#define IE_ANALYSIS		13
#define SDFL			17

//MoleculeInfo *Preprocess(	char *Filename);

MoleculeInfo *MolecularCageAnalysis(char *Filename, int ExecMode);

void MultipleMolecularCageAnalysis (	char *Filename,
										vector <MoleculeInfo *> *MMCA_Vector,
										int ExecMode);

MoleculeInfo *PerformMCA (Complex *Chemical, int ExecMode);


#endif /* PREPROCESS_HH_ */
