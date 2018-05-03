//	Software for single molecule analysis
//
//	windows.hh
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : January 25th 2016


#ifndef _CGAL_INTERFACE_H_
#define _CGAL_INTERFACE_H_

#define WIN_MAX_ITER	100

///		Maximum number of iterations before checking if all windows were computed (moving alphas methods)
#define MAX_ITER_CHECK	100

///		Identifiers for window computation method selection
#define F_ALPHA_W0		 	2
#define	VORO_GRID			3
#define	MV_ALPHA_W			5
#define	MV_ALPHA_UW			7
#define GEN_MESH			11
#define UW_ALPHA_OPT		13
#define ID_WC_ALL			2*3*5*7*11*13

// SMA Libraries
//#include "windows_aux.hh"


//#include "../SMA_Core/molecule_info.hh"
#include "cell_grid.hh"

//#include "../SMA_Interface/cgal_stdio.hh"

using namespace std;

using namespace boost;

/* 	ComputeWindows
 *		Calculates the windows of the molecule through its alpha shape, storing the information into the MoleculeInfo object pointed by the only argument of the function.
 *
 *	Input:
 *		Molecule - Pointer to MoleculeInfo object containing complete information about the molecule (chemical, geometrical, graphical). Window info is added in this stage.
 *	Return:
 *		(In Molecule object) - Complex with added information: holes (nodes, internal and external) and channels between holes.
 */
bool ComputeWindows( MoleculeInfo *Molecule);

Complex *GetExpandedVoronoi (	Complex *VoroGraph);


/* 	VoroEntryPaths
 *		Compute the number of entry paths to the molecule based on the Voronoi edges and internal/external Voronoi nodes identification based on alpha shapes (computation performed at internal_cells.cc).
 *
 *	Input:
 *		Molecule - Pointer to molecule's chemical information.
 *	Return:
 *		(In Molecule object) - Complex with added information: holes (nodes, internal and external) and channels between holes.
 */
int VoroEntryPaths( Complex *MolVoronoi);


#endif
