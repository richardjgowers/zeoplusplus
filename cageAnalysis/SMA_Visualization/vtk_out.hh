/*
 * vtk_out.hh
 *
 *  Created on: Mar 6, 2017
 *      Author: ismael.gomez
 */

#ifndef VTK_OUT_HH_
#define VTK_OUT_HH_

#include "../SMA_Core/molecule_info.hh"

void PrintVoronoiGraphMolecule(char *FileBaseName);

void PrintChemical (char *FileBaseName);

void PrintChemicalSDF (char *Location);

void PrintWindowsToFile(char *Basename);

void PrintAlphaShape (char *Basename);

#endif /* VTK_OUT_HH_ */
