/*
 * preprocess.cc
 *
 *  Created on: Mar 3, 2016
 *      Author: ismael.gomez
 */

#include "preprocess.hh"

/////////////////


MoleculeInfo *PerformMCA (Complex *Chemical, int ExecMode)
{

	/////	ERROR CHECK
	///		If Chemical is null, stop

	if (Chemical == NULL) {
		return NULL;
	}

	/////	CREATE MOLECULE_INFO OBJECT
	///		Empty object
	//		Add pruned chemical to it

	MoleculeInfo *Molecule = new MoleculeInfo();

	/////	COMPUTE CHEMICAL INFORMATION
	///		Prune chemical (previously read from file)
	//		Insert into molecule


	if (ExecMode%NO_PRUNE == 0)
	{
		//cout << "Not pruned" << endl;
		Complex *ChemCopy = Chemical->ComplexCopy();
		Molecule->SetChemical(ChemCopy);
	}
	else
	{
		//cout << "Pruned" << endl;

		Complex *PrunedChemical = Chemical->PruneComplex();

		if (PrunedChemical == NULL)
		{
			delete Molecule;
			return NULL;
		}

		Molecule->SetChemical(PrunedChemical);
	}


	/////	COMPUTE VORONOI GRAPH
	///		Use molecule's chemical (after processing)
	//		Insert into molecule info
	Complex *MolVoronoi = ComputeVoronoi( Molecule->GetChemical() );

	Molecule->SetVoronoiGraph(MolVoronoi);

	Complex *VoroCH = MolVoronoi->ComplexCopy();

	Molecule->SetVoroConvexHull(VoroCH);


	/////	MOLECULE'S CHARACTERISTIC
	///
	//

	///		Test (exploratory mode)
	bool Test = false;
	if (ExecMode%EXPLORATORY == 0)
	{
		Test = CenterShadowTest(Molecule);
		return Molecule;
	}
	//

	//Test = CenterShadowTest(Molecule);

	/////	COMPUTE INTERNAL CELLS
	///		Stop processing if internal cells computing failed
	//

	/*if( !MoleculeInternalCells(Molecule) )
	{
		delete Molecule;
		return NULL;
	}*/

	/////	COMPUTE SHADOW CHARACTERISTIC FOR BEST VORONOI NODE
	///		Only if previous test failed

	/*if (ExecMode%EXPLORATORY == 0 && !Test)
	{
		Test = MaxIntSphereShadowTest(Molecule);
		if (!Test)
		{
			return Molecule;
		}
	}*/


	CenterShadowTest(Molecule);


	int SRC = ProcessShadow(Molecule);

	//	If no points were reclassified as internal, finish
	if (SRC == 0)
	{
		//delete Molecule;
		//return NULL;
		return Molecule;
	}

	Complex *VoroCluster = ClusterizeComplex(Molecule->GetVoronoiGraph(), Molecule->GetChemical(), ATOMIC_CLUSTER_THRESHOLD);
	delete Molecule->GetVoronoiGraph();
	Molecule->SetVoronoiGraph(VoroCluster);


	/////	COMPUTE MOLECULE WINDOWS
	///		Call function and check if it succesfully computed
	//

	bool CWFlag = true;
	if ( ExecMode%WINDOW_COMP == 0)
	{
		//CWFlag = ComputeWindows( Molecule, ID_WC_ALL);
		//CWFlag = ComputeWindows( Molecule, UW_ALPHA_OPT);

		CWFlag = ComputeWindowsPER( Molecule);

		//if (!CWFlag) cout << "Warning: Windows were not successfully computed" << endl;
		Molecule->SetWinCorrect(CWFlag);
	}

	///// 	RETURN MOLECULE_INFO OBJECT
	///		Containing:
	//		- Chemical info (pruned)
	//		- Voronoi tessellation
	//		- Window list and window number

	return Molecule;

}

MoleculeInfo *MolecularCageAnalysis(char *Filename, int ExecMode)
{

	// TODO: Parse filename

	/////	READ CHEMICAL INFORMATION FROM FILE
	///		Read it from file (.mol).
	//		Store it in Complex object.

	Complex *Chemical = new Complex();

	if (!ReadMolFile(Chemical, Filename))
	{
		//cout << "Error reading file" << endl;
		delete Chemical;
		return NULL;
	}

	// TODO check extension of the file and load in depends
	/*if (!ReadXYZFile(Chemical, Filename))
	{
		cout << "Error reading file" << endl;
		return NULL;
	}*/

	/////	PERFORM MOLECULAR CAGE ANALYSIS
	///		Retrieve molecular relevant information in a MoleculeInfo object.
	//		Allocation performed at a lower layer.

	MoleculeInfo *Molecule = PerformMCA(Chemical, ExecMode);

	delete Chemical;

	return Molecule;

}

void MultipleMolecularCageAnalysis (	char *Filename,
										vector <MoleculeInfo *> *MMCA_Vector,
										int ExecMode)
{

	// TODO: Check that the extension is appropriate

	/////	LOAD ALL MOLECULES FROM FILE
	///		Molecules stored in vector of Complex *
	vector <Complex *> *MoleculeVector = new vector <Complex *>();

	//ReadFrameinfoFile(MoleculeVector, Filename);
	ReadSDF(MoleculeVector, Filename);

	//cout << "Molecules read: " << MoleculeVector->size() << endl;

	/////	FOR EACH MOLECULE, PROCESS
	///		Perform all computations allowed by the software
	//

	vector <Complex *>::iterator MolIt;

	int Count = 0;

	//cout << "Mols read: " << MoleculeVector->size() << endl;

	for (MolIt = MoleculeVector->begin(); MolIt != MoleculeVector->end(); MolIt++)
	{
		Complex *Chemical = *(MolIt);

		//char *CName = Chemical->GetCName();

		//cout << "Processing " << CName << endl;

		//	TODO - RESTORE (with ExecMode argument)
		MoleculeInfo *Molecule = PerformMCA(Chemical, ExecMode);

		if (Molecule != NULL)
		{
			MMCA_Vector->push_back(Molecule);
			Count++;
		}
		//else
		//{
			//cout << "Molecule failed to be read" << endl;
		//}

	}

	///		ERASE MOLECULES (After analysis)
	for (MolIt = MoleculeVector->begin(); MolIt != MoleculeVector->end(); MolIt++)
	{
		Complex *Chemical = *(MolIt);

		delete Chemical;
	}

	delete MoleculeVector;


}

/************
 *
 */






