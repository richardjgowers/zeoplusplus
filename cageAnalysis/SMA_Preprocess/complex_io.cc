//	Software for single molecule analysis
//
//	complex_io.cc
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : January 25th 2016

#include "complex_io.hh"


/***********************************/
/****** PRINT COMPLEX TO FILE ******/
/***********************************/

bool PrintToXYZFile (		Complex *TargetComplex,
							char *Filename)
{

	FILE *f = fopen(Filename, "w");

	if (!f) return false;	

	fprintf(f, "\t%d\n\n", TargetComplex->GetNPoints());

	for (int i = 0; i < TargetComplex->GetNPoints(); i++)
	{
		Point *PAux = TargetComplex->GetPointList()[i];
		
		//fprintf(f, "O %f %f %f\n", PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2) );
		fprintf(f, "C\t%f\t%f\t%f\n", PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2) );

	}

	fclose(f);

	return true;


}

bool PrintToMolFile (	Complex *TargetComplex,
						char *Filename)
{

	FILE *f = fopen(Filename, "w");

	if (!f) return false;

	/////	MOL FILE HEADER
	//

	fprintf(f, "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n");
	fprintf(f, "M  V30 BEGIN CTAB\nM  V30 COUNTS %d %d 0 0 0\nM  V30 BEGIN ATOM\n", TargetComplex->GetNPoints(), TargetComplex->GetNEdges());

	/////	PRINT POINTS
	//

	for (int i = 0; i < TargetComplex->GetNPoints(); i++)
	{
		Point *PAux = TargetComplex->GetPointList()[i];

		// TODO Print the adequate element
		char PT = '\0';

		switch(PAux->GetType())
		{
			case internal_point:
				PT = 'O';
				break;
			case external_point:
				PT = 'C';
				break;
			case boundary_point:
				PT = 'H';
				break;
			case undefined_point:
				PT = 'N';
				break;
			case atom_point:
				PT = 'F';
				break;
			case hybrid_point:
				PT = 'S';
				break;
			default:
				break;
		}

		fprintf(f, "M  V30 %d %c %lf %lf %lf 0\n", PAux->GetPointId()+1, PT ,PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2) );

	}

	/////	PRINT EDGES
	//
	
	fprintf(f, "M  V30 END ATOM\nM  V30 BEGIN BOND\n");
	
	for (int i = 0; i < TargetComplex->GetNEdges(); i++)
	{
		Edge *EAux = TargetComplex->GetEdgeList()[i];
		
		fprintf(f, "M  V30 %d 1 %d %d\n", i+1, EAux->GetPSource()->GetPointId()+1, EAux->GetPTarget()->GetPointId()+1);
	}

	fprintf(f, "M  V30 END BOND\nM  V30 END CTAB\nM  END\n");

	

	fclose(f);

	return true;

}

bool PrintToMolFile (	vector <Complex *> &TargetComplexes,
						char *Filename)
{

	FILE *f = fopen(Filename, "w");

	if (!f) return false;

	/////	MOL FILE HEADER
	//

	fprintf(f, "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n");

	int 	TotalPoints = 1,
			TotalEdges = 1;

	vector <int> MolBaseInd_Points, MolBaseInd_Edges;

	MolBaseInd_Points.push_back(1);
	MolBaseInd_Edges.push_back(1);

	for (int j = 0; j < TargetComplexes.size(); j++)
	{
		TotalPoints += TargetComplexes.at(j)->GetNPoints();
		TotalEdges += TargetComplexes.at(j)->GetNEdges();

		MolBaseInd_Points.push_back( TotalPoints );
		MolBaseInd_Edges.push_back(TotalEdges);

	}

	fprintf(f, "M  V30 BEGIN CTAB\nM  V30 COUNTS %d %d 0 0 0\nM  V30 BEGIN ATOM\n", TotalPoints, TotalEdges);

	/////	PRINT POINTS
	//

	for (int j = 0; j < TargetComplexes.size(); j++)
	{
		Complex *TargetComplex = TargetComplexes.at(j);

		for (int i = 0; i < TargetComplex->GetNPoints(); i++)
		{
			Point *PAux = TargetComplex->GetPointList()[i];

			// TODO Print the adequate element
			char PT = '\0';

			switch(PAux->GetType())
			{
				case internal_point:
					PT = 'O';
					break;
				case external_point:
					PT = 'C';
					break;
				case boundary_point:
					PT = 'H';
					break;
				case undefined_point:
					PT = 'N';
					break;
				case atom_point:
					PT = 'F';
					break;
				case hybrid_point:
					PT = 'S';
					break;
				default:
					break;
			}

			fprintf(f, "M  V30 %d %c %lf %lf %lf 0\n", PAux->GetPointId() + MolBaseInd_Points.at(j) , PT ,PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2) );

		}
	}

	/////	PRINT EDGES
	//

	fprintf(f, "M  V30 END ATOM\nM  V30 BEGIN BOND\n");

	for (int j = 0; j < TargetComplexes.size(); j++)
	{
		Complex *TargetComplex = TargetComplexes.at(j);

		for (int i = 0; i < TargetComplex->GetNEdges(); i++)
		{
			Edge *EAux = TargetComplex->GetEdgeList()[i];

			fprintf(f, "M  V30 %d 1 %d %d\n", i + MolBaseInd_Edges.at(j), EAux->GetPSource()->GetPointId() + MolBaseInd_Points.at(j), EAux->GetPTarget()->GetPointId() + MolBaseInd_Points.at(j));
		}
	}


	fprintf(f, "M  V30 END BOND\nM  V30 END CTAB\nM  END\n");



	fclose(f);

	return true;

}


bool PrintToMolFile_Molecule (	Complex *Chemical,
								char *Filename)
{

	FILE *f = fopen(Filename, "w");

	if (!f) return false;

	/////	MOL FILE HEADER
	//

	fprintf(f, "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n");
	fprintf(f, "M  V30 BEGIN CTAB\nM  V30 COUNTS %d %d 0 0 0\nM  V30 BEGIN ATOM\n", Chemical->GetNPoints(), Chemical->GetNEdges());

	/////	PRINT POINTS
	//

	for (int i = 0; i < Chemical->GetNPoints(); i++)
	{
		Point *PAux = Chemical->GetPointList()[i];

		// TODO Print the adequate element
		char PT[3] = {'\0','\0','\0'};

		double Radius = PAux->GetRadius();

		if (Radius == H_RAD)
		{
			PT[0] = 'H';
		}
		else if (Radius == C_RAD)
		{
			PT[0] = 'C';
		}
		else if (Radius == O_RAD)
		{
			PT[0] = 'O';
		}
		else if (Radius == N_RAD)
		{
			PT[0] = 'N';
		}
		else if (Radius == P_RAD)
		{
			PT[0] = 'P';
		}
		else if (Radius == F_RAD)
		{
			PT[0] = 'F';
		}
		else if (Radius == B_RAD)
		{
			PT[0] = 'B';
		}
		else if (Radius == CL_RAD)
		{
			PT[0] = 'C';
			PT[1] = 'l';
		}
		else if (Radius == SI_RAD)
		{
			PT[0] = 'S';
			PT[1] = 'i';
		}
		else if (Radius == BR_RAD)
		{
			PT[0] = 'B';
			PT[1] = 'r';
		}
		else
		{
			cout << "ERROR: Type of element not considered!! - Radius: " << Radius << endl;
		}


		fprintf(f, "M  V30 %d %s %lf %lf %lf 0\n", PAux->GetPointId()+1, PT ,PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2) );

	}

	/////	PRINT EDGES
	//

	fprintf(f, "M  V30 END ATOM\nM  V30 BEGIN BOND\n");

	for (int i = 0; i < Chemical->GetNEdges(); i++)
	{
		Edge *EAux = Chemical->GetEdgeList()[i];

		fprintf(f, "M  V30 %d 1 %d %d\n", i+1, EAux->GetPSource()->GetPointId()+1, EAux->GetPTarget()->GetPointId()+1);
	}

	fprintf(f, "M  V30 END BOND\nM  V30 END CTAB\nM  END\n");



	fclose(f);

	return true;

}


bool PrintToVTK (	Complex *TargetComplex,
					char *Basename)
{

	FILE *f;

	char Filename[FILENAME_MAX];

	strncpy(Filename, Basename, FILENAME_MAX);
	strncat(Filename, ".vtk", FILENAME_MAX);

	f = fopen(Filename, "w");

	///		Write data with VTK format

	///		Header

	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "vtk data for file %s\n", Basename);
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET POLYDATA\n");


	///		Points
	int NPoints = TargetComplex->GetNPoints();

	fprintf(f, "POINTS %d double\n", NPoints);

	for (int i = 0; i < NPoints; i++)
	{
		///	Read the coordinates of the point
		double Coordinates[DIM_SMA];
		TargetComplex->GetPointList()[i]->GetCoordinates(Coordinates);

		///	Print the coordinates
		fprintf(f, "%.3f %.3f %.3f\n", Coordinates[X_COORD], Coordinates[Y_COORD], Coordinates[Z_COORD]);


	}

	///		Edges
	int NEdges = TargetComplex->GetNEdges();

	fprintf(f, "LINES %d %d\n", NEdges, 3*NEdges);

	for (int i = 0; i < NEdges; i++)
	{
		///	Read the indices of the edges
		int 	Id1 = TargetComplex->GetEdgeList()[i]->GetPSource()->GetPointId(),
				Id2 = TargetComplex->GetEdgeList()[i]->GetPTarget()->GetPointId();

		///	Print the coordinates
		fprintf(f, "2 %d %d\n", Id1, Id2);

	}


	///		Close file and return

	fclose(f);

	return true;

}


bool PrintWindowToVTK (	Window *TWin,
						char *Basename)
{

	FILE *f;

	char Filename[FILENAME_MAX];

	strncpy(Filename, Basename, FILENAME_MAX);
	strncat(Filename, ".vtk", FILENAME_MAX);

	f = fopen(Filename, "w");

	if (!f)
	{
		cout << "Error writting windows\n" << endl;
		return false;
	}


	///		Normalize window just before printing.
	TWin->Normalize();

	///		Write data with VTK format

	///		Header

	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "vtk data for file %s\n", Basename);
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET POLYDATA\n");

	///		Points
	int NPoints = TWin->GetNPoints();

	fprintf(f, "POINTS %d double\n", NPoints);

	for (int i = 0; i < NPoints; i++)
	{
		///	Read the coordinates of the point
		double Coordinates[DIM_SMA];
		TWin->GetPointList()[i]->GetCoordinates(Coordinates);

		///	Print the coordinates
		fprintf(f, "%.3f %.3f %.3f\n", Coordinates[X_COORD], Coordinates[Y_COORD], Coordinates[Z_COORD]);


	}

	///		Triangles

	int NTriangles = TWin->GetSimplices()->size();

	//cout << "N triangles" << NTriangles << endl;

	list <Simplex *>::iterator SIt;

	fprintf(f, "TRIANGLE_STRIPS %d %d\n", NTriangles, 4*NTriangles);

	for (SIt = TWin->GetSimplices()->begin(); SIt != TWin->GetSimplices()->end(); SIt++)
	{
		Simplex *Triangle = *(SIt);

		Point 	*P0 = Triangle->GetSimplexPoints()->at(0),
				*P1 = Triangle->GetSimplexPoints()->at(1),
				*P2 = Triangle->GetSimplexPoints()->at(2);

		double 	Coords0[DIM_SMA],
				Coords1[DIM_SMA],
				Coords2[DIM_SMA];

		P0->GetCoordinates(Coords0);
		P1->GetCoordinates(Coords1);
		P2->GetCoordinates(Coords2);

		int 	Id0 = TWin->GetPointByCoordinates(Coords0)->GetPointId(),
				Id1 = TWin->GetPointByCoordinates(Coords1)->GetPointId(),
				Id2 = TWin->GetPointByCoordinates(Coords2)->GetPointId();

		//cout << "Ids: " << Id0 << " " << Id1 << " " << Id2 << endl;

		fprintf(f, "3 %d %d %d\n", Id0, Id1, Id2);

	}

	//cout << "Triangles printed" << endl;

	///		Close file and return

	fclose(f);

	return true;



}


bool ReadMolFile(	Complex *MoleculeChem,	// Complex already initialized
					char *Filename)
{

	FILE *f = fopen(Filename, "r");

	if (!f) {
		cout << "Error reading file " << Filename << endl;
		return false;
	}

	/////	MOL FILE HEADER
	//

	int NPoints = 0, NEdges = 0;

	fscanf(f, "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n");
	fscanf(f, "M  V30 BEGIN CTAB\nM  V30 COUNTS %d %d 0 0 0\nM  V30 BEGIN ATOM\n", &NPoints, &NEdges);

	/*MoleculeChem->SetNPoints(NPoints);
	MoleculeChem->SetNEdges(NEdges);*/

	
	/////	READ ATOMS
	//

	if (NPoints == 0 || NEdges == 0)
	{
		fclose(f);
		return false;
	}

	for (int i = 0; i < NPoints; i++)
	{
		//Point *PAux = TargetComplex->GetPointList()[i];
		double Coordinates[DIM_SMA] = {0., 0., 0.};
		int PointId = -1;
		//char Element = ' ';
		char Element[3] = {' ', ' ', '\0'};
		
		fscanf(f, "M  V30 %d %s %lf %lf %lf 0\n", &PointId, Element, &(Coordinates[X_COORD]), &(Coordinates[Y_COORD]), &(Coordinates[Z_COORD]) );

		//cout << "Element: " << Element << " - Point Id " << PointId << " - Point: " << Coordinates[X_COORD] << " " << Coordinates[Y_COORD] << " " << Coordinates[Z_COORD] <<  endl;

		double Radius = AtomRadius(Element);

		//cout << "Radius: " << Radius << endl;

		// TODO - FIX this!! Point id must be read from file
		MoleculeChem->InsertPointByValue(i+1, Coordinates, Radius, atom_point);

	}

	/////	READ BONDS
	//
	
	fscanf(f, "M  V30 END ATOM\nM  V30 BEGIN BOND\n");
	
	for (int i = 0; i < NEdges; i++)
	{
		//Edge *EAux = TargetComplex->GetEdgeList()[i];
		int PSourceId = 0, PTargetId = 0;
		int Iter = 0, BT = 0;
		
		fscanf(f, "M  V30 %d %d %d %d\n", &Iter, &BT, &PSourceId, &PTargetId);

		Point *PSource = MoleculeChem->GetPointById(PSourceId);
		Point *PTarget = MoleculeChem->GetPointById(PTargetId);

		MoleculeChem->InsertEdgeByValue(PSource, PTarget);
	
	}


	fclose(f);

	return true;

}

bool ReadXYZFile(	Complex *MoleculeChem,
					char *Filename)
{

	FILE *f = fopen(Filename, "r");

	if (!f) return false;

	/////	XYZ FILE HEADER
	//

	int NPoints = 0;

	fscanf(f, "\n%d\n", &NPoints);

	/////	READ POINTS
	//

	for (int i = 0; i < NPoints; i++)
	{
		//Point *PAux = TargetComplex->GetPointList()[i];
		double Coordinates[DIM_SMA] = {0., 0., 0.};
		char Element = ' ';
		
		fscanf(f, "%c %lf %lf %lf 0\n", &Element, &(Coordinates[X_COORD]), &(Coordinates[Y_COORD]), &(Coordinates[Z_COORD]) );

		if (Element == 'H') {
			MoleculeChem->InsertPointByValue(Coordinates, HIDROGEN_BOND_LEN, atom_point);
		}
		else {
			MoleculeChem->InsertPointByValue(Coordinates, GENERIC_BOND_LEN, atom_point);
		}

	}

	///// 	COMPUTE BONDS
	//

	for (int i = 0; i < NPoints; i++) 
	{
		Point *P1 = MoleculeChem->GetPointList()[i];

		for (int j = i+1; j < NPoints; j++)
		{
			Point *P2 = MoleculeChem->GetPointList()[j];

			if (	Distance (P1->GetCoordinates(), P2->GetCoordinates()) < BOND_THRESHOLD + P1->GetRadius() + P2->GetRadius() ) 
			{
				MoleculeChem->InsertEdgeByValue(P1, P2);
			}
		}
	}
		

	fclose(f);

	return true;
	

}


bool ReadXYZPoints (	list <Point *> *Points,
						char *Filename)
{

	FILE *f = fopen(Filename, "r");

	if (!f) return false;

	/////	XYZ FILE HEADER
	//

	int NPoints = 0;

	fscanf(f, "\n%d\n", &NPoints);

	/////	READ POINTS
	//

	for (int i = 0; i < NPoints; i++)
	{
		//Point *PAux = TargetComplex->GetPointList()[i];
		double Coordinates[DIM_SMA] = {0., 0., 0.};
		char Element = ' ';

		fscanf(f, "%c %lf %lf %lf 0\n", &Element, &(Coordinates[X_COORD]), &(Coordinates[Y_COORD]), &(Coordinates[Z_COORD]) );

		Point *P_xyz = new Point(i, Coordinates, 0.0, undefined_point);

		Points->push_back(P_xyz);

	}

	return true;


}


bool ReadFrameinfoFile(	vector <Complex *> &MoleculeChemList,
						char *Filename)
{
	FILE *f = fopen(Filename, "r");

	if (!f) return false;

	int MolId = -1, MolIdAux = -1;

	bool EofFlag = false;


	// Storage for point info (created in this context to solve "last atom" problems
	char Element[2];
	double Useless;
	double Radius, Coordinates[DIM_SMA];

	// Read atoms information until file is finished
	while (!EofFlag)
	{
		Complex *MoleculeChem = new Complex();

		// Insert first point of the molecule into it (read in previous iteration)
		if (MolId != -1)
		{
			MoleculeChem->InsertPointByValue(Coordinates, Radius/2, atom_point);
			//MoleculeChem->InsertPointByValue(Coordinates, Radius/2, undefined_point);
		}

		// READ MOLECULE'S POINTS
		while (MolId == MolIdAux)
		{


			int ScanOut = fscanf(	f,
									"%s    %lf   %lf   %lf   %lf   %lf   %lf  %lf   %lf  %lf  %lf    %d  \n",
									Element,
									&Useless,
									&Useless,
									&Useless,
									&Radius,
									&Useless,
									&Useless,
									&Useless,
									&(Coordinates[0]),
									&(Coordinates[1]),
									&(Coordinates[2]),
									&MolId);



			// Set the aux value when reading for the first time
			if (MolIdAux == -1)
			{
				MolIdAux = MolId;
			}

			if (ScanOut == EOF)
			{
				EofFlag = true;
				break;
			}

			if (MolIdAux != MolId)
			{
				break;
			}

			MoleculeChem->InsertPointByValue(Coordinates, Radius/2, atom_point);
			/*if (Element == 'H') {
				MoleculeChem->InsertPointByValue(Coordinates, HIDROGEN_BOND_LEN, atom_point);
			}
			else {
				MoleculeChem->InsertPointByValue(Coordinates, GENERIC_BOND_LEN, atom_point);
			}*/

		}

		///// COMPUTE MOLECULE'S BONDS
		for (int i = 0; i < MoleculeChem->GetNPoints(); i++)
		{
			Point *P1 = MoleculeChem->GetPointList()[i];

			for (int j = i+1; j < MoleculeChem->GetNPoints(); j++)
			{
				Point *P2 = MoleculeChem->GetPointList()[j];

				if ( Distance (P1->GetCoordinates(), P2->GetCoordinates()) < BOND_THRESHOLD + P1->GetRadius() + P2->GetRadius() )
				{
					MoleculeChem->InsertEdgeByValue(P1, P2);
				}
			}
		}

		///// INSERT MOLECULE INTO MOLECULE LIST
		MoleculeChemList.push_back(MoleculeChem);
		MolIdAux = MolId;


	}

	fclose(f);

	return true;
}


bool FrameidToMol (	char *FrameidFilename,
					char *MolFilename)
{

	/////		READ FRAMEID FILE
	///			Store the information into Complex
	//			Information about elements is stored in vector

	FILE *f = fopen(FrameidFilename, "r");

	if (!f) return false;

	int MolId = -1, MolIdAux = -1;

	bool EofFlag = false;


	vector <char *> ElementsVector;
	// Storage for point info (created in this context to solve "last atom" problems

	double Useless;
	double Radius, Coordinates[DIM_SMA];

	///		Read molecule information: only the first one is read

	Complex *MoleculeChem = new Complex();



	// READ MOLECULE'S POINTS
	while (MolId == MolIdAux)
	{

		char *Element = new char[2];

		int ScanOut = fscanf(	f,
								"%s    %lf   %lf   %lf   %lf   %lf   %lf  %lf   %lf  %lf  %lf    %d  \n",
								Element,
								&Useless,
								&Useless,
								&Useless,
								&Radius,
								&Useless,
								&Useless,
								&Useless,
								&(Coordinates[0]),
								&(Coordinates[1]),
								&(Coordinates[2]),
								&MolId);



		// Set the aux value when reading for the first time
		if (MolIdAux == -1)
		{
			MolIdAux = MolId;
		}

		if (ScanOut == EOF)
		{
			EofFlag = true;
			break;
		}

		if (MolIdAux != MolId)
		{
			break;
		}

		MoleculeChem->InsertPointByValue(Coordinates, Radius/2, atom_point);

		//cout << " " << Element;

		ElementsVector.push_back(Element);


		for (int i = 0; i < MoleculeChem->GetNPoints(); i++)
		{
			Point *P1 = MoleculeChem->GetPointList()[i];

			for (int j = i+1; j < MoleculeChem->GetNPoints(); j++)
			{
				Point *P2 = MoleculeChem->GetPointList()[j];

				if ( Distance (P1->GetCoordinates(), P2->GetCoordinates()) < BOND_THRESHOLD + P1->GetRadius() + P2->GetRadius() )
				{
					MoleculeChem->InsertEdgeByValue(P1, P2);
				}
			}
		}


	}

	//cout << endl;

	fclose(f);

	/// 	Compute bonds


	/////	PRINT TO MOL FILE
	///		Read the information about element from ElementsVector
	//		Graphical and geometrical information read from MoleculeChem

	FILE *f2 = fopen(MolFilename, "w");

	if (!f) return false;

	/////	MOL FILE HEADER
	//

	fprintf(f2, "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n");
	fprintf(f2, "M  V30 BEGIN CTAB\nM  V30 COUNTS %d %d 0 0 0\nM  V30 BEGIN ATOM\n", MoleculeChem->GetNPoints(), MoleculeChem->GetNEdges());

	/////	PRINT POINTS
	//

	for (int i = 0; i < MoleculeChem->GetNPoints(); i++)
	{
		Point *PAux = MoleculeChem->GetPointList()[i];

		//double Radius = PAux->GetRadius();

		char *Element = ElementsVector.at(i);

		//cout << " " << Element;


		fprintf(f, "M  V30 %d %s %lf %lf %lf 0\n", PAux->GetPointId()+1, Element ,PAux->GetCoordinate(0), PAux->GetCoordinate(1), PAux->GetCoordinate(2) );

	}

	//cout << endl;

	/////	PRINT EDGES
	//

	fprintf(f, "M  V30 END ATOM\nM  V30 BEGIN BOND\n");

	for (int i = 0; i < MoleculeChem->GetNEdges(); i++)
	{
		Edge *EAux = MoleculeChem->GetEdgeList()[i];

		fprintf(f2, "M  V30 %d 1 %d %d\n", i+1, EAux->GetPSource()->GetPointId()+1, EAux->GetPTarget()->GetPointId()+1);
	}

	fprintf(f2, "M  V30 END BOND\nM  V30 END CTAB\nM  END\n");



	fclose(f2);



	return true;

}


/*********************************************/
/************	READ MDI FILES	**************/
/*********************************************/

bool ReadMDI (	vector <Complex *> &Chemical,
				char *MDIName)
{

	int LinesRead = 1;

	FILE *f = fopen(MDIName, "r");

	while (LinesRead > 0)
	{
		Complex *MDIComplex = new Complex();

		LinesRead = fscanf(f, "%s\n", MDIComplex->GetCName());

		int NPoints = 0;

		fscanf(f, "%d\n", &NPoints);

		for (int i = 0; i < NPoints; i++)
		{

			double Coordinates[DIM_SMA];

			int PointId = -1,
				MolId = 0;

			char Atom = '\0';

			fscanf(f , "    %d     %c        %d  %lf %lf %lf\n", MolId, Atom, PointId, Coordinates[X_COORD], Coordinates[Y_COORD], Coordinates[Z_COORD]);

		}



	}

	return true;

}

/*********************************************/
/************	READ XYZ POINTS	**************/
/*********************************************/

bool ReadPoints( 	char *Filename,
					list <Point *> &PointList)
{

	FILE *f = fopen(Filename, "r");

	if (!f)
	{
		cout << "Read err" << endl;
		return false;
	}

	int NR = 1;

	int PointId = 0;

	while (NR != 0)
	{
		double Coords[DIM_SMA];

		NR = fscanf(f, "%lf %lf %lf", &(Coords[X_COORD]), &(Coords[Y_COORD]), &(Coords[Z_COORD]) );

		if (NR == 3)
		{
			Point *P = new Point();

			P->SetPointId(PointId++);

			P->SetCoordinate(X_COORD, Coords[X_COORD]);
			P->SetCoordinate(Y_COORD, Coords[Y_COORD]);
			P->SetCoordinate(Z_COORD, Coords[Z_COORD]);

			//cout << "Point read" << endl;

			//P->PrintPoint();

			PointList.push_back(P);

		}
		else
		{
			break;
		}
	}

	fclose(f);

	return true;

}


/*********************************************/
/************	READ SDF FILES	**************/
/*********************************************/


/*****		STRICT SDF FILES 	*****/


bool ReadSDF_MolName (	char *SDFName,
						FILE *SDFfile)
{

	int i = 0;

	char c = '\0';

	///		Read molecule name

	while (c != '\n')
	{
		c = fgetc(SDFfile);

		if (feof(SDFfile))
		{
			return false;
		}

		//cout << "[" << c << "]";

		SDFName[i++] = c;
	}

	//if (i == 0)

	SDFName[i-1] = '\0';

	///		Skip next two lines

	int Counter = 0;

	while (Counter != 2)
	{
		c = fgetc(SDFfile);

		if (c == '\n') {
			Counter++;
		}
	}

	return true;


}

bool ReadSDF_Header (	int &NPoints,
						int &NEdges,
						FILE *SDFfile)
{

	char c = '\0';

	char 	NPointsC[4] = {'\0', '\0', '\0', '\0'},
			NEdgesC[4] = {'\0', '\0', '\0', '\0'};

	int PCount = 0,
		ECount = 0;


	while (c != '\n')
	{
		c = fgetc(SDFfile);

		if (PCount >= 3 && ECount < 3)
		{
			//cout << PCount << " " << ECount << endl;

			NEdgesC[ECount] = c;
			ECount++;
		}

		if (PCount < 3)
		{
			NPointsC[PCount] = c;
			PCount++;
		}


	}

	//cout << "NP: " << NPointsC << " NE: " <<  NEdgesC << endl;

	NPoints = atoi(NPointsC);

	//cout << NPoints << endl;

	NEdges = atoi(NEdgesC);

	//cout << NEdges << endl;

	return true;



/*	char c = '\0';

	char 	NPointsC[8] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'},
			NEdgesC[8] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};

	bool Blank = false;

	int Count = 0, LocalCC = 0;


	/// Read the line
	while (c != '\n')
	{
		c = fgetc(SDFfile);


		if (c != ' ' && Blank == true) {
			Count++;
			Blank = false;
		}

		if (c == ' ')
		{
			LocalCC = 0;
			Blank = true;
		}

		///		Store coordinates values
		if (Count == 1 && Blank == false)
		{
			NPointsC[LocalCC++] = c;
		}

		if (Count == 2 && Blank == false)
		{
			NEdgesC[LocalCC++] = c;
		}
	}

	NPoints = atoi(NPointsC);
	NEdges = atoi(NEdgesC);

	return true;*/

}

/*
 *
 */

bool ReadSDF_Point (	double *Coordinates,
						char *AtomType,
						FILE *SDFfile)
{

	char c = '\0';

	char 	Coord0[8] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'},
			Coord1[8] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'},
			Coord2[8] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};

	int CoordsCount = 0, LocalCC = 0;

	bool Blank = true;

	/// Read the line
	while (c != '\n')
	{
		c = fgetc(SDFfile);

		if (c != ' ' && Blank == true) {
			CoordsCount++;
			Blank = false;
		}

		if (c == ' ')
		{
			LocalCC = 0;
			Blank = true;
		}

		///		Store coordinates values
		if (CoordsCount == 1 && Blank == false)
		{
			Coord0[LocalCC] = c;
			LocalCC++;
		}

		if (CoordsCount == 2 && Blank == false)
		{
			Coord1[LocalCC++] = c;
		}

		if (CoordsCount == 3 && Blank == false)
		{
			Coord2[LocalCC++] = c;
		}

		if (CoordsCount == 4 && Blank == false)
		{
			AtomType[LocalCC++] = c;
		}

	}

	Coordinates[X_COORD] = atof(Coord0);
	Coordinates[Y_COORD] = atof(Coord1);
	Coordinates[Z_COORD] = atof(Coord2);

	return true;

}

bool ReadSDF_Edge (  	int *EdgeIndices,
						FILE *SDFfile,
						int NPointsHead)
{

	char c = '\0';

	char 	Id0[4] = {'\0', '\0', '\0', '\0'},
			Id1[4] = {'\0', '\0', '\0', '\0'};

	int Id0Count = 0, Id1Count = 0;

	bool Blank = true;

	/// Read the line
	while (c != '\n')
	{
		c = fgetc(SDFfile);

		if (Id1Count < 3 && Id0Count >= 3)
		{
			Id1[Id1Count] = c;
			Id1Count++;
		}

		if (Id0Count < 3)
		{
			Id0[Id0Count] = c;
			Id0Count++;
		}


	}

	EdgeIndices[0] = atoi(Id0);
	EdgeIndices[1] = atoi(Id1);

	if (EdgeIndices[0] == 0 || EdgeIndices[0] > NPointsHead) {
		return false;
	}

	if (EdgeIndices[1] == 0 || EdgeIndices[1] > NPointsHead) {
		return false;
	}

	return true;

	/*char c = '\0';

	char 	Id0[5] = {'\0', '\0', '\0', '\0', '\0'},
			Id1[5] = {'\0', '\0', '\0', '\0', '\0'};

	int IdsCount = 0, LocalCC = 0;

	bool Blank = true;

	/// Read the line
	while (c != '\n')
	{

		c = fgetc(SDFfile);

		if (c != ' ' && Blank == true) {
			IdsCount++;
			Blank = false;
		}

		if (c == ' ')
		{
			LocalCC = 0;
			Blank = true;
		}

		///		Store coordinates values
		if (IdsCount == 1 && Blank == false)
		{
			Id0[LocalCC++] = c;
		}

		if (IdsCount == 2 && Blank == false)
		{
			Id1[LocalCC++] = c;
		}

	}

	EdgeIndices[0] = atoi(Id0);
	EdgeIndices[1] = atoi(Id1);

	if (EdgeIndices[0] == 0 || EdgeIndices[0] > NPointsHead) {
		return false;
	}

	if (EdgeIndices[1] == 0 || EdgeIndices[1] > NPointsHead) {
		return false;
	}

	return true;*/


}


/*
 *
 */
bool ReadSDF_NextMol ( FILE *SDFfile)
{

	char c = '\0';

	int Counter = 0;

	while (Counter < 4)
	{

		c = fgetc(SDFfile);

		if (c == '$') {
			Counter++;
		}
		else {
			Counter = 0;
		}

	}

	// Read the last \n
	if (feof(SDFfile) != 0) {
		return true;
	}

	c = fgetc(SDFfile);

	return false;

}


/*
 *
 */
int ReadSDF_Single ( 	Complex *Chemical,
						FILE *SDFfile)
{

	///		Read molecule id

	char *ChemId = new char[COMPLEX_NAME_SZ];
	ChemId[0] = '\0';

	//char Trash[100];

	int 	NPoints = 0,
			NEdges = 0;

	///		Allocate chemical
	//Chemical = new Complex();

	ReadSDF_MolName(ChemId, SDFfile);

	Chemical->SetCName(ChemId);

	if (feof(SDFfile)) {
		delete Chemical;
		return SDF_READ_ERR_EOF;
	}

	if (strlen(ChemId) < 1)
	{
		delete Chemical;
		return SDF_EOF;
	}

	//cout << Chemical->GetCName() << endl;

	//cout << ChemId << endl;

	ReadSDF_Header(NPoints, NEdges, SDFfile);
	if (feof(SDFfile)) {
		delete Chemical;
		return SDF_READ_ERR_EOF;
	}

	if (NPoints == 0 || NEdges == 0)
	{
		//cout << "Wrong number of points or edges in " << Chemical->GetCName() << endl;
		delete Chemical;
		ReadSDF_NextMol(SDFfile);
		if (feof(SDFfile)) {
			return SDF_READ_ERR_EOF;
		}
		return SDF_READ_ERR;
	}

	//cout << "NPoints: " << NPoints << " - NEdges: " << NEdges << endl;

	///		Read molecule's points

	for (int i = 0; i < NPoints; i++)
	{

		double Coordinates[DIM_SMA] = {0., 0., 0.};

		char AtomType[3] = {'\0', '\0', '\0'};

		///		If a point fails to get read, skip rest of molecule and return
		if (!ReadSDF_Point(Coordinates, AtomType, SDFfile))
		{
			//cout << "Wrong point in " << Chemical->GetCName() << endl;
			delete Chemical;
			ReadSDF_NextMol(SDFfile);
			if (feof(SDFfile)) {
				return SDF_READ_ERR_EOF;
			}
			return SDF_READ_ERR;
		}

		///		Store point after successful read
		Point *NewPoint = new Point();

		NewPoint->SetPointId(i+1);
		NewPoint->SetCoordinate(X_COORD, Coordinates[X_COORD]);
		NewPoint->SetCoordinate(Y_COORD, Coordinates[Y_COORD]);
		NewPoint->SetCoordinate(Z_COORD, Coordinates[Z_COORD]);

		NewPoint->SetType(atom_point);

		double Radius = AtomRadius(AtomType);

		NewPoint->SetRadius(Radius);

		Chemical->InsertPointByValue(NewPoint);

	}

	//cout << "Points read" << endl;

	///		Read molecule's edges
	for (int i = 0; i < NEdges; i++)
	{
		int Edges[2] = {-1, -1};

		///		If a point fails to get read, skip rest of molecule and return
		if (!ReadSDF_Edge(Edges, SDFfile, NPoints))
		{
			//cout << "Wrong edge in " << Chemical->GetCName() << endl;
			delete Chemical;
			ReadSDF_NextMol(SDFfile);
			if (feof(SDFfile)) {
				return SDF_READ_ERR_EOF;
			}
			return SDF_READ_ERR;
		}

		//cout << Edges[0] << " " << Edges[1] << endl;

		Point 	*PSource = Chemical->GetPointById(Edges[0]),
				*PTarget = Chemical->GetPointById(Edges[1]);

		Chemical->InsertEdgeByValue(PSource, PTarget);

	}

	//cout << "Edges read" << endl;

	//Chemical->PrintComplex();

	///		Read the rest of characters until next molecule is found and finish

	bool Last = ReadSDF_NextMol(SDFfile);

	if (Last)
	{
		return SDF_EOF;
	}

	if (feof(SDFfile))
	{
		return SDF_EOF;
	}

	return SDF_OK;

}

/*
 *
 */

bool ReadSDF (	vector <Complex *> *Chemicals,
				char *Filename)
{

	FILE *f = fopen(Filename, "r");

	///		Read one by one the molecules until the end of the file
	//		TODO



	bool EndFile = false;

	while (!EndFile)
	{
		Complex *Chemical = new Complex();
		//Complex *Chemical = NULL;
		int ROut = ReadSDF_Single(Chemical, f);

		EndFile = ( (ROut == SDF_READ_ERR_EOF) || (ROut == SDF_EOF) );

		if (ROut != SDF_READ_ERR_EOF && ROut != SDF_READ_ERR)
		{
			//cout << Chemical->GetCName() << endl;
			Chemicals->push_back(Chemical);
		}
		/*else
		{
			delete Chemical;
		}*/

	}

	//cout << "Finshed reading" << endl;

	fclose(f);

	return true;

}

/*****		LOOSELY FORMATTED SDF FILES 	*****/

bool ReadSDFL_Single(	Complex *Chemical,
						FILE *SDFfile)
{


	///		Read the header

	char 	*MolName,
			Discard[COMPLEX_NAME_SZ];

	int NPoints, NEdges;

	MolName = new char[COMPLEX_NAME_SZ];

	int HeadSize = fscanf(SDFfile, "%s\n%s\n  %d %d 0     0  0  0  0  0  0999 V2000\n", MolName, Discard, &NPoints, &NEdges);

	///		Empty header ==> End of file
	if (HeadSize <= 0)
	{
		return false;
	}

	//Chemical = new Complex();

	Chemical->SetCName(MolName);

	///		Read the points

	for (int i = 0; i < NPoints; i++)
	{
		double Coordinates[DIM_SMA];

		char AtomType[3] = {' ', ' ', '\0'};

		HeadSize = fscanf(SDFfile, "%lf  %lf  %lf  %s   0  0  0  0  0  0  0  0  0  0  0  0\n", &(Coordinates[X_COORD]), &(Coordinates[Y_COORD]), &(Coordinates[Z_COORD]), AtomType );

		//cout << "Element: " << AtomType << " - Point: " << Coordinates[X_COORD] << " " << Coordinates[Y_COORD] << " " << Coordinates[Z_COORD] <<  endl;

		double Radius = AtomRadius(AtomType);

		Point *P_i = new Point(i, Coordinates, Radius, atom_point);

		Chemical->InsertPointByValue(P_i);

	}

	///		Read the edges

	for (int i = 0; i < NEdges; i++)
	{
		int Source, Target;
		int BondType;

		HeadSize = fscanf(SDFfile, "%d  %d %d 0 0 0 0\n", &Source, &Target, &BondType);

		Point 	*PSource = Chemical->GetPointById(Source),
				*PTarget = Chemical->GetPointById(Target);

		//Edge *E_i = new Edge(PSource, PTarget, 0.);

		Chemical->InsertEdgeByValue(PSource, PTarget);

	}

	///		Read the end of the molecule
	HeadSize = fscanf(SDFfile, "M  END\n$$$$\n");

	return true;

}

bool ReadSDFL (	vector <Complex *> *Chemicals,
				char *Filename)
{

	FILE *f = fopen(Filename, "r");

	//cout << "File open: " << endl;

	bool Go = true;

	//Complex *Chemical;
	//ReadSDFL_Single(Chemical, f);

	while (Go)
	{
		Complex *Chemical = new Complex();

		Go = ReadSDFL_Single(Chemical, f);

		if (Go)
		{
			//cout << "Pushing" << endl;

			//Chemical->PrintComplex();

			Chemicals->push_back(Chemical);
		}
		else
		{
			delete Chemical;
		}

	}

	//cout << "Number of chemicals " << Chemicals->size() << endl;

	fclose(f);

	return true;


}






// COMPLEX_IO.CC
