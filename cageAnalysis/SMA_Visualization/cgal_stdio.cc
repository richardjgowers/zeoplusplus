//	Software for single molecule analysis
//
//	cgal_stdio.cc
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : January 25th 2016

#include "cgal_stdio.hh"

void PrintClassType (Classification_type CType)
{
	switch (CType)
	{
		case Alpha_shape_3::INTERIOR:
			cout << "INTERIOR";
			break;
		case Alpha_shape_3::EXTERIOR:
			cout << "EXTERIOR";
			break;
		case Alpha_shape_3::SINGULAR:
			cout << "SINGULAR";
			break;
		case Alpha_shape_3::REGULAR:
			cout << "REGULAR";
			break;
		default:
			cout << "ALL";
			break;
	}


}

void PrintAlphaCells (Alpha_shape_3 &MoleculeAlphaShape, Classification_type CType, double Alpha)
{


	list<Cell_handle> CellList;
	list<Cell_handle>::const_iterator CIt;

	// Get the list of cells of the alpha shape
	MoleculeAlphaShape.get_alpha_shape_cells(back_inserter(CellList), CType, Alpha);


	cout << "--- PRINTING CELLS: ";
	PrintClassType(CType);
	cout << " --- " << endl;

	if (CellList.size() == 0)
		cout << "(none)" << endl;

	// Go though the list of cells, printing vertices belonging to each cell
	for (CIt = CellList.begin(); CIt != CellList.end(); CIt++)
	{

		Cell_handle Cell = *CIt;

		for (int i = 0; i < CELL_N_V; i++)
		{
			cout << "(" << Cell->vertex(i)->point().x() << " ," << Cell->vertex(i)->point().y() << " ,"<< Cell->vertex(i)->point().z() << "); "; 
		}

		cout << endl;

	}

}

void PrintAlphaFacets (Alpha_shape_3 &MoleculeAlphaShape, Classification_type CType, double Alpha)
{


	list<Facet> FacetList;
	list<Facet>::const_iterator FacetListIt;
	
	//Facet_iterator FIt;

	// Get the list of facets of the alpha shape
	MoleculeAlphaShape.get_alpha_shape_facets(back_inserter(FacetList), CType, Alpha);

	cout << "--- PRINTING FACETS: ";
	PrintClassType(CType);
	cout << " --- " << endl;

	if (FacetList.size() == 0)
		cout << "(none)" << endl;

	// Go though the list of facets, printing vertices belonging to each cell
	for (FacetListIt = FacetList.begin(); FacetListIt != FacetList.end(); FacetListIt++)
	{

		Facet TargetF = *FacetListIt;
		int OppositeVertId = TargetF.second;

		// Run around the cell, checking that the index for the vertices remain to the facet
		for (int i = 0; i < CELL_N_V; i++)
		{
			if (i != OppositeVertId)
				cout << "(" << TargetF.first->vertex(i)->point().x() << ", " << TargetF.first->vertex(i)->point().y() << ", " << TargetF.first->vertex(i)->point().z() << " );";
		}
		cout << endl;
		
	 }

	cout << "There are " << FacetList.size() << " facets" << endl;

}

void PrintAlphaEdges (Alpha_shape_3 &MoleculeAlphaShape, Classification_type CType, double Alpha)
{


	list<CG_edge> EdgeList;
	list<CG_edge>::const_iterator EIt;

	// Get the list of edges of the alpha shape
	MoleculeAlphaShape.get_alpha_shape_edges(back_inserter(EdgeList), CType, Alpha);

	// Go though the list of edges, printing vertices belonging to each cell

	cout << "--- PRINTING EDGES: ";
	PrintClassType(CType);
	cout << " --- " << endl;

	if (EdgeList.size() == 0)
		cout << "(none)" << endl;

	for (EIt = EdgeList.begin(); EIt != EdgeList.end(); ++EIt)
	{

		CG_edge TargetEdge = *EIt;

		int Vert1 = TargetEdge.second, Vert2 = TargetEdge.third;

		for (int i = 0; i < FACET_N_V; i++)
			for (int j = i+1; j < FACET_N_V; j++)
			{
				if ( (Vert1 == i && Vert2 == j) || (Vert1 == j && Vert2 == i) )
				{
					cout << "(" << TargetEdge.first->vertex(i)->point().x() << " ," << TargetEdge.first->vertex(i)->point().y() <<  " ," << TargetEdge.first->vertex(i)->point().z() << "); "; 
					cout << "(" << TargetEdge.first->vertex(j)->point().x() << " ," << TargetEdge.first->vertex(j)->point().y() <<  " ," << TargetEdge.first->vertex(j)->point().z() << "); ";
					cout << endl; 
					break;
				}
			}

		

	}

}


void PrintAlphaVertices (Alpha_shape_3 &MoleculeAlphaShape, Classification_type CType, double Alpha)
{


	list<Vertex_handle> VertexList;
	list<Vertex_handle>::const_iterator VertexListIt;
	
	// Get the list of vertices of the alpha shape
	MoleculeAlphaShape.get_alpha_shape_vertices(back_inserter(VertexList), CType, Alpha);

	

	cout << "--- PRINTING VERTICES: ";
	PrintClassType(CType);
	cout << " --- " << endl;

	if (VertexList.size() == 0)
		cout << "(none)" << endl;

	// Go though the list of facets, printing vertices belonging to each cell

	for (VertexListIt = VertexList.begin(); VertexListIt != VertexList.end(); VertexListIt++)
	{

		Vertex_handle TargetF = *VertexListIt;

		cout << "(" << TargetF->point().x() << ", " << TargetF->point().y() << ", " << TargetF->point().z() << " );";
		cout << endl;
		
	 }

}
