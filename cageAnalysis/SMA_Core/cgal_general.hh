//	Software for single molecule analysis
//
//	cgal_general.hh
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : February 1st 2016


//#include "../SMA_Core/geometry.hh"
//#include "../SMA_Core/molecule_info.hh"

// CGAL Libraries
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>

#include <CGAL/Dynamic_matrix.h>
#include <CGAL/Linear_algebraCd.h>




////// CGAL typedefs

typedef CGAL::Exact_predicates_inexact_constructions_kernel					CG_Kernel;

// Types for unweighted alpha shapes

//typedef CGAL::Simple_cartesian<double>               						CG_Kernel;
/*typedef CGAL::Alpha_shape_vertex_base_3<CG_Kernel>               			Vb;
typedef CGAL::Alpha_shape_cell_base_3<CG_Kernel>                 			Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>      					Tds;
typedef CGAL::Delaunay_triangulation_3<CG_Kernel,Tds,CGAL::Fast_location>  	Delaunay_triangulation;
typedef CGAL::Alpha_shape_3<Delaunay_triangulation>                    		Alpha_shape_3;*/



// Types for weighted alpha shapes
typedef CGAL::Regular_triangulation_euclidean_traits_3<CG_Kernel> 			Gt;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>         						Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           						Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> 						Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       						Triangulation_3;

typedef Gt::Weighted_point                  								Weighted_point;
typedef Gt::Bare_point                      								Bare_point;

typedef CGAL::Alpha_shape_3<Triangulation_3>        						Alpha_shape_3;

// Alpha shape (general) types
typedef Alpha_shape_3::Alpha_iterator               						Alpha_iterator;
typedef Alpha_shape_3::NT                           						NT;

typedef Alpha_shape_3::Cell_handle 											Cell_handle;
typedef Alpha_shape_3::Edge 												CG_edge;
typedef Alpha_shape_3::Facet 												Facet;
typedef Alpha_shape_3::Vertex_handle										Vertex_handle;

typedef Alpha_shape_3::Classification_type									Classification_type;

typedef CGAL::Polyhedron_3<CG_Kernel> 										Polyhedron;
typedef Polyhedron::Halfedge_around_facet_circulator 						Halfedge_facet_circulator;
typedef Polyhedron::Facet_iterator 											Facet_iterator;

typedef CG_Kernel::Point_3                                  				CG_point;

// Type for cross product
typedef CGAL::Vector_3<CG_Kernel>											CG_vector_3;



//typedef CGAL::Linear_algebraHd<CG_Kernel>									Linear_algebra;


// Constants

#define CELL_N_V	4
#define FACET_N_V	3
#define EDGE_N_V	2

//void RetrievePointFromFacet (	Facet TargetFacet,
//								int FPIndex,
//								double *Point);



