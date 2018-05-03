# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Tell make that these are phony targets
#.PHONY: all clean
.PHONY:

include config.mk

# Object names
objs=		geometry.o \
			point.o \
			edge.o \
			complex.o \
			window.o \
			molecule_info.o \
			simplex.o \
			cell_grid.o \
			dijkstra_path.o \
			dijkstra_aux.o \
			internal_cells.o \
			complex_io.o \
			volume.o \
			window_comp.o \
			node_characteristic.o \
			voro_interface.o \
			window_PER.o \
			alpha_shape.o \
			clustering.o \
			preprocess.o \
			chem_info.o \
			internal_external_class.o \
			surface_area.o \
			SMA.o
			
			#voroedge_classificator.o \
			windows_aux.o \
			cgal_stdio.o \
			cgal_general.o \
			
			
			

# Paths to source files
objs_path=	SMA_Core/geometry.o \
			SMA_Core/point.o \
			SMA_Core/edge.o \
			SMA_Core/complex.o \
			SMA_Core/window.o \
			SMA_Core/molecule_info.o \
			SMA_Core/simplex.o \
			SMA_Preprocess/alpha_shape.o \
			SMA_Preprocess/window_comp.o \
			SMA_Preprocess/node_characteristic.o \
			SMA_Preprocess/complex_io.o \
			SMA_Preprocess/chem_info.o \
			SMA_Preprocess/preprocess.o \
			SMA_Preprocess/internal_cells.o \
			SMA_Preprocess/clustering.o \
			SMA_Preprocess/cell_grid.o \
			SMA_Preprocess/voro_interface.o \
			SMA_Preprocess/window_PER.o \
			SMA_Postprocess/dijkstra_path.o \
			SMA_Postprocess/volume.o \
			SMA_Postprocess/dijkstra_aux.o \
			SMA_Postprocess/internal_external_class.o \
			SMA_Postprocess/surface_area.o \
			SMA.o
			#SMA_Preprocess/voroedge_classificator.o \
					SMA_Analysis/windows_aux.o \
						SMA_Interface/cgal_stdio.o \
						SMA_GeoGraphComp/cgal_general.o \
						
			

# Build all of the executable files
all: single_mol_analysis1

# Build the libraries
single_mol_analysis1: $(objs_path) single_mol_analysis.cc
	$(CXX) $(CFLAGS) $(objs) $(E_LIB) single_mol_analysis.cc -o single_mol_analysis

%.o: %.cc
	$(CXX) $(CFLAGS) -c $<

# Clean up the executable files
clean:
	rm *.o *~ single_mol_analysis