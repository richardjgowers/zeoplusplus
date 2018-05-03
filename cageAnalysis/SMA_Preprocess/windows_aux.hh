//	Software for single molecule analysis
//
//	windows_aux.hh
//
// 	Author   : Ismael Gomez Garcia
// 	Email    : 
//	Date     : February 5th 2016

#include "../SMA_Core/complex.hh"

bool IsWindow(Complex *Window);

//double FacetArea (	double Point1[DIM_SMA],
//					double Point2[DIM_SMA],
//					double Point3[DIM_SMA]);

Complex *CompressToSphere (	Complex *TargetComplex);
