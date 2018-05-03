/*
 * window.hh
 *
 *  Created on: Mar 7, 2016
 *      Author: ismael.gomez
 */

#ifndef WINDOW_HH_
#define WINDOW_HH_

#include "complex.hh"

#include "simplex.hh"

#include <list>


class Window : public Complex
{

public:

	Window(Point **PointList, Edge **EdgeList, int NPoints, int NEdges, double Surface, list <Simplex *> *Simplices);
	Window();

	~Window();

	double GetSurface();
	list <Simplex *> *GetSimplices();

	void InsertSimplex(Simplex *NewSimplex);

	list <Simplex *>::iterator EraseSimplex(list <Simplex *>::iterator Target);

	void SetSurface(double Surface);

	void PrintWindow();

	bool SameWindow(Window *WindowT);

	void SurfaceAreaUpdate();


	Complex *JointSimplices();


private:

	double Surface;

	list <Simplex *> *Simplices;

};


#endif /* WINDOW_HH_ */
