// Voro++, a 3D cell-based Voronoi library
//
// dikjstra_aux.hh
//
// Routines for modified dikjstra analysis
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 20th 2015


#ifndef _DIKJSTRA_AUX_H
#define _DIKJSTRA_AUX_H

//#include "../SMA_Core/complex.hh"

#include "../SMA_Core/molecule_info.hh"

#include <stdlib.h>
#include <stdio.h>


#define QUEUE_SIZE	1000

/***** VISITED LIST STRUCTURE AND RELATED FUNCTIONS DECLARATION *****/

#define VISITED_SIZE 1000

typedef struct VisitedList_
{
	int VisitedElements[VISITED_SIZE];
	int Length;

} VisitedList;

VisitedList *NewVisitedList ();

void InsertInVisited(	int Element,
			VisitedList *VList);

int CheckIfVisited(	int Element,
			VisitedList *VList);

void DestroyVisitedList(VisitedList *VList);

/***** WEIGHT LIST STRUCTURES AND RELATED FUNCTIONS DECLARATION *****/

#define WLIST_SIZE	10000
#define MAX_WEIGHT	1000.

typedef struct WListElement_
{
	int IdPoint;
	float Weight;

} WListElement;

typedef struct WeightList_
{
	WListElement **ElementList;
	int Length;

} WeightList;


WListElement *NewWListElement( 	int IdPoint, 
				float Weight);

WeightList *NewWeightList();

void AddWeight (	int IdPoint,
			float Weight,
			WeightList *WList);

void UpdateWeight(	int IdPoint, 
			float NewWeight,
			WeightList *WList);

float MaximumWeight(	WeightList *WList);

int PointWithMaxWeight(	WeightList *WList,
			VisitedList *VList);

void DestroyWeightList(	WeightList *WList);

/***** PRIORITY QUEUE STRUCTURE AND MANAGEMENT FUNCTIONS *****/
typedef struct PriorityQueue_
{
	PriorityQueue_ *Next;
	int PointId;

} PriorityQueue;

PriorityQueue *ExtractFromQueue(	int TargetPoint,
					PriorityQueue *PQList);

PriorityQueue *AddElementToPriorityQueue (	int PointId,
						PriorityQueue *PQList);





#endif

