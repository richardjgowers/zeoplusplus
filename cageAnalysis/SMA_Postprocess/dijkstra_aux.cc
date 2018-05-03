// Voro++, a 3D cell-based Voronoi library
//
// dikjstra_aux.c
//
// Routine implementation for modified dikjstra analysis
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 20th 2015

#include "dijkstra_aux.hh"

/***** WEIGHT LIST MANAGEMENT FUNCTIONS *****/

WListElement *NewWListElement( 	int IdPoint, 
				float Weight)
{

	WListElement *Element = (WListElement *) malloc (sizeof(WListElement));

	Element->IdPoint = IdPoint;
	Element->Weight = Weight;

	return Element;
}

WeightList *NewWeightList()
{
	WeightList *WList = (WeightList *) malloc(sizeof(WeightList));

	// Allocate array of elements
	WList->ElementList = (WListElement **) malloc (WLIST_SIZE*sizeof(WListElement *));

	WList->Length = 0;

	return WList;

}

void AddWeight (	int IdPoint,
			float Weight,
			WeightList *WList)
{

	int LastElementIndex = WList->Length;

	// Check errors
	if (LastElementIndex >= WLIST_SIZE)
	{
		printf("Warning: Maximum size of weightlist has been surpassed\n");
		return;
	}


	WListElement *NewElement = NewWListElement(IdPoint, Weight);

	// Store new weight entry in table
	WList->ElementList[LastElementIndex] = NewElement;

	// Update length
	WList->Length++;

}

void UpdateWeight(	int IdPoint, 
			float NewWeight,
			WeightList *WList)
{

	int i, flag = 0;

	// Look for the point being updated, and update weight if new value is bigger than previous one
	for (i = 0; i<WList->Length; i++)
		if (WList->ElementList[i]->IdPoint == IdPoint)
		{
			flag = 1;
			if (WList->ElementList[i]->Weight < NewWeight)
			{

				//printf("New weight of point %d is %f\n", IdPoint, NewWeight);
				WList->ElementList[i]->Weight = NewWeight;
				return;
			}
		}

	// Warning if no element was found
	if (flag == 0)
		printf("Warning: UpdateWeight couldn't find point %d\n", IdPoint);

}

float MaximumWeight(	WeightList *WList)
{

	int i;
	float Max = 0.0;


	for (i = 0; i < WList->Length; i++)
		if (Max < WList->ElementList[i]->Weight)
				Max = WList->ElementList[i]->Weight;

	return Max;
}

int PointWithMaxWeight(	WeightList *WList,
			VisitedList *VList)
{
	int i;
	int IdMax = -1;
	float MaxWeight = 0.0;

	for (i = 0; i < WList->Length; i++)
		if (MaxWeight < WList->ElementList[i]->Weight)
			if(CheckIfVisited(WList->ElementList[i]->IdPoint, VList) != 0)
			{
				MaxWeight = WList->ElementList[i]->Weight;
				IdMax = WList->ElementList[i]->IdPoint;
			}

	return IdMax;


}

void DestroyWeightList(WeightList *WList)
{
	int i;

	// Free elements one by one
	for (i = 0; i < WList->Length; i++)
		free(WList->ElementList[i]);

	// Free element list
	free(WList->ElementList);

	// Free weight list
	free(WList);
}

/***** PRIORITY QUEUE STRUCTURE AND MANAGEMENT FUNCTIONS *****/

PriorityQueue *AddElementToPriorityQueue (	int PointId,
						PriorityQueue *PQList)
{

	// Auxiliary pointer to the first element of the list
	PriorityQueue *PQAux = PQList;
	PriorityQueue *PQNew = NULL;

//printf("\tAdding %d to PQList\n", PointId);

	// Check if the queue is empty
	if (PQList == NULL)
	{
		PQNew = (PriorityQueue *) calloc(sizeof(PriorityQueue), 1);
		PQNew->PointId = PointId;
		return PQNew;
	}

	// Check if a node with the same id is already in the queue
	while (PQAux != NULL)
		if (PQAux->PointId == PointId)
		{
//printf("\t\tNot added -> Already in queue\n");
			return PQList;
		}
		else
			PQAux = PQAux->Next;

	// Element to be inserted
	PQNew = (PriorityQueue *) calloc(sizeof(PriorityQueue), 1);
	PQNew->PointId = PointId;

	// Instert element at the beggining
	PQNew->Next = PQList;

	return PQNew;

}

PriorityQueue *ExtractFromQueue(	int TargetPoint,
					PriorityQueue *PQList)
{
	PriorityQueue *PQAux = PQList, *PQAuxPrev = NULL;

	while (PQAux->PointId != TargetPoint)
	{
		PQAuxPrev = PQAux;		
		PQAux = PQAux->Next;
		if (PQAux == NULL)
		{
			printf("Warning: ExtractFromQueue couldn't find point %d\n", TargetPoint);
			return PQList;
		}
	}

	// If the element is not the first
	if (PQAuxPrev != NULL)
	{
		PQAuxPrev->Next = PQAux->Next;
		free(PQAux);
		return PQList;
	}

	// If the element is exactly the first
	if (PQAux != NULL)
	{
		PQAuxPrev = PQAux->Next;
		free(PQAux);
		return PQAuxPrev;
	}

	// If the element is not in the list
	return PQList;	

	


}	

// TODO: Destroy priority queue

/***** VISITED LIST RELATED FUNCTIONS *****/

VisitedList *NewVisitedList ()
{
	VisitedList *VList = (VisitedList *) malloc (sizeof(VisitedList));

	VList->Length = 0;

	return VList;
}

void InsertInVisited(	int Element,
			VisitedList *VList)
{

	// TODO: Error control

	VList->VisitedElements[VList->Length] = Element;
	VList->Length++;

}

int CheckIfVisited(	int Element,
			VisitedList *VList)
{

	int i;

	for (i = 0; i < VList->Length; i++)
		if (VList->VisitedElements[i] == Element)
			return 0;

	return -1;

}

void DestroyVisitedList(VisitedList *VList)
{
	if (VList != NULL)
		free(VList);
}


