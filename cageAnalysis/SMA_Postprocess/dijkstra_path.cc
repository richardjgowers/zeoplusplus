// Voro++, a 3D cell-based Voronoi library
//
// dijkstra_path.cc
//
// Routines for single molecule analysis
//
// Author   : Ismael Gomez Garcia
// Email    : 
// Date     : November 16th 2015

#include "dijkstra_path.hh"


/*********************************************************/
/************* POINT CLASSIFICATION FUNCTIONS ************/
/*********************************************************/

/* ClasifyPoins
 *
 * Checks the given complex (with the information of one molecule trapped in a cage).
 * Updates the information about the types of points deciding if they are internal or external.
 */
void ClassifyPoints(	Complex *MoleculeC)
{

	// TODO: Apply Delaunay triangulation to the complex and calculate internal and external edges

	// Convex closure of the object is "inside"
	for (int i=0; i < MoleculeC->GetNPoints(); i++)
		if (MoleculeC->GetPointList()[i]->GetType() == undefined_point)
			MoleculeC->GetPointList()[i]->SetType(internal_point);
}


bool RestrictionCriterion ( 	Complex *MoleculeC, 
								Point *TargetPoint)
{

	Edge **EdgeList;
	int Size = 0;

	double ChannelRadius, OppositeRadius, TargetRadius = TargetPoint->GetRadius();

	TargetPoint->PrintPoint();


	// Get the edges that connect the point with a boundary point
	EdgeList = MoleculeC->GetEdgesByPointOppositeRestricted(TargetPoint, boundary_point, Size);

	// If there were no edges like that, the point is internal
	if (Size == 0)
		return false;

	// For each edge, check if the restriction criterion is satisfied; if one edge satisfies it, return true
	for (int i = 0; i < Size; i++)
	{
		ChannelRadius = EdgeList[i]->GetChannelRadius();
		OppositeRadius = EdgeList[i]->GetOpposite(TargetPoint)->GetRadius();

		EdgeList[i]->PrintEdge();

		if (OppositeRadius > ChannelRadius && ChannelRadius > TargetRadius)
			return true;
		
	}
	
	return false;


}

/* ClassifyPointsRestrictionCriterion
 *
 * Checks the given complex (with the information of one molecule trapped in a cage).
 * Updates the information about the types of points deciding if they are internal or external.
 * Classification criterion: Checks restriction through edge, deciding if the point is internal or external in depends on that. 
 * - Check undefined points connected to boundary points.
 * - If there is a restriction due to the edge, it considers the point to be internal. 
 * - If the edge does not strongly restrict (i.e., the cavity of the point restricts more), the point is considered external.
 */

void ClassifyPointsRestrictionCriterion (Complex *MoleculeC)
{

	// Expected to have 10% points of non-boundary points (maybe more)
	int ExpectedUndefPoints = MoleculeC->GetNPoints()/10;
	
	vector <Point *> UndefPointsVector(ExpectedUndefPoints);

	// Create the list of undefined points
	for (int i = 0; i < MoleculeC->GetNPoints(); i++)
	{
		if (MoleculeC->GetPointList()[i]->GetType() == undefined_point)
		{
			if ( RestrictionCriterion(MoleculeC, MoleculeC->GetPointList()[i]) )
			{
				MoleculeC->GetPointList()[i]->SetType(external_point);
			}
			else
			{
				MoleculeC->GetPointList()[i]->SetType(internal_point);
			}
		}
	}


}


/********** FUNCTIONS FOR CALCULATION OF CELL ACCESIBLE SIZE **********/


/* GetChannelSize
 *
 * Function that calculates the size of a channel formed by two nodes and the edge that connects them.
 * Returns the radius of the maximum sphere capable of crossing that channel (i.e., the minimum radius of the channel).
 */
double GetChannelSize(	Complex *MoleculeC, 
			int SourcePointId, 
			int TargetPointId)
{

	double SourceSize, ChannelSize, TargetSize;

	Edge *Channel = MoleculeC->GetEdgeByIds(SourcePointId, TargetPointId);

	ChannelSize = Channel->GetChannelRadius();
	SourceSize = Channel->GetPSource()->GetRadius();
	TargetSize = Channel->GetPTarget()->GetRadius();
	
	// Return the minimum of these three values
	if (SourceSize <= ChannelSize)
	{
		if (SourceSize <= TargetSize)
			return SourceSize;
		else 
			return TargetSize;
	}
	else if (TargetSize <= ChannelSize)
		return TargetSize;

	return ChannelSize;

}


/* VisitNode
 *
 * Visits a node as part of the Modified Dijkstra Algorithm, marking it as visited, adding its neighbors to the queue (if necessary)
 * and computing weights.
 * 
 */
void VisitNode(	int PointId,
		Complex *MoleculeC,
		PriorityQueue **PQList,
		VisitedList *VList,
		WeightList *WList)
{

	int i;
	double Weight;

	int *NeighborIdList = new int[MoleculeC->GetNEdges()];
	int NNeigh;

	PointType Type;
	Point *TargetPoint = MoleculeC->GetPointById(PointId);

	//printf("Visiting node %d\n", PointId);

	// Insert node in visited list
	if (CheckIfVisited(PointId, VList) != 0)
		InsertInVisited(PointId, VList);
	else
		return; // TODO: Err msg
	
	// Calculate node neighbor list
	NNeigh = MoleculeC->CalculateNeighborIdList(TargetPoint, NeighborIdList);
	
	// Update neighbor weights and insert unvisited neighbors in list
	for (i = 0; i < NNeigh; i++)
	{
		// Get the size of the channel of interest and compare with current weight
		Weight = GetChannelSize(MoleculeC, PointId, NeighborIdList[i]);
		if (WList->ElementList[PointId]->Weight < Weight)
			Weight = WList->ElementList[PointId]->Weight;

		UpdateWeight(NeighborIdList[i], Weight, WList);

		// Check if node is external
		Type = MoleculeC->GetPointById(NeighborIdList[i])->GetType();
		// Insert unvisited nodes in queue (as long as they're internal)
		if (CheckIfVisited(NeighborIdList[i], VList) != 0 && Type == internal_point)
			*PQList = AddElementToPriorityQueue(NeighborIdList[i], *PQList);
			//InsertInQueue(NeighborIdList[i], PriorityQueue);

	}	
	
	//TODO Destroy(NeighborIdList);
}


/* CellAccessibleSize
 *
 * Checks a particular cell inside the given complex (with the information of one molecule trapped in a cage).
 * Returns the radius of the maximum sphere capable of escaping from that cell to the outside.
 */
double CellAccessibleSize (	Complex *MoleculeC,
							Point *TargetCell)
{
	int i;
	int PointId;

	// To calculate maximum
	int AuxId;
	PointType AuxType;

	double MaxWeight = 0.0;

	//Queue *PriorityQueue = NewQueue();
	VisitedList *VList = NewVisitedList();
	WeightList *WList = NewWeightList();
	PriorityQueue *PQList = NULL;

	// Error control
	if (TargetCell->GetType() != internal_point)
		return FLOAT_ERR;

	/***** Modified Dijkstra Algorithm *****/

	// Create weight list (all weights are 0 except the initial one)
	for (i = 0; i < MoleculeC->GetNPoints(); i++)
	{
		PointId = MoleculeC->GetPointList()[i]->GetPointId(); 
		if ( PointId != TargetCell->GetPointId())
			AddWeight(PointId, 0., WList);
		else
			AddWeight(PointId, TargetCell->GetRadius(), WList );
	}

	// Insert external and boundary nodes in visited list (so they're never visited)
	for (i = 0; i < MoleculeC->GetNPoints(); i++)
	{
		AuxType = MoleculeC->GetPointList()[i]->GetType();
		if (AuxType == boundary_point || AuxType == external_point)
			InsertInVisited(MoleculeC->GetPointList()[i]->GetPointId(), VList);
	}		
			
	// Visit initial node (insert in visited list, update neighbors)
	VisitNode(TargetCell->GetPointId(), MoleculeC, &PQList, VList, WList);

	// While queue is not empty, keep visiting elements
	while (PQList != NULL)
	{
		PointId = PointWithMaxWeight(WList, VList);
		PQList = ExtractFromQueue(PointId, PQList);
		VisitNode(PointId, MoleculeC, &PQList, VList, WList);
	}

	// TODO Change so it takes only external nodes
	for (i = 0; i < WList->Length; i++)
	{
		AuxId = WList->ElementList[i]->IdPoint; 
		AuxType = MoleculeC->GetPointById(AuxId)->GetType();
		if (WList->ElementList[i]->Weight > MaxWeight)
			if (AuxType == external_point || AuxType == boundary_point)
				MaxWeight = WList->ElementList[i]->Weight;
	}
	
	/*****/
	// Finish algorithm


	//DestroyQueue(PriorityQueue);
	DestroyVisitedList(VList);
	DestroyWeightList(WList);

	return MaxWeight;
}


/*
 *
 */
double MoleculeAccessibleSize(Complex *MolVoronoi)
{
	double MaxSize = 0.0;

	for (int i = 0; i < MolVoronoi->GetNPoints(); i++)
	{
		Point *PAux = MolVoronoi->GetPointList()[i];

		if (PAux->GetType() == internal_point)
		{
			double NodeSize = CellAccessibleSize(MolVoronoi, PAux);

			if (NodeSize > MaxSize)
			{
				MaxSize = NodeSize;
			}
		}

	}

	return MaxSize;

}

/********** FUNCTIONS FOR CALCULATION OF CELL PATHS **********/



/* CellEscapePath
 *
 * Checks a particular cell inside the given complex (with the information of one molecule trapped in a cage).
 * Returns the escape path for the biggest molecule trapped in that cell.
 */
Point **CellEscapePath (	Complex *MoleculeC,
							Point *TargetCell)
{

	std::vector<Complex *> PathList;

	return NULL;
	
	//return PathList;
}

/* CellProbeEscapePath
 *
 * Checks a particular cell inside the given complex (with the information of one molecule trapped in a cage).
 * Returns the escape path for a probe of given size trapped in that cell (Size_probe < Size_cell).
 */
Point **CellProbeEscapePath (	Complex *MoleculeC,
								Point *TargetCell)
{

	Point **Path = NULL;
	
	return Path;
}





