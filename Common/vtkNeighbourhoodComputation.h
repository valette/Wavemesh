/***************************************************************************
vtkNeighbourhoodComputation.h  -  description
-------------------
begin                : January 4 2006
copyright            : (C) 2006 by Sebastien Valette
email                : 
***************************************************************************/
/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

#ifndef _VTKNEIGHBOURHOODCOMPUTATION_H_
#define _VTKNEIGHBOURHOODCOMPUTATION_H_

#include "vtkCommand.h"
#include "vtkObjectFactory.h"
#include "vtkSurface.h"

/// A Class to compute Neighbourhoods on meshes
/// The class is designed for efficient and possibly multithreaded computation
class VTK_EXPORT vtkNeighbourhoodComputation : public vtkObject
{
public:
	/// the public constructor
	static vtkNeighbourhoodComputation* New();

	/// Method to initialize the class for a specific Input
	void SetInput(vtkSurface *Mesh);

	/// Compute the NRing around the Cell.
	void ComputeNRingCells(vtkIdType Cell,int RingSize,vtkIdList *FList);

	/// Compute the Cells around the input cell within the input Distance.
	void ComputeDistanceRingCells (vtkIdType Cell,double Distance, vtkIdList *FList);

	// Defines the type of the origin cells (0=faces 1=Vertices)
	void SetCellType (int Type)
	{ this->CellType=Type;};

	/// Returns the Input mesh
	vtkSurface *GetInput(){return (this->Input);};

protected:

	vtkNeighbourhoodComputation();
	~vtkNeighbourhoodComputation();

private:

	// Type of the origin cells (0=faces 1=Vertices)
	int CellType;
	
	// this parameter stores the number of times the methods 
	// ComputeNRingCells() or ComputeDistanceRingCells() were called.
	int Time;
	
	// this method increases the "Time" class member, and check wether there is an overflow.
	// In case of overflow, Time is reset to 0 and the visited fields are also reset.
	// Note : overflow will very unlikely happen , unless you call for Neighborhood computation 
	// a lot of times (more INT_MAX)
	void IncreaseTime();
	
	// Reset the arrays defining whether an item was already visited
	void InitArrays();

	// The input mesh
	vtkSurface *Input;

	// those arrays define which elements have already been visited
	vtkIntArray *VisitedCells;
	vtkIntArray *VisitedEdges;
	vtkIntArray *VisitedVertices;

	// IdLists statically created to speed up the neighborhood computation
	vtkIdList *VList;
	vtkIdList *VList2;
	vtkIdList *EList;
};
#endif
