/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

// .NAME vtkRandomTriangulation - Creates meshes with random triangulation

#ifndef __VTKRANDOMTRIANGULATION
#define __VTKRANDOMTRIANGULATION

#include "vtkSurface.h"

class VTK_EXPORT vtkRandomTriangulation : public vtkObject
{
public:

	/// Description
	/// Builds a random triangulation
	/// the entry Type gives the distribution :
	/// 0 :uniform plane 
	/// 1 : non uniform plane 
	/// 2: plane made of 4 different regions with different densities
	/// 3 : half pipe
	/// 4 : half pipe cut on one corner
	static vtkRandomTriangulation *New();

	static vtkSurface * BuildRandomTriangulation (int NumberOfPoints, int Type);


protected:
  vtkRandomTriangulation();
  ~vtkRandomTriangulation();
};

#endif
