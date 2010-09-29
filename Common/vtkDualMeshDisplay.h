/*=========================================================================

  Program:   vtkDualMeshDisplay
  Module:    vtkSurface
  Language:  C++
  Date:      2008/03
  Author:    Sebastien VALETTE

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  
#ifndef __vtkDualMeshDisplay_h
#define __vtkDualMeshDisplay_h

#include "vtkSurface.h"
#include "RenderWindow.h"

class VTK_EXPORT vtkDualMeshDisplay : public vtkObject
{

public:

	/// The Constructor vtkSurfaceToDual::New();
	static vtkDualMeshDisplay *New();
	vtkTypeMacro(vtkDualMeshDisplay,vtkObject);

	void SetInput(vtkSurface *Input);

	void SetColors(vtkIntArray *Colors);

	vtkSurface* GetOutput();

	RenderWindow *GetRenderWindow();
	
	void Update();

protected:

	void CreateDual();

	RenderWindow *Window;

	vtkIntArray *Colors;
	
	vtkSurface *Input;

	vtkSurface *Dual;

	/// constructor
	vtkDualMeshDisplay();

	/// desctructor
	~vtkDualMeshDisplay();
};
#endif
