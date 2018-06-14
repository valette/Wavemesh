/*=========================================================================

Program:   viewer : a simple mesh visualization example
Module:    vtkSurface
Language:  C++
Date:      2003/07
Author:   Sebastien Valette

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

// .NAME meshviewer 
// .SECTION Description

#include "vtkSurface.h"
#include "RenderWindow.h"
#include "vtkOBJExporter.h"

/// A simple example of mesh visualization
/// Usage : viewer mesh1 mesh2 .... meshN
/// where mesh1, mesh2, ... meshN are all mesh files. They will be displayed in the same window

int main( int argc, char *argv[] )
{
	RenderWindow *Window=RenderWindow::New();

	vtkSurface *Mesh;
	for (int i=0;i<argc-1;i++)
	{
		// Load the mesh and create the vtkSurface data structure
		Mesh=vtkSurface::New();
		cout <<"load : "<<argv[i+1]<<endl;
		Mesh->CreateFromFile(argv[i+1]);

		// prints to standard output the mesh caracteristics
		Mesh->DisplayMeshProperties();

		// Create a renderwindow
		if (i==0)
			Window->SetInputData(Mesh);
		else
			Window->AddPolyData(Mesh);

		Window->SetWindowName(argv[i+1]);
		Mesh->Delete();
	}

	// start display and interaction
	Window->Render();
	Window->Interact();
	
	vtkOBJExporter *Export=vtkOBJExporter::New();
	Export->SetInput(Window->GetvtkRenderWindow());
	Export->SetFilePrefix("export");
	Export->Write();
	
	
	Window->Delete();
	return (0);
}
