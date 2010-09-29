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

/// A simple example of mesh visualization
/// Usage : viewer mesh1 mesh2 .... meshN
/// where mesh1, mesh2, ... meshN are all mesh files. They will be displayed in different windows,
/// but the windows will have the same viewports

int main( int argc, char *argv[] )
{
	RenderWindow **Window=new RenderWindow*[argc-1];
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
		Window[i]=RenderWindow::New();
		Window[i]->SetInput(Mesh);

		// synchronize the viewport with the first created window		
		if (i>0)
			Window[i]->AttachToRenderWindow(Window[0]);

		Window[i]->Render();
		Window[i]->SetWindowName(argv[i+1]);
		Mesh->Delete();
	}

	// start interaction
	Window[0]->Interact();
	
	// Delete objects before exit
	for (int i=0;i<argc-1;i++)
		Window[i]->Delete();
	delete [] Window;

	return (0);
}
