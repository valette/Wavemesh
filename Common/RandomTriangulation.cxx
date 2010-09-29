/*=========================================================================

Program:   Random Triangulation generator
Module:    wavemesh.cxx
Language:  C++
Date:      2006/12
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
#include <vtkPLYWriter.h>
#include "vtkSurface.h"
#include "RenderWindow.h"
#include "vtkRandomTriangulation.h"

int main( int argc, char *argv[] )
{
	if (argc<3)
	{
		cout<<"Random Triangulation generator"<<endl;
		cout<<"usage : \"RandomTriangulation numberofvertices type\" "<<endl;
		cout<<"available types :"<<endl;
		cout<<"0 : uniform plane "<<endl;
		cout<<"1 : non uniform plane"<<endl;
		cout<<"2 : plane made of 4 different regions with different densities"<<endl;
		cout<<"3 : half pipe"<<endl;
		cout<<"4 : half pipe cut on one corner"<<endl;
		return (1);
	}
	
	vtkRandomTriangulation *Triangulation=vtkRandomTriangulation::New();
	vtkSurface *Mesh=Triangulation->BuildRandomTriangulation(atoi(argv[1]),atoi(argv[2]));
	
	Mesh->DisplayMeshProperties();
	
	RenderWindow *Window=RenderWindow::New();
	Window->SetInput(Mesh);
	Window->Render();
	Window->Interact();
	
	vtkPLYWriter *Writer=vtkPLYWriter::New();
	Writer->SetInput(Mesh);
	Writer->SetFileName("mesh.ply");
	Writer->Write();
	Writer->Delete();
	
	Mesh->Delete();
	Window->Delete();	
	return (0);
}	
	
