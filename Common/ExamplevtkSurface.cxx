/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */     

#include <vtkCommand.h>
#include "vtkSurface.h"
#include "RenderWindow.h"

// Here is a small sample program showing elementary methods to use vtkSurface

int main( int argc, char *argv[] )
{
	vtkIdType i;			// variable de boucle 
	vtkIdType v1,v2;		// stockage de l'identite de points

	// points and cells definitions
	static double x[6][3]={{0,0,1},{0.707,0.707,0},{0.707,-0.707,0},{-0.707,-0.707,0},{-0.707,0.707,0},{0,0,-1}};
	static int pts[8][3]={{0,1,2},{0,2,3},{0,3,4},{0,4,1},{5,1,2},{5,2,3},{5,3,4},{5,4,1} }; 

	vtkSurface *test=vtkSurface::New();

	// create vertices
	for (i=0;i<6;i++) test->AddVertex(x[i]);

	// create triangles
	for (i=0;i<6;i++) 
		test->AddFace(pts[i][0],pts[i][1],pts[i][2]);

	// display points coordinates
	for (int i=0;i<test->GetNumberOfPoints();i++)
	{
		double Point[3];
		test->GetPoint(i,Point);
		cout<<"Coordinates for vertex "<<i<<" : "<<Point[0]<<" "<<Point[1]<<" "<<Point[2]<<endl;
	}

	// loop on edges
	for (i=0;i<test->GetNumberOfEdges();i++)
	{
		test->GetEdgeVertices(i,v1,v2);
		cout<<"Edge "<<i<<" has vertices "<<v1<<" and "<<v2<<endl;
	}

	// display mesh

	RenderWindow *Window=RenderWindow::New();
	Window->SetInput(test);
	Window->Render();
	Window->Interact();

	// add more triangles
	for (i=6;i<8;i++) 
		test->AddFace(pts[i][0],pts[i][1],pts[i][2]);
	
	//render again
	Window->Render();
	Window->Interact();

	//delete objects
	Window->Delete();
	test->Delete();

	return(0);
}
