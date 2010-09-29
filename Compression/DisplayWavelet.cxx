/*=========================================================================

  Program:  DisplayWavelet
  Module:    DisplayWavelet.cxx
  Language:  C++
  Date:      2010/09
  Auteur:    Sebastien Valette
  This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME DisplayWavelet 
// .SECTION Description

#include "RenderWindow.h"
#include "vtkWaveletSubdivisionFilter.h"

int main( int argc, char *argv[] )
{

	vtkSurface *Mesh = vtkSurface::New();
	Mesh->AddVertex(-1,-0.8,0);
	Mesh->AddVertex(1,-0.8,0);
	Mesh->AddVertex(0,1,0);
	Mesh->AddFace(0,1,2);

	int NumberOfFilters=5;
	vtkWaveletSubdivisionFilter *Filters[NumberOfFilters];

	for (int i=0;i<NumberOfFilters;i++)
		Filters[i]=vtkWaveletSubdivisionFilter::New();

	int Arithmetics=1;

	// Parse optionnal arguments
	int ArgumentsIndex=1;
	while (ArgumentsIndex<argc)
	{
		if (strcmp(argv[ArgumentsIndex],"-a")==0)
		{
			Arithmetics=atoi(argv[ArgumentsIndex+1]);
			cout<<"Arithmetics (0=floats (zerotree) ; 1=integers (original wavemesh)): "<<Arithmetics<<endl;
				for (int i=0;i<NumberOfFilters;i++)
					Filters[i]->SetArithmeticType(Arithmetics);
		}

		if (strcmp(argv[ArgumentsIndex],"-l")==0)
		{
			int Lifting=atoi(argv[ArgumentsIndex+1]);
			if (Lifting==-1)
			{
				cout<<"Lifting scheme deactivated"<<endl;
				for (int i=0;i<NumberOfFilters;i++)
					Filters[i]->SetLifting(0);
			}
			else
			{
				for (int i=0;i<NumberOfFilters;i++)
				{
					Filters[i]->SetLifting(1);
					Filters[i]->SetLiftingRadius(Lifting);
				}
				cout<<"Lifting radius : "<<Lifting<<endl;
			}
		}

	ArgumentsIndex+=2;
	}

	Filters[0]->SetInput(Mesh);
	for (int i=0;i<NumberOfFilters;i++)
	{
		Filters[i]->Subdivide();
		if (i!=NumberOfFilters-1)
			Filters[i+1]->SetInput(Filters[i]->GetOutput());
	}


	Filters[NumberOfFilters-1]->DisplayWavelet();
	
	RenderWindow *Window=RenderWindow::New();
	Window->SetInput(Filters[NumberOfFilters-1]->GetOutput());
	
	Window->Render();
	Window->Interact();

}

