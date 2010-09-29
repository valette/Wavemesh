/*=========================================================================

  Program:   Wavemesh : a progressive lossless compression scheme for 3D triagular meshes
  Module:    wavemesh.cxx
  Language:  C++
  Date:      2008/08
  Auteur:    Sebastien Valette
    This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME wavemesh 
// .SECTION Description

#include "vtkMultiresolutionIO.h"

int main( int argc, char *argv[] )
{
	cout << "Wavemesh, a wavelet based 3D mesh compression"<<endl;
	cout << "Copyright (C) 2008 S. Valette @ CREATIS, France "<<endl;
    cout << "(http://www-creatis.insa-lyon.fr/~valette)"<<endl;
	cout << "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY."<<endl;
	cout << "This program uses the VTK Library (www.kitware.com) for visualization"<<endl;
	cout << "This program uses the Range encoder (www.compressconsult.com) for coding"<<endl<<endl;

	if(argc<2)
	{
		cout<<"Usage :"<<endl;
		cout<<"Compression : wavemesh c file [options]"<<endl;
		cout<<"Decompression : wavemesh d file [options]"<<endl;
		cout<<endl;
		cout<<"Optionnal arguments :"<<endl;
		cout<<"-d 0/1/2 changes display type (default : 0):"<<endl;
		cout<<"              0: no display"<<endl;
		cout<<"              1: display all levels"<<endl;
		cout<<"              2: display with stop at each level"<<endl;
		cout<<"-a 0/1 changes arithmetics : 0=floats (zerotree coding)"<<endl;
		cout<<"                             1=integers (original wavemesh coding)"<<endl;
		cout<<"-q n : defines the coordinates quantization to n bits (default : 12 bits)"<<endl;
		cout<<"-b n : defines the number of bitplanes for zerotree coding (default : 12)"<<endl;
		cout<<"-l r : defines the lifting radius (-1=no lifting. Default : -1)"<<endl;
		cout<<"-g 0/1 : enables/disable geometric constraints for simplification"<<endl;
		cout<<"-et threshold : defines the edge angle threshold (default : 0.3)"<<endl;
		cout<<"-wt threshold : defines the wavelet ratio threshold (default : 0.25)"<<endl;
		cout<<"-o filename : defines the compressed output file name (default : out.ddd)"<<endl;
		
		return (0);
	}
	else
	{
		cout<<"For mor help on wavemesh, execute the program without any argument"<<endl;
	}

	vtkMultiresolutionIO *MIO=vtkMultiresolutionIO::New();
	int Arithmetics=1;

	// Parse optionnal arguments
	int ArgumentsIndex=3;
	while (ArgumentsIndex<argc)
	{
		if (strcmp(argv[ArgumentsIndex],"-d")==0)
		{
			int Display=atoi(argv[ArgumentsIndex+1]);
            MIO->SetDisplay(Display);
			cout<<"Display="<<Display<<endl;
		}

		if (strcmp(argv[ArgumentsIndex],"-a")==0)
		{
			Arithmetics=atoi(argv[ArgumentsIndex+1]);
			cout<<"Arithmetics (0=floats (zerotree) ; 1=integers (original wavemesh)): "<<Arithmetics<<endl;
			MIO->SetArithmeticType(Arithmetics);
		}

		if (strcmp(argv[ArgumentsIndex],"-q")==0)
		{
			int Quantization=atoi(argv[ArgumentsIndex+1]);
			cout<<"Coordinates Quantization :"<<Quantization<<endl;
			MIO->SetQuantization(Quantization);
		}

		if (strcmp(argv[ArgumentsIndex],"-b")==0)
		{
			int NumberOfBitPlanes=atoi(argv[ArgumentsIndex+1]);
			cout<<"Number Of BitPlanes :"<<NumberOfBitPlanes<<endl;
			MIO->SetNumberOfBitPlanes(NumberOfBitPlanes);
		}

		if (strcmp(argv[ArgumentsIndex],"-l")==0)
		{
			int Lifting=atoi(argv[ArgumentsIndex+1]);
			if (Lifting==-1)
			{
				cout<<"Lifting scheme deactivated"<<endl;
				MIO->SetLifting(0);
			}
			else
			{
				MIO->SetLifting(1);
				MIO->SetLiftingRadius(Lifting);
				cout<<"Lifting radius : "<<Lifting<<endl;
			}
		}

		if (strcmp(argv[ArgumentsIndex],"-g")==0)
		{
			int Geometry=atoi(argv[ArgumentsIndex+1]);
			cout<<"Geometric constraints for simplification :"<<Geometry<<endl;
			MIO->SetGeometricalConstraint(Geometry);
		}
		
		if (strcmp(argv[ArgumentsIndex],"-et")==0)
		{
			double Threshold=atof(argv[ArgumentsIndex+1]);
			cout<<"Edge angle threshold :"<<Threshold<<endl;
			MIO->SetEdgeAngleThreshold(Threshold);
		}

		if (strcmp(argv[ArgumentsIndex],"-wt")==0)
		{
			double Threshold=atof(argv[ArgumentsIndex+1]);
			cout<<"Wavelet ratio threshold :"<<Threshold<<endl;
			MIO->SetWGC(Threshold);
		}

		if (strcmp(argv[ArgumentsIndex],"-o")==0)
		{
			cout<<"Output File : "<<argv[ArgumentsIndex+1]<<endl;
			MIO->SetFileName(argv[ArgumentsIndex+1]);
		}					
		ArgumentsIndex+=2;
	}

	if (strcmp(argv[1],"c")==0)
	{
		vtkSurface *Mesh = vtkSurface::New();
		cout<<"Load : "<<argv[2]<<endl;
		Mesh->CreateFromFile(argv[2]);
		Mesh->DisplayMeshProperties();

		if (Arithmetics==1)
			Mesh->QuantizeCoordinates(MIO->GetQuantization());

		MIO->SetInput(Mesh);

		MIO->Analyse();
		MIO->Synthetize();
		MIO->Approximate();
		MIO->Write();
		Mesh->Delete();
	}
	else
	{
		MIO->SetFileName(argv[2]);
		MIO->Read();
	}

	MIO->Delete();
}

