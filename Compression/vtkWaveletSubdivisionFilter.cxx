//  This software is governed by the GPL license (see License.txt)

#include <vector>
#include <math.h>

#include <vtkPolyDataMapper.h>
#include <vtkCellData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLookupTable.h>                                                     
#include <vtkActor.h>                                                           
#include <vtkRenderer.h>                                                        
#include <vtkCamera.h>

#include "vtkWaveletSubdivisionFilter.h"
#include "RenderWindow.h"


// first include the required header files for the vtk classes we are using

//--------------------------------------------------------------------------
vtkWaveletSubdivisionFilter* vtkWaveletSubdivisionFilter::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkWaveletSubdivisionFilter");
	if(ret)
	{
		return (vtkWaveletSubdivisionFilter*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return new vtkWaveletSubdivisionFilter;
}

void vtkWaveletSubdivisionFilter::SaveWavelets(const char *Filename)
{
	int i;
	fstream	WavF;
	WavF.open (Filename, ofstream::out | ofstream::trunc);
	WavF<<"Input: "<<this->SubdivisionInput->GetNumberOfPoints()<<" points"<<endl;
	WavF<<"Output: "<<this->SubdivisionOutput->GetNumberOfPoints()<<" points"<<endl;
	for	(i=0;i<this->IntegerWavelets->GetNumberOfTuples();i++)
	{
		double	value1,value2,value3,Wav[3];
		this->IntegerWavelets->GetTuple(i,Wav);
		value1=Wav[0];
		value2=Wav[1];
		value3=Wav[2];
		WavF<<value1<<" "<<value2<<" "<<value3<<endl;
	}
	WavF.close();
}

void vtkWaveletSubdivisionFilter::WriteCoefficients()
{
	qsmodel QsmGeometry,QsmGeometry2,QsmGeometry3;
	int i,j;
	if (this->ArithmeticType==0)
		return;

	double Wavelet[3];
	int ch;
	short max=-1000,min=1000,size=0;
	char sign;
	for (i=0;i<this->IntegerWavelets->GetNumberOfTuples();i++)
	{
		this->IntegerWavelets->GetTuple(i,Wavelet);
		for (j=0;j<3;j++)
		{
			if (Wavelet[j]<min) min=(short) Wavelet[j];
			if (Wavelet[j]>max) max=(short) Wavelet[j];
		}
	}
	size=max-min+1;
	std::vector<int> table;
	table.resize(size);

	for (i=0;i<size;i++) table[i]=0;

	for (i=0;i<this->IntegerWavelets->GetNumberOfTuples();i++)
	{			
		this->IntegerWavelets->GetTuple(i,Wavelet);				
		for (j=0;j<3;j++)
			table[(int)Wavelet[j]-min]++;
	}

	double entropy=0,p,denom;

	denom=this->IntegerWavelets->GetNumberOfTuples()*3;

	for (i=0;i<size;i++)
	{
		if (table[i]>0)
		{
			p=(double) table[i]/denom;
			entropy+=-p*log(p)/log(2.0);
		}
	}

	this->EstimatedGeometryBits=entropy*3*this->IntegerWavelets->GetNumberOfTuples();
	table.empty();

	QsmGeometry.initqsmodel(size,15,1000,NULL,1);

	if (this->CoordinatesCoupling==1)
	{
		QsmGeometry2.initqsmodel(size,15,1000,NULL,1);
		QsmGeometry3.initqsmodel(size,15,1000,NULL,1);
	}

	if (min<0)
	{
		min=-min;
		sign=1;
	}
	else sign=0;

	this->ArithmeticCoder->EncodeByte(sign);
	this->ArithmeticCoder->EncodeWord(min);
	this->ArithmeticCoder->EncodeWord(size);

	if (sign==1) min=-min;

	if (this->CoordinatesCoupling==0)
	{
		for (i=0;i<this->IntegerWavelets->GetNumberOfTuples();i++)
		{
			this->IntegerWavelets->GetTuple(i,Wavelet);
			for (j=0;j<3;j++)
			{
				Wavelet[j]-=min;
				ch=(int)Wavelet[j];
				this->ArithmeticCoder->Encode(ch,&QsmGeometry);
			}
		}
	}
	else
	{
		for (i=0;i<this->IntegerWavelets->GetNumberOfTuples();i++)
		{
			this->IntegerWavelets->GetTuple(i,Wavelet);
			Wavelet[0]-=min;
			ch=(int)Wavelet[0];
			this->ArithmeticCoder->Encode(ch,&QsmGeometry);

			Wavelet[1]-=min;
			ch=(int)Wavelet[1];
			this->ArithmeticCoder->Encode(ch,&QsmGeometry2);

			Wavelet[2]-=min;
			ch=(int)Wavelet[2];
			this->ArithmeticCoder->Encode(ch,&QsmGeometry3);
		}
	}
}

void vtkWaveletSubdivisionFilter::ReadCoefficients()
{
	qsmodel QsmGeometry,QsmGeometry2,QsmGeometry3;

	vtkSurface *Out=this->SubdivisionOutput;
	vtkSurface *Input=this->SubdivisionInput;
	int i,j;

	if (this->ArithmeticType==1)
	{
		double Wavelet[3];
		int ch;
		short min,max;
		char sign;


		sign=this->ArithmeticCoder->DecodeByte();
		min=this->ArithmeticCoder->DecodeWord();
		max=this->ArithmeticCoder->DecodeWord();

		if (sign==1) min=-min;

		QsmGeometry.initqsmodel(max,15,1000,NULL,0);

		if(this->CoordinatesCoupling==1)
		{
			QsmGeometry2.initqsmodel(max,15,1000,NULL,0);
			QsmGeometry3.initqsmodel(max,15,1000,NULL,0);
		}

		if (!this->IntegerWavelets)
		{
			this->IntegerWavelets=vtkIntArray::New();
			this->IntegerWavelets->SetNumberOfComponents(3);
			this->IntegerWavelets->SetNumberOfTuples(Out->GetNumberOfPoints()-Input->GetNumberOfPoints());
		}

		if (this->CoordinatesCoupling==0)
		{
			for (i=0;i<Out->GetNumberOfPoints()-Input->GetNumberOfPoints();i++)
			{
				for (j=0;j<3;j++)
				{
					ch=this->ArithmeticCoder->Decode(&QsmGeometry);
					Wavelet[j]=ch+min;
				}
				this->IntegerWavelets->SetTuple(i,Wavelet);
			}
		}
		else
		{
			for (i=0;i<Out->GetNumberOfPoints()-Input->GetNumberOfPoints();i++)
			{
				ch=this->ArithmeticCoder->Decode(&QsmGeometry);
				Wavelet[0]=ch+min;
				ch=this->ArithmeticCoder->Decode(&QsmGeometry2);
				Wavelet[1]=ch+min;
				ch=this->ArithmeticCoder->Decode(&QsmGeometry3);
				Wavelet[2]=ch+min;
				this->IntegerWavelets->SetTuple(i,Wavelet);
			}
		}
	}
}

void vtkWaveletSubdivisionFilter::SetInput(vtkSurface *input)
{
	this->vtkProcessObject::SetNthInput(0, input);
	this->SubdivisionInput=input;
}

void vtkWaveletSubdivisionFilter::SolveInverseProblem (vtkIdType initial_face)
{

	vtkLookupTable *lut=0;
	RenderWindow *Window=0;

	double TransparencyForConqueredFaces=0.4;
	double TransparencyForConqueredVertices=0.4;

	if (this->DisplayEfficiency==1)
	{
		Window=RenderWindow::New();
		double x[4];
		int i;
		this->CellTypes->SetNumberOfValues(this->MergeInput->GetNumberOfCells());
		for (i=0;i<MergeInput->GetNumberOfCells();i++)
			this->CellTypes->SetValue(i,0);
		lut=vtkLookupTable::New();
		lut->SetNumberOfTableValues(10);

		// The value 0 is for faces that are not part of the solution
		// they are the problem...
		x[0]=0;
		x[1]=0;
		x[2]=0;
		x[3]=1;
		lut->SetTableValue (0,x);

		// The value 1 is for faces that remain intact 
		x[0]=1;
		x[1]=1;
		x[2]=1;
		x[3]=TransparencyForConqueredFaces;
		lut->SetTableValue (1,x);

		// The value 2 is for faces that are merged two by two
		x[0]=0;
		x[1]=0;
		x[2]=1;
		x[3]=TransparencyForConqueredFaces;
		lut->SetTableValue (2,x);

		// The value 3 is for faces that are merged 3 by 3
		x[0]=0;
		x[1]=1;
		x[2]=1;
		x[3]=TransparencyForConqueredFaces;
		lut->SetTableValue (3,x);

		// The value 4 is for faces that are merged 4 by 4
		x[0]=0;
		x[1]=1;
		x[2]=0;
		x[3]=TransparencyForConqueredFaces;
		lut->SetTableValue (4,x);

		// The value 5 is for faces merged with an edge flip (1st kind)
		x[0]=1;
		x[1]=0;
		x[2]=0;
		x[3]=TransparencyForConqueredFaces;
		lut->SetTableValue (5,x);

		// The value 6 is for faces merged with an edge flip (2nd kind)
		x[0]=1;
		x[1]=0;
		x[2]=0;
		x[3]=TransparencyForConqueredFaces;
		lut->SetTableValue (6,x);

		// The value 7 is for parent vertices
		x[0]=0;
		x[1]=0;
		x[2]=1;
		x[3]=TransparencyForConqueredVertices;
		lut->SetTableValue (7,x);

		// The value 8 is for child vertices
		x[0]=1;
		x[1]=1;
		x[2]=1;
		x[3]=TransparencyForConqueredVertices;
		lut->SetTableValue (8,x);

		lut->SetRange(0.0,9.0);
	}

	this->SolveInverseSubdivision4(initial_face);
	if (this->DisplayEfficiency==1)
	{
		vtkPolyData *CopyMesh=vtkPolyData::New();
		CopyMesh->ShallowCopy(this->MergeInput);
		CopyMesh->GetCellData()->SetScalars(this->CellTypes);
		vtkIdType i;
		vtkIdList *VList=vtkIdList::New();
		VList->SetNumberOfIds(1);

		vtkPolyData *VerticesMesh=vtkPolyData::New();
		VerticesMesh->Allocate(this->MergeInput->GetNumberOfPoints());
		VerticesMesh->SetPoints(CopyMesh->GetPoints());

		vtkIntArray *VerticesTypes=vtkIntArray::New();
		VerticesTypes->SetNumberOfValues(MergeInput->GetNumberOfPoints());
		for (i=0;i<MergeInput->GetNumberOfPoints();i++)
		{
			if (this->vertices[i].type==PARENT)
				VerticesTypes->SetValue(i,7);
			else
				VerticesTypes->SetValue(i,8);
		}
		VerticesMesh->GetCellData()->SetScalars(VerticesTypes);

		for (i=0;i<this->MergeInput->GetNumberOfPoints();i++)
		{
			VerticesMesh->InsertNextCell(VTK_VERTEX,1,&i);
		}

		Window->SetInput(CopyMesh);
		Window->SetLookupTable(lut);
		Window->AddPolyData(VerticesMesh);
		Window->Render();
		Window->SetWindowName("Simplification Efficiency");
		Window->Interact();
		Window->Delete();
		lut->Delete();
		CopyMesh->Delete();
		VerticesMesh->Delete();
		VerticesTypes->Delete();
	}
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
vtkPolyData *vtkWaveletSubdivisionFilter::GetInput()
{
	if (this->NumberOfInputs < 1) return NULL;

	return (vtkSurface *)(this->Inputs[0]);
}

void vtkWaveletSubdivisionFilter::ComputeFastLiftingCoeff()
{
	double N1,N2,R;
	N1=this->SubdivisionInput->GetNumberOfCells();
	N2=this->SubdivisionOutput->GetNumberOfCells();
	R=N2/N1;
	N1=floor(0.5+21.0-9.0*(exp(2.0*R-8.0)));
	this->LiftNum=(int) N1;
}

void vtkWaveletSubdivisionFilter::ComputeVertexContribution(vtkIdType p1,vtkIdType p2, vtkFloatingPointType *V1)
{
	vtkIdType Edge,v1,v2,f1,f2;
	vtkFloatingPointType Vertex[3];
	int Valence;

	vtkPoints *Points=this->SubdivisionInput->GetPoints();

	Points->GetPoint(p1, Vertex);

	if (this->SubdivisionInput->GetNumberOfBoundaries(p1)>0)
	{
		// TO DO!!!
		return;
	}


	Valence=this->SubdivisionInput->GetValence(p1);
	if (Valence==3)
	{
		V1[0]+=2.0*Vertex[0]/12.0;
		V1[1]+=2.0*Vertex[1]/12.0;
		V1[2]+=2.0*Vertex[2]/12.0;

		Edge=this->SubdivisionInput->IsEdge(p1,p2);
		this->SubdivisionInput->GetEdgeFaces(Edge,f1,f2);

		v1=this->SubdivisionInput->GetThirdPoint(f1,p1,p2);
		Points->GetPoint(v1, Vertex);
		V1[0]-=Vertex[0]/12.0;
		V1[1]-=Vertex[1]/12.0;
		V1[2]-=Vertex[2]/12.0;

		v1=this->SubdivisionInput->GetThirdPoint(f2,p1,p2);
		Points->GetPoint(v1, Vertex);	
		V1[0]-=Vertex[0]/12.0;
		V1[1]-=Vertex[1]/12.0;
		V1[2]-=Vertex[2]/12.0;

		return;
	}
	if (Valence==4)
	{
		V1[0]+=1.0*Vertex[0]/8.0;
		V1[1]+=1.0*Vertex[1]/8.0;
		V1[2]+=1.0*Vertex[2]/8.0;

		Edge=this->SubdivisionInput->IsEdge(p1,p2);
		this->SubdivisionInput->GetEdgeFaces(Edge,f1,f2);
		v1=this->SubdivisionInput->GetThirdPoint(f1,p1,p2);
		this->SubdivisionInput->Conquer(f1,p1,v1,f2,v2);	
		Points->GetPoint(v2, Vertex);	
		V1[0]-=Vertex[0]/8.0;
		V1[1]-=Vertex[1]/8.0;
		V1[2]-=Vertex[2]/8.0;


		return;

	}

	//	if (Valence==6)
	{

		double ratio=16.0;
		Edge=this->SubdivisionInput->IsEdge(p1,p2);
		this->SubdivisionInput->GetEdgeFaces(Edge,f1,f2);
		v1=this->SubdivisionInput->GetThirdPoint(f1,p1,p2);
		Points->GetPoint(v1, Vertex);
		V1[0]+=Vertex[0]/ratio;
		V1[1]+=Vertex[1]/ratio;
		V1[2]+=Vertex[2]/ratio;

		this->SubdivisionInput->Conquer(f1,p1,v1,f2,v2);
		Points->GetPoint(v2, Vertex);
		V1[0]-=Vertex[0]/ratio;
		V1[1]-=Vertex[1]/ratio;
		V1[2]-=Vertex[2]/ratio;

		this->SubdivisionInput->GetEdgeFaces(Edge,f1,f2);
		v1=this->SubdivisionInput->GetThirdPoint(f2,p1,p2);

		Points->GetPoint(v1, Vertex);
		V1[0]+=Vertex[0]/ratio;
		V1[1]+=Vertex[1]/ratio;
		V1[2]+=Vertex[2]/ratio;

		this->SubdivisionInput->Conquer(f2,p1,v1,f1,v2);
		Points->GetPoint(v2, Vertex);
		V1[0]-=Vertex[0]/16.0;
		V1[1]-=Vertex[1]/16.0;
		V1[2]-=Vertex[2]/16.0;

		return;
	}

}
void vtkWaveletSubdivisionFilter::ComputeVertexContribution2(vtkIdType p1,vtkIdType p2, double *V1)
{
	vtkIdType Edge,v1,v2,f1,f2;
	vtkFloatingPointType Vertex[3];	
	int Valence;

	vtkPoints *Points=this->SubdivisionInput->GetPoints();
	vtkSurface *MInput=this->SubdivisionInput;

	Points->GetPoint(p1, Vertex);	

	if (MInput->GetNumberOfBoundaries(p1)>0)
	{
		// TO DO!!!
		return;
	}


	Valence=MInput->GetValence(p1);
	if (Valence==3)
	{
		Points->GetPoint(p1, Vertex);	
		V1[0]+=3.0*Vertex[0]/12.0;
		V1[1]+=3.0*Vertex[1]/12.0;
		V1[2]+=3.0*Vertex[2]/12.0;

		Points->GetPoint(p2, Vertex);	
		V1[0]-=Vertex[0]/12.0;
		V1[1]-=Vertex[1]/12.0;
		V1[2]-=Vertex[2]/12.0;

		Edge=MInput->IsEdge(p1,p2);
		MInput->GetEdgeFaces(Edge,f1,f2);
		v1=MInput->GetThirdPoint(f1,p1,p2);
		Points->GetPoint(v1, Vertex);	
		V1[0]-=Vertex[0]/12.0;
		V1[1]-=Vertex[1]/12.0;
		V1[2]-=Vertex[2]/12.0;


		v1=MInput->GetThirdPoint(f2,p1,p2);
		Points->GetPoint(v1, Vertex);	
		V1[0]-=Vertex[0]/12.0;
		V1[1]-=Vertex[1]/12.0;
		V1[2]-=Vertex[2]/12.0;


		return;

	}
	if (Valence==4)
	{
		Points->GetPoint(p1, Vertex);	
		V1[0]+=2.0*Vertex[0]/8.0;
		V1[1]+=2.0*Vertex[1]/8.0;
		V1[2]+=2.0*Vertex[2]/8.0;

		Points->GetPoint(p2, Vertex);	
		V1[0]-=Vertex[0]/8.0;
		V1[1]-=Vertex[1]/8.0;
		V1[2]-=Vertex[2]/8.0;

		Edge=MInput->IsEdge(p1,p2);
		MInput->GetEdgeFaces(Edge,f1,f2);
		v1=MInput->GetThirdPoint(f1,p1,p2);
		MInput->Conquer(f1,p1,v1,f2,v2);
		Points->GetPoint(v2, Vertex);	

		V1[0]-=Vertex[0]/8.0;
		V1[1]-=Vertex[1]/8.0;
		V1[2]-=Vertex[2]/8.0;

		return;

	}

	{

		int j;
		double Sum,Sj;
		double pi=4.0*atan(1.0);

		Edge=MInput->IsEdge(p1,p2);
		MInput->GetEdgeFaces(Edge,f1,f2);
		v1=MInput->GetThirdPoint(f1,p1,p2);	
		Points->GetPoint(p2, Vertex);	
		Sj=5.0/4.0+1.0/(2*(double)Valence)-0.5;

		Sum=Sj;
		V1[0]+=Vertex[0]*Sj;
		V1[1]+=Vertex[1]*Sj;
		V1[2]+=Vertex[2]*Sj;
		for (j=1;j<Valence;j++)
		{
			Points->GetPoint(v1, Vertex);	
			Sj=0.25+cos(2.0*pi*(double)j/(double)Valence)
				+0.5*cos(4.0*pi*(double)j/(double)Valence)/(double)Valence;
			MInput->Conquer(f1,p1,v1,f2,v2);
			v1=v2;
			f1=f2;
			Sum+=Sj;
			V1[0]+=Vertex[0]*Sj;
			V1[1]+=Vertex[1]*Sj;
			V1[2]+=Vertex[2]*Sj;
		}
		Points->GetPoint(p1, Vertex);	
		V1[0]-=Vertex[0]*Sum;
		V1[1]-=Vertex[1]*Sum;
		V1[2]-=Vertex[2]*Sum;

		return;
	}
}

void vtkWaveletSubdivisionFilter::ComputeNewVertexCoordinates(vtkIdType p1,vtkIdType p2, vtkFloatingPointType *V1)
{
	vtkFloatingPointType P1[3], P2[3], Vertex[3];		

	int Val1,Val2,Val;
	int approach=1;
	vtkIdType Edge, v1,v2,f1,f2;

	vtkSurface *MInput=this->SubdivisionInput;
	vtkPoints *Points=MInput->GetPoints();

	if (approach==1)
	{
		V1[0]=0;
		V1[1]=0;
		V1[2]=0;

		if ((MInput->GetNumberOfBoundaries(p1)==0)&&
			(MInput->GetNumberOfBoundaries(p2)==0))
		{
			Val1=MInput->GetValence(p1);
			if (Val1==6)
			{
				Val2=MInput->GetValence(p2);			
				if (Val2==6)	
				{
					// here Valence(p1)=6 and Valence (p2=6)

					Edge=MInput->IsEdge(p1,p2);
					MInput->GetEdgeFaces(Edge,f1,f2);
					v1=MInput->GetThirdPoint(f1,p1,p2);				
					Points->GetPoint(v1, Vertex);	
					V1[0]+=Vertex[0]/8.0;
					V1[1]+=Vertex[1]/8.0;
					V1[2]+=Vertex[2]/8.0;


					MInput->Conquer(f1,p1,v1,f2,v2);
					Points->GetPoint(v2, Vertex);	
					V1[0]-=Vertex[0]/16.0;
					V1[1]-=Vertex[1]/16.0;
					V1[2]-=Vertex[2]/16.0;

					MInput->Conquer(f1,p2,v1,f2,v2);
					Points->GetPoint(v2, Vertex);	
					V1[0]-=Vertex[0]/16.0;
					V1[1]-=Vertex[1]/16.0;
					V1[2]-=Vertex[2]/16.0;

					MInput->GetEdgeFaces(Edge,f1,f2);
					v1=MInput->GetThirdPoint(f2,p1,p2);
					Points->GetPoint(v1, Vertex);	
					V1[0]+=Vertex[0]/8.0;
					V1[1]+=Vertex[1]/8.0;
					V1[2]+=Vertex[2]/8.0;


					MInput->Conquer(f2,p1,v1,f1,v2);
					Points->GetPoint(v2, Vertex);	
					V1[0]-=Vertex[0]/16.0;
					V1[1]-=Vertex[1]/16.0;
					V1[2]-=Vertex[2]/16.0;

					MInput->Conquer(f2,p2,v1,f1,v2);
					Points->GetPoint(v2, Vertex);	
					V1[0]-=Vertex[0]/16.0;
					V1[1]-=Vertex[1]/16.0;
					V1[2]-=Vertex[2]/16.0;
					return;
				}
				else
				{
					v1=p1;
					p1=p2;
					p2=v1;
					Val=Val1;
					Val1=Val2;
					Val2=Val;
				}
			}
			else
			{
				Val2=MInput->GetValence(p2);
				if (Val2!=6)
				{
					this->ComputeVertexContribution(p1,p2,V1);
					this->ComputeVertexContribution(p2,p1,V1);
					return;
				}
			}

			// here Valence(p1)!=6 and Valence(p2)=6

			vtkIdType Test;
			Test=MInput->IsEdge(p1,p2);
			this->ComputeVertexContribution(p1,p2,V1);
			this->ComputeVertexContribution(p2,p1,V1);
			V1[0]=V1[0]/2.0;
			V1[1]=V1[1]/2.0;
			V1[2]=V1[2]/2.0;
			return;
		}
	}
	else 
	{
		this->SubdivisionOutput->GetPoints()->GetPoint(p1, P1);	
		this->SubdivisionOutput->GetPoints()->GetPoint(p2, P2);	

		if (this->GeometryPrediction==0)
		{
			V1[0]=0;
			V1[1]=0;
			V1[2]=0;
		}
		else
		{
			V1[0]=0;
			V1[1]=0;
			V1[2]=0;

			this->ComputeVertexContribution(p1,p2,V1);
			this->ComputeVertexContribution(p2,p1,V1);
		}
	}
}


void vtkWaveletSubdivisionFilter::Reconstruct()
{
	vtkIdType PtId,p1,p2,NNewPoints,NOldPoints;
	double *Wavelet;
	vtkFloatingPointType P1[3], P2[3], Po1[3], Po2[3], V1[3];		
	vtkFloatingPointType V[3];
	int *IntWavelet;

	vtkPoints *NewPoints=this->SubdivisionOutput->GetPoints();
	vtkPoints *OldPoints=this->SubdivisionInput->GetPoints();
	vtkSparseMatrix::NonZeroElement *NonZero;


	NNewPoints=this->SubdivisionOutput->GetNumberOfPoints();
	NOldPoints=this->SubdivisionInput->GetNumberOfPoints();

	vtkDoubleArray *Wavelets=this->Wavelets;

	if ((this->GeometryPrediction==1)&&(this->ArithmeticType==0))
	{
		Wavelets=vtkDoubleArray::New();
		Wavelets->DeepCopy(this->Wavelets);
		vtkFloatingPointType Vertex[3];

		for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
		{
			p1=this->vertices[PtId].parent1;
			p2=this->vertices[PtId].parent2;
			this->ComputeNewVertexCoordinates(p1,p2,Vertex);

			Wavelet=Wavelets->GetPointer((PtId-NOldPoints)*3);

			Wavelet[0]+=Vertex[0];
			Wavelet[1]+=Vertex[1];
			Wavelet[2]+=Vertex[2];
		}
	}
	else
	{
		Wavelets=this->Wavelets;
	}

	if (this->ProcessLifting==0)
	{
		for (PtId=0;PtId<NOldPoints;PtId++)
		{
			OldPoints->GetPoint(PtId, P1);	
			P2[0]=P1[0];
			P2[1]=P1[1];
			P2[2]=P1[2];
			NewPoints->SetPoint(PtId, P2);	
		}
	}


	// Inversion de l'approximation (lifting classique):
	if (this->ProcessLifting==1)
	{
		if (this->ArithmeticType==0)
		{
			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				OldPoints->GetPoint(PtId, P1);				
				NewPoints->GetPoint(PtId, P2);	
				P2[0]=P1[0];
				P2[1]=P1[1];
				P2[2]=P1[2];

				NonZero=this->Alpha->FirstHor[PtId];

				while (NonZero)
				{
					Wavelet=Wavelets->GetPointer((NonZero->j)*3);

					P2[0]-=NonZero->Value*Wavelet[0];
					P2[1]-=NonZero->Value*Wavelet[1];
					P2[2]-=NonZero->Value*Wavelet[2];

					NonZero=NonZero->NextHor;
				}
				NewPoints->SetPoint(PtId, P2);	
			}
		}
		else	// if (this->ArithmeticType==1)
		{
			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				V[0]=0;
				V[1]=0;
				V[2]=0;

				NonZero=this->Alpha->FirstHor[PtId];

				while (NonZero)
				{
					IntWavelet=IntegerWavelets->GetPointer((NonZero->j)*3);

					V[0]+=NonZero->Value*IntWavelet[0];
					V[1]+=NonZero->Value*IntWavelet[1];
					V[2]+=NonZero->Value*IntWavelet[2];

					NonZero=NonZero->NextHor;
				}			
				OldPoints->GetPoint(PtId, P1);				
				NewPoints->GetPoint(PtId, P2);					
				P2[0]=P1[0]-floor(V[0]+0.5);
				P2[1]=P1[1]-floor(V[1]+0.5);
				P2[2]=P1[2]-floor(V[2]+0.5);
				NewPoints->SetPoint(PtId, P2);	
			}
		}
	}
	// Inversion de l'approximation (lifting rapide):
	if (this->ProcessLifting==2)
	{
		this->ComputeFastLiftingCoeff();

		if (this->ArithmeticType==0)
		{
			double k,den;
			k=this->LiftNum;
			den=this->LiftDen;
			k=k/den;
			V[0]=0;
			V[1]=0;
			V[2]=0;

			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				NewPoints->SetPoint(PtId,V);
			}

			for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
			{
				p1=this->vertices[PtId].parent1;
				p2=this->vertices[PtId].parent2;

				Wavelet=Wavelets->GetPointer((PtId-NOldPoints)*3);
				NewPoints->GetPoint(p1, Po1);				
				NewPoints->GetPoint(p2, Po2);	

				Po1[0]-=(vtkFloatingPointType)Wavelet[0]*k;
				Po1[1]-=(vtkFloatingPointType)Wavelet[1]*k;
				Po1[2]-=(vtkFloatingPointType)Wavelet[2]*k;	

				Po2[0]-=(vtkFloatingPointType)Wavelet[0]*k;
				Po2[1]-=(vtkFloatingPointType)Wavelet[1]*k;
				Po2[2]-=(vtkFloatingPointType)Wavelet[2]*k;

				NewPoints->SetPoint(p1, Po1);				
				NewPoints->SetPoint(p2, Po2);	

			}

			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				OldPoints->GetPoint(PtId, P1);				
				NewPoints->GetPoint(PtId, Po1);	

				Po1[0]=Po1[0]+P1[0];
				Po1[1]=Po1[1]+P1[1];
				Po1[2]=Po1[2]+P1[2];

				NewPoints->SetPoint(PtId, Po1);	
			}
		}
		else	// if (this->ArithmeticType==1)
		{
			double k,den;
			k=this->LiftNum;
			den=this->LiftDen;
			k=k/den;
			V[0]=0;
			V[1]=0;
			V[2]=0;

			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				NewPoints->SetPoint(PtId,V);
			}

			for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
			{
				p1=this->vertices[PtId].parent1;
				p2=this->vertices[PtId].parent2;

				IntWavelet=IntegerWavelets->GetPointer((PtId-NOldPoints)*3);
				NewPoints->GetPoint(p1, Po1);				
				NewPoints->GetPoint(p2, Po2);	

				Po1[0]-=(vtkFloatingPointType)IntWavelet[0]*k;
				Po1[1]-=(vtkFloatingPointType)IntWavelet[1]*k;
				Po1[2]-=(vtkFloatingPointType)IntWavelet[2]*k;	

				Po2[0]-=(vtkFloatingPointType)IntWavelet[0]*k;
				Po2[1]-=(vtkFloatingPointType)IntWavelet[1]*k;
				Po2[2]-=(vtkFloatingPointType)IntWavelet[2]*k;

				NewPoints->SetPoint(p1, Po1);				
				NewPoints->SetPoint(p2, Po2);	

			}

			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				OldPoints->GetPoint(PtId, P1);					
				NewPoints->GetPoint(PtId, Po1);		

				Po1[0]=(vtkFloatingPointType)ceil(Po1[0])+P1[0];
				Po1[1]=(vtkFloatingPointType)ceil(Po1[1])+P1[1];
				Po1[2]=(vtkFloatingPointType)ceil(Po1[2])+P1[2];

				NewPoints->SetPoint(PtId, Po1);		

			}
		}
	}

	// Calcul des nouveaux points:
	if (this->ArithmeticType==0)
	{		
		for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
		{
			p1=this->vertices[PtId].parent1;
			p2=this->vertices[PtId].parent2;

			NewPoints->GetPoint(p1, P1);	
			NewPoints->GetPoint(p2, P2);	
			NewPoints->GetPoint(PtId, V1);	


			V1[0]=0.5*(P1[0]+P2[0]);
			V1[1]=0.5*(P1[1]+P2[1]);
			V1[2]=0.5*(P1[2]+P2[2]);

			Wavelet=Wavelets->GetPointer((PtId-NOldPoints)*3);

			V1[0]+=Wavelet[0];
			V1[1]+=Wavelet[1];
			V1[2]+=Wavelet[2];

			NewPoints->SetPoint(PtId, V1);	

		}
	}
	else	// if (this->ArithmeticType==1)
	{
		for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
		{
			p1=this->vertices[PtId].parent1;
			p2=this->vertices[PtId].parent2;

			NewPoints->GetPoint(p1, P1);	
			NewPoints->GetPoint(p2, P2);	

			V1[0]=0.5*(P1[0]+P2[0]);
			V1[1]=0.5*(P1[1]+P2[1]);
			V1[2]=0.5*(P1[2]+P2[2]);

			IntWavelet=IntegerWavelets->GetPointer((PtId-NOldPoints)*3);

			V1[0]=ceil(V1[0]+IntWavelet[0]);
			V1[1]=ceil(V1[1]+IntWavelet[1]);
			V1[2]=ceil(V1[2]+IntWavelet[2]);

			NewPoints->SetPoint(PtId, V1);	
		}
	}

	if ((this->GeometryPrediction==1)&&(this->ArithmeticType==0))
	{
		Wavelets->Delete();
	}

	this->SubdivisionOutput->Modified();                   
}

void vtkWaveletSubdivisionFilter::Approximate()
{
	vtkIdType PtId,p1,p2,NNewPoints,NOldPoints;
	double *Wavelet;
	int *IntWavelet;
	double P1[3], P2[3], Po1[3], Po2[3], V1[3];		
	double V[3];

	vtkPoints *NewPoints=this->SubdivisionOutput->GetPoints();
	vtkPoints *OldPoints=this->SubdivisionInput->GetPoints();

	vtkSparseMatrix::NonZeroElement *NonZero;

	NNewPoints=this->SubdivisionOutput->GetNumberOfPoints();
	NOldPoints=this->SubdivisionInput->GetNumberOfPoints();

	if (this->ArithmeticType==0)
	{
		if (!this->Wavelets)
		{
			this->Wavelets=vtkDoubleArray::New();
			this->Wavelets->SetNumberOfComponents(3);
			this->Wavelets->SetNumberOfTuples(NNewPoints-NOldPoints);
		}
	}
	else	//(this->ArithmeticType==1)
	{
		if (!this->IntegerWavelets)
		{
			this->IntegerWavelets=vtkIntArray::New();
			this->IntegerWavelets->SetNumberOfComponents(3);
			this->IntegerWavelets->SetNumberOfTuples(NNewPoints-NOldPoints);
		}
	}


	// Calcul des coefficients d'ondelette

	if (this->ArithmeticType==0)
	{
		for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
		{
			p1=this->vertices[PtId].parent1;
			p2=this->vertices[PtId].parent2;

			NewPoints->GetPoint(PtId, V1);	
			NewPoints->GetPoint(p1, P1);	
			NewPoints->GetPoint(p2, P2);	

			Wavelet=Wavelets->GetPointer((PtId-NOldPoints)*3);

			Wavelet[0]=V1[0]-0.5*(P1[0]+P2[0]);
			Wavelet[1]=V1[1]-0.5*(P1[1]+P2[1]);
			Wavelet[2]=V1[2]-0.5*(P1[2]+P2[2]);
		}
	}
	else
	{
		for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
		{
			p1=this->vertices[PtId].parent1;
			p2=this->vertices[PtId].parent2;

			NewPoints->GetPoint(PtId, V1);	
			NewPoints->GetPoint(p1, P1);	
			NewPoints->GetPoint(p2, P2);	

			IntWavelet=IntegerWavelets->GetPointer((PtId-NOldPoints)*3);

			IntWavelet[0]=(int) floor(V1[0]-0.5*(P1[0]+P2[0]));
			IntWavelet[1]=(int) floor(V1[1]-0.5*(P1[1]+P2[1]));
			IntWavelet[2]=(int) floor(V1[2]-0.5*(P1[2]+P2[2]));
		}
	}

	// Calcul de l'approximation (no lifting):
	if (this->ProcessLifting==0)
	{
		for (PtId=0;PtId<NOldPoints;PtId++)
		{
			NewPoints->GetPoint(PtId, P1);	
			OldPoints->SetPoint(PtId, P1);	
		}
	}

	// Calcul de l'approximation (lifting classique):
	if (this->ProcessLifting==1)
	{
		if (this->ArithmeticType==0)
		{
			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				NewPoints->GetPoint(PtId, P2);	

				NonZero=this->Alpha->FirstHor[PtId];

				while (NonZero)
				{
					Wavelet=Wavelets->GetPointer((NonZero->j)*3);

					P2[0]+=NonZero->Value*Wavelet[0];
					P2[1]+=NonZero->Value*Wavelet[1];
					P2[2]+=NonZero->Value*Wavelet[2];

					NonZero=NonZero->NextHor;
				}

				OldPoints->SetPoint(PtId, P2);	
			}
		}
		else	// if (this->ArithmeticType==1)
		{
			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				V[0]=0;
				V[1]=0;
				V[2]=0;

				NonZero=this->Alpha->FirstHor[PtId];

				while (NonZero)
				{
					IntWavelet=IntegerWavelets->GetPointer((NonZero->j)*3);

					V[0]+=NonZero->Value*IntWavelet[0];
					V[1]+=NonZero->Value*IntWavelet[1];
					V[2]+=NonZero->Value*IntWavelet[2];

					NonZero=NonZero->NextHor;
				}

				NewPoints->GetPoint(PtId, P1);	

				P2[0]=P1[0]+floor(V[0]+0.5);
				P2[1]=P1[1]+floor(V[1]+0.5);
				P2[2]=P1[2]+floor(V[2]+0.5);

				OldPoints->SetPoint(PtId, P2);	
			}
		}
	}

	// Calcul de l'approximation (lifting rapide):
	if (this->ProcessLifting==2)
	{
		this->ComputeFastLiftingCoeff();
		if (this->ArithmeticType==1)
		{	
			double k,den;
			k=this->LiftNum;
			den=this->LiftDen;
			k=k/den;
			V[0]=0;
			V[1]=0;
			V[2]=0;

			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				OldPoints->SetPoint(PtId,V);
			}

			for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
			{
				p1=this->vertices[PtId].parent1;
				p2=this->vertices[PtId].parent2;

				IntWavelet=IntegerWavelets->GetPointer((PtId-NOldPoints)*3);

				OldPoints->GetPoint(p1, Po1);	
				OldPoints->GetPoint(p2, Po2);	

				Po1[0]+=(double)(IntWavelet[0]*k);
				Po1[1]+=(double)(IntWavelet[1]*k);
				Po1[2]+=(double)(IntWavelet[2]*k);

				Po2[0]+=(double)(IntWavelet[0]*k);
				Po2[1]+=(double)(IntWavelet[1]*k);
				Po2[2]+=(double)(IntWavelet[2]*k);

				OldPoints->SetPoint(p1, Po1);	
				OldPoints->SetPoint(p2, Po2);	
			}

			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				NewPoints->GetPoint(PtId, P1);	
				OldPoints->GetPoint(PtId, Po1);	

				Po1[0]=(vtkFloatingPointType)floor(Po1[0])+P1[0];
				Po1[1]=(vtkFloatingPointType)floor(Po1[1])+P1[1];
				Po1[2]=(vtkFloatingPointType)floor(Po1[2])+P1[2];

				NewPoints->SetPoint(PtId, P1);	
				OldPoints->SetPoint(PtId, Po1);	
			}
		}
		else	// if (this->ArithmeticType==0)
		{
			double k,den;
			k=this->LiftNum;
			den=this->LiftDen;
			k=k/den;
			V[0]=0;
			V[1]=0;
			V[2]=0;

			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				OldPoints->SetPoint(PtId,V);
			}

			for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
			{
				p1=this->vertices[PtId].parent1;
				p2=this->vertices[PtId].parent2;

				Wavelet=Wavelets->GetPointer((PtId-NOldPoints)*3);

				OldPoints->GetPoint(p1, Po1);	
				OldPoints->GetPoint(p2, Po2);	

				Po1[0]+=(double)(Wavelet[0]*k);
				Po1[1]+=(double)(Wavelet[1]*k);
				Po1[2]+=(double)(Wavelet[2]*k);

				Po2[0]+=(double)(Wavelet[0]*k);
				Po2[1]+=(double)(Wavelet[1]*k);
				Po2[2]+=(double)(Wavelet[2]*k);

				OldPoints->SetPoint(p1, Po1);	
				OldPoints->SetPoint(p2, Po2);	
			}

			for (PtId=0;PtId<NOldPoints;PtId++)
			{
				NewPoints->GetPoint(PtId, P1);	
				OldPoints->GetPoint(PtId, Po1);	

				Po1[0]=Po1[0]+P1[0];
				Po1[1]=Po1[1]+P1[1];
				Po1[2]=Po1[2]+P1[2];

				OldPoints->SetPoint(PtId, Po1);	
			}
		}
	}

	// Lifting dual (prédiction des coefficients d'ondelette)
	if (this->GeometryPrediction==1)
	{
		vtkFloatingPointType Vertex[3];
		if (this->ArithmeticType==0)
		{
			for (PtId=NOldPoints;PtId<NNewPoints;PtId++)
			{
				Wavelet=Wavelets->GetPointer((PtId-NOldPoints)*3);

				p1=this->vertices[PtId].parent1;
				p2=this->vertices[PtId].parent2;

				this->ComputeNewVertexCoordinates(p1,p2,Vertex);

				Wavelet[0]-=Vertex[0];
				Wavelet[1]-=Vertex[1];
				Wavelet[2]-=Vertex[2];
			}
		}
	}

	this->SubdivisionInput->Modified();
}


void vtkWaveletSubdivisionFilter::EncodeMidpoint(int code)
{
	if (this->RegularSubdivisionFlag==1) return;

	this->ArithmeticCoder->EncodeMidpoint(code);
}

int vtkWaveletSubdivisionFilter::DecodeMidpoint()
{
	if (this->RegularSubdivisionFlag==1)
		return (1);

	return (this->ArithmeticCoder->DecodeMidpoint());
}

void vtkWaveletSubdivisionFilter::EncodeFaceSwap(int code,int pos)
{
	this->ArithmeticCoder->EncodeFaceSwap(code);
}

int vtkWaveletSubdivisionFilter::DecodeFaceSwap(int pos)
{
	return (this->ArithmeticCoder->DecodeFaceSwap());
}

void vtkWaveletSubdivisionFilter::EncodeEdgeSwap(int code,int pos)
{
	this->EncodeFaceSwap(code,pos);

}

int vtkWaveletSubdivisionFilter::DecodeEdgeSwap(int pos)
{
	return(this->DecodeFaceSwap(pos));
}

void vtkWaveletSubdivisionFilter::DisplayWavelet()
{
	int number_of_vertices,number_of_old_vertices;

	vtkSurface *ptr=this->SubdivisionInput;
	vtkSurface *ptr2=this->SubdivisionOutput;

	vtkIdType i,v1,v2,v3;

	number_of_old_vertices=ptr->GetNumberOfPoints();
	number_of_vertices=ptr2->GetNumberOfPoints();

	vtkSparseMatrix::NonZeroElement *NonZero1,*NonZero2;

	double factor=1;

	double p[3],Distance,Nearest=1000000,xmin,ymin,xmax,ymax,x,y;

	v1=0;

	xmin=10000;
	xmax=-10000;
	ymin=10000;
	ymax=-10000;

	for (i=0;i<number_of_vertices;i++)
	{
		ptr2->GetPoints()->GetPoint(i,p);
		if (p[0]<xmin) xmin=p[0];
		if (p[0]>xmax) xmax=p[0];
		if (p[1]<ymin) ymin=p[1];
		if (p[1]>ymax) ymax=p[1];
	}

	y=0.5*(ymin+ymax);
	x=0.5*(xmin+xmax);

	for (i=number_of_old_vertices;i<number_of_vertices;i++)
	{
		ptr2->GetPoints()->GetPoint(i,p);
		Distance=(p[0]-x)*(p[0]-x)+(p[1]-y)*(p[1]-y);
		if ((Distance<Nearest)&&(p[2]>-0.1))
		{
			Nearest=Distance;
			v1=i;
		}
	}

	for( v3=0;v3<ptr2->GetNumberOfPoints();v3++)
	{
		ptr2->GetPoints()->GetPoint(v3,p);
		if (p[2]>-0.2)
		{
			p[2]=0;
			ptr2->GetPoints()->SetPoint(v3,p);
		}
	}

	ptr2->GetPoints()->GetPoint(v1,p);
	p[2]=factor;
	ptr2->GetPoints()->SetPoint(v1,p);

	if (this->ProcessLifting==1)
	{
		NonZero1=this->Alpha->FirstVer[v1-number_of_old_vertices];		

		while (NonZero1!=0)
		{
			v2=NonZero1->i;
			NonZero2=this->P->FirstVer[v2];
			while (NonZero2!=0)
			{
				v3=NonZero2->i;
				ptr2->GetPoints()->GetPoint(v3,p);
				p[2]-=NonZero1->Value*NonZero2->Value*factor;
				NonZero2=NonZero2->NextVer;
				ptr2->GetPoints()->SetPoint(v3,p);	
			}
			NonZero1=NonZero1->NextVer;
		}
	}

	this->Modified();
	ptr2->Modified();
}

void vtkWaveletSubdivisionFilter::DisplayScalingFunction()
{
	int number_of_vertices,number_of_old_vertices;

	vtkSurface *ptr=this->SubdivisionInput;
	vtkSurface *ptr2=this->SubdivisionOutput;

	vtkIdType i,v1;

	number_of_old_vertices=ptr->GetNumberOfPoints();
	number_of_vertices=ptr2->GetNumberOfPoints();

	double factor=1;

	double p[3],Distance,Nearest=1000000,xmin,ymin,xmax,ymax,x,y;

	v1=0;

	xmin=10000;
	xmax=-10000;
	ymin=10000;
	ymax=-10000;

	for (i=0;i<number_of_vertices;i++)
	{
		ptr2->GetPoints()->GetPoint(i,p);
		if (p[0]<xmin) xmin=p[0];
		if (p[0]>xmax) xmax=p[0];
		if (p[1]<ymin) ymin=p[1];
		if (p[1]>ymax) ymax=p[1];
	}

	y=0.5*(ymin+ymax);
	x=0.5*(xmin+xmax);

	for (i=0;i<number_of_old_vertices;i++)
	{
		ptr2->GetPoints()->GetPoint(i,p);
		Distance=(p[0]-x)*(p[0]-x)+(p[1]-y)*(p[1]-y);
		if ((Distance<Nearest)&&(p[2]>-0.1))
		{
			Nearest=Distance;
			v1=i;
		}
	}

	ptr2->GetPoints()->GetPoint(v1,p);
	p[2]=factor;
	ptr2->GetPoints()->SetPoint(v1,p);
	ptr->GetPoints()->SetPoint(v1,p);

	this->Modified();
	ptr2->Modified();

}

void vtkWaveletSubdivisionFilter::ClearWavelets()
{
	int i;
	double w[3];
	w[0]=0;
	w[1]=0;
	w[2]=0;

	for (i=0;i<this->Wavelets->GetNumberOfTuples();i++)
		this->Wavelets->SetTuple(i,w);
}

void vtkWaveletSubdivisionFilter::MakePMatrix(vtkSurface *ptr,vtkSurface *ptr2)
{
	int number_of_vertices,number_of_old_vertices;
	vtkIdType t,v1,v2;

	number_of_old_vertices=ptr->GetNumberOfPoints();
	number_of_vertices=ptr2->GetNumberOfPoints();

	this->P->Init(number_of_vertices,number_of_old_vertices);

	for (t=0;t<number_of_old_vertices;t++)
		this->P->AddValue(t,t,1);

	for (t=number_of_old_vertices;t<number_of_vertices;t++)
	{
		v1=this->vertices[t].parent1;
		v2=this->vertices[t].parent2;
		this->P->AddValue(t,v1,0.5);
		this->P->AddValue(t,v2,0.5);
	}
}

void vtkWaveletSubdivisionFilter::MakeInMatrixWeighted(vtkSurface *ptr,vtkSparseMatrix *In)
{
	int i,number_of_vertices,count2;
	vtkIdType t,edge,f1,f2,v1,v2;
	vtkIdList *Edges=vtkIdList::New();

	number_of_vertices=ptr->GetNumberOfPoints();

	In->Init(number_of_vertices,number_of_vertices);

	for (t=0;t<number_of_vertices;t++)
	{
		count2=0;

		ptr->GetVertexNeighbourEdges(t,Edges);

		for (i=0;i<Edges->GetNumberOfIds();i++)
		{
			edge=Edges->GetId(i);
			ptr->GetEdgeFaces(edge,f1,f2);
			f1=this->faces[f1].lowrestype;
			if (f2<0)
				f2=0;
			else
				f2=this->faces[f2].lowrestype;

			ptr->GetEdgeVertices(edge,v1,v2);
			if (v1==t) v1=v2;

			count2+=f1+f2;
			In->AddValue(v1,t,(double)(f1+f2)/24);
		}
		In->AddValue(t,t,(double)count2/24);
	}

	Edges->Delete();
}

void vtkWaveletSubdivisionFilter::MakeInMatrixUnweighted(vtkSurface *ptr, vtkSparseMatrix *In)
{
	int i,number_of_vertices,count;
	vtkIdType t,v1,v2,f1,f2;
	vtkIdList *Edges=vtkIdList::New();

	number_of_vertices=ptr->GetNumberOfPoints();

	In->Init(number_of_vertices,number_of_vertices);

	for (t=0;t<number_of_vertices;t++)
	{
		count=0;

		ptr->GetVertexNeighbourEdges(t,Edges);

		for (i=0;i<Edges->GetNumberOfIds();i++)
		{
			ptr->GetEdgeVertices(Edges->GetId(i),v1,v2);

			if (v1==t) v1=v2;

			ptr->GetEdgeFaces(Edges->GetId(i),f1,f2);

			if (f2<0)
			{
				In->AddValue(v1,t,(double)1/24);
				count++;
			}
			else
			{
				In->AddValue(v1,t,(double)1/12);
				count+=2;
			}
		}
		In->AddValue(t,t,(double)count/24);
	}
	Edges->Delete();
}

void vtkWaveletSubdivisionFilter::Orthogonalize ()
{
	if (this->ProcessLifting!=1)
		return;

	vtkSurface *ptr1=this->SubdivisionInput;
	vtkSurface *ptr2=this->SubdivisionOutput;
	int radius=this->LiftingRadius;

	vtkSparseMatrix *In1,*In2;
	vtkSparseMatrix::NonZeroElement *NextNonZero1,*NextNonZero2;

	vtkIdType i,j,k,PtId,PtId2,nPtId;
	vtkIdList *HorizontalList,*VerticalList,*TempList,*SortList;
	int offset;

	double *tuple,*tuple2,sum,valeur;
	vtkDoubleArray *Matrix, *Matrix2,*Second,*Diag,*Vect2,*Sol;
	int Height,Width,MaxHeight=100,MaxWidth=100,update;

	Matrix=vtkDoubleArray::New();
	Matrix->SetNumberOfComponents(MaxWidth);
	Matrix->SetNumberOfTuples(MaxHeight);

	Second=vtkDoubleArray::New();
	Second->SetNumberOfValues(MaxHeight);

	Vect2=vtkDoubleArray::New();
	Vect2->SetNumberOfValues(MaxWidth);

	Matrix2=vtkDoubleArray::New();
	Matrix2->SetNumberOfComponents(MaxWidth);
	Matrix2->SetNumberOfTuples(MaxWidth);

	Diag=vtkDoubleArray::New();
	Diag->SetNumberOfValues(MaxWidth);

	Sol=vtkDoubleArray::New();
	Sol->SetNumberOfValues(MaxWidth);

	update=0;

	In1=vtkSparseMatrix::New();
	In2=vtkSparseMatrix::New();

	this->MakeInMatrixWeighted(ptr1,In1);
	this->MakeInMatrixUnweighted(ptr2,In2);
	this->MakePMatrix(ptr1,ptr2);
	this->Alpha->Init(ptr1->GetNumberOfPoints(),ptr2->GetNumberOfPoints()-ptr1->GetNumberOfPoints());

	HorizontalList=vtkIdList::New();
	VerticalList=vtkIdList::New();
	SortList=vtkIdList::New();

	double Ex=0,Ex2=0,n=0;

	for (PtId=ptr1->GetNumberOfPoints();PtId<ptr2->GetNumberOfPoints();PtId++)
	{

		HorizontalList->SetNumberOfIds(2);
		HorizontalList->SetId(0,this->vertices[PtId].parent1);
		HorizontalList->SetId(1,this->vertices[PtId].parent2);

		if (radius>0)
		{
			for (i=0;i<radius;i++)
			{
				ptr1->GetNeighbours(HorizontalList,VerticalList);
				TempList=VerticalList;
				VerticalList=HorizontalList;
				HorizontalList=TempList;
			}
		}
		ptr1->GetNeighbours(HorizontalList,VerticalList);


		// Tri des Indices des sommets de la liste VerticalList
		SortList->SetNumberOfIds(0);
		Height=VerticalList->GetNumberOfIds();

		for (i=0;i<Height;i++)
		{
			nPtId=-1;
			k=-1;
			for (j=0;j<Height;j++)
			{
				PtId2=VerticalList->GetId(j);
				if (PtId2>nPtId)
				{
					nPtId=PtId2;
					k=j;
				}
			}
			SortList->InsertNextId(nPtId);
			VerticalList->SetId(k,-1);
		}

		VerticalList->SetNumberOfIds(Height);

		for (i=0;i<Height;i++)
			VerticalList->SetId(i,SortList->GetId(Height-i-1));

		// Tri des Indices des sommets de la liste HorizontalList
		SortList->SetNumberOfIds(0);
		Width=HorizontalList->GetNumberOfIds();

		for (i=0;i<Width;i++)
		{
			nPtId=-1;
			k=-1;
			for (j=0;j<Width;j++)
			{
				PtId2=HorizontalList->GetId(j);
				if (PtId2>nPtId)
				{
					nPtId=PtId2;
					k=j;
				}
			}
			SortList->InsertNextId(nPtId);
			HorizontalList->SetId(k,-1);
		}

		HorizontalList->SetNumberOfIds(Width);

		for (i=0;i<Width;i++)
			HorizontalList->SetId(i,SortList->GetId(Width-i-1));

		if (Height>MaxHeight)
		{
			MaxHeight=2*Height;
			update=1;
		}

		if (Width>MaxWidth)
		{
			MaxWidth=2*Width;
			update=1;
		}

		if (update==1)
		{
			Matrix->Delete();
			Second->Delete();
			Vect2->Delete();
			Matrix2->Delete();
			Diag->Delete();
			Sol->Delete();

			update=0;

			Matrix=vtkDoubleArray::New();
			Matrix->SetNumberOfComponents(MaxWidth);
			Matrix->SetNumberOfTuples(MaxHeight);

			Second=vtkDoubleArray::New();
			Second->SetNumberOfValues(MaxHeight);

			Vect2=vtkDoubleArray::New();
			Vect2->SetNumberOfValues(MaxWidth);


			Matrix2=vtkDoubleArray::New();
			Matrix2->SetNumberOfComponents(MaxWidth);
			Matrix2->SetNumberOfTuples(MaxWidth);

			Diag=vtkDoubleArray::New();
			Diag->SetNumberOfValues(MaxWidth);

			Sol=vtkDoubleArray::New();
			Sol->SetNumberOfValues(MaxWidth);

		}


		//calcul de la matrice
		tuple=Matrix->GetPointer(0);

		for (i=0;i<Height;i++,tuple+=MaxWidth)
		{
			NextNonZero1=In1->FirstHor[VerticalList->GetId(i)];
			j=0;

			while (j<Width)
			{
				if (NextNonZero1==0)
				{
					for (;j<Width;j++)
						tuple[j]=0;
				}
				else
				{
					PtId2=NextNonZero1->j;
					if (PtId2<HorizontalList->GetId(j))
						NextNonZero1=NextNonZero1->NextHor;
					else
					{
						if (PtId2>HorizontalList->GetId(j))
						{
							tuple[j]=0;
							j++;
						}
						else
						{
							tuple[j]=NextNonZero1->Value;
							j++;
							NextNonZero1=NextNonZero1->NextHor;
						}
					}
				}
			}
		}

		//calcul du second membre
		for (i=0;i<Height;i++)
		{
			PtId2=VerticalList->GetId(i);
			sum=0;

			NextNonZero1=In2->FirstVer[PtId];
			NextNonZero2=this->P->FirstVer[PtId2];
			while ((NextNonZero2!=0)&&(NextNonZero1!=0))
			{
				if (NextNonZero1->i>NextNonZero2->i)
				{
					NextNonZero2=NextNonZero2->NextVer;
				}
				else
				{
					if (NextNonZero1->i<NextNonZero2->i)
					{
						NextNonZero1=NextNonZero1->NextVer;
					}
					else
					{
						sum+=NextNonZero1->Value*NextNonZero2->Value;
						NextNonZero2=NextNonZero2->NextVer;
						NextNonZero1=NextNonZero1->NextVer;
					}
				}
			}
			Second->SetValue(i,sum);
		}

		// calcule tA*A  dans Matrix2
		tuple2=Matrix2->GetPointer(0);
		for (i=0;i<Width;i++,tuple2+=MaxWidth)
		{
			for (j=i;j<Width;j++)
			{
				sum=0.0;
				tuple=Matrix->GetPointer(0);
				for (k=0;k<Height;k++,tuple+=MaxWidth)
				{
					sum=sum+tuple[i]*tuple[j];
				}
				tuple2[j]=sum;
			}
		}

		//calcule tA*B dans Vect2 
		for (i=0;i<Width;i++)
		{
			sum=0;
			tuple=Matrix->GetPointer(0);
			for (j=0;j<Height;j++,tuple+=MaxWidth)
			{
				sum+=tuple[i]*Second->GetValue(j);
			}
			Vect2->SetValue(i,sum);
		}

		//	decomposition de choleski
		tuple=Matrix2->GetPointer(0);
		for (i=0;i<Width;i++,tuple+=MaxWidth)
		{
			tuple2=tuple;
			for (j=i;j<Width;j++,tuple2+=MaxWidth)
			{

				for (sum=tuple[j],k=i-1;k>=0;k--)
					sum-=tuple[k]*tuple2[k];
				if (i==j)
				{
					if (sum<=0.0)
						cout<<"erreur de zero : sum="<<sum<<endl;

					Diag->SetValue(i,sqrt(sum));
				}
				else
					tuple2[i]=sum/Diag->GetValue(i);
			}
		}

		// résolution de choleski:
		tuple=Matrix2->GetPointer(0);
		for (i=0;i<Width;i++,tuple+=MaxWidth)
		{

			for (sum=Vect2->GetValue(i),k=i-1;k>=0;k--)
				sum-=tuple[k]*Sol->GetValue(k);

			Sol->SetValue(i,valeur=sum/Diag->GetValue(i));
		}
		for (i=Width-1;i>=0;i--,tuple-=MaxWidth)
		{
			tuple2=tuple;
			for (sum=Sol->GetValue(i),k=i+1;k<Width;k++,tuple2+=MaxWidth)
				sum-=tuple2[i]*Sol->GetValue(k);

			Sol->SetValue(i,sum/Diag->GetValue(i));
		}

		//ajout à alpha

		offset=ptr1->GetNumberOfPoints();

		for(i=0;i<HorizontalList->GetNumberOfIds();i++)
		{
			Ex+=Sol->GetValue(i);
			Ex2+=Sol->GetValue(i)*Sol->GetValue(i);
			n++;
			this->Alpha->AddValue(HorizontalList->GetId(i),PtId-offset,Sol->GetValue(i));
		}
	}

	double moyenne,variance,ecart;

	moyenne=Ex/n;
	variance=(moyenne-(Ex2/n)*(Ex2/n));
	ecart=sqrt(variance);

	In1->Delete();
	In2->Delete();

	Matrix->Delete();
	Second->Delete();
	Vect2->Delete();
	Matrix2->Delete();
	Diag->Delete();
	Sol->Delete();

	HorizontalList->Delete();
	VerticalList->Delete();
	SortList->Delete();
}


int vtkWaveletSubdivisionFilter::IsFaceSwitch(vtkIdType v1,vtkIdType v2,vtkIdType v3)
{
	vtkIdType v4,v5,v6,v7,f1,f2;

	if (this->IOType==2)
	{
		v4=this->DecodeFaceSwap(1);
		if (v4==1)
			return (1);
		else
			return (-1);
	}

	v4=this->vertices[this->PointsIds->GetId(v1)].old_to_new;
	v5=this->vertices[this->PointsIds->GetId(v2)].old_to_new;
	v6=this->vertices[this->PointsIds->GetId(v3)].old_to_new;

	this->Bits1to2++;

	this->MergeOutput->GetEdgeFaces(v7=this->MergeOutput->IsEdge(v4,v5),f1,f2);

	v7=this->MergeOutput->GetThirdPoint(f1,v4,v5);
	if (v7!=v6)
		f1=f2;

	if (this->faces[f1].vertex>=0)
	{
		if (this->IOType==1)
			this->EncodeFaceSwap(1,1);
		return (1);
	}
	else
	{
		if (this->IOType==1)
			this->EncodeFaceSwap(0,1);
		return (-1);
	}
}

vtkIdType vtkWaveletSubdivisionFilter::GetEdgeSwap(vtkIdType v1,vtkIdType v2,vtkIdType v3)
{
	vtkIdType v4,v5,v6,v7,f1,f2;

	if (this->IOType==2)
	{
		v4=this->DecodeEdgeSwap(2);
		if (v4==1)
			return (v1);
		else
			return (v2);
	}

	v4=this->vertices[PointsIds->GetId(v1)].old_to_new;
	v5=this->vertices[PointsIds->GetId(v2)].old_to_new;
	v6=this->vertices[PointsIds->GetId(v3)].old_to_new;

	this->MergeOutput->GetEdgeFaces(this->MergeOutput->IsEdge(v4,v5),f1,f2);

	this->BitsSwap++;

	v7=this->MergeOutput->GetThirdPoint(f1,v4,v5);

	if (v7!=v6)
		f1=f2;

	if (this->faces[f1].vertex==this->vertices[v4].new_to_old)
	{
		if (this->IOType==1)
			this->EncodeEdgeSwap(1,2);
		return (v1);
	}
	else
	{
		if (this->IOType==1)
			this->EncodeEdgeSwap(0,2);
		return (v2);
	}
}
vtkIdType vtkWaveletSubdivisionFilter::GetEdgeMidPoint(vtkIdType v1,vtkIdType v2)
{
	vtkIdType e1,v3,v4,v1p,v2p;
	double p1[3],p2[3],p3[3];

	vtkSurface *Mesh=this->SubdivisionOutput;
	switch (this->SubdivisionType)
	{
	case 0:

		Mesh->GetPoint(v1,p1);
		Mesh->GetPoint(v2,p2);
		p3[0]=0.5*p1[0]+0.5*p2[0];
		p3[1]=0.5*p1[1]+0.5*p2[1];
		p3[2]=0.5*p1[2]+0.5*p2[2];
		v4=Mesh->AddVertex(p3);
		this->vertices[v4].parent1=v1;
		this->vertices[v4].parent2=v2;
		return (v4);

	case 1:
		switch(this->IOType)
		{
		case 0:
			this->BitsMidPoints++;
			v1p=this->vertices[this->PointsIds->GetId(v1)].old_to_new;
			v2p=this->vertices[this->PointsIds->GetId(v2)].old_to_new;
			e1=this->MergeOutput->IsEdge(v1p,v2p);
			if (e1<0)
				cout<<endl<<"probleme d'arete"<<endl;
			v3=this->edgesvector[e1].child;
			if (v3>=0)
			{
				this->BitsMidPoints1++;
				if (this->ConnectivityOnly==0)
				{
					this->MergeInput->GetPoints()->GetPoint(v3,p3);	
					v4=Mesh->AddVertex(p3);
				}
				else
				{
					Mesh->GetPoint(v1,p1);
					Mesh->GetPoint(v2,p2);
					p3[0]=0.5*p1[0]+0.5*p2[0];
					p3[1]=0.5*p1[1]+0.5*p2[1];
					p3[2]=0.5*p1[2]+0.5*p2[2];
					v4=Mesh->AddVertex(p3);
				}
				PointsIds->SetId(v4,v3);
				this->vertices[v4].parent1=v1;
				this->vertices[v4].parent2=v2;

				return (v4);
			}
			else
			{
				return (-1);
			}
		case 1:

			this->BitsMidPoints++;
			v1p=this->vertices[this->PointsIds->GetId(v1)].old_to_new;
			v2p=this->vertices[this->PointsIds->GetId(v2)].old_to_new;
			e1=this->MergeOutput->IsEdge(v1p,v2p);
			if (e1<0)
				cout<<endl<<"problème d'arête"<<endl;
			v3=this->edgesvector[e1].child;

			e1=this->SubdivisionInput->IsEdge(v1,v2);
			if (e1<0)
				cout<<endl<<"problème d'arête"<<endl;
			v4=this->EdgeMidPoints->GetId(e1);
			if (v4>=0)
			{
				this->BitsMidPoints1++;
				PointsIds->SetId(v4,v3);
				this->EncodeMidpoint(1);
				return (v4);
			}
			else
			{
				this->EncodeMidpoint(0);
				return (-1);
			}
		case 2:

			v3=this->DecodeMidpoint();

			if (v3>0)
			{
				Mesh->GetPoint(v1,p1);
				Mesh->GetPoint(v2,p2);
				p3[0]=0.5*p1[0]+0.5*p2[0];
				p3[1]=0.5*p1[1]+0.5*p2[1];
				p3[2]=0.5*p1[2]+0.5*p2[2];
				v4=Mesh->AddVertex(p3);
				this->vertices[v4].parent1=v1;
				this->vertices[v4].parent2=v2;
				return (v4);
			}
			else
			{
				return (-1);
			}
		}
	default:
		return (-1);

	}
}

void vtkWaveletSubdivisionFilter::Subdivide()
{
	this->Subdivision();
}

void vtkWaveletSubdivisionFilter::AddChildEdgeToTree(vtkIdType Child, vtkIdType Parent)
{
	if ((Child==-1)||(Parent==-1))
		cout<<"Probleme fils-pere!!!!"<<endl;

	if (this->TreeParentEdge->GetId(Child)==-1)
	{
		this->TreeParentEdge->SetId(Child,Parent);
		this->TreeNextChildEdge->SetId(Child,this->TreeFirstChildEdge->GetId(Parent));
		this->TreeFirstChildEdge->SetId(Parent,Child);
	}
	return;	
}

void vtkWaveletSubdivisionFilter::InitTree()
{
	int i;
	this->TreeFirstChildEdge->SetNumberOfIds(this->SubdivisionInput->GetNumberOfEdges());
	for (i=0;i<this->SubdivisionInput->GetNumberOfEdges();i++)
		this->TreeFirstChildEdge->SetId(i,-1);
	this->TreeNextChildEdge->SetNumberOfIds(this->SubdivisionInput->GetNumberOfEdges()*4);
	this->TreeParentEdge->SetNumberOfIds(this->SubdivisionInput->GetNumberOfEdges()*4);
	for (i=0;i<this->SubdivisionInput->GetNumberOfEdges()*4;i++)
	{
		this->TreeNextChildEdge->SetId(i,-1);
		this->TreeParentEdge->SetId(i,-1);
	}
}

void vtkWaveletSubdivisionFilter::Subdivision()
{
	vtkIdType i,e1,e2,e3,v1,v2,v3,v4,v5,v6,f1,f2,vp1,vp2,vp3,edge1;
	if ((this->SubdivisionType==0)&&(this->IOType==1))
		return;

	vtkSurface *Out=this->SubdivisionOutput;
	vtkSurface *Input=this->SubdivisionInput;

	vtkBitArray *FacesSwapPerformed=vtkBitArray::New();
	vtkBitArray *EdgesSwap=vtkBitArray::New();

	vtkIdList *EdgeMidPoints2=vtkIdList::New();

	int LeftCheck;
	int RightCheck;

	this->BitsMidPoints=0;
	this->BitsMidPoints1=0;
	this->Bits1to3=0;
	this->Bits1to2=0;
	this->BitsSwap=0;

	if (this->IOType!=1)
	{
		this->InitTree();

		EdgeMidPoints->SetNumberOfIds(Input->GetNumberOfEdges());
		Out->Init(Input->GetNumberOfPoints()+Input->GetNumberOfEdges()
			,Input->GetNumberOfCells()*4
			,Input->GetNumberOfEdges()*2+Input->GetNumberOfCells()*3);
		for (i=0;i<Input->GetNumberOfPoints();i++)
		{
			Out->AddVertex(Input->GetPoint(i));
		}
	}
	EdgeMidPoints2->SetNumberOfIds(Input->GetNumberOfEdges());

	if (this->SubdivisionType==0)
	{
		this->vertices.resize(Input->GetNumberOfPoints()+Input->GetNumberOfEdges());
		this->faces.resize(Input->GetNumberOfCells()*4);
		for (i=0;i<Input->GetNumberOfEdges();i++)
		{
			EdgeMidPoints2->SetId(i,-2);
		}
	}
	else
	{
		if (this->IOType<2)
		{
			for (i=0;i<Input->GetNumberOfPoints();i++)
			{
				v1=PointsIds->GetId(i);
				vp1=this->vertices[v1].new_to_old;	
				PointsIds->SetId(i,this->vertices[v1].new_to_old);
			}
		}
		else
		{
			this->vertices.resize(Input->GetNumberOfPoints()+Input->GetNumberOfEdges());
			this->faces.resize(Input->GetNumberOfCells()*4);
		}

		EdgesSwap->Allocate(Input->GetNumberOfEdges());
		FacesSwapPerformed->Allocate(Input->GetNumberOfCells());
		for (i=0;i<Input->GetNumberOfCells();i++)
			FacesSwapPerformed->SetValue(i,0);
		for (i=0;i<Input->GetNumberOfEdges();i++)
		{
			EdgesSwap->SetValue(i,0);
			EdgeMidPoints2->SetId(i,-2);
		}
	}

	vtkIdType Vertices[4];
	int j;
	for (i=0;i<Input->GetNumberOfCells();i++)
	{
		Input->GetFaceVertices(i,v1,v2,v3);
		Input->GetFaceVertices(i,Vertices[0],Vertices[1],Vertices[2]);
		Vertices[3]=Vertices[0];
		for (j=0;j<3;j++)
		{
			e1=Input->IsEdge(Vertices[j],Vertices[j+1]);
			if (EdgeMidPoints2->GetId(e1)<-1)
			{
				EdgeMidPoints2->SetId(e1,1);
				EdgeMidPoints->SetId(e1,this->GetEdgeMidPoint(Vertices[j],Vertices[j+1]));
			}
		}
	}

	for (i=0;i<Input->GetNumberOfCells();i++)
	{
		Input->GetFaceVertices(i,v1,v2,v3);

		e1=Input->IsEdge(v1,v2);
		e2=Input->IsEdge(v1,v3);
		e3=Input->IsEdge(v2,v3);

		v4=EdgeMidPoints->GetId(e1);
		v5=EdgeMidPoints->GetId(e2);
		v6=EdgeMidPoints->GetId(e3);

		if (v4<0)
		{
			if (v5<0)
			{
				if (v6<0)
				{
					this->faces[i].lowrestype=1;
				}
				else
				{
					this->faces[i].lowrestype=2;
				}
			}
			else
			{
				if (v6<0)
				{
					this->faces[i].lowrestype=2;
				}
				else
				{
					this->faces[i].lowrestype=3;

				}
			}
		}
		else
		{
			if (v5<0)
			{
				if (v6<0)
				{
					this->faces[i].lowrestype=2;
				}
				else
				{
					this->faces[i].lowrestype=3;
				}
			}
			else
			{
				if (v6<0)
				{
					this->faces[i].lowrestype=3;
				}
				else
				{

					this->faces[i].lowrestype=4;
				}
			}
		}

		if (this->SubdivisionType==1)
		{
			if ((this->faces[i].lowrestype==3)||(this->faces[i].lowrestype==4))
				FacesSwapPerformed->SetValue(i,1);
			else
				FacesSwapPerformed->SetValue(i,0);
		}
	}

	for (i=0;i<Input->GetNumberOfCells();i++)
	{
		Input->GetFaceVertices(i,v1,v2,v3);

		e1=Input->IsEdge(v1,v2);
		e2=Input->IsEdge(v1,v3);
		e3=Input->IsEdge(v2,v3);

		v4=EdgeMidPoints->GetId(e1);
		v5=EdgeMidPoints->GetId(e2);
		v6=EdgeMidPoints->GetId(e3);

		if (v4<0)
		{
			if (v5<0)
			{
				if (v6<0)
				{
				}
				else
				{
					v4=v6;
					v5=v1;
					v1=v2;
					v2=v3;
					v3=v5;
				}
			}
			else
			{
				if (v6<0)
				{
					v4=v5;
					v5=v2;
					v2=v1;
					v1=v3;
					v3=v5;
				}
				else
				{

				}
			}
		}
		else
		{
			if (v5<0)
			{
				if (v6<0)
				{
				}
				else
				{
					v5=v1;
					v1=v3;
					v3=v2;
					v2=v5;
					v5=v6;
					v6=v4;
				}
			}
			else
			{
				if (v6<0)
				{

					v6=v5;
					v5=v4;
					v4=v1;
					v1=v2;
					v2=v3;
					v3=v4;
				}
				else
				{
				}
			}
		}

		vtkIdType Edge1=-1,Edge2=-1,Edge3;
		switch(this->faces[i].lowrestype)
		{
		case 1:
			if (this->IOType!=1)
			{
				Out->AddFace(v1,v2,v3);
				this->AddChildEdgeToTree(Out->IsEdge(v1,v2),
					Input->IsEdge(v1,v2));
				this->AddChildEdgeToTree(Out->IsEdge(v1,v3),
					Input->IsEdge(v1,v3));
				this->AddChildEdgeToTree(Out->IsEdge(v3,v2),
					Input->IsEdge(v3,v2));
			}
			break;

		case 2:
			if (this->IOType!=1)
			{
				Out->AddFace(v1,v4,v3);
				Out->AddFace(v4,v2,v3);
				Edge1=Input->IsEdge(v1,v2);
				this->AddChildEdgeToTree(Out->IsEdge(v1,v4),Edge1);
				this->AddChildEdgeToTree(Out->IsEdge(v4,v2),Edge1);
				this->AddChildEdgeToTree(Out->IsEdge(v4,v3),Edge1);
				this->AddChildEdgeToTree(Out->IsEdge(v1,v3),
					Input->IsEdge(v1,v3));
				this->AddChildEdgeToTree(Out->IsEdge(v3,v2),
					Input->IsEdge(v3,v2));
			}
			if (FacesSwapPerformed->GetValue(i)==0)
			{

				Input->Conquer(i,v1,v3,f1,v5);
				Input->Conquer(i,v2,v3,f2,v6);

				if (f1<0)
				{
					LeftCheck=0;
					f1=i;
				}
				else
				{
					if ((FacesSwapPerformed->GetValue(f1)==0)
						&&(EdgeMidPoints->GetId(Input->IsEdge(v3,v5))<0))
						LeftCheck=1;
					else
						LeftCheck=0;
				}

				if (f2<0)
				{
					f2=i;
					RightCheck=0;
				}
				else
				{
					if ((FacesSwapPerformed->GetValue(f2)==0)
						&&(EdgeMidPoints->GetId(Input->IsEdge(v3,v6))<0))
						RightCheck=1;
					else
						RightCheck=0;
				}

				if (LeftCheck==1)
				{
					if (RightCheck==1)
					{
						if (this->IsFaceSwitch(v1,v2,v3)==1)
						{

							v5=this->GetEdgeSwap(v1,v2,v3);
							if (EdgesSwap->GetValue(Input->IsEdge(v3,v5))==1)
								cout<<"double permutation!!"<<endl;
							EdgesSwap->SetValue(Input->IsEdge(v3,v5),1);

							if (v5==v1)
							{
								FacesSwapPerformed->SetValue(i,1);
								FacesSwapPerformed->SetValue(f1,1);
							}
							else
							{
								FacesSwapPerformed->SetValue(i,1);
								FacesSwapPerformed->SetValue(f2,1);
							}
						}												
					}
					else
					{
						if (this->IsFaceSwitch(v1,v2,v3)==1)
						{
							if (EdgesSwap->GetValue(Input->IsEdge(v3,v1))==1)
								cout<<"double permutation!!"<<endl;
							EdgesSwap->SetValue(Input->IsEdge(v3,v1),1);
							FacesSwapPerformed->SetValue(i,1);
							FacesSwapPerformed->SetValue(f1,1);
						}
					}
				}
				else
				{
					if (RightCheck==1)
					{
						if (this->IsFaceSwitch(v1,v2,v3)==1)
						{
							if (EdgesSwap->GetValue(Input->IsEdge(v3,v2))==1)
								cout<<"double permutation!!"<<endl;
							EdgesSwap->SetValue(Input->IsEdge(v3,v2),1);
							FacesSwapPerformed->SetValue(i,1);
							FacesSwapPerformed->SetValue(f2,1);
						}
					}
				}

			}
			break;

		case 3:
			if (this->IOType!=1)
			{
				Edge1=Input->IsEdge(v1,v3);
				Edge2=Input->IsEdge(v2,v3);

				Out->AddFace(v3,v5,v6);
				this->AddChildEdgeToTree(Out->IsEdge(v3,v5),Edge1);
				this->AddChildEdgeToTree(Out->IsEdge(v3,v6),Edge2);
			}

			if (IOType==2)
			{
				if (this->DecodeFaceSwap(4)==0)
					v4=v5;
				else
					v4=v6;
			}
			else
			{
				vp1=this->vertices[PointsIds->GetId(v1)].old_to_new;
				vp2=this->vertices[PointsIds->GetId(v2)].old_to_new;
				vp3=this->vertices[PointsIds->GetId(v3)].old_to_new;
				edge1=this->MergeOutput->IsEdge(vp1,vp2);
				this->MergeOutput->GetEdgeFaces(edge1,f1,f2);
				if (f2>=0)
				{
					if (this->MergeOutput->GetThirdPoint(f2,vp1,vp2)==vp3)
						f1=f2;
				}
				v4=this->faces[f1].vertex;

				if (v4==PointsIds->GetId(v5))
					v4=v5;
				else
					v4=v6;	

				this->Bits1to3++;
				if (this->IOType==1)
				{
					if (v4==v5)
						this->EncodeFaceSwap(0,4);
					else
						this->EncodeFaceSwap(1,4);
				}
			}

			if (this->IOType!=1)
			{
				if (v4==v5)
				{
					Out->AddFace(v1,v2,v5);
					Out->AddFace(v6,v5,v2);
					this->AddChildEdgeToTree(Out->IsEdge(v6,v2),Edge2);
					this->AddChildEdgeToTree(Out->IsEdge(v1,v5),Edge1);
					this->AddChildEdgeToTree(Out->IsEdge(v1,v2),
						Input->IsEdge(v1,v2));
					this->AddChildEdgeToTree(Out->IsEdge(v2,v5),Edge1);
					this->AddChildEdgeToTree(Out->IsEdge(v5,v6),Edge2);
				}
				else
				{
					Out->AddFace(v6,v5,v1);
					Out->AddFace(v2,v6,v1);

					this->AddChildEdgeToTree(Out->IsEdge(v6,v2),Edge2);
					this->AddChildEdgeToTree(Out->IsEdge(v1,v5),Edge1);
					this->AddChildEdgeToTree(Out->IsEdge(v1,v2),
						Input->IsEdge(v1,v2));
					this->AddChildEdgeToTree(Out->IsEdge(v6,v1),Edge2);
					this->AddChildEdgeToTree(Out->IsEdge(v5,v6),Edge1);
				}
			}
			break;

		case 4:
			if (this->IOType!=1)
			{
				Edge1=Input->IsEdge(v1,v2);
				Edge2=Input->IsEdge(v1,v3);
				Edge3=Input->IsEdge(v2,v3);
				Out->AddFace(v1,v4,v5);
				Out->AddFace(v6,v5,v4);
				Out->AddFace(v2,v6,v4);
				Out->AddFace(v3,v5,v6);
				this->AddChildEdgeToTree(Out->IsEdge(v1,v4),Edge1);
				this->AddChildEdgeToTree(Out->IsEdge(v2,v4),Edge1);
				this->AddChildEdgeToTree(Out->IsEdge(v4,v6),Edge1);

				this->AddChildEdgeToTree(Out->IsEdge(v1,v5),Edge2);
				this->AddChildEdgeToTree(Out->IsEdge(v3,v5),Edge2);
				this->AddChildEdgeToTree(Out->IsEdge(v4,v5),Edge2);

				this->AddChildEdgeToTree(Out->IsEdge(v2,v6),Edge3);
				this->AddChildEdgeToTree(Out->IsEdge(v3,v6),Edge3);
				this->AddChildEdgeToTree(Out->IsEdge(v5,v6),Edge3);
			}

		default:
			break;
		}
	}

	if ((this->ProcessLifting==1)&&(this->IOType!=1))
	{
		this->Orthogonalize();
	}

	if ((this->SubdivisionType==1)&&(this->IOType!=1))
	{
		for (i=0;i<Input->GetNumberOfEdges();i++)
		{
			if (EdgesSwap->GetValue(i)==1)
			{
				Input->GetEdgeVertices(i,v1,v2);
				Out->FlipEdge(Out->IsEdge(v1,v2));
			}
		}
	}

	EdgeMidPoints2->Delete();
	FacesSwapPerformed->Delete();
	EdgesSwap->Delete();

	Out->Squeeze();
}

void vtkWaveletSubdivisionFilter::ConquerEdge(vtkIdType v1,vtkIdType v2,vtkIdType & face,vtkIdType &vertex)
{
	vtkIdType f1;
	vtkIdType edge=this->MergeInput->IsEdge(v1,v2);

	this->MergeInput->GetEdgeFaces(edge,face,f1);

	if (faces[face].hirestype==1)
	{
		face=f1;
	}

	vertex=this->MergeInput->GetThirdPoint(face,v1,v2);

	return;
}

int vtkWaveletSubdivisionFilter::AddFace(vtkIdType v1,vtkIdType v2,vtkIdType v3,vtkIdType& edge1,vtkIdType& edge2,vtkIdType& edge3)
{
	vtkFloatingPointType point[3];
	vtkIdType face,nv1,nv2,nv3;

	this->MergeInput->GetPoints()->GetPoint(v1,point);
	nv1=vertices[v1].old_to_new;
	if (nv1==-1)
	{
		nv1=this->MergeOutput->AddVertex(point);
		vertices[v1].old_to_new=nv1;
		vertices[nv1].new_to_old=v1;
	}

	this->MergeInput->GetPoints()->GetPoint(v2,point);
	nv2=vertices[v2].old_to_new;
	if (nv2==-1)
	{
		nv2=this->MergeOutput->AddVertex(point);
		vertices[v2].old_to_new=nv2;
		vertices[nv2].new_to_old=v2;
	}

	this->MergeInput->GetPoints()->GetPoint(v3,point);
	nv3=vertices[v3].old_to_new;
	if (nv3==-1)
	{
		nv3=this->MergeOutput->AddVertex(point);
		vertices[v3].old_to_new=nv3;
		vertices[nv3].new_to_old=v3;
	}

	face=this->MergeOutput->AddFace(nv1,nv2,nv3);
	edge1=this->MergeOutput->IsEdge(nv1,nv2);
	edge2=this->MergeOutput->IsEdge(nv1,nv3);
	edge3=this->MergeOutput->IsEdge(nv2,nv3);

	return(face);
}

void vtkWaveletSubdivisionFilter::switchcontour(vtkIdType edge)
{
	vtkIdType v1,v2,v3,f1,f2,edge2;
	contour_element *next,*previous;

	this->MergeOutput->GetEdgeVertices(edge,v1,v2);
	v1=this->vertices[v1].new_to_old;
	v2=this->vertices[v2].new_to_old;

	if ((v3=edgesvector[edge].child)<0)
	{
		edge2=this->MergeInput->IsEdge(v1,v2);
		this->MergeInput->GetEdgeFaces(edge2,f1,f2);
		if (f2<0) return;
	}
	else
	{
		edge2=this->MergeInput->IsEdge(v1,v3);
		this->MergeInput->GetEdgeFaces(edge2,f1,f2);
		if (f2<0) return;
	}

	if (edgesvector[edge].contour==0)
	{
		if (v3<0)
		{
			edge2=this->MergeInput->IsEdge(v1,v2);
			this->MergeInput->GetEdgeFaces(edge2,f1,f2);
			if ((this->faces[f1].hirestype>=0)&&(this->faces[f2].hirestype>=0)) return;
		}
		else
		{
			edge2=this->MergeInput->IsEdge(v1,v3);
			this->MergeInput->GetEdgeFaces(edge2,f1,f2);
			if ((this->faces[f1].hirestype>=0)&&(this->faces[f2].hirestype>=0)) return;
		}

		edgesvector[edge].contour=new contour_element;
		edgesvector[edge].contour->edge=edge;

		if (LastElement==0)
		{
			FirstElement=edgesvector[edge].contour;
			LastElement=edgesvector[edge].contour;
			edgesvector[edge].contour->NextElement=0;
			edgesvector[edge].contour->PreviousElement=0;
		}
		else
		{
			edgesvector[edge].contour->NextElement=0;
			edgesvector[edge].contour->PreviousElement=LastElement;
			LastElement->NextElement=edgesvector[edge].contour;
			LastElement=edgesvector[edge].contour;
		}
	}
	else
	{
		next=edgesvector[edge].contour->NextElement;
		previous=edgesvector[edge].contour->PreviousElement;

		if ((next==0)&&(previous==0))
		{
			LastElement=0;
			FirstElement=0;
			delete edgesvector[edge].contour;
			edgesvector[edge].contour=0;
			return;
		}
		else
		{
			if (next==0)
			{
				previous->NextElement=0;
				LastElement=previous;
				delete edgesvector[edge].contour;
				edgesvector[edge].contour=0;
				return;
			}
			else
			{
				if (previous==0)
				{
					FirstElement=next;
					next->PreviousElement=0;
					delete edgesvector[edge].contour;
					edgesvector[edge].contour=0;
					return;
				}
				else
				{
					next->PreviousElement=previous;
					previous->NextElement=next;
					delete edgesvector[edge].contour;
					edgesvector[edge].contour=0;
					return;		
				}
			}
		}
	}
}

int vtkWaveletSubdivisionFilter::MatchParents(vtkIdType v1,vtkIdType p1,vtkIdType p2,vtkIdType f1,
											  vtkIdType f2,int incidentfaces)
{
	int Test,i;
	double s1[3],s2[3],s3[3],t1[3],t2[3],t3[3],r;

	Test=MatchParents2(v1,p1,p2,f1,f2,incidentfaces);

	if (Test==0)
		return (0);
	else
	{
		if ((this->GeometryCriterion==1)&&(this->MergeInput->GetNumberOfPoints()>200))
		{
			vtkPoints *points=this->MergeInput->GetPoints();

			if (this->MergeInput->IsSharpVertex(v1)==1)
			{
				points->GetPoint(v1,s1);
				points->GetPoint(p1,s2);
				points->GetPoint(p2,s3);
				for (i=0;i<3;i++)
				{
					t1[i]=s1[i]-0.5*(s2[i]+s3[i]);
					t2[i]=s2[i]-s3[i];
				}
				t3[0]=t1[1]*t2[2]-t1[2]*t2[1];
				t3[1]=t1[2]*t2[0]-t1[0]*t2[2];
				t3[2]=t1[0]*t2[1]-t1[1]*t2[0];	

				r=sqrt(t3[1]*t3[1]+t3[2]*t3[2]+t3[0]*t3[0])
					/(t2[1]*t2[1]+t2[2]*t2[2]+t2[0]*t2[0]);

				if (r<this->WaveletTreshold)
					return (1);
				else
					return (0);
			}
			/*

			if (this->MergeInput->IsSharpVertex(v1)==1)
			{
			double Normal[3];
			vtkDoubleArray *Areas=this->MergeInput->GetTrianglesAreas();
			vtkDoubleArray *Normals=this->MergeInput->GetTrianglesNormals();
			int Face;
			vtkIdList *Faces=vtkIdList::New();
			this->MergeInput->GetVertexNeighbourFaces(v1,Faces);
			double SVolume=0,SArea=0,Area;

			points->GetPoint(v1,s1);
			points->GetPoint(p1,s2);
			points->GetPoint(p2,s3);
			for (i=0;i<3;i++)
			{
			t1[i]=s1[i]-0.5*(s2[i]+s3[i]);
			}

			for (i=0;i<Faces->GetNumberOfIds();i++)
			{
			Face=Faces->GetId(i);
			Area=Areas->GetValue(Face);
			SArea+=Area;
			Normals->GetTuple(Face,Normal);
			SVolume+=fabs(t1[0]*Normal[0]+t1[1]*Normal[1]+t1[2]*Normal[2])*Area/3.0;
			}
			Faces->Delete();
			if (pow(SVolume,0.33333)/(sqrt(SArea))>this->WaveletTreshold)
			return (0);
			else
			return (1);
			}
			*/
			else
				return (1);


		}
		else
			return (1);
	}
}

int vtkWaveletSubdivisionFilter::MatchParents2(vtkIdType v1,vtkIdType p1,vtkIdType p2,vtkIdType f1,
											   vtkIdType f2,int incidentfaces)
{
	vtkIdType edge1,nv1,np1,np2,f3,f4,f5,f6,v2,v3,v4,v5,v6;

	if (this->vertices[v1].type==PARENT) return (0);

	if (p1==p2) return (0);

	if (this->MergeInput->GetNumberOfBoundaries(v1)>0)
	{
		edge1=this->MergeInput->IsEdge(v1,p1);
		this->MergeInput->GetEdgeFaces(edge1,f3,f4);
		if (f4>=0) return (0);
		edge1=this->MergeInput->IsEdge(v1,p2);
		this->MergeInput->GetEdgeFaces(edge1,f3,f4);
		if (f4>=0) return (0);
	}	

	if (this->vertices[v1].type==CHILD)
	{
		if (((this->vertices[v1].parent1!=p1)&&(this->vertices[v1].parent1!=p2))||
			((this->vertices[v1].parent2!=p1)&&(this->vertices[v1].parent2!=p2))) 
			return (0);
		else
			return (1);
	}

	if ((this->vertices[v1].valence-incidentfaces>3)||
		(this->vertices[v1].valence-incidentfaces<0)||
		(this->vertices[v1].valence-incidentfaces==1))
	{
		return 0;
	}

	edge1=this->MergeInput->IsEdge(p1,p2);
	if (edge1>=0) return (0);

	np1=this->vertices[p1].old_to_new;
	np2=this->vertices[p2].old_to_new;

	if ((np1>=0)&&(np2>=0))
	{
		edge1=this->MergeOutput->IsEdge(np1,np2);
		if (edge1>=0)
		{
			return(0);
			nv1=this->edgesvector[edge1].child;
			if ((nv1<0)||((nv1>=0)&&(nv1!=v1))) return (0);
		}
	}

	this->MergeInput->Conquer(f1,v1,p1,f3,v2);
	this->MergeInput->Conquer(f2,v1,p2,f4,v3);

	if ((f3==-1)&&(f4==-1))
	{
		return (1);
	}
	else if ((f3==-1)||(f4==-1))
	{
		return (0);
	}


	if (v2==v3)
	{
		return (1);
		if (this->PossibleParent(v2)==1)
			return (1);
		else
			return (0);
	}

	this->MergeInput->Conquer(f3,v1,v2,f5,v3);

	if ((this->vertices[v2].type==CHILD)&&(this->vertices[v3].type==CHILD))
	{
		return (1);
	}

	if (this->vertices[v2].type==CHILD)
	{
		if ((p1!=this->vertices[v2].parent1)&&(p1!=this->vertices[v2].parent2))
			return (0);

		if ((v3!=this->vertices[v2].parent1)&&(v3!=this->vertices[v2].parent2))
		{
			this->MergeInput->Conquer(f5,v2,v3,f6,v4);

			if (this->MergeInput->IsEdge(p1,v3)>=0) 
			{
				if (this->MatchParents(v3,p2,v4,f4,f6,3)!=1)
					return (0);
			}
			if ((this->vertices[p1].old_to_new>=0)&&(this->vertices[v3].old_to_new>=0))
			{
				if (this->MergeOutput->IsEdge(this->vertices[p1].old_to_new,this->vertices[v3].old_to_new)>=0)
				{
					if (this->MatchParents(v3,p2,v4,f4,f6,3)!=1)
						return (0);
				}
			}
		}
	}

	if (this->vertices[v3].type==CHILD)
	{
		if ((p2!=this->vertices[v3].parent1)&&(p2!=this->vertices[v3].parent2))
			return (0);

		if ((v2!=this->vertices[v3].parent1)&&(v2!=this->vertices[v3].parent2))
		{
			this->MergeInput->Conquer(f5,v2,v3,f6,v4);

			if (this->MergeInput->IsEdge(p2,v2)>=0)
			{
				if (this->MatchParents(v2,p1,v4,f3,f6,3)!=1)
					return (0);
			}
			if ((this->vertices[p2].old_to_new>=0)&&(this->vertices[v2].old_to_new>=0))
			{
				if (this->MergeOutput->IsEdge(this->vertices[p2].old_to_new,this->vertices[v2].old_to_new)>=0)
				{
					if (this->MatchParents(v2,p1,v4,f3,f6,3)!=1)
						return (0);
				}
			}
		}
	}


	int NO1=0,NO2=0;

	if (this->MergeInput->IsEdge(p2,v2)>=0) NO1=1;
	if ((this->vertices[p2].old_to_new>=0)&&(this->vertices[v2].old_to_new>=0))
	{
		if (this->MergeOutput->IsEdge(this->vertices[p2].old_to_new,this->vertices[v2].old_to_new)>=0)
			NO1=1;
	}

	if (this->MergeInput->IsEdge(p1,v3)>=0)	NO2=1;
	if ((this->vertices[p1].old_to_new>=0)&&(this->vertices[v3].old_to_new>=0))
	{
		if (this->MergeOutput->IsEdge(this->vertices[p1].old_to_new,this->vertices[v3].old_to_new)>=0)
			NO2=1;
	}

	if (NO1+NO2==2)
	{
		return (0);
	}

	this->MergeInput->GetVertexNeighbours(p1,this->vlist);
	int i=this->vlist->GetNumberOfIds()-1;
	for (;i>=0;i--)
	{
		v2=this->vlist->GetId(i);
		if (this->vertices[v2].type==CHILD)
		{
			v3=this->vertices[v2].parent1;
			v4=this->vertices[v2].parent2;

			if ((v3==p2)||(v4==p2))
			{
				if (p2==v4)
				{
					v6=v4;
					v4=v3;
					v3=v6;
				}

				edge1=this->MergeInput->IsEdge(v2,p2);
				this->MergeInput->GetEdgeFaces(edge1,f3,f4);

				if (this->faces[f3].hirestype<0)
				{
					v5=this->MergeInput->GetThirdPoint(f3,p2,v2);

					if ((this->MergeInput->IsFace(v5,p1,v2)>=0)&&
						(this->MergeInput->IsFace(v4,p1,v2)>=0)&&
						(v5!=v1))
					{
						if (this->MergeInput->IsEdge(v4,v5)>=0)
							return (0);
						if (this->vertices[v5].old_to_new>=0)
						{
							if (this->MergeOutput->IsEdge(this->vertices[v4].old_to_new,this->vertices[v5].old_to_new)>=0)
								return (0);
						}
					}
				}
				else
				{
					if (f4>=0)
					{
						f3=f4;
						if (this->faces[f3].hirestype<0)
						{
							v5=this->MergeInput->GetThirdPoint(f3,p2,v2);
							if ((this->MergeInput->IsFace(v5,p1,v2)>=0)&&
								(this->MergeInput->IsFace(v4,p1,v2)>=0)&&
								(v5!=v1))
							{
								if (this->MergeInput->IsEdge(v4,v5)>=0)
									return (0);
								if (this->vertices[v5].old_to_new>=0)
								{
									if (this->MergeOutput->IsEdge(this->vertices[v4].old_to_new,this->vertices[v5].old_to_new)>=0)
										return (0);
								}
							}
						}
					}
				}
			}

			if ((v3=p1)||(v4==p1))
			{
				if (p1==v4)
				{
					v6=v4;
					v4=v3;
					v3=v6;
				}

				edge1=this->MergeInput->IsEdge(v2,p1);
				this->MergeInput->GetEdgeFaces(edge1,f3,f4);

				if (this->faces[f3].hirestype<0)
				{
					v5=this->MergeInput->GetThirdPoint(f3,p1,v2);
					if ((this->MergeInput->IsFace(v5,p2,v2)>=0)&&
						(this->MergeInput->IsFace(v4,p2,v2)>=0)&&
						(v5!=v1))
					{
						if (this->MergeInput->IsEdge(v4,v5)>=0)
							return (0);
						if (this->vertices[v5].old_to_new>=0)
						{
							if (this->MergeOutput->IsEdge(this->vertices[v4].old_to_new,this->vertices[v5].old_to_new)>=0)
								return (0);
						}
					}
				}
				else
				{
					if (f4>=0)
					{
						f3=f4;
						if (this->faces[f3].hirestype<0)
						{
							v5=this->MergeInput->GetThirdPoint(f3,p1,v2);
							if ((this->MergeInput->IsFace(v5,p2,v2)>=0)&&
								(this->MergeInput->IsFace(v4,p2,v2)>=0)&&
								(v5!=v1))
							{
								if (this->MergeInput->IsEdge(v4,v5)>=0)
									return (0);
								if (this->vertices[v5].old_to_new>=0)
								{
									if (this->MergeOutput->IsEdge(this->vertices[v4].old_to_new,this->vertices[v5].old_to_new)>=0)
										return (0);
								}
							}
						}
					}
				}
			}
		}
	}

	return (1);

}

int vtkWaveletSubdivisionFilter::PossibleParent(vtkIdType v1)
{
	vtkIdType v2,v3,v4,f1,f2,e1,f[2];
	int i,j;

	if (vertices[v1].type==CHILD) return (0);
	if (vertices[v1].type==PARENT) return (1);

	return (1); // A.G : i might be wrong but the code should exit here whatever

	this->MergeInput->GetVertexNeighbours(v1,vlist);

	for (i=0;i<this->vlist->GetNumberOfIds();i++)
	{
		v2=this->vlist->GetId(i);
		if (this->vertices[v2].type==CHILD)
		{
			if ((v1!=this->vertices[v2].parent1)&&
				(v1!=this->vertices[v2].parent2))
			{
				e1=this->MergeInput->IsEdge(v1,v2);
				this->MergeInput->GetEdgeFaces(e1,f[0],f[1]);
				for (j=0;j<2;j++)
				{
					f1=f[j];
					v3=this->MergeInput->GetThirdPoint(f1,v1,v2);
					if ((v3!=this->vertices[v2].parent1)&&
						(v3!=this->vertices[v2].parent2))
					{
						this->MergeInput->Conquer(f1,v2,v3,f2,v4);
						if (this->MatchParents(v3,v1,v4,f1,f2,2)!=1) 
							return(0);
					}
				}
			}
		}
	}
	return (1);
}

void vtkWaveletSubdivisionFilter::Merge2To1(vtkIdType v1,vtkIdType v2, vtkIdType v3,
											vtkIdType v4,vtkIdType f1,vtkIdType f2)
{

	vtkIdType f5,edge1,edge2,edge3;
	f5=this->AddFace(v1,v2,v3,edge1,edge2,edge3);

	faces[f5].lowrestype=2;
	faces[f5].vertex=-1;
	faces[f1].hirestype=1;
	faces[f2].hirestype=1;

	edgesvector[edge1].child=v4;
	edgesvector[edge2].child=-1;
	edgesvector[edge3].child=-1;

	vertices[v4].type=CHILD;
	vertices[v1].type=PARENT;
	vertices[v2].type=PARENT;
	vertices[v3].type=PARENT;

	vertices[v4].parent1=v1;
	vertices[v4].parent2=v2;

	if (this->SolveProcess<2)
	{
		switchcontour(edge1);
		switchcontour(edge2);
		switchcontour(edge3);
	}
	else
	{
		switchcontour2(edge1);
		switchcontour2(edge2);
		switchcontour2(edge3);
	}

	this->vertices[v3].valence--;
	this->remaining_faces-=2;
	this->M2to1++;

	if (this->DisplayEfficiency==1)
	{
		this->CellTypes->SetValue(f1,2);
		this->CellTypes->SetValue(f2,2);
	}
}

void vtkWaveletSubdivisionFilter::MergeAndSwap1(vtkIdType v1,vtkIdType v2, vtkIdType v3,
												vtkIdType v4,vtkIdType v5,vtkIdType f1,vtkIdType f2,vtkIdType f3)
{

	vtkIdType f5,edge1,edge2,edge3,edge4,edge5,edge6;

	this->faces[f1].hirestype=1;
	this->faces[f2].hirestype=1;
	this->faces[f3].hirestype=1;

	f5=this->AddFace(v1,v2,v4,edge1,edge2,edge3);

	this->faces[f5].lowrestype=2;
	this->faces[f5].vertex=v1;

	this->edgesvector[edge1].child=v5;
	this->edgesvector[edge2].child=-1;
	this->edgesvector[edge3].child=-1;

	this->vertices[v5].type=CHILD;
	this->vertices[v1].type=PARENT;
	this->vertices[v2].type=PARENT;
	this->vertices[v3].type=PARENT;
	this->vertices[v4].type=PARENT;

	this->vertices[v5].parent1=v1;
	this->vertices[v5].parent2=v2;

	this->M2to1++;
	this->Swap1++;

	f5=this->AddFace(v1,v3,v4,edge4,edge5,edge6);

	this->faces[f5].lowrestype=1;
	this->faces[f5].vertex=-1;

	this->edgesvector[edge4].child=-1;
	this->edgesvector[edge5].child=-1;
	this->edgesvector[edge6].child=-1;

	if (this->SolveProcess<2)
	{
		switchcontour(edge1);
		switchcontour(edge3);
		switchcontour(edge4);
		switchcontour(edge6);
	}
	else
	{
		switchcontour2(edge1);
		switchcontour2(edge3);
		switchcontour2(edge4);
		switchcontour2(edge6);
	}

	this->remaining_faces-=3;
	this->M1to1++;

	this->vertices[v1].valence++;
	this->vertices[v4].valence++;
	this->vertices[v3].valence--;

	if (this->DisplayEfficiency==1)
	{
		this->CellTypes->SetValue(f1,5);
		this->CellTypes->SetValue(f2,5);
		this->CellTypes->SetValue(f3,5);
	}
}

void vtkWaveletSubdivisionFilter::MergeAndSwap2(vtkIdType v1,vtkIdType v2, vtkIdType v3,
												vtkIdType v4,vtkIdType v5,vtkIdType v6,vtkIdType f1,vtkIdType f2,vtkIdType f3, vtkIdType f4)
{
	vtkIdType f5,edge1,edge2,edge3,edge4,edge5,edge6;

	faces[f1].hirestype=1;
	faces[f2].hirestype=1;
	faces[f3].hirestype=1;
	faces[f4].hirestype=1;

	f5=this->AddFace(v1,v3,v6,edge1,edge2,edge3);

	faces[f5].lowrestype=2;
	faces[f5].vertex=v3;

	edgesvector[edge1].child=v4;
	edgesvector[edge2].child=-1;
	edgesvector[edge3].child=-1;

	vertices[v1].type=PARENT;
	vertices[v2].type=PARENT;
	vertices[v3].type=PARENT;
	vertices[v4].type=CHILD;
	vertices[v5].type=CHILD;
	vertices[v6].type=PARENT;

	vertices[v5].parent1=v2;
	vertices[v5].parent2=v3;
	vertices[v4].parent1=v1;
	vertices[v4].parent2=v3;

	this->M2to1++;
	this->Swap2++;

	f5=this->AddFace(v2,v3,v6,edge4,edge5,edge6);

	faces[f5].lowrestype=2;
	faces[f5].vertex=v3;

	edgesvector[edge4].child=v5;
	edgesvector[edge5].child=-1;
	edgesvector[edge6].child=-1;

	if (this->SolveProcess<2)
	{
		switchcontour(edge1);
		switchcontour(edge2);
		switchcontour(edge4);
		switchcontour(edge5);
	}
	else
	{
		switchcontour2(edge1);
		switchcontour2(edge2);
		switchcontour2(edge4);
		switchcontour2(edge5);
	}

	this->remaining_faces-=4;
	this->M2to1++;
	this->vertices[v3].valence++;
	this->vertices[v6].valence++;

	if (this->DisplayEfficiency==1)
	{
		this->CellTypes->SetValue(f1,6);
		this->CellTypes->SetValue(f2,6);
		this->CellTypes->SetValue(f3,6);
		this->CellTypes->SetValue(f4,6);
	}
}

void vtkWaveletSubdivisionFilter::Merge1To1(vtkIdType v1,vtkIdType v2, vtkIdType v3,vtkIdType f1)
{

	vtkIdType f5,edge1,edge2,edge3;
	f5=this->AddFace(v1,v2,v3,edge1,edge2,edge3);

	faces[f5].lowrestype=1;
	faces[f5].vertex=-1;
	faces[f1].hirestype=1;

	edgesvector[edge1].child=-1;
	edgesvector[edge2].child=-1;
	edgesvector[edge3].child=-1;

	vertices[v1].type=PARENT;
	vertices[v2].type=PARENT;
	vertices[v3].type=PARENT;

	if (this->SolveProcess<2)
	{
		switchcontour(edge1);
		switchcontour(edge2);
		switchcontour(edge3);
	}
	else
	{
		switchcontour2(edge1);
		switchcontour2(edge2);
		switchcontour2(edge3);
	}

	this->remaining_faces-=1;

	if (this->DisplayEfficiency==1)
	{
		this->CellTypes->SetValue(f1,1);
	}

	this->M1to1++;
}

void vtkWaveletSubdivisionFilter::Merge3To1(vtkIdType v1,vtkIdType v2, vtkIdType v3,
											vtkIdType v4, vtkIdType v5,vtkIdType f1,vtkIdType f2,vtkIdType f3)
{
	vtkIdType f5,edge1,edge2,edge3;
	f5=this->AddFace(v1,v2,v3,edge1,edge2,edge3);

	faces[f5].lowrestype=3;
	faces[f5].vertex=v5;
	faces[f1].hirestype=1;
	faces[f2].hirestype=1;
	faces[f3].hirestype=1;

	edgesvector[edge1].child=v4;
	edgesvector[edge2].child=v5;
	edgesvector[edge3].child=-1;

	vertices[v4].type=CHILD;
	vertices[v5].type=CHILD;
	vertices[v1].type=PARENT;
	vertices[v2].type=PARENT;
	vertices[v3].type=PARENT;

	vertices[v4].parent1=v1;
	vertices[v4].parent2=v2;
	vertices[v5].parent1=v1;
	vertices[v5].parent2=v3;

	if (this->SolveProcess<2)
	{
		switchcontour(edge1);
		switchcontour(edge2);
		switchcontour(edge3);
	}
	else
	{
		switchcontour2(edge1);
		switchcontour2(edge2);
		switchcontour2(edge3);
	}

	this->vertices[v2].valence--;
	this->remaining_faces-=3;
	this->M3to1++;

	if (this->DisplayEfficiency==1)
	{
		this->CellTypes->SetValue(f1,3);
		this->CellTypes->SetValue(f3,3);
		this->CellTypes->SetValue(f2,3);
	}
}




void vtkWaveletSubdivisionFilter::SolveInverseSubdivision4(vtkIdType initial_face)
{
	int	numPoints=this->MergeInput->GetNumberOfPoints();
	int numCells=this->MergeInput->GetNumberOfCells();
	int numEdges=0;

	vtkIdType i,ptId0,ptId1,ptId2,ptId3,ptId4=-1,ptId5,ptId6=-1,ptId7=-1,ptId8=-1,edge1,edge2,edge3;

	int found,iteration=0;
	int TestLeft,TestRight,MergeType,BestValence=-1,TestValence;

	vtkIdType f1,f2,f3=-1,f4,f5=-1,f6=-1,f7=-1;
	vtkIdType ptId;

	int val1,val2,val3;

	if (this->GeometryCriterion==1)
	{
		this->MergeInput->ComputeSharpVertices(this->CurvatureTreshold);
	}

	this->SolveProcess=3;

	this->SubdivisionType=1;

	numEdges=this->MergeInput->GetNumberOfEdges();

	this->edgesvector.resize(numEdges);
	this->faces.resize(numCells);
	this->vertices.resize(numPoints);

	for (i=0;i<numEdges;i++)
	{
		this->edgesvector[i].contour=0;
		this->edgesvector[i].child=-1;
	}

	for (i=0;i<numCells;i++)
	{
		this->faces[i].hirestype=-1;
		this->faces[i].lowrestype=-1;
	}

	for (i=0;i<numPoints;i++)
	{
		this->vertices[i].type=-1;
		this->vertices[i].parent1=-1;
		this->vertices[i].parent2=-1;
		this->vertices[i].new_to_old=-1;
		this->vertices[i].old_to_new=-1;
	}

	//////////////////////////////////////////////////////////////////////

	this->MergeOutput->Init(this->MergeInput);

	for (ptId=0;ptId<numPoints;ptId++)
	{
		this->MergeInput->GetVertexNeighbourFaces(ptId,this->vlist);

		int essai=this->vertices[ptId].valence=this->vlist->GetNumberOfIds();

		if (essai>6) this->vertices[ptId].type=PARENT;

		ptId1=this->MergeInput->GetNumberOfBoundaries(ptId);

		if (ptId1>1) this->vertices[ptId].type=PARENT;
		if ((ptId1==1)&&(essai<2)) this->vertices[ptId].type=PARENT;
		if ((ptId1==1)&&(essai>4)) this->vertices[ptId].type=PARENT;
		if ((ptId1==0)&&(essai<4)) this->vertices[ptId].type=PARENT;
	}

	this->remaining_faces=this->MergeInput->GetNumberOfCells();

	this->M1to1=0;
	this->M2to1=0;
	this->M3to1=0;
	this->M4to1=0;
	this->Swap1=0;
	this->Swap2=0;

	found=0;

	for (f1=0;f1<this->MergeInput->GetNumberOfCells();f1++)
	{
		if (this->faces[f1].hirestype<0)
		{
			this->MergeInput->GetFaceVertices(f1,ptId0,ptId1,ptId2);
			this->MergeInput->Conquer(f1,ptId0,ptId1,f2,ptId3);
			this->MergeInput->Conquer(f1,ptId0,ptId2,f3,ptId4);
			this->MergeInput->Conquer(f1,ptId1,ptId2,f4,ptId5);

			if ((f2>=0)&&(f3>=0)&&(f4>=0))
			{
				val1=this->vertices[ptId3].valence;
				val2=this->vertices[ptId4].valence;
				val3=this->vertices[ptId5].valence;

				if ((this->faces[f2].hirestype<0)&&
					(this->faces[f3].hirestype<0)&&
					(this->faces[f4].hirestype<0)
					&&(((val1!=6)&&(this->MergeInput->GetNumberOfBoundaries(ptId3)==0))
					||((val2!=6)&&(this->MergeInput->GetNumberOfBoundaries(ptId4)==0))
					||((val3!=6)&&(this->MergeInput->GetNumberOfBoundaries(ptId5)==0)))
					&&(this->vertices[ptId3].type!=CHILD)
					&&(this->vertices[ptId4].type!=CHILD)
					&&(this->vertices[ptId5].type!=CHILD))
				{
					if (this->MatchParents(ptId0,ptId3,ptId4,f2,f3,3)==1)
					{
						if (this->MatchParents(ptId1,ptId3,ptId5,f2,f4,3)==1)
						{
							if (this->MatchParents(ptId2,ptId4,ptId5,f3,f4,3)==1)
							{
								f5=this->AddFace(ptId3,ptId5,ptId4,edge2,edge1,edge3);	

								faces[f5].lowrestype=4;

								faces[f1].hirestype=1;
								faces[f2].hirestype=1;
								faces[f3].hirestype=1;
								faces[f4].hirestype=1;

								edgesvector[edge1].child=ptId0;
								edgesvector[edge2].child=ptId1;
								edgesvector[edge3].child=ptId2;

								vertices[ptId0].type=CHILD;
								vertices[ptId1].type=CHILD;
								vertices[ptId2].type=CHILD;

								vertices[ptId3].type=PARENT;
								vertices[ptId4].type=PARENT;
								vertices[ptId5].type=PARENT;

								vertices[ptId0].parent1=ptId3;
								vertices[ptId0].parent2=ptId4;

								vertices[ptId1].parent1=ptId3;
								vertices[ptId1].parent2=ptId5;

								vertices[ptId2].parent1=ptId4;
								vertices[ptId2].parent2=ptId5;

								switchcontour2(edge1);
								switchcontour2(edge2);
								switchcontour2(edge3);

								found=1;

								this->M4to1++;

								this->remaining_faces-=4;

								if (this->DisplayEfficiency==1)
								{
									this->CellTypes->SetValue(f1,4);
									this->CellTypes->SetValue(f2,4);
									this->CellTypes->SetValue(f3,4);
									this->CellTypes->SetValue(f4,4);
								}
								this->ConquerRegularSubdivisions();
							}
						}	
					}
				}
			}
		}
	}

	for (f1=0;f1<this->MergeInput->GetNumberOfCells();f1++)
	{
		if (this->faces[f1].hirestype<0)
		{
			this->MergeInput->GetFaceVertices(f1,ptId0,ptId1,ptId2);
			this->MergeInput->Conquer(f1,ptId0,ptId1,f2,ptId3);
			this->MergeInput->Conquer(f1,ptId0,ptId2,f3,ptId4);
			this->MergeInput->Conquer(f1,ptId1,ptId2,f4,ptId5);

			val1=this->vertices[ptId0].valence;
			val2=this->vertices[ptId1].valence;
			val3=this->vertices[ptId2].valence;

			if ((f2>=0)&&(f3>=0)&&(f4>=0))
			{
				val1=this->vertices[ptId3].valence;
				val2=this->vertices[ptId4].valence;
				val3=this->vertices[ptId5].valence;		
				if ((this->faces[f2].hirestype<0)&&
					(this->faces[f3].hirestype<0)&&
					(this->faces[f4].hirestype<0)
					&&(this->vertices[ptId3].type!=CHILD)
					&&(this->vertices[ptId4].type!=CHILD)
					&&(this->vertices[ptId5].type!=CHILD))
				{
					if (this->MatchParents(ptId0,ptId3,ptId4,f2,f3,3)==1)
					{
						if (this->MatchParents(ptId1,ptId3,ptId5,f2,f4,3)==1)
						{
							if (this->MatchParents(ptId2,ptId4,ptId5,f3,f4,3)==1)
							{
								f5=this->AddFace(ptId3,ptId5,ptId4,edge2,edge1,edge3);	

								faces[f5].lowrestype=4;

								faces[f1].hirestype=1;
								faces[f2].hirestype=1;
								faces[f3].hirestype=1;
								faces[f4].hirestype=1;

								edgesvector[edge1].child=ptId0;
								edgesvector[edge2].child=ptId1;
								edgesvector[edge3].child=ptId2;

								vertices[ptId0].type=CHILD;
								vertices[ptId1].type=CHILD;
								vertices[ptId2].type=CHILD;

								vertices[ptId3].type=PARENT;
								vertices[ptId4].type=PARENT;
								vertices[ptId5].type=PARENT;

								vertices[ptId0].parent1=ptId3;
								vertices[ptId0].parent2=ptId4;

								vertices[ptId1].parent1=ptId3;
								vertices[ptId1].parent2=ptId5;

								vertices[ptId2].parent1=ptId4;
								vertices[ptId2].parent2=ptId5;

								switchcontour2(edge1);
								switchcontour2(edge2);
								switchcontour2(edge3);

								found=1;

								this->M4to1++;

								this->remaining_faces-=4;	

								if (this->DisplayEfficiency==1)
								{
									this->CellTypes->SetValue(f1,4);
									this->CellTypes->SetValue(f2,4);
									this->CellTypes->SetValue(f3,4);
									this->CellTypes->SetValue(f4,4);
								}
								this->ConquerRegularSubdivisions();
							}
						}	
					}
				}
			}
		}
	}
	if (found==0)
	{
		f1=initial_face;
		this->MergeInput->GetFaceVertices(f1,ptId0,ptId1,ptId2);
		this->Merge1To1(ptId0,ptId1,ptId2,f1);
	}


	while ((this->FirstElement!=0)||(this->FirstElementRegular!=0))
	{
		this->ConquerRegularSubdivisions();

		if (this->FirstElement==0)
		{
			if (this->GeometryCriterion==1)
			{
				this->MergeInput->DeleteSharpVertices();
			}

			this->MergeOutput->Squeeze();
			return;
		}

		found=0;

		edge1=FirstElement->edge;

		ptId1=edgesvector[edge1].child;

		if (ptId1<0) 
			goto failchild;

		this->MergeOutput->GetEdgeVertices(edge1,ptId0,ptId2);

		ptId0=this->vertices[ptId0].new_to_old;
		ptId2=this->vertices[ptId2].new_to_old;

		this->ConquerEdge(ptId0,ptId1,f1,ptId3);

		this->MergeInput->Conquer(f1,ptId1,ptId3,f2,ptId4);

		if (ptId4==ptId2)
			goto fail1;


		this->MergeInput->Conquer(f2,ptId4,ptId1,f3,ptId2);

		this->MergeInput->Conquer(f2,ptId3,ptId4,f4,ptId5);

		goto fail4;
		//  fin fusion 4->1

fail1:
		TestLeft=0;
		TestRight=0;

		this->MergeInput->Conquer(f1,ptId0,ptId3,f3,ptId4);
		if (f3==-1)
			goto FailLeft1;

		if (faces[f3].hirestype>0)
			goto FailLeft1;

		if ((PossibleParent(ptId4)!=1)||(MatchParents(ptId3,ptId2,ptId4,f2,f3,3)!=1))
			goto FailLeft1;

		TestLeft=1;

FailLeft1:
		this->MergeInput->Conquer(f2,ptId2,ptId3,f4,ptId5);

		if (f4==-1)
			goto FailRight1;

		if (this->faces[f4].hirestype>0)
			goto FailRight1;

		if ((this->PossibleParent(ptId5)!=1)||
			(this->MatchParents(ptId3,ptId0,ptId5,f1,f4,3)!=1))
			goto FailRight1;

		TestRight=1;

FailRight1:
		if (TestLeft==1)
		{
			if (TestRight==1)
			{
				if (this->vertices[ptId0].valence>this->vertices[ptId2].valence)
				{
					Merge3To1(ptId2,ptId0,ptId4,ptId1,ptId3,f1,f2,f3);
					goto Success;
				}
				else
				{
					Merge3To1(ptId0,ptId2,ptId5,ptId1,ptId3,f1,f2,f4);
					goto Success;
				}
			}
			else
			{
				Merge3To1(ptId2,ptId0,ptId4,ptId1,ptId3,f1,f2,f3);
				goto Success;
			}
		}
		else
		{
			if (TestRight==1)
			{
				Merge3To1(ptId0,ptId2,ptId5,ptId1,ptId3,f1,f2,f4);
				goto Success;
			}
		}


		if (this->PossibleParent(ptId3)!=1)
			goto fail;

		Merge2To1(ptId0,ptId2,ptId3,ptId1,f1,f2);

		goto Success;

fail4:
		TestRight=0;
		TestLeft=0;

		if ((PossibleParent(ptId4)==1)&&
			(MatchParents(ptId3,ptId0,ptId4,f1,f2,2)==1))
			TestLeft=1;

		if ((PossibleParent(ptId3)==1)&&
			(MatchParents(ptId4,ptId2,ptId3,f3,f2,2)==1))
			TestRight=1;

		if (TestLeft==1)
		{
			if (TestRight==1)
			{
				if (this->vertices[ptId4].valence>this->vertices[ptId3].valence)
				{
					Merge3To1(ptId0,ptId4,ptId2,ptId3,ptId1,f1,f2,f3);
					goto Success;
				}
				else
				{
					Merge3To1(ptId2,ptId3,ptId0,ptId4,ptId1,f1,f2,f3);
					goto Success;
				}
			}
			else
			{
				Merge3To1(ptId0,ptId4,ptId2,ptId3,ptId1,f1,f2,f3);
				goto Success;
			}
		}
		else
		{
			if (TestRight==1)
			{
				Merge3To1(ptId2,ptId3,ptId0,ptId4,ptId1,f1,f2,f3);
				goto Success;
			}
		}

		///ici interviennent les permutations d'arête....

		TestLeft=0;
		TestRight=0;

		this->MergeInput->Conquer(f2,ptId3,ptId4,f4,ptId5);

		if (f4<0)
			goto failperm2;

		if (this->faces[f4].hirestype>=0)
			goto failperm2;

		if (this->PossibleParent(ptId5)!=1)
			goto failperm2;

		if (this->PossibleParent(ptId4)!=1)
			goto FailLeft4;

		if (this->MatchParents(ptId3,ptId0,ptId5,f1,f4,3)!=1) 
			goto FailLeft4;

		if (this->MergeInput->IsEdge(ptId0,ptId4)>=0) 
			goto FailLeft4;

		if ((this->vertices[ptId0].old_to_new>=0)&&
			(this->vertices[ptId4].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId0].old_to_new,this->vertices[ptId4].old_to_new)>=0)
				goto FailLeft4;
		}

		TestLeft=1;

FailLeft4:

		if (this->PossibleParent(ptId3)!=1) 
			goto FailRight4;

		if (this->MatchParents(ptId4,ptId2,ptId5,f3,f4,3)!=1) 
			goto FailRight4;

		if (this->MergeInput->IsEdge(ptId2,ptId3)>=0) 
			goto FailRight4;

		if ((this->vertices[ptId2].old_to_new>=0)&&
			(this->vertices[ptId3].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId2].old_to_new,this->vertices[ptId3].old_to_new)>=0) 
				goto FailRight4;
		}

		TestRight=1;

FailRight4:
		if (TestLeft==1)
		{
			if (TestRight==1)
			{
				if (this->vertices[ptId4].valence+this->vertices[ptId0].valence
					<this->vertices[ptId3].valence+this->vertices[ptId2].valence)
				{
					this->MergeAndSwap2(ptId5,ptId2,ptId0,ptId3,ptId1,ptId4,f1,f2,f3,f4);
					goto Success;
				}
				else
				{
					this->MergeAndSwap2(ptId0,ptId5,ptId2,ptId1,ptId4,ptId3,f1,f2,f3,f4);
					goto Success;
				}
			}
			else
			{
				this->MergeAndSwap2(ptId5,ptId2,ptId0,ptId3,ptId1,ptId4,f1,f2,f3,f4);
				goto Success;
			}
		}
		else
		{
			if (TestRight==1)
			{
				this->MergeAndSwap2(ptId0,ptId5,ptId2,ptId1,ptId4,ptId3,f1,f2,f3,f4);
				goto Success;
			}
		}

failperm2:
		if (this->PossibleParent(ptId3)!=1) 
			goto fail;

		if (this->PossibleParent(ptId4)!=1) 
			goto fail;

		TestLeft=0;
		TestRight=0;

		if (this->MergeInput->IsEdge(ptId2,ptId3)>=0) 
			goto FailLeft5;

		if ((this->vertices[ptId2].old_to_new>=0)&&
			(this->vertices[ptId3].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId2].old_to_new,this->vertices[ptId3].old_to_new)>=0) 
				goto FailLeft5;
		}

		TestLeft=1;

FailLeft5:
		if (this->MergeInput->IsEdge(ptId0,ptId4)>=0) 
			goto FailRight5;

		if ((this->vertices[ptId0].old_to_new>=0)&&
			(this->vertices[ptId4].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId0].old_to_new,this->vertices[ptId4].old_to_new)>=0) 
				goto FailRight5;
		}

		TestRight=1;

FailRight5:
		if (TestLeft==1)
		{
			if (TestRight==1)
			{
				if (this->vertices[ptId4].valence+this->vertices[ptId0].valence
					-this->vertices[ptId3].valence>this->vertices[ptId3].valence
					+this->vertices[ptId2].valence-this->vertices[ptId4].valence)
				{
					this->MergeAndSwap1(ptId2,ptId0,ptId4,ptId3,ptId1,f1,f2,f3);
					goto Success;
				}
				else
				{
					this->MergeAndSwap1(ptId0,ptId2,ptId3,ptId4,ptId1,f1,f2,f3);
					goto Success;
				}
			}
			else
			{
				this->MergeAndSwap1(ptId2,ptId0,ptId4,ptId3,ptId1,f1,f2,f3);
				goto Success;
			}
		}
		else
		{
			if (TestRight==1)
			{
				this->MergeAndSwap1(ptId0,ptId2,ptId3,ptId4,ptId1,f1,f2,f3);
				goto Success;
			}
		}

		goto fail;
		/////////////////////////////////////////////

failchild:
		this->MergeOutput->GetEdgeVertices(edge1,ptId0,ptId1);

		ptId0=this->vertices[ptId0].new_to_old;
		ptId1=this->vertices[ptId1].new_to_old;

		this->ConquerEdge(ptId0,ptId1,f1,ptId2);

		TestLeft=0;
		TestRight=0;

		this->MergeInput->Conquer(f1,ptId0,ptId2,f2,ptId3);

		if (f2<0)
			goto FailLeft2;

		if (this->faces[f2].hirestype>0) 
			goto FailLeft2;

		this->MergeInput->Conquer(f2,ptId2,ptId3,f3,ptId4);

		if (f3<0) 
			goto FailLeft2;

		if (this->faces[f3].hirestype>=0) 
			goto FailLeft2;

		if (this->MatchParents(ptId2,ptId1,ptId4,f1,f3,3)!=1) 
			goto FailLeft2;

		if (this->MatchParents(ptId3,ptId0,ptId4,f2,f3,2)!=1) 
			goto FailLeft2;

		if (this->PossibleParent(ptId4)!=1) 
			goto FailLeft2;

		TestLeft=1;

FailLeft2:
		this->MergeInput->Conquer(f1,ptId1,ptId2,f4,ptId5);

		if (f4<0) 
			goto FailRight2;

		if (this->faces[f4].hirestype>0) 
			goto FailRight2;

		this->MergeInput->Conquer(f4,ptId2,ptId5,f5,ptId6);

		if (f5<0) 
			goto FailRight2;

		if (this->faces[f5].hirestype>0) 
			goto FailRight2;

		if (this->MatchParents(ptId5,ptId6,ptId1,f5,f4,2)!=1) 
			goto FailRight2;

		if (this->MatchParents(ptId2,ptId6,ptId0,f5,f1,3)!=1) 
			goto FailRight2;

		if (this->PossibleParent(ptId6)!=1) 
			goto FailRight2;

		TestRight=1;

FailRight2:
		if (TestLeft==1)
		{
			if (TestRight==1)
			{
				if (this->vertices[ptId0].valence>this->vertices[ptId1].valence)
				{
					this->Merge3To1(ptId4,ptId0,ptId1,ptId3,ptId2,f1,f2,f3);
					goto Success;
				}
				else
				{
					this->Merge3To1(ptId6,ptId1,ptId0,ptId5,ptId2,f1,f4,f5);
					goto Success;
				}
			}
			else
			{
				this->Merge3To1(ptId4,ptId0,ptId1,ptId3,ptId2,f1,f2,f3);
				goto Success;
			}
		}
		else
		{
			if (TestRight==1)
			{
				this->Merge3To1(ptId6,ptId1,ptId0,ptId5,ptId2,f1,f4,f5);
				goto Success;
			}
		}

		TestLeft=0;
		TestRight=0;

		this->MergeInput->Conquer(f1,ptId0,ptId2,f2,ptId3);

		if (f2<0) 
			goto FailLeft3;

		if (this->faces[f2].hirestype>0) 
			goto FailLeft3;

		if (this->PossibleParent(ptId3)!=1) 
			goto FailLeft3;

		if (this->MatchParents(ptId2,ptId1,ptId3,f1,f2,2)!=1) 
			goto FailLeft3;

		TestLeft=1;

FailLeft3:
		this->MergeInput->Conquer(f1,ptId1,ptId2,f3,ptId4);

		if (f3<0) 
			goto FailRight3;

		if (this->faces[f3].hirestype>0) 
			goto FailRight3;

		if (this->PossibleParent(ptId4)!=1) 
			goto FailRight3;

		if (this->MatchParents(ptId2,ptId0,ptId4,f1,f3,2)!=1) 
			goto FailRight3;

		TestRight=1;

FailRight3:
		if (TestLeft==1)
		{
			if (TestRight==1)
			{
				if (this->vertices[ptId0].valence>this->vertices[ptId1].valence)
				{
					this->Merge2To1(ptId1,ptId3,ptId0,ptId2,f1,f2);
					goto Success;
				}
				else
				{
					this->Merge2To1(ptId0,ptId4,ptId1,ptId2,f1,f3);
					goto Success;
				}
			}
			else
			{
				this->Merge2To1(ptId1,ptId3,ptId0,ptId2,f1,f2);
				goto Success;
			}
		}
		else
		{
			if (TestRight==1)
			{
				this->Merge2To1(ptId0,ptId4,ptId1,ptId2,f1,f3);
				goto Success;
			}
		}

		//	ici interviennent les permutations d'arête

		TestLeft=0;
		TestRight=0;

		this->MergeInput->Conquer(f1,ptId0,ptId2,f2,ptId3);

		if (f2<0) 
			goto FailLeft6;

		if (this->faces[f2].hirestype>=0) 
			goto FailLeft6;

		this->MergeInput->Conquer(f2,ptId0,ptId3,f3,ptId4);

		if (f3<0) 
			goto FailLeft6;

		if (this->faces[f3].hirestype>=0) 
			goto FailLeft6;

		//
		this->MergeInput->Conquer(f2,ptId3,ptId2,f4,ptId5);
		if (f4<0) 
			goto FailLeft6;

		if (this->faces[f4].hirestype>=0) 
			goto FailLeft6;

		if (this->PossibleParent(ptId4)!=1) 
			goto FailLeft6;

		if (this->PossibleParent(ptId5)!=1) 
			goto FailLeft6;

		if (this->MatchParents(ptId3,ptId4,ptId5,f3,f4,3)!=1)
			goto FailLeft6;

		if (this->MatchParents(ptId2,ptId1,ptId5,f1,f4,3)!=1)
			goto FailLeft6;

		if (this->MergeInput->IsEdge(ptId0,ptId5)>=0) 
			goto FailLeft6;

		if ((this->vertices[ptId0].old_to_new>=0)&&
			(this->vertices[ptId5].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId0].old_to_new,this->vertices[ptId5].old_to_new)>=0) 
				goto FailLeft6;
		}

		TestLeft=1;

FailLeft6:
		this->MergeInput->Conquer(f1,ptId1,ptId2,f5,ptId6);

		if (f5<0) 
			goto FailRight6;

		if (this->faces[f5].hirestype>=0) 
			goto FailRight6;

		this->MergeInput->Conquer(f5,ptId1,ptId6,f6,ptId7);
		if (f6<0) 
			goto FailRight6;

		if (this->faces[f6].hirestype>=0) 
			goto FailRight6;
		//
		this->MergeInput->Conquer(f5,ptId6,ptId2,f7,ptId8);

		if (f7<0) 
			goto FailRight6;

		if (this->faces[f7].hirestype>=0) 
			goto FailRight6;


		if (this->PossibleParent(ptId7)!=1) 
			goto FailRight6;

		if (this->PossibleParent(ptId8)!=1) 
			goto FailRight6;

		if (this->MatchParents(ptId6,ptId7,ptId8,f6,f7,3)!=1) 
			goto FailRight6;

		if (this->MatchParents(ptId2,ptId0,ptId8,f1,f7,3)!=1)
			goto FailRight6;

		if (this->MergeInput->IsEdge(ptId1,ptId8)>=0)
			goto FailRight6;

		if ((this->vertices[ptId1].old_to_new>=0)&&(this->vertices[ptId8].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId1].old_to_new,this->vertices[ptId8].old_to_new)>=0) 
				goto FailRight6;
		}

		TestRight=1;

FailRight6:
		if (TestLeft==1)
		{
			if (TestRight==1)
			{
				if (this->vertices[ptId0].valence+this->vertices[ptId5].valence<
					this->vertices[ptId1].valence+this->vertices[ptId8].valence)
				{
					this->MergeAndSwap2(ptId1,ptId4,ptId5,ptId2,ptId3,ptId0,f1,f2,f3,f4);
					goto Success;
				}
				else
				{
					this->MergeAndSwap2(ptId0,ptId7,ptId8,ptId2,ptId6,ptId1,f1,f5,f6,f7);
					goto Success;
				}
			}
			else
			{
				this->MergeAndSwap2(ptId1,ptId4,ptId5,ptId2,ptId3,ptId0,f1,f2,f3,f4);
				goto Success;
			}
		}
		else
		{
			if (TestRight==1)
			{
				this->MergeAndSwap2(ptId0,ptId7,ptId8,ptId2,ptId6,ptId1,f1,f5,f6,f7);
				goto Success;
			}
		}


		this->MergeInput->Conquer(f1,ptId1,ptId2,f2,ptId3);

		MergeType=0;

		if (f2<0) 
			goto FailRight7;

		if (this->faces[f2].hirestype>=0) 
			goto FailRight7;

		this->MergeInput->Conquer(f1,ptId0,ptId2,f3,ptId4);

		if (f3<0) 
			goto FailRight7;

		if (this->faces[f3].hirestype>=0) 
			goto FailRight7;

		if (this->PossibleParent(ptId3)!=1) 
			goto FailRight7;

		if (this->PossibleParent(ptId4)!=1) 
			goto FailRight7;

		if (this->MatchParents(ptId2,ptId3,ptId4,f2,f3,3)!=1) 
			goto FailRight7;

		TestLeft=0;
		TestRight=0;

		if (this->MergeInput->IsEdge(ptId1,ptId4)>=0) 
			goto FailLeft7;

		if ((this->vertices[ptId1].old_to_new>=0)&&
			(this->vertices[ptId4].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId1].old_to_new,this->vertices[ptId4].old_to_new)>=0) 
				goto FailLeft7;
		}

		MergeType=1;

		BestValence=this->vertices[ptId1].valence+this->vertices[ptId4].valence-this->vertices[ptId0].valence;


FailLeft7:
		if (this->MergeInput->IsEdge(ptId0,ptId3)>=0) 
			goto FailRight7;

		if ((this->vertices[ptId0].old_to_new>=0)&&
			(this->vertices[ptId3].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId0].old_to_new,this->vertices[ptId3].old_to_new)>=0) 
				goto FailRight7;
		}

		TestRight=1;

		TestValence=this->vertices[ptId0].valence+this->vertices[ptId3].valence-this->vertices[ptId1].valence;

		if (MergeType==0)
		{
			MergeType=2;
			BestValence=TestValence;
		}
		else
		{
			if (TestValence<BestValence)
			{
				MergeType=2;
				BestValence=TestValence;
			}
		}

FailRight7:
		this->MergeInput->Conquer(f1,ptId0,ptId2,f3,ptId4);

		if (f3<0) 
			goto FailRight8;

		if (this->faces[f3].hirestype>=0) 
			goto FailRight8;

		this->MergeInput->Conquer(f3,ptId2,ptId4,f5,ptId6);

		if (f5<0) 
			goto FailRight8;

		if (this->faces[f5].hirestype>=0) 
			goto FailRight8;

		if (this->PossibleParent(ptId4)!=1) 
			goto FailRight8;

		if (this->PossibleParent(ptId6)!=1) 
			goto FailRight8;

		if (this->MatchParents(ptId2,ptId1,ptId6,f1,f5,3)!=1) 
			goto FailRight8;

		if (this->MergeInput->IsEdge(ptId0,ptId6)>=0) 
			goto FailLeft8;

		if ((this->vertices[ptId0].old_to_new>=0)&&
			(this->vertices[ptId6].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId0].old_to_new,this->vertices[ptId6].old_to_new)>=0) 
				goto FailLeft8;
		}

		TestValence=this->vertices[ptId0].valence+this->vertices[ptId6].valence-this->vertices[ptId4].valence;

		if (MergeType==0)
		{
			MergeType=3;
			BestValence=TestValence;
		}
		else
		{
			if (TestValence<BestValence)
			{
				MergeType=3;
				BestValence=TestValence;
			}
		}

FailLeft8:
		if (this->MergeInput->IsEdge(ptId1,ptId4)>=0) 
			goto FailRight8;

		if ((this->vertices[ptId1].old_to_new>=0)&&
			(this->vertices[ptId4].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId1].old_to_new,this->vertices[ptId4].old_to_new)>=0) 
				goto FailRight8;
		}

		TestValence=this->vertices[ptId1].valence+this->vertices[ptId4].valence-this->vertices[ptId0].valence;

		if (MergeType==0)
		{
			MergeType=4;
			BestValence=TestValence;
		}
		else
		{
			if (TestValence<BestValence)
			{
				MergeType=4;
				BestValence=TestValence;
			}
		}

FailRight8:
		this->MergeInput->Conquer(f1,ptId1,ptId2,f2,ptId3);

		if (f2<0) 
			goto FailRight9;

		if (this->faces[f2].hirestype>=0) 
			goto FailRight9;

		this->MergeInput->Conquer(f2,ptId2,ptId3,f4,ptId5);

		if (f4<0) 
			goto FailRight9;

		if (this->faces[f4].hirestype>=0)
			goto FailRight9;

		if (this->PossibleParent(ptId3)!=1) 
			goto FailRight9;

		if (this->PossibleParent(ptId5)!=1) 
			goto FailRight9;

		if (this->MatchParents(ptId2,ptId0,ptId5,f1,f4,3)!=1) 
			goto FailRight9;

		if (this->MergeInput->IsEdge(ptId0,ptId3)>=0) 
			goto FailLeft9;

		if ((this->vertices[ptId0].old_to_new>=0)&&
			(this->vertices[ptId3].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId0].old_to_new,this->vertices[ptId3].old_to_new)>=0) 
				goto FailLeft9;
		}

		TestValence=this->vertices[ptId0].valence+this->vertices[ptId3].valence-this->vertices[ptId1].valence;

		if (MergeType==0)
		{
			MergeType=5;
			BestValence=TestValence;
		}
		else
		{
			if (TestValence<BestValence)
			{
				MergeType=5;
				BestValence=TestValence;
			}
		}


FailLeft9:
		if (this->MergeInput->IsEdge(ptId1,ptId5)>=0) 
			goto FailRight9;

		if ((this->vertices[ptId1].old_to_new>=0)&&
			(this->vertices[ptId5].old_to_new>=0))
		{
			if (this->MergeOutput->IsEdge(this->vertices[ptId1].old_to_new,this->vertices[ptId5].old_to_new)>=0)
				goto FailRight9;
		}

		TestValence=this->vertices[ptId1].valence+this->vertices[ptId5].valence-this->vertices[ptId3].valence;

		if (MergeType==0)
		{
			MergeType=6;
			BestValence=TestValence;
		}
		else
		{
			if (TestValence<BestValence)
			{
				MergeType=6;
				BestValence=TestValence;
			}
		}

FailRight9:
		switch (MergeType)
		{
		case 1:
			{
				this->MergeAndSwap1(ptId4,ptId3,ptId0,ptId1,ptId2,f1,f2,f3);
				goto Success;
			}
		case 2:
			{
				this->MergeAndSwap1(ptId3,ptId4,ptId1,ptId0,ptId2,f1,f2,f3);
				goto Success;
			}
		case 3:
			{
				this->MergeAndSwap1(ptId6,ptId1,ptId4,ptId0,ptId2,f1,f3,f5);
				goto Success;
			}
		case 4:
			{
				this->MergeAndSwap1(ptId1,ptId6,ptId0,ptId4,ptId2,f1,f3,f5);
				goto Success;
			}
		case 5:
			{
				this->MergeAndSwap1(ptId0,ptId5,ptId1,ptId3,ptId2,f1,f2,f4);
				goto Success;
			}
		case 6:
			{
				this->MergeAndSwap1(ptId5,ptId0,ptId3,ptId1,ptId2,f1,f2,f4);
				goto Success;
			}

		case 0:
			{
				goto faillast;
			}
		}


faillast:	
		if (this->PossibleParent(ptId2)!=1)
			goto fail;

		this->Merge1To1(ptId0,ptId1,ptId2,f1);

Success:
		found=1;	

fail:
		if (found==0)
		{
			edge1=FirstElement->edge;
			switchcontour(edge1);
		}
		iteration++;
	}

	if (this->GeometryCriterion==1)
		this->MergeInput->DeleteSharpVertices();

	this->MergeOutput->Squeeze();

}

void vtkWaveletSubdivisionFilter::switchcontour2(vtkIdType edge)
{
	vtkIdType v1,v2,v3,f1,f2,edge2;
	contour_element *next,*previous;

	this->MergeOutput->GetEdgeVertices(edge,v1,v2);

	v1=this->vertices[v1].new_to_old;
	v2=this->vertices[v2].new_to_old;

	if ((v3=edgesvector[edge].child)<0)
	{
		edge2=this->MergeInput->IsEdge(v1,v2);
		this->MergeInput->GetEdgeFaces(edge2,f1,f2);

		if (f2<0) 
			return;
	}
	else
	{
		edge2=this->MergeInput->IsEdge(v1,v3);
		this->MergeInput->GetEdgeFaces(edge2,f1,f2);

		if (f2<0) 
			return;
	}

	if (edgesvector[edge].contour==0)
	{
		if (v3<0)
		{
			edge2=this->MergeInput->IsEdge(v1,v2);
			this->MergeInput->GetEdgeFaces(edge2,f1,f2);

			if ((this->faces[f1].hirestype>=0)&&
				(this->faces[f2].hirestype>=0))
				return;

			edgesvector[edge].contour=new contour_element;
			edgesvector[edge].contour->edge=edge;

			if (LastElement==0)
			{
				FirstElement=edgesvector[edge].contour;
				LastElement=edgesvector[edge].contour;

				edgesvector[edge].contour->NextElement=0;
				edgesvector[edge].contour->PreviousElement=0;
			}
			else
			{
				edgesvector[edge].contour->NextElement=0;
				edgesvector[edge].contour->PreviousElement=LastElement;

				LastElement->NextElement=edgesvector[edge].contour;
				LastElement=edgesvector[edge].contour;
			}
		}
		else
		{
			edge2=this->MergeInput->IsEdge(v1,v3);
			this->MergeInput->GetEdgeFaces(edge2,f1,f2);

			if ((this->faces[f1].hirestype>=0)&&
				(this->faces[f2].hirestype>=0)) 
				return;

			edgesvector[edge].contour=new contour_element;
			edgesvector[edge].contour->edge=edge;

			if (LastElementRegular==0)
			{
				FirstElementRegular=edgesvector[edge].contour;
				LastElementRegular=edgesvector[edge].contour;

				edgesvector[edge].contour->NextElement=0;
				edgesvector[edge].contour->PreviousElement=0;
			}
			else
			{
				edgesvector[edge].contour->NextElement=0;
				edgesvector[edge].contour->PreviousElement=LastElementRegular;

				LastElementRegular->NextElement=edgesvector[edge].contour;
				LastElementRegular=edgesvector[edge].contour;
			}
		}
	}
	else	// if (edgesvector[edge].contour!=0)
	{
		next=this->edgesvector[edge].contour->NextElement;
		previous=this->edgesvector[edge].contour->PreviousElement;

		if ((next==0)&&(previous==0))
		{
			if (this->FirstElement==this->edgesvector[edge].contour)
			{
				this->LastElement=0;
				this->FirstElement=0;

				delete this->edgesvector[edge].contour;

				this->edgesvector[edge].contour=0;

				return;
			}
			else
			{
				this->LastElementRegular=0;
				this->FirstElementRegular=0;

				delete this->edgesvector[edge].contour;

				this->edgesvector[edge].contour=0;

				return;
			}
		}
		else
		{
			if (next==0)
			{
				previous->NextElement=0;

				if (this->LastElement==this->edgesvector[edge].contour)
				{
					this->LastElement=previous;
				}
				else
				{
					this->LastElementRegular=previous;
				}

				delete this->edgesvector[edge].contour;

				this->edgesvector[edge].contour=0;

				return;
			}
			else
			{
				if (previous==0)
				{
					if (this->FirstElement==this->edgesvector[edge].contour)
					{
						this->FirstElement=next;
					}
					else
					{
						this->FirstElementRegular=next;
					}

					next->PreviousElement=0;

					delete this->edgesvector[edge].contour;

					this->edgesvector[edge].contour=0;

					return;
				}
				else
				{
					next->PreviousElement=previous;
					previous->NextElement=next;

					delete this->edgesvector[edge].contour;

					this->edgesvector[edge].contour=0;

					return;		
				}
			}
		}
	}
}

void vtkWaveletSubdivisionFilter::ConquerRegularSubdivisions()
{
	vtkIdType edge1,edge2,edge3,ptId0,ptId1,ptId2,ptId3,ptId4,ptId5,f1,f2,f3,f4,f5;
	int found;
	contour_element *next,*previous;

	while ((this->FirstElementRegular!=0))
	{
		while (this->FirstElementRegular!=0)
		{
			found=0;
			edge1=FirstElementRegular->edge;
			ptId1=edgesvector[edge1].child;

			if (ptId1<0)
				goto failregular;

			this->MergeOutput->GetEdgeVertices(edge1,ptId0,ptId2);

			ptId0=this->vertices[ptId0].new_to_old;
			ptId2=this->vertices[ptId2].new_to_old;

			this->ConquerEdge(ptId0,ptId1,f1,ptId3);
			this->MergeInput->Conquer(f1,ptId1,ptId3,f2,ptId4);

			if (ptId4==ptId2)
				goto failregular;

			this->MergeInput->Conquer(f2,ptId4,ptId1,f3,ptId2);
			this->MergeInput->Conquer(f2,ptId3,ptId4,f4,ptId5);

			if (f4<0)
				goto failregular;

			if (faces[f4].hirestype>0)
				goto failregular;	

			if (((vertices[ptId3].type==PARENT)||(vertices[ptId4].type==PARENT))
				||(vertices[ptId5].type==CHILD))
				goto failregular;

			if (MatchParents(ptId3,ptId0,ptId5,f1,f4,3)==0)
				goto failregular;

			if (MatchParents(ptId4,ptId2,ptId5,f3,f4,3)==0)
				goto failregular;

			f5=AddFace(ptId0,ptId2,ptId5,edge1,edge2,edge3);

			faces[f5].lowrestype=4;

			faces[f1].hirestype=1;
			faces[f2].hirestype=1;
			faces[f3].hirestype=1;
			faces[f4].hirestype=1;

			edgesvector[edge1].child=ptId1;
			edgesvector[edge2].child=ptId3;
			edgesvector[edge3].child=ptId4;

			vertices[ptId1].type=CHILD;
			vertices[ptId3].type=CHILD;
			vertices[ptId4].type=CHILD;

			vertices[ptId0].type=PARENT;
			vertices[ptId2].type=PARENT;
			vertices[ptId5].type=PARENT;

			vertices[ptId1].parent1=ptId0;
			vertices[ptId1].parent2=ptId2;

			vertices[ptId3].parent1=ptId0;
			vertices[ptId3].parent2=ptId5;

			vertices[ptId4].parent1=ptId2;
			vertices[ptId4].parent2=ptId5;

			switchcontour2(edge1);
			switchcontour2(edge2);
			switchcontour2(edge3);

			this->M4to1++;

			this->remaining_faces-=4;

			found=1;

			if (this->DisplayEfficiency==1)
			{
				this->CellTypes->SetValue(f1,4);
				this->CellTypes->SetValue(f2,4);
				this->CellTypes->SetValue(f3,4);
				this->CellTypes->SetValue(f4,4);
			}

failregular:
			if (found==0)
			{
				next=this->edgesvector[edge1].contour->NextElement;
				previous=this->edgesvector[edge1].contour->PreviousElement;

				if ((next==0)&&(previous==0))
				{
					this->LastElementRegular=0;
					this->FirstElementRegular=0;
				}
				else
				{
					if (next==0)
					{
						previous->NextElement=0;
						this->LastElementRegular=previous;
					}
					else
					{
						if (previous==0)
						{
							this->FirstElementRegular=next;
							next->PreviousElement=0;
						}
						else
						{
							next->PreviousElement=previous;
							previous->NextElement=next;
						}
					}
				}

				if (this->LastElement==0)
				{
					this->FirstElement=edgesvector[edge1].contour;
					this->LastElement=edgesvector[edge1].contour;
					edgesvector[edge1].contour->NextElement=0;
					edgesvector[edge1].contour->PreviousElement=0;
				}
				else
				{
					this->edgesvector[edge1].contour->NextElement=0;
					this->edgesvector[edge1].contour->PreviousElement=this->LastElement;
					this->LastElement->NextElement=edgesvector[edge1].contour;
					this->LastElement=edgesvector[edge1].contour;
				}
			}

		}	//End of "while (this->FirstElementRegular!=0)"
	}// End of "while ((this->FirstElementRegular!=0))"
}
