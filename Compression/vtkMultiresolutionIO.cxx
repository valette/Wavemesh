/*=========================================================================

Program:   Mailleur 3D multi-résolution (Creatis 1997 ~)
Module:    vtkMultiresolutionIO.cxx
Language:  C++
Date:      2003/05
Auteurs:   Sebastien Valette
  This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME vtkMultiresolutionIO
// .SECTION Description

#include <list>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkIdTypeArray.h>
#include <vtkObjectFactory.h>
#include <vtkLookupTable.h>
#include <vtkWindowToImageFilter.h>
#include <vtkBMPWriter.h>
#include <vtkPLYWriter.h>
#include <vtkCamera.h>

#include "vtkSurface.h"
#include "vtkMultiresolutionIO.h"


void vtkMultiresolutionIO::PrintInsignificantCoeffs()
{
	int coord=0;
	std::list<WaveletCoefficient>::iterator CoeffIter;

	cout<<"InsignificantCoeffs:"<<endl;
	for (CoeffIter=this->InsignificantCoeffs[coord].begin();
		CoeffIter!=this->InsignificantCoeffs[coord].end();
		CoeffIter++)
	{
		cout<<"Filter "<<CoeffIter->FilterId<<", Wavelet "<<CoeffIter->WaveletId<<endl;
	}
}
void vtkMultiresolutionIO::PrintWavelets()
{
	int i,j,count=1;
	for	(i=0;i<this->NumberOfFilters;i++)
	{
		for (j=0;j<this->Filters[i]->Wavelets->GetNumberOfTuples();j++)
		{
			double	value1,value2,value3,Wav[3];
			this->Filters[i]->Wavelets->GetTuple(j,Wav);
			value1=Wav[0];
			value2=Wav[1];
			value3=Wav[2];
			cout<<count++<<" : Filter "<<i<<", wavelet "<<j<<" : "<<value1<<" "<<value2<<" "<<value3<<endl;
		}
	}
}
void vtkMultiresolutionIO::SaveWavelets()
{
	int i,j;
	fstream	WavF;
	WavF.open ("wavelets.dat", std::ofstream::out | std::ofstream::trunc|ios::binary);
	for	(i=0;i<this->NumberOfFilters;i++)
	{
		for (j=0;j<this->Filters[i]->Wavelets->GetNumberOfTuples();j++)
		{
			double	value1,value2,value3,Wav[3];
			this->Filters[i]->Wavelets->GetTuple(j,Wav);
			value1=Wav[0];
			value2=Wav[1];
			value3=Wav[2];
			WavF.write((char*) &value1,sizeof(double));
			WavF.write((char*) &value2,sizeof(double));
			WavF.write((char*) &value3,sizeof(double));
		}
	}
	WavF.close();
}

void vtkMultiresolutionIO::LoadWavelets()
{
	fstream WavF;
	WavF.open ("wavelets.dat", std::ofstream::in|ios::binary);
	int i,j;

	for	(i=0;i<this->NumberOfFilters;i++)
	{
		for (j=0;j<this->Filters[i]->Wavelets->GetNumberOfTuples();j++)
		{
			double	value1,value2,value3,Wav[3];
			WavF.read((char*) &value1,sizeof	(double));
			WavF.read((char*) &value2,sizeof	(double));
			WavF.read((char*) &value3,sizeof	(double));
			Wav[0]=value1;
			Wav[1]=value2;
			Wav[2]=value3;
			this->Filters[i]->Wavelets->SetTuple(j,Wav);
		}
	}
	WavF.close();
}

void vtkMultiresolutionIO::EncodeLiftingProperties()
{
	if (this->Lifting==0)
		this->ArithmeticCoder->EncodeByte(1);
	else
	{
		if (this->Lifting==2)
			this->ArithmeticCoder->EncodeByte(0);
		else
			this->ArithmeticCoder->EncodeByte(this->LiftingRadius+2);
	}
	this->ArithmeticCoder->EncodeBit(this->GeometryPrediction!=0);
}
void vtkMultiresolutionIO::DecodeLiftingProperties()
{
	int LiftingTest=this->ArithmeticCoder->DecodeByte();
	LiftingTest-=2;
	if (LiftingTest>=0)
	{
		this->Lifting=1;
		this->LiftingRadius=LiftingTest;
	}
	else
	{
		if (LiftingTest==-1)
		{
			this->Lifting=0;
			this->LiftingRadius=0;			
		}
		else
		{
			this->Lifting=2;
			this->LiftingRadius=0;
		}
	}
	this->GeometryPrediction=(this->ArithmeticCoder->DecodeBit()==true);

}

void vtkMultiresolutionIO::EncodeMeshConnectivity(vtkSurface *Mesh)
{
	vtkIdType i;
	vtkIdType n,v1,v2,v3;

	this->ArithmeticCoder->EncodeWord(Mesh->GetNumberOfPoints());
	this->ArithmeticCoder->EncodeWord(Mesh->GetNumberOfCells());
	qsmodel QsmConnectivity;
	QsmConnectivity.initqsmodel(Mesh->GetNumberOfPoints(),15,1000,NULL,1);

	n=Mesh->GetNumberOfCells();
	for (i=0;i<n;i++)
	{
		Mesh->GetFaceVertices(i,v1,v2,v3);
		this->ArithmeticCoder->Encode(v1,&QsmConnectivity);
		this->ArithmeticCoder->Encode(v2,&QsmConnectivity);
		this->ArithmeticCoder->Encode(v3,&QsmConnectivity);
	}
}
vtkSurface *vtkMultiresolutionIO::DecodeMeshConnectivity()
{
	vtkIdType i;
	int v1,v2,v3;
	int NumberOfPoints,NumberOfFaces;
	double P[3];
	P[0]=0;
	P[1]=0;
	P[2]=0;

	vtkSurface *Mesh=vtkSurface::New();
	NumberOfPoints=this->ArithmeticCoder->DecodeWord();
	NumberOfFaces=this->ArithmeticCoder->DecodeWord();
	Mesh->Init(NumberOfPoints,NumberOfFaces,NumberOfPoints+NumberOfFaces+1000);


	for (i=0;i<NumberOfPoints;i++)
		Mesh->AddVertex(P);

	qsmodel QsmConnectivity;
	QsmConnectivity.initqsmodel(NumberOfPoints,15,1000,NULL,0);


	for (i=0;i<NumberOfFaces;i++)
	{
		v1=this->ArithmeticCoder->Decode(&QsmConnectivity);
		v2=this->ArithmeticCoder->Decode(&QsmConnectivity);
		v3=this->ArithmeticCoder->Decode(&QsmConnectivity);
		Mesh->AddFace(v1,v2,v3);
	}

	vtkSurface *Mesh2=vtkSurface::New();
	Mesh2->CreateFromPolyData(Mesh);
	Mesh->Delete();
	return (Mesh2);


}

void vtkMultiresolutionIO::EncodeMeshGeometry(vtkSurface *Mesh)
{

	vtkIdType i,j;
	double P[3];
	int NumberOfPoints;

	NumberOfPoints=Mesh->GetNumberOfPoints();

	if (this->ArithmeticType==0)
	{
		for (i=0;i<NumberOfPoints;i++)
		{
			Mesh->GetPointCoordinates(i,P);
			for (j=0;j<3;j++)
			{
				this->ArithmeticCoder->EncodeFloat(P[j]);
			}
		}
	}
	else
	{
		int Value;
		for (i=0;i<NumberOfPoints;i++)
		{
			Mesh->GetPointCoordinates(i,P);
			for (j=0;j<3;j++)
			{
				Value=(int) P[j];
				this->ArithmeticCoder->EncodeWord(Value+16384);
			}
		}
	}
}
void vtkMultiresolutionIO::DecodeMeshGeometry(vtkSurface *Mesh)
{

	vtkIdType i,j;
	double P[3];
	int NumberOfPoints;

	NumberOfPoints=Mesh->GetNumberOfPoints();
	if (this->ArithmeticType==0)
	{
		for (i=0;i<NumberOfPoints;i++)
		{
			for (j=0;j<3;j++)
			{		
				P[j]=this->ArithmeticCoder->DecodeFloat();
			}
			Mesh->SetPointCoordinates(i,P);
		}
	}
	else
	{
		for (i=0;i<NumberOfPoints;i++)
		{
			for (j=0;j<3;j++)
			{		
				P[j]=(this->ArithmeticCoder->DecodeWord()-16384);
			}
			Mesh->SetPointCoordinates(i,P);
		}
	}
}

void vtkMultiresolutionIO::EncodeScalingFactors(vtkSurface *Mesh)
{
	double Factor, Tx, Ty, Tz;
	Mesh->GetScalingFactors(Factor,Tx,Ty,Tz);
	this->ArithmeticCoder->EncodeFloat((float) Factor);
	this->ArithmeticCoder->EncodeFloat((float) Tx);
	this->ArithmeticCoder->EncodeFloat((float) Ty);
	this->ArithmeticCoder->EncodeFloat((float) Tz);
}
void vtkMultiresolutionIO::DecodeScalingFactors(vtkSurface *Mesh)
{
	double Factor, Tx, Ty, Tz;
	Factor=this->ArithmeticCoder->DecodeFloat();
	Tx=this->ArithmeticCoder->DecodeFloat();
	Ty=this->ArithmeticCoder->DecodeFloat();
	Tz=this->ArithmeticCoder->DecodeFloat();
	Mesh->SetScalingFactors(Factor,Tx,Ty,Tz);
}


void vtkMultiresolutionIO::SetLifting(int lifting)
{
	int i;
	this->Lifting=lifting;

	if (this->NumberOfFilters>0)
		for (i=0;i<this->NumberOfFilters;i++)
			this->Filters[i]->SetLifting(lifting);
}
void vtkMultiresolutionIO::SetArithmeticType (int Type)
{
	int i;
	this->ArithmeticType=Type;
	if (this->NumberOfFilters>0)
		for (i=0;i<this->NumberOfFilters;i++)
			this->Filters[i]->SetArithmeticType(Type);
};

void vtkMultiresolutionIO::SetLiftingRadius(int radius)
{
	int i;
	this->LiftingRadius=radius;
	if (this->NumberOfFilters>0)
		for (i=0;i<this->NumberOfFilters;i++)
			this->Filters[i]->SetLiftingRadius(radius);
}


void vtkMultiresolutionIO::InsertNextFilter(vtkWaveletSubdivisionFilter *Filter)
{
	int i;

	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		this->Filters[i+1]=this->Filters[i];
		this->SynthesisMeshes[i+2]=this->Filters[i+1]->GetSubdivisionInput();
	}

	this->SynthesisMeshes[1]=Filter->GetSubdivisionInput();
	this->Filters[0]=Filter;
	this->SynthesisMeshes[0]=Filter->GetOutput();

	this->NumberOfFilters++;
	Filter->SetPointsIds(this->PointsIds);
}

void vtkMultiresolutionIO::EncodeSignificance(int S,int Coord)
{
	this->ArithmeticCoder->EncodeBit(S!=0);
	return;
	this->ArithmeticCoder->EncodeSignificance(S,Coord,0);
	this->SignificanceCodes++;
};

void vtkMultiresolutionIO::EncodeCoeffRefinement(int Coeff)
{
	this->ArithmeticCoder->EncodeSignOrRefinement(Coeff);
};

void vtkMultiresolutionIO::EncodeSign(int Sign)
{
	this->ArithmeticCoder->EncodeSignOrRefinement(Sign);
};

int vtkMultiresolutionIO::DecodeSignificance(int Coord)
{
	return (this->ArithmeticCoder->DecodeBit());
	int S; 
	S=this->ArithmeticCoder->DecodeSignificance(Coord,0);
	this->SignificanceCodes++;
	return (S);
};

int vtkMultiresolutionIO::DecodeCoeffRefinement()
{
	int Coeff; 
	Coeff=this->ArithmeticCoder->DecodeSignOrRefinement();
	return (Coeff);
};

int vtkMultiresolutionIO::DecodeSign()
{
	int Sign;
	Sign=this->ArithmeticCoder->DecodeSignOrRefinement();
	return (Sign);
};


void vtkMultiresolutionIO::GetFirstGoodEdge(vtkIdType Edge,vtkIdType FilterId, vtkIdType &Edge2, vtkIdType &FilterId2, vtkIdType &Vertex)
{

	Edge2=Edge;
	FilterId2=FilterId;

	Vertex=this->EdgeMidPoints[FilterId2]->GetId(Edge);
	if (Vertex>=0)
	{
		return;
	}
	while ((Vertex<0)&&(FilterId2>0))
	{
		Edge2=this->Filters[FilterId2]->TreeFirstChildEdge->GetId(Edge2);
		FilterId2--;
		Vertex=this->EdgeMidPoints[FilterId2]->GetId(Edge2);		
	}
	return;
}

void vtkMultiresolutionIO::ComputeWaveletTree ()
{
	int i,j,k,NumberOfEdges;
	double *Wavelet,*Max,*Max2;
	double AbsoluteWavelet;
	int NumberOfVertices;
	vtkIdType Vertex;
	for (i=0;i<this->NumberOfFilters;i++)
	{
		this->Filters[i]->TreeFirstChildEdge=this->Filters[i]->TreeFirstChildEdge;
		this->TreeNextChildEdge[i]=this->Filters[i]->TreeNextChildEdge;
		this->TreeParentEdge[i]=this->Filters[i]->TreeParentEdge;
		this->Wavelets[i]=this->Filters[i]->Wavelets;
		this->EdgeMidPoints[i]=this->Filters[i]->EdgeMidPoints;

		NumberOfEdges=this->Filters[i]->GetSubdivisionInput()->GetNumberOfEdges();
		this->Maxima[i]=vtkDoubleArray::New();
		this->Maxima[i]->SetNumberOfValues(NumberOfEdges*3);
		for (j=0;j<NumberOfEdges;j++)
		{
			Max=this->Maxima[i]->GetPointer(j*3);
			Max[0]=0;
			Max[1]=0;
			Max[2]=0;
		}
	}

	for (i=0;i<this->NumberOfFilters;i++)
	{
		NumberOfVertices=this->Filters[i]->GetSubdivisionInput()->GetNumberOfPoints();
		NumberOfEdges=this->Filters[i]->GetSubdivisionInput()->GetNumberOfEdges();


		if (i<this->NumberOfFilters-1)
		{
			for (j=0;j<NumberOfEdges;j++)
			{
				Vertex=this->EdgeMidPoints[i]->GetId(j)-NumberOfVertices;
				if (Vertex>=0)
				{
					Max=this->Maxima[i+1]->GetPointer(this->TreeParentEdge[i+1]->GetId(j)*3);
					Wavelet=this->Wavelets[i]->GetPointer(Vertex*3);
					for (k=0;k<3;k++)
					{
						AbsoluteWavelet=Wavelet[k];
						if (fabs(Max[k])<fabs(AbsoluteWavelet))
							Max[k]=AbsoluteWavelet;
					}
				}
			}
		}
		if (i>0)
		{
			NumberOfEdges=this->Filters[i]->GetOutput()->GetNumberOfEdges();	

			for (j=0;j<NumberOfEdges;j++)
			{
				Max=this->Maxima[i]->GetPointer(this->TreeParentEdge[i]->GetId(j)*3);
				Max2=this->Maxima[i-1]->GetPointer(j*3);
				for (k=0;k<3;k++)
				{
					if (fabs(Max[k])<fabs(Max2[k]))
						Max[k]=Max2[k];
				}
			}
		}
	}
}

void vtkMultiresolutionIO::EncodeProgressivePrecision()
{
	this->ComputeWaveletTree();

	vtkIdType i,j,BitPlane;
	WaveletCoefficientSet Set;
	WaveletCoefficient Coeff;

	double Max[3],*Wav,AbsoluteCoeff;
	vtkIdType Vertex,Edge;
	vtkIdType NumberOfVertices,FilterId;
	int Count=0; 

	for (i=0;i<3;i++)
	{
		this->InsignificantCoeffs[i].clear();
		this->InsignificantSets[i].clear();
		this->SignificantCoeffs[i].clear();
	}
	this->ArithmeticCoder->StartCoding();
	this->ArithmeticCoder->EncodeByte((this->NumberOfBitPlanes<<4)+this->NumberOfStartBitPlanes);
	this->ArithmeticCoder->InitConnectivityQSModels(1);

	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		if (this->Filters[i]->GetSubdivisionType()==0)
		{
			this->ArithmeticCoder->EncodeBit(0);
		}
		else
		{
			if (this->Filters[i]->GetSubdivisionType()==1)
			{
				this->ArithmeticCoder->EncodeBit(1);
				this->ArithmeticCoder->EncodeBit(0);
			}
			else
			{
				this->ArithmeticCoder->EncodeBit(1);
				this->ArithmeticCoder->EncodeBit(1);
			}
		}

		this->Filters[i]->SetIOType(1);
		this->Filters[i]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[i]->Subdivide();
	}
	
	this->ConnectivityBytes=this->ArithmeticCoder->StopCoding();
	this->ArithmeticCoder->StartCoding();

	Max[0]=0;
	Max[1]=0;
	Max[2]=0;

	NumberOfVertices=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints();

	for (i=0;i<this->SynthesisMeshes[this->NumberOfFilters]->GetNumberOfEdges();i++)
	{
		Wav=this->Maxima[this->NumberOfFilters-1]->GetPointer(i*3);
		for (j=0;j<3;j++)
		{
			AbsoluteCoeff=Wav[j];
			if (fabs(Max[j])<fabs(AbsoluteCoeff))
				Max[j]=AbsoluteCoeff;
		}
		this->GetFirstGoodEdge(i,this->NumberOfFilters-1,Coeff.Edge,Coeff.FilterId,Vertex);

		if (Vertex>=0)
		{
			Vertex-=this->Filters[Coeff.FilterId]->GetSubdivisionInput()->GetNumberOfPoints();
			Coeff.WaveletId=Vertex;
			// Add entry to LIP
			Wav=this->Wavelets[Coeff.FilterId]->GetPointer(Vertex*3);
			for (j=0;j<3;j++)
			{
				AbsoluteCoeff=Wav[j];
				if (fabs(Max[j])<fabs(AbsoluteCoeff))
					Max[j]=AbsoluteCoeff;
				Coeff.Wavelet=AbsoluteCoeff;
				Coeff.Count=Count++;
				Coeff.Sn=2;
				InsignificantCoeffs[j].push_back(Coeff);
			}

			// Add entry to LIS
			Set.Edge=Coeff.Edge;
			Set.FilterId=Coeff.FilterId;
			Set.Type=1;
			Set.Sn=2;
			Set.Count=Count++;

			Wav=this->Maxima[Set.FilterId]->GetPointer(Set.Edge*3);
			if (Set.FilterId>0)
				for (j=0;j<3;j++)
				{
					Set.Max=Wav[j];
					InsignificantSets[j].push_back(Set);
				}
		}
	}

	AbsoluteCoeff=Max[0];
	if (fabs(AbsoluteCoeff)<fabs(Max[1]))
		AbsoluteCoeff=Max[1];
	if (fabs(AbsoluteCoeff)<fabs(Max[2]))
		AbsoluteCoeff=Max[2];

	int N,Descendants,Sn;
	vtkIdType Edge2,FilterId2;
	double T,MaxCoeff,TestCoeff;
	std::list<WaveletCoefficientSet>::iterator SetIter;
	std::list<WaveletCoefficient>::iterator CoeffIter;
	std::list<WaveletCoefficientSet>::iterator SetIter2;
	std::list<WaveletCoefficient>::iterator CoeffIter2;

	N=(int) floor(log(fabs(AbsoluteCoeff))/log(2.00));
	for (i=0;i<3;i++)
	{
		for (CoeffIter=InsignificantCoeffs[i].begin();
			CoeffIter!=InsignificantCoeffs[i].end();CoeffIter++)
			CoeffIter->N=N;
	}

	for (i=0;i<this->NumberOfFilters;i++)
	{
		this->SauvWavelets[i]=vtkDoubleArray::New();
		this->SauvWavelets[i]->DeepCopy(Wavelets[i]);
	}

	this->ArithmeticCoder->EncodeByte(N+128);
	this->ArithmeticCoder->EncodeByte(this->NumberOfBitPlanes);
	this->ArithmeticCoder->EncodeByte(this->NumberOfStartBitPlanes);
	T=pow(2.0,N);

	for (BitPlane=0;BitPlane<this->NumberOfBitPlanes;BitPlane++)
	{
		if (BitPlane>=this->NumberOfStartBitPlanes)
			this->ArithmeticCoder->StartCoding();
		this->ArithmeticCoder->InitZerotreeQSModels(1);
		for (i=0;i<3;i++)
		{
			int number=0;

			//2.1 For Each Entry (i,j) in the LIP
			for (CoeffIter=InsignificantCoeffs[i].begin();
				CoeffIter!=InsignificantCoeffs[i].end();)
			{
				number++;

				Sn=(fabs(CoeffIter->Wavelet)>=T);
				this->EncodeSignificance(Sn,i);

				//2.1.1 Output Sn(i,j)
				if (Sn==0)
					CoeffIter++;
				else
				{
					// If Sn(i,j)=1 then move (i,j) to the LSP and output the sign of cij
					Sn=(CoeffIter->Wavelet<0);
					this->EncodeSign(Sn);

					Coeff.Edge=CoeffIter->Edge;
					Coeff.FilterId=CoeffIter->FilterId;
					Coeff.Wavelet=CoeffIter->Wavelet;
					Coeff.WaveletId=CoeffIter->WaveletId;
					Coeff.N=N;
					SignificantCoeffs[i].push_back(Coeff);

					CoeffIter=InsignificantCoeffs[i].erase(CoeffIter);
				}
			}

			//2.2 for each entry (i,j) in the LIS do:
			for (SetIter=InsignificantSets[i].begin();
				SetIter!=InsignificantSets[i].end();)
			{
				if (SetIter->Type==1)
				{
					// 2.2.1 if the entry is of type A
					// Output Sn(D(i,j))

					Sn=(fabs(SetIter->Max)>=T);
					this->EncodeSignificance(Sn,i);

					if (Sn==0)
						SetIter++;
					else
					{
						// is Sn(i,j)=1 then:
						Descendants=0;
						MaxCoeff=0;
						FilterId=SetIter->FilterId-1;

						// for each (k,l) in O(i,j)
						Edge=this->Filters[SetIter->FilterId]->TreeFirstChildEdge->GetId(SetIter->Edge);
						while (Edge>=0)
						{
							this->GetFirstGoodEdge(Edge,FilterId,Edge2,FilterId2,Vertex);
							TestCoeff=this->Maxima[FilterId2]->GetPointer(Edge2*3)[i];
							if (fabs(MaxCoeff)<fabs(TestCoeff))
								MaxCoeff=TestCoeff;

							if (Vertex>=0)
							{
								if (FilterId2>0)
									Descendants++;
								Coeff.WaveletId=Vertex-this->Filters[FilterId2]->GetSubdivisionInput()->GetNumberOfPoints();
								Coeff.Edge=Edge2;
								Coeff.FilterId=FilterId2;
								Coeff.Wavelet=this->Wavelets[FilterId2]->GetPointer(Coeff.WaveletId*3)[i];
								Coeff.N=N;
								Coeff.Count=Count;
								Coeff.Sn=2;

								//Output Sn(k,l)
								Sn=(fabs(Coeff.Wavelet)>=T);

								this->EncodeSignificance(Sn,i);
								if (Sn==0)
									// if Sn(k,l)=0 then add (k,l) to the end of the LIP
									InsignificantCoeffs[i].push_back(Coeff);
								else
								{
									// if Sn(k,l)=1 then add (k,l) to the end of the LSP and output the sign of cij
									Sn=Coeff.Wavelet<0;
									this->EncodeSign(Sn);
									SignificantCoeffs[i].push_back(Coeff);
								}

							}
							if (FilterId>=0)
								Edge=this->TreeNextChildEdge[SetIter->FilterId]->GetId(Edge);
							else
								Edge=-1;
						}
						Count++;
						// if L(i,j) non empty
						if ((Descendants>0)&&(FilterId>0))
						{

							SetIter->Max=MaxCoeff;
							SetIter->Type=0;
						}
						else
						{
							// else remove (i,j) from the LIS
							SetIter=InsignificantSets[i].erase(SetIter);
						}
					}
				}

				// 2.2.2 if the entry is of type B then
				if (SetIter!=InsignificantSets[i].end())
				{
					if (SetIter->Type==0)
					{
						Sn=(fabs(SetIter->Max)>=T);
						// Output Sn(L(i,j))
						this->EncodeSignificance(Sn,i);

						if (Sn==0)
							SetIter++;
						else
						{
							// if Sn(L(i,j))=1 then
							MaxCoeff=0;
							Edge=this->Filters[SetIter->FilterId]->TreeFirstChildEdge->GetId(SetIter->Edge);
							FilterId=SetIter->FilterId-1;
							// add each (k,l) In O(i,j) to the end of the LIS as an entry of type A
							while (Edge>=0)
							{
								this->GetFirstGoodEdge(Edge,FilterId,Edge2,FilterId2,Vertex);
								if ((Vertex>=0)&&(FilterId2>0))
								{
									Set.Edge=Edge2;
									Set.FilterId=FilterId2;
									Set.Max=this->Maxima[FilterId2]->GetPointer(Edge2*3)[i];
									Set.Type=1;
									Set.Sn=2;
									Set.Count=Count;
									InsignificantSets[i].push_back(Set);
								}
								Edge=this->TreeNextChildEdge[SetIter->FilterId]->GetId(Edge);
							}
							// remove (i,j) from the LIS;
							SetIter=InsignificantSets[i].erase(SetIter);
							Count++;
						}
					}
				}
			}		

			// Refinement pass
			int newnumber=0,k,Ptest;
			number=0;
			for(CoeffIter=SignificantCoeffs[i].begin();
				CoeffIter!=SignificantCoeffs[i].end();CoeffIter++)
			{
				number++;
				if (CoeffIter->N!=N)
				{
					newnumber++;
					TestCoeff=floor(fabs(CoeffIter->Wavelet)/T);
					Ptest=(int) TestCoeff;
					if (Ptest%2==1)
					{
						this->EncodeCoeffRefinement(1);
					}
					else
					{
						this->EncodeCoeffRefinement(0);
					}

				}
			}

			number=0;
			for (k=0;k<this->NumberOfFilters;k++)
			{
				for (j=0;j<this->Filters[k]->GetOutput()->GetNumberOfPoints()
					-this->Filters[k]->GetSubdivisionInput()->GetNumberOfPoints();j++)
				{
					Wav=this->Wavelets[k]->GetPointer(j*3);
					if (fabs(Wav[i])>=T)
						number++;
				}
			}

		}
		N--;
		T=T/2;
		if (BitPlane>=this->NumberOfStartBitPlanes-1)
			this->DataSize[BitPlane]=this->ArithmeticCoder->StopCoding();
		else
			this->DataSize[BitPlane]=0;
	}

	this->ArithmeticCoder->CloseFile();
	if (this->WriteRepport==1)
	{

		std::ofstream Repport;
		Repport.open ("repport.txt", std::ofstream::out | std::ofstream::trunc);

		double duration,ratio;
#if ( (VTK_MAJOR_VERSION >= 5))
		duration = this->Timer->GetUniversalTime()- this->StartTime;
#else
		duration = this->Timer->GetCurrentTime()- this->StartTime;
#endif
		ratio= ((double)this->Filters[0]->GetOutput()->GetNumberOfCells())/duration;

		int Sum;
		int NumberOfVertices=this->SynthesisMeshes[0]->GetNumberOfPoints();
		Repport<<"Zerotree Encoding of a mesh with "<<this->SynthesisMeshes[0]->GetNumberOfCells()<<
			" faces and "<<NumberOfVertices<<" vertices, "<<this->NumberOfFilters<<" resolution levels"<<endl;

		Repport<<"Coding of mesh connectivity+ base mesh geometry : "<<this->ConnectivityBytes<<
			" bytes ( "<<this->ConnectivityBytes*8<<" bits, "<<
			8.0*(double)this->ConnectivityBytes/(double)NumberOfVertices<<" bits/vertex)"<<endl;
		Repport<<"Encoding time : "<<duration<<" seconds ( "<<ratio<<" faces/s)"<<endl;

		Sum=this->ConnectivityBytes+this->DataSize[this->NumberOfStartBitPlanes-1];

		Repport<<"Bitplanes 0 to "<<this->NumberOfStartBitPlanes-1<<" : "
			<<this->DataSize[this->NumberOfStartBitPlanes-1]<<" Bytes, Total Data : "<<
			Sum<<" Bytes ( "<<Sum*8<<" bits, "<<8.0*(double)Sum/(double)NumberOfVertices<<" bits/Vertex)"<<endl;

		for (BitPlane=this->NumberOfStartBitPlanes;BitPlane<this->NumberOfBitPlanes;BitPlane++)
		{
			Sum+=this->DataSize[BitPlane];
			Repport<<"Bitplane "<<BitPlane<<" : "<<this->DataSize[BitPlane]<<" Bytes, Total Data : "<<
				Sum<<" Bytes ( "<<Sum*8<<" bits, "<<8.0*(double)Sum/(double)NumberOfVertices<<" bits/Vertex)"<<endl;
		}
		Repport.close();
	}
}

void vtkMultiresolutionIO::DecodeProgressivePrecision()
{
	int i,j,BitPlane;
	WaveletCoefficientSet Set;
	WaveletCoefficient Coeff;

	vtkIdType Vertex,Edge;
	int NumberOfVertices,FilterId,Count=0;

	for (i=0;i<3;i++)
	{
		this->InsignificantCoeffs[i].clear();
		this->InsignificantSets[i].clear();
		this->SignificantCoeffs[i].clear();
	}

	this->ArithmeticCoder->StartDecoding();

	int SubdivisionType;
	int BitPlanesCode=this->ArithmeticCoder->DecodeByte();
	this->NumberOfBitPlanes=BitPlanesCode>>4;
	this->NumberOfStartBitPlanes=BitPlanesCode-(this->NumberOfBitPlanes<<4);

	this->ArithmeticCoder->InitConnectivityQSModels(0);

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		if (this->ArithmeticCoder->DecodeBit()==0)
			SubdivisionType=0;
		else
		{
			if (this->ArithmeticCoder->DecodeBit()==0)
				SubdivisionType=1;
			else
				SubdivisionType=2;
		}
		
		this->Filters[j]=this->NewFilter(SubdivisionType);

		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);

		this->Filters[j]->SetSubdivisionType(SubdivisionType);
		this->Filters[j]->SetIOType(2);	
		this->Filters[j]->SetLifting(this->Lifting);
		this->Filters[j]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[j]->GeometryPrediction=this->GeometryPrediction;
		this->Filters[j]->SetArithmeticType(0);
		this->EdgeMidPoints[j]=this->Filters[j]->EdgeMidPoints;
		this->TreeNextChildEdge[j]=this->Filters[j]->TreeNextChildEdge;

		this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[j]->Subdivide();
		this->Filters[j]->Wavelets=vtkDoubleArray::New();
		this->Filters[j]->Wavelets->SetNumberOfComponents(3);
		this->Filters[j]->Wavelets->SetNumberOfTuples(this->Filters[j]->GetOutput()->GetNumberOfPoints()-
			this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints());

		this->Wavelets[j]=this->Filters[j]->Wavelets;
		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();

		if (this->DisplayText)
			cout<<"Level "<<this->NumberOfFilters-j<<": "<<this->SynthesisMeshes[j]->GetNumberOfCells()
			<<" faces, "<<this->SynthesisMeshes[j]->GetNumberOfPoints()<<" vertices "<<endl;

	}

	int k;
	double Wav[3];
	for (k=0;k<this->NumberOfFilters;k++)
	{
		for (j=0;j<this->Filters[k]->GetOutput()->GetNumberOfPoints()
			-this->Filters[k]->GetSubdivisionInput()->GetNumberOfPoints();j++)
		{
			this->Wavelets[k]->GetTuple(j,Wav);
			Wav[0]=0;
			Wav[1]=0;
			Wav[2]=0;
			this->Wavelets[k]->SetTuple(j,Wav);
		}
	}
	this->Reconstruct();
	if (this->DisplayText)
		cout<<"Final mesh: "<<this->SynthesisMeshes[0]->GetNumberOfCells()<<" faces"<<endl;
	if (this->Display!=0)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
		this->MeshWindow->Render();
		this->MeshWindow->SetWindowName("Progressive precision reconstruction");
		cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
		this->MeshWindow->Interact();

	}
	this->ArithmeticCoder->StopDecoding();

	this->ArithmeticCoder->StartDecoding();

	NumberOfVertices=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints();
	for (i=0;i<this->SynthesisMeshes[this->NumberOfFilters]->GetNumberOfEdges();i++)
	{
		this->GetFirstGoodEdge(i,this->NumberOfFilters-1,Coeff.Edge,Coeff.FilterId,Vertex);
		if (Vertex>=0)
		{
			Vertex-=this->Filters[Coeff.FilterId]->GetSubdivisionInput()->GetNumberOfPoints();
			Coeff.WaveletId=Vertex;
			Coeff.Count=Count++;
			Coeff.Sn=2;
			// Add entry to LIP
			for (j=0;j<3;j++)
				InsignificantCoeffs[j].push_back(Coeff);

			// Add entry to LIS
			Set.Edge=Coeff.Edge;
			Set.FilterId=Coeff.FilterId;
			Set.Type=1;
			Set.Count=Count++;
			Set.Sn=2;

			if (Set.FilterId>0)
				for (j=0;j<3;j++)
					InsignificantSets[j].push_back(Set);
		}
		else
		{
//			cout<<"Bizarre!!!"<<endl;
		}

	}


	int N,Descendants,Sn;
	vtkIdType Edge2,FilterId2;
	double T,MaxCoeff;
	std::list<WaveletCoefficientSet>::iterator SetIter;
	std::list<WaveletCoefficient>::iterator CoeffIter;

	N=this->ArithmeticCoder->DecodeByte()-128;
	this->NumberOfBitPlanes=this->ArithmeticCoder->DecodeByte();
	this->NumberOfStartBitPlanes=this->ArithmeticCoder->DecodeByte();

	T=pow(2.0,N);
	for (i=0;i<3;i++)
	{
		for (CoeffIter=InsignificantCoeffs[i].begin();
			CoeffIter!=InsignificantCoeffs[i].end();CoeffIter++)
			CoeffIter->N=N;
	}


	for (BitPlane=0;BitPlane<this->NumberOfBitPlanes;BitPlane++)
	{
		if (BitPlane>=this->NumberOfStartBitPlanes)
			this->ArithmeticCoder->StartDecoding();
		this->ArithmeticCoder->InitZerotreeQSModels(0);

		for (i=0;i<3;i++)
		{
			int number=0;
			//2.1 For Each Entry (i,j) in the LIP
			for (CoeffIter=InsignificantCoeffs[i].begin();
				CoeffIter!=InsignificantCoeffs[i].end();)
			{
				number++;
				//2.1.1 Input Sn(i,j)
				Sn=this->DecodeSignificance(i);

				if (Sn==0)
					CoeffIter++;
				else
				{

					// If Sn(i,j)=1 then move (i,j) to the LSP and input the sign of cij
					Sn=this->DecodeSign();
					if (Sn==1)
						this->Wavelets[CoeffIter->FilterId]->GetPointer(CoeffIter->WaveletId*3)[i]=-1.5*pow(2.0,N);
					else
						this->Wavelets[CoeffIter->FilterId]->GetPointer(CoeffIter->WaveletId*3)[i]=1.5*pow(2.0,N);

					Coeff.WaveletId=CoeffIter->WaveletId;
					Coeff.Edge=CoeffIter->Edge;
					Coeff.FilterId=CoeffIter->FilterId;
					Coeff.N=N;
					SignificantCoeffs[i].push_back(Coeff);
					CoeffIter=InsignificantCoeffs[i].erase(CoeffIter);
				}
			}

			//2.2 for each entry (i,j) in the LIS do:
			for (SetIter=InsignificantSets[i].begin();
				SetIter!=InsignificantSets[i].end();)
			{
				if (SetIter->Type==1)
				{
					// 2.2.1 if the entry is of type A
					// Output Sn(D(i,j))
					Sn=this->DecodeSignificance(i);


					if (Sn==0)
						SetIter++;
					else
					{
						// is Sn(i,j)=1 then:
						Descendants=0;
						MaxCoeff=0;
						FilterId=SetIter->FilterId-1;


						// for each (k,l) in O(i,j)
						Edge=this->Filters[SetIter->FilterId]->TreeFirstChildEdge->GetId(SetIter->Edge);
						while (Edge>=0)
						{
							this->GetFirstGoodEdge(Edge,FilterId,Edge2,FilterId2,Vertex);
							if (Vertex>=0)
							{
								if (FilterId2>0)
									Descendants++;

								Coeff.WaveletId=Vertex-this->Filters[FilterId2]->GetSubdivisionInput()->GetNumberOfPoints();
								Coeff.Edge=Edge2;
								Coeff.FilterId=FilterId2;
								Coeff.N=N;


								//Input Sn(k,l)
								Sn=this->DecodeSignificance(i);

								if (Sn==0)
								{
									// if Sn(k,l)=0 then add (k,l) to the end of the LIP
									InsignificantCoeffs[i].push_back(Coeff);
								}
								else
								{
									// if Sn(k,l)=1 then add (k,l) to the end of the LSP and output the sign of cij
									Sn=this->DecodeSign();
									if (Sn==1)
										this->Wavelets[Coeff.FilterId]->GetPointer(Coeff.WaveletId*3)[i]=-1.5*pow(2.0,N);
									else
										this->Wavelets[Coeff.FilterId]->GetPointer(Coeff.WaveletId*3)[i]=1.5*pow(2.0,N);

									SignificantCoeffs[i].push_back(Coeff);
								}
							}
							if (FilterId>=0)
								Edge=this->TreeNextChildEdge[SetIter->FilterId]->GetId(Edge);
							else
								Edge=-1;
						}
						// if L(i,j) non empty
						if ((Descendants>0)&&(FilterId>0))
						{
							SetIter->Type=0;
						}
						else
						{
							// else remove (i,j) from the LIS
							SetIter=InsignificantSets[i].erase(SetIter);
						}
					}
				}

				// 2.2.2 if the entry is of type B then
				if (SetIter!=InsignificantSets[i].end())
				{
					if (SetIter->Type==0)
					{
						Sn=this->DecodeSignificance(i);

						// Input Sn(L(i,j))
						if (Sn==0)
							SetIter++;
						else
						{
							// if Sn(L(i,j))=1 then
							MaxCoeff=0;
							Edge=this->Filters[SetIter->FilterId]->TreeFirstChildEdge->GetId(SetIter->Edge);
							FilterId=SetIter->FilterId-1;
							// add each (k,l) In O(i,j) to the end of the LIS as an entry of type A
							while (Edge>=0)
							{
								this->GetFirstGoodEdge(Edge,FilterId,Edge2,FilterId2,Vertex);
								if ((Vertex>=0)&&(FilterId2>0))
								{
									Set.Edge=Edge2;
									Set.FilterId=FilterId2;
									Set.Type=1;
									InsignificantSets[i].push_back(Set);
								}
								Edge=this->TreeNextChildEdge[SetIter->FilterId]->GetId(Edge);
							}
							// remove (i,j) from the LIS;
							SetIter=InsignificantSets[i].erase(SetIter);
						}
					}
				}
			}		

			// Refinement pass
			int newnumber=0,ref;
			double *WaveletPointer;
			number=0;
			for(CoeffIter=SignificantCoeffs[i].begin();
				CoeffIter!=SignificantCoeffs[i].end();CoeffIter++)
			{
				number++;
				if (CoeffIter->N!=N)
				{
					newnumber++;
					WaveletPointer=this->Wavelets[CoeffIter->FilterId]->GetPointer(CoeffIter->WaveletId*3);
					ref=this->DecodeCoeffRefinement();

					if (WaveletPointer[i]<0)
					{
						if (ref==1)
							WaveletPointer[i]-=pow(2.0, N-1);
						else
							WaveletPointer[i]+=pow(2.0, N-1);
					}
					else
					{
						if (ref==1)
							WaveletPointer[i]+=pow(2.0, N-1);
						else
							WaveletPointer[i]-=pow(2.0,N-1);
					}
				}
			}
		}

		N--;

		cout<<"Bitplane :"<<BitPlane<<", Threshold="<<T<<endl;
		T=T/2;

		this->SynthesisMeshes[0]->Modified();
		if (BitPlane>=this->NumberOfStartBitPlanes-1)
		{
			this->ArithmeticCoder->StopDecoding();
			this->Reconstruct();
			if (this->Display!=0)
			{
				this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
				this->MeshWindow->Render();
			}
			if ((this->Display==2)||((BitPlane==this->NumberOfStartBitPlanes-1)&&(this->Display!=0)))
				this->MeshWindow->Interact();
			if (this->WriteOutput==1)
			{
				if (this->FileType==0)
				{
					std::stringstream strfile;
					strfile<<"Mesh"<<BitPlane<<".iv";
					this->Filters[0]->GetOutput()->WriteInventor(strfile.str().c_str());
				}
				else
				{
					std::stringstream strfile;
					strfile<<"Mesh"<<BitPlane<<".ply";
					vtkPLYWriter *Writer=vtkPLYWriter::New();
					Writer->SetInputData(this->Filters[0]->GetOutput());
					Writer->SetFileName(strfile.str().c_str());
					Writer->Write();
					Writer->Delete();
				}
			}
			if (this->Capture==1)
			{
				std::stringstream strfile;
				strfile<<"Mesh"<<BitPlane<<".bmp";
				this->MeshWindow->Capture(strfile.str().c_str());
			}
		}
	}

	if (this->Display!=0)
	{
		this->MeshWindow->Render();
		this->MeshWindow->Interact();
	}
	this->ArithmeticCoder->CloseFile();
	this->Output=this->SynthesisMeshes[0];
}

void vtkMultiresolutionIO::Execute()
{

}

void vtkMultiresolutionIO::Approximate()
{
	int i;
	for (i=0;i<this->NumberOfFilters;i++)
	{
		this->Filters[i]->SetGeometryPrediction(this->GeometryPrediction);
		this->Filters[i]->SetLifting(this->Lifting);
		this->Filters[i]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[i]->Approximate();

	}
}

void vtkMultiresolutionIO::Reconstruct()
{
	int i;
	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		this->Filters[i]->SetGeometryPrediction(this->GeometryPrediction);
		this->Filters[i]->SetLifting(this->Lifting);
		this->Filters[i]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[i]->Reconstruct();
		this->SynthesisMeshes[i]=this->Filters[i]->GetOutput();
	}
	this->SynthesisMeshes[0]->Modified();
	this->SynthesisMeshes[0]->GetPoints()->Modified();
	this->Output=this->SynthesisMeshes[0];
	this->Output->Register(this);
}
void vtkMultiresolutionIO::Read()
{
	cout<<"Reconstruction:"<<endl;

	this->ArithmeticCoder=vtkArithmeticCoder::New();
	this->ArithmeticCoder->OpenFile(this->FileName,0);
	this->ArithmeticCoder->InitConnectivityQSModels(0);
	this->ArithmeticCoder->StartDecoding();
	this->ArithmeticType=this->ArithmeticCoder->DecodeBit();

	vtkSurface *BaseMesh=this->DecodeMeshConnectivity();
	this->DecodeMeshGeometry(BaseMesh);
	this->DecodeLiftingProperties();
	this->NumberOfFilters=this->ArithmeticCoder->DecodeByte();

	int Type,i;
	int bit1,bit2;
	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		bit1=this->ArithmeticCoder->DecodeBit();
		bit2=this->ArithmeticCoder->DecodeBit();
		Type=bit1*2+bit2;

		this->Filters[i]=this->NewFilter(Type);
		this->Filters[i]->SetSubdivisionType(Type);
	}


	this->ArithmeticCoder->StopDecoding();


	if (this->DisplayText)
		cout<<"Level 0 : "<<BaseMesh->GetNumberOfCells()
		<<" faces, "<<BaseMesh->GetNumberOfPoints()<<" vertices and "
		<<BaseMesh->GetNumberOfEdges()<<" edges "<<endl;

	this->SynthesisMeshes[this->NumberOfFilters]=BaseMesh;
	if (this->ArithmeticType!=0)
        this->DecodeProgressiveResolution();
	else
		this->DecodeProgressivePrecision();
}

void vtkMultiresolutionIO::DecodeProgressiveResolution()
{
	int j;
	double start2,finish2;

	if (this->Display!=0)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[this->NumberOfFilters]);
		this->MeshWindow->Render();
		this->MeshWindow->SetWindowName("Progressive resolution reconstruction");
#if ( (VTK_MAJOR_VERSION >= 5))
		start2=this->Timer->GetUniversalTime();
		finish2=this->Timer->GetUniversalTime();
		while (finish2 - start2<this->Time)
			finish2=this->Timer->GetUniversalTime();

#else
		start2=this->Timer->GetCurrentTime();
		finish2=this->Timer->GetCurrentTime();
		while (finish2 - start2<this->Time)
			finish2=this->Timer->GetCurrentTime();
#endif			
		this->MeshWindow->Interact();
	}

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);
		this->Filters[j]->SetIOType(2);	
		this->Filters[j]->SetLifting(this->Lifting);
		this->Filters[j]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[j]->SetArithmeticType(1);
		this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[j]->Quantization=this->Quantization;

		this->ArithmeticCoder->StartDecoding();
		if (j==this->NumberOfFilters-1)
			this->DecodeScalingFactors(this->Filters[j]->GetSubdivisionInput());


		this->Filters[j]->Subdivide();
		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();
		this->ArithmeticCoder->StopDecoding();
		this->ArithmeticCoder->StartDecoding();
		this->Filters[j]->ReadCoefficients();
		this->ArithmeticCoder->StopDecoding();
		this->Filters[j]->Reconstruct();

		if (this->DisplayText)
			cout<<"Level "<<this->NumberOfFilters-j<<": "<<this->SynthesisMeshes[j]->GetNumberOfCells()
			<<" faces, "<<this->SynthesisMeshes[j]->GetNumberOfPoints()<<" vertices and "
			<<this->SynthesisMeshes[j]->GetNumberOfEdges()<<" edges "<<endl;


		if (this->Display!=0)
		{
			this->MeshWindow->SetInputData(this->Filters[j]->GetOutput());
			this->MeshWindow->Render();

#if ( (VTK_MAJOR_VERSION >= 5))
			start2=this->Timer->GetUniversalTime();
			finish2=this->Timer->GetUniversalTime();
			while (finish2 - start2<this->Time)
				finish2=this->Timer->GetUniversalTime();

#else
			start2=this->Timer->GetCurrentTime();
			finish2=this->Timer->GetCurrentTime();
			while (finish2 - start2<this->Time)
				finish2=this->Timer->GetCurrentTime();
#endif
			if (this->Capture==1)
			{
				std::stringstream strfile;
				strfile<<"Mesh"<<this->NumberOfFilters-j<<".bmp";
				this->MeshWindow->Capture(strfile.str().c_str());
			}

			if (this->Display==2)
			{
				cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
				this->MeshWindow->Interact();
			}
		}

	}

	if (this->DisplayText)
		cout<<"Final mesh: "<<this->SynthesisMeshes[0]->GetNumberOfCells()<<" faces"<<endl;
	this->ArithmeticCoder->CloseFile();

	this->SynthesisMeshes[0]->Modified();
	if (this->Display!=0)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
		this->MeshWindow->Render();
		cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
		this->MeshWindow->Interact();
	}
	if (this->WriteOutput==1)
	{
		this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->UnQuantizeCoordinates();
		if (this->FileType==0)
			this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->WriteInventor("Mesh0.iv");
		else
		{
			vtkPLYWriter *Writer=vtkPLYWriter::New();
			Writer->SetFileName("Mesh0.ply");
			Writer->SetInputData(this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput());
			Writer->Write();
			Writer->Delete();
		}
		double Factor,Tx,Ty,Tz;
		this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetScalingFactors(Factor,Tx,Ty,Tz);
		for (j=0;j<this->NumberOfFilters;j++)
		{
			std::stringstream strfile;
			this->Filters[j]->GetOutput()->SetScalingFactors(Factor,Tx,Ty,Tz);
			this->Filters[j]->GetOutput()->UnQuantizeCoordinates();
			if (this->FileType==0)
			{
				strfile<<"Mesh"<<this->NumberOfFilters-j<<".iv";
				this->Filters[j]->GetOutput()->WriteInventor(strfile.str().c_str());
			}
			else
			{
				strfile<<"Mesh"<<this->NumberOfFilters-j<<".ply";
				vtkPLYWriter *Writer=vtkPLYWriter::New();
				Writer->SetFileName(strfile.str().c_str());
				Writer->SetInputData(this->Filters[j]->GetOutput());
				Writer->Write();
				Writer->Delete();
			}
		}
	}
}
void vtkMultiresolutionIO::Write()
{
	vtkIdType Id;
	vtkSurface *BaseMesh;

	this->PointsIds->SetNumberOfIds(this->Filters[0]->GetOutput()->GetNumberOfPoints());
	for (Id=0;Id<this->Filters[0]->GetOutput()->GetNumberOfPoints();Id++)
	{
		this->PointsIds->SetId(Id,Id);
	}
	BaseMesh=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput();

	this->ArithmeticCoder=vtkArithmeticCoder::New();
	this->ArithmeticCoder->OpenFile(this->FileName,1);
	this->ArithmeticCoder->InitConnectivityQSModels(1);
	this->ArithmeticCoder->StartCoding();
	this->ArithmeticCoder->EncodeBit(this->ArithmeticType!=0);
	this->EncodeMeshConnectivity(BaseMesh);
	this->EncodeMeshGeometry(BaseMesh);
	this->EncodeLiftingProperties();
	this->ArithmeticCoder->EncodeByte(this->NumberOfFilters);

	int Type,i;
	int bit1,bit2;
	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		Type=this->Filters[i]->GetSubdivisionType();
		bit1=Type&2;
		bit2=Type&1;
		this->ArithmeticCoder->EncodeBit(bit1!=0);
		this->ArithmeticCoder->EncodeBit(bit2!=0);
	}


	this->BaseMeshBitRate=this->ArithmeticCoder->StopCoding()*8;
	if (this->ArithmeticType==1)
		this->EncodeProgressiveResolution();
	else
		this->EncodeProgressivePrecision();
}
void vtkMultiresolutionIO::EncodeProgressiveResolution()
{
	int j;
	int GeometryBitrate[1000];
	int ConnectivityBitrate[1000];


	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		this->ArithmeticCoder->StartCoding();
		this->Filters[j]->SetIOType(1);
		this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[j]->SetPointsIds(this->PointsIds);
		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);

		
		if (j==this->NumberOfFilters-1)
		{
			if (this->Filters[0]->GetMergeInput()!=0)
                this->EncodeScalingFactors(this->Filters[0]->GetMergeInput());
			else
				this->EncodeScalingFactors(this->Filters[0]->GetOutput());

		}

		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();

		this->Filters[j]->Subdivide();
		ConnectivityBitrate[j]=this->ArithmeticCoder->StopCoding()*8;
		this->ArithmeticCoder->StartCoding();
		this->Filters[j]->WriteCoefficients();
		this->Filters[j]->Reconstruct();	
		GeometryBitrate[j]=this->ArithmeticCoder->StopCoding()*8;
	}
	this->ArithmeticCoder->CloseFile();

	if (this->WriteRepport==1)
	{
		int connectivitybits=0;
		double encodedbits=0,duration,ratio;

#if ( (VTK_MAJOR_VERSION >= 5))
		duration = this->Timer->GetUniversalTime()- this->StartTime;
#else
		duration = this->Timer->GetCurrentTime()- this->StartTime;
#endif
		ratio= (this->Filters[0]->GetOutput()->GetNumberOfCells())/duration;

		std::ofstream Repport;
		Repport.open ("repport.txt", std::ofstream::out | std::ofstream::trunc);
		Repport<<"Filename: "<<this->FileName<<endl;

		Repport<<"Quantization : "<<this->Quantization<<" bits"<<endl;
		if (this->Lifting==0)
			Repport<<"No lifting"<<endl;
		else
		{
			if (this->Lifting==2)
				Repport<<"Fast 0-ring Lifting"<<endl;
			else
				Repport<<"Lifting radius: "<<this->LiftingRadius<<endl;
		}

		if (this->GeometricConstraint==0)
			Repport<<"No Wavelet Geometrical Criterion"<<endl;
		else
		{
			Repport<<"Geometry threshold: "<<this->EdgeAngleThreshold<<endl;
			Repport<<"Wavelet threshold: "<<this->WGC<<endl;
		}
		Repport<<"Total execution time : "<<duration<<" seconds : "<<ratio<<" faces/s"<<endl;


		double numberoffaces=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfCells();
		double numberofvertices=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints();
		double cost=32.0+ceil(log(numberofvertices)/log(2.0))*3.0*numberoffaces;
		connectivitybits=(int) cost;
		encodedbits=this->BaseMeshBitRate;

		Repport<<"Level 0 : "<<this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfCells()
			<<"f, "<<this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints()
			<<"v, total data: "<<(int) this->BaseMeshBitRate<<" bits (connectivity: "<<connectivitybits<<
			"bits)"<<endl;
		for (j=this->NumberOfFilters-1;j>=0;j--)
		{
			encodedbits+=ConnectivityBitrate[j]+GeometryBitrate[j];
			
			double RelativeBitrate=encodedbits/(double) this->Filters[0]->GetOutput()->GetNumberOfPoints();

			connectivitybits+=ConnectivityBitrate[j];
			Repport<<"Level "<<this->NumberOfFilters-j
				<<": "<<this->Filters[j]->GetOutput()->GetNumberOfCells()
				<<"f, "<<this->Filters[j]->GetOutput()->GetNumberOfPoints()
				<<"v, valence entropy= "<<this->Filters[j]->GetOutput()->GetValenceEntropy()
				<<", total data: "<<RelativeBitrate
				<<" bits/v (connectivity: "<<connectivitybits
				<<"bits, "<<
				(float)((float)ConnectivityBitrate[j])
				/((float)this->Filters[j]->GetOutput()->GetNumberOfPoints()-(float)this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints())
				<<" bits/vertex for this level)"<<endl;
		};

		Repport<<"Global coding: "<<(double) encodedbits/this->Filters[0]->GetOutput()->GetNumberOfPoints()
			<<" bits/vertex, connectivity : "<<
			(double) connectivitybits/this->Filters[0]->GetOutput()->GetNumberOfPoints()<<
			" bits/vertex, geometry : "<<
			(double) (encodedbits-connectivitybits)/this->Filters[0]->GetOutput()->GetNumberOfPoints()
			<<"bits/vertex"<<endl;
		Repport<<"File size: "<<((int) encodedbits)/8<<"bytes"<<endl;

		Repport.close();
	}
}

void vtkMultiresolutionIO::Analyse()
{
	double start2, finish2;
	int i,j,revert;

	this->AnalysisMeshes[0]=this->Input;
	this->MeshWindow->SetInputData(this->AnalysisMeshes[0]);
	if (this->Display!=0)
	{
		this->MeshWindow->Render();
		cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
		this->MeshWindow->Interact();
	}

#if ( (VTK_MAJOR_VERSION >= 5))
	this->StartTime = this->Timer->GetUniversalTime();
#else
	this->StartTime = this->Timer->GetCurrentTime();
#endif

	i=0;
	int merged=1;
	while (merged==1)
	{
		if (i+1>=this->MaxNumberOfLevels)
		{
			merged=0;
		}
		else
		{
			this->Filters[i]=vtkWaveletSubdivisionFilter::New();
			this->Filters[i]->SetGeometryPrediction(this->GeometryPrediction);
			this->Filters[i]->SetLifting(this->Lifting);
			this->Filters[i]->SetLiftingRadius(this->LiftingRadius);
			this->Filters[i]->SetArithmeticType(this->ArithmeticType);
			this->Filters[i]->Quantization=this->Quantization;
			this->Filters[i]->SetMergeInput(this->AnalysisMeshes[i]);
			this->Filters[i]->SetGeometryCriterion(this->GeometricConstraint);
			this->Filters[i]->GeometryPrediction=this->GeometryPrediction;
			this->Filters[i]->SetCurvatureTreshold(this->EdgeAngleThreshold);
			this->Filters[i]->SetWaveletTreshold(this->WGC);
			this->Filters[i]->SetDisplayEfficiency(this->DisplayEfficiency);

			this->Filters[i]->SolveInverseProblem(0);
			this->AnalysisMeshes[i+1]=this->Filters[i]->GetMergeOutput();

			j=0;
			while (((this->AnalysisMeshes[i+1]->GetNumberOfPoints()==this->AnalysisMeshes[i]->GetNumberOfPoints())
				||(this->Filters[i]->remaining_faces!=0))&&(j<this->AnalysisMeshes[i]->GetNumberOfCells()-1))
			{
//				cout<<"Fail "<<j<<", "<<this->Filters[i]->remaining_faces<<" remaining faces"<<endl;
				this->Filters[i]->Delete();
				this->Filters[i]=vtkWaveletSubdivisionFilter::New();
				this->Filters[i]->SetMergeInput(this->AnalysisMeshes[i]);
				this->Filters[i]->SetGeometryPrediction(this->GeometryPrediction);
				this->Filters[i]->SetLifting(this->Lifting);
				this->Filters[i]->SetLiftingRadius(this->LiftingRadius);
				this->Filters[i]->Quantization=this->Quantization;
				this->Filters[i]->SetArithmeticType(this->ArithmeticType);
				this->Filters[i]->SetGeometryCriterion(this->GeometricConstraint);
				this->Filters[i]->SetCurvatureTreshold(this->EdgeAngleThreshold);
				this->Filters[i]->SetWaveletTreshold(this->WGC);
				this->Filters[i]->SetDisplayEfficiency(this->DisplayEfficiency);

				j++;
//				this->AnalysisMeshes[i+1]->Delete();
				this->Filters[i]->SolveInverseProblem(j);
				this->AnalysisMeshes[i+1]=this->Filters[i]->GetMergeOutput();
			}

			if (this->DisplayText)
			{
				cout<<"i="<<i<<" ,original : "<<this->Filters[i]->GetMergeInput()->GetNumberOfCells()<<" faces, "<<this->Filters[i]->GetMergeInput()->GetPoints()->GetNumberOfPoints()<<" vertices";
				cout<<", simplified : "<<this->Filters[i]->GetMergeOutput()->GetNumberOfCells()<<" faces, "<<this->Filters[i]->GetMergeOutput()->GetPoints()->GetNumberOfPoints()<<" vertices, "<<this->Filters[i]->GetMergeOutput()->GetNumberOfEdges()<<" edges"<<endl;
			}

			if (this->Filters[i]->remaining_faces!=0)
				revert=0;

			if (this->Display!=0)
			{
				if (this->GeometricConstraint==1)
				{	
					vtkLookupTable *lut=vtkLookupTable::New();
					this->MeshWindow->SetLookupTable(lut);
					this->AnalysisMeshes[i+1]->ComputeSharpVertices(this->EdgeAngleThreshold);
					lut->SetHueRange(0.667,0.0);
				}
				this->MeshWindow->SetInputData(this->AnalysisMeshes[i+1]);
				this->MeshWindow->Render();
#if ( (VTK_MAJOR_VERSION >= 5))
				start2=this->Timer->GetUniversalTime();
				finish2=this->Timer->GetUniversalTime();
				while (finish2 - start2<this->Time)
					finish2=this->Timer->GetUniversalTime();

#else
				start2=this->Timer->GetCurrentTime();
				finish2=this->Timer->GetCurrentTime();
				while (finish2 - start2<this->Time)
					finish2=this->Timer->GetCurrentTime();
#endif
				if (this->Display==2)
				{
					cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
					this->MeshWindow->Interact();
				}
			}

			if ((this->AnalysisMeshes[i]->GetNumberOfCells()>0)
				&&(this->Filters[i]->GetMergeOutput()->GetNumberOfPoints()<this->Filters[i]->GetMergeInput()->GetNumberOfPoints()))
				merged=1;
			else
				merged=0;
			if (this->Filters[i]->remaining_faces!=0)
				merged=0;
		}
		i++;

	}
	this->NumberOfFilters=i-1;


	if (this->Display!=0)
		this->MeshWindow->Interact();
}
void vtkMultiresolutionIO::Synthetize()
{
	double start2, finish2;
	vtkIdType Id;
	int j;

	this->PointsIds->SetNumberOfIds(this->AnalysisMeshes[0]->GetNumberOfPoints());
	for (Id=0;Id<this->AnalysisMeshes[0]->GetNumberOfPoints();Id++)
	{
		this->PointsIds->SetId(Id,Id);
	}

	this->SynthesisMeshes[this->NumberOfFilters]=vtkSurface::New();
	this->SynthesisMeshes[this->NumberOfFilters]->CreateFromPolyData(this->AnalysisMeshes[this->NumberOfFilters]);

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		this->Filters[j]->SetPointsIds(this->PointsIds);
		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);			
		this->Filters[j]->Subdivide();

		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();
		if (this->Display!=0)
		{
			this->MeshWindow->SetInputData(this->Filters[j]->GetOutput());
			this->MeshWindow->Render();
#if ( (VTK_MAJOR_VERSION >= 5))
			start2=this->Timer->GetUniversalTime();
			finish2=this->Timer->GetUniversalTime();
			while (finish2 - start2<this->Time)
				finish2=this->Timer->GetUniversalTime();

#else
			start2=this->Timer->GetCurrentTime();
			finish2=this->Timer->GetCurrentTime();
			while (finish2 - start2<this->Time)
				finish2=this->Timer->GetCurrentTime();
#endif
		}
	}
	if (this->Display!=0)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
		this->MeshWindow->Render();
		this->MeshWindow->Interact();
	}
}

void vtkMultiresolutionIO::DisplayHires()
{
	this->MeshWindow->SetInputData(this->SynthesisMeshes[this->NumberOfFilters-1]);
	this->MeshWindow->Render();
	this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
	this->MeshWindow->Render();
	this->MeshWindow->Interact();
}	

vtkMultiresolutionIO* vtkMultiresolutionIO::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMultiresolutionIO");
	if(ret)
	{
		return (vtkMultiresolutionIO*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkMultiresolutionIO);
}

vtkMultiresolutionIO::vtkMultiresolutionIO()
{
	this->EdgeAngleThreshold=0.3;
	this->GeometricConstraint=0;
	this->NumberOfFilters=0;
	this->WGC=0.25;
	this->Quantization=12;
	this->ArithmeticType=1;
	this->Lifting=0;
	this->LiftingRadius=1;
	this->DisplayEfficiency=0;
	this->Display=0;
	this->DisplayText=0;
	this->WriteRepport=1;
	this->NumberOfBitPlanes=12;
	this->NumberOfStartBitPlanes=1;
	this->WriteOutput=1;
	this->GeometryPrediction=1;
	this->MaxNumberOfLevels=99;
	this->Input=0;
	this->Output=0;
	this->FileType=1;

	this->MeshWindow=RenderWindow::New();
	this->PointsIds=vtkIdList::New();
	this->SignificanceCodes=0;
	this->Timer=vtkTimerLog::New();
	this->SetFileName("out.ddd");
}

vtkMultiresolutionIO::~vtkMultiresolutionIO()
{
	this->PointsIds->Delete();
	int i;
	for (i=0;i<this->NumberOfFilters;i++)
		this->Filters[i]->Delete();
	
	if (this->NumberOfFilters!=0)
		this->SynthesisMeshes[this->NumberOfFilters]->Delete();
		
	if (this->Input)
		this->Input->UnRegister(this);

	if (this->Output)
		this->Output->UnRegister(this);

	this->MeshWindow->Delete();
	this->Timer->Delete();
}
