/*==================================================================================================================
  Program:   3D mesh sequence coding using the combination of spatial and temporal wavelet analysis (Creatis 1997 ~)
  Module:    vtkMultiresolutionIOSeqII.cpp
  Language:  C++
  Date:      2006/10
  Auteurs:   Jae-Won Cho
  Base codeur: vtkMultiresolutionIO.cxx (Sebastien Valette), vtkMultiresolutionIOSeq(Min-Su Kim and Jae-Won Cho)
==================================================================================================================*/

// .NAME vtkMultiresolutionIOSeqII
// .SECTION Description

#include <vtkPLYWriter.h>

#include "vtkMultiresolutionIOSeqII.h"

#define minimum(X, Y)  ((X) < (Y) ? (X) : (Y))
#define maximum(X, Y)  ((X) > (Y) ? (X) : (Y))

vtkMultiresolutionIOSeqII* vtkMultiresolutionIOSeqII::New()
{
  // First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMultiresolutionIOSeqII");
	if(ret)
	{
		return (vtkMultiresolutionIOSeqII*)ret;
	}
  // If the factory was unable to create the object, then create it here.
  return (new vtkMultiresolutionIOSeqII);
}

vtkMultiresolutionIOSeqII::vtkMultiresolutionIOSeqII():vtkMultiresolutionIO()
{
	int i;

	// initialization of member variables
	globalConnectivityBits=0;
	globalGeometryBits=0;
	globalTotalBits=0;
	NumOfPoints=0;
	FileSize=0;

	for(i = 0; i < 100; i++)
	{
		localConnectivityBits[i] = 0;
		localGeometryBits[i] = 0;
		localTotalBits[i] = 0;
	}
}

vtkMultiresolutionIOSeqII::~vtkMultiresolutionIOSeqII()
{
}

void vtkMultiresolutionIOSeqII::CopyForEncoding(int frameindex, int mode)
{
	int j;

// strcpy(inputfilename, prefixOfFilename);
//	strcat(inputfilename, _itoa(frameindex, buf, 10));
//	strcat(inputfilename, ".ply");

	std::stringstream strfile;
	strfile<<prefixOfFilename<<frameindex<<".ply";
	strcpy(inputfilename,strfile.str().c_str());


	vtkSurface *Mesh = vtkSurface::New();
	Mesh->CreateFromFile(inputfilename);

	if (mode == 0)		// for the first frame
	{
		Mesh->QuantizeCoordinates(this->Quantization);
		Mesh->GetScalingFactors(Factor, Tx, Ty, Tz);

		this->SetInput(Mesh);

		this->Analyse();		//vtkSurface(Mesh)->(MergeInput(0)->MergeOutput(0+1))->(MergeInput(0+1)((==MergeOutput(0+1))->MergeOutput)->(MergeInput->MergeOutput)
		this->Synthetize();
		this->Approximate();		
	}
	else if (mode == 1)	// for the other frames
	{
		Mesh->QuantizeCoordinates(Factor, Tx, Ty, Tz);  // we will use same factors with 1st Frame
		this->UpdateGeometry(Mesh);                      // we only update geometry of the other frame
		this->Approximate();
	}
	else;

	BaseMesh=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput();		// base mesh read
	
	// copy of base mesh for temporal wavelet transform in encoder side
	BaseMeshEncodingBuffer[frameindex] = vtkPoints::New();
   BaseMeshEncodingBuffer[frameindex]->DeepCopy(BaseMesh->GetPoints());

	// copy of spatial wavelet coefficients for temporal wavelet transform in encoder side
	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		WaveletsEncodingBuffer[frameindex][j] = vtkIntArray::New();		
		WaveletsEncodingBuffer[frameindex][j]->DeepCopy(this->Filters[j]->IntegerWavelets);
	}
}

void vtkMultiresolutionIOSeqII::CopyForDecoding(int frameindex, int mode)
{
   vtkIdType j;

	if(mode == 0)
   {
		this->Filters[0]=vtkWaveletSubdivisionFilter::New();

		this->ArithmeticCoder->StartDecoding();
		BaseMesh=this->DecodeMeshConnectivity();
		this->DecodeMeshGeometry(BaseMesh);
		this->DecodeScalingFactors(BaseMesh);

		int LiftingTest=this->ArithmeticCoder->DecodeByte();
		LiftingTest-=2;
		if (LiftingTest>=0)
		{
			this->Lifting=1;
			this->Filters[0]->SetLifting(1);
			this->LiftingRadius=LiftingTest;
			this->Filters[0]->SetLiftingRadius(LiftingTest);
		}
		else
		{
			if (LiftingTest==-1)
			{
				this->Lifting=0;
				this->Filters[0]->SetLifting(0);
			}
			else
			{
				this->Lifting=2;
				this->Filters[0]->SetLifting(2);
			}
		}
		this->NumberOfFilters=this->ArithmeticCoder->DecodeByte();
		this->ArithmeticCoder->StopDecoding();

		// copy of base mesh for temporal wavelet transform in decoder side
		BaseMeshDecodingBuffer[frameindex] = vtkPoints::New();
		BaseMeshDecodingBuffer[frameindex]->DeepCopy(BaseMesh->GetPoints());

		this->Filters[0]->SetInput(BaseMesh);
		this->Filters[this->NumberOfFilters-1]=this->Filters[0];	

		for (j=this->NumberOfFilters-1;j>=0;j--)
		{
			//1. Set BaseMesh
			if (j!=this->NumberOfFilters-1)
			{
				this->Filters[j]=vtkWaveletSubdivisionFilter::New();
				this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);
			}
			this->Filters[j]->SetIOType(2);	
			this->Filters[j]->SetLifting(this->Lifting);
			this->Filters[j]->SetLiftingRadius(this->LiftingRadius);
			this->Filters[j]->SetArithmeticType(1);
			this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
			this->Filters[j]->Quantization=this->Quantization;
			this->Filters[j]->SetSubdivisionType(1);

			//2. Decode connectivity of each resolution meshes
			this->ArithmeticCoder->StartDecoding();		
			this->Filters[j]->Subdivide();
			this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();
			this->ArithmeticCoder->StopDecoding();
			///////////////////////////////////////

			//3. Decode new points of each resolution meshes
			this->ArithmeticCoder->StartDecoding();
			this->Filters[j]->ReadCoefficients();
			this->ArithmeticCoder->StopDecoding();

			// copy of spatial wavelet coefficients for temporal wavelet transform in decoder side
			WaveletsDecodingBuffer[frameindex][j] = vtkIntArray::New();
			WaveletsDecodingBuffer[frameindex][j]->DeepCopy(this->Filters[j]->IntegerWavelets);

			//4. Reconstruction using connectivity from 1 and geometry from 2.
			this->Filters[j]->Reconstruct();
		}

		this->Output=this->SynthesisMeshes[0];

		BaseMesh->GetScalingFactors(Factor,Tx,Ty,Tz);
		this->Output->SetScalingFactors(Factor,Tx,Ty,Tz);
   }
   
	else
	{
		this->ArithmeticCoder->StartDecoding();	
		this->DecodeMeshGeometry(this->BaseMesh);
		this->ArithmeticCoder->StopDecoding();

		// copy of base mesh for temporal wavelet transform in decoder side
		BaseMeshDecodingBuffer[frameindex] = vtkPoints::New();
		BaseMeshDecodingBuffer[frameindex]->DeepCopy(BaseMesh->GetPoints());
		/////////

		this->Filters[this->NumberOfFilters-1]->SetInput(BaseMesh);

		for (j=this->NumberOfFilters-1;j>=0;j--)
		{
			this->ArithmeticCoder->StartDecoding();
			this->Filters[j]->ReadCoefficients();
			this->ArithmeticCoder->StopDecoding();

			// copy of spatial wavelet coefficients for temporal wavelet transform in decoder side
			WaveletsDecodingBuffer[frameindex][j] = vtkIntArray::New();
			WaveletsDecodingBuffer[frameindex][j]->DeepCopy(this->Filters[j]->IntegerWavelets);

   		this->Filters[j]->Reconstruct();
		}	

		this->Output=this->SynthesisMeshes[0];

		BaseMesh->GetScalingFactors(Factor,Tx,Ty,Tz);
		this->Output->SetScalingFactors(Factor,Tx,Ty,Tz);
	}
}

void vtkMultiresolutionIOSeqII::WriteIFrame(int frameindex)
{
   vtkIdType Id;
	int i, j;
	int temp;
	double pt[3];
   
	double numberoffaces=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfCells();	// for calculating the bitrate
	double numberofvertices=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints();	// for calculating the bitrate
   double cost=32.0+ceil(log(numberofvertices)/log(2.0))*3.0*numberoffaces;	// for calculating the bitrate
	
   globalConnectivityBits = (int) cost;	// initialize encoding cost
   localConnectivityBits[this->NumberOfFilters] = (int) cost;
	NumOfPoints = this->Filters[0]->GetOutput()->GetNumberOfPoints();	// for calculating the bitrate

	// number of vertex of the highest resolution level(original mesh)
	for (Id=0;Id<this->AnalysisMeshes[0]->GetNumberOfPoints();Id++)
	{
		this->PointsIds->SetId(Id,Id);		
	}

	BaseMesh=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput();		// base mesh read
	
	// exchange basemesh by BaseMeshEncodingBuffer
	for(i=0; i<BaseMesh->GetNumberOfPoints(); i++)
	{
		BaseMeshEncodingBuffer[frameindex]->GetPoint(i,pt);	// the first frame BaseMesh
		BaseMesh->SetPointCoordinates(i,pt);
	}
	
	//////////Begin of Encode BaseMesh///////////////////////////
	// connectivity coding of BaseMesh
 	this->ArithmeticCoder->StartCoding();
	this->EncodeMeshConnectivity(BaseMesh);
	this->EncodeMeshGeometry(BaseMesh);
	this->EncodeScalingFactors(this->Filters[0]->GetMergeInput()); //this->Filters[0]->GetMergeInput() Input of vtksurface->QuantizeCoordinates()

	if (this->Lifting==0)
		this->ArithmeticCoder->EncodeByte(1);
	else
	{
		if (this->Lifting==2)
			this->ArithmeticCoder->EncodeByte(0);
		else //
			this->ArithmeticCoder->EncodeByte(this->LiftingRadius+2);
	}

	this->ArithmeticCoder->EncodeByte(this->NumberOfFilters);

   temp = this->ArithmeticCoder->StopCoding()*8;
   localTotalBits[this->NumberOfFilters] = temp;  // initialize local totalbits.
	globalTotalBits = temp;	// for calculating the bitrate
	/////////End of Encode BaseMesh/////////////////////////////////

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		/////////////Begin of Encode connectivity 
		this->ArithmeticCoder->StartCoding();
		this->Filters[j]->SetIOType(1);		// 1: write, 0: read

		this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[j]->SetPointsIds(this->PointsIds);
		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);		// BaseMesh
		this->Filters[j]->Subdivide();

		// exchange spatial wavelet coefficients by WaveletsEncodingBuffer
		this->Filters[j]->IntegerWavelets->DeepCopy(WaveletsEncodingBuffer[frameindex][j]);
		/////////////////////

		//temp = 0;	// for calculating the bitrate
		temp = this->ArithmeticCoder->StopCoding()*8;	// for calculating the bitrate
		globalConnectivityBits = globalConnectivityBits + temp;	// for calculating the bitrate
		globalTotalBits = globalTotalBits +temp;	// for calculating the bitrate
      localConnectivityBits[j] = localConnectivityBits[j] +temp;
      localTotalBits[j] = localTotalBits[j] + temp;
		/////////////End of Encode connectivity

		/////////////Begin of Encode geometry of wavelet coefficients
		this->ArithmeticCoder->StartCoding();
		this->Filters[j]->WriteCoefficients();
		this->Filters[j]->Reconstruct();		

      temp = this->ArithmeticCoder->StopCoding()*8;
      globalTotalBits = globalTotalBits + temp;	// for calculating the bitrate
      localTotalBits[j] = localTotalBits[j] + temp;
		/////////////End of Encode geometry of wavelet coefficients
	}
}

void vtkMultiresolutionIOSeqII::WritePFrame(int frameindex)
{
	vtkIdType i, j;
   int temp;
	double pt[3];
	
	BaseMesh=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput();		// base mesh read

	// exchange basemesh by BaseMeshEncodingBuffer
	for(i=0; i<BaseMesh->GetNumberOfPoints(); i++)
	{
		BaseMeshEncodingBuffer[frameindex]->GetPoint(i,pt);	// 0th frame BaseMesh
		BaseMesh->SetPointCoordinates(i,pt);
	}

	//////////Begin of Encode BaseMesh///////////////////////////
	this->ArithmeticCoder->StartCoding();
	this->EncodeMeshGeometry(BaseMesh);			/// MOTION VECTOR OF GEOMETRY 
   temp = this->ArithmeticCoder->StopCoding()*8;
   localTotalBits[this->NumberOfFilters] = localTotalBits[this->NumberOfFilters] + temp;
   globalTotalBits = globalTotalBits + temp;
	/////////End of Encode BaseMesh/////////////////////////////////

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		// exchange spatial wavelet coefficients by WaveletsEncodingBuffer
		this->Filters[j]->IntegerWavelets->DeepCopy(WaveletsEncodingBuffer[frameindex][j]);
		/////////////////////

		////////////Begin of Encode geometry of wavelet coefficients
		this->ArithmeticCoder->StartCoding();
      this->Filters[j]->ArithmeticCoder = this->ArithmeticCoder;
		this->Filters[j]->WriteCoefficients();		// give motion vector here
		this->Filters[j]->Reconstruct();		// why?	
      temp = this->ArithmeticCoder->StopCoding()*8;
		localTotalBits[j] = localTotalBits[j] + temp;
      globalTotalBits = globalTotalBits + temp;

      /////////////End of Encode geometry of wavelet coefficients
	}
}

void vtkMultiresolutionIOSeqII::ReadIFrame(int frameindex)
{
	int i, j;
	double pt[3];
	
	this->Filters[0]=vtkWaveletSubdivisionFilter::New();

	this->ArithmeticCoder->StartDecoding();
	BaseMesh=this->DecodeMeshConnectivity();
	this->DecodeMeshGeometry(BaseMesh);
	this->DecodeScalingFactors(BaseMesh);

	int LiftingTest=this->ArithmeticCoder->DecodeByte();
	LiftingTest-=2;
	if (LiftingTest>=0)
	{
		this->Lifting=1;
		this->Filters[0]->SetLifting(1);
		this->LiftingRadius=LiftingTest;
		this->Filters[0]->SetLiftingRadius(LiftingTest);
	}
	else
	{
		if (LiftingTest==-1)
		{
			this->Lifting=0;
			this->Filters[0]->SetLifting(0);
		}
		else
		{
			this->Lifting=2;
			this->Filters[0]->SetLifting(2);
		}
	}
	this->NumberOfFilters=this->ArithmeticCoder->DecodeByte();
	this->ArithmeticCoder->StopDecoding();

	// exchange basemesh by BaseMeshDecodingBuffer
	for(i=0; i<BaseMesh->GetNumberOfPoints(); i++)
	{
		BaseMeshDecodingBuffer[frameindex]->GetPoint(i,pt);	// 0th frame BaseMesh
		BaseMesh->SetPointCoordinates(i,pt);
	}
	/////	

	this->Filters[0]->SetInput(BaseMesh);
	this->Filters[this->NumberOfFilters-1]=this->Filters[0];	

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		// set BaseMesh
		if (j!=this->NumberOfFilters-1)
		{
			this->Filters[j]=vtkWaveletSubdivisionFilter::New();
			this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);
		}
		this->Filters[j]->SetIOType(2);	
		this->Filters[j]->SetLifting(this->Lifting);
		this->Filters[j]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[j]->SetArithmeticType(1);
		this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[j]->Quantization=this->Quantization;
		this->Filters[j]->SetSubdivisionType(1);

		// decode connectivity of each resolution meshes
		this->ArithmeticCoder->StartDecoding();		
		this->Filters[j]->Subdivide();
		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();
		this->ArithmeticCoder->StopDecoding();
		///////////////////////////////////////

		// decode new points of each resolution meshes
		this->ArithmeticCoder->StartDecoding();
		this->Filters[j]->ReadCoefficients();
		this->ArithmeticCoder->StopDecoding();

		// exchange spatial wavelet coefficients by WaveletsDecodingBuffer
		this->Filters[j]->IntegerWavelets->DeepCopy(WaveletsDecodingBuffer[frameindex][j]);
		/////////////////////

		// reconstruction using connectivity
		this->Filters[j]->Reconstruct();
	}

	this->Output=this->SynthesisMeshes[0];

   BaseMesh->GetScalingFactors(Factor,Tx,Ty,Tz);
   this->Output->SetScalingFactors(Factor,Tx,Ty,Tz);
}

void vtkMultiresolutionIOSeqII::ReadPFrame(int frameindex)
{
	vtkIdType i, j;
	double pt[3];

	this->ArithmeticCoder->StartDecoding();	
	this->DecodeMeshGeometry(this->BaseMesh);
	this->ArithmeticCoder->StopDecoding();

	// exchange basemesh by BaseMeshDecodingBuffer
	for(i=0; i<BaseMesh->GetNumberOfPoints(); i++)
	{
		BaseMeshDecodingBuffer[frameindex]->GetPoint(i,pt);	// 0th frame BaseMesh
		BaseMesh->SetPointCoordinates(i,pt);
	}
	/////

	this->Filters[this->NumberOfFilters-1]->SetInput(BaseMesh);

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		// decode wavelet coefficients
		this->ArithmeticCoder->StartDecoding();
		this->Filters[j]->ReadCoefficients();
		this->ArithmeticCoder->StopDecoding();

		// exchange spatial wavelet coefficients by WaveletsDecodingBuffer
		this->Filters[j]->IntegerWavelets->DeepCopy(WaveletsDecodingBuffer[frameindex][j]);
		/////////////////////

		// reconstruction using connectivity
   	this->Filters[j]->Reconstruct();
	}	

	this->Output=this->SynthesisMeshes[0];

   BaseMesh->GetScalingFactors(Factor,Tx,Ty,Tz);
   this->Output->SetScalingFactors(Factor,Tx,Ty,Tz);
}

void vtkMultiresolutionIOSeqII::UpdateGeometry(vtkSurface *Mesh)
{
    double pt[3];
    vtkIdType Id;

    //Update the geometry. for the highest resolution level only.
    //With this method, we can avoid the connectivity encode and decode.
    for (Id=0;Id<this->SynthesisMeshes[0]->GetNumberOfPoints();Id++)
	{
        Mesh->GetPointCoordinates(this->PointsIds->GetId(Id),pt);
        this->SynthesisMeshes[0]->SetPointCoordinates(Id, pt);
	}
}

void vtkMultiresolutionIOSeqII::EncodeSequence()
{
	vtkIdType i;

	/* Encode sequence Part I : Store Temporal Wavelet Coefficients */
	// Storing spatial wavelet coefficients into BaseMeshEncodingBuffer and WaveletsEncodingBuffer
	for(i = first; i <= last; i++)
	{
		if(i == first)
		{
			CopyForEncoding(first, 0);		// for the first
		}

		else if(i != first)
		{
			CopyForEncoding(i, 1);	// for the other frames
		}
	}

// for print wavelet coefficients
/*	double pt[3];
	for(vtkIdType l = first; l <= last; l++)
	{
		BaseMeshEncodingBuffer[l]->GetPoint(0, pt);
		cout << pt[0] << endl;
	}
*/

/*	double pt[3];
	for(vtkIdType l = first; l <= last; l++)
	{
		WaveletsEncodingBuffer[l][5]->GetTuple(0, pt);
		cout << pt[0] << endl;
	}
*/
	// Temporal wavelet transform
	TemporalWaveletDecomposition();

	// for print wavelet coefficients
/*	double pt[3];
	for(vtkIdType l = first; l <= last; l++)
	{
		BaseMeshEncodingBuffer[l]->GetPoint(0, pt);
		cout << pt[0] << endl;
	}
*/

	double pt[3];
	for(vtkIdType l = first; l <= last; l++)
	{
		WaveletsEncodingBuffer[l][5]->GetTuple(0, pt);
		cout << pt[0] << endl;
	}

	/* Encode sequence Part II : Arithmetic coding by using exchange the spatial wavelet coefficients by temporal wavelet coefficients */
   this->ArithmeticCoder = vtkArithmeticCoder::New();

   char outbitstream[1000];


	std::stringstream strfile;

	if(tempoWaveDecompLevel == 0)
	{
		strfile<<prefixOfFilename<<"_F"<<first<<"_L"<<last<<"_Q"
			<<this->Quantization<<"_TWDL"<<tempoWaveDecompLevel<<".wmc";


//		strcpy(outbitstream, prefixOfFilename);   
//		strcat(outbitstream, "_F");
//		strcat(outbitstream, itoa(first, buf, 10));
//		strcat(outbitstream, "_L");
//		strcat(outbitstream, itoa(last, buf, 10));
//		strcat(outbitstream, "_Q");
//		strcat(outbitstream, itoa(this->Quantization, buf, 10));
//		strcat(outbitstream, "_TWDL");
//		strcat(outbitstream, itoa(tempoWaveDecompLevel, buf, 10)); 
//		strcat(outbitstream, ".wmc");
	}

	else
	{
		strfile<<prefixOfFilename<<"_F"<<first<<"_L"<<last<<"_Q"
			<<this->Quantization<<"_TWDL"<<tempoWaveDecompLevel<<"_TWM"<<tempoWaveMode<<"_TWF"
			<<tempoWaveFilter<<".wmc";

//		strcpy(outbitstream, prefixOfFilename);   
//		strcat(outbitstream, "_F");
//		strcat(outbitstream, itoa(first, buf, 10));
//		strcat(outbitstream, "_L");
//		strcat(outbitstream, itoa(last, buf, 10));
//		strcat(outbitstream, "_Q");
//		strcat(outbitstream, itoa(this->Quantization, buf, 10));
//		strcat(outbitstream, "_TWDL");
//		strcat(outbitstream, itoa(tempoWaveDecompLevel, buf, 10)); 
//		strcat(outbitstream, "_TWM");
//		strcat(outbitstream, itoa(tempoWaveMode, buf, 10)); 
//		strcat(outbitstream, "_TWF");
//		strcat(outbitstream, itoa(tempoWaveFilter, buf, 10)); 
//		strcat(outbitstream, ".wmc");
	}
   
 //  this->ArithmeticCoder->OpenOutput(outbitstream,1);
	strcpy(outbitstream,strfile.str().c_str());
   this->ArithmeticCoder->OpenFile(outbitstream,1);
	this->ArithmeticCoder->InitConnectivityQSModels(1);

   for(i = first ; i <= last ; i++)
	{
      if(i==first)
      {
         Encode(first, 0);	// for the first frame
      }

		else
		{
			Encode(i, 1);	// for the other frames
		} 
	}

//   this->ArithmeticCoder->CloseOutput();
   this->ArithmeticCoder->CloseFile();
	
   GetGlobalBitrate();
   GetLocalBitrate();
}

void vtkMultiresolutionIOSeqII::DecodeSequence()
{
   vtkIdType i;

	char inbitstream[1000];
	std::stringstream strfile;

	if(tempoWaveDecompLevel == 0)
	{
		strfile<<prefixOfFilename<<"_F"<<first<<"_L"<<last<<"_Q"
			<<this->Quantization<<"_TWDL"<<tempoWaveDecompLevel<<".wmc";

//		strcpy(inbitstream, prefixOfFilename);   
//		strcat(inbitstream, "_F");
//		strcat(inbitstream, itoa(first, buf, 10));
//		strcat(inbitstream, "_L");
//		strcat(inbitstream, itoa(last, buf, 10));
//		strcat(inbitstream, "_Q");
//		strcat(inbitstream, itoa(this->Quantization, buf, 10));
//		strcat(inbitstream, "_TWDL");
//		strcat(inbitstream, itoa(tempoWaveDecompLevel, buf, 10));
//		strcat(inbitstream, ".wmc");
	}

	else
	{
		strfile<<prefixOfFilename<<"_F"<<first<<"_L"<<last<<"_Q"
			<<this->Quantization<<"_TWDL"<<tempoWaveDecompLevel<<"_TWM"<<tempoWaveMode<<"_TWF"
			<<tempoWaveFilter<<".wmc";


//		strcpy(inbitstream, prefixOfFilename);   
//		strcat(inbitstream, "_F");
//		strcat(inbitstream, itoa(first, buf, 10));
//		strcat(inbitstream, "_L");
//		strcat(inbitstream, itoa(last, buf, 10));
//		strcat(inbitstream, "_Q");
//		strcat(inbitstream, itoa(this->Quantization, buf, 10));
//		strcat(inbitstream, "_TWDL");
//		strcat(inbitstream, itoa(tempoWaveDecompLevel, buf, 10));
//		strcat(inbitstream, "_TWM");
//		strcat(inbitstream, itoa(tempoWaveMode, buf, 10)); 
//		strcat(inbitstream, "_TWF");
//		strcat(inbitstream, itoa(tempoWaveFilter, buf, 10)); 
//		strcat(inbitstream, ".wmc");
	}


	strcpy(inbitstream,strfile.str().c_str());
 

	/* Decode sequence Part I : Store Temporal Wavelet Coefficients */
	// Storing temporal wavelet coefficients into BaseMeshDecodingBuffer and WaveletsDecodingBuffer
   this->ArithmeticCoder = vtkArithmeticCoder::New();

//   this->ArithmeticCoder->OpenOutput(inbitstream,0);
   this->ArithmeticCoder->OpenFile(inbitstream,0);
   this->ArithmeticCoder->InitConnectivityQSModels(0);

	for(i = first; i <= last; i++)
	{
		if(i == first)
		{
			CopyForDecoding(first, 0);		// I-Frame
		}

		else if(i != first)
		{
			CopyForDecoding(i, 1);	// For P-Frame
		}
	}

//	this->ArithmeticCoder->CloseOutput();
	this->ArithmeticCoder->CloseFile();

	// Temporal inverse wavelet transform
	TemporalWaveletReconstruction();

   /* Decode sequence Part II : Arithmetic decoding by using exchange temporal wavelet coefficients by reconstructed spatial wavelet coefficients */
   this->ArithmeticCoder = vtkArithmeticCoder::New();

 //  this->ArithmeticCoder->OpenOutput(inbitstream,0);
   this->ArithmeticCoder->OpenFile(inbitstream,0);
   this->ArithmeticCoder->InitConnectivityQSModels(0);

   for(i = first ; i <= last; i++)
	{
		if(i==first)
      {
         Decode(first, 0);	// decode I-Frame
      }
		else
		{
			Decode(i, 1);	// decode P-Frame
		}

	}

//   this->ArithmeticCoder->CloseOutput();
   this->ArithmeticCoder->CloseFile();
}

void vtkMultiresolutionIOSeqII::Encode(int frameindex, int mode)
{
//   strcpy(inputfilename, prefixOfFilename);
//	strcat(inputfilename, itoa(frameindex, buf, 10));
//	strcat(inputfilename, ".ply");

	std::stringstream strfile;
	strfile<<prefixOfFilename<<frameindex<<".ply";
	strcpy(inputfilename,strfile.str().c_str());

	vtkSurface *Mesh = vtkSurface::New();
	Mesh->CreateFromFile(inputfilename);

	if (mode == 0)		// for the first frame
	{
		Mesh->QuantizeCoordinates(this->Quantization);
		Mesh->GetScalingFactors(Factor, Tx, Ty, Tz);

		this->SetInput(Mesh);

		this->Analyse();		//vtkSurface(Mesh)->(MergeInput(0)->MergeOutput(0+1))->(MergeInput(0+1)((==MergeOutput(0+1))->MergeOutput)->(MergeInput->MergeOutput)
		this->Synthetize();
		this->Approximate();
		this->WriteIFrame(frameindex);
	}
	else if (mode == 1)	// for the other frames
	{
		Mesh->QuantizeCoordinates(Factor, Tx, Ty, Tz);  // we will use same factors with 1st Frame
		this->UpdateGeometry(Mesh);                      // we only update geometry of 2nd frame.
		this->Approximate();
		this->WritePFrame(frameindex);
	}
	else;
}

void vtkMultiresolutionIOSeqII::Decode(int frameindex, int mode)
{
	std::stringstream strfile;

	

	if(tempoWaveDecompLevel == 0)
	{
		strfile<<"C:/cvsroot/WaveMSeq/datas/"<<prefixOfFilename<<"/PreviousWorks/M3DMC/Q"
			<<this->Quantization<<"/TWDL"<<tempoWaveDecompLevel<<"/j/"<<"Recon_Q"
			<<this->Quantization<<"_TWDL"<<tempoWaveDecompLevel<<"_"<<prefixOfFilename<<frameindex<<".ply";

//		strcpy(outputfilename,"C:/cvsroot/WaveMSeq/datas/");
//		strcat(outputfilename,prefixOfFilename);
//		strcat(outputfilename,"/PreviousWorks/M3DMC/Q");	// M3DMC : Modified 3DMC
//		strcat(outputfilename, itoa(this->Quantization, buf, 10));
//		strcat(outputfilename, "/TWDL");
//		strcat(outputfilename, itoa(tempoWaveDecompLevel, buf, 10));
//		strcat(outputfilename, "/j/");
//		strcat(outputfilename, "Recon_Q");        // set prefix for output filename.
//		strcat(outputfilename, itoa(this->Quantization, buf, 10));
//		strcat(outputfilename, "_TWDL");
//		strcat(outputfilename, itoa(tempoWaveDecompLevel, buf, 10));
//		strcat(outputfilename, "_");
//		strcat(outputfilename, prefixOfFilename);    
//		strcat(outputfilename, itoa(frameindex, buf, 10));
//		strcat(outputfilename, ".ply");
	}

	else
	{
		strfile<<"C:/cvsroot/WaveMSeq/datas/"<<prefixOfFilename<<"/STW/Q"
			<<this->Quantization<<"/TWDL"<<tempoWaveDecompLevel<<"/TWM"<<tempoWaveMode
			<<"/TWF"<<tempoWaveFilter<<"/j/"<<"Recon_Q"<<this->Quantization
			<<"_TWDL"<<tempoWaveDecompLevel<<"_TWM"<<tempoWaveMode<<"_TWF"<<tempoWaveFilter
			<<"_"<<prefixOfFilename<<frameindex<<".ply";


//		strcpy(outputfilename,"C:/cvsroot/WaveMSeq/datas/");
//		strcat(outputfilename,prefixOfFilename);
//		strcat(outputfilename,"/STW/Q");	// STW : Spatial-Time Wavelet
//		strcat(outputfilename, itoa(this->Quantization, buf, 10));
//		strcat(outputfilename, "/TWDL");
//		strcat(outputfilename, itoa(tempoWaveDecompLevel, buf, 10));
//		strcat(outputfilename, "/TWM");
//		strcat(outputfilename, itoa(tempoWaveMode, buf, 10)); 
//		strcat(outputfilename, "/TWF");
//		strcat(outputfilename, itoa(tempoWaveFilter, buf, 10));
//		strcat(outputfilename, "/j/");
//		strcat(outputfilename, "Recon_Q");        // set prefix for output filename.
//		strcat(outputfilename, itoa(this->Quantization, buf, 10));
//		strcat(outputfilename, "_TWDL");
//		strcat(outputfilename, itoa(tempoWaveDecompLevel, buf, 10));
//		strcat(outputfilename, "_TWM");
//		strcat(outputfilename, itoa(tempoWaveMode, buf, 10)); 
//		strcat(outputfilename, "_TWF");
//		strcat(outputfilename, itoa(tempoWaveFilter, buf, 10)); 
//		strcat(outputfilename, "_");
//		strcat(outputfilename, prefixOfFilename);    
//		strcat(outputfilename, itoa(frameindex, buf, 10));
//		strcat(outputfilename, ".ply");
	}


		strcpy(outputfilename,strfile.str().c_str());
 

	if (mode == 0)		// for the first frame
	{
		ReadIFrame(frameindex);
		this->GetOutput()->GetScalingFactors(Factor, Tx, Ty, Tz);            
	}
	else if (mode == 1)	// for the other frames
	{
		ReadPFrame(frameindex);
		this->GetOutput()->SetScalingFactors(Factor, Tx, Ty, Tz);            
	}
	else;

	this->GetOutput()->UnQuantizeCoordinates();

	vtkPLYWriter *Writer=vtkPLYWriter::New();
	Writer->SetFileName(outputfilename);
	Writer->SetFileTypeToASCII();
	Writer->SetInputData(this->GetOutput());
	Writer->Write();
	Writer->Delete();
   
   for (vtkIdType j=1; j<5; j++)
   {
      char outputfilename[1000];
		std::stringstream strfile2;


		if(tempoWaveDecompLevel == 0)
		{
			strfile2<<"C:/cvsroot/WaveMSeq/datas/"<<prefixOfFilename<<"/PreviousWorks/M3DMC/Q"
				<<this->Quantization<<"/TWDL"<<tempoWaveDecompLevel<<"/j/"<<"Recon_Q"
				<<this->Quantization<<"_TWDL"<<tempoWaveDecompLevel<<"_"<<prefixOfFilename<<frameindex<<".ply";


//			strcpy(outputfilename,"C:/cvsroot/WaveMSeq/datas/");
//			strcat(outputfilename,prefixOfFilename);
//			strcat(outputfilename,"/PreviousWorks/M3DMC/Q");
//			strcat(outputfilename, itoa(this->Quantization, buf, 10));
//			strcat(outputfilename, "/TWDL");
//			strcat(outputfilename, itoa(tempoWaveDecompLevel, buf, 10));
//			strcat(outputfilename, "/j");
//			strcat(outputfilename, itoa(this->NumberOfFilters - j, buf, 10));
//			strcat(outputfilename, "/");
//			strcat(outputfilename, "Recon_Q");        // set prefix for output filename.
//			strcat(outputfilename, itoa(this->Quantization, buf, 10));
//			strcat(outputfilename, "_TWDL");
//			strcat(outputfilename, itoa(tempoWaveDecompLevel, buf, 10));
//			strcat(outputfilename, "_");
//			strcat(outputfilename, prefixOfFilename);
//			strcat(outputfilename, itoa(frameindex, buf, 10));
//			strcat(outputfilename, "_j");
//			strcat(outputfilename, itoa(this->NumberOfFilters - j, buf, 10));
//			strcat(outputfilename, ".ply");
		}

		else
		{
			strfile2<<"C:/cvsroot/WaveMSeq/datas/"<<prefixOfFilename<<"/STW/Q"
				<<this->Quantization<<"/TWDL"<<tempoWaveDecompLevel<<"/TWM"<<tempoWaveMode
				<<"/TWF"<<tempoWaveFilter<<"/j/"<<"Recon_Q"<<this->Quantization
				<<"_TWDL"<<tempoWaveDecompLevel<<"_TWM"<<tempoWaveMode<<"_TWF"<<tempoWaveFilter
				<<"_"<<prefixOfFilename<<frameindex<<".ply";


//			strcpy(outputfilename,"C:/cvsroot/WaveMSeq/datas/");
//			strcat(outputfilename,prefixOfFilename);
//			strcat(outputfilename,"/STW/Q");
//			strcat(outputfilename, itoa(this->Quantization, buf, 10));
//			strcat(outputfilename, "/TWDL");
//			strcat(outputfilename, itoa(tempoWaveDecompLevel, buf, 10));
//			strcat(outputfilename, "/TWM");
//			strcat(outputfilename, itoa(tempoWaveMode, buf, 10)); 
//			strcat(outputfilename, "/TWF");
//			strcat(outputfilename, itoa(tempoWaveFilter, buf, 10));
//			strcat(outputfilename, "/j");
//			strcat(outputfilename, itoa(this->NumberOfFilters - j, buf, 10));
//			strcat(outputfilename, "/");
//			strcat(outputfilename, "Recon_Q");        // set prefix for output filename.
//			strcat(outputfilename, itoa(this->Quantization, buf, 10));
//			strcat(outputfilename, "_TWDL");
//			strcat(outputfilename, itoa(tempoWaveDecompLevel, buf, 10));
//			strcat(outputfilename, "_TWM");
//			strcat(outputfilename, itoa(tempoWaveMode, buf, 10)); 
//			strcat(outputfilename, "_TWF");
//			strcat(outputfilename, itoa(tempoWaveFilter, buf, 10));
//			strcat(outputfilename, "_");
//			strcat(outputfilename, prefixOfFilename);
//			strcat(outputfilename, itoa(frameindex, buf, 10));
//			strcat(outputfilename, "_j");
//			strcat(outputfilename, itoa(this->NumberOfFilters - j, buf, 10));
//			strcat(outputfilename, ".ply");
		}

		strcpy(outputfilename,strfile2.str().c_str());
 
      this->Filters[j]->GetOutput()->SetScalingFactors(Factor, Tx, Ty, Tz);
      this->Filters[j]->GetOutput()->UnQuantizeCoordinates();

      vtkPLYWriter *Writer=vtkPLYWriter::New();
	   Writer->SetFileName(outputfilename);
      Writer->SetFileTypeToASCII();
	   Writer->SetInputData(this->Filters[j]->GetOutput());
	   Writer->Write();
	   Writer->Delete();
   }
}

void vtkMultiresolutionIOSeqII::TemporalWaveletDecomposition()
{
	vtkIdType h, i, j, k;
	double pt[3];

	// temporal wavelet transform for base mesh
	for(i = 0; i < BaseMeshEncodingBuffer[first]->GetNumberOfPoints(); i++)
	{
		for(j = 0; j < 3; j++)
		{
			double *lineWaveletBuffer = new double [last-first+1];

			for(k = first; k <= last; k++)
			{
				BaseMeshEncodingBuffer[k]->GetPoint(i, pt);
				lineWaveletBuffer[k-first] = pt[j];
			}

			if(tempoWaveMode == 0)		// tempoWaveMode : 0 = dyadic, 1 = packet
			{ 
				oneD_DWT_Dyadic(lineWaveletBuffer, last-first+1, tempoWaveDecompLevel);
			}

			else
			{
				oneD_DWT_Packet(lineWaveletBuffer, last-first+1, tempoWaveDecompLevel);
			}

			for(k = first; k <= last; k++)
			{
				BaseMeshEncodingBuffer[k]->GetPoint(i, pt);
				pt[j] = lineWaveletBuffer[k-first];
				BaseMeshEncodingBuffer[k]->SetPoint(i, pt);
			}

			delete [] lineWaveletBuffer;
		}
	}

	// temporal wavelet transform for spatial wavelet coefficients
	for(h = this->NumberOfFilters-1; h >= 0; h--)
	{
		for(i = 0; i < this->Filters[h]->IntegerWavelets->GetNumberOfTuples(); i++)
		{
			for(j = 0; j < 3; j++)
			{
				double *lineWaveletBuffer = new double [last-first+1];

				for(k = first; k <= last; k++)
				{
					WaveletsEncodingBuffer[k][h]->GetTuple(i, pt);
					lineWaveletBuffer[k-first] = pt[j];
				}

				if(tempoWaveMode == 0)		// tempoWaveMode : 0 = dyadic, 1 = packet
				{
					oneD_DWT_Dyadic(lineWaveletBuffer, last-first+1, tempoWaveDecompLevel);
				}

				else
				{
					oneD_DWT_Packet(lineWaveletBuffer, last-first+1, tempoWaveDecompLevel);
				}

				for(k = first; k <= last; k++)
				{
					WaveletsEncodingBuffer[k][h]->GetTuple(i, pt);
					pt[j] = lineWaveletBuffer[k-first];
					WaveletsEncodingBuffer[k][h]->SetTuple(i, pt);	
				}

				delete [] lineWaveletBuffer;
			}
		}
	}
}

void vtkMultiresolutionIOSeqII::TemporalWaveletReconstruction()
{
	vtkIdType h, i, j, k;
	double pt[3];

	// temporal wavelet transform for base mesh
	for(i = 0; i < BaseMeshDecodingBuffer[first]->GetNumberOfPoints(); i++)
	{
		for(j = 0; j < 3; j++)
		{
			double *lineWaveletBuffer = new double [last-first+1];

			for(k = first; k <= last; k++)
			{
				BaseMeshDecodingBuffer[k]->GetPoint(i, pt);
				lineWaveletBuffer[k-first] = pt[j];
			}

			if(tempoWaveMode == 0)		// tempoWaveMode : 0 = dyadic, 1 = packet
			{
				oneD_IDWT_Dyadic(lineWaveletBuffer, last-first+1, tempoWaveDecompLevel);
			}

			else
			{
				oneD_IDWT_Packet(lineWaveletBuffer, last-first+1, tempoWaveDecompLevel);
			}

			for(k = first; k <= last; k++)
			{
				BaseMeshDecodingBuffer[k]->GetPoint(i, pt);
				pt[j] = lineWaveletBuffer[k-first];
				BaseMeshDecodingBuffer[k]->SetPoint(i, pt);
			}

			delete [] lineWaveletBuffer;
		}
	}

	// temporal wavelet transform for spatial wavelet coefficients
	for(h = this->NumberOfFilters-1; h >= 0; h--)
	{
		for(i = 0; i < this->Filters[h]->IntegerWavelets->GetNumberOfTuples(); i++)
		{
			for(j = 0; j < 3; j++)
			{
				double *lineWaveletBuffer = new double [last-first+1];

				for(k = first; k <= last; k++)
				{
					WaveletsDecodingBuffer[k][h]->GetTuple(i, pt);
					lineWaveletBuffer[k-first] = pt[j];
				}

				if(tempoWaveMode == 0)		// tempoWaveMode : 0 = dyadic, 1 = packet
				{
					oneD_IDWT_Dyadic(lineWaveletBuffer, last-first+1, tempoWaveDecompLevel);
				}

				else
				{
					oneD_IDWT_Packet(lineWaveletBuffer, last-first+1, tempoWaveDecompLevel);
				}

				for(k = first; k <= last; k++)
				{
					WaveletsDecodingBuffer[k][h]->GetTuple(i, pt);
					pt[j] = lineWaveletBuffer[k-first];
					WaveletsDecodingBuffer[k][h]->SetTuple(i, pt);	
				}

				delete [] lineWaveletBuffer;
			}
		}
	}
}

/* Dyadic */
/* 1-D Discrete Wavelet Transform */
void vtkMultiresolutionIOSeqII::oneD_DWT_Dyadic(double *InputData, vtkIdType data_leng, vtkIdType level)
{
	vtkIdType i = 0, j = 0, index;
	vtkIdType sub_data_leng;

	double *data;

	for(j = 0; j < level; j++)
	{
		sub_data_leng = data_leng / (vtkIdType)pow(2, j);
		data = new double [sub_data_leng];
		
		for(i = 0; i < sub_data_leng; i++)
		{
			data[i] = InputData[i];
		}

		if(tempoWaveFilter == 0)		// tempoWaveFilter : 0 = harr, 1 = 5/3 tap lifting, 2 = 9/7 tap lifting
		{
			fHarr(data, sub_data_leng);
		}
		
		else if(tempoWaveFilter == 1)
		{
			forward_lifting_53(data, sub_data_leng);
		}

		else
		{
			forward_lifting_97(data, sub_data_leng);
		}
				
		index = 0;

		for(i = 0; i < sub_data_leng; i = i + 2)
		{
			InputData[index++] = data[i];
		}
		
		for(i = 1; i < sub_data_leng; i = i + 2)
		{
			InputData[index++] = data[i];
		}

		delete [] data;
	}
}

/* 1-D Inverse Discrete Wavelet Transform */
void vtkMultiresolutionIOSeqII::oneD_IDWT_Dyadic(double *InputData, vtkIdType data_leng, vtkIdType level)
{
	vtkIdType i = 0, j = 0, index;
	vtkIdType sub_data_leng;

	double *data;

	for(j = level - 1; j >= 0; j--)
	{
		sub_data_leng = data_leng / (vtkIdType)pow(2, j);
		data = new double [sub_data_leng];

		index = 0;

		for(i = 0; i < sub_data_leng; i = i + 2)
		{
			data[i] = InputData[index++];
		}

		for(i = 1; i < sub_data_leng; i = i + 2)
		{
			data[i] = InputData[index++];
		}

		if(tempoWaveFilter == 0)		// tempoWaveFilter : 0 = harr, 1 = 5/3 tap lifting, 2 = 9/7 tap lifting
		{
			iHarr(data, sub_data_leng);
		}
		
		else if(tempoWaveFilter == 1)
		{
			inverse_lifting_53(data, sub_data_leng);
		}

		else
		{
			inverse_lifting_97(data, sub_data_leng);
		}

		for(i = 0; i < sub_data_leng; i++)
		{
			InputData[i] = data[i];
		}

		delete [] data;
	}
}

/* packet */
/* 1-D Discrete Wavelet Transform */
void vtkMultiresolutionIOSeqII::oneD_DWT_Packet(double *InputData, vtkIdType data_leng, vtkIdType level)
{
	vtkIdType i = 0, j = 0, k = 0, index1, index2;
	vtkIdType sub_data_leng;

	double *data;

	for(j = 0; j < level; j++)
	{
		index1 = 0;
		index2 = 0;

		sub_data_leng = data_leng / (vtkIdType)pow(2, j);
		data = new double [sub_data_leng];
		
		for(k = 1; k <= (vtkIdType)pow(2, j); k++)
		{	
			for(i = index2; i < index2 + sub_data_leng; i++)
			{
				data[i - index2] = InputData[i];
			}

			if(tempoWaveFilter == 0)		// tempoWaveFilter : 0 = harr, 1 = 5/3 tap lifting, 2 = 9/7 tap lifting
			{
				fHarr(data, sub_data_leng);
			}
			
			else if(tempoWaveFilter == 1)
			{
				forward_lifting_53(data, sub_data_leng);
			}
			
			else
			{
				forward_lifting_97(data, sub_data_leng);
			}

			for(i = index2; i < index2 + sub_data_leng; i = i + 2)
			{
				InputData[index1++] = data[i - index2];
			}
			
			for(i = index2 + 1; i < index2 + sub_data_leng; i = i + 2)
			{
				InputData[index1++] = data[i - index2];
			}

			index2 += sub_data_leng;
		}

		delete [] data;
	}
}

/* 1-D Inverse Discrete Wavelet Transform */
void vtkMultiresolutionIOSeqII::oneD_IDWT_Packet(double *InputData, vtkIdType data_leng, vtkIdType level)
{
	vtkIdType i = 0, j = 0, k = 0, index1, index2;
	vtkIdType sub_data_leng;

	double *data;

	for(j = level - 1; j >= 0; j--)
	{
		index1 = 0;
		index2 = 0;

		sub_data_leng = data_leng / (long)pow(2, j);
		data = new double [sub_data_leng];

		for(k = 1; k <= (vtkIdType)pow(2, j); k++)
		{
			for(i = index2; i < index2 + sub_data_leng; i = i + 2)
			{
				data[i - index2] = InputData[index1++];
			}

			for(i = index2 + 1; i < index2 + sub_data_leng; i = i + 2)
			{
				data[i - index2] = InputData[index1++];
			}

			if(tempoWaveFilter == 0)		// tempoWaveFilter : 0 = harr, 1 = 5/3 tap lifting, 2 = 9/7 tap lifting
			{
				iHarr(data, sub_data_leng);
			}
			
			else if(tempoWaveFilter == 1)
			{
				inverse_lifting_53(data, sub_data_leng);
			}

			else
			{
				inverse_lifting_97(data, sub_data_leng);
			}
			
			for(i = index2; i < index2 + sub_data_leng; i++)
			{
				InputData[i] = data[i - index2];
			}

			index2 += sub_data_leng;
		}

		delete [] data;
	}
}

/* forward lifting based on 5/3 tap */
void vtkMultiresolutionIOSeqII::forward_lifting_53(double *InputData, vtkIdType sub_data_leng)
{
	vtkIdType i;
	double *OutputData = new double [sub_data_leng];

	// step1 //
	for(i = 0; i < sub_data_leng; i++)
	{
		if(i == 0)
		{
			OutputData[i] = -1. / 4. * InputData[i + 2] + 1. / 2. * InputData[i + 1] + 6. / 8. * InputData[i];
		}

		else if(i == sub_data_leng - 1)
		{
			OutputData[i] = InputData[i] - InputData[i - 1];
		}

		else if(i != sub_data_leng - 1 && i % 2 == 1)
		{
			OutputData[i] = InputData[i] - ( (InputData[i - 1] + InputData[i + 1]) / 2. );
		}
	}

	// step2 //
	for(i = 0; i < sub_data_leng; i++)
	{
		if(i != 0 && i % 2 == 0)
		{
			OutputData[i] = InputData[i] + ( (OutputData[i - 1] + OutputData[i + 1]) / 4. );
		}
	}

	for(i = 0; i < sub_data_leng; i++)
	{
		InputData[i] = OutputData[i];
	}

	delete [] OutputData;
}

/* inverse lifting based on 5/3 tap */
void vtkMultiresolutionIOSeqII::inverse_lifting_53(double *InputData, vtkIdType sub_data_leng)
{
	vtkIdType i;
	double *OutputData = new double [sub_data_leng];

	// step1 //
	for(i = 0; i < sub_data_leng; i++)
	{
		if(i == 0)
		{
			OutputData[i] = -1. / 2. * InputData[i + 1] + InputData[i];
		}

		if(i == sub_data_leng - 1)
		{
			OutputData[i] = -1. / 4. * InputData[i - 2] + InputData[i - 1] + 6. / 8. * InputData[i];
		}

		else if(i != 0 && i % 2 == 0)
		{
			OutputData[i] = InputData[i] - ( (InputData[i - 1] + InputData[i + 1]) / 4. );
		}
	}

	// step2 //
	for(i = 0; i < sub_data_leng; i++)
	{
		if(i != sub_data_leng - 1 && i % 2 == 1)
		{
			OutputData[i] = InputData[i] + ( (OutputData[i - 1] + OutputData[i + 1]) / 2. );
		}
	}

	for(i = 0; i < sub_data_leng; i++)
	{
		InputData[i] = OutputData[i];
	}

	delete [] OutputData;
}

/* forward lifting based on 9/7 tap */
void vtkMultiresolutionIOSeqII::forward_lifting_97(double *InputData, vtkIdType sub_data_leng)
{
	vtkIdType i, position = 0;
	double *OutputData = new double [sub_data_leng];

	double alpha = -1.586134342059924;
	double beta = -0.052980118572961;
	double gamma = 0.882911075530934;
	double delta = 0.443506852043971;
	double K = 1.230174104914001;

	// step1 //
	for(i = 0; i < sub_data_leng; i++)
	{
		// boundary processing
		if(i == sub_data_leng - 1)
		{
			position = minimum( abs((i + 1) % (2 * (sub_data_leng - 1))) , (2 * (sub_data_leng - 1)) - abs(((i + 1) % (2 * (sub_data_leng - 1)))) );
			OutputData[i] = InputData[i] + alpha * (InputData[i - 1] + InputData[position]);
		}

		// general processing
		else if(i != sub_data_leng - 1 && i % 2 == 1)
		{
			OutputData[i] = InputData[i] + alpha * (InputData[i - 1] + InputData[i + 1]);
		}
	}

	// step2 //
	for(i = 0; i < sub_data_leng; i++)
	{
		// boundary processing
		if(i == 0)
		{
			position = minimum( abs((i - 1) % (2 * (sub_data_leng - 1))) , (2 * (sub_data_leng - 1)) - abs(((i - 1) % (2 * (sub_data_leng - 1)))) );
			OutputData[i] = InputData[i] + beta * (OutputData[position] + OutputData[i + 1]);
		}

		// general processing
		else if(i != 0 && i % 2 == 0)
		{
			OutputData[i] = InputData[i] + beta * (OutputData[i - 1] + OutputData[i + 1]);
		}
	}

	// step3 //
	for(i = 0; i < sub_data_leng; i++)
	{
		// boundary processing
		if(i == sub_data_leng - 1)
		{
			position = minimum( abs((i + 1) % (2 * (sub_data_leng - 1))) , (2 * (sub_data_leng - 1)) - abs(((i + 1) % (2 * (sub_data_leng - 1)))) );
			OutputData[i] = OutputData[i] + gamma * (OutputData[i - 1] + OutputData[position]);
		}

		// general processing
		else if(i != sub_data_leng - 1 && i % 2 == 1)
		{
			OutputData[i] = OutputData[i] + gamma * (OutputData[i - 1] + OutputData[i + 1]);
		}
	}

	// step4 //
	for(i = 0; i < sub_data_leng; i++)
	{
		// boundary processing
		if(i == 0)
		{
			position = minimum( abs((i - 1) % (2 * (sub_data_leng - 1))) , (2 * (sub_data_leng - 1)) - abs(((i - 1) % (2 * (sub_data_leng - 1)))) );
			OutputData[i] = OutputData[i] + delta * (OutputData[position] + OutputData[i + 1]);
		}

		// general processing
		else if(i != 0 && i % 2 == 0)
		{
			OutputData[i] = OutputData[i] + delta * (OutputData[i - 1] + OutputData[i + 1]);
		}
	}

	// step5 //
	for(i = 0; i < sub_data_leng; i++)
	{
		if(i % 2 == 1)
		{
			OutputData[i] = K * OutputData[i];
		}
	}

	// step6 //
	for(i = 0; i < sub_data_leng; i++)
	{
		if(i % 2 == 0)
		{
			OutputData[i] = OutputData[i]/K ;
		}
	}

	for(i = 0; i < sub_data_leng; i++)
	{
		InputData[i] = OutputData[i];
	}

	delete [] OutputData;
}

/* inverse lifting based on 9/7 tap */
void vtkMultiresolutionIOSeqII::inverse_lifting_97(double *InputData, vtkIdType sub_data_leng)
{
	vtkIdType i, position = 0;
	double *OutputData = new double [sub_data_leng];

	double alpha = -1.586134342059924;
	double beta = -0.052980118572961;
	double gamma = 0.882911075530934;
	double delta = 0.443506852043971;
	double K = 1.230174104914001;

	// step1 //
	for(i = 0; i < sub_data_leng; i++)
	{
		if(i % 2 == 1)
		{
			OutputData[i] = InputData[i] / K ;
		}
	}

	// step2 //
	for(i = 0; i < sub_data_leng; i++)
	{
		if(i % 2 == 0)
		{
			OutputData[i] = K * InputData[i];
		}
	}

	// step3 //
	for(i = 0; i < sub_data_leng; i++)
	{
		// boundary processing
		if(i == 0)
		{
			position = minimum( abs((i - 1)%(2 * (sub_data_leng - 1))) , (2 * (sub_data_leng - 1)) - abs(((i - 1) % (2 * (sub_data_leng - 1)))) );
			OutputData[i] = OutputData[i] - delta * (OutputData[position] + OutputData[i + 1]);
		}

		// general processing
		else if(i != 0 && i % 2 == 0)
		{
			OutputData[i] = OutputData[i] - delta * (OutputData[i - 1] + OutputData[i + 1]);
		}
	}

	// step4 //
	for(i = 0; i < sub_data_leng; i++)
	{
		// boundary processing
		if(i == sub_data_leng - 1)
		{
			position = minimum( abs((i + 1) % (2 * (sub_data_leng - 1))) , (2 * (sub_data_leng - 1)) - abs(((i + 1) % (2 * (sub_data_leng - 1)))) );
			OutputData[i] = OutputData[i] - gamma * (OutputData[i - 1] + OutputData[position]);
		}

		// general processing
		else if(i != sub_data_leng - 1 && i % 2 == 1)
		{
			OutputData[i] = OutputData[i] - gamma * (OutputData[i - 1] + OutputData[i + 1]);
		}
	}

	// step5 //
	for(i = 0; i < sub_data_leng; i++)
	{
		// boundary processing
		if(i == 0)
		{
			position = minimum( abs((i - 1) % (2 * (sub_data_leng - 1))) , (2 * (sub_data_leng - 1)) - abs(((i - 1) % (2 * (sub_data_leng - 1)))) );
			OutputData[i] = OutputData[i] - beta * (OutputData[position] + OutputData[i + 1]);
		}

		// general processing
		else if(i != 0 && i % 2 == 0)
		{
			OutputData[i] = OutputData[i] - beta * (OutputData[i - 1] + OutputData[i + 1]);
		}
	}

	// step6 //
	for(i = 0; i < sub_data_leng; i++)
	{
		// boundary processing
		if(i == sub_data_leng - 1)
		{
			position = minimum( abs((i + 1) % (2 * (sub_data_leng - 1))) , (2 * (sub_data_leng - 1)) - abs(((i + 1) % (2 * (sub_data_leng - 1)))) );
			OutputData[i] = OutputData[i] - alpha * (OutputData[i - 1] + OutputData[position]);
		}

		// general processing
		else if(i != sub_data_leng - 1 && i % 2 == 1)
		{
			OutputData[i] = OutputData[i] - alpha * (OutputData[i - 1] + OutputData[i + 1]);
		}
	}

	for(i = 0; i < sub_data_leng; i++)
	{
		InputData[i] = OutputData[i];
	}

	delete [] OutputData;
}

/* Harr 2/2 tap filter */
void vtkMultiresolutionIOSeqII::fHarr(double *InputData, vtkIdType sub_data_leng)
{
	vtkIdType i;
	double *OutputData = new double [sub_data_leng];

	for(i = 0; i < sub_data_leng/2; i++)
	{
		OutputData[2*i] = 1.0/sqrt(2.0) * (InputData[2*i] + InputData[2*i + 1]);
		OutputData[2*i + 1] = 1.0/sqrt(2.0) * (InputData[2*i] - InputData[2*i + 1]);
	}

	for(i = 0; i < sub_data_leng; i++)
	{
		InputData[i] = OutputData[i];
	}

	delete [] OutputData;
}

void vtkMultiresolutionIOSeqII::iHarr(double *InputData, vtkIdType sub_data_leng)
{
	vtkIdType i;
	double *OutputData = new double [sub_data_leng];

	for(i = 0; i < sub_data_leng/2; i++)
	{
		OutputData[2*i] = 1.0/sqrt(2.0) * (InputData[2*i] + InputData[2*i + 1]);
		OutputData[2*i + 1] = 1.0/sqrt(2.0) * (InputData[2*i] - InputData[2*i + 1]);
	}

	for(i = 0; i < sub_data_leng; i++)
	{
		InputData[i] = OutputData[i];
	}

	delete [] OutputData;
}

void vtkMultiresolutionIOSeqII::GetGlobalBitrate()
{
   int NbOfFrames = last + 1 - first;
   double globalConnectivityBitrate;
   double globalGeometryBitrate;
   double globalTotalBitrate;

   cout << "Filename : " << prefixOfFilename << endl;
   cout << "First frame subfix : " << first << endl;
   cout << "Final frame subfix : " << last << endl;
	cout << "Number of Frames : " << NbOfFrames << endl;
   cout << "Temporal Wavelet Decomposition Level : " << tempoWaveDecompLevel << endl;
	cout << "Temporal Wavelet Mode : " << tempoWaveMode << endl;
	cout << "Temporal Wavelet Filter : " << tempoWaveFilter << endl;
	cout << "Quantization level : " << this->Quantization << "bits/coordinate" << endl;

	// Total Bit Rate
	globalTotalBitrate = (double)globalTotalBits/(double)NumOfPoints/(double)NbOfFrames;
	cout << "Total Bit Rate = " << globalTotalBitrate << "(bits/vertex/frame)" << endl;

	// Connectivity Bit Rate
	globalConnectivityBitrate = (double)globalConnectivityBits/(double)NumOfPoints/(double)NbOfFrames;
	cout << "Connectivity Bit Rate = " << globalConnectivityBitrate << "(bits/vertex/frame)" << endl;

	// Geometry Bit Rate
	globalGeometryBits = globalTotalBits - globalConnectivityBits;
	globalGeometryBitrate = (double)globalGeometryBits/(double)NumOfPoints/(double)NbOfFrames;
	cout << "Geometry Bit Rate = " << globalGeometryBitrate << "(bits/vertex/frame)" << endl;

	// File Size
	FileSize = globalTotalBits/8;
	cout << "File Size = " << FileSize << "(bytes)" << endl;
}

//Bits/|v^j|/|frames|
void vtkMultiresolutionIOSeqII::GetLocalBitrate()
{
   int NbOfFrames = last + 1 - first;

	double TotalBitrate;
   double ConnectivityBitrate;
   double GeometryBitrate;
   int GeometryBits;
   int TotalBits = 0;
   int ConnectivityBits = 0;

   vtkIdType i;
 
   vtkIdType NumberOfPoints, NumberOfCells;
   
   NumberOfPoints = this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints();
   NumberOfCells = this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfCells();

   TotalBits = localTotalBits[this->NumberOfFilters];
   ConnectivityBits = localConnectivityBits[this->NumberOfFilters];
   GeometryBits = TotalBits - ConnectivityBits;

   TotalBitrate = (double)(localTotalBits[this->NumberOfFilters]) / (double)(NumberOfPoints)/(double)NbOfFrames;
   ConnectivityBitrate = (double)(localConnectivityBits[this->NumberOfFilters]) / (double)(NumberOfPoints)/(double)NbOfFrames;
   GeometryBitrate = (double)GeometryBits / (double)(NumberOfPoints)/(double)NbOfFrames;

   cout << "j = " << 0 << ", ";
   cout << NumberOfCells << "f," << NumberOfPoints << "v, ";
   cout << "\tTotal = " << TotalBits << " bits, ";
   cout << "C = " << ConnectivityBits << " bits, ";
   cout << "G = " << GeometryBits << " bits" << endl;

   cout << "\t\tTotal = " << TotalBitrate << "b/v/f, ";
   cout << "C = " << ConnectivityBitrate << "b/v/f, ";
   cout << "G = " << GeometryBitrate << "b/v/f" << endl;


   for (i=this->NumberOfFilters - 1; i>=0; i--)
   {
      NumberOfPoints = this->Filters[i]->GetOutput()->GetNumberOfPoints();
      NumberOfCells = this->Filters[i]->GetOutput()->GetNumberOfCells();
      
      TotalBits += localTotalBits[i];
      ConnectivityBits += localConnectivityBits[i];
      GeometryBits = TotalBits - ConnectivityBits;
      
      TotalBitrate = (double) TotalBits / (double)(NumberOfPoints)/(double)NbOfFrames;
      ConnectivityBitrate = ConnectivityBits / (double)(NumberOfPoints)/(double)NbOfFrames;
      GeometryBitrate = (double)GeometryBits / (double)(NumberOfPoints)/(double)NbOfFrames;

      cout << "j = " << this->NumberOfFilters - i << ", ";
      cout << NumberOfCells << "f," << NumberOfPoints << "v, ";
      
      cout << "\tTotal = " << TotalBits << " bits, ";
      cout << "C = " << ConnectivityBits << " bits, ";
      cout << "G = " << GeometryBits << " bits" << endl;

      cout << "\t\t Total = " << TotalBitrate << "b/v/f, ";
      cout << "C = " << ConnectivityBitrate << "b/v/f, ";
      cout << "G = " << GeometryBitrate << "b/v/f" << endl;
   }
}
