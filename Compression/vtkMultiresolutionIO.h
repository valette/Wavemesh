/*=========================================================================

Program:   Mailleur 3D multi-résolution (Creatis 1997 ~)
Module:    vtkMultiresolutionIO.h
Language:  C++
Date:      2003/05
Auteurs:   Sebastien Valette
  This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME vtkMultiresolutionIO
// .SECTION Description


#ifndef __vtkMultiresolutionIO_h
#define __vtkMultiresolutionIO_h

#include <list>

#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkTimerLog.h>

#include "vtkSurface.h"
#include "vtkWaveletSubdivisionFilter.h"
#include "vtkArithmeticCoder.h"
#include "RenderWindow.h"

#define VTK_FILE_BYTE_ORDER_BIG_ENDIAN 0
#define VTK_FILE_BYTE_ORDER_LITTLE_ENDIAN 1

class VTK_EXPORT vtkMultiresolutionIO : public vtkObject
{

	struct WaveletCoefficientSet
	{
		vtkIdType FilterId;
		vtkIdType Edge;
		int Type; //=1 if the set contains descendants of edges
				  //=0 if the set contains descendants of edges except direct children

		double Max;
		int Count;
		char Sn;
		
		WaveletCoefficientSet()
		{FilterId=-1;Edge=-1;Type=-1;Max=0;Count=-1;Sn=0;};
	};

	struct WaveletCoefficient
	{
		vtkIdType FilterId,N;
		vtkIdType Edge;
		vtkIdType WaveletId;
		double Wavelet;
		int Count;
		char Sn;
	};
public:
	static vtkMultiresolutionIO *New();
	vtkTypeMacro(vtkMultiresolutionIO,vtkObject);


	// Description:
	// Returns the ith filter of the entire pyramid (no range check is performed)
	// WARNING : filters are sorted from the finest to the coarsest one
	// Exemple : GetFilter(0) will return the highest resolution filter.
	vtkWaveletSubdivisionFilter *GetFilter(int i) {return (this->Filters[i]);};

	// Description:
	// Specify file name.
	void SetFileName(const char* File) {this->FileName=File;};

	// Description:
	// Specify the mesh to encode.
	void SetInput(vtkSurface *Input) {this->Input=Input;this->Input->Register(this);};

	// Description:
	// Returns the decoded mesh.
	vtkSurface *GetOutput() {return (this->Output);};


	// Description:
	// Inserts a new Filter into the hierarchy
	void InsertNextFilter(vtkWaveletSubdivisionFilter *Filter);


	// Description:
	// Sets The Quantization
	// Only for Coding (integer Wavelets)
	void SetQuantization(int q) {this->Quantization=q;};
	
	int GetQuantization () {return this->Quantization;};

	// Description:
	// Sets The Number of bit planes
	// Only for Zerotree Coding (float wavelets)
	void SetNumberOfBitPlanes(int n) {this->NumberOfBitPlanes=n;};

	// Description:
	// Sets The Number of first bitplanes encoded altogether
	void SetNumberOfStartBitPlanes(int n) {this->NumberOfStartBitPlanes=n;};


	// Description:
	// Sets The Geometrical Constraint mode On/Off(0 : Off; 1:On)
	// Only for Coding
	void SetGeometricalConstraint(int M) {this->GeometricConstraint=M;};

	// Description:
	// Sets The Edge Angle Threshold (usually Et=0.3 OR 0<Et<1)
	// Only for Coding
	void SetEdgeAngleThreshold(double Et) {this->EdgeAngleThreshold=Et;};

	// Description:
	// Sets The Wavelet Geometrical Criterion (WGC) (Usually WGC=0.25 OR 0<WGC<1)
	// Only For Coding
	void SetWGC(double WGC) {this->WGC=WGC;};

	// Description:
	// Sets The Lifting (0=No Lifting; 1: Lifting 2:Fast 0-ring Lifting)
	// Only For Coding
	void SetLifting(int lifting);

	// Description:
	// Sets The Lifting Radius
	// Only For Coding
	void SetLiftingRadius(int radius);

	// Description:
	// Sets The Display Time for non iterative Display
	void SetDisplayTime(double Time) {this->Time=Time;};

	// Sets the arithmetic used for the wavelet decomposition.
	// 0 : float arithmetics, 1 : integer arithmetics
	// Note : for Compression, you must use integer arithmetics
	void SetArithmeticType (int Type);

	// Description:
	// Sets the capture of the reconstructed meshes to BMP file (0:off 1 :on)
	void SetCapture (int s) {this->Capture=s;};


	// Description:
	// Sets The Display (0:off 1:on 2:iterative display)
	void SetDisplay(int type) {this->Display=type;};

	// Description:
	// Sets the Display of the simplification efficiency (0: on 1: on)
	void SetDisplayEfficiency(int D) {this->DisplayEfficiency=D;};

	RenderWindow *GetRenderWindow(){return (this->MeshWindow);};

	// Description:
	// Process Analysis
	void Analyse();

	// Description:
	// Process Synthesis
	void Synthetize();

	// Description:
	// Process Geometry Approximations (Entire Pyramid)
	void Approximate();

	// Description:
	// Process Geometry Reconstructions (Entire Pyramid)
	void Reconstruct();

	// Description:
	// Writes Data to Disk. Depending on the arithmetic type, it will encode with zerotree or entropy coding
	void Write();

	// Description:
	// Reads Data from disk. Depending on the arithmetic type of the file, zerotree or entropy coding will be used
	void Read();

	// Description:
	// Sets On/off the writing of reconstructed meshes
	void SetWriteOutput(int s){this->WriteOutput=s;};

	// Description:
	// Sets the file format for meshes output : 0 : .iv 1 :.ply
	void SetFileType(int type){this->FileType=type;};

	// Description:
	// Sets On/off the writing of compression repport (0: Off 1: On)
	void SetWriteRepport(int s){this->WriteRepport=s;};

	// Description:
	// Sets On/off the on screen display of text about multiresolution analysis(0: Off 1: On)
	void SetDisplayText(int s) {this->DisplayText=s;};

	// Description:
	// Computes the tree of Abs(WaveletCoefficients), needed for zerotree coding
	void ComputeWaveletTree ();

	// Description:
	// Returns the number of filters embedded in the class
	vtkIdType GetNumberOfFilters(){return this->NumberOfFilters;};

	// Description:
	// Returns the number of filters embedded in the class
	vtkSurface *GetOutput(int Level){return this->SynthesisMeshes[Level];};

	// Set the maximum number of resolution levels
	// default value=99
	void SetMaxNumberOfLevels(int N)
	{ this->MaxNumberOfLevels=N;};

	// Defines wether butterfly dual lifting will be used
	void SetGeometryPrediction(int p)
	{this->GeometryPrediction=p;};

protected:
	vtkMultiresolutionIO();
	~vtkMultiresolutionIO();

	virtual vtkWaveletSubdivisionFilter *NewFilter(int Type)
	{ return vtkWaveletSubdivisionFilter::New();};	


	void Execute();
	int ConnectivityBytes;
	const char *FileName;
	vtkArithmeticCoder *ArithmeticCoder;

	int NumberOfBitPlanes;
	int NumberOfStartBitPlanes;
	void DisplayHires();


	void EncodeMeshConnectivity(vtkSurface *Mesh);
	vtkSurface *DecodeMeshConnectivity();

	void EncodeMeshGeometry(vtkSurface *Mesh);
	void DecodeMeshGeometry(vtkSurface *Mesh);

	void EncodeScalingFactors(vtkSurface *Mesh);
	void DecodeScalingFactors(vtkSurface *Mesh);

	void EncodeLiftingProperties();
	void DecodeLiftingProperties();

	vtkWaveletSubdivisionFilter *Filters[100];
	vtkSurface *AnalysisMeshes[100];
	vtkSurface *SynthesisMeshes[100];

	int Lifting;

	int LiftingRadius;
	int GeometryPrediction;
	vtkIdType NumberOfFilters;
	int Quantization;
	int GeometricConstraint;
	int Display;
	int DisplayEfficiency;

	int WriteOutput;
	int WriteRepport;
	int DisplayText;
	int Capture;
	vtkIdList *PointsIds;
	vtkSurface *Input;
	vtkSurface *Output;

	// the maximum of levels
	int MaxNumberOfLevels;

	void EncodeProgressiveResolution();
	void DecodeProgressiveResolution();
	void EncodeProgressivePrecision();
	void DecodeProgressivePrecision();

private:

	//BTX

	std::list<WaveletCoefficientSet> InsignificantSets[3];
	std::list<WaveletCoefficient> InsignificantCoeffs[3];
	std::list<WaveletCoefficient> SignificantCoeffs[3];


	void GetFirstGoodEdge(vtkIdType Edge,vtkIdType FilterId, vtkIdType &Edge2, vtkIdType &FilterId2, vtkIdType &Vertex);

	void EncodeSignificance(int S,int Coord);

	void EncodeCoeffRefinement(int Coeff);
	void EncodeSign(int Sign);

	int DecodeSignificance(int Coord);
	int DecodeCoeffRefinement();
	int DecodeSign();


	double StartTime;
	vtkTimerLog *Timer;

	int FileType; //0 :.iv output 1 : .ply output
	int ArithmeticType;

	// Sets The Lifting (0=No Lifting; 1: Lifting 2:Fast 0-ring Lifting)

	double Time;

	double EdgeAngleThreshold;
	double WGC;


	//  vtkMultiresolutionIO(const vtkMultiresolutionIO&);  // Not implemented.
	//  void operator=(const vtkMultiresolutionIO&);  // Not implemented.

	vtkDoubleArray *Maxima[100];
	//	vtkIdList *TreeFirstChildEdge[100];
	vtkIdList *TreeNextChildEdge[100];
	vtkIdList *TreeParentEdge[100];
	vtkIdList *EdgeMidPoints[100];
	vtkDoubleArray *Wavelets[100];
	int DataSize[100];


	RenderWindow *MeshWindow;

	vtkDoubleArray *SauvWavelets[100];

	//ETX

	void SaveWavelets();
	void LoadWavelets();
	void PrintWavelets();

	int SignificanceCodes;
	void PrintInsignificantCoeffs();

	int BaseMeshBitRate;

};



#endif
