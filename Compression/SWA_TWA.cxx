/*==================================================================================================================================

  Program:   Test module for 3D mesh sequence coding using the combination of spatial and temporal wavelet analysis (Creatis 1997 ~)
  Module:    SWA_TWA.cxx (SWA : Spatial Wavelet Transform + TWA : Temporal Wavelet Analysis)
  Language:  C++
  Date:      2006/10
  Auteurs:   Jae-Won CHO
==================================================================================================================================*/

// .NAME wavemesh 
// .SECTION Description
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"

#include "vtkSurface.h"
#include "vtkWaveletSubdivisionFilter.h"
#include "vtkMultiresolutionIOSeqII.h"
#include "math.h"


//***************************************************************************************************
//test module for 3D mesh sequence coding using the combination of spatial and temporal wavelet analysis.
//It uses the class, "vtkMultiresolutionIOSeqII".

//Input arguments :
//arg1 = the prefix of 3D dynamic mesh file name
//arg2 = index of the first frame
//arg3 = index of the last frame
//arg4 = quantization level
//arg5 = temporal wavelet decomposition level
//arg6 = temporal wavelet analysis/synthesis mode (0 : dyadic, 1 : packet)
//arg7 = temporal wavelet analysis/synthesis  filter (0 : Harr(2/2 tap), 1 : Le Gall(5/3 tap), 2 : Daubechies(9/7 tap))
//
//Example : cowheavy 0 127 12 5 0 2
//*******************************************************************************************************

char inputfilename[100];
char prefixOfFilename[100];
char outputfilename[100];
char buf[1000];

int quantize;						// quantization level
double Factor, Tx, Ty, Tz;		// To use the same scaling factors for "quantize" and "unquantize"
vtkIdType i, j, k;

vtkIdType first;					// index of the first frame
vtkIdType last;					// index of the last frame

vtkIdType tempoWaveDecompLevel;	// temporal wavelet decomposition level
vtkIdType tempoWaveMode;			// temporal wavelet analysis/synthesis mode : 0 = dyadic, 1 = packet
vtkIdType tempoWaveFilter;			// temporal wavelet analysis/synthesis filter : 0 = Haar(2/2 tap), 1 = Le Gall(5/3 tap), 2 = Daubechies(9/7 tap)

int main( int argc, char *argv[] )
{
	strcpy(prefixOfFilename,argv[1]);		//load prefix only.
	first = atoi(argv[2]);
   last = atol(argv[3]);
	quantize = atoi(argv[4]);	
	tempoWaveDecompLevel = atoi(argv[5]);

	if(tempoWaveDecompLevel != 0)				// If temporal wavelet decomposition level is 0, it means that we use only spatial wavelet analysis for each frame
	{
		tempoWaveMode = atoi(argv[6]);
		tempoWaveFilter = atoi(argv[7]);
	}

	/////Begin of Compression///////////////////////////////////////////////////////
	vtkMultiresolutionIOSeqII *MWriter=vtkMultiresolutionIOSeqII::New();
	MWriter->SetFileType(1);	// 0 : IV format; 1 : PLY format.
	MWriter->SetDisplay(0);	// 0 = No display ; 1 = Display auto ; 2 = Display interactif
	MWriter->SetWriteRepport(0);
	MWriter->SetWriteOutput(1);	// 1 = Write output mesh of each level, 0 = Do not write output mesh
	MWriter->SetDisplayText(0);
	MWriter->SetDisplayTime(0.5);
	MWriter->SetCapture(0);		// IF on decompresse, ecrit un fichier BMP pour chaque maillage
   MWriter->SetGeometryPrediction(0);	// Use prediction based on Buterfly scheme
	MWriter->SetGeometricalConstraint(0);	// 1 = use A geometry based criterion to constraint coarsening
	MWriter->SetEdgeAngleThreshold(0.3);	// IF geometry, define feature edges
	MWriter->SetWGC(0.1);	// IF geometry, define the threshold on Wavelet coefficient to decide to remove or not a point.
   MWriter->SetNumberOfBitPlanes(12);	// IF arithm = 0; number of maximum bit planes to transfert
	
   MWriter->SetArithmeticType(1);	// 0 = flottants;  1 = entiers;
	MWriter->SetLifting(2);	// 0 = No Lifting; 1 = Lifting; 2 = Fast Lifting;
	MWriter->SetLiftingRadius(1);	// IF lifting == 1; Radius of the lifting computation
	
	MWriter->SetQuantization(quantize);	
	
   MWriter->SetFileNamePrefix(prefixOfFilename);
   MWriter->SetFirstFrameId(first);
   MWriter->SetLastFrameId(last);
   MWriter->SetTempoWaveDecompLevel(tempoWaveDecompLevel);
	MWriter->SetTempoWaveMode(tempoWaveMode);
	MWriter->SetTempoWaveFilter(tempoWaveFilter);

   MWriter->EncodeSequence();
   MWriter->Delete();
	/////End of Compression///////////////////////////////////////////////////////

	/////Begin of Decompression///////////////////////////////////////////////////	
   vtkMultiresolutionIOSeqII *MReader=vtkMultiresolutionIOSeqII::New();
   MReader->SetFileNamePrefix(prefixOfFilename);
   MReader->SetFirstFrameId(first);
   MReader->SetLastFrameId(last);
   MReader->SetTempoWaveDecompLevel(tempoWaveDecompLevel);
	MReader->SetTempoWaveMode(tempoWaveMode);
	MReader->SetTempoWaveFilter(tempoWaveFilter);
   MReader->SetQuantization(quantize);
   MReader->DecodeSequence();
   MReader->Delete();	
   //////End of Decompression///////////////////////////////////////////////////

	return(0);
}
