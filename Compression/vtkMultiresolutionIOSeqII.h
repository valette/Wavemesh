/*==================================================================================================================
  Program:   3D mesh sequence coding using the combination of spatial and temporal wavelet analysis (Creatis 1997 ~)
  Module:    vtkMultiresolutionIOSeqII.h
  Language:  C++
  Date:      2006/10
  Auteurs:   Jae-Won Cho
  Base codeur: vtkMultiresolutionIO.h (Sebastien Valette), vtkMultiresolutionIOSeq.h (Min-Su Kim and Jae-Won Cho)
==================================================================================================================*/
// .NAME vtkMultiresolutionIOSeqII
// .SECTION Description

#ifndef __vtkMultiresolutionIOSeqII_h
#define __vtkMultiresolutionIOSeqII_h

#include "vtkMultiresolutionIO.h"

class VTK_EXPORT vtkMultiresolutionIOSeqII : public vtkMultiresolutionIO
{
public:

   vtkTypeMacro(vtkMultiresolutionIOSeqII,vtkMultiresolutionIO);

   static vtkMultiresolutionIOSeqII *New();
   
   vtkMultiresolutionIOSeqII(); 
   ~vtkMultiresolutionIOSeqII();
   
   void SetFileNamePrefix(char *prefix){ strcpy(this->prefixOfFilename, prefix);};	// set prefix of filename

   void SetFirstFrameId(int i) { this->first = i; };		// set the first frame index

   void SetLastFrameId(int i) { this->last = i; };			// set the last frame index

   void SetTempoWaveDecompLevel(int i){ this->tempoWaveDecompLevel = i; };		// set the temporal wavelet decomposition level

	void SetTempoWaveMode(int i){ this->tempoWaveMode = i; };						// set the temporal wavelet analysis/synthesis mode (0 : dyadic, 1 : packet)

	void SetTempoWaveFilter(int i){ this->tempoWaveFilter = i; };					// set the temporal wavelet filter bank (0 : Haar, 1 : Le Gall, 2 : Daubechies)

   void WriteIFrame(int frameindex);				// write the first framme

   void WritePFrame(int frameindex);				// write the other frames

   void ReadIFrame(int frameindex);					// read the first frame

   void ReadPFrame(int frameindex);					// read the other frames

   void UpdateGeometry(vtkSurface *Mesh);			// for the other frame coding

   void CopyForEncoding(int frameindex, int mode);					// for storing spatial wavelet coefficients into buffer for temporal wavelet transform(when encoding)

	void CopyForDecoding(int frameindex, int mode);					// for storing temporal wavelet coefficients into buffer for temporal inverse wavelet transform(when decoding)

   void EncodeSequence();								// encoding sequence

   void DecodeSequence();								// decoding sequence

   void GetGlobalBitrate();							// calculate global bitrate

   void GetLocalBitrate();								// calculate local bitrate

	void TemporalWaveletDecomposition();			// temporal wavelet decomposition

	void TemporalWaveletReconstruction();			// temporal wavelet reconstruction

	void oneD_DWT_Dyadic(double *InputData, vtkIdType data_leng, vtkIdType level);		// dyadic wavelet

	void oneD_IDWT_Dyadic(double *InputData, vtkIdType data_leng, vtkIdType level);		// dyadic inverse wavelet

	void oneD_DWT_Packet(double *InputData, vtkIdType data_leng, vtkIdType level);		// packet wavelet

	void oneD_IDWT_Packet(double *InputData, vtkIdType data_leng, vtkIdType level);		// packet inverse wavelet

	void fHarr(double *InputData, vtkIdType sub_data_leng);										// wavelet filter 0 : haar

	void iHarr(double *InputData, vtkIdType sub_data_leng);										// wavelet filter 0 : haar

	void forward_lifting_53(double *InputData, vtkIdType sub_data_leng);						// wavelet filter 1 : 5/3 tap lifting

	void inverse_lifting_53(double *InputData, vtkIdType sub_data_leng);						// wavelet filter 1 : 5/3 tap lifting

	void forward_lifting_97(double *InputData, vtkIdType sub_data_leng);						// wavelet filter 2 : 9/7 tap lifting

	void inverse_lifting_97(double *InputData, vtkIdType sub_data_leng);						// wavelet filter 2 : 9/7 tap lifting

private:

   void Encode(int frameindex, int mode);				// encode frame

   void Decode(int frameindex, int mode);				// decode frame

	vtkPoints *BaseMeshEncodingBuffer[8192];			// base mesh buffer for encoding (maximum frame length : 8192)

	vtkPoints *BaseMeshDecodingBuffer[8192];			// base mesh buffer for decoding

   vtkSurface *BaseMesh;

	vtkIntArray *WaveletsEncodingBuffer[8192][100];	// wavelet coefficients buffer for encoding (maximum frame length : 8192, maximum spatial resolution level : 100)

	vtkIntArray *WaveletsDecodingBuffer[8192][100];	// wavelet coefficients buffer for decoding

   char inputfilename[1000];								
   char prefixOfFilename[1000];
   char outputfilename[1000];
   char buf[1000];

   vtkIdType first;					// first index of the input sequence
   vtkIdType last;					// last index of the input sequence
   int tempoWaveDecompLevel;		// temporal wavelet decomposition level
	int tempoWaveMode;				// temporal wavelet mode : 0 = dyadic, 1 = packet
	int tempoWaveFilter;				// temporal wavelet filter : 0 = harr, 1 = 5/3 tap lifting, 2 = 9/7 tap lifting
   int subfixOfFilename;
   double Factor;
   double Tx;
   double Ty;
   double Tz;

   // for evaluation //
   int globalConnectivityBits;
   int globalGeometryBits;
   int globalTotalBits;

   int localConnectivityBits[100];
   int localGeometryBits[100];
   int localTotalBits[100];

   int NumOfPoints;
   int FileSize;

};

#endif
