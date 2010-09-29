/*=========================================================================

Program:   Mailleur 3D multi-résolution (Creatis 1997 ~)
Module:    vtkArithmeticCoder.h
Language:  C++
Date:      2003/06
Auteurs:   Sebastien Valette
  This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME vtkArithmeticCoder
// .SECTION Description

#ifndef __vtkArithmeticCoder_h
#define __vtkArithmeticCoder_h
#include <vtkObject.h>

// RangeCoder
#include "vtkArithmeticCoderBase.h"

class VTK_EXPORT vtkArithmeticCoder : public vtkArithmeticCoderBase
{
public:

	static vtkArithmeticCoder *New();

	//////////////////////////////////////////////////////////////////////////////
	//  Wavemesh specific methods
	//////////////////////////////////////////////////////////////////////////////

	/// initialize the models for Zerotree Encoding
	void InitZerotreeQSModels(int Direction);

	/// initialize the models for the Connectivity coding
	void InitConnectivityQSModels(int Direction);

	/// Encode a bit of significance
	void EncodeSignificance(int S,int Coord,int Size);

	/// Encode a bit of sign or refinment
	void EncodeSignOrRefinement(int R);

	/// Decode a significance bit
	int DecodeSignificance(int Coord, int Size);

	/// decode a sign or refinement bit
	int DecodeSignOrRefinement();

	/// encode a bit for midpoint detection
	void EncodeMidpoint(int code);

	/// decode a bit for midpoint detection
	int DecodeMidpoint();

	/// encode a face swap bit
	void EncodeFaceSwap(int code);

	/// decode a face swap bit
	int DecodeFaceSwap();

protected:

	vtkArithmeticCoder();
	~vtkArithmeticCoder();

	////////////////////////////////////////////////////////////////////////////
	//  Wavemesh specific members
	////////////////////////////////////////////////////////////////////////////

	qsmodel QsmSignAndRefinement;
	qsmodel QsmSignificance[3][4];

	qsmodel QsmMidpoints;
	qsmodel QsmSwaps;
};

#endif
