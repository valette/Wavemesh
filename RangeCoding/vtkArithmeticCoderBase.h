/*=========================================================================

Program:   Mailleur 3D multi-résolution (Creatis 1997 ~)
Module:    vtkArithmeticCoderBase.h
Language:  C++
Date:      2003/06
Auteurs:   Sebastien Valette
  This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME vtkArithmeticCoderBase
// .SECTION Description

#ifndef __vtkArithmeticCoderBase_h
#define __vtkArithmeticCoderBase_h
#include <iostream>
#include <fstream>
#include <vtkObject.h>

// RangeCoder
#include "port.h"
#include "QuasiStaticModel.h"
#include "RangeEncoder.h"
#include "RangeDecoder2.h"

//////////////////////////////////////////////////////////////////////////////////////
// This class is a general purpose arithmetic coder. It uses the RangeCoder coded by
// Michael Schindler (www.compressconsult.com). He kindly allowed us to freely use his
// code inside our project
//
///////////////////////////////////////////////////////////////////////////////////////
class VTK_EXPORT vtkArithmeticCoderBase : public vtkObject
{
public:
	static vtkArithmeticCoderBase *New();

	//Description:
	// Opens the file before starting IO access
	//Direction=0 for reading, 1 for writing
	void OpenFile(const char* file,int direction);

	/// closes the file
	void CloseFile();

	/// Begin Coding
	void StartCoding();

	/// Stop Coding
	void StartDecoding();

	/// stops the coding. Returns the number of bytes written
	int StopCoding();

	/// stops the decoding
	void StopDecoding();

	////////////////////////////////////////////////////////////////////////////////////
	// The following methods are usefull for encoding values without a model
	// Thus, these methods do not compress the data. They just store them...
	////////////////////////////////////////////////////////////////////////////////////

	/// encode a byte
	void EncodeByte(unsigned char B);

	/// decode a byte and returns it
	unsigned char DecodeByte();

	/// encode a word
	void EncodeWord(int W);

	/// decode a word and returns it
	int DecodeWord();
	
	/// encode an int
	void EncodeInt(int I);

	/// decode an int and returns it
	int DecodeInt();

	/// encode a bit
	void EncodeBit(bool b);

	/// decode 
	bool DecodeBit();

	/// Encode a floating point value (lossless)
	void EncodeFloat(float F);

	/// Decode a floating point value (lossless)
	float DecodeFloat();

	//////////////////////////////////////////////////////////////////////////////
	// To achieve real compression, use these methods :)
	//////////////////////////////////////////////////////////////////////////////

	/// Encode a code with a specific QuasiStaticModel
	void Encode(int code,qsmodel *Model);

	/// Decode a code with a specific QuasiStaticModel
	int Decode(qsmodel *Model);
	
	rangedecoder2 *GetDecoder()
		{return &this->Rd;};

protected:

	vtkArithmeticCoderBase();
	~vtkArithmeticCoderBase();

	/// the direction of file IO
	/// 0: read 1: write
	int Direction;

	/**
	** Rangecoder from Michael Schindler
	**/
	rangeencoder Rc;
	rangedecoder2 Rd;
};

#endif
