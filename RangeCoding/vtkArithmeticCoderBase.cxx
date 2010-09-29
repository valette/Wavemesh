/*=========================================================================

Program:   Mailleur 3D multi-résolution (Creatis 1997 ~)
Module:    vtkArithmeticCoderBase.cxx
Language:  C++
Date:      2003/06
Auteurs:   Sebastien Valette
  This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME vtkArithmeticCoderBase
// .SECTION Description

#include "vtkArithmeticCoderBase.h"
#include "vtkObjectFactory.h"

vtkArithmeticCoderBase* vtkArithmeticCoderBase::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkArithmeticCoder");
	if(ret)
	{
		return (vtkArithmeticCoderBase*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkArithmeticCoderBase);
}

void vtkArithmeticCoderBase::Encode(int code,qsmodel *Model)
{
	int syfreq,ltfreq;
	Model->qsgetfreq(code,&syfreq,&ltfreq);
	this->Rc.encode_shift(syfreq,ltfreq,15);
	Model->qsupdate(code);
}

int vtkArithmeticCoderBase::Decode(qsmodel *Model)
{
	int syfreq,ltfreq,ch;
	ltfreq = this->Rd.decode_culshift(15);
	ch = Model->qsgetsym( ltfreq);
	Model->qsgetfreq(ch,&syfreq,&ltfreq);
	this->Rd.decode_update(syfreq, ltfreq, 1<<15);
	Model->qsupdate(ch);
	return (ch);
}

void vtkArithmeticCoderBase::OpenFile(const char* file,int direction)
{
	this->Direction=direction;
	if (direction==1)
		this->Rc.open_file(file);
	else
		this->Rd.open_file(file);
};

void vtkArithmeticCoderBase::CloseFile()
{
	if (this->Direction==1)
		this->Rc.close_file();
	else
		this->Rd.close_file();
};

vtkArithmeticCoderBase::vtkArithmeticCoderBase()
{
}

vtkArithmeticCoderBase::~vtkArithmeticCoderBase()
{
}

void vtkArithmeticCoderBase::StartCoding()
{
	this->Rc.start_encoding(0,0);
}

void vtkArithmeticCoderBase::StartDecoding()
{
	this->Rd.start_decoding();
}

void vtkArithmeticCoderBase::EncodeByte( unsigned char B)
{
	this->Rc.encode_byte(B);
}

unsigned char vtkArithmeticCoderBase::DecodeByte()
{
	return(this->Rd.decode_byte());
}

void vtkArithmeticCoderBase::EncodeWord(int W)
{
	this->Rc.encode_short(W);
}

int vtkArithmeticCoderBase::DecodeWord()
{
	return(this->Rd.decode_short());
}

void vtkArithmeticCoderBase::EncodeInt(int I)
{
	int Buff;
	Buff=I>>16;
	this->Rc.encode_short(Buff);
	Buff=I&((1<<16)-1);
	this->Rc.encode_short(Buff);
}

int vtkArithmeticCoderBase::DecodeInt()
{
	int number1=this->Rd.decode_short();
	int number2=this->Rd.decode_short();
	return((number1<<16)+number2);
}

void vtkArithmeticCoderBase::EncodeBit(bool b)
{
	this->Rc.encode_shift((freq)1,(freq)(b),(freq)1);
}

bool vtkArithmeticCoderBase::DecodeBit()
{
	bool tmp = (this->Rd.decode_culshift(1)!=0);
	this->Rd.decode_update(1,tmp,(freq)1<<1);
	return ( tmp);
}

void vtkArithmeticCoderBase::EncodeFloat( float F)
{
	unsigned char *Bytes;
	Bytes=(unsigned char *) &F;

	for (unsigned int i=0;i<sizeof (float);i++)
		this->Rc.encode_byte(Bytes[i]);
}

float vtkArithmeticCoderBase::DecodeFloat()
{
	float F;
	unsigned char *Bytes;
	Bytes=(unsigned char *) &F;

	for (unsigned int i=0;i<sizeof (float);i++)
		Bytes[i]=this->Rd.decode_byte();
	return (F);
}

int vtkArithmeticCoderBase::StopCoding()
{
	return (this->Rc.done_encoding());
}

void vtkArithmeticCoderBase::StopDecoding()
{
	this->Rd.done_decoding();
}

