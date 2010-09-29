/*=========================================================================

Program:   Mailleur 3D multi-résolution (Creatis 1997 ~)
Module:    vtkArithmeticCoder.cxx
Language:  C++
Date:      2003/06
Auteurs:   Sebastien Valette
  This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME vtkArithmeticCoder
// .SECTION Description


#include "vtkArithmeticCoder.h"
#include "vtkObjectFactory.h"


vtkArithmeticCoder* vtkArithmeticCoder::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkArithmeticCoder");
	if(ret)
	{
		return (vtkArithmeticCoder*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkArithmeticCoder);
}

vtkArithmeticCoder::vtkArithmeticCoder()
{
}

vtkArithmeticCoder::~vtkArithmeticCoder()
{
}

void vtkArithmeticCoder::InitZerotreeQSModels(int Direction)
{
	int i,j;

	for (i=0;i<3;i++)
		for (j=0;j<4;j++)
			this->QsmSignificance[i][j].initqsmodel(2,10,100,NULL,Direction);

	this->QsmSignAndRefinement.initqsmodel(2,10,100,NULL,Direction);
}

void vtkArithmeticCoder::InitConnectivityQSModels(int Direction)
{
	this->QsmMidpoints.initqsmodel(2,10,100,NULL,Direction);
	this->QsmSwaps.initqsmodel(2,10,100,NULL,Direction);
}

void vtkArithmeticCoder::EncodeSignOrRefinement(int R)
{
	int syfreq,ltfreq;
	this->QsmSignAndRefinement.qsgetfreq(R,&syfreq,&ltfreq);
	this->Rc.encode_shift(syfreq,ltfreq,10);
	this->QsmSignAndRefinement.qsupdate(R);

}

int vtkArithmeticCoder::DecodeSignOrRefinement()
{
	int syfreq,ltfreq,ch;
	ltfreq = this->Rd.decode_culshift(10);
	ch = this->QsmSignAndRefinement.qsgetsym( ltfreq);
	this->QsmSignAndRefinement.qsgetfreq(ch,&syfreq,&ltfreq);
	this->Rd.decode_update(syfreq, ltfreq, 1<<10);
	this->QsmSignAndRefinement.qsupdate(ch);
	return (ch);

}
void vtkArithmeticCoder::EncodeSignificance(int S,int Coord,int Size)
{
	int syfreq,ltfreq;

	this->QsmSignificance[Coord][Size].qsgetfreq(S,&syfreq,&ltfreq);
	this->Rc.encode_shift(syfreq,ltfreq,10);
	this->QsmSignificance[Coord][Size].qsupdate(S);

}

int vtkArithmeticCoder::DecodeSignificance(int Coord,int Size)
{
	int syfreq,ltfreq,ch;

	ltfreq = this->Rd.decode_culshift(10);
	ch = this->QsmSignificance[Coord][Size].qsgetsym( ltfreq);
	this->QsmSignificance[Coord][Size].qsgetfreq(ch,&syfreq,&ltfreq);
	this->Rd.decode_update(syfreq, ltfreq, 1<<10);
	this->QsmSignificance[Coord][Size].qsupdate(ch);
	return (ch);

}

void vtkArithmeticCoder::EncodeMidpoint(int code)
{
	int syfreq,ltfreq;

	this->QsmMidpoints.qsgetfreq(code,&syfreq,&ltfreq);
	this->Rc.encode_shift(syfreq,ltfreq,10);
	this->QsmMidpoints.qsupdate(code);

}
int vtkArithmeticCoder::DecodeMidpoint()
{
	int syfreq,ltfreq,ch;

	ltfreq = this->Rd.decode_culshift(10);
	ch = this->QsmMidpoints.qsgetsym( ltfreq);
	this->QsmMidpoints.qsgetfreq(ch,&syfreq,&ltfreq);
	this->Rd.decode_update( syfreq, ltfreq, 1<<10);
	this->QsmMidpoints.qsupdate(ch);
	return (ch);

}

void vtkArithmeticCoder::EncodeFaceSwap(int code)
{
	int syfreq,ltfreq;

	this->QsmSwaps.qsgetfreq(code,&syfreq,&ltfreq);
	this->Rc.encode_shift(syfreq,ltfreq,10);
	this->QsmSwaps.qsupdate(code);


}
int vtkArithmeticCoder::DecodeFaceSwap()
{
	int syfreq,ltfreq,ch;
	ltfreq = this->Rd.decode_culshift(10);
	ch = this->QsmSwaps.qsgetsym(ltfreq);
	this->QsmSwaps.qsgetfreq(ch,&syfreq,&ltfreq);
	this->Rd.decode_update( syfreq, ltfreq, 1<<10);
	this->QsmSwaps.qsupdate(ch);
	return (ch);
}
