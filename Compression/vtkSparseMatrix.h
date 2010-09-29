/*=========================================================================

  Program:   Classe pour matrices creuses (Creatis 1997/2002)
  Module:    vtkSparseMatrix.h
  Language:  C++
  Date:      2002/03
  Auteurs:   Sebastien Valette
  This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME vtkSparseMatrix 
// .SECTION Description
#include <vector>
#include <vtkCommand.h>
#include <vtkPoints.h>
#include <vtkObjectFactory.h>
#include <vtkDataObject.h>

#ifndef __vtkSparseMatrix_h
#define __vtkSPArseMatrix_h
//----------------------------------------------------------------------
/**
 ** DECLARATION
 **
 **/

class VTK_EXPORT vtkSparseMatrix : public vtkObject
//class vtkSparseMatrix :public vtkObject
{
public:

	class NonZeroElement
	{
	public:

		double Value;
		int i,j;
		NonZeroElement *NextHor, *NextVer;

		NonZeroElement () {this->NextHor=0; this->NextVer=0;};

	};

	static vtkSparseMatrix *New();
    vtkTypeMacro(vtkSparseMatrix,vtkObject);


	// Description:
	// Create a similar type object.
    //	vtkDataObject *MakeObject() {return vtkSparseMatrix::New();};

	//BTX
	std::vector <NonZeroElement*> FirstHor;
	std::vector <NonZeroElement*> FirstVer;
	//ETX

	void Init(int i,int j);
	void AddValue(int i, int j, double Value);
	void Multiply(vtkPoints *in, vtkPoints *out) {};
	void Multiply(vtkSparseMatrix *in, vtkSparseMatrix *out) {};

	vtkSparseMatrix(const vtkSparseMatrix&) {};
	void operator=(const vtkSparseMatrix&) {};

protected:

	vtkSparseMatrix() {}; //constructeur
	~vtkSparseMatrix(); //Destructeur 

private:


};

#endif

