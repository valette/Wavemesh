/*=========================================================================

  Program:   SMF Writer
  Module:    vtkSMFWriter.cxx
  Language:  C++
  Date:      2003/05
  Auteurs:   Sebastien Valette
=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

// .NAME vtkSMFWriter - write SMF files
// .SECTION Description
// vtkSMFWriter writes SMF files (.smf) files.

#ifndef __vtkSMFWriter_h
#define __vtkSMFWriter_h

#include <vtkPolyDataWriter.h>

class VTK_EXPORT vtkSMFWriter : public vtkPolyDataWriter
{
public:
	static vtkSMFWriter *New();
	vtkTypeMacro(vtkSMFWriter,vtkPolyDataWriter);
	//virtual void PrintSelf(ostream& os, vtkIndent indent);

protected:
	vtkSMFWriter() {};
	~vtkSMFWriter() {};

	void WriteData()
	{
		std::ofstream File;
		File.open (this->FileName, ofstream::out | ofstream::trunc);

		vtkIdType i;
		vtkIdType nverts=this->GetInput()->GetNumberOfPoints();
		vtkIdType nfaces=this->GetInput()->GetNumberOfCells();

		File << "# Generated from vtkSMFWriter" << endl;
		File<< "# " <<  nverts<< " vertices" << endl;
		File << "# " << nfaces << " faces" << endl;

		double P[3];
		for(i=0; i<nverts; i++)
		{
			this->GetInput()->GetPoint(i,P);
			File << "v "<< P[0] << " " << P[1] << " " << P[2] << endl;
		}

		vtkIdType  *Vertices,NumberOfVertices;
		vtkIdType  j;
		for(i=0; i<nfaces; i++)
		{
			this->GetInput()->GetCellPoints(i,NumberOfVertices,Vertices);
			File << "f ";
			for (j=0;j<NumberOfVertices;j++)
				File<< Vertices[j]+1 << " ";
			File<< endl;
		}
		File.close();
	}

	private:
	vtkSMFWriter(const vtkSMFWriter&);  // Not implemented.
	void operator=(const vtkSMFWriter&);  // Not implemented.
};
vtkStandardNewMacro(vtkSMFWriter);
#endif

