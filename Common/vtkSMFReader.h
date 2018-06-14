/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSMFReader.h,v $
  Language:  C++

  Copyright (c) 2003 Arnaud Gelas, Min-Su KIM 
  All rights reserved.

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Arnaud Gelas, Min-Su KIM
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

// .NAME vtkSMFReader - read Qsilm .smf files
// .SECTION Description
// vtkSMFReader is a source object that reads Qsilm .smf
// files. The output of this source object is polygonal data.
// .SECTION See Also
// vtkSMFImporter

#ifndef __vtkSMFReader_h
#define __vtkSMFReader_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>

class VTK_EXPORT vtkSMFReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkSMFReader *New();
  
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of Qsilm .smf file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  vtkSMFReader();
  ~vtkSMFReader();
  
  void Execute();

  char *FileName;
private:
  vtkSMFReader(const vtkSMFReader&);  // Not implemented.
  void operator=(const vtkSMFReader&);  // Not implemented.

	int AddNormalVector(char *line, vtkFloatArray *NormalVector);
	int AddTextureCoordinate(char *line, vtkFloatArray *TextureTuple);
	int AddColorComponent(char *line,vtkFloatArray *colorTuple);
	int AddFace(char *line, vtkCellArray *polys);
	int AddVertex(char *line, vtkPoints *points);
};

#endif


