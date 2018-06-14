/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOFFReader.h,v $
  Language:  C++

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

// .NAME vtkOFFReader - Read an ASCII file formated to OFF format.
// .SECTION Description

#ifndef __vtkOFFReader_h
#define __vtkOFFReader_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>

#define VTK_FILE_BYTE_ORDER_BIG_ENDIAN 0
#define VTK_FILE_BYTE_ORDER_LITTLE_ENDIAN 1

class VTK_EXPORT vtkOFFReader : public vtkPolyDataAlgorithm
{
public:
  static vtkOFFReader *New();
  vtkTypeMacro(vtkOFFReader,vtkPolyDataAlgorithm);

  // Description:
  // Specify file name.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  void SetInfoOnCellsOn() {this->InfoOnCells = 1;};
  void SetInfoOnCellsOff() {this->InfoOnCells = 0;};

protected:
  vtkOFFReader();
  ~vtkOFFReader();

  char *FileName;
  int SwapBytes;

  int NumberOfPoints;
  int NumberOfCells;
  
  void ExecuteInformation();
  void Execute();

private:
  vtkOFFReader(const vtkOFFReader&);  // Not implemented.
  void operator=(const vtkOFFReader&);  // Not implemented.

  int InfoOnCells;

};

#endif
