/*=========================================================================

  Program:   vtkTag
  Module:    vtkSurface
  Language:  C++
  Date:      2008/04
  Auteur:    Sebastien VALETTE

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

#ifndef __vtkTag_h
#define __vtkTag_h

#include <limits.h>
#include "vtkObjectFactory.h"

class VTK_EXPORT vtkTag : public vtkObject
{

public:

	/// The Constructor vtkRecons::New();
	static vtkTag *New()
	{
		// First try to create the object from the vtkObjectFactory
		vtkObject *ret = vtkObjectFactory::CreateInstance ("vtkTag");
		if (ret)
		{
			return (vtkTag *) ret;
		}
		// If the factory was unable to create the object, then create it here.
		return (new vtkTag);
	}

	vtkTypeMacro(vtkTag,vtkObject);

	void SetNumberOfItems(int NumberOfItems)
	{
		this->LastTag->SetNumberOfValues(NumberOfItems);
		this->ResetCounter();
	}
	
	void Reset()
	{
		if (this->Time==INT_MAX)
			this->ResetCounter();
		else
			this->Time++;
	}
	
	bool IsTagged(vtkIdType Item)
	{
		return (this->LastTag->GetValue(Item)==this->Time);
	}
	
	void Tag(vtkIdType Item)
	{
		this->LastTag->SetValue(Item,this->Time);
	}

	void UnTag(vtkIdType Item)
	{
		this->LastTag->SetValue(Item,this->Time-1);
	}	

protected:

	void ResetCounter()
	{
		int NumberOfItems=this->LastTag->GetSize();
		this->Time=0;
		for (int i=0;i<NumberOfItems;i++)
			this->LastTag->SetValue(i,-1);
	}

	vtkIntArray *LastTag;
	
	int Time;

	/// constructor
	vtkTag()
	{
		this->LastTag=vtkIntArray::New();
	}

	/// desctructor
	virtual ~vtkTag()
	{
		this->LastTag->Delete();
	}
};

class VTK_EXPORT vtkTagWithList : public vtkTag
{

public:

	static vtkTagWithList *New()
	/// The Constructor vtkRecons::New();
	{
		// First try to create the object from the vtkObjectFactory
		vtkObject *ret = vtkObjectFactory::CreateInstance ("vtkTag");
		if (ret)
		{
			return (vtkTagWithList *) ret;
		}
		// If the factory was unable to create the object, then create it here.
		return (new vtkTagWithList);
	}

	vtkTypeMacro(vtkTagWithList,vtkObject);
	vtkGetObjectMacro(TaggedItems, vtkIdList)

	void Reset()
	{
		this->TaggedItems->Reset();
		this->vtkTag::Reset();
	}
	
	bool IsTagged(vtkIdType Item)
	{
		return (this->LastTag->GetValue(Item)==this->Time);
	}
	
	void Tag(vtkIdType Item)
	{
		if (!this->IsTagged(Item))
			this->TaggedItems->InsertNextId(Item);
		this->vtkTag::Tag(Item);
	}

protected:

	vtkIdList *TaggedItems;
	
	/// constructor
	vtkTagWithList()
	{
		this->TaggedItems=vtkIdList::New();
	}

	/// desctructor
	virtual ~vtkTagWithList()
	{
		this->TaggedItems->Delete();
	}
};

#endif
