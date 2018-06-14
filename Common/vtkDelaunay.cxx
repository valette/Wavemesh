/*=========================================================================

  Program:   Mesh to Delaunay Mesh
  Module:    vtkDelaunay.cxx
  Language:  C++
  Date:      2008/02
  Auteurs:   Sebastien VALETTE
             Arnaud Gelas

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette, Arnaud Gelas
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

#include <vtkObjectFactory.h>
#include "vtkDelaunay.h"

bool vtkDelaunay::IsEdgeFlippable(vtkIdType Edge)
{
    vtkIdType v1, v2, v3, v4, f1, f2;
    this->Mesh->GetEdgeFaces(Edge, f1, f2);
    if (f2<0)
        return false;
    this->Mesh->GetEdgeVertices(Edge, v1, v2);
    v3 = this->Mesh->GetThirdPoint(f1, v1, v2);
    v4 = this->Mesh->GetThirdPoint(f2, v1, v2);
    if (this->Mesh->IsEdge(v3,v4)) {
        return false;
    } else {
        return true;
    }
}

void vtkDelaunay::SetInputData(vtkSurface *Input)
{
	if (this->Mesh)
		this->Mesh->UnRegister(this);
	
	Input->Register(this);
	this->Mesh=Input;
}

void vtkDelaunay::DelaunayConform()
{
	// Empty the queue
	this->NonDelaunayEdges->Reset();

	// push all edges
	for( vtkIdType Edge = 0; Edge < this->Mesh->GetNumberOfEdges( ); Edge++)
	{
		if(this->Mesh->IsEdgeActive(Edge))
		{
			this->PushEdge(Edge);
		}
	}
	this->Update();
}

void vtkDelaunay::PushEdge(vtkIdType Edge)
{
	if (!this->IsEdgeFlippable(Edge))
		return;

	if(this->Mesh->IsEdgeManifold(Edge))
	{
		double t = vtkDelaunay::ComputeDelaunayConformingCriterion( this->Mesh, Edge );
		this->NonDelaunayEdges->DeleteId(Edge);
		if (t>0.00001)
			this->NonDelaunayEdges->Insert( -t, Edge );
	}	
}

void vtkDelaunay::Update()
{
	vtkIdType f1, f2;
	double t;
	
	this->NumberOfPerformedFlips=0;

	// flip the edges in the queue
	while( NonDelaunayEdges->GetNumberOfItems( ) != 0 )
	{
		vtkIdType Edge = NonDelaunayEdges->Pop( 0, t );

		vtkIdType v1, v2, v3, v4;
		this->Mesh->GetEdgeVertices( Edge, v1, v2 );
		this->Mesh->GetEdgeFaces( Edge, f1, f2 );
		v3 = this->Mesh->GetThirdPoint( f1, v1, v2 );
		v4 = this->Mesh->GetThirdPoint( f2, v1, v2 );

		if (this->Mesh->IsEdge(v3,v4) < 0)
		{
			this->Mesh->FlipEdge( Edge );
			this->NumberOfPerformedFlips++;
			this->UpdateEdge(Edge);

			this->PushEdge(this->Mesh->IsEdge( v3, v1 ));
			this->PushEdge(this->Mesh->IsEdge( v3, v4 ));
			this->PushEdge(this->Mesh->IsEdge( v3, v2 ));
			this->PushEdge(this->Mesh->IsEdge( v4, v1 ));
			this->PushEdge(this->Mesh->IsEdge( v4, v2 ));
		}
	}
}

vtkDelaunay *vtkDelaunay::New ()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject *ret = vtkObjectFactory::CreateInstance ("vtkDelaunay");
	if (ret)
	{
		return (vtkDelaunay *) ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkDelaunay);

}

vtkDelaunay::vtkDelaunay()
{
	this->Mesh=0;
	this->NonDelaunayEdges=vtkPriorityQueue::New();
}

vtkDelaunay::~vtkDelaunay()
{
	this->NonDelaunayEdges->Delete();

	if (this->Mesh)
		this->Mesh->UnRegister(this);
}

