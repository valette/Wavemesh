//  This software is governed by the GPL license (see License.txt)
#include "vtkSparseMatrix.h"

vtkSparseMatrix::~vtkSparseMatrix()	
{
	int i;
	NonZeroElement *next,*last;
	for (i=0;i<(int) this->FirstHor.size();i++)
	{
		next=this->FirstHor[i];
		while (next!=0)
		{
			last=next;
			next=next->NextHor;
			delete last;
		}
	}
};

vtkSparseMatrix* vtkSparseMatrix::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkSparseMatrix");
	if(ret)
	{
		return (vtkSparseMatrix*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkSparseMatrix);
}

void vtkSparseMatrix::Init(int i,int j)
{
	this->FirstHor.resize(i);
	this->FirstVer.resize(j);
}

void vtkSparseMatrix::AddValue(int i, int j, double Value)
{
	
int found_place,new_place;
NonZeroElement *new_non_zero,*last,*next;

	found_place=0;
	new_place=1;
	last=0;
	next=this->FirstHor[i];
	if (next==0)
	{
		new_non_zero=new NonZeroElement;
		this->FirstHor[i]=new_non_zero;
		new_non_zero->NextHor=0;
		new_non_zero->i=i;
		new_non_zero->j=j;
		new_non_zero->Value=Value;
	}
	else
	{
		while (found_place==0)
		{
			if (next->j>j)
			{
				new_non_zero=new NonZeroElement;
				if (new_non_zero==0) 
					cout<<"probleme mémoire"<<endl;
				new_non_zero->Value=Value;
				new_non_zero->i=i;
				new_non_zero->j=j;
				new_non_zero->NextHor=next;
				found_place=1;
				if (last==0)
					this->FirstHor[i]=new_non_zero;
				else
					last->NextHor=new_non_zero;
			}
			else
			{
				if (next->j==j)
				{
					next->Value+=Value;
					found_place=1;
					new_place=0;
				}
				else
				{
					last=next;
					next=next->NextHor;
					if (next==0)
					{
						new_non_zero=new  NonZeroElement;
						if (new_non_zero==0)
							cout<<"problème mémoire"<<endl;;
						new_non_zero->j=j;
						new_non_zero->i=i;
						new_non_zero->Value=Value;
						new_non_zero->NextHor=0;
						last->NextHor=new_non_zero;
						found_place=1;
					}
				}
			}
		}
	}
	if (new_place==1)
	{
		last=0;
		found_place=0;
		next=this->FirstVer[j];
		if (next==0)
		{
			this->FirstVer[j]=new_non_zero;
			new_non_zero->NextVer=0;
		}
		else
		{
			while (found_place==0)
			{
				if (next->i>i)
				{
					new_non_zero->NextVer=next;
					found_place=1;
					if (last==0)
						this->FirstVer[j]=new_non_zero;
					else
						last->NextVer=new_non_zero;
				}
				else
				{
					last=next;
					next=next->NextVer;
					if (next==0)
					{
						new_non_zero->NextVer=0;
						last->NextVer=new_non_zero;
						found_place=1;
					}
				}
			}
		}
	}
}
