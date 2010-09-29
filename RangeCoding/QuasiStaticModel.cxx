/*
qsmodel.c     headerfile for quasistatic probability model

(c) Michael Schindler
1997, 1998, 2000
http://www.compressconsult.com/
michael@compressconsult.com

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.  It may be that this
program violates local patents in your country, however it is
belived (NO WARRANTY!) to be patent-free here in Austria.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA.

Qsmodel is a quasistatic probability model that periodically
(at chooseable intervals) updates probabilities of symbols;
it also allows to initialize probabilities. Updating is done more
frequent in the beginning, so it adapts very fast even without
initialisation.

it provides function for creation, deletion, query for probabilities
and symbols and model updating.

for usage see example.c
*/
#include "QuasiStaticModel.h"

/* default tablesize 1<<TBLSHIFT */
#define TBLSHIFT 7

/* rescale frequency counts */
void qsmodel::dorescale()
{
	int i, cf, missing;
	if (this->nextleft)  /* we have some more before actual rescaling */
	{
		this->incr++;
		this->left = this->nextleft;
		this->nextleft = 0;
		return;
	}

	if (this->rescale < this->targetrescale)  /* double rescale interval if needed */
	{
		this->rescale <<= 1;
		if (this->rescale > this->targetrescale)
			this->rescale = this->targetrescale;
	}

	cf = missing = this->cf[this->n];  /* do actual rescaling */

	for (i=this->n-1; i; i--)
	{
		int tmp = this->newf[i];
		cf -= tmp;
		this->cf[i] = cf;
		tmp = tmp>>1 | 1;
		missing -= tmp;
		this->newf[i] = tmp;
	}

	if (cf!=this->newf[0])
	{   
		this->deleteqsmodel();
		std::cout<<"BUG: rescaling left %d total frequency"<<std::endl;
	}

	this->newf[0] = this->newf[0]>>1 | 1;
	missing -= this->newf[0];
	this->incr = missing / this->rescale;
	this->nextleft = missing % this->rescale;
	this->left = this->rescale - this->nextleft;
	if (this->search != NULL)
	{
		i=this->n;
		while (i)
		{   int start, end;
			end = (this->cf[i]-1) >> this->searchshift;
			i--;
			start = this->cf[i] >> this->searchshift;
			while (start<=end)
			{
				this->search[start] = i;
				start++;
			}
		}
	}
}

/* initialisation of qsmodel                           */
/* m   qsmodel to be initialized                       */
/* n   number of symbols in that model                 */
/* lg_totf  base2 log of total frequency count         */
/* rescale  desired rescaling interval, should be < 1<<(lg_totf+1) */
/* init  array of int's to be used for initialisation (NULL ok) */
/* compress  set to 1 on compression, 0 on decompression */
void qsmodel::initqsmodel(int n, int lg_totf, int rescale, int *init, int compress )
{  
	this->deleteqsmodel();
	this->n = n;
	this->targetrescale = rescale;
	this->searchshift = lg_totf - TBLSHIFT;
	if (this->searchshift < 0)
		this->searchshift = 0;
	
	this->cf = new uint2[n+1];//malloc(test);
	this->newf = new int[n+1]; //)*sizeof(uint2));
	this->cf[n] = 1<<lg_totf;
	this->cf[0] = 0;
	if (compress)
		this->search = NULL;
	else
	{
		this->search = new int[(1<<TBLSHIFT)+1];
		this->search[1<<TBLSHIFT] = n-1;
	}
	this->resetqsmodel(init);
}

/* reinitialisation of qsmodel                         */
/* m   qsmodel to be initialized                       */
/* init  array of int's to be used for initialisation (NULL ok) */
void qsmodel::resetqsmodel(int *init)
{
	int i, end, initval;
	this->rescale = this->n>>4 | 2;
	this->nextleft = 0;
	if (init == NULL)
	{
		initval = this->cf[this->n] / this->n;
		end = this->cf[this->n] % this->n;
		for (i=0; i<end; i++)
			this->newf[i] = initval+1;
		for (; i<this->n; i++)
			this->newf[i] = initval;
	}
	else
		for(i=0; i<this->n; i++)
			this->newf[i] = init[i];
	this->dorescale();
}

/* deletion of qsmodel m                               */
void qsmodel::deleteqsmodel()
{
	if (this->cf==0)
		return;

	delete [] this->cf;
	delete [] this->newf;
	if (this->search != NULL)
		delete [] this->search;
}

/* retrieval of estimated frequencies for a symbol     */
/* m   qsmodel to be questioned                        */
/* sym  symbol for which data is desired; must be <n   */
/* sy_f frequency of that symbol                       */
/* lt_f frequency of all smaller symbols together      */
/* the total frequency is 1<<lg_totf                   */
void qsmodel::qsgetfreq(int sym, int *sy_f, int *lt_f )
{ 
	*sy_f = this->cf[sym+1] - (*lt_f = this->cf[sym]);
}	

/* find out symbol for a given cumulative frequency    */
/* m   qsmodel to be questioned                        */
/* lt_f  cumulative frequency                          */
int qsmodel::qsgetsym(int lt_f )
{
	int lo, hi;
	int *tmp;
	tmp = &this->search[(lt_f>>this->searchshift)];
	lo = *tmp;
	hi = *(tmp+1) + 1;
	while (lo+1 < hi )
	{
		int mid = (lo+hi)>>1;
		if (lt_f < this->cf[mid])
			hi = mid;
		else
			lo = mid;
	}
	return lo;
}

/* update model                                        */
/* m   qsmodel to be updated                           */
/* sym  symbol that occurred (must be <n from init)    */
void qsmodel::qsupdate(int sym )
{
	if (this->left <= 0)
		dorescale();
	this->left--;
	this->newf[sym] += this->incr;
}

qsmodel::qsmodel()
{
	this->cf=0;
}

qsmodel::~qsmodel()
{
	this->deleteqsmodel();
}

