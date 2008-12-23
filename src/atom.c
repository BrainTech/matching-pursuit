/***************************************************************************
 *   Copyright (C) 2006 by Piotr J. Durka Dobieslaw Ircha, Rafal Kus, Marek Matysiak   *
 *   durka@fuw.edu.pl, rircha@fuw.edu.pl, rkus@fuw.edu.pl				     	*
 *   Department of Biomedical Physics at Warsaw University			     		*
 *   http://brain.fuw.edu.pl, http://eeg.pl						     		*
 *												     		*
 *   This program is free software; you can redistribute it and/or modify	     		*
 *   it under the terms of the GNU General Public License as published by	     		*
 *   the Free Software Foundation; either version 2 of the License, or 		     	*
 *   (at your option) any later version.							     		*
 *												     		*
 *   This program is distributed in the hope that it will be useful,		     		*
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of	     	*
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 		*
 *   GNU General Public License for more details.					     		*
 *												     		*
 *   You should have received a copy of the GNU General Public License		     	*
 *   along with this program; if not, write to the					     		*
 *   Free Software Foundation, Inc.,							     		*
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.			     	*
 ***************************************************************************/

#define _GNU_SOURCE

#include<stdlib.h>
#include<math.h>
#include"atom.h"
#include"types.h"
#include"vector.h"

#ifdef __MINGW32__
	#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

Atom* allocateAtom(unsigned short int numberOfAllocatedChannels, unsigned char mp5Type)
{
	Atom *atom = (Atom *)malloc(sizeof(Atom));
	
	if(mp5Type & MMP1) // in case of mmp1 algorithm, the additional, last channel is assign to keep intermediate values, which helps to speed up mmp1 algotihm many times
	{
		atom->RS      = fVectorAllocate(numberOfAllocatedChannels + 1);
		atom->RC      = fVectorAllocate(numberOfAllocatedChannels + 1);
	}
	else
	{
		atom->RS      = fVectorAllocate(numberOfAllocatedChannels);
		atom->RC      = fVectorAllocate(numberOfAllocatedChannels);
	}
	
	atom->phase   = fVectorAllocate(numberOfAllocatedChannels);
	atom->feature = 0x0;

	return atom;
}

void freeAtom(Atom *atom)
{
	fVectorFree(atom->RS);
	fVectorFree(atom->RC);
	fVectorFree(atom->phase);

	free(atom);    
}

void allocateAtomElements(Atom *atom, unsigned short int numberOfAllocatedChannels, unsigned char mp5Type)
{
	if(mp5Type & MMP1) // in case of mmp1 algorithm, the additional, last channel is assign to keep intermediate values, which helps to speed up mmp1 algotihm many times
	{
		atom->RS      = fVectorAllocate(numberOfAllocatedChannels + 1);
		atom->RC      = fVectorAllocate(numberOfAllocatedChannels + 1);
	}
	else
	{
		atom->RS      = fVectorAllocate(numberOfAllocatedChannels);
		atom->RC      = fVectorAllocate(numberOfAllocatedChannels);
	}
	
	atom->phase   = fVectorAllocate(numberOfAllocatedChannels);
	atom->feature = 0x0;
}

void freeAtomElements(Atom *atom)
{
	fVectorFree(atom->RS);
	fVectorFree(atom->RC);
	fVectorFree(atom->phase);
}

void copyAtom(const Atom *sourceAtom, Atom *copyAtom, unsigned short int numberOfAllocatedChannels)
{
	unsigned short int channel;

	copyAtom->scaleIndex             = sourceAtom->scaleIndex;  
	copyAtom->position               = sourceAtom->position;
	copyAtom->rifling                = sourceAtom->rifling;
	copyAtom->randomShiftInFrequency = sourceAtom->randomShiftInFrequency;
	copyAtom->KS                     = sourceAtom->KS;
	copyAtom->KC                     = sourceAtom->KC;
	copyAtom->KM                     = sourceAtom->KM;
	
	for(channel=0;channel<numberOfAllocatedChannels;channel++)
	{
		*(copyAtom->RS + channel)    = *(sourceAtom->RS    + channel);
		*(copyAtom->RC + channel)    = *(sourceAtom->RC    + channel);
		*(copyAtom->phase + channel) = *(sourceAtom->phase + channel);
	}

	copyAtom->feature  = sourceAtom->feature;
}
