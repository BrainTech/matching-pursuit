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

#ifndef _ATOM_H_

	#define _ATOM_H_

	#include"types.h"

	Atom* allocateAtom(unsigned short int numberOfAllocatedChannels, unsigned char mp5Type);
	void freeAtom(Atom *atom);
	void allocateAtomElements(Atom *atom, unsigned short int numberOfAllocatedChannels, unsigned char mp5Type);
	void freeAtomElements(Atom *atom);
	void copyAtom(const Atom *sourceAtom, Atom *copyAtom, unsigned short int numberOfAllocatedChannels);
	void printFitedAtoms(const MP5Parameters *mp5Parameters, const Dictionary *dictionary, Atom *atom, char where, unsigned short int orderIndex);

#endif

