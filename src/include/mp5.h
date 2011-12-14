/*************************************************************************************
 *   Copyright (C) 2006 by Piotr J. Durka Dobieslaw Ircha, Rafal Kus, Marek Matysiak *
 *   durka@fuw.edu.pl, rircha@fuw.edu.pl, rkus@fuw.edu.pl				     	     *
 *   Department of Biomedical Physics at Warsaw University			     		     *
 *   http://brain.fuw.edu.pl, http://eeg.pl						     		         *
 *												    								 *
 *   This program is free software; you can redistribute it and/or modify			 *
 *   it under the terms of the GNU General Public License as published by			 *
 *   the Free Software Foundation; either version 2 of the License, or				 *
 *   (at your option) any later version.											 *
 *												     								 *
 *   This program is distributed in the hope that it will be useful,	     		 *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of	     			 *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 					 *
 *   GNU General Public License for more details.					   		   		 *
 *												     								 *
 *   You should have received a copy of the GNU General Public License		     	 *
 *   along with this program; if not, write to the					     			 *
 *   Free Software Foundation, Inc.,							    				 *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.			 			 *
 *************************************************************************************/

#ifndef _MP5_H_

	#define _MP5_H_

	#include"types.h"


	void setNumberOfAnalysedChannelsAndNumberOfResultsFiles(MP5Parameters *mp5Parameters);
	void setMP5Parameters(const Dictionary *dictionary, MP5Parameters *mp5Parameters);
	void freeMP5Parameters(const Dictionary *dictionary, MP5Parameters *mp5Parameters);
	void makeSinCosExpTable(const Dictionary *dictionary, MP5Parameters *mp5Parameters);
	void normAtomTable(const Dictionary *dictionary, const MP5Parameters *mp5Parameters, Atom *atom);
	void makeSinCosExpAtomTable(const Dictionary *dictionary, MP5Parameters *mp5Parameters, const Atom *atom);
	void makeAtomTable(MP5Parameters *mp5Parameters, const Atom *atom, unsigned short int channelNumber);
	double findSignalEnergy(const double *signalTable, unsigned int offsetExpandedDimension);
	void findResidue(double *residueTable, const double *atomTable, double modulus, unsigned int offsetExpandedDimension);
	STATUS returndAmplitudeAndModulus(const Atom *atom, float *amplitude, float *modulus, unsigned short int channelNumber, char *info);
	void findAtomDataDotProduct(const Dictionary *dictionary, const MP5Parameters *mp5Parameters, Atom *currentAtom, const double *dataTable, unsigned int channelNumber, unsigned char mode);
	STATUS findUnknowPhaseDI(Atom *currentAtom, double *modulus, unsigned int channelNumber);
	STATUS findUnknowPhaseAM(Atom *currentAtom, double *modulus, unsigned int numberOfAnalysedChannels);
	void findGaborDataDotProductFFT(const Dictionary *dictionary, const MP5Parameters *mp5Parameters, Atom *currentGabor, const double *dataTable, unsigned short int channelNumber, unsigned char mode);
	double findLambda(double residueEnergy);

#endif
