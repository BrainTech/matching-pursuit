/***************************************************************************
 *   Copyright (C) 2006 by Piotr J. Durka Dobieslaw Ircha, Rafal Kus       *
 *   durka@fuw.edu.pl, rircha@fuw.edu.pl, rkus@fuw.edu.pl                  *
 *   Department of Biomedical Physics at Warsaw University                 *
 *   http://brain.fuw.edu.pl, http://eeg.pl                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef _MP5_H_

	#define _MP5_H_
	
	#include"types.h"
				      
	void setNumberOfAnalysedChannelsAndNumberOfResultsFiles(MP5Parameters *mp5Parameters, DataParameters *dataParameters);
	void setMP5Parameters(MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary);
	void freeMP5Parameters(MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary);		  
	void makeSinCosExpTable(MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary);
	void normSinCosGaborTable(const MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary, Gabor *gabor);
	void makeSinCosGaborTable(MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary, const Gabor  *gabor);
        void makeGaborTable(MP5Parameters *mp5Parameters, const Gabor *gabor, unsigned short int channelNumber);
	double findSignalEnergy(const double *signalTable, unsigned short int dimOffset);
	void findResidue(double *residueTable, const double *gaborTable, double modulus, unsigned int dimExpand);
	STATUS returndAmplitudeAndModulus(const Gabor *gabor, float *amplitude, float *modulus, unsigned short int channelNumber, char *info);
	void findGaborDataDotProduct(const MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary, Gabor *currentGabor, const double *dataTable, unsigned short int channelNumber, unsigned char mode);
	STATUS findUnknowPhaseDI(Gabor *currentGabor, double *modulus, unsigned short int channelNumber);
	STATUS findUnknowPhaseAM(Gabor *currentGabor, double *modulus, unsigned short int numberOfAnalysedChannels);

#endif
