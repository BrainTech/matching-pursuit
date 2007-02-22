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

#define _GNU_SOURCE

#include<stdlib.h>
#include<math.h>
#include"include/gabor.h"
#include"include/types.h"
#include"include/vector.h"

#ifdef __MINGW32__
#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

Gabor* allocateGabor(unsigned short int numberOfAnalysedChannels)
{
	Gabor *gabor = (Gabor *)malloc(sizeof(Gabor));

	gabor->RS      = fVectorAllocate(numberOfAnalysedChannels);
	gabor->RC      = fVectorAllocate(numberOfAnalysedChannels);
	gabor->phase   = fVectorAllocate(numberOfAnalysedChannels);
	gabor->feature = 0x0;

	return gabor;
}

void freeGabor(Gabor *gabor)
{
	fVectorFree(gabor->RS);
	fVectorFree(gabor->RC);
	fVectorFree(gabor->phase);

	free(gabor);    
}

void allocateGaborElements(Gabor *gabor, unsigned short int numberOfAnalysedChannels)
{
	gabor->RS      = fVectorAllocate(numberOfAnalysedChannels);
	gabor->RC      = fVectorAllocate(numberOfAnalysedChannels);
	gabor->phase   = fVectorAllocate(numberOfAnalysedChannels);
	gabor->feature = 0x0;
}

void freeGaborElements(Gabor *gabor)
{
	fVectorFree(gabor->RS);
	fVectorFree(gabor->RC);
	fVectorFree(gabor->phase);
}

void copyGabor(const Gabor *sourceGabor, Gabor *copyGabor, unsigned short int numberOfAnalysedChannels)
{
	unsigned short int channel;

	copyGabor->scaleIndex = sourceGabor->scaleIndex;  
	copyGabor->position   = sourceGabor->position;
	copyGabor->rifling    = sourceGabor->rifling;
	copyGabor->KS         = sourceGabor->KS;
	copyGabor->KC         = sourceGabor->KC;
	copyGabor->KM         = sourceGabor->KM;
    
	for(channel=0;channel<numberOfAnalysedChannels;channel++)
	{
		*(copyGabor->RS + channel)    = *(sourceGabor->RS    + channel);
		*(copyGabor->RC + channel)    = *(sourceGabor->RC    + channel);
		*(copyGabor->phase + channel) = *(sourceGabor->phase + channel);
	}

	copyGabor->feature  = sourceGabor->feature;
}

void printFitedGabors(const DataParameters *dataParameters, const MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary, Gabor *gabor, char where, unsigned short int orderIndex)
{

	unsigned long int sample;
	unsigned short int channel;
	
	double sinPhase, cosPhase;
	double sinPart, cosPart, sinCosPart;
	double phase, amplitude, modulus;

	double KS, KC, KM;
	double RS, RC;
	
	FILE *fileOut = NULL;

	if(where == 'M')  // print gabors on monitor screen
	    fileOut = stdout;
	else if(where == 'F')
	    fileOut = dataParameters->fitedGaborsFile;

	fprintf(fileOut,"\n IN ORDER OF FITED: %hu \n",orderIndex);

	fprintf(fileOut,"  position:        %hu\n",gabor->position);

	if(gabor->feature & DIRACDELTA)
	    fprintf(fileOut,"  scale:           0\n");
	else if(gabor->feature & FFTWAVE)
	    fprintf(fileOut,"  scale:           %u\n",(unsigned int)(mp5Parameters->dimOffset + 1));
	else
    	    fprintf(fileOut,"  scaleIndex:      %hu, scale: %hu\n",gabor->scaleIndex,*(gaborDictionary->tableOfScalesInOptimalDictionary + gabor->scaleIndex));

//	fprintf(fileOut,"  frequency:       %lf\n",gabor->frequency);
	fprintf(fileOut,"  basicFrequency:  %lf\n",*(gaborDictionary->tableOfFrequenciesInOptimalDictionary + gabor->scaleIndex));
	fprintf(fileOut,"  basicPeriod:     %u\n",*(gaborDictionary->tableOfPeriodsInOptimalDictionary + gabor->scaleIndex));
	
	for(sample=0;sample<*(gaborDictionary->tableOfPeriodsInOptimalDictionary + gabor->scaleIndex);sample++)
		fprintf(fileOut,"  %lf ",*(*(mp5Parameters->sinTable + gabor->scaleIndex)+ sample));
	fprintf(fileOut,"\n");
	
	for(sample=0;sample<*(gaborDictionary->tableOfPeriodsInOptimalDictionary + gabor->scaleIndex);sample++)
		fprintf(fileOut,"  %lf ",*(*(mp5Parameters->cosTable + gabor->scaleIndex)+ sample));
	fprintf(fileOut,"\n");

	fprintf(fileOut,"  rifling:         %u\n",gabor->rifling);
	fprintf(fileOut,"  KS:      %1.16lf  \n",(gabor->KS));
	fprintf(fileOut,"  KC:      %1.16lf  \n",(gabor->KC));
	fprintf(fileOut,"  KM:      %1.16lf  \n",(gabor->KM));
	
	for (channel = 0;channel<mp5Parameters->numberOfAnalysedChannels;channel++)
	{
		KS = gabor->KS;
		KC = gabor->KC;
		KM = gabor->KM;
	    
		RS = *(gabor->RS    + channel);
		RC = *(gabor->RC    + channel);

		phase = *(gabor->phase + channel);

		sincos(phase,&sinPhase,&cosPhase);
		sinPart    = KS*(sinPhase*sinPhase);
		cosPart    = KC*(cosPhase*cosPhase);
		sinCosPart = 2.0*KM*sinPhase*cosPhase;

		amplitude = sqrt((sinPart + cosPart ) - sinCosPart);
		modulus   = (RC*cosPhase - RS*sinPhase)/amplitude;

		fprintf(fileOut,"  channel:    %-hu      \n",channel);
		fprintf(fileOut,"  RS:         %-1.16lf  \n",RS);
		fprintf(fileOut,"  RC:         %-1.16lf  \n",RC);
		fprintf(fileOut,"  phase:      %-1.16lf  \n",phase);
		fprintf(fileOut,"  modulus:    %-23.16lf \n",modulus);
		fprintf(fileOut,"  amplitude:  %-23.16lf \n",amplitude);
	}
}
