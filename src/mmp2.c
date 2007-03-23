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

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include"include/gabor.h"
#include"include/mmp2.h"
#include"include/mp5.h"
#include"include/queue.h"
#include"include/tools.h"
#include"include/types.h"

#ifdef __MINGW32__
	#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

/* MMP1 - MULTICHANNEL MATCHING PURSUIT */

void firstIterationMMP2(MP5Parameters *mp5Parameters, DataParameters *dataParameters, GaborDictionary *gaborDictionary)
{
    unsigned long int stepInToolbar, step = 0;
    unsigned short int sample,channel;
    unsigned short int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
    unsigned       int gaborsCounter;
    unsigned       int indexGabor = 0;
    int tmpNumberOfAnalysedChannels = (int)numberOfAnalysedChannels;
    int inc = 1;

    double **multiChannelSignalTable  = mp5Parameters->multiChannelSignalTable;
    double **multiChannelResidueTable = mp5Parameters->multiChannelResidueTable = mp5Parameters->multiChannelSignalTable;

    double *prevGaborTable;    
    double *signalTable, *residueTable;

    unsigned short int dimOffset   = mp5Parameters->dimOffset;
    unsigned       int dimExpand   = mp5Parameters->dimExpand;

    double t1, t2;

    double modulusesTable[numberOfAnalysedChannels];
    double sumOfModuluses;
    double bestSumOfModuluses;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
    double tmpBestModulusesTable[numberOfAnalysedChannels];

    Gabor *gabor;
    Gabor *bestGabor = allocateGabor(mp5Parameters->numberOfAnalysedChannels);

    mp5Parameters->totalSignalEnergy = 0.0;

	double tmpDataValue;
	
	for(sample=0;sample<dimOffset;sample++)
	{
		tmpDataValue = 0.0;
		
		for(channel=0;channel<numberOfChosenChannels;channel++)
			tmpDataValue+= *(*(multiChannelSignalTable + channel) + sample);

		*(meanSignalTable + sample) = tmpDataValue/numberOfAnalysedChannels;
		mp5Parameters->meanSignalEnergy+= (*(meanSignalTable + sample))*(*(meanSignalTable + sample));
	}
	
    stepInToolbar = gaborDictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;

    t1 = Clock();

    bestSumOfModuluses = 0.0;
    mp5Parameters->totalResidueEnergy = 0.0;

	signalTable = meanSignalTable;

    for(gaborsCounter=0;gaborsCounter<gaborDictionary->sizeOfDictionary;gaborsCounter++)
    {
	gabor = gaborDictionary->gaborsTable + gaborsCounter;

	if(!(gabor->feature & INCORRECTGABOR))
	{
		normSinCosGaborTable(mp5Parameters,gaborDictionary,gabor);

		findGaborDataDotProduct(mp5Parameters,gaborDictionary,gabor,signalTable,channel,FIRST_ITERATION);
    
		if(findUnknowPhaseDI(gabor,(modulusesTable + channel),channel) == ERROR)
		    continue;
	    }

    	if(findUnknowPhaseAM(gabor,modulusesTable,numberOfAnalysedChannels)==ERROR)
		continue; 

	    sumOfModuluses = 0.0;

	    for (channel= 0;channel<numberOfAnalysedChannels;channel++) 
		sumOfModuluses += *(modulusesTable+channel);  

	    if(sumOfModuluses>bestSumOfModuluses)
	    {
		bestSumOfModuluses = sumOfModuluses;
		indexGabor = gaborsCounter;
	    	memcpy((void *)tmpBestModulusesTable,(void *)modulusesTable,numberOfAnalysedChannels*sizeof(double));
	    }

	    if(dataParameters->verbose & VERBOSE_PRINT_PROCESSBAR)
		if(gaborsCounter%stepInToolbar==0)
		{
		    toolbar(step);
		    step++;
		}
	}
    }

    gabor = gaborDictionary->gaborsTable + indexGabor;
    mp5Parameters->previousGabor = gabor;
    gabor->feature|=GABOR_WAS_HIT;
    copyGabor(gabor,bestGabor,mp5Parameters->numberOfAnalysedChannels);
    memcpy((void *)bestModulusesTable,(void *)tmpBestModulusesTable,numberOfAnalysedChannels*sizeof(double));

    makeSinCosGaborTable(mp5Parameters,gaborDictionary,bestGabor);

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
        makeGaborTable(mp5Parameters,bestGabor,channel);

        residueTable   = *(multiChannelResidueTable      + channel);
        prevGaborTable = *(mp5Parameters->prevGaborTable + channel);

        findResidue(residueTable,prevGaborTable,*(bestModulusesTable + channel),dimExpand);
    }

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
        residueTable   = *(multiChannelResidueTable + channel);
	*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,dimOffset);
	mp5Parameters->totalResidueEnergy+= (*(mp5Parameters->residueEnergyInEachChannel + channel));
    }

    addNode(mp5Parameters->fitted,(void *)bestGabor);

    t2 = Clock();

    printf(" ATOM: [%3d], SEC.: %6.2lf, SIG: %6.2lf, MOD: %6.2lf, RES: %6.2lf, RES/SIG: %6.2lf",1,(t2-t1),mp5Parameters->totalSignalEnergy,ddot(&tmpNumberOfAnalysedChannels,bestModulusesTable,&inc,bestModulusesTable,&inc),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

    if(gabor->feature & DIRACDELTA)	
        printf(" D\n");
    else if(gabor->feature & GABORWAVE)	
        printf(" G\n");
    else if(gabor->feature & FFTWAVE)	
        printf(" F\n");
    
    printf("\n");
    fflush(stdout);

    if(dataParameters->verbose & VERBOSE_PRINT_FITED_GABORS)
	printFitedGabors(dataParameters,mp5Parameters,gaborDictionary,bestGabor,'F',1);

}

void nextIterationMMP2(MP5Parameters *mp5Parameters, DataParameters *dataParameters, GaborDictionary *gaborDictionary)
{
    unsigned long  int stepInToolbar, step = 0;
    unsigned short int channel;
    unsigned short int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
    unsigned short int iterationCounter;
    unsigned       int gaborsCounter;
    unsigned       int indexGabor = 0;
    int tmpNumberOfAnalysedChannels = (int)numberOfAnalysedChannels;
    int inc = 1;

    double **multiChannelResidueTable = mp5Parameters->multiChannelResidueTable;

    double *prevGaborTable;    
    double *residueTable;

    unsigned short int dimOffset   = mp5Parameters->dimOffset;
    unsigned       int dimExpand   = mp5Parameters->dimExpand;

    double t1, t2;

    double modulusesTable[numberOfAnalysedChannels];
    double sumOfModuluses;
    double bestSumOfModuluses;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
    double tmpBestModulusesTable[numberOfAnalysedChannels];

    Gabor *gabor;
    Gabor *bestGabor = NULL;

    stepInToolbar = gaborDictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;

    iterationCounter = 1;

    while(mp5Parameters->totalResidueEnergy>(mp5Parameters->totalSignalEnergy - mp5Parameters->totalSignalEnergy*mp5Parameters->energyPercent/100.0) && iterationCounter<mp5Parameters->maxNumberOfIterations)
    {
	t1 = Clock();

	bestSumOfModuluses = 0.0;
	mp5Parameters->totalResidueEnergy = 0.0;
	step = 0;
	bestGabor = allocateGabor(mp5Parameters->numberOfAnalysedChannels);

	for(gaborsCounter=0;gaborsCounter<gaborDictionary->sizeOfDictionary;gaborsCounter++)
	{
	    gabor = gaborDictionary->gaborsTable + gaborsCounter;

	    if(!(gabor->feature & GABOR_WAS_HIT) && !(gabor->feature & INCORRECTGABOR))
	    {
		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
		    prevGaborTable = *(mp5Parameters->prevGaborTable + channel);
			findGaborDataDotProduct(mp5Parameters,gaborDictionary,gabor,prevGaborTable,channel,NEXT_ITERATION);
		    if(findUnknowPhaseDI(gabor,(modulusesTable + channel),channel) == ERROR)
			continue;
		}

	        if(findUnknowPhaseAM(gabor,modulusesTable,numberOfAnalysedChannels)==ERROR)
		    continue; 

		sumOfModuluses = 0.0;
	    
		for (channel= 0;channel<numberOfAnalysedChannels;channel++) 
		    sumOfModuluses += *(modulusesTable+channel);  
	    
		if(sumOfModuluses>bestSumOfModuluses)
		{
		    bestSumOfModuluses = sumOfModuluses;
		    indexGabor = gaborsCounter;
		    memcpy((void *)tmpBestModulusesTable,(void *)modulusesTable,numberOfAnalysedChannels*sizeof(double));
		}

		if(dataParameters->verbose & VERBOSE_PRINT_PROCESSBAR)
		    if(gaborsCounter%stepInToolbar==0)
		    {
			toolbar(step);
			step++;
		    }
	    }
	}

	gabor = gaborDictionary->gaborsTable + indexGabor;
	mp5Parameters->previousGabor = gabor;
	gabor->feature|=GABOR_WAS_HIT;
	copyGabor(gabor,bestGabor,mp5Parameters->numberOfAnalysedChannels);
        memcpy((void *)bestModulusesTable,(void *)tmpBestModulusesTable,numberOfAnalysedChannels*sizeof(double));

	makeSinCosGaborTable(mp5Parameters,gaborDictionary,bestGabor);

	for(channel=0;channel<numberOfAnalysedChannels;channel++)
	{
	    makeGaborTable(mp5Parameters,bestGabor,channel);

	    prevGaborTable = *(mp5Parameters->prevGaborTable + channel);
	    residueTable   = *(multiChannelResidueTable      + channel);
	    findResidue(residueTable,prevGaborTable,*(bestModulusesTable + channel),dimExpand);
	}

	for(channel=0;channel<numberOfAnalysedChannels;channel++)
	{
	    residueTable   = *(multiChannelResidueTable + channel);
	    *(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,dimOffset);
	    mp5Parameters->totalResidueEnergy+= (*(mp5Parameters->residueEnergyInEachChannel + channel));
	}

	addNode(mp5Parameters->fitted,(void *)bestGabor);

	iterationCounter++;

	t2 = Clock();

	printf(" ATOM: [%3d], SEC.: %6.2lf, SIG: %6.2lf, MOD: %6.2lf, RES: %6.2lf, RES/SIG: %6.2lf",iterationCounter,(t2-t1),mp5Parameters->totalSignalEnergy,ddot(&tmpNumberOfAnalysedChannels,bestModulusesTable,&inc,bestModulusesTable,&inc),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

	if(gabor->feature & DIRACDELTA)	
	    printf(" D\n");
	else if(gabor->feature & GABORWAVE)	
	    printf(" G\n");
	else if(gabor->feature & FFTWAVE)	
	    printf(" F\n");
    
	printf("\n");
	fflush(stdout);

        if(dataParameters->verbose & VERBOSE_PRINT_FITED_GABORS)
	    printFitedGabors(dataParameters,mp5Parameters,gaborDictionary,bestGabor,'F',iterationCounter);
    }
}
