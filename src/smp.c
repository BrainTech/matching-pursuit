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

#include<math.h>
#include<stdlib.h>
#include"include/def.h"
#include"include/gabor.h"
#include"include/mp5.h"
#include"include/queue.h"
#include"include/smp.h"
#include"include/tools.h"
#include"include/types.h"

/* SMP - SINGLE CHANNEL MATCHING PURSUIT */

void firstIterationSMP(MP5Parameters *mp5Parameters, const DataParameters *dataParameters, GaborDictionary *gaborDictionary)
{
    unsigned long int stepInToolbar, step = 0;
    unsigned int gaborsCounter;
    unsigned int bestGaborIndex = 0;

    double *signalTable    = mp5Parameters->singleChannelSignalTable;
    double *residueTable   = mp5Parameters->singleChannelResidueTable = mp5Parameters->singleChannelSignalTable;

    unsigned short int dimOffset   = mp5Parameters->dimOffset;
    unsigned       int dimExpand   = mp5Parameters->dimExpand;

    double t1, t2;

    double modulus, absOfModulus;
    double bestModulus = 0, absOfBestModulus;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;

    Gabor *gabor;

    Gabor *bestGabor = allocateGabor(mp5Parameters->numberOfAnalysedChannels);

    t1 = Clock();

    *(mp5Parameters->signalEnergyInEachChannel) = findSignalEnergy(signalTable,dimOffset);

    absOfBestModulus = 0.0;

    stepInToolbar = gaborDictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;

    for(gaborsCounter=0;gaborsCounter<gaborDictionary->sizeOfDictionary;gaborsCounter++)
    {
	gabor = gaborDictionary->gaborsTable + gaborsCounter;

	if(!(gabor->feature & INCORRECTGABOR))
	{

	    normSinCosGaborTable(mp5Parameters,gaborDictionary,gabor);
	    findGaborDataDotProduct(mp5Parameters,gaborDictionary,gabor,signalTable,0,FIRST_ITERATION);

	    if(findUnknowPhaseDI(gabor,&modulus,0) == ERROR)
		continue;

	    absOfModulus = fabs(modulus);

	    if(absOfModulus>absOfBestModulus)
	    {
	        bestModulus      = modulus;
	        absOfBestModulus = fabs(bestModulus);
	        bestGaborIndex = gaborsCounter;
	    }
	}
	
	if(dataParameters->verbose & VERBOSE_PRINT_PROCESSBAR)
	    if(gaborsCounter%stepInToolbar==0)
	    {
		toolbar(step);
		step++;
	    }
    }

    gabor = gaborDictionary->gaborsTable + bestGaborIndex;    
    mp5Parameters->previousGabor = gabor;
    *bestModulusesTable = bestModulus;

    gabor->feature|=GABOR_WAS_HIT;

    copyGabor(gabor,bestGabor,mp5Parameters->numberOfAnalysedChannels);

    makeSinCosGaborTable(mp5Parameters,gaborDictionary,bestGabor);

    makeGaborTable(mp5Parameters,bestGabor,0);

    findResidue(residueTable,*(mp5Parameters->prevGaborTable),*bestModulusesTable,dimExpand);
    *(mp5Parameters->residueEnergyInEachChannel) = findSignalEnergy(signalTable,dimOffset);

    addNode(mp5Parameters->fitted,(void *)bestGabor);

    t2 = Clock();
	
    printf(" ATOM: [%-3d], SEC.: %-6.2lf, SIG: %-6.2lf, MOD: %-6.2lf, RES: %-6.2lf, RES/SIG: %-6.2lf %% ",1,(t2-t1),*(mp5Parameters->signalEnergyInEachChannel),pow(*(mp5Parameters->bestModulusesTable),2.0),*(mp5Parameters->residueEnergyInEachChannel),((*(mp5Parameters->residueEnergyInEachChannel))/(*(mp5Parameters->signalEnergyInEachChannel)))*100.0);

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

void nextIterationSMP(MP5Parameters *mp5Parameters, const DataParameters *dataParameters, GaborDictionary *gaborDictionary)
{
    unsigned long  int stepInToolbar, step = 0;
    unsigned short int iterationsCounter;
    unsigned       int gaborsCounter;
    unsigned int bestGaborIndex = 0;
  
    double *residueTable = mp5Parameters->singleChannelResidueTable;

    unsigned short int dimOffset   = mp5Parameters->dimOffset;
    unsigned       int dimExpand   = mp5Parameters->dimExpand;

    double t1, t2;

    double modulus, absOfModulus;

    double bestModulus = 0, absOfBestModulus;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;

    Gabor *gabor, *bestGabor;

    stepInToolbar = gaborDictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;

    iterationsCounter = 1;

    while((*(mp5Parameters->residueEnergyInEachChannel)>(*(mp5Parameters->signalEnergyInEachChannel) - (*(mp5Parameters->signalEnergyInEachChannel))*mp5Parameters->energyPercent/100.0)) && iterationsCounter<mp5Parameters->maxNumberOfIterations)
    {
	step = 0;
	t1 = Clock();
	absOfBestModulus = 0.0;
	bestGabor = allocateGabor(mp5Parameters->numberOfAnalysedChannels);

	for(gaborsCounter=0;gaborsCounter<gaborDictionary->sizeOfDictionary;gaborsCounter++)
	{	
	    gabor = gaborDictionary->gaborsTable + gaborsCounter;

	    if(!(gabor->feature & GABOR_WAS_HIT) && !(gabor->feature & INCORRECTGABOR)) 
	    {

		findGaborDataDotProduct(mp5Parameters,gaborDictionary,gabor,*(mp5Parameters->prevGaborTable),0,NEXT_ITERATION);

		if(findUnknowPhaseDI(gabor,&modulus,0) == ERROR)
		    continue;

		absOfModulus = fabs(modulus);

		if(absOfModulus>absOfBestModulus)
		{
		    bestModulus      = modulus;
		    absOfBestModulus = fabs(bestModulus);
		    bestGaborIndex = gaborsCounter;
		}
	    }

	    if(dataParameters->verbose & VERBOSE_PRINT_PROCESSBAR)
	        if(gaborsCounter%stepInToolbar==0)
		{
		    toolbar(step);
		    step++;
		}
	}

	gabor = gaborDictionary->gaborsTable + bestGaborIndex;
	mp5Parameters->previousGabor = gabor;
	*bestModulusesTable = bestModulus;

	gabor->feature|=GABOR_WAS_HIT;

	copyGabor(gabor,bestGabor,mp5Parameters->numberOfAnalysedChannels);

	makeSinCosGaborTable(mp5Parameters,gaborDictionary,bestGabor);

	makeGaborTable(mp5Parameters,bestGabor,0);

	findResidue(residueTable,*(mp5Parameters->prevGaborTable),*bestModulusesTable,dimExpand);
	
	*(mp5Parameters->residueEnergyInEachChannel) = findSignalEnergy(residueTable,dimOffset);

	addNode(mp5Parameters->fitted,(void *)bestGabor);

	iterationsCounter++;

	t2 = Clock();
	
    
	printf(" ATOM: [%-3d], SEC.: %-6.2lf, SIG: %-6.2lf, MOD: %-6.2lf, RES: %-6.2lf, RES/SIG: %-6.2lf %% ",iterationsCounter,(t2-t1),*(mp5Parameters->signalEnergyInEachChannel),pow((*(mp5Parameters->bestModulusesTable)),2.0),*(mp5Parameters->residueEnergyInEachChannel),((*(mp5Parameters->residueEnergyInEachChannel))/(*(mp5Parameters->signalEnergyInEachChannel))*100.0));

	if(gabor->feature & DIRACDELTA)	
	    printf(" D\n");
	else if(gabor->feature & GABORWAVE)	
	    printf(" G\n");
	else if(gabor->feature & FFTWAVE)	
	    printf(" F\n");
	printf("\n");

	fflush(stdout);
	
	if(dataParameters->verbose & VERBOSE_PRINT_FITED_GABORS)
	    printFitedGabors(dataParameters,mp5Parameters,gaborDictionary,bestGabor,'F',(unsigned short int)iterationsCounter);
	

    }
}
