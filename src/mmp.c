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


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include"atom.h"
#include"mp5.h"
#include"queue.h"
#include"tools.h"
#include"types.h"
#include"vector.h"

extern unsigned char applicationMode;

/* MMP - MULTICHANNEL MATCHING PURSUIT */

static BOOLEAN getBestModulusesTable(unsigned char MMPType, 
									 Atom *atom,
									 double *modulusesTable,
									 double *tmpBestModulusesTable,
									 unsigned short int numberOfAnalysedChannels,
									 double *bestSumOfModuluses)
{
	int inc = 1;
	double sumOfModuluses = 0.0; // sum of moduluses or sum of sqr moduluses

	if(MMPType & MMP1)
	{		
		if(findUnknowPhaseAM(atom,modulusesTable,numberOfAnalysedChannels)==ERROR)
			return FALSE;

		sumOfModuluses = dasum(numberOfAnalysedChannels,modulusesTable,inc);
	}
	else if(MMPType & MMP3)
		sumOfModuluses = ddot(numberOfAnalysedChannels,modulusesTable,inc,modulusesTable,inc);
	else if(MMPType & MMP4)
		sumOfModuluses = *(modulusesTable + 0);

	if(sumOfModuluses>*bestSumOfModuluses)
	{
		*bestSumOfModuluses = sumOfModuluses;
		memcpy((void *)tmpBestModulusesTable,(void *)modulusesTable,numberOfAnalysedChannels*sizeof(double));
		return TRUE;
	}
	return FALSE;
}

void firstIterationMMP(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned long int stepInToolbar, step = 1;
    unsigned short int channel;
    unsigned short int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
    unsigned       int atomsCounter;
    unsigned       int indexAtom = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;   
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
	const int inc = 1;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    double **multiChannelSignalTable  = mp5Parameters->multiChannelSignalTable;
    double **multiChannelResidueTable = mp5Parameters->multiChannelResidueTable = mp5Parameters->multiChannelSignalTable;

    double *prevAtomTable;    
    double *signalTable, *residueTable;

    const unsigned int offsetExpandedDimension = mp5Parameters->offsetExpandedDimension;
	const unsigned int numberOfNonFFTAtoms     = dictionary->numberOfNonFFTAtoms;
	
    double t1, t2;

    double modulusesTable[numberOfAnalysedChannels];
    double bestSumOfModuluses;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
    float  *bestPhasesTable    = mp5Parameters->bestPhasesTable;
 
    double tmpBestModulusesTable[numberOfAnalysedChannels];

    Atom *atom;
    Atom *bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

    mp5Parameters->totalSignalEnergy = 0.0;

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
		signalTable = *(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);
		*(mp5Parameters->signalEnergyInEachChannel + channel) = findSignalEnergy(signalTable,offsetExpandedDimension);
		mp5Parameters->totalSignalEnergy = mp5Parameters->totalSignalEnergy + (*(mp5Parameters->signalEnergyInEachChannel + channel));
    }

    stepInToolbar = dictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;

    t1 = Clock();
	
	if(applicationMode & PROCESS_SERVER_MODE)
	{
		printf("ATOM\t%3u\t%3u\t%6.2lf\t%6.2lf\n",0,dictionary->sizeOfDictionary,0.0,0.0);
		fflush(stdout);
	}

    bestSumOfModuluses = 0.0;
    mp5Parameters->totalResidueEnergy = 0.0;

	atomsCounter = 0;
	atom = dictionary->atomsTable;
	
	if(mp5Parameters->FFT)
	{
		if(dictionary->diracInDictionary || dictionary->gaussInDictionary)
		{
			for(atomsCounter=0;atomsCounter<numberOfNonFFTAtoms;atomsCounter++)
			{
				if(!(atom->feature & INCORRECTGABOR))
				{
					normAtomTable(dictionary,mp5Parameters,atom);

					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{
						signalTable = *(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);

						findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
		
						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
							continue;
					} 

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					if((atomsCounter+1)%stepInToolbar==0)
					{
						if(applicationMode & PROCESS_USER_MODE)
						{
							if(mp5Parameters->progressBar)
							{
								toolbar(step);
								step++;
							}
						}
						else
						{
							printf("TESTED %u\n",atomsCounter);
							fflush(stdout);
						}
					}
				}
				atom++;
			}
		}
		
		/* from now the code is FFT for SIN/COS waves */
		if(dictionary->sinCosInDictionary)
		{
			for(channel=0;channel<numberOfAnalysedChannels;channel++)
			{
				signalTable = *(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);
				findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
			}

			numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

			for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
			{
				if(!(atom->feature & INCORRECTGABOR))
				{
					normAtomTable(dictionary,mp5Parameters,atom);

					for(channel=0;channel<numberOfAnalysedChannels;channel++)
						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
							continue;

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					if((atomsCounter+1)%stepInToolbar==0)
					{
						if(applicationMode & PROCESS_USER_MODE)
						{
							if(mp5Parameters->progressBar)
							{
								toolbar(step);
								step++;
							}
						}
						else
						{
							printf("TESTED %u\n",atomsCounter);
							fflush(stdout);
						}
					}
				}
				atom++;
				atomsCounter++;
			}
		}
		/* end of FFT code for SIN/COS waves */

		for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++)
		{
			numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);
			
			for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
			{
				for(channel=0;channel<numberOfAnalysedChannels;channel++)
				{
					signalTable = *(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);
					findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
				}

				numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

				for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
				{
					if(!(atom->feature & INCORRECTGABOR))
					{
						normAtomTable(dictionary,mp5Parameters,atom);

						for(channel=0;channel<numberOfAnalysedChannels;channel++)
							if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
								continue;

						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						if((atomsCounter+1)%stepInToolbar==0)
						{
							if(applicationMode & PROCESS_USER_MODE)
							{
								if(mp5Parameters->progressBar)
								{
									toolbar(step);
									step++;
								}
							}
							else
							{
								printf("TESTED %u\n",atomsCounter);
								fflush(stdout);
							}
						}
					}
					atom++;
					atomsCounter++;
				}
			}
		}
	}
	else
	{
		for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
		{
			atom = dictionary->atomsTable + atomsCounter;

			if(!(atom->feature & INCORRECTGABOR))
			{
				normAtomTable(dictionary,mp5Parameters,atom);

				for(channel=0;channel<numberOfAnalysedChannels;channel++)
				{
					signalTable = *(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);

					findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
		
					if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
						continue;
				} 

				if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
					indexAtom = atomsCounter;

				if((atomsCounter+1)%stepInToolbar==0)
				{
					if(applicationMode & PROCESS_USER_MODE)
					{
						if(mp5Parameters->progressBar)
						{
							toolbar(step);
							step++;
						}
					}
					else
					{
						printf("TESTED %u\n",atomsCounter);
						fflush(stdout);
					}
				}
			}
		}
    }

    atom = dictionary->atomsTable + indexAtom;
    mp5Parameters->previousAtom = atom;
    atom->feature|=GABOR_WAS_HIT;
    copyAtom(atom,bestAtom,mp5Parameters->numberOfAnalysedChannels);
    memcpy((void *)bestModulusesTable,(void *)tmpBestModulusesTable,numberOfAnalysedChannels*sizeof(double));
    memcpy((void *)bestPhasesTable,(void *)(bestAtom->phase),numberOfAnalysedChannels*sizeof(float));

    makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
        makeAtomTable(mp5Parameters,bestAtom,channel);
		residueTable   =  *(multiChannelResidueTable + mp5Parameters->chosenChannels[channel] - 1);
        prevAtomTable = *(mp5Parameters->prevAtomTable + channel);

        findResidue(residueTable,prevAtomTable,*(bestModulusesTable + channel),offsetExpandedDimension);
    }

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
		residueTable = *(multiChannelResidueTable + mp5Parameters->chosenChannels[channel] - 1);
		*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,offsetExpandedDimension);
		mp5Parameters->totalResidueEnergy+= (*(mp5Parameters->residueEnergyInEachChannel + channel));
    }

    addNode(mp5Parameters->fitted,(void *)bestAtom);

    t2 = Clock();

	if(applicationMode & PROCESS_USER_MODE)
	{
		printf(" ATOM: [%3d], SEC.: %6.2lf, SIG: %6.2lf, MOD: %6.2lf, RES: %6.2lf, RES/SIG: %6.2lf",1,(t2-t1),mp5Parameters->totalSignalEnergy,ddot(numberOfAnalysedChannels,bestModulusesTable,inc,bestModulusesTable,inc),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

		if(atom->feature & DIRACDELTA)	
			printf(" D\n");
		else if(atom->feature & GAUSSFUNCTION)
			printf(" N\n");
		else if(atom->feature & SINCOSWAVE)	
			printf(" H\n");
		else if(atom->feature & GABORWAVE)	
			printf(" G\n");

		printf("\n");
		fflush(stdout);
	}
}

void nextIterationMMP(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned long  int stepInToolbar, step = 1;
    unsigned short int channel = 0;
    unsigned short int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
	unsigned short int lastChannel = numberOfAnalysedChannels; 
    unsigned short int iterationCounter;
    unsigned       int atomsCounter;
    unsigned       int indexAtom = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;   
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
    const int inc = 1;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    double **multiChannelResidueTable = mp5Parameters->multiChannelResidueTable;

    double *prevAtomTable;    
    double *residueTable;

    const unsigned       int offsetExpandedDimension = mp5Parameters->offsetExpandedDimension;
	const unsigned       int numberOfNonFFTAtoms     = dictionary->numberOfNonFFTAtoms;
	
    double t1, t2;
	double energyProgress, iterationProgress;

    double modulusesTable[numberOfAnalysedChannels];
    double bestSumOfModuluses;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
    float  *bestPhasesTable    = mp5Parameters->bestPhasesTable;
    double tmpBestModulusesTable[numberOfAnalysedChannels];
	double RC = 0.0;
	double RS = 0.0;
	
    Atom *atom;
    Atom *bestAtom = NULL;
	
    stepInToolbar = dictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;

    iterationCounter = 1;

    while(mp5Parameters->totalResidueEnergy>(mp5Parameters->totalSignalEnergy - mp5Parameters->totalSignalEnergy*mp5Parameters->energyPercent/100.0) && iterationCounter<mp5Parameters->maximalNumberOfIterations)
    {
		if(applicationMode & PROCESS_SERVER_MODE)
		{
			energyProgress    = (((mp5Parameters->totalSignalEnergy - mp5Parameters->totalResidueEnergy))/mp5Parameters->totalSignalEnergy)*100.0;
			iterationProgress = (iterationCounter/mp5Parameters->maximalNumberOfIterations)*100.0;
			printf("ATOM\t%3u\t%3u\t%6.2lf\t%6.2lf \n",iterationCounter,dictionary->sizeOfDictionary,energyProgress,iterationProgress);			
			fflush(stdout);
		}

		t1 = Clock();

		step = 1;
		atomsCounter = 0;
		bestSumOfModuluses = 0.0;
		mp5Parameters->totalResidueEnergy = 0.0;
		bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

		atom = dictionary->atomsTable;
		
		if(mp5Parameters->FFT)
		{
			if(dictionary->diracInDictionary || dictionary->gaussInDictionary)
			{
				for(atomsCounter=0;atomsCounter<numberOfNonFFTAtoms;atomsCounter++)
				{
					if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
					{
						if(mp5Parameters->MPType & MMP1)
						{
							prevAtomTable = *mp5Parameters->prevAtomTable;
							findAtomDataDotProduct(dictionary,mp5Parameters,atom,prevAtomTable,lastChannel,MMP1_NEXT_ITERATION);
				
							RS = *(atom->RS + lastChannel);
							RC = *(atom->RC + lastChannel);
						}
						else
						{
							for(channel=0;channel<numberOfAnalysedChannels;channel++)
							{	
								prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
								findAtomDataDotProduct(dictionary,mp5Parameters,atom,prevAtomTable,channel,NEXT_ITERATION);
							}
						}
						
						for(channel=0;channel<numberOfAnalysedChannels;channel++)
						{	
							if(mp5Parameters->MPType & MMP1)
							{
								if(((*(bestPhasesTable))*(*(bestPhasesTable + channel))) > 0 )
								{
									*(atom->RS + channel) = (*(atom->RS + channel)) - (*(bestModulusesTable + channel))*RS;
									*(atom->RC + channel) = (*(atom->RC + channel)) - (*(bestModulusesTable + channel))*RC;
								}
								else
								{
									*(atom->RS + channel) = (*(atom->RS + channel)) + (*(bestModulusesTable + channel))*RS;
									*(atom->RC + channel) = (*(atom->RC + channel)) + (*(bestModulusesTable + channel))*RC;
								}
							}
							
							if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
								continue;
						}

						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						if((atomsCounter+1)%stepInToolbar==0)
						{
							if(applicationMode & PROCESS_USER_MODE)
							{
								if(mp5Parameters->progressBar)
								{
									toolbar(step);
									step++;
								}
							}
							else
							{
								printf("TESTED %u\n",atomsCounter);
								fflush(stdout);
							}
						}	
					}
					atom++;
				}
			}
			
			/* from now the code is FFT for SIN/COS waves */
			if(dictionary->sinCosInDictionary)
			{
				if(mp5Parameters->MPType & MMP1)
				{
					prevAtomTable = *mp5Parameters->prevAtomTable;
					findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,prevAtomTable,lastChannel,MMP1_NEXT_ITERATION);
				}
				else
				{
					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{
						prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
						findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,prevAtomTable,channel,NEXT_ITERATION);
					}
				}
				
				numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

				for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
				{
					if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
					{
						if(mp5Parameters->MPType & MMP1)
						{
							RS = *(atom->RS + lastChannel);
							RC = *(atom->RC + lastChannel);
						}
						
						for(channel=0;channel<numberOfAnalysedChannels;channel++)
						{	
							if(mp5Parameters->MPType & MMP1)
							{
								if(((*(bestPhasesTable))*(*(bestPhasesTable + channel))) > 0 )
								{
									*(atom->RS + channel) = (*(atom->RS + channel)) - (*(bestModulusesTable + channel))*RS;
									*(atom->RC + channel) = (*(atom->RC + channel)) - (*(bestModulusesTable + channel))*RC;
								}
								else
								{
									*(atom->RS + channel) = (*(atom->RS + channel)) + (*(bestModulusesTable + channel))*RS;
									*(atom->RC + channel) = (*(atom->RC + channel)) + (*(bestModulusesTable + channel))*RC;
								}
							}

							if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
								continue;
						}
	
						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						if((atomsCounter+1)%stepInToolbar==0)
						{
							if(applicationMode & PROCESS_USER_MODE)
							{
								if(mp5Parameters->progressBar)
								{
									toolbar(step);
									step++;
								}
							}
							else
							{
								printf("TESTED %u\n",atomsCounter);
								fflush(stdout);
							}
						}	
					}
					atom++;
					atomsCounter++;
				}
			}
			/* end of FFT code for SIN/COS waves */
			
			for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++)
			{
				numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);
				
				for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
				{
					if(mp5Parameters->MPType & MMP1)
					{
						prevAtomTable = *mp5Parameters->prevAtomTable;
						findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,prevAtomTable,lastChannel,MMP1_NEXT_ITERATION);
					}
					else
					{
						for(channel=0;channel<numberOfAnalysedChannels;channel++)
						{
							prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
							findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,prevAtomTable,channel,NEXT_ITERATION);
						}
					}
					
					numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

					for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
					{
						if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
						{
							if(mp5Parameters->MPType & MMP1)
							{
								RS = *(atom->RS + lastChannel);
								RC = *(atom->RC + lastChannel);
							}
							
							for(channel=0;channel<numberOfAnalysedChannels;channel++)
							{	
								if(mp5Parameters->MPType & MMP1)
								{
									if(((*(bestPhasesTable))*(*(bestPhasesTable + channel))) > 0 )
									{
										*(atom->RS + channel) = (*(atom->RS + channel)) - (*(bestModulusesTable + channel))*RS;
										*(atom->RC + channel) = (*(atom->RC + channel)) - (*(bestModulusesTable + channel))*RC;
									}
									else
									{
										*(atom->RS + channel) = (*(atom->RS + channel)) + (*(bestModulusesTable + channel))*RS;
										*(atom->RC + channel) = (*(atom->RC + channel)) + (*(bestModulusesTable + channel))*RC;
									}
								}
								
								if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
									continue;
							}
	
							if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
								indexAtom = atomsCounter;
								
							if((atomsCounter+1)%stepInToolbar==0)
							{
								if(applicationMode & PROCESS_USER_MODE)
								{
									if(mp5Parameters->progressBar)
									{
										toolbar(step);
										step++;
									}
								}
								else
								{
									printf("TESTED %u\n",atomsCounter);
									fflush(stdout);
								}
							}	
						}
						atom++;
						atomsCounter++;
					}
				}
			}
		}
		else
		{		
			for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
			{
				atom = dictionary->atomsTable + atomsCounter;

				if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
				{			
					if(mp5Parameters->MPType & MMP1)
					{
						prevAtomTable = *mp5Parameters->prevAtomTable;
						findAtomDataDotProduct(dictionary,mp5Parameters,atom,prevAtomTable,lastChannel,MMP1_NEXT_ITERATION);
				
						RS = *(atom->RS + lastChannel);
						RC = *(atom->RC + lastChannel);
					}
					else
					{
						for(channel=0;channel<numberOfAnalysedChannels;channel++)
						{
							prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
							findAtomDataDotProduct(dictionary,mp5Parameters,atom,prevAtomTable,channel,NEXT_ITERATION);
						}
					}
						
					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{	
						if(mp5Parameters->MPType & MMP1)
						{
							if(((*(bestPhasesTable))*(*(bestPhasesTable + channel))) > 0)
							{
								*(atom->RS + channel) = (*(atom->RS + channel)) - (*(bestModulusesTable + channel))*RS;
								*(atom->RC + channel) = (*(atom->RC + channel)) - (*(bestModulusesTable + channel))*RC;
							}
							else
							{
								*(atom->RS + channel) = (*(atom->RS + channel)) + (*(bestModulusesTable + channel))*RS;
								*(atom->RC + channel) = (*(atom->RC + channel)) + (*(bestModulusesTable + channel))*RC;
							}
						}
						
						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
							continue;
					}

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;
				
					if((atomsCounter+1)%stepInToolbar==0)
					{
						if(applicationMode & PROCESS_USER_MODE)
						{
							if(mp5Parameters->progressBar)
							{
								toolbar(step);
								step++;
							}
						}
						else
						{
							printf("TESTED %u\n",atomsCounter);
							fflush(stdout);
						}
					}	
				}
			}
		}

		atom = dictionary->atomsTable + indexAtom;
		mp5Parameters->previousAtom = atom;
		atom->feature|=GABOR_WAS_HIT;
		copyAtom(atom,bestAtom,mp5Parameters->numberOfAnalysedChannels);
       	memcpy((void *)bestModulusesTable,(void *)tmpBestModulusesTable,numberOfAnalysedChannels*sizeof(double));
		memcpy((void *)bestPhasesTable,(void *)(bestAtom->phase),numberOfAnalysedChannels*sizeof(float));

		makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);

		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
			makeAtomTable(mp5Parameters,bestAtom,channel);

			prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
			residueTable   =  *(multiChannelResidueTable + mp5Parameters->chosenChannels[channel] - 1);
			findResidue(residueTable,prevAtomTable,*(bestModulusesTable + channel),offsetExpandedDimension);
		}

		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
			residueTable = *(multiChannelResidueTable + mp5Parameters->chosenChannels[channel] - 1);
			*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,offsetExpandedDimension);
			mp5Parameters->totalResidueEnergy+= (*(mp5Parameters->residueEnergyInEachChannel + channel));
		}

		addNode(mp5Parameters->fitted,(void *)bestAtom);
	
		iterationCounter++;

		t2 = Clock();

		if(applicationMode & PROCESS_USER_MODE)
		{
			printf(" ATOM: [%3d], SEC.: %6.2lf, SIG: %6.2lf, MOD: %6.2lf, RES: %6.2lf, RES/SIG: %6.2lf",iterationCounter,(t2-t1),mp5Parameters->totalSignalEnergy,ddot(numberOfAnalysedChannels,bestModulusesTable,inc,bestModulusesTable,inc),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

			if(atom->feature & DIRACDELTA)	
				printf(" D\n");
			else if(atom->feature & GAUSSFUNCTION)
				printf(" N\n");
			else if(atom->feature & SINCOSWAVE)	
				printf(" H\n");
			else if(atom->feature & GABORWAVE)	
				printf(" G\n");
		
			printf("\n");
			fflush(stdout);
		}
    }
}

// New MMP4 algorithm
static void countMeanSignalOrResidumOverChannels(MP5Parameters *mp5Parameters, double **multiChannelSignalTable)
{
	const unsigned int offsetExpandedDimension = mp5Parameters->offsetExpandedDimension;

	unsigned int sample;
	double tmpDataValue;
	double *meanSignalTable = mp5Parameters->meanSignalTable;
	int channel;

	for(sample=0;sample<offsetExpandedDimension;sample++)
	{
		tmpDataValue = 0.0;
		
		for(channel=0;channel<mp5Parameters->numberOfAnalysedChannels;channel++)
		{
			tmpDataValue+= (*(*(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1) + sample));
		}
		*(meanSignalTable + sample) = tmpDataValue/mp5Parameters->numberOfAnalysedChannels;
	}
}

void firstAndNextIterationMMP4(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned long  int stepInToolbar, step = 1;
    unsigned short int channel = 0;
    unsigned short int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
    unsigned short int iterationCounter;
    unsigned       int atomsCounter;
    unsigned       int indexAtom = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;   
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
    const int inc = 1;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;
    const unsigned  int offsetExpandedDimension = mp5Parameters->offsetExpandedDimension;
	const unsigned  int numberOfNonFFTAtoms     = dictionary->numberOfNonFFTAtoms;

    double *signalTable  = NULL;
    double *residueTable = NULL;
    double *prevAtomTable;    
	double **multiChannelSignalTable  = mp5Parameters->multiChannelSignalTable;
    double **multiChannelResidueTable = mp5Parameters->multiChannelResidueTable = mp5Parameters->multiChannelSignalTable;

	countMeanSignalOrResidumOverChannels(mp5Parameters,multiChannelSignalTable);

    mp5Parameters->totalSignalEnergy = 0.0;
	
	for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
		signalTable = *(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);
		*(mp5Parameters->signalEnergyInEachChannel + channel) = findSignalEnergy(signalTable,offsetExpandedDimension);
		mp5Parameters->totalSignalEnergy = mp5Parameters->totalSignalEnergy + (*(mp5Parameters->signalEnergyInEachChannel + channel));
    }
	
	mp5Parameters->totalResidueEnergy = mp5Parameters->totalSignalEnergy;
	
    double t1, t2;
	double energyProgress, iterationProgress;

    double modulusesTable[numberOfAnalysedChannels];
    double bestSumOfModuluses;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
    double tmpBestModulusesTable[numberOfAnalysedChannels];
	
    Atom *atom;
    Atom *bestAtom = NULL;
	
    stepInToolbar = dictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;

    iterationCounter = 0;

    while(mp5Parameters->totalResidueEnergy>(mp5Parameters->totalSignalEnergy - mp5Parameters->totalSignalEnergy*mp5Parameters->energyPercent/100.0) && iterationCounter<mp5Parameters->maximalNumberOfIterations)
    {
		if(applicationMode & PROCESS_SERVER_MODE)
		{
			energyProgress    = (((mp5Parameters->totalSignalEnergy - mp5Parameters->totalResidueEnergy))/mp5Parameters->totalSignalEnergy)*100.0;
			iterationProgress = (iterationCounter/mp5Parameters->maximalNumberOfIterations)*100.0;
			printf("ATOM\t%3u\t%3u\t%6.2lf\t%6.2lf \n",iterationCounter,dictionary->sizeOfDictionary,energyProgress,iterationProgress);			
			fflush(stdout);
		}

		t1 = Clock();

		step = 1;
		atomsCounter = 0;
		bestSumOfModuluses = 0.0;
		mp5Parameters->totalResidueEnergy = 0.0;
		bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

		atom = dictionary->atomsTable;
		residueTable = mp5Parameters->meanSignalTable;
	
		if(mp5Parameters->FFT)
		{
			if(dictionary->diracInDictionary || dictionary->gaussInDictionary)
			{
				for(atomsCounter=0;atomsCounter<numberOfNonFFTAtoms;atomsCounter++)
				{
					if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
					{
						if(iterationCounter==0)
							normAtomTable(dictionary,mp5Parameters,atom);
						
						findAtomDataDotProduct(dictionary,mp5Parameters,atom,residueTable,0,FIRST_ITERATION);
							
						if(findUnknowPhaseDI(atom,(modulusesTable + 0),0) == ERROR)
							continue;
					}

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,1,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					if((atomsCounter+1)%stepInToolbar==0)
					{
						if(applicationMode & PROCESS_USER_MODE)
						{
							if(mp5Parameters->progressBar)
							{
								toolbar(step);
								step++;
							}
						}
						else
						{
							printf("TESTED %u\n",atomsCounter);
							fflush(stdout);
						}
					}
					atom++;
				}
			}
			
			/* from now the code is FFT for SIN/COS waves */
			if(dictionary->sinCosInDictionary)
			{
				findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,residueTable,0,FIRST_ITERATION);
				numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

				for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
				{
					if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
					{
						if(iterationCounter==0)
							normAtomTable(dictionary,mp5Parameters,atom);

						if(findUnknowPhaseDI(atom,(modulusesTable + 0),0) == ERROR)
							continue;
	
						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						if((atomsCounter+1)%stepInToolbar==0)
						{
							if(applicationMode & PROCESS_USER_MODE)
							{
								if(mp5Parameters->progressBar)
								{
									toolbar(step);
									step++;
								}
							}
							else
							{
								printf("TESTED %u\n",atomsCounter);
								fflush(stdout);
							}
						}	
					}
					atom++;
					atomsCounter++;
				}
				/* end of FFT code for SIN/COS waves */
			}
			
			for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++)
			{
				numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);
				
				for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
				{
					findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,residueTable,0,FIRST_ITERATION);

					numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

					for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
					{
						if(iterationCounter==0)
							normAtomTable(dictionary,mp5Parameters,atom);

						if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
						{								
							if(findUnknowPhaseDI(atom,(modulusesTable + 0),0) == ERROR)
								continue;
						}
	
						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,1,&bestSumOfModuluses))
							indexAtom = atomsCounter;
								
						if((atomsCounter+1)%stepInToolbar==0)
						{
							if(applicationMode & PROCESS_USER_MODE)
							{
								if(mp5Parameters->progressBar)
								{
									toolbar(step);
									step++;
								}
							}
							else
							{
								printf("TESTED %u\n",atomsCounter);
								fflush(stdout);
							}
						}	
						atom++;
						atomsCounter++;
					}
				}
			}
		}
		else
		{			
			for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
			{
				atom = dictionary->atomsTable + atomsCounter;

				if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
				{			
					if(iterationCounter==0)
						normAtomTable(dictionary,mp5Parameters,atom);

					findAtomDataDotProduct(dictionary,mp5Parameters,atom,residueTable,0,FIRST_ITERATION);
						
					if(findUnknowPhaseDI(atom,(modulusesTable + 0),0) == ERROR)
						continue;
				}

				if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,1,&bestSumOfModuluses))
					indexAtom = atomsCounter;
				
				if((atomsCounter+1)%stepInToolbar==0)
				{
					if(applicationMode & PROCESS_USER_MODE)
					{
						if(mp5Parameters->progressBar)
						{
							toolbar(step);
							step++;
						}
					}
					else
					{
						printf("TESTED %u\n",atomsCounter);
						fflush(stdout);
					}
				}	
			}
		}

		atom = dictionary->atomsTable + indexAtom;
		mp5Parameters->previousAtom = atom;
		atom->feature|=GABOR_WAS_HIT;
		
	    for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
			signalTable = *(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);
			findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
			findUnknowPhaseDI(atom,(modulusesTable + channel),channel);
		}
		
		copyAtom(atom,bestAtom,mp5Parameters->numberOfAnalysedChannels); 
       	memcpy((void *)bestModulusesTable,(void *)modulusesTable,numberOfAnalysedChannels*sizeof(double));

		makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);

		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
			makeAtomTable(mp5Parameters,bestAtom,channel);
			prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
			residueTable  = *(multiChannelResidueTable + mp5Parameters->chosenChannels[channel] - 1);
			findResidue(residueTable,prevAtomTable,*(bestModulusesTable + channel),offsetExpandedDimension);
		}
		
		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
			residueTable = *(multiChannelResidueTable + mp5Parameters->chosenChannels[channel] - 1);
			*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,offsetExpandedDimension);
			mp5Parameters->totalResidueEnergy+= (*(mp5Parameters->residueEnergyInEachChannel + channel));
		}

		countMeanSignalOrResidumOverChannels(mp5Parameters,multiChannelResidueTable);
		
		addNode(mp5Parameters->fitted,(void *)bestAtom);
	
		iterationCounter++;

		t2 = Clock();

		if(applicationMode & PROCESS_USER_MODE)
		{
			printf(" ATOM: [%3d], SEC.: %6.2lf, SIG: %6.2lf, MOD: %6.2lf, RES: %6.2lf, RES/SIG: %6.2lf",iterationCounter,(t2-t1),mp5Parameters->totalSignalEnergy,ddot(numberOfAnalysedChannels,bestModulusesTable,inc,bestModulusesTable,inc),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

			if(atom->feature & DIRACDELTA)	
				printf(" D\n");
			else if(atom->feature & GAUSSFUNCTION)
				printf(" N\n");
			else if(atom->feature & SINCOSWAVE)	
				printf(" H\n");
			else if(atom->feature & GABORWAVE)	
				printf(" G\n");
		
			printf("\n");
			fflush(stdout);
		}
    }
}
