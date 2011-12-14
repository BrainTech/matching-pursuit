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


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include"atom.h"
#include"dic.h"
#include"io_mp5.h"
#include"matrix.h"
#include"mp5.h"
#include"queue.h"
#include"tools.h"
#include"types.h"
#include"vector.h"

extern unsigned char applicationMode;

/* MMP - MULTICHANNEL MATCHING PURSUIT */


static BOOLEAN getBestModulusesTable(unsigned short int MMPType,
									 Atom *atom,
									 double *modulusesTable,
									 double *tmpBestModulusesTable,
									 unsigned int numberOfAnalysedChannels,
									 double *bestSumOfModuluses)
{
	int inc = 1;
	double sumOfModuluses = 0.0; // sum of moduluses or sum of sqr moduluses

	if((MMPType & MMP1) || (MMPType & MMP11) || (MMPType & MMP12) || (MMPType & MMP21))
	{
		if(findUnknowPhaseAM(atom,modulusesTable,numberOfAnalysedChannels)==ERROR)
			return FALSE;

		sumOfModuluses = dasum(numberOfAnalysedChannels,modulusesTable,inc);

	}
	else if((MMPType & MMP3) || (MMPType & MMP33) || (MMPType & MMP23) || (MMPType & MMP32))
		sumOfModuluses = ddot(numberOfAnalysedChannels,modulusesTable,inc,modulusesTable,inc);

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
    unsigned       int channel;
    unsigned       int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
    unsigned       int atomsCounter;
    unsigned       int tmpAtomsCounter;
    unsigned       int indexAtom = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
	const int inc = 1;
	Progress       progress;
	progress.applicationMode = applicationMode;
    progress.stepInToolbar   = dictionary->initialNumberOfAtoms/NUMBER_OF_STEPS_IN_TOOLBAR;
    progress.step            = 1;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    double **multiChannelSignalTable  = mp5Parameters->multiChannelSignalTable;
    double **multiChannelResidueTable = mp5Parameters->multiChannelResidueTable = mp5Parameters->multiChannelSignalTable;

    double *prevAtomTable;
    double *signalTable, *residueTable;

    const unsigned int epochExpandedSize   = mp5Parameters->epochExpandedSize;

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
		signalTable = *(multiChannelSignalTable + channel);
		*(mp5Parameters->signalEnergyInEachChannel + channel) = findSignalEnergy(signalTable,epochExpandedSize);
		mp5Parameters->totalSignalEnergy = mp5Parameters->totalSignalEnergy + (*(mp5Parameters->signalEnergyInEachChannel + channel));
    }

    t1 = Clock();

	if(applicationMode & PROCESS_SERVER_MODE)
	{
		printf("ATOM\t%3u\t%3u\t%6.2f\t%6.2f\n",0,dictionary->finalNumberOfAtoms,0.0,0.0);
		fflush(stdout);
	}

    bestSumOfModuluses = 0.0;
    mp5Parameters->totalResidueEnergy = 0.0;

	atomsCounter    = 0;
	tmpAtomsCounter = 0;

	if(mp5Parameters->FFT)
	{
		if(dictionary->diracInDictionary)
		{
			atom = dictionary->diracAtomsTable;

			for(atomsCounter=0;atomsCounter<dictionary->numberOfFinalDiracFunctions;atomsCounter++)
			{
				if((atom->feature & INCORRECTGABOR)==0)
				{
					normAtomTable(dictionary,mp5Parameters,atom);

					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{
						signalTable = *(multiChannelSignalTable + channel);

						findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);

						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
						{
							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
							break;						
						}							
					}

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
				}
				atom++;
			}
		}

		if(dictionary->gaussInDictionary)
		{
			atom = dictionary->gaussAtomsTable;

			for(tmpAtomsCounter=0;tmpAtomsCounter<dictionary->numberOfFinalGaussFunctions;tmpAtomsCounter++)
			{
				if((atom->feature & INCORRECTGABOR)==0)
				{
					normAtomTable(dictionary,mp5Parameters,atom);

					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{
						signalTable = *(multiChannelSignalTable + channel);

						findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);

						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
						{
							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
							break;						
						}							
					}

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);

				}
				atom++;
				atomsCounter++;
			}
		}

		/* from now the code is FFT for SIN/COS waves */
		if(dictionary->sinCosInDictionary)
		{
			atom = dictionary->sinCosAtomsTable;

			for(channel=0;channel<numberOfAnalysedChannels;channel++)
			{
				signalTable = *(multiChannelSignalTable + channel);
				findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
			}

			numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

			for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
			{
				if((atom->feature & INCORRECTGABOR)==0)
				{
					normAtomTable(dictionary,mp5Parameters,atom);

					for(channel=0;channel<numberOfAnalysedChannels;channel++)
						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
						{
							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
							break;						
						}							

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
	
				}
				atom++;
				atomsCounter++;
			}
		}
		/* end of FFT code for SIN/COS waves */

		if(dictionary->gaborInDictionary)
		{
			atom = dictionary->gaborAtomsTable;

			for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
			{
				numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);

				for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
				{
					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{
						signalTable = *(multiChannelSignalTable + channel);
						findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
					}

					numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

					for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
					{
						if((atom->feature & INCORRECTGABOR)==0)
						{
							normAtomTable(dictionary,mp5Parameters,atom);

							for(channel=0;channel<numberOfAnalysedChannels;channel++)
							{
								if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
								{
									printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
									break;						
								}							
							}
							if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
								indexAtom = atomsCounter;

							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
						}
						atom++;
						atomsCounter++;
					}
				}
			}
		}
	}
	else
	{
		for(atomsCounter=0;atomsCounter<dictionary->initialNumberOfAtoms;atomsCounter++)
		{
			atom = getAtom(dictionary,atomsCounter);

			if((atom->feature & INCORRECTGABOR) || ((atom->feature & STOCHASTIC_ATOM)==0))
			{
				printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
				continue;
			}

			normAtomTable(dictionary,mp5Parameters,atom);

			for(channel=0;channel<numberOfAnalysedChannels;channel++)
			{
				signalTable = *(multiChannelSignalTable + channel);

				findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);

				if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
				{
					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					break;						
				}
			}

			if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
				indexAtom = atomsCounter;

			printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
		}
	}
	
	atom = getAtom(dictionary,indexAtom);

    mp5Parameters->previousAtom = atom;
    atom->feature|=ATOM_WAS_SELECTED;
    copyAtom(atom,bestAtom,mp5Parameters->numberOfAnalysedChannels);
    memcpy((void *)bestModulusesTable,(void *)tmpBestModulusesTable,numberOfAnalysedChannels*sizeof(double));
    memcpy((void *)bestPhasesTable,(void *)(bestAtom->phase),numberOfAnalysedChannels*sizeof(float));

    makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
        makeAtomTable(mp5Parameters,bestAtom,channel);
		residueTable  = *(multiChannelResidueTable + channel);
        prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
        findResidue(residueTable,prevAtomTable,*(bestModulusesTable + channel),epochExpandedSize);        
    }

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
		residueTable = *(multiChannelResidueTable + channel);
		*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,epochExpandedSize);
		mp5Parameters->totalResidueEnergy+= (*(mp5Parameters->residueEnergyInEachChannel + channel));
    }

    addNode(mp5Parameters->fitted,(void *)bestAtom);

    t2 = Clock();

	if(applicationMode & PROCESS_USER_MODE)
	{
		printf(" ATOM: [%3d], SEC.: %6.2f, SIG: %6.2f, MOD: %6.2f, RES: %6.2f, RES/SIG: %6.2f",1,(t2-t1),mp5Parameters->totalSignalEnergy,ddot(numberOfAnalysedChannels,bestModulusesTable,inc,bestModulusesTable,inc),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

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
    unsigned       int channel = 0;
    unsigned       int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
	unsigned       int lastChannel = numberOfAnalysedChannels;
    unsigned short int iterationCounter;
    unsigned       int atomsCounter;
    unsigned       int tmpAtomsCounter;
    unsigned       int indexAtom = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
    const int inc = 1;
	Progress       progress;
	progress.applicationMode = applicationMode;
    progress.stepInToolbar   = dictionary->initialNumberOfAtoms/NUMBER_OF_STEPS_IN_TOOLBAR;
    progress.step            = 1;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    double **multiChannelResidueTable = mp5Parameters->multiChannelResidueTable;

    double *prevAtomTable;
    double *residueTable;

    const unsigned       int epochExpandedSize   = mp5Parameters->epochExpandedSize;

    double t1, t2;
	double energyProgress, iterationProgress;

    double modulusesTable[numberOfAnalysedChannels];
    double bestSumOfModuluses;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
    float  *bestPhasesTable    = mp5Parameters->bestPhasesTable;
    double tmpBestModulusesTable[numberOfAnalysedChannels];
    double energyStopCondition = mp5Parameters->totalSignalEnergy *(1 - mp5Parameters->energyPercent/100.0);
	double RC = 0.0;
	double RS = 0.0;

    Atom *atom;
    Atom *bestAtom = NULL;

    iterationCounter = 1;

    while((mp5Parameters->totalResidueEnergy > energyStopCondition) && (iterationCounter<mp5Parameters->maximalNumberOfIterations))
    {
		if(applicationMode & PROCESS_SERVER_MODE)
		{
			energyProgress    = (((mp5Parameters->totalSignalEnergy - mp5Parameters->totalResidueEnergy))/mp5Parameters->totalSignalEnergy)*100.0;
			iterationProgress = (iterationCounter/mp5Parameters->maximalNumberOfIterations)*100.0;
			printf("ATOM\t%3u\t%3u\t%6.2f\t%6.2f \n",iterationCounter,dictionary->finalNumberOfAtoms,energyProgress,iterationProgress);
			fflush(stdout);
		}

		t1 = Clock();

	    progress.step   = 1;
		atomsCounter    = 0;
		tmpAtomsCounter = 0;
		bestSumOfModuluses = 0.0;
		mp5Parameters->totalResidueEnergy = 0.0;
		bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

		if(mp5Parameters->FFT)
		{
			if(dictionary->diracInDictionary)
			{
				atom = dictionary->diracAtomsTable;

				for(atomsCounter=0;atomsCounter<dictionary->numberOfFinalDiracFunctions;atomsCounter++)
				{
					if(((atom->feature & ATOM_WAS_SELECTED)==0) && ((atom->feature & INCORRECTGABOR)==0))
					{
						if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							{
								printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
								break;
							}
						}

						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					}

					atom++;
				}
			}

			if(dictionary->gaussInDictionary)
			{
				atom = dictionary->gaussAtomsTable;

				for(tmpAtomsCounter=0;tmpAtomsCounter<dictionary->numberOfFinalGaussFunctions;tmpAtomsCounter++)
				{
					if(((atom->feature & ATOM_WAS_SELECTED)==0) && ((atom->feature & INCORRECTGABOR)==0))
					{
						if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							{
								printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
								break;
							}
						}

						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);

					}
					atom++;
					atomsCounter++;
				}
			}

			/* from now the code is FFT for SIN/COS waves */
			if(dictionary->sinCosInDictionary)
			{
				atom = dictionary->sinCosAtomsTable;

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
					if(((atom->feature & ATOM_WAS_SELECTED)==0) && ((atom->feature & INCORRECTGABOR)==0))
					{
						if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
						{
							RS = *(atom->RS + lastChannel);
							RC = *(atom->RC + lastChannel);
						}

						for(channel=0;channel<numberOfAnalysedChannels;channel++)
						{
							if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							{
								printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
								break;
							}
						}

						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
						
					}

					atom++;
					atomsCounter++;
				}
			}
			/* end of FFT code for SIN/COS waves */

			if(dictionary->gaborInDictionary)
			{
				atom = dictionary->gaborAtomsTable;

				for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
				{
					numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);

					for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
					{
						if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							if(((atom->feature & ATOM_WAS_SELECTED)==0) && ((atom->feature & INCORRECTGABOR)==0))
							{
								if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
								{
									RS = *(atom->RS + lastChannel);
									RC = *(atom->RC + lastChannel);
								}

								for(channel=0;channel<numberOfAnalysedChannels;channel++)
								{
									if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
									{
										printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
										break;
									}
								}

								if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
									indexAtom = atomsCounter;

								printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
							}
							atom++;
							atomsCounter++;
						}
					}
				}
			}
		}
		else
		{
			for(atomsCounter=0;atomsCounter<dictionary->initialNumberOfAtoms;atomsCounter++)
			{
	 			atom = getAtom(dictionary,atomsCounter);
	 			
				if((atom->feature & ATOM_WAS_SELECTED) || (atom->feature & INCORRECTGABOR) || ((atom->feature & STOCHASTIC_ATOM)==0))
				{
					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					continue;
				}
	 			
				if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
					if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
					{
						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
						break;
					}
				}

				if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
					indexAtom = atomsCounter;

				printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
			}
		}

		atom = getAtom(dictionary,indexAtom);		

		mp5Parameters->previousAtom = atom;
		atom->feature|=ATOM_WAS_SELECTED;
		copyAtom(atom,bestAtom,mp5Parameters->numberOfAnalysedChannels);
	   	memcpy((void *)bestModulusesTable,(void *)tmpBestModulusesTable,numberOfAnalysedChannels*sizeof(double));
		memcpy((void *)bestPhasesTable,(void *)(bestAtom->phase),numberOfAnalysedChannels*sizeof(float));

		makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);

		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
			makeAtomTable(mp5Parameters,bestAtom,channel);

			prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
			residueTable  = *(multiChannelResidueTable + channel);
			findResidue(residueTable,prevAtomTable,*(bestModulusesTable + channel),epochExpandedSize);
		}

		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
			residueTable = *(multiChannelResidueTable + channel);
			*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,epochExpandedSize);
			mp5Parameters->totalResidueEnergy+= (*(mp5Parameters->residueEnergyInEachChannel + channel));
		}

		addNode(mp5Parameters->fitted,(void *)bestAtom);

		iterationCounter++;

		t2 = Clock();

		if(applicationMode & PROCESS_USER_MODE)
		{
			printf(" ATOM: [%3d], SEC.: %6.2f, SIG: %6.2f, MOD: %6.2f, RES: %6.2f, RES/SIG: %6.2f",iterationCounter,(t2-t1),mp5Parameters->totalSignalEnergy,ddot(numberOfAnalysedChannels,bestModulusesTable,inc,bestModulusesTable,inc),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

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

void firstIterationMultiChannelMultiTrial(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
	unsigned short int MPType = mp5Parameters->MPType;
    unsigned       int channel;
    unsigned       int numberOfAnalysedChannels      = mp5Parameters->numberOfAnalysedChannels;
    unsigned       int numberOfReadChannelsAndEpochs = mp5Parameters->numberOfReadChannelsAndEpochs;
    unsigned       int atomsCounter;
    unsigned       int tmpAtomsCounter;
    unsigned       int indexAtom = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
	const int inc = 1;
	Progress       progress;
	progress.applicationMode = applicationMode;
    progress.stepInToolbar   = dictionary->initialNumberOfAtoms/NUMBER_OF_STEPS_IN_TOOLBAR;
    progress.step            = 1;

//	float tmpAmplitude = 0;
//	float tmpModulus   = 0;
//	double totalMMP2Modulus = 0.0;
//	double totalMMP2Residue = 0.0;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    double **multiChannelSignalTable  = mp5Parameters->multiChannelSignalTable;
	double **meanSignalTable          = mp5Parameters->meanSignalTable;
	double **meanResidueTable         = mp5Parameters->meanResidueTable = mp5Parameters->meanSignalTable;

    double *prevAtomTable;
    double *signalTable, *residueTable;

    const unsigned int epochExpandedSize   = mp5Parameters->epochExpandedSize;

    double t1, t2;

    double modulusesTable[numberOfAnalysedChannels];
    double bestSumOfModuluses;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
    float  *bestPhasesTable    = mp5Parameters->bestPhasesTable;

    double tmpBestModulusesTable[numberOfAnalysedChannels];

    Atom *atom;
    Atom *bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

    mp5Parameters->totalSignalEnergy = 0.0;

    for(channel=0;channel<numberOfReadChannelsAndEpochs;channel++)
    {
		signalTable = *(multiChannelSignalTable + channel);
		*(mp5Parameters->signalEnergyInEachChannel + channel) = findSignalEnergy(signalTable,epochExpandedSize);
		mp5Parameters->totalSignalEnergy = mp5Parameters->totalSignalEnergy + (*(mp5Parameters->signalEnergyInEachChannel + channel));
    }

	if((MPType & MMP12) || (MPType & MMP32))
		countMeanSignalOrResidumOverEpochs(mp5Parameters,multiChannelSignalTable);
	else // MMP21 MMP23
		countMeanSignalOrResidumOverChannels(mp5Parameters,multiChannelSignalTable);

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
		signalTable = *(meanSignalTable + channel);
		*(mp5Parameters->meanSignalEnergyInEachChannel + channel) = findSignalEnergy(signalTable,epochExpandedSize);
		mp5Parameters->oneChannelSignalEnergy = mp5Parameters->oneChannelSignalEnergy + (*(mp5Parameters->meanSignalEnergyInEachChannel + channel));
    }

    t1 = Clock();

	if(applicationMode & PROCESS_SERVER_MODE)
	{
		printf("ATOM\t%3u\t%3u\t%6.2f\t%6.2f\n",0,dictionary->finalNumberOfAtoms,0.0,0.0);
		fflush(stdout);
	}

    bestSumOfModuluses = 0.0;
    mp5Parameters->totalResidueEnergy = 0.0;

	atomsCounter    = 0;
	tmpAtomsCounter = 0;

	if(mp5Parameters->FFT)
	{
		if(dictionary->diracInDictionary)
		{
			atom = dictionary->diracAtomsTable;

			for(atomsCounter=0;atomsCounter<dictionary->numberOfFinalDiracFunctions;atomsCounter++)
			{
				if((atom->feature & INCORRECTGABOR)==0)
				{
					normAtomTable(dictionary,mp5Parameters,atom);

					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{
						signalTable = *(meanSignalTable + channel);

						findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);

						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
						{
							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
							break;						
						}							
					}

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					
				}
				atom++;
			}
		}

		if(dictionary->gaussInDictionary)
		{
			atom = dictionary->gaussAtomsTable;

			for(tmpAtomsCounter=0;tmpAtomsCounter<dictionary->numberOfFinalGaussFunctions;tmpAtomsCounter++)
			{
				if((atom->feature & INCORRECTGABOR)==0)
				{
					normAtomTable(dictionary,mp5Parameters,atom);

					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{
						signalTable = *(meanSignalTable + channel);

						findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);

						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
						{
							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
							break;						
						}							
					}

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
		
				}
				atom++;
				atomsCounter++;
			}
		}

		/* from now the code is FFT for SIN/COS waves */
		if(dictionary->sinCosInDictionary)
		{
			atom = dictionary->sinCosAtomsTable;

			for(channel=0;channel<numberOfAnalysedChannels;channel++)
			{
				signalTable = *(meanSignalTable + channel);
				findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
			}

			numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

			for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
			{
				if((atom->feature & INCORRECTGABOR)==0)
				{
					normAtomTable(dictionary,mp5Parameters,atom);

					for(channel=0;channel<numberOfAnalysedChannels;channel++)
						if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
						{
							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
							break;						
						}							

					if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
						indexAtom = atomsCounter;

					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);

				}
				atom++;
				atomsCounter++;
			}
		}
		/* end of FFT code for SIN/COS waves */

		if(dictionary->gaborInDictionary)
		{
			atom = dictionary->gaborAtomsTable;

			for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
			{
				numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);

				for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
				{
					for(channel=0;channel<numberOfAnalysedChannels;channel++)
					{
						signalTable = *(meanSignalTable + channel);
						findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);
					}

					numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

					for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
					{
						if((atom->feature & INCORRECTGABOR)==0)
						{
							normAtomTable(dictionary,mp5Parameters,atom);

							for(channel=0;channel<numberOfAnalysedChannels;channel++)
							{
								if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
								{
									printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
									break;						
								}							
							}
							if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
								indexAtom = atomsCounter;

							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					
						}
						atom++;
						atomsCounter++;
					}
				}
			}
		}
	}
	else
	{
		for(atomsCounter=0;atomsCounter<dictionary->initialNumberOfAtoms;atomsCounter++)
		{
			atom = getAtom(dictionary,atomsCounter);
			
			if((atom->feature & INCORRECTGABOR) || ((atom->feature & STOCHASTIC_ATOM)==0))
			{
				printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
				continue;
			}							
			
			normAtomTable(dictionary,mp5Parameters,atom);

			for(channel=0;channel<numberOfAnalysedChannels;channel++)
			{
				signalTable = *(meanSignalTable + channel);

				findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,channel,FIRST_ITERATION);

				if(findUnknowPhaseDI(atom,(modulusesTable + channel),channel) == ERROR)
				{
					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					break;						
				}							
			}

			if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
				indexAtom = atomsCounter;

			printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);								

		}
    }

	atom = getAtom(dictionary,indexAtom);
    mp5Parameters->previousAtom = atom;
    atom->feature|=ATOM_WAS_SELECTED;
    copyAtom(atom,bestAtom,mp5Parameters->numberOfAnalysedChannels);
    memcpy((void *)bestModulusesTable,(void *)tmpBestModulusesTable,numberOfAnalysedChannels*sizeof(double));
    memcpy((void *)bestPhasesTable,(void *)(bestAtom->phase),numberOfAnalysedChannels*sizeof(float));

    makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
        makeAtomTable(mp5Parameters,bestAtom,channel);
		residueTable  = *(meanResidueTable + channel);
        prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
        findResidue(residueTable,prevAtomTable,*(bestModulusesTable + channel),epochExpandedSize);
    }

    for(channel=0;channel<numberOfAnalysedChannels;channel++)
    {
		residueTable = *(meanResidueTable + channel);
		*(mp5Parameters->meanResidueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,epochExpandedSize);
		mp5Parameters->oneChannelResidueEnergy+=(*(mp5Parameters->meanResidueEnergyInEachChannel + channel));
	}

    addNode(mp5Parameters->fitted,(void *)bestAtom);

	/*
		for(channel=0;channel<numberOfReadChannelsAndEpochs;channel++)
		{
			returnAmplitudeAndModulusForMMP2DI(mp5Parameters,dictionary,bestAtom,&tmpAmplitude,&tmpModulus,channel);
			totalMMP2Modulus = totalMMP2Modulus + tmpModulus*tmpModulus;
			signalTable = *(mp5Parameters->multiChannelSignalTable + channel);
			*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(signalTable,epochExpandedSize);
			totalMMP2Residue = totalMMP2Residue + (*(mp5Parameters->residueEnergyInEachChannel + channel));
		}
	*/

    t2 = Clock();

	if(applicationMode & PROCESS_USER_MODE)
	{
		printf(" ATOM: [%3d], SEC.: %6.2f, SIG: %6.2f, MOD: %6.2f, RES: %6.2f, RES/SIG: %6.2f",1,(t2-t1),mp5Parameters->oneChannelSignalEnergy,ddot(numberOfAnalysedChannels,bestModulusesTable,inc,bestModulusesTable,inc),mp5Parameters->oneChannelResidueEnergy,(mp5Parameters->oneChannelResidueEnergy/mp5Parameters->oneChannelSignalEnergy)*100.0);
//		printf(" SIG: %-6.2f, MOD: %-6.2f, RES: %-6.2f, RES/SIG: %-6.2f %% ",mp5Parameters->totalSignalEnergy,totalMMP2Modulus,totalMMP2Residue,(totalMMP2Residue/mp5Parameters->totalSignalEnergy)*100.0);

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

void nextIterationMultiChannelMultiTrial(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned       int channel = 0;
    unsigned       int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
//    unsigned       int numberOfReadChannelsAndEpochs = mp5Parameters->numberOfReadChannelsAndEpochs;
	unsigned       int lastChannel = numberOfAnalysedChannels;
    unsigned short int iterationCounter;
    unsigned       int atomsCounter;
    unsigned       int tmpAtomsCounter;
    unsigned       int indexAtom = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
    const int inc = 1;
	Progress       progress;
	progress.applicationMode = applicationMode;
    progress.stepInToolbar   = dictionary->initialNumberOfAtoms/NUMBER_OF_STEPS_IN_TOOLBAR;
    progress.step            = 1;   

//	float tmpAmplitude = 0;
//	float tmpModulus   = 0;
//	double totalMMP2Modulus = 0.0;
//	double totalMMP2Residue = 0.0;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

	double **meanResidueTable = mp5Parameters->meanResidueTable = mp5Parameters->meanSignalTable;

    double *prevAtomTable;
	double /* *signalTable,*/ *residueTable;

    const unsigned       int epochExpandedSize   = mp5Parameters->epochExpandedSize;

    double t1, t2;
	double energyProgress, iterationProgress;

    double modulusesTable[numberOfAnalysedChannels];
    double bestSumOfModuluses;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
    float  *bestPhasesTable    = mp5Parameters->bestPhasesTable;
    double tmpBestModulusesTable[numberOfAnalysedChannels];
    double energyStopCondition = mp5Parameters->totalSignalEnergy *(1 - mp5Parameters->energyPercent/100.0);    
	double RC = 0.0;
	double RS = 0.0;

    Atom *atom;
    Atom *bestAtom = NULL;

    iterationCounter = 1;

    while((mp5Parameters->oneChannelResidueEnergy>energyStopCondition) && (iterationCounter<mp5Parameters->maximalNumberOfIterations))
    {
		if(applicationMode & PROCESS_SERVER_MODE)
		{
			energyProgress    = (((mp5Parameters->totalSignalEnergy - mp5Parameters->totalResidueEnergy))/mp5Parameters->totalSignalEnergy)*100.0;
			iterationProgress = (iterationCounter/mp5Parameters->maximalNumberOfIterations)*100.0;
			printf("ATOM\t%3u\t%3u\t%6.2f\t%6.2f \n",iterationCounter,dictionary->finalNumberOfAtoms,energyProgress,iterationProgress);
			fflush(stdout);
		}

		t1 = Clock();

	    progress.step   = 1;   
		atomsCounter    = 0;
		tmpAtomsCounter = 0;
		bestSumOfModuluses = 0.0;
		mp5Parameters->totalResidueEnergy = 0.0;
		mp5Parameters->oneChannelResidueEnergy = 0.0;
		bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

//		tmpAmplitude = 0.0;
//		tmpModulus   = 0.0;
//		totalMMP2Modulus = 0.0;
//		totalMMP2Residue = 0.0;

		if(mp5Parameters->FFT)
		{
			if(dictionary->diracInDictionary)
			{
				atom = dictionary->diracAtomsTable;

				for(atomsCounter=0;atomsCounter<dictionary->numberOfFinalDiracFunctions;atomsCounter++)
				{
					if(((atom->feature & ATOM_WAS_SELECTED)==0) && ((atom->feature & INCORRECTGABOR)==0))
					{
						if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							{
								printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
								break;						
							}							
						}

						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					}
					atom++;
				}
			}

			if(dictionary->gaussInDictionary)
			{
				atom = dictionary->gaussAtomsTable;

				for(tmpAtomsCounter=0;tmpAtomsCounter<dictionary->numberOfFinalGaussFunctions;tmpAtomsCounter++)
				{
					if(!(atom->feature & ATOM_WAS_SELECTED) && !(atom->feature & INCORRECTGABOR))
					{
						if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							{
								printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
								break;						
							}							
						}

						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
		
					}
					atom++;
					atomsCounter++;
				}
			}

			/* from now the code is FFT for SIN/COS waves */
			if(dictionary->sinCosInDictionary)
			{
				atom = dictionary->sinCosAtomsTable;

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
					if(((atom->feature & ATOM_WAS_SELECTED)==0) && ((atom->feature & INCORRECTGABOR)==0))
					{
						if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
						{
							RS = *(atom->RS + lastChannel);
							RC = *(atom->RC + lastChannel);
						}

						for(channel=0;channel<numberOfAnalysedChannels;channel++)
						{
							if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							{
								printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
								break;						
							}							
						}

						if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
							indexAtom = atomsCounter;

						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					}
					atom++;
					atomsCounter++;
				}
			}
			/* end of FFT code for SIN/COS waves */

			if(dictionary->gaborInDictionary)
			{
				atom = dictionary->gaborAtomsTable;

				for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
				{
					numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);

					for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
					{
						if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
							if(!(atom->feature & ATOM_WAS_SELECTED) && !(atom->feature & INCORRECTGABOR))
							{
								if(mp5Parameters->MPType & MMP1)
								{
									RS = *(atom->RS + lastChannel);
									RC = *(atom->RC + lastChannel);
								}

								for(channel=0;channel<numberOfAnalysedChannels;channel++)
								{
									if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
									{
										printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
										break;						
									}							
								}

								if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
									indexAtom = atomsCounter;

								printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);

							}
							atom++;
							atomsCounter++;
						}
					}
				}
			}
		}
		else
		{
			for(atomsCounter=0;atomsCounter<dictionary->initialNumberOfAtoms;atomsCounter++)
			{
				atom = getAtom(dictionary,atomsCounter);

				if((atom->feature & ATOM_WAS_SELECTED) || (atom->feature & INCORRECTGABOR) || ((atom->feature & STOCHASTIC_ATOM)==0))
				{
					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					continue;
				}							

				if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
					if((mp5Parameters->MPType & MMP1) || (mp5Parameters->MPType & MMP11))
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
					{
						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
						break;
					}							
						
				}

				if(getBestModulusesTable(mp5Parameters->MPType,atom,modulusesTable,tmpBestModulusesTable,numberOfAnalysedChannels,&bestSumOfModuluses))
					indexAtom = atomsCounter;

				printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
				
			}
		}
		
		atom = getAtom(dictionary,indexAtom);
		mp5Parameters->previousAtom = atom;
		atom->feature|=ATOM_WAS_SELECTED;
		copyAtom(atom,bestAtom,mp5Parameters->numberOfAnalysedChannels);
	   	memcpy((void *)bestModulusesTable,(void *)tmpBestModulusesTable,numberOfAnalysedChannels*sizeof(double));
		memcpy((void *)bestPhasesTable,(void *)(bestAtom->phase),numberOfAnalysedChannels*sizeof(float));
	
		makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);
	
		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{

			makeAtomTable(mp5Parameters,bestAtom,channel);
	
			prevAtomTable = *(mp5Parameters->prevAtomTable + channel);
			residueTable  = *(meanResidueTable + channel);
			findResidue(residueTable,prevAtomTable,*(bestModulusesTable + channel),epochExpandedSize);
		}
	
		for(channel=0;channel<numberOfAnalysedChannels;channel++)
		{
			residueTable = *(meanResidueTable + channel);
			*(mp5Parameters->meanResidueEnergyInEachChannel + channel) = findSignalEnergy(residueTable,epochExpandedSize);
			mp5Parameters->oneChannelResidueEnergy+= (*(mp5Parameters->meanResidueEnergyInEachChannel + channel));
		}
	
		addNode(mp5Parameters->fitted,(void *)bestAtom);
	
		iterationCounter++;
	
		t2 = Clock();
	
		if(applicationMode & PROCESS_USER_MODE)
		{
			printf(" ATOM: [%3d], SEC.: %6.2f, SIG: %6.2f, MOD: %6.2f, RES: %6.2f, RES/SIG: %6.2f",iterationCounter,(t2-t1),mp5Parameters->oneChannelSignalEnergy,ddot(numberOfAnalysedChannels,bestModulusesTable,inc,bestModulusesTable,inc),mp5Parameters->oneChannelResidueEnergy,(mp5Parameters->oneChannelResidueEnergy/mp5Parameters->oneChannelSignalEnergy)*100.0);
	//			printf(" SIG: %-6.2f, MOD: %-6.2f, RES: %-6.2f, RES/SIG: %-6.2f %% ",mp5Parameters->totalSignalEnergy,totalMMP2Modulus,totalMMP2Residue,(totalMMP2Residue/mp5Parameters->totalSignalEnergy)*100.0);
	
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
