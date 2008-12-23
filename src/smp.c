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

#include<math.h>
#include<stdlib.h>
#include"def.h"
#include"atom.h"
#include"io_mp5.h"
#include"mp5.h"
#include"queue.h"
#include"smp.h"
#include"tools.h"
#include"types.h"

/* SMP - SINGLE CHANNEL MATCHING PURSUIT */

extern unsigned char applicationMode;

static unsigned int stepInToolbar = 0;
static unsigned int step = 1;
static double   bestModulus = 0.0;
static double   absOfBestModulus = 0.0;

static BOOLEAN bestModulusUpdate(double modulus)
{
	double absOfModulus = fabs(modulus);

	if(absOfModulus>absOfBestModulus)
	{
		bestModulus      = modulus;
		absOfBestModulus = absOfModulus;
		return TRUE;
	}
	else
		return FALSE;
}

static void printProgress(unsigned int atomsCounter, unsigned char progressBar)
{		
	if((atomsCounter+1)%stepInToolbar==0)
	{
		if(applicationMode & PROCESS_USER_MODE)
		{
			if(progressBar)
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

void firstIterationSMPMMP2(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned int atomsCounter = 0;
    unsigned int bestAtomIndex = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;   
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;

	unsigned short int channel;
	
	float tmpAmplitude = 0;
	float tmpModulus   = 0;
	double totalMMP2Modulus = 0.0;
	double totalMMP2Residue = 0.0;
	
	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;
	
    double *signalTable    = mp5Parameters->singleChannelSignalTable;
    double *residueTable   = mp5Parameters->singleChannelResidueTable = mp5Parameters->singleChannelSignalTable;

    unsigned       int offsetExpandedDimension = mp5Parameters->offsetExpandedDimension;
	unsigned       int numberOfNonFFTAtoms     = dictionary->numberOfNonFFTAtoms;
	/* Atoms assocciated with FFT are sin/cos waves and Gabor functions. The dot product of Dirac Deltas and Gauss functions with signal can not be astimated by means of FFT algorithm. */
	
    double t1, t2;

	double modulus;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
	bestModulus        = 0.0;  
	absOfBestModulus   = 0.0;
	
    Atom *atom;
			
    Atom *bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

	if(mp5Parameters->MPType & SMP)
	{
		mp5Parameters->oneChannelSignalEnergy = findSignalEnergy(signalTable,offsetExpandedDimension);
	}
	else // MMP2
	{
		mp5Parameters->totalSignalEnergy = 0.0;
		
		double *tmpSignalTable = NULL;
		mp5Parameters->oneChannelSignalEnergy = findSignalEnergy(signalTable,offsetExpandedDimension);
		
		for(channel=0;channel<mp5Parameters->numberOfAnalysedChannels;channel++)
		{
			tmpSignalTable = *(mp5Parameters->multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);
			*(mp5Parameters->signalEnergyInEachChannel + channel) = findSignalEnergy(tmpSignalTable,offsetExpandedDimension);
			mp5Parameters->totalSignalEnergy = mp5Parameters->totalSignalEnergy + (*(mp5Parameters->signalEnergyInEachChannel + channel));
		}
	}
	
    t1 = Clock();

	if(applicationMode & PROCESS_SERVER_MODE)
	{
		printf("ATOM\t%3u\t%3u\t%6.2lf\t%6.2lf\n",0,dictionary->sizeOfDictionary,0.0,0.0);
		fflush(stdout);
	}

    stepInToolbar = dictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;
	step = 1;
	
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
					findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);

					if(findUnknowPhaseDI(atom,&modulus,0))
					{	
						if(bestModulusUpdate(modulus))
							bestAtomIndex = atomsCounter;
					}
				}
				atom++;
				printProgress(atomsCounter,mp5Parameters->progressBar);				
			}
		}

		if(dictionary->sinCosInDictionary)
		{
			/* from now the code is FFT for SIN/COS waves */
			findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);
			numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

			for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
			{
				if(!(atom->feature & INCORRECTGABOR))
				{
					normAtomTable(dictionary,mp5Parameters,atom);
					if(findUnknowPhaseDI(atom,&modulus,0))
					{	
						if(bestModulusUpdate(modulus))
							bestAtomIndex = atomsCounter;
					}
				}
				atom++;
				atomsCounter++;
				printProgress(atomsCounter,mp5Parameters->progressBar);	
			}
			/* end of FFT code for SIN/COS waves */
		}

		for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++)
		{
			numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);
		
			for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
			{
				findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);
				numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

				for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
				{
					if(!(atom->feature & INCORRECTGABOR))
					{
						normAtomTable(dictionary,mp5Parameters,atom);

						if(findUnknowPhaseDI(atom,&modulus,0))
						{							
							if(bestModulusUpdate(modulus))
								bestAtomIndex = atomsCounter;		
						}
					}
					atom++;
					atomsCounter++;
					printProgress(atomsCounter,mp5Parameters->progressBar);
				}
			}
		}
	}
	else
	{
		for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
		{
			if(!(atom->feature & INCORRECTGABOR))
			{
				normAtomTable(dictionary,mp5Parameters,atom);

				findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);
				findUnknowPhaseDI(atom,&modulus,0);

				if(findUnknowPhaseDI(atom,&modulus,0))
				{	
					if(bestModulusUpdate(modulus))
						bestAtomIndex = atomsCounter;
				}
			}
			atom++;
			printProgress(atomsCounter,mp5Parameters->progressBar);
		}
	}

	atom = dictionary->atomsTable + bestAtomIndex;
    atom->feature|=GABOR_WAS_HIT;
    mp5Parameters->previousAtom = atom;
	*bestModulusesTable = bestModulus;
	
    copyAtom(atom,bestAtom,mp5Parameters->numberOfAllocatedChannels);

    makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);
    makeAtomTable(mp5Parameters,bestAtom,0);
	
    findResidue(residueTable,*(mp5Parameters->prevAtomTable),*bestModulusesTable,offsetExpandedDimension);
    mp5Parameters->oneChannelResidueEnergy = findSignalEnergy(signalTable,offsetExpandedDimension);
		
    addNode(mp5Parameters->fitted,(void *)bestAtom);

	findLambda(mp5Parameters->oneChannelResidueEnergy);
	
	if(mp5Parameters->MPType & MMP2)
	{
		double *signalInParticularChannel = NULL;
		for(channel=0;channel<mp5Parameters->numberOfAnalysedChannels;channel++)
		{
			returnAmplitudeAndModulusForMMP2DI(mp5Parameters,dictionary,bestAtom,&tmpAmplitude,&tmpModulus,channel);
			totalMMP2Modulus = totalMMP2Modulus + tmpModulus*tmpModulus;
			signalInParticularChannel = *(mp5Parameters->multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);
			*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(signalInParticularChannel,offsetExpandedDimension);
			totalMMP2Residue = totalMMP2Residue + (*(mp5Parameters->residueEnergyInEachChannel + channel));
		}
	}
	
    t2 = Clock();

	if(applicationMode & PROCESS_USER_MODE)
	{
		printf(" ATOM: [%-3d], SEC.: %-6.2lf, SIG: %-6.2lf, MOD: %-6.2lf, RES: %-6.2lf, RES/SIG: %-6.2lf %% LAMBDA %6.2lf",1,(t2-t1),mp5Parameters->oneChannelSignalEnergy,pow(*(mp5Parameters->bestModulusesTable),2.0),mp5Parameters->oneChannelResidueEnergy,(mp5Parameters->oneChannelResidueEnergy/mp5Parameters->oneChannelSignalEnergy)*100.0,0.0);
		if(mp5Parameters->MPType & MMP2)
			printf(" SIG: %-6.2lf, MOD: %-6.2lf, RES: %-6.2lf, RES/SIG: %-6.2lf %% ",mp5Parameters->totalSignalEnergy,totalMMP2Modulus,totalMMP2Residue,(totalMMP2Residue/mp5Parameters->totalSignalEnergy)*100.0);			
		
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

void nextIterationSMPMMP2(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned int atomsCounter = 0;
    unsigned int bestAtomIndex = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;   
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
    unsigned short int iterationCounter;

	unsigned short int channel;
	
	float tmpAmplitude = 0;
	float tmpModulus   = 0;
	double totalMMP2Modulus = 0.0;
	double totalMMP2Residue = 0.0;
	
	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    double *residueTable = mp5Parameters->singleChannelResidueTable;

    unsigned       int offsetExpandedDimension = mp5Parameters->offsetExpandedDimension;
	unsigned       int numberOfNonFFTAtoms     = dictionary->numberOfNonFFTAtoms;

    double t1, t2;
	double energyProgress, iterationProgress;
	double lambda;

    double *bestModulusesTable = mp5Parameters->bestModulusesTable;
	double modulus     = 0.0;
	bestModulus        = 0.0;  
	absOfBestModulus   = 0.0;

    Atom *atom, *bestAtom;

    stepInToolbar = dictionary->sizeOfDictionary/NUMBER_OF_STEPS_IN_TOOLBAR;

    iterationCounter = 1;

    while((mp5Parameters->oneChannelSignalEnergy>(mp5Parameters->oneChannelSignalEnergy - mp5Parameters->oneChannelSignalEnergy*mp5Parameters->energyPercent/100.0)) && iterationCounter<mp5Parameters->maximalNumberOfIterations)
    {
		if(applicationMode & PROCESS_SERVER_MODE)
		{
			energyProgress    = ((mp5Parameters->oneChannelSignalEnergy - mp5Parameters->oneChannelResidueEnergy)/mp5Parameters->oneChannelSignalEnergy)*100.0;
			iterationProgress = (iterationCounter/mp5Parameters->maximalNumberOfIterations)*100.0;
			printf("ATOM\t%3u\t%3u\t%6.2lf\t%6.2lf\n",iterationCounter,dictionary->sizeOfDictionary,energyProgress,iterationProgress);
		} 

		t1 = Clock();

		step             = 1;
		modulus          = 0.0;
		bestModulus      = 0.0;  
		absOfBestModulus = 0.0;
		atomsCounter     = 0;

		bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

		atom = dictionary->atomsTable;

		tmpAmplitude = 0;
		tmpModulus   = 0;
		totalMMP2Modulus = 0.0;
		totalMMP2Residue = 0.0;
		
		if(mp5Parameters->FFT)
		{
			if(dictionary->diracInDictionary || dictionary->gaussInDictionary)
			{
				for(atomsCounter=0;atomsCounter<numberOfNonFFTAtoms;atomsCounter++)
				{
					if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
					{
						findAtomDataDotProduct(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);

						if(findUnknowPhaseDI(atom,&modulus,0))
						{	
							if(bestModulusUpdate(modulus))
								bestAtomIndex = atomsCounter;
						}
					}
					atom++;
					printProgress(atomsCounter,mp5Parameters->progressBar);				
				}
			}
		
			if(dictionary->sinCosInDictionary)
			{
				findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);
				numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

				for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
				{
					if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
					{
						if(findUnknowPhaseDI(atom,&modulus,0))
						{
							if(bestModulusUpdate(modulus))
								bestAtomIndex = atomsCounter;
						}
					}
					atom++;
					atomsCounter++;
					printProgress(atomsCounter,mp5Parameters->progressBar);
				}
			}
			/* end of FFT code for SIN/COS waves */

			for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++)
			{
				numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);
				
				for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
				{
					findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);
					numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

					for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
					{
						if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
						{
							if(findUnknowPhaseDI(atom,&modulus,0))
							{
								if(bestModulusUpdate(modulus))
									bestAtomIndex = atomsCounter;
							}
						}
						
						atom++;
						atomsCounter++;
						printProgress(atomsCounter,mp5Parameters->progressBar);
					}
				}
			}
		}
		else
		{
			for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
			{
				if(!(atom->feature & GABOR_WAS_HIT) && !(atom->feature & INCORRECTGABOR))
				{
					findAtomDataDotProduct(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);

					if(findUnknowPhaseDI(atom,&modulus,0))
					{
						if(bestModulusUpdate(modulus))
							bestAtomIndex = atomsCounter;
					}
				}
				atom++;
				printProgress(atomsCounter,mp5Parameters->progressBar);
			}
		}
	
		atom = dictionary->atomsTable + bestAtomIndex;

		mp5Parameters->previousAtom = atom;
		*bestModulusesTable = bestModulus;
		atom->feature|=GABOR_WAS_HIT;
		
		copyAtom(atom,bestAtom,mp5Parameters->numberOfAllocatedChannels);

		makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);

		makeAtomTable(mp5Parameters,bestAtom,0);

		findResidue(residueTable,*(mp5Parameters->prevAtomTable),*bestModulusesTable,offsetExpandedDimension);
		
		mp5Parameters->oneChannelResidueEnergy = findSignalEnergy(residueTable,offsetExpandedDimension);

		addNode(mp5Parameters->fitted,(void *)bestAtom);
		
		iterationCounter++;

		lambda = findLambda(mp5Parameters->oneChannelResidueEnergy);
	
		if(mp5Parameters->MPType & MMP2)
		{
			double *signalInParticularChannel = NULL;
			for(channel=0;channel<mp5Parameters->numberOfAnalysedChannels;channel++)
			{
				returnAmplitudeAndModulusForMMP2DI(mp5Parameters,dictionary,bestAtom,&tmpAmplitude,&tmpModulus,channel);
				totalMMP2Modulus = totalMMP2Modulus + tmpModulus*tmpModulus;
				signalInParticularChannel = *(mp5Parameters->multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1);
				*(mp5Parameters->residueEnergyInEachChannel + channel) = findSignalEnergy(signalInParticularChannel,offsetExpandedDimension);
				totalMMP2Residue = totalMMP2Residue + (*(mp5Parameters->residueEnergyInEachChannel + channel));
			}
		}

		t2 = Clock();

		if(applicationMode & PROCESS_USER_MODE)
		{
			printf(" ATOM: [%-3d], SEC.: %-6.2lf, SIG: %-6.2lf, MOD: %-6.2lf, RES: %-6.2lf, RES/SIG: %-6.2lf %% LAMBDA %6.2lf",iterationCounter,(t2-t1),mp5Parameters->oneChannelSignalEnergy,pow((*(mp5Parameters->bestModulusesTable)),2.0),mp5Parameters->oneChannelResidueEnergy,(mp5Parameters->oneChannelResidueEnergy/mp5Parameters->oneChannelSignalEnergy)*100.0,lambda);
			if(mp5Parameters->MPType & MMP2)
				printf(" SIG: %-6.2lf, MOD: %-6.2lf, RES: %-6.2lf, RES/SIG: %-6.2lf %% ",mp5Parameters->totalSignalEnergy,totalMMP2Modulus,totalMMP2Residue,(totalMMP2Residue/mp5Parameters->totalSignalEnergy)*100.0);			

			if(atom->feature & DIRACDELTA)	
				printf(" D\n");
			else if(atom->feature & GAUSSFUNCTION)
				printf(" N\n");
			else if(atom->feature & SINCOSWAVE)	
				printf(" H\n");
			else if(atom->feature & GABORWAVE)
				printf(" G\n");
		
			printf("\n");
		}
		fflush(stdout);
	}
}
