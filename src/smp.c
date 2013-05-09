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
#define _GNU_SOURCE

#include<float.h>
#include<math.h>
#include<stdlib.h>
#include"def.h"
#include"dic.h"
#include"atom.h"
#include"io_mp5.h"
#include"mp5.h"
#include"queue.h"
#include"smp.h"
#include"tools.h"
#include"types.h"

/* SMP - SINGLE CHANNEL MATCHING PURSUIT */

extern unsigned char applicationMode;

static BOOLEAN bestModulusUpdateL1(double bestModulusL1, double modulusL1)
{
	return modulusL1<bestModulusL1 ? TRUE : FALSE;
}

static BOOLEAN bestModulusUpdateL2(double bestModulusL2, double modulusL2)
{
	double absOfModulusL2     = fabs(modulusL2);
	double absOfBestModulusL2 = fabs(bestModulusL2);
	return absOfModulusL2>absOfBestModulusL2 ? TRUE : FALSE;
}

void firstIterationSMPMMP2(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned int atomsCounter    = 0;
    unsigned int tmpAtomsCounter = 0;
    unsigned int bestAtomIndex   = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
	Progress       progress;
	progress.applicationMode = applicationMode;
    progress.stepInToolbar   = dictionary->initialNumberOfAtoms/NUMBER_OF_STEPS_IN_TOOLBAR;
    progress.step            = 1;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    double *signalTable  = mp5Parameters->singleChannelSignalTable;
    double *residueTable = signalTable;
    mp5Parameters->singleChannelResidueTable = residueTable;

    unsigned int epochExpandedSize   = mp5Parameters->epochExpandedSize;
	/* Atoms assocciated with FFT are sin/cos waves and Gabor functions. The dot product of Dirac Deltas and Gauss functions with signal can not be astimated by means of FFT algorithm. */

    double t1, t2;

	double modulusL1     = 0.0;
	double modulusL2     = 0.0;
	double bestModulusL1 = FLT_MAX;
	double bestModulusL2 = 0.0;

    double *bestModulusesTableL2 = mp5Parameters->bestModulusesTableL2;

    Atom *atom;

    Atom *bestAtom = allocateAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

	mp5Parameters->totalSignalEnergy  = findSignalEnergy(signalTable,epochExpandedSize);
	mp5Parameters->totalResidueEnergy = mp5Parameters->totalSignalEnergy;

    t1 = Clock();

	if(applicationMode & PROCESS_SERVER_MODE)
	{
		printf("ATOM\t%3u\t%3u\t%6.2f\t%6.2f\n",0,dictionary->finalNumberOfAtoms,0.0,0.0);
		fflush(stdout);
	}

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
					findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);

					if(findUnknowPhaseDI(atom,&modulusL2,0))
					{
						if(mp5Parameters->normType & L1)
						{
							modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

							if(bestModulusUpdateL1(bestModulusL1,modulusL1))
							{
								bestModulusL1 = modulusL1;
								bestModulusL2 = modulusL2;
								bestAtomIndex = atomsCounter;
							}
						}
						else
						{
							if(bestModulusUpdateL2(bestModulusL2,modulusL2))
							{
								bestModulusL2 = modulusL2;
								bestAtomIndex = atomsCounter;
							}
						}
					}
				}
				atom++;
				printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
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
					findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);

					if(findUnknowPhaseDI(atom,&modulusL2,0))
					{
						if(mp5Parameters->normType & L1)
						{
							modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

							if(bestModulusUpdateL1(bestModulusL1,modulusL1))
							{
								bestModulusL1 = modulusL1;
								bestModulusL2 = modulusL2;
								bestAtomIndex = atomsCounter;
  							}
						}
						else
						{
							if(bestModulusUpdateL2(bestModulusL2,modulusL2))
							{
								bestModulusL2 = modulusL2;
								bestAtomIndex = atomsCounter;
							}
						}
					}
				}
				atom++;
				atomsCounter++;
				printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
			}
		}

		if(dictionary->sinCosInDictionary)
		{
			atom = dictionary->sinCosAtomsTable;

			/* from now the code is FFT for SIN/COS waves */
			findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);
			numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

			for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
			{
				if((atom->feature & INCORRECTGABOR)==0)
				{
					normAtomTable(dictionary,mp5Parameters,atom);
					if(findUnknowPhaseDI(atom,&modulusL2,0))
					{
						if(mp5Parameters->normType & L1)
						{
							modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

							if(bestModulusUpdateL1(bestModulusL1,modulusL1))
							{
								bestModulusL1 = modulusL1;
								bestModulusL2 = modulusL2;
								bestAtomIndex = atomsCounter;
 								}
						}
						else
						{
							if(bestModulusUpdateL2(bestModulusL2,modulusL2))
							{
								bestModulusL2 = modulusL2;
								bestAtomIndex = atomsCounter;
							}
						}
					}
				}
				atom++;
				atomsCounter++;
				printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
			}
			/* end of FFT code for SIN/COS waves */
		}

		if(dictionary->gaborInDictionary)
		{
			atom = dictionary->gaborAtomsTable;

			for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
			{
				numberOfStepsInPosition = *(numberOfStepsInPositionAtParticularScale + scaleIndex);

				for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
				{

					findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);
					numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

					for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
					{
						if((atom->feature & INCORRECTGABOR)==0)
						{
							normAtomTable(dictionary,mp5Parameters,atom);

							if(findUnknowPhaseDI(atom,&modulusL2,0))
							{									
								if(mp5Parameters->normType & L1)
								{
									modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

									if(bestModulusUpdateL1(bestModulusL1,modulusL1))
									{
										bestModulusL1 = modulusL1;
										bestModulusL2 = modulusL2;
										bestAtomIndex = atomsCounter;
  									}
								}
								else
								{
									if(bestModulusUpdateL2(bestModulusL2,modulusL2))
									{
										bestModulusL2 = modulusL2;
										bestAtomIndex = atomsCounter;
									}
								}
							}
						}
						atom++;
						atomsCounter++;
						printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
					}
				}
			}
		}
	}
	else // only for STOCHASTIC DICTIOANARY
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

			findAtomDataDotProduct(dictionary,mp5Parameters,atom,signalTable,0,FIRST_ITERATION);
			findUnknowPhaseDI(atom,&modulusL2,0);

			if(findUnknowPhaseDI(atom,&modulusL2,0))
			{
				if(mp5Parameters->normType & L1)
				{
					modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

					if(bestModulusUpdateL1(bestModulusL1,modulusL1))
					{
						bestModulusL1 = modulusL1;
						bestModulusL2 = modulusL2;
						bestAtomIndex = atomsCounter;
					}
				}
				else
				{
					if(bestModulusUpdateL2(bestModulusL2,modulusL2))
					{
						bestModulusL2 = modulusL2;
						bestAtomIndex = atomsCounter;
					}
				}
			}
			printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
		}
	}

    atom = getAtom(dictionary,bestAtomIndex);

    atom->feature|=ATOM_WAS_SELECTED;
    mp5Parameters->previousAtom = atom;
    *bestModulusesTableL2 = bestModulusL2;
    copyAtom(atom,bestAtom,mp5Parameters->numberOfAllocatedChannels);

    makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);
	makeAtomTable(*(mp5Parameters->prevAtomTable),mp5Parameters,bestAtom,0);
	findResidue(residueTable,*(mp5Parameters->prevAtomTable),*bestModulusesTableL2,epochExpandedSize);
	mp5Parameters->totalResidueEnergy = findSignalEnergy(residueTable,epochExpandedSize);

    addNode(mp5Parameters->fitted,(void *)bestAtom);

    t2 = Clock();

	if(applicationMode & PROCESS_USER_MODE)
	{
		printf(" ATOM: [%-3d], SEC.: %-6.2f, SIG: %-6.2f, MOD: %-6.2f, RES: %-6.2f, RES/SIG: %-6.2f %%",1,(t2-t1),mp5Parameters->totalSignalEnergy,pow(*(mp5Parameters->bestModulusesTableL2),2.0),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

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
    unsigned int atomsCounter    = 0;
    unsigned int tmpAtomsCounter = 0;
    unsigned int bestAtomIndex = 0;
	unsigned short int scaleIndex;
	unsigned       int positionIndex;
	unsigned 	   int frequencyIndex;
	unsigned       int numberOfStepsInFrequency;
	unsigned       int numberOfStepsInPosition;
    unsigned short int iterationCounter;

	Progress       progress;
	progress.applicationMode = applicationMode;
    progress.stepInToolbar   = dictionary->initialNumberOfAtoms/NUMBER_OF_STEPS_IN_TOOLBAR;
    progress.step            = 1;

	const unsigned short int numberOfStepsInScale                       = dictionary->numberOfStepsInScale;
	const unsigned       int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
	const unsigned       int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    double *residueTable       = mp5Parameters->singleChannelResidueTable;
	double energyStopCondition = mp5Parameters->totalSignalEnergy *(1 - mp5Parameters->energyPercent/100.0);    

    unsigned int epochExpandedSize   = mp5Parameters->epochExpandedSize;

    double t1, t2;
	double energyProgress, iterationProgress;

    double *bestModulusesTableL2 = mp5Parameters->bestModulusesTableL2;
	double modulusL1     = 0.0;
	double modulusL2     = 0.0;
	double bestModulusL1 = FLT_MAX;
	double bestModulusL2 = 0.0;

    Atom *atom, *bestAtom;

    iterationCounter = 1;

    while((mp5Parameters->totalResidueEnergy>energyStopCondition) && (iterationCounter<mp5Parameters->maximalNumberOfIterations))
    {
		if(applicationMode & PROCESS_SERVER_MODE)
		{
			energyProgress    = ((mp5Parameters->totalSignalEnergy - mp5Parameters->totalResidueEnergy)/mp5Parameters->totalSignalEnergy)*100.0;
			iterationProgress = (iterationCounter/mp5Parameters->maximalNumberOfIterations)*100.0;
			printf("ATOM\t%3u\t%3u\t%6.2f\t%6.2f\n",iterationCounter,dictionary->finalNumberOfAtoms,energyProgress,iterationProgress);
		}

		t1 = Clock();

	    progress.step    = 1;
		modulusL1        = 0.0;
		bestModulusL1    = FLT_MAX;
		modulusL2        = 0.0;
		bestModulusL2    = 0.0;

		atomsCounter     = 0;
		tmpAtomsCounter  = 0;

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
						findAtomDataDotProduct(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);

						if(findUnknowPhaseDI(atom,&modulusL2,0))
						{
							if(mp5Parameters->normType & L1)
							{
								modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

								if(bestModulusUpdateL1(bestModulusL1,modulusL1))
								{
									bestModulusL1 = modulusL1;
									bestModulusL2 = modulusL2;
									bestAtomIndex = atomsCounter;
  								}
							}
							else
							{
								if(bestModulusUpdateL2(bestModulusL2,modulusL2))
								{
									bestModulusL2 = modulusL2;
									bestAtomIndex = atomsCounter;
								}
							}
						}
					}
					atom++;
					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
				}
			}

			if(dictionary->gaussInDictionary)
			{
				atom = dictionary->gaussAtomsTable;

				for(tmpAtomsCounter=0;tmpAtomsCounter<dictionary->numberOfFinalGaussFunctions;tmpAtomsCounter++)
				{
					if(!(atom->feature & ATOM_WAS_SELECTED) && !(atom->feature & INCORRECTGABOR))
					{
						findAtomDataDotProduct(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);

						if(findUnknowPhaseDI(atom,&modulusL2,0))
						{
							if(mp5Parameters->normType & L1)
							{
								modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

								if(bestModulusUpdateL1(bestModulusL1,modulusL1));
								{
									bestModulusL1 = modulusL1;
									bestModulusL2 = modulusL2;
									bestAtomIndex = atomsCounter;
  								}
							}
							else
							{
								if(bestModulusUpdateL2(bestModulusL2,modulusL2))
								{
									bestModulusL2 = modulusL2;
									bestAtomIndex = atomsCounter;
								}
							}
						}
					}
					atom++;
					atomsCounter++;
					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
				}
			}

			if(dictionary->sinCosInDictionary)
			{
				atom = dictionary->sinCosAtomsTable;

				findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);
				numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);

				for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
				{
					if(!(atom->feature & ATOM_WAS_SELECTED) && !(atom->feature & INCORRECTGABOR))
					{
						if(findUnknowPhaseDI(atom,&modulusL2,0))
						{
							if(mp5Parameters->normType & L1)
							{
								modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

								if(bestModulusUpdateL1(bestModulusL1,modulusL1))
								{
									bestModulusL1 = modulusL1;
									bestModulusL2 = modulusL2;
									bestAtomIndex = atomsCounter;
  								}
							}
							else
							{
								if(bestModulusUpdateL2(bestModulusL2,modulusL2))
								{
									bestModulusL2 = modulusL2;
									bestAtomIndex = atomsCounter;
								}
							}
						}
					}
					atom++;
					atomsCounter++;
					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
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
						findGaborDataDotProductFFT(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);
						numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);

						for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
						{
							if(((atom->feature & ATOM_WAS_SELECTED)==0) && ((atom->feature & INCORRECTGABOR)==0))
							{
								if(findUnknowPhaseDI(atom,&modulusL2,0))
								{
									if(mp5Parameters->normType & L1)
									{
										modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

										if(bestModulusUpdateL1(bestModulusL1,modulusL1))
										{
											bestModulusL1 = modulusL1;
											bestModulusL2 = modulusL2;
											bestAtomIndex = atomsCounter;
  										}
									}
									else
									{
										if(bestModulusUpdateL2(bestModulusL2,modulusL2))
										{
											bestModulusL2 = modulusL2;
											bestAtomIndex = atomsCounter;
										}
									}
								}
							}

							atom++;
							atomsCounter++;
							printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
						}
					}
				}
			}
		}
		else // only for STOCHASTIC DICTIOANARY
		{
			for(atomsCounter=0;atomsCounter<dictionary->initialNumberOfAtoms;atomsCounter++)
			{
				atom = getAtom(dictionary,atomsCounter);

				if((atom->feature & ATOM_WAS_SELECTED) || (atom->feature & INCORRECTGABOR) || ((atom->feature & STOCHASTIC_ATOM)==0))
				{

/*					printf(" atom wybrany: %hu\n",atom->feature & ATOM_WAS_SELECTED);
					printf(" atom zepsuty: %hu\n",(atom->feature & INCORRECTGABOR));
					printf(" atom losowy: %hu\n",(atom->feature & STOCHASTIC_ATOM));
					printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);*/
					continue;
				}							

				findAtomDataDotProduct(dictionary,mp5Parameters,atom,*(mp5Parameters->prevAtomTable),0,NEXT_ITERATION);

				if(findUnknowPhaseDI(atom,&modulusL2,0))
				{
					if(mp5Parameters->normType & L1)
					{
						modulusL1 = findModulusL1(dictionary,mp5Parameters,atom);

						if(bestModulusUpdateL1(bestModulusL1,modulusL1))
						{
							bestModulusL1 = modulusL1;
							bestModulusL2 = modulusL2;
							bestAtomIndex = atomsCounter;
						}
					}
					else
					{
						if(bestModulusUpdateL2(bestModulusL2,modulusL2))
						{
							bestModulusL2 = modulusL2;
							bestAtomIndex = atomsCounter;
						}
					}
				}
				printInformationAboutProgress(mp5Parameters,&progress,atomsCounter);
			}
		}



		atom = getAtom(dictionary,bestAtomIndex);
		mp5Parameters->previousAtom = atom;
		*bestModulusesTableL2 = bestModulusL2;
		atom->feature|=ATOM_WAS_SELECTED;
		copyAtom(atom,bestAtom,mp5Parameters->numberOfAllocatedChannels);
		makeSinCosExpAtomTable(dictionary,mp5Parameters,bestAtom);
		makeAtomTable(*(mp5Parameters->prevAtomTable),mp5Parameters,bestAtom,0);
        findResidue(residueTable,*(mp5Parameters->prevAtomTable),*bestModulusesTableL2,epochExpandedSize);
		mp5Parameters->totalResidueEnergy = findSignalEnergy(residueTable,epochExpandedSize);
		addNode(mp5Parameters->fitted,(void *)bestAtom);
		iterationCounter++;

		t2 = Clock();

		if(applicationMode & PROCESS_USER_MODE)
		{
			printf(" ATOM: [%-3d], SEC.: %-6.2f, SIG: %-6.2f, MOD: %-6.2f, RES: %-6.2f, RES/SIG: %-6.2f %% ",iterationCounter,(t2-t1),mp5Parameters->totalSignalEnergy,pow((*(mp5Parameters->bestModulusesTableL2)),2.0),mp5Parameters->totalResidueEnergy,(mp5Parameters->totalResidueEnergy/mp5Parameters->totalSignalEnergy)*100.0);

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
