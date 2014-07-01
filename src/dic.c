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
 
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include"atom.h"
#include"def.h"
#include"dic.h"
#include"matrix.h"
#include"mp5.h"
#include"random.h"
#include"types.h"
#include"vector.h"

#define kiloByte 1024.0
#define megaByte 1048576.0
#define dOmegaCorrection (0.606)
#define dUCorrection     (0.678)

extern unsigned char applicationMode;

static void setPermuteTable(unsigned int *permuteTable, unsigned int numberOfAtoms)
{
	unsigned int atomNumber;

	for(atomNumber = 0;atomNumber<numberOfAtoms;atomNumber++)
		permuteTable[atomNumber] = atomNumber;
}

static double estimateDilationFactor(double energyError)
{
	double a;
	const double sqrEnergy = energyError*energyError;
	const double firstTerm  = (2.0 - sqrEnergy);
	const double secondTerm = (sqrEnergy*sqrEnergy - 2.0*sqrEnergy + 2.0);
	const double thirdterm  = (1.0 - sqrEnergy)*(1.0 - sqrEnergy);

	a = (1.0 + energyError*sqrt(firstTerm*secondTerm))/thirdterm;

	return a;
}

/*static double estimateTInOptimalDictionary(double gaborScale, double dilationFactor)
{
	const double constA = log(0.5*(dilationFactor + 1.0/dilationFactor));
	const double dOmega = dOmegaCorrection*sqrt((4*M_PI/(gaborScale*gaborScale))*constA);
	return (unsigned int)(((2*M_PI)/dOmega) + 0.5);
}*/

static double estimateTInOptimalDictionary(double gaborScale, double energyError)
{
	const double constA = sqrt(-8.0*M_PI*log(1.0 - energyError*energyError));
	const double dOmega = constA/gaborScale;
	return (unsigned int)(((2.0*M_PI)/dOmega) + 0.5);
}

/*
static double estimateDUInOptimalDictionary(double gaborScale, double dilationFactor)
{
	const double constA = log(0.5*(dilationFactor + 1.0/dilationFactor));
	const double  dU = dUCorrection*sqrt(((gaborScale*gaborScale)/M_PI)*constA);
	return dU < 1.0 ? 1.0: dU;
}
*/

static double estimateDUInOptimalDictionary(double gaborScale, double energyError)
{
	const double constA = sqrt(-2.0*log(1.0 - 2.0*energyError*energyError)/M_PI);
	const double  dU = gaborScale*constA;
	return dU < 1.0 ? 1.0: dU;
}

static unsigned int estimateNumberOfStepsInFrequency(double dOmega)
{
	return (unsigned int)((M_PI - dOmega)/dOmega + 0.5);
}

static unsigned int estimateNumberOfStepsInPosition(unsigned int epochSize,double dU)
{
	return (unsigned int)((epochSize - dU)/dU + 0.5 );
}

static unsigned int findInterval(double gaborScale, unsigned char mode)
{
	unsigned int intervalRange;

	if(mode & GAUSS_NON_GAUSS)
		intervalRange    = (unsigned int)(sqrt(-LOG_EPS_DOT_PRODUCT/M_PI)*gaborScale + 1.5);
	else if(mode & GAUSS_GAUSS)
		intervalRange  = (unsigned int)(sqrt(-LOG_EPS_DOT_PRODUCT/M_2PI)*gaborScale + 1.5);

	return intervalRange;
}

void setMP5Sizes(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
	/* in this place one has to estimate the size of the tables, required by mp5 algorithms.
	    The gabor function presented in following way:
	    g = K*exp(-pi*(t-u)^2/s^2)*exp(-2piwt)
	    has scale = s/sqrt(pi)
	    If we want to get accurate results, the decomposition should be procceded in the range at least 4*scale
 	    In case of FFT optimalisation the FFT frequencies must be equal to optimal dictionary frequencies:
		A = ln(0.5*(a + a^-1))
		DF_fft = 2pi/N
		DF = (2/s)sqrt(pi*A)
		DF = DF_fft
		2pi/N = (2/s)sqrt(pi*A)
		N = s*sqrt(pi/A)

		On the other hand N = s*sqrt(pi/A) is a basic period for the particular scale.
		If one reads this file, he founds, that the smaller stepi in period (and the in frequency = 2pi/period)
		is equal to sqrt(pi/A). Periods for scales are equal to s*sqrt(pi/A).
		To sum it up, if one has gabor period, this is  at the same time a size of FFT at which
		one gets frequencies for optimalal dictionary.
	*/

	// maximal scale is set to maximal scale o gabor
	double maximalScale  = dictionary->tableOfScalesInOptimalDictionary[dictionary->numberOfStepsInScale - 1];
	double maximalPeriod = dictionary->tableOfPeriodsInOptimalDictionary[dictionary->numberOfStepsInScale - 1];

	/* Now we check, whether the FFT size is at lest equal to or larger then range, for which dot product has required accurancy */

	unsigned int intervalRange;

	intervalRange = findInterval(maximalScale,GAUSS_NON_GAUSS);
    double maximalFFTSizeForGaborAtoms = 0.0;

    if(mp5Parameters->FFT)
    {
        const double factor   = (2.0*intervalRange + 1)/maximalPeriod; // the dot product is calculated in the range:  position of atom +-interval -> 2*interval + 1 points*/

        if(factor>1.0)
            maximalFFTSizeForGaborAtoms = ((unsigned int)(factor + 1.0))*maximalPeriod;
        else
            maximalFFTSizeForGaborAtoms = maximalPeriod;
    }
    else
        maximalFFTSizeForGaborAtoms = 2.0*intervalRange + 1;

    mp5Parameters->marginalSize       = maximalFFTSizeForGaborAtoms/2 + 1;
    mp5Parameters->exponensTableSize  = mp5Parameters->epochSize + mp5Parameters->marginalSize;
    mp5Parameters->epochExpandedSize  = mp5Parameters->epochSize + 2*mp5Parameters->marginalSize;
    mp5Parameters->fftTableSize       = mp5Parameters->epochExpandedSize;

/*    printf(" period: %f %u\n",maximalPeriod,intervalRange);
	printf(" mp5Parameters->marginalSize: %u\n",mp5Parameters->marginalSize);
	printf(" mp5Parameters->exponensTableSize:  %u\n",mp5Parameters->exponensTableSize);
	printf(" mp5Parameters->marginalSize:       %u\n",mp5Parameters->marginalSize);
	printf(" mp5Parameters->epochExpandedSize:  %u\n",mp5Parameters->epochExpandedSize);
	printf(" mp5Parameters->fftTableSize:       %u\n",mp5Parameters->fftTableSize);
	printf(" calkowity rozmiar tablicy wynosi %u \n",mp5Parameters->epochExpandedSize);
	fflush(stdout);*/
}

static void findIntegerScalesForEnergyErrorParameter(MP5Parameters *mp5Parameters, Dictionary *dictionary)
{
	const double       numerator = log(mp5Parameters->epochSize);
	const double       dilationFactor = estimateDilationFactor(dictionary->energyError);
	unsigned short int scaleIndex;
	double             denominator;
	unsigned short int numberOfStepsInScale;

    denominator                      = log(dilationFactor);
    numberOfStepsInScale             = (unsigned short int)(numerator/denominator + 0.5);	
	dictionary->numberOfStepsInScale = numberOfStepsInScale;

	dictionary->tableOfScalesInOptimalDictionary       = (double *)dVectorAllocate(numberOfStepsInScale);
	dictionary->tableOfPeriodsInOptimalDictionary      = (unsigned int *)uiVectorAllocate(numberOfStepsInScale);
	dictionary->tableOfFrequenciesInOptimalDictionary  = (double *)dVectorAllocate(numberOfStepsInScale);

	dictionary->numberOfStepsInFrequencyAtParticularScale = uiVectorAllocate(numberOfStepsInScale);
	dictionary->numberOfStepsInPositionAtParticularScale  = uiVectorAllocate(numberOfStepsInScale);

	for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
        dictionary->tableOfScalesInOptimalDictionary[scaleIndex] = pow(dilationFactor,(scaleIndex + 1.0));

}

static BOOLEAN testScaleToPeriodFactorAtom(const Dictionary *dictionary, unsigned short int scaleIndex, unsigned int rifling, double scaleToPeriodFactor)
{

	const unsigned int period = (unsigned int)((*(dictionary->tableOfPeriodsInOptimalDictionary + scaleIndex))/((double)rifling) + 0.5);
	const double       scale  = (*(dictionary->tableOfScalesInOptimalDictionary  + scaleIndex))/sqrt(M_PI);

	if(scale<=scaleToPeriodFactor*period)
		return FALSE;

	return TRUE;
} 

/*
#ifdef DICTIONARY_DUMPING
	static void atomToString(const MP5Parameters *mp5Parameters, const Dictionary *dictionary, const Atom *atom, FILE* file)
	{
		fprintf(file," scale:           %f\n",dictionary->tableOfScalesInOptimalDictionary[atom->scaleIndex]/sqrt(M_PI));
		fprintf(file," scale:           %f\n",dictionary->tableOfScalesInOptimalDictionary[atom->scaleIndex]);		
		fprintf(file," position:        %u\n",atom->position);
		fprintf(file," rifling:         %u\n",atom->rifling);
		fprintf(file," frequency:       %f\n",(atom->rifling + 1)*M_2PI/(dictionary->tableOfPeriodsInOptimalDictionary[atom->scaleIndex]));
		fprintf(file," basicFrequency:  %f\n",dictionary->tableOfFrequenciesInOptimalDictionary[atom->scaleIndex]);
		fprintf(file," basicPeriod:     %u\n",dictionary->tableOfPeriodsInOptimalDictionary[atom->scaleIndex]);
		fprintf(file," \n");
	}
#endif
*/
static void printSizeOfDictionaryAndSizeOfSinCosExpTables(const Dictionary *dictionary, const MP5Parameters *mp5Parameters)
{
    unsigned short int scaleIndex;
    unsigned long int numberOfPointsInSinCosTables = 0;
    unsigned long int numberOfPointsInExpTables    = 0;
    unsigned long int totalNumberOfPoints          = 0;

    unsigned int sizeOfOneAtom = estimateSizeOfAtom(mp5Parameters->numberOfAllocatedChannels,mp5Parameters->MPType);

    for(scaleIndex=0;scaleIndex<dictionary->numberOfStepsInScale;scaleIndex++)
		numberOfPointsInSinCosTables+= (*(dictionary->tableOfPeriodsInOptimalDictionary + scaleIndex));

    numberOfPointsInExpTables = mp5Parameters->epochExpandedSize*dictionary->numberOfStepsInScale;

    totalNumberOfPoints = 2*numberOfPointsInSinCosTables + numberOfPointsInExpTables;

    printf(" INFORMATION ABOUT INITIAL DICTIONARY: \n");
    printf(" \n");
    printf(" NUMBER OF DIRAC'S DELTAS:                      %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfInitialDiracFunctions,dictionary->numberOfInitialDiracFunctions*sizeOfOneAtom/megaByte);
    printf(" NUMBER OF GAUSS FUNCTIONS:                     %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfInitialGaussFunctions,dictionary->numberOfInitialGaussFunctions*sizeOfOneAtom/megaByte);
    printf(" NUMBER OF SIN/COS:                             %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfInitialSinCosFunctions,dictionary->numberOfInitialSinCosFunctions*sizeOfOneAtom/megaByte);
    printf(" NUMBER OF ALL GABORS:                          %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfInitialGaborFunctions,dictionary->numberOfInitialGaborFunctions*sizeOfOneAtom/megaByte);
    printf(" NUMBER OF CORRECT GABORS:                      %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfInitialCorrectGabors,dictionary->numberOfInitialCorrectGabors*sizeOfOneAtom/megaByte);
    printf(" TOTAL NUMBER OF ATOMS (WITH INCORRECT GABORS): %12u, SIZE: %10.3lf (MB) \n",dictionary->initialNumberOfAtoms,(1.0*dictionary->initialNumberOfAtoms*sizeOfOneAtom)/megaByte);
    printf(" \n");

    printf(" INFORMATION ABOUT FINAL DICTIONARY: \n");
    printf(" \n");
    printf(" NUMBER OF DIRAC'S DELTAS:                      %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfFinalDiracFunctions,dictionary->numberOfFinalDiracFunctions*sizeOfOneAtom/megaByte);
    printf(" NUMBER OF GAUSS FUNCTIONS:                     %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfFinalGaussFunctions,dictionary->numberOfFinalGaussFunctions*sizeOfOneAtom/megaByte);
    printf(" NUMBER OF SIN/COS:                             %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfFinalSinCosFunctions,dictionary->numberOfFinalSinCosFunctions*sizeOfOneAtom/megaByte);
    printf(" NUMBER OF ALL GABORS:                          %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfFinalGaborFunctions,dictionary->numberOfFinalGaborFunctions*sizeOfOneAtom/megaByte);
//    printf(" NUMBER OF CORRECT GABORS:                      %12u, SIZE: %10.3lf (MB) \n",dictionary->numberOfFinalCorrectGabors,dictionary->numberOfFinalCorrectGabors*sizeOfOneAtom/megaByte);
    printf(" TOTAL NUMBER OF ATOMS (WITH INCORRECT GABORS): %12u, SIZE: %10.3lf (MB) \n",dictionary->finalNumberOfAtoms,(1.0*dictionary->finalNumberOfAtoms*sizeOfOneAtom)/megaByte);
    printf(" \n");
    printf(" SIZE OF SIN TABLES:                            %10.3lf (MB)\n",(1.0*numberOfPointsInSinCosTables*sizeof(double))/megaByte);
    printf(" SIZE OF COS TABLES:                            %10.3lf (MB)\n",(1.0*numberOfPointsInSinCosTables*sizeof(double))/megaByte);
    printf(" SIZE OF EXP TABLES:                            %10.3lf (MB)\n",(1.0*numberOfPointsInExpTables*sizeof(double))/megaByte);
    printf(" SIZE OF ALOCATED SIN/COS/EXP TABLES:           %10.3lf (MB)\n",(1.0*totalNumberOfPoints*sizeof(double))/megaByte);
    printf(" \n");
    printf(" TOTAL MEMORY USAGE BY DICTIONARY:              %10.3lf (MB)\n",(1.0*dictionary->initialNumberOfAtoms*sizeOfOneAtom + 1.0*totalNumberOfPoints*sizeof(double))/megaByte);

    fflush(stdout);
}

void analyseDictionarySizeAndType(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
	if(applicationMode & PROCESS_USER_MODE)
		printf("\n ANALYSIS OF DICTIONARY TYPE AND DICTIONARY SIZE \n");

	double gaborScale;

	unsigned int epochSize   = mp5Parameters->epochSize;
    const double energyError = dictionary->energyError;

	unsigned int gaborsCounter = 0U;

    findIntegerScalesForEnergyErrorParameter(mp5Parameters,dictionary);

	const unsigned short int numberOfStepsInScale = dictionary->numberOfStepsInScale;

	double       *tableOfScalesInOptimalDictionary      = dictionary->tableOfScalesInOptimalDictionary;
    unsigned int *tableOfPeriodsInOptimalDictionary     = dictionary->tableOfPeriodsInOptimalDictionary;
	double       *tableOfFrequenciesInOptimalDictionary = dictionary->tableOfFrequenciesInOptimalDictionary;

    unsigned int numberOfStepsInFrequency = 0;
	unsigned int numberOfStepsInPosition  = 0;

	unsigned short int scaleIndex  = 0;
	double             omega       = 0.0;
    double             dOmega      = 0.0;
	unsigned int       dT          = 0;
	double             dU          = 0.0;

	unsigned int frequencyStepNumber;

	for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
	{
		gaborScale = tableOfScalesInOptimalDictionary[scaleIndex];
		dT         = *(tableOfPeriodsInOptimalDictionary     + scaleIndex) = estimateTInOptimalDictionary(gaborScale,energyError);
		*(tableOfFrequenciesInOptimalDictionary + scaleIndex) = M_2PI/dT;

		dU = estimateDUInOptimalDictionary(gaborScale,energyError);

		*(dictionary->numberOfStepsInFrequencyAtParticularScale + scaleIndex) = estimateNumberOfStepsInFrequency(M_2PI/dT);
		*(dictionary->numberOfStepsInPositionAtParticularScale + scaleIndex)  = estimateNumberOfStepsInPosition(epochSize,dU);

	}

	if(dictionary->diracInDictionary)
    {
        dictionary->numberOfInitialDiracFunctions = epochSize;

        if((dictionary->typeOfDictionary & OCTAVE_FIXED))
        	dictionary->numberOfFinalDiracFunctions = dictionary->numberOfInitialDiracFunctions;
        else // stochastic dictionary
        	dictionary->numberOfFinalDiracFunctions = round(dictionary->numberOfInitialDiracFunctions*dictionary->stochasticDictionaryReductionCoefficient);
    }
	else
    {
		dictionary->numberOfInitialDiracFunctions = 0;
		dictionary->numberOfFinalDiracFunctions   = 0;
    }

	if(dictionary->gaussInDictionary)
    {
		dictionary->numberOfInitialGaussFunctions = numberOfStepsInScale*epochSize;

        if((dictionary->typeOfDictionary & OCTAVE_FIXED))
        	dictionary->numberOfFinalGaussFunctions = dictionary->numberOfInitialGaussFunctions;
        else // stochastic dictionary
        	dictionary->numberOfFinalGaussFunctions = round(dictionary->numberOfInitialGaussFunctions*dictionary->stochasticDictionaryReductionCoefficient);
    }
	else
    {
		dictionary->numberOfInitialGaussFunctions = 0;
		dictionary->numberOfFinalGaussFunctions   = 0;
    }

	if(dictionary->sinCosInDictionary)
	{
		dOmega     = *(tableOfFrequenciesInOptimalDictionary + numberOfStepsInScale - 1);
		numberOfStepsInFrequency = estimateNumberOfStepsInFrequency(dOmega);

		dictionary->numberOfInitialSinCosFunctions = numberOfStepsInFrequency;

        if((dictionary->typeOfDictionary & OCTAVE_FIXED))
        	dictionary->numberOfFinalSinCosFunctions = dictionary->numberOfInitialSinCosFunctions;
        else // stochastic dictionary
        	dictionary->numberOfFinalSinCosFunctions = round(dictionary->numberOfInitialSinCosFunctions*dictionary->stochasticDictionaryReductionCoefficient);

	}
	else
    {
		dictionary->numberOfInitialSinCosFunctions = 0;
		dictionary->numberOfFinalSinCosFunctions   = 0;
    }

	if(dictionary->gaborInDictionary)
    {
		for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
		{
			gaborScale = tableOfScalesInOptimalDictionary[scaleIndex];
			dOmega = *(tableOfFrequenciesInOptimalDictionary + scaleIndex);
			numberOfStepsInFrequency = estimateNumberOfStepsInFrequency(dOmega);

			for(frequencyStepNumber = 0;frequencyStepNumber < numberOfStepsInFrequency;frequencyStepNumber++)
			{
				omega = (frequencyStepNumber + 1)*dOmega;
				dU = estimateDUInOptimalDictionary(gaborScale,energyError);
				numberOfStepsInPosition = estimateNumberOfStepsInPosition(epochSize,dU);
				gaborsCounter+=numberOfStepsInPosition;
			}
		}
	}

	if(dictionary->gaborInDictionary)
    {
		dictionary->numberOfInitialGaborFunctions = gaborsCounter;

        if((dictionary->typeOfDictionary & OCTAVE_FIXED))
        	dictionary->numberOfFinalGaborFunctions = dictionary->numberOfInitialGaborFunctions;
        else // stochastic dictionary
        	dictionary->numberOfFinalGaborFunctions = round(dictionary->numberOfInitialGaborFunctions*dictionary->stochasticDictionaryReductionCoefficient);

    }
	else
    {
		dictionary->numberOfInitialGaborFunctions = 0;
		dictionary->numberOfFinalGaborFunctions = 0;
    }
		
	dictionary->initialNumberOfAtoms = dictionary->numberOfInitialDiracFunctions  + /* number of Dirac Delta      */
							   	       dictionary->numberOfInitialGaussFunctions  + /* number of  Gauss Functions */
								       dictionary->numberOfInitialSinCosFunctions + /* number of  Sin/Cos Waves   */
								       dictionary->numberOfInitialGaborFunctions;   /* number of Gabor Waves      */

	dictionary->finalNumberOfAtoms = dictionary->numberOfFinalDiracFunctions  + /* number of Dirac Delta      */
							   	     dictionary->numberOfFinalGaussFunctions  + /* number of  Gauss Functions */
								     dictionary->numberOfFinalSinCosFunctions + /* number of  Sin/Cos Waves   */
								     dictionary->numberOfFinalGaborFunctions;   /* number of Gabor Waves      */


	setMP5Sizes(dictionary,mp5Parameters);

}

void allocateDictionary(Dictionary *dictionary, const MP5Parameters *mp5Parameters)
{
    unsigned int             atomsCounter = 0;
    const unsigned short int numberOfAllocatedChannels = mp5Parameters->numberOfAllocatedChannels;
    Atom *atomPointer = NULL;

    dictionary->atomsTable              = (Atom **)malloc(4*sizeof(Atom*));

    dictionary->atomsTable[0] = (Atom *)malloc(dictionary->numberOfInitialDiracFunctions*sizeof(Atom));
    dictionary->atomsTable[1] = (Atom *)malloc(dictionary->numberOfInitialGaussFunctions*sizeof(Atom));
    dictionary->atomsTable[2] = (Atom *)malloc(dictionary->numberOfInitialSinCosFunctions*sizeof(Atom));
    dictionary->atomsTable[3] = (Atom *)malloc(dictionary->numberOfInitialGaborFunctions*sizeof(Atom));

	dictionary->tmpRandomTable = uiVectorAllocate(dictionary->initialNumberOfAtoms);

    dictionary->diracAtomsTable  = dictionary->atomsTable[0];
	dictionary->gaussAtomsTable  = dictionary->atomsTable[1];
	dictionary->sinCosAtomsTable = dictionary->atomsTable[2];
	dictionary->gaborAtomsTable  = dictionary->atomsTable[3];

    for(atomsCounter=0;atomsCounter<dictionary->numberOfInitialDiracFunctions;atomsCounter++)
    {
		atomPointer = &(dictionary->diracAtomsTable[atomsCounter]);
		allocateAtomElements(atomPointer,numberOfAllocatedChannels,mp5Parameters->MPType);
    }

    for(atomsCounter=0;atomsCounter<dictionary->numberOfInitialGaussFunctions;atomsCounter++)
    {
		atomPointer = &(dictionary->gaussAtomsTable[atomsCounter]);
		allocateAtomElements(atomPointer,numberOfAllocatedChannels,mp5Parameters->MPType);
    }

    for(atomsCounter=0;atomsCounter<dictionary->numberOfInitialSinCosFunctions;atomsCounter++)
    {
		atomPointer = &(dictionary->sinCosAtomsTable[atomsCounter]);
		allocateAtomElements(atomPointer,numberOfAllocatedChannels,mp5Parameters->MPType);
    }

    for(atomsCounter=0;atomsCounter<dictionary->numberOfInitialGaborFunctions;atomsCounter++)
    {
		atomPointer = &(dictionary->gaborAtomsTable[atomsCounter]);
		allocateAtomElements(atomPointer,numberOfAllocatedChannels,mp5Parameters->MPType);
    }
}

void freeDictionary(Dictionary *dictionary)
{
    unsigned int atomsCounter;

	if(dictionary->tableOfScalesInOptimalDictionary!=NULL)
		dVectorFree(dictionary->tableOfScalesInOptimalDictionary);
	if(dictionary->numberOfStepsInFrequencyAtParticularScale!=NULL)
		uiVectorFree(dictionary->numberOfStepsInFrequencyAtParticularScale);
	if(dictionary->numberOfStepsInPositionAtParticularScale!=NULL)
		uiVectorFree(dictionary->numberOfStepsInPositionAtParticularScale);
	if(dictionary->tableOfPeriodsInOptimalDictionary!=NULL)
		uiVectorFree(dictionary->tableOfPeriodsInOptimalDictionary);
	if(dictionary->tableOfFrequenciesInOptimalDictionary!=NULL)
		dVectorFree(dictionary->tableOfFrequenciesInOptimalDictionary);
	if(dictionary->tmpRandomTable!=NULL)
		uiVectorFree(dictionary->tmpRandomTable);

	if(dictionary->atomsTable!=NULL)
	{
	    for(atomsCounter=0;atomsCounter<dictionary->numberOfInitialDiracFunctions;atomsCounter++)
			freeAtomElements(&(dictionary->diracAtomsTable[atomsCounter]));

	    for(atomsCounter=0;atomsCounter<dictionary->numberOfInitialGaussFunctions;atomsCounter++)
			freeAtomElements(&(dictionary->gaussAtomsTable[atomsCounter]));

	    for(atomsCounter=0;atomsCounter<dictionary->numberOfInitialSinCosFunctions;atomsCounter++)
			freeAtomElements(&(dictionary->sinCosAtomsTable[atomsCounter]));

	    for(atomsCounter=0;atomsCounter<dictionary->numberOfInitialGaborFunctions;atomsCounter++)
			freeAtomElements(&(dictionary->gaborAtomsTable[atomsCounter]));

	    free(dictionary->atomsTable);
	}

}

void makeDictionary(Dictionary *dictionary, const MP5Parameters *mp5Parameters)
{
    unsigned int numberOfCorrectGabors = 0;
	unsigned int atomNumber            = 0;
    Atom        *atomPointer = NULL;

	double scaleToPeriodFactor = dictionary->scaleToPeriodFactor;
    const unsigned int epochSize = mp5Parameters->epochSize;
          unsigned int leftSize = 0, rightSize = 0;

	double energyError = dictionary->energyError;

	double *tableOfScalesInOptimalDictionary      = dictionary->tableOfScalesInOptimalDictionary;
	double *tableOfFrequenciesInOptimalDictionary = dictionary->tableOfFrequenciesInOptimalDictionary;

    const unsigned short int numberOfStepsInScale = dictionary->numberOfStepsInScale;

    unsigned int numberOfStepsInFrequency = 0;
    unsigned int numberOfStepsInPosition  = 0;

    unsigned short int scaleIndex     = 0;
	double             gaborScale     = 0;
    unsigned       int frequencyIndex = 0;
    unsigned       int positionIndex  = 0;
    unsigned       int rifling        = 0;

	double omega  = 0.0;
	double dOmega = 0.0;
    double dU     = 0.0;

	if(dictionary->typeOfDictionary & OCTAVE_STOCH)
	{
		if(dictionary->randomSeed == AUTO_RANDOM_SEED)
			mt_seed();
		else
			mt_seed32new(dictionary->randomSeed);
	}

	if(applicationMode & PROCESS_USER_MODE)
		printf("\n START DICTIONARY GENERAITING \n\n");

    /* set Dirac Delta */
	if(dictionary->diracInDictionary)
	{
	    atomPointer = dictionary->diracAtomsTable;
		atomPointer->feature = 0x0000;

		for(positionIndex=0;positionIndex<epochSize;positionIndex++)
		{
			atomPointer->scaleIndex  = 0; // although scaleIndex = 0 means scale = 1, we set this value,
										  // because this fild cand not be empty. However in case of Dirac Delta
									      // this fild is not used, so it will not cause any disturbance

			atomPointer->rifling  = 0;
			atomPointer->position = (unsigned int)positionIndex;
			atomPointer->feature|=DIRACDELTA;

			leftSize  = atomPointer->position;
			rightSize = (unsigned int)(epochSize - atomPointer->position - 1);

			if(leftSize<=rightSize)
				atomPointer->feature|=LEFT_SIDE_POSITION_IN_EPOCH;

			atomPointer++;
						
		}

        if((dictionary->typeOfDictionary) & OCTAVE_STOCH)
        {
            atomPointer  = dictionary->diracAtomsTable;

			for(atomNumber = 0;atomNumber<dictionary->numberOfInitialDiracFunctions;atomNumber++)
				((atomPointer + atomNumber)->feature)&= (~STOCHASTIC_ATOM);

            atomPointer  = dictionary->diracAtomsTable;

			setPermuteTable(dictionary->tmpRandomTable,dictionary->numberOfInitialDiracFunctions);
			permuteTable(dictionary->tmpRandomTable,dictionary->numberOfInitialDiracFunctions);

			for(atomNumber = 0;atomNumber<dictionary->numberOfFinalDiracFunctions;atomNumber++)
				((atomPointer + dictionary->tmpRandomTable[atomNumber])->feature)|= STOCHASTIC_ATOM;
        }
		else
		{
            atomPointer = dictionary->diracAtomsTable;

			for(atomNumber = 0;atomNumber<dictionary->numberOfInitialDiracFunctions;atomNumber++)
				((atomPointer + atomNumber)->feature) & ~(STOCHASTIC_ATOM);
		}
	}

	if(dictionary->gaussInDictionary)
	{
	    atomPointer = dictionary->gaussAtomsTable;
		atomPointer->feature = 0x0000;

		/* set Gauss Functions */
		for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++) // scaleindex form 0 to numberOfStepsInScale-2 is for gauss/gabor function
		{																   // scale for scaleIndex= numberOfStepsInScale-1 us for sin/cos wave
            for(positionIndex=0;positionIndex<epochSize;positionIndex++)
            {
                atomPointer->scaleIndex = scaleIndex;
				atomPointer->position   = (unsigned int)(dU*(positionIndex + 1));

				atomPointer->feature|=GAUSSFUNCTION;
				leftSize  = atomPointer->position;
				rightSize = (unsigned int)(epochSize - atomPointer->position - 1);

				if(leftSize<=rightSize)
					atomPointer->feature|=LEFT_SIDE_POSITION_IN_EPOCH;

				atomPointer++;
			}
		}

        if((dictionary->typeOfDictionary) & OCTAVE_STOCH)
        {
            atomPointer  = dictionary->gaussAtomsTable;

			for(atomNumber = 0;atomNumber<dictionary->numberOfInitialGaussFunctions;atomNumber++)
				((atomPointer + atomNumber)->feature)&= (~STOCHASTIC_ATOM);

            atomPointer  = dictionary->gaussAtomsTable;

			setPermuteTable(dictionary->tmpRandomTable,dictionary->numberOfInitialGaussFunctions);
			permuteTable(dictionary->tmpRandomTable,dictionary->numberOfInitialGaussFunctions);

			for(atomNumber = 0;atomNumber<dictionary->numberOfFinalGaussFunctions;atomNumber++)
				((atomPointer + dictionary->tmpRandomTable[atomNumber])->feature)|= STOCHASTIC_ATOM;
        }
		else
		{
            atomPointer = dictionary->gaussAtomsTable;

			for(atomNumber = 0;atomNumber<dictionary->numberOfInitialGaussFunctions;atomNumber++)
				((atomPointer + atomNumber)->feature) & ~(STOCHASTIC_ATOM);
		}
	}

	/* set Sin/Cos Waves */
	if(dictionary->sinCosInDictionary)
	{
	    atomPointer = dictionary->sinCosAtomsTable;	    
		atomPointer->feature = 0x0000;

		dOmega     = tableOfFrequenciesInOptimalDictionary[numberOfStepsInScale-1];
		numberOfStepsInFrequency = estimateNumberOfStepsInFrequency(dOmega);

		for(frequencyIndex=0,rifling=1;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++,rifling++)
		{
			atomPointer->scaleIndex = (unsigned short int)(numberOfStepsInScale - 1);  // although scaleIndex = numberOfScale - 1 means scale = epochSize, we set this value,
																					   // because this fild can not be empty. However in case of FFT WAVE
																				       // this fild is not used, so it will not cause any disturbance

            atomPointer->rifling    = rifling;

            if((epochSize%2)==0)
                atomPointer->position = (epochSize - 1)/2;
            else
                atomPointer->position = epochSize/2;

            atomPointer->feature|=SINCOSWAVE;

            leftSize  = atomPointer->position;
            rightSize = (unsigned int)(epochSize - atomPointer->position - 1);

            if(leftSize<=rightSize)
                atomPointer->feature|=LEFT_SIDE_POSITION_IN_EPOCH;

            atomPointer++;
		}

        if((dictionary->typeOfDictionary) & OCTAVE_STOCH)
        {
            atomPointer  = dictionary->sinCosAtomsTable;

			for(atomNumber = 0;atomNumber<dictionary->numberOfInitialSinCosFunctions;atomNumber++)
				((atomPointer + atomNumber)->feature)&= (~STOCHASTIC_ATOM);

            atomPointer  = dictionary->sinCosAtomsTable;

			setPermuteTable(dictionary->tmpRandomTable,dictionary->numberOfInitialSinCosFunctions);
			permuteTable(dictionary->tmpRandomTable,dictionary->numberOfInitialSinCosFunctions);

			for(atomNumber = 0;atomNumber<dictionary->numberOfFinalSinCosFunctions;atomNumber++)
				((atomPointer + dictionary->tmpRandomTable[atomNumber])->feature)|= STOCHASTIC_ATOM;
        }
		else
		{
            atomPointer = dictionary->sinCosAtomsTable;

			for(atomNumber = 0;atomNumber<dictionary->numberOfInitialSinCosFunctions;atomNumber++)
				((atomPointer + atomNumber)->feature) & ~(STOCHASTIC_ATOM);
		}
	}

	/* set Gabor waves */
	if(dictionary->gaborInDictionary)
	{
		atomPointer = dictionary->gaborAtomsTable;
		atomPointer->feature = 0x0000;

	    for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
		{

			gaborScale = tableOfScalesInOptimalDictionary[scaleIndex];

			dU = estimateDUInOptimalDictionary(gaborScale,energyError);
			numberOfStepsInPosition = estimateNumberOfStepsInPosition(epochSize,dU);

			for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
			{
				dOmega = tableOfFrequenciesInOptimalDictionary[scaleIndex];
    	        numberOfStepsInFrequency = estimateNumberOfStepsInFrequency(dOmega);

				for(frequencyIndex=0,rifling=1;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++,rifling=rifling+1)
				{
					omega = (frequencyIndex + 1)*dOmega;

				    atomPointer->scaleIndex = scaleIndex;
					atomPointer->rifling    = rifling;
					atomPointer->position   = (unsigned int)(dU*(positionIndex + 1));

					atomPointer->feature|=GABORWAVE;
					leftSize  = atomPointer->position;
					rightSize = (unsigned int)(epochSize - atomPointer->position - 1);

					if(leftSize<=rightSize)
						atomPointer->feature|=LEFT_SIDE_POSITION_IN_EPOCH;

					if(testScaleToPeriodFactorAtom(dictionary,scaleIndex,rifling,scaleToPeriodFactor))
						numberOfCorrectGabors++;
					else
                        atomPointer->feature|=INCORRECTGABOR;

                    atomPointer++;
				}
			}
	    }

        if((dictionary->typeOfDictionary) & OCTAVE_STOCH)
        {
            atomPointer = dictionary->gaborAtomsTable;

			for(atomNumber = 0;atomNumber<dictionary->numberOfInitialGaborFunctions;atomNumber++)
				((atomPointer + atomNumber)->feature)&= (~STOCHASTIC_ATOM);

            atomPointer                 = dictionary->gaborAtomsTable;

			setPermuteTable(dictionary->tmpRandomTable,dictionary->numberOfInitialGaborFunctions);
			permuteTable(dictionary->tmpRandomTable,dictionary->numberOfInitialGaborFunctions);

			for(atomNumber = 0;atomNumber<dictionary->numberOfFinalGaborFunctions;atomNumber++)
				((atomPointer + dictionary->tmpRandomTable[atomNumber])->feature)|= STOCHASTIC_ATOM;

        }
		else
		{
            atomPointer = dictionary->gaborAtomsTable;

			for(atomNumber = 0;atomNumber<dictionary->numberOfInitialGaborFunctions;atomNumber++)
				((atomPointer + atomNumber)->feature)|=STOCHASTIC_ATOM;
		}

		dictionary->numberOfInitialCorrectGabors = numberOfCorrectGabors;

/*		if((dictionary->typeOfDictionary) & OCTAVE_STOCH)
			dictionary->numberOfFinalCorrectGabors = numberOfRandomSelectedAtoms;
		else
			dictionary->numberOfFinalCorrectGabors = dictionary->numberOfInitialCorrectGabors;		*/
	}
		
/*
	atomPointer = dictionary->atomsTable;
	unsigned int atomsCounter;
	for(atomsCounter=0;atomsCounter<dictionary->initialNumberOfAtoms;atomsCounter++)
	{
		scaleIndex = atomPointer->scaleIndex;

		printf("%u %hu %f %hu %u %hu %hu %f %f\n",atomsCounter,
												 scaleIndex,
												 *(dictionary->tableOfScalesInOptimalDictionary + atomPointer->scaleIndex),
												 atomPointer->position,
												 atomPointer->rifling,
												 (atomPointer->feature & LEFT_SIDE_POSITION_IN_EPOCH),
												 (atomPointer->feature & DIRACDELTA)    |
												 (atomPointer->feature & GAUSSFUNCTION) |
												 (atomPointer->feature & SINCOSWAVE)    |
												 (atomPointer->feature & GABORWAVE),
												 *(tableOfFrequenciesInOptimalDictionary + scaleIndex),
												 *(tableOfPositionsInOptimalDictionary   + scaleIndex));

		getchar();

		atomPointer++;
	}
*/

/*
	if(applicationMode & PROCESS_USER_MODE)
		printSizeOfDictionaryAndSizeOfSinCosExpTables(dictionary,mp5Parameters);

	#ifdef DICTIONARY_DUMPING
		unsigned int atomNumber = 0;
		char nameOfResultsFile[LENGTH_OF_NAME_OF_RESULTS_FILE];
        unsigned short int charCounter;
        unsigned short int lengthOfDataFileWithOutExpand;

        for(charCounter=0;charCounter<LENGTH_OF_NAME_OF_RESULTS_FILE;charCounter++)
            *(nameOfResultsFile + charCounter) = '\0';

        char *dot = strrchr(mp5Parameters->nameOfDataFile,'.');
        int len;

        char tmpString[LENGTH_OF_TMP_STRING];

        bzero((void *)tmpString,LENGTH_OF_TMP_STRING);

        if(dot==NULL)
        {
            strcpy(nameOfResultsFile,mp5Parameters->nameOfDataFile);
			sprintf(nameOfResultsFile,"_dictionary.dump");
		}
        else
        {
            len = strlen(dot);
            strncat(nameOfResultsFile,mp5Parameters->nameOfDataFile,strlen(mp5Parameters->nameOfDataFile)-len);
            lengthOfDataFileWithOutExpand = (unsigned short int)strlen(nameOfResultsFile);
    		sprintf((nameOfResultsFile + lengthOfDataFileWithOutExpand),"_dictionary.dump");
        }

		FILE *file = fopen(nameOfResultsFile,"wt");

		if(dictionary->diracInDictionary)
		{
			for(atomNumber=0;atomNumber<dictionary->numberOfInitialDiracFunctions;atomNumber++)
            {
                fprintf(file,"%u\n",atomNumber);
				atomToString(mp5Parameters,dictionary,&dictionary->atomsDiracTable[atomNumber],file);
            }
		}

		if(dictionary->gaussInDictionary)
		{
			for(atomNumber=0;atomNumber<dictionary->numberOfInitialGaussFunctions;atomNumber++)
            {
                fprintf(file,"%u\n",atomNumber); 
            	atomToString(mp5Parameters,dictionary,&dictionary->atomsGaussTable[atomNumber],file);
            }
		}

		if(dictionary->sinCosInDictionary)
		{
            fprintf(file,"%u\n",atomNumber);
			for(atomNumber=0;atomNumber<dictionary->numberOfInitialSinCosFunctions;atomNumber++)
            {
                fprintf(file,"%u\n",atomNumber);
				atomToString(mp5Parameters,dictionary,&dictionary->atomsSinCosTable[atomNumber],file);
            }
		}

		if(dictionary->gaborInDictionary)
		{
			for(atomNumber=0;atomNumber<dictionary->numberOfInitialGaborFunctions;atomNumber++)
            {
                fprintf(file,"%u\n",atomNumber);
				atomToString(mp5Parameters,dictionary,&dictionary->atomsGaborTable[atomNumber],file);
            }
		}
		fflush(file);
		fclose(file);
	#endif
*/
}

void reinitDictionary(Dictionary *dictionary, const MP5Parameters *mp5Parameters)
{
    unsigned int atomNumber;

	if(dictionary->diracInDictionary)
	{
		for(atomNumber=0;atomNumber<dictionary->numberOfInitialDiracFunctions;atomNumber++)
			(dictionary->diracAtomsTable[atomNumber]).feature  = 0x00;
	}

	if(dictionary->gaussInDictionary)
	{
		for(atomNumber=0;atomNumber<dictionary->numberOfInitialGaussFunctions;atomNumber++)
			(dictionary->gaussAtomsTable[atomNumber]).feature    = 0x00;
	}

	if(dictionary->sinCosInDictionary)
	{
		for(atomNumber=0;atomNumber<dictionary->numberOfInitialSinCosFunctions;atomNumber++)
			(dictionary->sinCosAtomsTable[atomNumber]).feature    = 0x00;
	}

	if(dictionary->gaborInDictionary)
	{
		for(atomNumber=0;atomNumber<dictionary->numberOfInitialGaborFunctions;atomNumber++)
			(dictionary->gaborAtomsTable[atomNumber]).feature    = 0x00;
	}

    makeDictionary(dictionary,mp5Parameters);

}

/*
void testAtomFeature(Dictionary *dictionary)
{
    Atom     *atomPointer = NULL;
    unsigned int correctAtomsCounter = 0;
    unsigned int atomsCounter = 0;

	dictionary->numberOfCorrectGabors = 0;

	if(applicationMode & PROCESS_USER_MODE)
	{
		printf(" THE INCORRECT GABORS ARE BEING REMOVED FROM DICTIONARY \n");
		printf(" UNCORRECT GABORS ARE GABORS WHICH SCALE IS SMALLER THEN K*PERIOD WAVE IN GABOR\n");
		fflush(stdout);
	}

    for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
    {
		atomPointer = dictionary->atomsTable + atomsCounter;
		if(atomPointer->feature & GABORWAVE)
			if(testScaleToPeriodFactorAtom(dictionary,atomPointer,dictionary->scaleToPeriodFactor))
				dictionary->numberOfCorrectGabors++;
    }

    correctAtomsCounter = 0;

    for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
    {
		atomPointer = dictionary->atomsTable + atomsCounter;

		if(!(atomPointer->feature & INCORRECTGABOR))
			correctAtomsCounter++;
	}

	if(applicationMode & PROCESS_USER_MODE)
	{
		printf("\n AFTER PREPROCESSING, THE NUMBER OF CORRECT GABORS IS: %u\n",correctAtomsCounter);
		fflush(stdout);
	}
}
*/
void resetDictionary(Dictionary *dictionary)
{
    unsigned int atomNumber;

	if(dictionary->diracInDictionary)
	{
		for(atomNumber=0;atomNumber<dictionary->numberOfInitialDiracFunctions;atomNumber++)
			(dictionary->diracAtomsTable[atomNumber]).feature&=~ATOM_WAS_SELECTED;
	}

	if(dictionary->gaussInDictionary)
	{
		for(atomNumber=0;atomNumber<dictionary->numberOfInitialGaussFunctions;atomNumber++)
			(dictionary->gaussAtomsTable[atomNumber]).feature&=~ATOM_WAS_SELECTED;
	}

	if(dictionary->sinCosInDictionary)
	{
		for(atomNumber=0;atomNumber<dictionary->numberOfInitialSinCosFunctions;atomNumber++)
			(dictionary->sinCosAtomsTable[atomNumber]).feature&=~ATOM_WAS_SELECTED;
	}

	if(dictionary->gaborInDictionary)
	{
		for(atomNumber=0;atomNumber<dictionary->numberOfInitialGaborFunctions;atomNumber++)
			(dictionary->gaborAtomsTable[atomNumber]).feature&=~ATOM_WAS_SELECTED;
	}
}

inline Atom* getAtom(const Dictionary *dictionary, unsigned int atomIndex)
{
	Atom *atom = NULL;
	static char flag = NO;
	static unsigned int fromDirac;
	static unsigned int fromDiracToGauss;
	static unsigned int fromDiracToSinCos;

	if(!flag)
	{
		fromDirac         = dictionary->numberOfInitialDiracFunctions;
		fromDiracToGauss  = dictionary->numberOfInitialDiracFunctions + dictionary->numberOfInitialGaussFunctions;
		fromDiracToSinCos = fromDiracToGauss + dictionary->numberOfInitialSinCosFunctions;
		flag = YES;
	}

	if((dictionary->diracInDictionary) && (atomIndex<fromDirac))
		atom = &dictionary->diracAtomsTable[atomIndex];
	else if((dictionary->gaussInDictionary) && (atomIndex<fromDiracToGauss))
		atom = &dictionary->gaussAtomsTable[atomIndex - fromDirac];
	else if((dictionary->sinCosInDictionary) && (atomIndex<fromDiracToSinCos))
		atom = &dictionary->sinCosAtomsTable[atomIndex - fromDiracToGauss];
	else if((dictionary->gaborInDictionary))
		atom = &dictionary->gaborAtomsTable[atomIndex - fromDiracToSinCos];

	return atom;

}
