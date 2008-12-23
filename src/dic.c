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

#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include"atom.h"
#include"def.h"
#include"dic.h"
#include"matrix.h"
#include"mp5.h"
#include"r250.h"
#include"types.h"
#include"vector.h"

#define kiloByte 1024.0
#define megaByte 1048576.0

extern unsigned char applicationMode;

static unsigned int findInterval(double gaborScale, unsigned char mode)
{
	unsigned int intervalRange;

	if(mode & GAUSS_NON_GAUSS)
		intervalRange    = (unsigned int)(sqrt(-LOG_EPS_DOT_PRODUCT/M_PI)*gaborScale + 1.5);
	else if(mode & GAUSS_GAUSS)
		intervalRange  = (unsigned int)(sqrt(-LOG_EPS_DOT_PRODUCT/M_2PI)*gaborScale + 1.5);

	return intervalRange;
}


void setMP5Dimensions(Dictionary *dictionary, MP5Parameters *mp5Parameters)
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
		To sum it up, if one has gabor period, this is  at the same time a dimension of FFT at which
		one gets frequencies for optimalal dictionary.	
	*/

	// maximal scale is set to maximal scale o gabor
	double maximalScale  = *(dictionary->tableOfScalesInOptimalDictionary  + dictionary->numberOfStepsInScale - 2);
	double maximalPeriod = *(dictionary->tableOfPeriodsInOptimalDictionary + dictionary->numberOfStepsInScale - 2);
	/* Now we check, whether the FFT dimension is at lest equal to or larger then range, for which dot product has required accurancy */

	unsigned int intervalRange;
	
	intervalRange = findInterval(maximalScale,GAUSS_NON_GAUSS);

	if(mp5Parameters->accuracy & FULL)
	{
		double maximalFFTDimensionForGaborAtoms = 0.0;						 

		if(mp5Parameters->FFT)
		{
			const double factor   = (2.0*intervalRange + 1)/maximalPeriod; // the dot product is calculated in the range:  position of atom +-interval -> 2*interval + 1 points*/

			if(factor>1.0)
				maximalFFTDimensionForGaborAtoms = ((unsigned int)(factor + 1.0))*maximalPeriod;
			else
				maximalFFTDimensionForGaborAtoms = maximalPeriod;
		}
		else
			maximalFFTDimensionForGaborAtoms = 8000;//maximalFFTDimensionForGaborAtoms = 2.0*intervalRange + 1;
	
		mp5Parameters->marginalDimension       = maximalFFTDimensionForGaborAtoms/2 + 1;	
		mp5Parameters->exponensTableDimension  = mp5Parameters->offsetDimension + mp5Parameters->marginalDimension;
		mp5Parameters->offsetExpandedDimension = mp5Parameters->offsetDimension + 2*mp5Parameters->marginalDimension;
		mp5Parameters->fftTableDimension       = mp5Parameters->offsetExpandedDimension;
	
		// printf(" period: %lf %u\n",maximalPeriod,intervalRange);
		// printf(" mp5Parameters->marginalDimension: %u\n",mp5Parameters->marginalDimension);
		// printf(" mp5Parameters->exponensTableDimension:  %u\n",mp5Parameters->exponensTableDimension);
		// printf(" mp5Parameters->marginalDimension:       %u\n",mp5Parameters->marginalDimension);
		// printf(" mp5Parameters->offsetExpandedDimension: %u\n",mp5Parameters->offsetExpandedDimension);
		// printf(" mp5Parameters->fftTableDimension:       %u\n",mp5Parameters->fftTableDimension);
		// printf(" calkowity rozmiar tablicy wynosi %u \n",mp5Parameters->offsetExpandedDimension);
	}
	else if(mp5Parameters->accuracy & FLOATING)
	{
		mp5Parameters->marginalDimension       = (maximalPeriod/2 + 1) > mp5Parameters->offsetDimension ? (maximalPeriod/2 + 1) : mp5Parameters->offsetDimension;
		mp5Parameters->exponensTableDimension  = mp5Parameters->offsetDimension + mp5Parameters->marginalDimension;
		mp5Parameters->offsetExpandedDimension = mp5Parameters->offsetDimension + 2*mp5Parameters->marginalDimension;
		mp5Parameters->fftTableDimension       = mp5Parameters->offsetExpandedDimension;
	}
}

static void findIntegerScalesForDilationFactorParameters(MP5Parameters *mp5Parameters, Dictionary *dictionary)
{
	const double       numerator            = log(mp5Parameters->offsetDimension);
	const double       denominator          = log(dictionary->dilationFactor);
	unsigned short int numberOfStepsInScale = (unsigned short int)(numerator/denominator + 0.5) + 1; // the last scale is for sin/cos wave.
	unsigned short int scaleIndex = 0;
	
	if(mp5Parameters->maxGaborScale>0.0)
	{
		unsigned short int tmpNumberOfScales = 0;
		unsigned short int counter;

		for(counter=0;counter<numberOfStepsInScale;counter++)
		{
			if(pow(dictionary->dilationFactor,(counter + 1.0))>mp5Parameters->maxGaborScale)
				break;
			tmpNumberOfScales++;
		}
		numberOfStepsInScale = tmpNumberOfScales + 1;
	}
	
	dictionary->numberOfStepsInScale = numberOfStepsInScale; 
	
	dictionary->tableOfScalesInOptimalDictionary                 = (double *)dVectorAllocate(numberOfStepsInScale);
	dictionary->tableOfPeriodsInOptimalDictionary                = (unsigned int *)uiVectorAllocate(numberOfStepsInScale);
	dictionary->tableOfFrequenciesInOptimalDictionary            = (double *)dVectorAllocate(numberOfStepsInScale);
	
	dictionary->tableOfPositionsInOptimalDictionary              = (double *)dVectorAllocate(numberOfStepsInScale);

	if(dictionary->typeOfDictionary & OCTAVE_FIXED)
	{                               
		for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++)
			*(dictionary->tableOfScalesInOptimalDictionary + scaleIndex) = pow(dictionary->dilationFactor,(scaleIndex + 1.0));
		
		*(dictionary->tableOfScalesInOptimalDictionary + numberOfStepsInScale - 1) = *(dictionary->tableOfScalesInOptimalDictionary + numberOfStepsInScale - 2);
	}
	else if(dictionary->typeOfDictionary & OCTAVE_STOCH)
	{
		time_t   seedTime;
		long int seed;
		double   exponent = 0.0;
		
		if(dictionary->randomSeed==AUTO_RANDOM_SEED)
			seed = time(&seedTime);
		else
			seed = dictionary->randomSeed;

		r250_init(seed);	

		for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++)
		{
			exponent = scaleIndex+(dr250() + 0.5);	
			*(dictionary->tableOfScalesInOptimalDictionary + scaleIndex) = pow(dictionary->dilationFactor,exponent);
		}			
		
		*(dictionary->tableOfScalesInOptimalDictionary + numberOfStepsInScale - 1) = *(dictionary->tableOfScalesInOptimalDictionary + numberOfStepsInScale - 2);
	}

	// for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
		// printf(" Ala ma kota %lf \n",*(dictionary->tableOfScalesInOptimalDictionary + scaleIndex));
	// printf("\n");
	dictionary->numberOfStepsInFrequencyAtParticularScale = uiVectorAllocate((unsigned int)numberOfStepsInScale);
	dictionary->numberOfStepsInPositionAtParticularScale  = uiVectorAllocate(numberOfStepsInScale);

}

static BOOLEAN testScaleToPeriodFactorAtom(const Dictionary *dictionary, unsigned short int scaleIndex, unsigned int rifling, double scaleToPeriodFactor)
{

	const unsigned int period = (unsigned int)((*(dictionary->tableOfPeriodsInOptimalDictionary + scaleIndex))/((double)rifling) + 0.5);
	const double       scale  = (*(dictionary->tableOfScalesInOptimalDictionary  + scaleIndex))/sqrt(M_PI);

	if(scale<scaleToPeriodFactor*period)
		return FALSE;
	
	return TRUE;
}

static void atomToString(const MP5Parameters *mp5Parameters, const Dictionary *dictionary, const Atom *atom, FILE* file)
{
	fprintf(file," scale:           %lf\n",(*(dictionary->tableOfScalesInOptimalDictionary + atom->scaleIndex))/sqrt(M_PI));
	fprintf(file," position:        %u\n",atom->position);
	fprintf(file," rifling:         %u\n",atom->rifling);
	fprintf(file," frequency:       %lf\n",(atom->rifling + 1)*M_2PI/(*(dictionary->tableOfPeriodsInOptimalDictionary + atom->scaleIndex)));
	fprintf(file," basicFrequency:  %lf\n",*(dictionary->tableOfFrequenciesInOptimalDictionary + atom->scaleIndex));
	fprintf(file," basicPeriod:     %u\n",*(dictionary->tableOfPeriodsInOptimalDictionary + atom->scaleIndex));
	fprintf(file," \n");
}

static void printSizeOfDictionaryAndSizeOfSinCosExpTables(const Dictionary *dictionary, const MP5Parameters *mp5Parameters)
{
    unsigned short int scaleIndex;
    unsigned long int numberOfPointsInSinCosTables = 0;
    unsigned long int numberOfPointsInExpTables    = 0;
    unsigned long int totalNumberOfPoints          = 0;
	unsigned int sizeOfDictionary = dictionary->numberOfDiracFunctions  +
									dictionary->numberOfGaussFunctions  +
									dictionary->numberOfSinCosFunctions + 
									dictionary->numberOfCorrectGabors;

    unsigned int sizeOfRSRCPhaseTablesOfOneAtom = 3*sizeof(float)*mp5Parameters->numberOfAllocatedChannels; // 3x because we have tables OF RS RC and phase

    for(scaleIndex=0;scaleIndex<dictionary->numberOfStepsInScale;scaleIndex++)
		numberOfPointsInSinCosTables+= (*(dictionary->tableOfPeriodsInOptimalDictionary + scaleIndex));

    numberOfPointsInExpTables = mp5Parameters->offsetExpandedDimension*dictionary->numberOfStepsInScale;
    
    totalNumberOfPoints = 2*numberOfPointsInSinCosTables + numberOfPointsInExpTables;

    printf(" INFORMATION ABOUT DICTIONARY: \n");
    printf(" \n");
    printf(" NUMBER OF DIRAC'S DELTAS:  %12hu\n",dictionary->numberOfDiracFunctions);
    printf(" NUMBER OF GAUSS FUNCTIONS: %12u\n",dictionary->numberOfGaussFunctions);
    printf(" NUMBER OF SIN/COS:         %12u\n",dictionary->numberOfSinCosFunctions);
    printf(" NUMBER OF CORRECT GABORS:  %12u\n",dictionary->numberOfCorrectGabors);
    printf(" \n");
    printf(" TOTAL NUMBER OF ATOMS:     %12u      %10.3lf (MB)\n",sizeOfDictionary,(1.0*dictionary->sizeOfDictionary*(sizeof(Atom) + sizeOfRSRCPhaseTablesOfOneAtom))/megaByte);
    printf(" \n");
    printf(" SIZE OF SIN TABLES:                          %10.3lf (MB)\n",(1.0*numberOfPointsInSinCosTables*sizeof(double))/megaByte);
    printf(" SIZE OF COS TABLES:                          %10.3lf (MB)\n",(1.0*numberOfPointsInSinCosTables*sizeof(double))/megaByte);
    printf(" SIZE OF EXP TABLES:                          %10.3lf (MB)\n",(1.0*numberOfPointsInExpTables*sizeof(double))/megaByte);
    printf(" SIZE OF ALOCATED SIN/COS/EXP TABLES:         %10.3lf (MB)\n",(1.0*totalNumberOfPoints*sizeof(double))/megaByte); 
    printf(" \n");
    printf(" TOTAL MEMORY USAGE BY DICTIONARY:            %10.3lf (MB)\n",(1.0*dictionary->sizeOfDictionary*(sizeof(Atom) + sizeOfRSRCPhaseTablesOfOneAtom) + 1.0*totalNumberOfPoints*sizeof(double))/megaByte);
    fflush(stdout);
}

void analyseDictionarySizeAndType(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
	if(applicationMode & PROCESS_USER_MODE)
		printf("\n ANALYSIS OF DICTIONARY TYPE AND DICTIONARY SIZE \n");

	unsigned int offsetDimension = mp5Parameters->offsetDimension;
    unsigned int gaussCounter;
	unsigned int gaborsCounter;
    
    dictionary->basicStepInPositionInSignal  = 1.0;
    dictionary->basicStepInFrequencyInSignal = 0.0; //?; now we don't now how much
    	
    if((dictionary->typeOfDictionary & OCTAVE_FIXED) || (dictionary->typeOfDictionary & OCTAVE_STOCH))
    {

        findIntegerScalesForDilationFactorParameters(mp5Parameters,dictionary);

		const unsigned short int numberOfStepsInScale = dictionary->numberOfStepsInScale;
        double					 basicStepInFrequencyInOptimalDictionary = 0.0;
		unsigned int             basicStepInPeriodInOptimalDictionary    = 0;
		double                   basicStepInPositionInOptimalDictionary  = 0;

		double	     *tableOfScalesInOptimalDictionary      = NULL;
        unsigned int *tableOfPeriodsInOptimalDictionary     = NULL;
		double	     *tableOfFrequenciesInOptimalDictionary = NULL;
		double	     *tableOfPositionsInOptimalDictionary   = NULL;
        unsigned int numberOfStepsInFrequency = 0;
		unsigned int numberOfStepsInPosition  = 0;
	
		unsigned short int scaleIndex = 0;
		unsigned       int DT         = 0;
        double             DF         = 0.0;	
		double             DU         = 0.0;
		
		double stepInFrequency = 0.0;

        double constA = 0.0;

		gaussCounter  = 0U;
		gaborsCounter = 0U;

        constA = log(0.5*(dictionary->dilationFactor + 1.0/dictionary->dilationFactor));

		if(dictionary->typeOfDictionary & OCTAVE_FIXED)
			basicStepInPeriodInOptimalDictionary = (unsigned int)(sqrt(M_PI/constA)); // ((unsigned int)(M_2PI/(2.0*sqrt(M_PI*constA))));
                                                // basicStepInFrequencyInOptimalDictionary is in radians degree
                                                // if we want to get period in samples, we have to use the following, simple equation:
                                                // T = 2pi/w
                                                // because of saving of computer's memory,
                                                // basicStepInPeriodInOptimalDictionary should be integer.
                                                // thanks to it, we will be able to generate all sin/cos
                                                // functions from one sin/cos
												
		else if(dictionary->typeOfDictionary & OCTAVE_STOCH)
		{
			if(mp5Parameters->FFT)
				basicStepInPeriodInOptimalDictionary = ((unsigned int)(sqrt(M_PI/constA)));
			else
				basicStepInPeriodInOptimalDictionary = ((unsigned int)(sqrt(M_PI/constA)))*dictionary->periodDensity;		
		}
			
        basicStepInFrequencyInOptimalDictionary = M_2PI/basicStepInPeriodInOptimalDictionary;
   
        dictionary->basicStepInFrequencyInOptimalDictionary = basicStepInFrequencyInOptimalDictionary;
		dictionary->basicStepInPeriodInOptimalDictionary    = basicStepInPeriodInOptimalDictionary;
  
		basicStepInPositionInOptimalDictionary = sqrt(constA/M_PI);

		tableOfScalesInOptimalDictionary       = dictionary->tableOfScalesInOptimalDictionary;
        tableOfPeriodsInOptimalDictionary      = dictionary->tableOfPeriodsInOptimalDictionary;
		tableOfFrequenciesInOptimalDictionary  = dictionary->tableOfFrequenciesInOptimalDictionary;
        tableOfPositionsInOptimalDictionary    = dictionary->tableOfPositionsInOptimalDictionary;

		for(scaleIndex=0;scaleIndex<(numberOfStepsInScale - 1);scaleIndex++)
		{				
			DT = *(tableOfPeriodsInOptimalDictionary     + scaleIndex) = (unsigned int)(basicStepInPeriodInOptimalDictionary*(*(tableOfScalesInOptimalDictionary + scaleIndex)) + 0.5);
			DF = *(tableOfFrequenciesInOptimalDictionary + scaleIndex) = M_2PI/DT;
			DU = (*(tableOfScalesInOptimalDictionary + scaleIndex))*basicStepInPositionInOptimalDictionary;
				
			 // printf(" %u %lf\n",scaleIndex,*(tableOfScalesInOptimalDictionary + scaleIndex));
			 // printf(" DT %u\n",DT);
			 // printf(" DF %lf\n",DF);
			 // printf(" DU %lf \n\n",DU);

			fflush(stdout);
				
			if(DU<=(double)(dictionary->basicStepInPositionInSignal))
				DU = (double)(dictionary->basicStepInPositionInSignal);
				
			*(tableOfPositionsInOptimalDictionary + scaleIndex) = DU;
	    
			if(dictionary->typeOfDictionary & OCTAVE_STOCH)
			{
				if(mp5Parameters->FFT)
					stepInFrequency = DF;
				else
					stepInFrequency = DF*dictionary->periodDensity;
			}
			else
				stepInFrequency = DF;
				
			numberOfStepsInFrequency = (unsigned int)((M_PI - 1.5*stepInFrequency)/stepInFrequency + 1);
			*(dictionary->numberOfStepsInFrequencyAtParticularScale + scaleIndex) = numberOfStepsInFrequency;

			numberOfStepsInPosition  = (unsigned int)((offsetDimension - 1.5*DU)/DU + 1);		
			*(dictionary->numberOfStepsInPositionAtParticularScale + scaleIndex) = numberOfStepsInPosition;
			
			if(dictionary->gaussInDictionary)
				gaussCounter +=numberOfStepsInPosition;
				
			gaborsCounter+=numberOfStepsInFrequency*numberOfStepsInPosition;
		
			if(scaleIndex==(numberOfStepsInScale-2)) // at this moment we know the longest period (it is for highest scale), so we have to set up
			{
				/* find the size of tables, which are required by algorithms 
				the size of tables depends on maximal scale of atoms */			
				setMP5Dimensions(dictionary,mp5Parameters);
			}
		}

		if(dictionary->diracInDictionary)
			dictionary->numberOfDiracFunctions = offsetDimension;

		dictionary->numberOfGaussFunctions = gaussCounter;
	
		if(dictionary->sinCosInDictionary)
		{						 
			DT = *(tableOfPeriodsInOptimalDictionary      + numberOfStepsInScale - 1) = *(tableOfPeriodsInOptimalDictionary      + numberOfStepsInScale - 2);
 			DF = *(tableOfFrequenciesInOptimalDictionary  + numberOfStepsInScale - 1) = M_2PI/DT;
									    			
			numberOfStepsInFrequency = *(dictionary->numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 2);
			*(dictionary->numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1) = numberOfStepsInFrequency;
		
			dictionary->numberOfSinCosFunctions = numberOfStepsInFrequency;
		}
				
		dictionary->sizeOfDictionary = dictionary->numberOfDiracFunctions  + /* number of Dirac Delta          */
									   dictionary->numberOfGaussFunctions  + /* number of  Gauss Functions  */
									   dictionary->numberOfSinCosFunctions + /* number of  Sin/Cos Waves   */
									   gaborsCounter;                        /* number of Gabor Waves */
	
		dictionary->numberOfNonFFTAtoms = dictionary->numberOfDiracFunctions + dictionary->numberOfGaussFunctions;
	}

}

void allocateDictionary(Dictionary *dictionary, const MP5Parameters *mp5Parameters)
{
    unsigned int             atomsCounter = 0;
    const unsigned short int numberOfAllocatedChannels = mp5Parameters->numberOfAllocatedChannels;
    Atom *atomPointer = NULL;

    dictionary->atomsTable = (Atom *)malloc(dictionary->sizeOfDictionary*sizeof(Atom));

    for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
    {
		atomPointer = (dictionary->atomsTable + atomsCounter);
		allocateAtomElements(atomPointer,numberOfAllocatedChannels, mp5Parameters->MPType);
    }
}

void freeDictionary(Dictionary *dictionary)
{
    unsigned int atom;

	if(dictionary->tableOfScalesInOptimalDictionary!=NULL)
		dVectorFree(dictionary->tableOfScalesInOptimalDictionary);
	if(dictionary->tableOfPeriodsInOptimalDictionary!=NULL)
		uiVectorFree(dictionary->tableOfPeriodsInOptimalDictionary);
	if(dictionary->tableOfFrequenciesInOptimalDictionary!=NULL)
		dVectorFree(dictionary->tableOfFrequenciesInOptimalDictionary);
	if(dictionary->tableOfPositionsInOptimalDictionary!=NULL)
		dVectorFree(dictionary->tableOfPositionsInOptimalDictionary);    
	if(dictionary->numberOfStepsInFrequencyAtParticularScale!=NULL)
		uiVectorFree(dictionary->numberOfStepsInFrequencyAtParticularScale);
	if(dictionary->numberOfStepsInPositionAtParticularScale!=NULL)
		uiVectorFree(dictionary->numberOfStepsInPositionAtParticularScale);

	if(dictionary->atomsTable!=NULL)
	{
		for(atom=0;atom<dictionary->sizeOfDictionary;atom++)
			if((dictionary->atomsTable + atom)!=NULL)
				freeAtomElements(dictionary->atomsTable + atom);
				
		free(dictionary->atomsTable);
	}
	
}

void makeDictionary(Dictionary *dictionary, const MP5Parameters *mp5Parameters)
{	
	unsigned int numberOfCorrectGabors = 0;
	
	time_t       seedTime = 0;
    long int     seed     = 0;
    Atom        *atomPointer = NULL;
	
	double scaleToPeriodFactor = dictionary->scaleToPeriodFactor;
    const unsigned int offsetDimension = mp5Parameters->offsetDimension;
          unsigned int leftSize = 0, rightSize = 0;

    const unsigned short int numberOfStepsInScale = dictionary->numberOfStepsInScale;

    const double *tableOfFrequenciesInOptimalDictionary  = dictionary->tableOfFrequenciesInOptimalDictionary;
    const double *tableOfPositionsInOptimalDictionary    = dictionary->tableOfPositionsInOptimalDictionary;

    const unsigned int *numberOfStepsInFrequencyAtParticularScale = dictionary->numberOfStepsInFrequencyAtParticularScale;
    const unsigned int *numberOfStepsInPositionAtParticularScale  = dictionary->numberOfStepsInPositionAtParticularScale;

    const unsigned int periodDensity = dictionary->periodDensity;

    unsigned int numberOfStepsInFrequency = 0;
    unsigned int numberOfStepsInPosition  = 0;
    
    unsigned short int scaleIndex     = 0;
    unsigned       int frequencyIndex = 0;
    unsigned       int positionIndex  = 0;
    unsigned       int rifling        = 0;
    double             frequency      = 0.0;

    double DF              = 0.0;
    double DU              = 0.0;
    double stepInFrequency = 0.0; 
    double        dr250Position  = 0.0;
    short  int    dn250Frequency = 0;
    double        dr250Frequency = 0.0;
	
    seed = time(&seedTime);

    r250_init(seed);	

    atomPointer = dictionary->atomsTable;
		
	if(applicationMode & PROCESS_USER_MODE)
		printf("\n START DICTIONARY GENERAITING \n\n");
		
    /* set Dirac Delta */
	if(dictionary->diracInDictionary)
	{
		for(positionIndex=0;positionIndex<(double)(offsetDimension);positionIndex++)
		{			
			atomPointer->scaleIndex  = 0; // although scaleIndex = 0 means scale = 1, we set this value, 
										  // because this fild cand not be empty. However in case of Dirac Delta 
									      // this fild is not used, so it will not cause any disturbance 

			atomPointer->rifling     = 0;
			atomPointer->position    = (unsigned int)positionIndex;
			atomPointer->feature|=DIRACDELTA;

			leftSize  = atomPointer->position;
			rightSize = (unsigned int)(offsetDimension - atomPointer->position - 1);

			if(leftSize<=rightSize)
				atomPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;
		
			atomPointer++;
		} 
	}
	
	if(dictionary->gaussInDictionary)
	{
		/* set Gauss Functions */
		for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++) // scaleindex form 0 to numberOfStepsInScale-2 is for gauss/gabor function
		{																   // scale for scaleIndex= numberOfStepsInScale-1 us for sin/cos wave
			DU = *(tableOfPositionsInOptimalDictionary   + scaleIndex);
			numberOfStepsInPosition  = *(numberOfStepsInPositionAtParticularScale + scaleIndex);	
			
			if(dictionary->typeOfDictionary & OCTAVE_FIXED)
			{        
				for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
				{
					atomPointer->scaleIndex = scaleIndex;
					atomPointer->position   = (unsigned int)(DU*(positionIndex + 1));	

					atomPointer->feature|=GAUSSFUNCTION;
					leftSize  = atomPointer->position;
					rightSize = (unsigned int)(offsetDimension - atomPointer->position - 1);
					
					if(leftSize<=rightSize)
						atomPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;

					atomPointer++;
				}
			}
			else if(dictionary->typeOfDictionary & OCTAVE_STOCH)
			{		
				for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
				{
					atomPointer->scaleIndex = scaleIndex;
					atomPointer->position   = (unsigned int)(DU*(positionIndex + 1) + DU*(dr250() - 0.5) + 0.5);	

					atomPointer->feature|=GAUSSFUNCTION;
					leftSize  = atomPointer->position;
					rightSize = (unsigned int)(offsetDimension - atomPointer->position - 1); 

					if(leftSize<=rightSize)
						atomPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;
					atomPointer++;
				}
			}
		}
	}

	/* set Sin/Cos Waves */
	if(dictionary->sinCosInDictionary)
	{
		DF = *(tableOfFrequenciesInOptimalDictionary + numberOfStepsInScale - 1);

		if(dictionary->typeOfDictionary & OCTAVE_STOCH)
		{
			if(mp5Parameters->FFT)
				stepInFrequency = DF;
			else
				stepInFrequency = DF*periodDensity;			
		}
		else
			stepInFrequency = DF;
			
		numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);	

		if(dictionary->typeOfDictionary & OCTAVE_STOCH)
		{
			if(mp5Parameters->FFT)
				dr250Frequency = DF*(r250() - 0.5);
			else
				dn250Frequency = r250n(periodDensity) - periodDensity/2;
		}
		
		for(frequencyIndex=0,rifling=1;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++,rifling=rifling+periodDensity)
		{
			atomPointer->scaleIndex = (unsigned short int)(numberOfStepsInScale - 1);  // although scaleIndex = numberOfScale - 1 means scale = offsetDimension, we set this value, 
																					   // because this fild can not be empty. However in case of FFT WAVE 
																				       // this fild is not used, so it will not cause any disturbance 
					
			if(dictionary->typeOfDictionary & OCTAVE_FIXED)
			{            
				frequency               = stepInFrequency*(frequencyIndex + 1);
				atomPointer->rifling    = rifling;
				atomPointer->randomShiftInFrequency = 0;

				if((offsetDimension%2)==0)
					atomPointer->position = (offsetDimension - 1)/2;
				else
					atomPointer->position = offsetDimension/2;

				atomPointer->feature|=SINCOSWAVE;

				leftSize  = atomPointer->position;
				rightSize = (unsigned int)(offsetDimension - atomPointer->position - 1);

				if(leftSize<=rightSize)
					atomPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;

				atomPointer++;
			}
			else if(dictionary->typeOfDictionary & OCTAVE_STOCH)
			{
				if(mp5Parameters->FFT)
					frequency = stepInFrequency*(frequencyIndex + 1);
				else
					frequency = stepInFrequency*(frequencyIndex + 1) + DF*(r250n(periodDensity) - periodDensity/2);					
                 
				atomPointer->rifling    = (unsigned int)(frequency/DF + 0.5);
				
				if(mp5Parameters->FFT)
					atomPointer->randomShiftInFrequency = dr250Frequency;
				else
					atomPointer->randomShiftInFrequency = 0;
				
				if((offsetDimension%2)==0)
					atomPointer->position = (offsetDimension - 1)/2;
				else
					atomPointer->position = offsetDimension/2;
				
				atomPointer->feature|=SINCOSWAVE;

				leftSize  = atomPointer->position;
				rightSize = (unsigned int)(offsetDimension - atomPointer->position - 1); 

				if(leftSize<=rightSize)
					atomPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;

				atomPointer++;		
			}
		}
	}
	
	/* set Gabor waves */
    for(scaleIndex=0;scaleIndex<(numberOfStepsInScale-1);scaleIndex++)
	{
		DF = *(tableOfFrequenciesInOptimalDictionary + scaleIndex);
		DU = *(tableOfPositionsInOptimalDictionary   + scaleIndex);

		if(dictionary->typeOfDictionary & OCTAVE_STOCH)
		{
			if(mp5Parameters->FFT)
				stepInFrequency = DF;
			else
				stepInFrequency = DF*periodDensity;			
		}
		else
			stepInFrequency = DF;
			
		numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);	
		numberOfStepsInPosition  = *(numberOfStepsInPositionAtParticularScale + scaleIndex);	

		if(dictionary->typeOfDictionary & OCTAVE_FIXED)
		{            
			for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
    	    {
				for(frequencyIndex=0,rifling=1;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++,rifling=rifling+1)
				{				
					atomPointer->scaleIndex = scaleIndex;
					frequency               = stepInFrequency*(frequencyIndex + 1);
					atomPointer->rifling    = rifling;
					atomPointer->randomShiftInFrequency = 0;
					atomPointer->position   = (unsigned int)(DU*(positionIndex + 1));	

					atomPointer->feature|=GABORWAVE;
					leftSize  = atomPointer->position;
					rightSize = (unsigned int)(offsetDimension - atomPointer->position - 1); 

					if(leftSize<=rightSize)
						atomPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;

					if(testScaleToPeriodFactorAtom(dictionary,scaleIndex,rifling,scaleToPeriodFactor))
						numberOfCorrectGabors++;
					else
						atomPointer->feature|=INCORRECTGABOR;
					
					atomPointer++;
				}  	    
			}
		}
		else if(dictionary->typeOfDictionary & OCTAVE_STOCH)
		{		
			for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
			{
				dr250Position = dr250();
				
				if(mp5Parameters->FFT)
					dr250Frequency = DF*(dr250() - 0.5);
				else
					dn250Frequency = r250n(periodDensity) - periodDensity/2;
				
				for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
				{
					if(mp5Parameters->FFT)
						frequency = stepInFrequency*(frequencyIndex + 1);                     
					else
						frequency = stepInFrequency*(frequencyIndex + 1) + DF*(r250n(periodDensity) - periodDensity/2);					
		
					atomPointer->rifling = (unsigned int)(frequency/DF + 0.5);
		
					if(mp5Parameters->FFT)
						atomPointer->randomShiftInFrequency = dr250Frequency;
					else
						atomPointer->randomShiftInFrequency = 0;
										
					atomPointer->scaleIndex = scaleIndex;
										
					if(mp5Parameters->FFT)
						atomPointer->position   = (unsigned int)(DU*(positionIndex + 1) + DU*(dr250Position -0.5) + 0.5);
					else
						atomPointer->position   = (unsigned int)(DU*(positionIndex + 1) + DU*(dr250() - 0.5) + 0.5);	
						
					atomPointer->feature|=GABORWAVE;
					leftSize  = atomPointer->position;
					rightSize = (unsigned int)(offsetDimension - atomPointer->position - 1);

					if(leftSize<=rightSize)
						atomPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;
						
					if(testScaleToPeriodFactorAtom(dictionary,scaleIndex,rifling,scaleToPeriodFactor))
						numberOfCorrectGabors++;
					else
						atomPointer->feature|=INCORRECTGABOR;
	
					atomPointer++;
				}
			}
		}
    }
	
	/*atomPointer = dictionary->atomsTable;
	unsigned int atomsCounter;
	for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
	{
		scaleIndex = atomPointer->scaleIndex;
		
		printf("%u %hu %lf %hu %u %hu %hu %lf %lf\n",atomsCounter,
												 scaleIndex,
												 *(dictionary->tableOfScalesInOptimalDictionary + atomPointer->scaleIndex),
												 atomPointer->position,
												 atomPointer->rifling,
												 (atomPointer->feature & LEFT_SIDE_POSITION_IN_OFFSET),
												 (atomPointer->feature & DIRACDELTA)    |
												 (atomPointer->feature & GAUSSFUNCTION) |
												 (atomPointer->feature & SINCOSWAVE)    |
												 (atomPointer->feature & GABORWAVE),			 
												 *(tableOfFrequenciesInOptimalDictionary + scaleIndex),
												 *(tableOfPositionsInOptimalDictionary   + scaleIndex));
		
		atomPointer++;
	}
	getchar();*/
	
	dictionary->numberOfCorrectGabors = numberOfCorrectGabors;

	if(applicationMode & PROCESS_USER_MODE)
		printSizeOfDictionaryAndSizeOfSinCosExpTables(dictionary,mp5Parameters);	

	atomPointer = dictionary->atomsTable;
	
	#ifdef DICTIONARY_DUMPING
		unsigned int atomNumber;
		FILE *file = fopen("dictionary.dump","wt");
		for(atomNumber=0;atomNumber<dictionary->sizeOfDictionary;atomNumber++)
		{
			atomToString(mp5Parameters,dictionary,atomPointer,file);
			atomPointer++;
		}
		fclose(file);
	#endif
	
}

void reinitDictionary(Dictionary *dictionary, const MP5Parameters *mp5Parameters)
{
    unsigned int atomsCounter;
    
    for(atomsCounter=0;atomsCounter<dictionary->sizeOfDictionary;atomsCounter++)
		(dictionary->atomsTable + atomsCounter)->feature = 0x0;

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
    unsigned int atomsCounter;
    const unsigned int sizeOfDictionary = dictionary->sizeOfDictionary;

    for(atomsCounter=0;atomsCounter<sizeOfDictionary;atomsCounter++)
		(dictionary->atomsTable + atomsCounter)->feature&=~GABOR_WAS_HIT;

}
