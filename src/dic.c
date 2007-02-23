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

#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"include/def.h"
#include"include/dic.h"
#include"include/gabor.h"
#include"include/r250.h"
#include"include/types.h"
#include"include/vector.h"

#define kiloByte 1024.0
#define megaByte 1048576.0


static void dictionaryToString(const GaborDictionary *gaborDictionary)
{
        unsigned short int scaleIndex;
 		
	printf("\n");
 	printf(" GaborDictionary\n");
        printf(" sizeOfDictionary:                        %u\n",gaborDictionary->sizeOfDictionary);
	if(gaborDictionary->typeOfDictionary & OCTAVE_FIXED)
	    printf(" typeOfDictionary:                        OCTAVE_FIXED\n");       
        else if(gaborDictionary->typeOfDictionary & OCTAVE_STOCH)
	    printf(" typeOfDictionary:                        OCTAVE_STOCH\n");  
                       
        printf(" scaleToPeriodFactor:                     %lf\n",gaborDictionary->scaleToPeriodFactor);
        printf(" dilationFactor:                          %lf\n",gaborDictionary->dilationFactor);
        printf(" numberOfstepsInScale                     %hu\n",gaborDictionary->numberOfStepsInScale);
        printf(" basicStepInPositionInSignal:             %lf\n",gaborDictionary->basicStepInPositionInSignal);
        printf(" basicStepInFrequencyInSignal:            %lf\n",gaborDictionary->basicStepInFrequencyInSignal);
        printf(" basicStepInFrequencyInOptimalDictionary: %lf\n",gaborDictionary->basicStepInFrequencyInOptimalDictionary);
        printf(" basicStepInPeriodInOptimalDictionary:    %u\n",gaborDictionary->basicStepInPeriodInOptimalDictionary);
        printf(" basicStepInPositionInOptimalDictionary:  %lf\n",gaborDictionary->basicStepInPositionInOptimalDictionary);
                
        for(scaleIndex=0;scaleIndex<gaborDictionary->numberOfStepsInScale;scaleIndex++)
        {
              printf(" [%hu]                                         \
                       tableOfScalesInOptimalDictionary:      %5hu   \
                       tableOfPeriodsInOptimalDictionary:     %12u     \
                       tableOfFrequenciesInOptimalDictionary: %7.6lf \
                       tableOfPositionsInOptimalDictionary:   %11.6lf \n", 
                       scaleIndex,
                       *(gaborDictionary->tableOfScalesInOptimalDictionary      + scaleIndex),     
                       *(gaborDictionary->tableOfPeriodsInOptimalDictionary     + scaleIndex), 
                       *(gaborDictionary->tableOfFrequenciesInOptimalDictionary + scaleIndex), 
                       *(gaborDictionary->tableOfPositionsInOptimalDictionary   + scaleIndex)); 
		       fflush(stdout);
        }
}

static void findIntegerScalesForDilationFactorParameters(const MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary)
{
	const unsigned short int dimOffset = mp5Parameters->dimOffset;
	const          double numerator    =  log(mp5Parameters->dimOffset);
	const          double denominator  =  log(gaborDictionary->dilationFactor);
	const unsigned short int numberOfStepsInScale = (unsigned short int)floor(numerator/denominator);
	unsigned short int numberOfDifferentScales = 1;
	unsigned short int scaleIndex,differentScalesCounter;
	unsigned short int currentScale  = 1;
	unsigned short int previousScale = 1;

	for(scaleIndex=1;scaleIndex<=numberOfStepsInScale;scaleIndex++)
	{
		currentScale = (unsigned short int)(pow(gaborDictionary->dilationFactor,(double)scaleIndex));
		
		if(currentScale==previousScale)
			continue;
		else
        	{
			numberOfDifferentScales++;
			previousScale=currentScale;
		}
	}

	if(currentScale!=dimOffset)
		numberOfDifferentScales++;

	gaborDictionary->numberOfStepsInScale = numberOfDifferentScales;
	
	gaborDictionary->tableOfScalesInOptimalDictionary      = (unsigned short int *)usiVectorAllocate(numberOfDifferentScales);
	gaborDictionary->tableOfPeriodsInOptimalDictionary     = (unsigned int *)uiVectorAllocate(numberOfDifferentScales);
	gaborDictionary->tableOfFrequenciesInOptimalDictionary = (double *)dVectorAllocate(numberOfDifferentScales);
	gaborDictionary->tableOfPositionsInOptimalDictionary   = (double *)dVectorAllocate(numberOfDifferentScales);

	if(gaborDictionary->typeOfDictionary & OCTAVE_FIXED)
	{
		currentScale  = 1;
		previousScale = 1;
		*(gaborDictionary->tableOfScalesInOptimalDictionary) = currentScale;
                                      
		differentScalesCounter = 1;
       
		for(scaleIndex=1;scaleIndex<=numberOfStepsInScale;scaleIndex++)
		{
			currentScale = (unsigned short int)(pow(gaborDictionary->dilationFactor,(double)scaleIndex));
	   
			if(currentScale==previousScale)
				continue;
			else
			{
				*(gaborDictionary->tableOfScalesInOptimalDictionary + differentScalesCounter) = currentScale;

				differentScalesCounter++;
				previousScale=currentScale;
			}
		}             
    
		*(gaborDictionary->tableOfScalesInOptimalDictionary + numberOfDifferentScales - 1) = dimOffset;

	}
	else if(gaborDictionary->typeOfDictionary & OCTAVE_STOCH)
	{
		time_t       seedTime;
		long int     seed;
		seed = time(&seedTime);

		r250_init(seed);	

		unsigned short int tmpTableOfPeriodsInOptimalDictionary[numberOfDifferentScales];
		unsigned short int firstScale  = 0;
		unsigned short int secondScale = 0;

		currentScale  = 1;
		previousScale = 1;
		*tmpTableOfPeriodsInOptimalDictionary = currentScale;
                                      
		differentScalesCounter = 1;
       
		for(scaleIndex=1;scaleIndex<=numberOfStepsInScale;scaleIndex++)
		{
			currentScale = (unsigned short int)(pow(gaborDictionary->dilationFactor,(double)scaleIndex));
	   
			if(currentScale==previousScale)
				continue;
			else
			{
				*(tmpTableOfPeriodsInOptimalDictionary + differentScalesCounter) = currentScale;

				differentScalesCounter++;
				previousScale=currentScale;
			}
		}             
    
		*(tmpTableOfPeriodsInOptimalDictionary + numberOfDifferentScales - 1) = dimOffset;

		for(scaleIndex=0;scaleIndex<numberOfDifferentScales-1;scaleIndex++)
		{
			firstScale  = *(tmpTableOfPeriodsInOptimalDictionary + scaleIndex);
			secondScale = *(tmpTableOfPeriodsInOptimalDictionary + scaleIndex + 1);
			*(gaborDictionary->tableOfScalesInOptimalDictionary + scaleIndex) = (unsigned short int)(firstScale + r250n((unsigned int)(secondScale - firstScale)));
		}
		
		*(gaborDictionary->tableOfScalesInOptimalDictionary + numberOfDifferentScales - 1) = dimOffset;
	
	}

	gaborDictionary->numberOfStepsInFrequencyAtParticularScale = uiVectorAllocate((unsigned int)numberOfDifferentScales);
	gaborDictionary->numberOfStepsInPositionAtParticularScale  = usiVectorAllocate(numberOfDifferentScales);
	
	gaborDictionary->memoryAllocated|=FIRST_PART_OF_TABLES_IN_DICTIONARY_ALLOCATED;

}

static void testScaleToPeriodFactorGabor(const GaborDictionary *gaborDictionary, Gabor *gabor, double scaleToPeriodFactor)
{
    const unsigned int period = (unsigned int)((*(gaborDictionary->tableOfPeriodsInOptimalDictionary +  gabor->scaleIndex))/((double)gabor->rifling) + 0.5);
    const double       scale  = (*(gaborDictionary->tableOfScalesInOptimalDictionary  +  gabor->scaleIndex))/sqrt(M_PI);

    if(scale<scaleToPeriodFactor*period)
	gabor->feature|=INCORRECTGABOR;
}


void printSizeOfDictionaryAndSizeOfSinCosExpTables(const MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary)
{
    unsigned short int scaleIndex;

    unsigned long int numberOfPointsInSinCosTables = 0;
    unsigned long int numberOfPointsInExpTables    = 0;
    unsigned long int totalNumberOfPoints          = 0;

    unsigned int sizeOfRSRCPhaseTablesOfOneGabor = 3*sizeof(float)*mp5Parameters->numberOfAnalysedChannels; // 3x because we have tables OF RS RC and phase

    for(scaleIndex=0;scaleIndex<gaborDictionary->numberOfStepsInScale;scaleIndex++)
	numberOfPointsInSinCosTables+= (*(gaborDictionary->tableOfPeriodsInOptimalDictionary + scaleIndex));

    numberOfPointsInExpTables = mp5Parameters->dimExpTable*gaborDictionary->numberOfStepsInScale;
    
    totalNumberOfPoints = 2*numberOfPointsInSinCosTables + numberOfPointsInExpTables;

    printf(" INFORMATION ABOUT DICTIONARY: \n");	
    printf(" \n");
    printf(" NUMBER OF DIRAC'S DELTAS:  %12hu\n",mp5Parameters->dimOffset);
    printf(" NUMBER OF GABORS:          %12u\n",gaborDictionary->sizeOfDictionary - mp5Parameters->dimOffset - gaborDictionary->numberSinCosFunctions + 1);
    printf(" NUMBER OF SIN/COS:         %12u\n",gaborDictionary->numberSinCosFunctions);
    printf(" \n");
    printf(" TOTAL NUMBER OF GABORS:    %12u      %10.3lf (MB)\n",gaborDictionary->sizeOfDictionary,(1.0*gaborDictionary->sizeOfDictionary*(sizeof(Gabor) + sizeOfRSRCPhaseTablesOfOneGabor))/megaByte);
    printf(" \n");
    printf(" SIZE OF SIN TABLES:                          %10.3lf (MB)\n",(1.0*numberOfPointsInSinCosTables*sizeof(double))/megaByte);
    printf(" SIZE OF COS TABLES:                          %10.3lf (MB)\n",(1.0*numberOfPointsInSinCosTables*sizeof(double))/megaByte);
    printf(" SIZE OF EXP TABLES:                          %10.3lf (MB)\n",(1.0*numberOfPointsInExpTables*sizeof(double))/megaByte);
    printf(" SIZE OF ALOCATED SIN/COS/EXP TABLES:         %10.3lf (MB)\n",(1.0*totalNumberOfPoints*sizeof(double))/megaByte); 
    printf(" \n");
    printf(" TOTAL MEMORY USAGE BY DICTIONARY:            %10.3lf (MB)\n",(1.0*gaborDictionary->sizeOfDictionary*(sizeof(Gabor) + sizeOfRSRCPhaseTablesOfOneGabor) + 1.0*totalNumberOfPoints*sizeof(double))/megaByte);
    fflush(stdout);
    
}


void analyseDictionarySizeAndType(MP5Parameters *mp5Parameters,GaborDictionary *gaborDictionary)
{
    printf("\n ANALYSIS OF DICTIONARY TYPE AND DICTIONARY SIZE \n");

    unsigned short int dimOffset = mp5Parameters->dimOffset;
    unsigned int gaborsCounter;

    gaborDictionary->basicStepInPositionInSignal  = 1.0;
    gaborDictionary->basicStepInFrequencyInSignal = 0.0; //?; now we don now how much
    
    if((gaborDictionary->typeOfDictionary & OCTAVE_FIXED) || (gaborDictionary->typeOfDictionary & OCTAVE_STOCH))
    {

        findIntegerScalesForDilationFactorParameters(mp5Parameters,gaborDictionary);

	const unsigned short int numberOfStepsInScale = gaborDictionary->numberOfStepsInScale;
        double       		 basicStepInFrequencyInOptimalDictionary;
	unsigned int basicStepInPeriodInOptimalDictionary = 0;

	unsigned short int *tableOfScalesInOptimalDictionary;
        unsigned       int *tableOfPeriodsInOptimalDictionary;
	double		   *tableOfFrequenciesInOptimalDictionary;
	double		   *tableOfPositionsInOptimalDictionary;
        unsigned       int numberOfStepsInFrequency = 0;
	unsigned short int numberOfStepsInPosition  = 0;
	
	unsigned short int scaleIndex;
	unsigned       int DT;
        double             DF;	
	double             DU;

	double stepInFrequency = 0.0;

        double constA;

	gaborsCounter = 0U;

        constA = log(0.5*(gaborDictionary->dilationFactor + 1.0/gaborDictionary->dilationFactor));

	if(gaborDictionary->typeOfDictionary & OCTAVE_FIXED)
	    basicStepInPeriodInOptimalDictionary = ((unsigned int)ceil(M_2PI/(2.0*sqrt(M_PI*constA))))*gaborDictionary->periodDensity;
                                                // basicStepInFrequencyInOptimalDictionary is in radians degree
                                                // if we want to get period in samples, we have to use the following, simple equation:
                                                // T = 2pi/w
                                                // because of saving of computer's memory,
                                                // basicStepInPeriodInOptimalDictionary should be integer.
                                                // thanks to it, we will be able to generate all sin/cos
                                                // functions from one sin/cos
	else if(gaborDictionary->typeOfDictionary & OCTAVE_STOCH)
	    basicStepInPeriodInOptimalDictionary = 2*((unsigned int)ceil(M_2PI/(2.0*sqrt(M_PI*constA))))*gaborDictionary->periodDensity;
						// in comparizon with equation above we multpily basicStepInPeriodInOptimalDictionary by 2
						// because the range between randomed gabors must be DF.
						// In case of OCTAVE_FIXED Dictionary gabors in frequency domain are positioned every DF
						// Now we want random gabor. To avoid situation when the range betwean two nearest gabors is
						// larger then DF, estimated from equations, we have divided it by 2, this mean we
						// multiply period by 2

	if(basicStepInPeriodInOptimalDictionary<4)
	    basicStepInPeriodInOptimalDictionary = 4;
	
	// Why minimal step in period is set to 4?
	// The smallest period in digital signal is 2.
	// If someone malicious set dilationFactor very small, he would get very large basicStepInPeriodInOptimalDictionary,
	// but on the other hand, if he set dilationFactor very large, he colud get basicStepInPeriodInOptimalDictionary smallest
	// then one sample !
	// to avoid such a situation, we have to give some minmal value. But why not 2?
	// Other periods are estimated as a multiplication of basicStepInPeriodInOptimalDictionary in a such way,
	// that frequencies related with these periods are not higher then PI - 2*PI/basicStepInPeriodInOptimalDictionary.
	// After some short modification of equeations one can find, that in order to get frequencies consisted
	// with the condition above, basicStepInPeriodInOptimalDictionary must be at least 4.
	    
	
	// basicStepInFrequencyInOptimalDictionary might have been non integer. 
	// We have found integer basicStepInPeriodInOptimalDictionary, and now we calculate
	// basicStepInFrequencyInOptimalDictionary for it.
  
        basicStepInFrequencyInOptimalDictionary = M_2PI/basicStepInPeriodInOptimalDictionary;
   
        gaborDictionary->basicStepInFrequencyInOptimalDictionary = basicStepInFrequencyInOptimalDictionary;
	gaborDictionary->basicStepInPeriodInOptimalDictionary    = basicStepInPeriodInOptimalDictionary;
  
	if(gaborDictionary->typeOfDictionary & OCTAVE_FIXED)
    	    gaborDictionary->basicStepInPositionInOptimalDictionary = sqrt(constA/M_PI);
	else if(gaborDictionary->typeOfDictionary & OCTAVE_STOCH)
    	    gaborDictionary->basicStepInPositionInOptimalDictionary = sqrt(constA/M_PI)/2.0; // why we do this, it is described above
	    
	if(gaborDictionary->basicStepInPositionInOptimalDictionary<=(double)(gaborDictionary->basicStepInPositionInSignal))
	    gaborDictionary->basicStepInPositionInOptimalDictionary = (double)(gaborDictionary->basicStepInPositionInSignal);
      
	tableOfScalesInOptimalDictionary      = gaborDictionary->tableOfScalesInOptimalDictionary;
        tableOfPeriodsInOptimalDictionary     = gaborDictionary->tableOfPeriodsInOptimalDictionary;
	tableOfFrequenciesInOptimalDictionary = gaborDictionary->tableOfFrequenciesInOptimalDictionary;
        tableOfPositionsInOptimalDictionary   = gaborDictionary->tableOfPositionsInOptimalDictionary;

        for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
	{
	    DT = *(tableOfPeriodsInOptimalDictionary     + scaleIndex) = basicStepInPeriodInOptimalDictionary*(*(tableOfScalesInOptimalDictionary + scaleIndex));
	    DF = *(tableOfFrequenciesInOptimalDictionary + scaleIndex) = M_2PI/DT;
	    DU = (*(tableOfScalesInOptimalDictionary     + scaleIndex))*sqrt(constA/M_PI);

	    if(DU<=gaborDictionary->basicStepInPositionInSignal)
		DU = gaborDictionary->basicStepInPositionInOptimalDictionary; // this means that DU = gaborDictionary->basicStepInPositionInOptimalDictionary = gaborDictionary->basicStepInPositionInSignal
         
	    *(tableOfPositionsInOptimalDictionary + scaleIndex) = DU;
	    
	    stepInFrequency = DF*gaborDictionary->periodDensity;

	    if(gaborDictionary->typeOfDictionary & OCTAVE_FIXED)
	    {
		numberOfStepsInFrequency = (unsigned int)((M_PI - 2.0*stepInFrequency)/stepInFrequency + 1);
		*(gaborDictionary->numberOfStepsInFrequencyAtParticularScale + scaleIndex) = numberOfStepsInFrequency;

		numberOfStepsInPosition  = (unsigned short int)((dimOffset - 2.0*DU)/DU + 1);		
		*(gaborDictionary->numberOfStepsInPositionAtParticularScale + scaleIndex)  = numberOfStepsInPosition;
	    }
	    else if(gaborDictionary->typeOfDictionary & OCTAVE_STOCH)
	    {
		numberOfStepsInFrequency = (unsigned int)((M_PI - 3.0*stepInFrequency)/stepInFrequency + 1);
		*(gaborDictionary->numberOfStepsInFrequencyAtParticularScale + scaleIndex) = numberOfStepsInFrequency;

		numberOfStepsInPosition  = (unsigned short int)((dimOffset - 3.0*DU)/DU + 1);		
		*(gaborDictionary->numberOfStepsInPositionAtParticularScale + scaleIndex)  = numberOfStepsInPosition;	    
	    }
	    
	    if(scaleIndex==(numberOfStepsInScale-1))
		gaborDictionary->numberSinCosFunctions = numberOfStepsInFrequency;

	    gaborsCounter+=numberOfStepsInFrequency*numberOfStepsInPosition;	
	}
        
	gaborDictionary->sizeOfDictionary = dimOffset     +                              /* number of Dirac Delta    */
					    gaborsCounter +                              /* number of Gabors Waves   */
					    gaborDictionary->numberSinCosFunctions;      /* number of  Sin/Cos Waves */

    }

}

void allocateDictionary(const MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary)
{
    unsigned int gaborsCounter;
    const unsigned short int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;
    Gabor *gaborPointer;

    gaborDictionary->gaborsTable = (Gabor *)malloc(gaborDictionary->sizeOfDictionary*sizeof(Gabor));

    for(gaborsCounter=0;gaborsCounter<gaborDictionary->sizeOfDictionary;gaborsCounter++)
    {
	gaborPointer = (gaborDictionary->gaborsTable + gaborsCounter);
	allocateGaborElements(gaborPointer,numberOfAnalysedChannels);
    }

    gaborDictionary->memoryAllocated|=SECOND_PART_OF_TABLES_IN_DICTIONARY_ALLOCATED;
}

void freeDictionary(GaborDictionary *gaborDictionary)
{
    unsigned int gabor;

    if(gaborDictionary->memoryAllocated & FIRST_PART_OF_TABLES_IN_DICTIONARY_ALLOCATED)
    {
	usiVectorFree(gaborDictionary->tableOfScalesInOptimalDictionary);
	uiVectorFree(gaborDictionary->tableOfPeriodsInOptimalDictionary);
	dVectorFree(gaborDictionary->tableOfFrequenciesInOptimalDictionary);
	dVectorFree(gaborDictionary->tableOfPositionsInOptimalDictionary);    
	uiVectorFree(gaborDictionary->numberOfStepsInFrequencyAtParticularScale);
	usiVectorFree(gaborDictionary->numberOfStepsInPositionAtParticularScale);
    }
    else if(gaborDictionary->memoryAllocated & SECOND_PART_OF_TABLES_IN_DICTIONARY_ALLOCATED)
    {
	for(gabor=0;gabor<gaborDictionary->sizeOfDictionary;gabor++)
	    freeGaborElements(gaborDictionary->gaborsTable + gabor);

	free(gaborDictionary->gaborsTable);
    }
}

void makeDictionary(const MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary)
{
    time_t       seedTime;
    long int     seed;
    Gabor        *gaborPointer = NULL;

    const unsigned short int dimOffset = mp5Parameters->dimOffset;
    unsigned short int leftSize, rightSize;

    const unsigned short int numberOfStepsInScale = gaborDictionary->numberOfStepsInScale;

    const double *tableOfFrequenciesInOptimalDictionary  = gaborDictionary->tableOfFrequenciesInOptimalDictionary;
    const double *tableOfPositionsInOptimalDictionary    = gaborDictionary->tableOfPositionsInOptimalDictionary;

    const unsigned       int *numberOfStepsInFrequencyAtParticularScale = gaborDictionary->numberOfStepsInFrequencyAtParticularScale;
    const unsigned short int *numberOfStepsInPositionAtParticularScale  = gaborDictionary->numberOfStepsInPositionAtParticularScale;

    const unsigned int periodDensity = gaborDictionary->periodDensity;

    unsigned       int numberOfStepsInFrequency;
    unsigned short int numberOfStepsInPosition;
    
    unsigned short int scaleIndex = 0;
    unsigned       int frequencyIndex;
    unsigned short int positionIndex;
    unsigned       int rifling    = 0;
    double             frequency  = 0.0;


    double DF              = 0.0;
    double DU              = 0.0;
    double stepInFrequency = 0.0; 
    
    seed = time(&seedTime);

    r250_init(seed);	

    printf("\n START DICTIONARY GENERAITING \n\n");

    gaborPointer = gaborDictionary->gaborsTable;

    unsigned counter = 0;

    /* set Dirac Delta */
    for(positionIndex=0;positionIndex<(double)(dimOffset);positionIndex++)
    {
	gaborPointer->scaleIndex  = 0; // although scaleIndex = 0 means scale = 1, we set this value, 
    				       // because this fild cand not be empty. However in case of Dirac Delta 
				       // this fild is not used, so it will not cause any disturbance 
	gaborPointer->rifling     = 0;
	gaborPointer->position    = (unsigned short int)positionIndex;
	gaborPointer->feature|=DIRACDELTA;

	leftSize  = gaborPointer->position;
	rightSize = (unsigned short int)(dimOffset - gaborPointer->position - 1);

	if(leftSize<=rightSize)
	    gaborPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;
		
	gaborPointer++;
	counter++;
    }

    /* set Gabor waves */
    for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
    {
	DF = *(tableOfFrequenciesInOptimalDictionary + scaleIndex);
        DU = *(tableOfPositionsInOptimalDictionary   + scaleIndex);  

	numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + scaleIndex);	
	numberOfStepsInPosition  = *(numberOfStepsInPositionAtParticularScale + scaleIndex);	

	stepInFrequency = DF*periodDensity;

	if(gaborDictionary->typeOfDictionary & OCTAVE_FIXED)
	{            
	    for(frequencyIndex=0,rifling=1;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++,rifling=rifling+periodDensity)
    	    {
		for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
		{
		    gaborPointer->scaleIndex = scaleIndex;
		    gaborPointer->rifling    = rifling;
		    gaborPointer->position   = (unsigned short int)(DU*(positionIndex + 1));	

		    gaborPointer->feature|=GABORWAVE;
		    leftSize  = gaborPointer->position;
		    rightSize = (unsigned short int)(dimOffset - gaborPointer->position - 1);

		    if(leftSize<=rightSize)
			gaborPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;

		    gaborPointer++;
		    counter++;
		}  	    
	    }
	}
	else if(gaborDictionary->typeOfDictionary & OCTAVE_STOCH)
	{
	    for(frequencyIndex=0;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++)
    	    {
		for(positionIndex=0;positionIndex<numberOfStepsInPosition;positionIndex++)
		{
		    gaborPointer->scaleIndex = scaleIndex;
		    frequency                = stepInFrequency*(frequencyIndex + 1) + DF*r250n(periodDensity);                     
		    gaborPointer->rifling    = (unsigned int)(frequency/DF + 0.5);

		    gaborPointer->position   = (unsigned short int)(DU*(positionIndex + 1) + DU*dr250());	

		    gaborPointer->feature|=GABORWAVE;
		    leftSize  = gaborPointer->position;
		    rightSize = (unsigned short int)(dimOffset - gaborPointer->position - 1);

		    if(leftSize<=rightSize)
			gaborPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;

		    gaborPointer++;
		    counter++;
		}
	    }
	}
    }

    DF = *(tableOfFrequenciesInOptimalDictionary + numberOfStepsInScale - 1);
    stepInFrequency = DF*periodDensity;
    numberOfStepsInFrequency = *(numberOfStepsInFrequencyAtParticularScale + numberOfStepsInScale - 1);	

    for(frequencyIndex=0,rifling=1;frequencyIndex<numberOfStepsInFrequency;frequencyIndex++,rifling=rifling+periodDensity)
    {
	gaborPointer->scaleIndex = (unsigned short int)(numberOfStepsInScale - 1);  // although scaleIndex = numberOfScale - 1 means scale = dimOffset, we set this value, 
							    			    // because this fild cand not be empty. However in case of FFT WAVE 
							    			    // this fild is not used, so it will not cause any disturbance 
	gaborPointer->rifling    = rifling;
	gaborPointer->position   = (unsigned short int)(dimOffset/2);
	gaborPointer->feature|=FFTWAVE;

	leftSize  = gaborPointer->position;
	rightSize = (unsigned short int)(dimOffset - gaborPointer->position - 1);

	if(leftSize<= rightSize)
	    gaborPointer->feature|=LEFT_SIDE_POSITION_IN_OFFSET;

	gaborPointer++;
	counter++;
    }
}

void printDictionaryToAsciFile(const DataParameters *dataParameters, const GaborDictionary *gaborDictionary)
{
	unsigned int gaborCounter;
	char typeOfGabor;
	char statusOfGabor;
	Gabor *gabor;
	double frequency;

	for(gaborCounter=0;gaborCounter<gaborDictionary->sizeOfDictionary;gaborCounter++)
	{
		gabor = gaborDictionary->gaborsTable + gaborCounter;

		if(gabor->feature & DIRACDELTA)
			typeOfGabor = 'D';
		else if(gabor->feature & GABORWAVE)
			typeOfGabor = 'G';
		else
			typeOfGabor = 'F';

		statusOfGabor = (char)(gabor->feature & INCORRECTGABOR ? 'I' : 'C');
		frequency = *(gaborDictionary->tableOfFrequenciesInOptimalDictionary + gabor->scaleIndex)*gabor->rifling;

		fprintf(dataParameters->dictionaryFile," %12u   %c   %5hu   %5hu   %7.6lf   %12u   %012.6lf   %5c \n",gaborCounter,typeOfGabor,gabor->position,*(gaborDictionary->tableOfScalesInOptimalDictionary + gabor->scaleIndex),frequency,gabor->rifling,frequency*dataParameters->samplingFrequency/M_2PI,statusOfGabor);
		fflush(stdout);
	}
}

void reinitDictionary(const MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary)
{
    unsigned int gaborsCounter;
    
    for(gaborsCounter=0;gaborsCounter<gaborDictionary->sizeOfDictionary;gaborsCounter++)
	(gaborDictionary->gaborsTable + gaborsCounter)->feature = 0x0;

    makeDictionary(mp5Parameters,gaborDictionary);

}

void testGaborFeature(GaborDictionary *gaborDictionary)
{
    Gabor     *gaborPointer = NULL;

    unsigned int correctGaborsCounter = 0;
    unsigned int gaborsCounter = 0;

    printf(" THE INCORRECT GABORS ARE BEING REMOVED FROM DICTIONARY \n");
    printf(" UNCORRECT GABORS ARE GABORS WHICH SCALE IS SMALLER THEN K*PERIOD WAVE IN GABOR\n");
    fflush(stdout);
    
    for(gaborsCounter=0;gaborsCounter<gaborDictionary->sizeOfDictionary;gaborsCounter++)
    {
	gaborPointer = gaborDictionary->gaborsTable + gaborsCounter;
	if(gaborPointer->feature & GABORWAVE)
	    testScaleToPeriodFactorGabor(gaborDictionary,gaborPointer,gaborDictionary->scaleToPeriodFactor);
    }

    correctGaborsCounter = 0;

    for(gaborsCounter=0;gaborsCounter<gaborDictionary->sizeOfDictionary;gaborsCounter++)
    {
	gaborPointer = gaborDictionary->gaborsTable + gaborsCounter;
	
	if(!(gaborPointer->feature & INCORRECTGABOR))
	    correctGaborsCounter++;
    }

    printf("\n AFTER PREPROCESSING, THE NUMBER OF CORRECT GABORS IS: %u \n",correctGaborsCounter);
    fflush(stdout);

}

void resetDictionary(GaborDictionary *gaborDictionary)
{
    unsigned int gaborsCounter;
    const unsigned int sizeOfDictionary = gaborDictionary->sizeOfDictionary;

    for(gaborsCounter=0;gaborsCounter<sizeOfDictionary;gaborsCounter++)
	(gaborDictionary->gaborsTable + gaborsCounter)->feature&= ~GABOR_WAS_HIT;

}
