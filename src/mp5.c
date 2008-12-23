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
#include<string.h>
#include<strings.h>
#include<unistd.h>
#include"atom.h"
#include"def.h"
#include"dic.h"
#include"fftw3.h"
#include"matrix.h"
#include"mp5.h"
#include"queue.h"
#include"tools.h"
#include"types.h"
#include"vector.h"

extern unsigned char applicationMode;

#ifdef __MINGW32__
    #define bzero(ptr,size) memset (ptr, 0, size);
	#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

#define findStartResample(sample,rifling,period)         ((sample*rifling)%period);
#define findNextResample(currentResample,rifling,period) ((currentResample+rifling)%period);

static unsigned int uiabs(unsigned int firstNumber, unsigned int secondNumber)
{
	return (firstNumber<=secondNumber) ? (secondNumber-firstNumber) : (firstNumber - secondNumber);
}

static void findFastlySinCos(const double *basicSin, const double *basicCos, double *newSin, double *newCos)
{
	const double oldCos = *newCos;
	const double oldSin = *newSin;
			
	*newCos = oldCos*(*basicCos) - oldSin*(*basicSin);
	*newSin = oldSin*(*basicCos) + oldCos*(*basicSin);	
}

static void findStartAndStopConditionsInFullRange(unsigned int position,
						  unsigned     int marginalDimension,
						  unsigned     int exponensTableDimension,
						  unsigned int *firstStart,
						  unsigned int *firstStop,
						  unsigned int *secondStart,
						  unsigned int *secondStop,
						  unsigned char feature)
{
/* In this implementation of MP algorithm calculation are performed in the range
 from 0 to N-1, where N is an offset's length.
 However, one can make MP estimation faster, if the fatures of atoms function will be
 used. 'Cosine-Atom' is a even function and 'Sine-Atom' is a odd function, therefore
 we don't have to evaluate dot products for the whole offset. Morover, Atom function deseaper
 very fast, so we can estimate dot products only there, where value of Atom Function is grater then
 some EPS. This function is looking for the beginning and the end of offset, where calculation
 should be permpormed. */

    const unsigned int startPosition = marginalDimension + position;
    const unsigned int endOfTable    = exponensTableDimension - position;
    
	if(feature & LEFT_SIDE_POSITION_IN_OFFSET)
    {
		*firstStart = 0;
		*firstStop  = startPosition;

		*secondStart = startPosition + 1;
		*secondStop  = endOfTable - 1;
    }
    else
    {
		*firstStart = 0;
		*firstStop  = endOfTable - 1;

		*secondStart = endOfTable;
		*secondStop  = startPosition;
    }
}

static void findStartAndStopConditionsInLimitedRange(const unsigned int atomsPosition,
													 unsigned int *startInterval,
													 unsigned int *stopInterval,
													 unsigned int *firstStart,
													 unsigned int *firstStop,
													 unsigned int *secondStart,
													 unsigned int *secondStop)
{
/* In this implementation of MP algorithm calculation are performed in the range
 from 0 to N-1, where N is an offset's length.
 However, one can make MP estimation faster, if the fatures of atoms function will be
 used. 'Cosine-Atom' is a even function and 'Sine-Atom' is a odd function, therefore
 we don't have to evaluate dot products for the whole offset. Morover, Atom function deseaper
 very fast, so we can estimate dot products only there, where value of Atom Function is grater then
 some EPS. This function is looking for the beginning and the end of offset, where calculation
 should be permpormed. */

	unsigned int newStartInterval, newStopInterval, tmpValue;

    const unsigned int tmpStartInterval = *startInterval;
    const unsigned int tmpStopInterval  = *stopInterval;

    newStartInterval = uiabs(*startInterval,atomsPosition);
    newStopInterval  = uiabs(*stopInterval,atomsPosition);

    if(newStartInterval>newStopInterval)
    {
		tmpValue         = newStartInterval;
		newStartInterval = newStopInterval;
		newStopInterval  = tmpValue;
    }

    if((*startInterval<atomsPosition && *stopInterval<atomsPosition) || (*startInterval>atomsPosition && *stopInterval>atomsPosition))
    {
		*firstStart = 1; /* these condition seem to be very strange, but thanks to them */
		*firstStop  = 0; /* some loops will not start, what is needed by us :-) */

		*secondStart = newStartInterval;
		*secondStop  = newStopInterval;
    }
    
	if(((*startInterval<=atomsPosition) && (*stopInterval>=atomsPosition)) || ((*startInterval>=atomsPosition) && (*stopInterval<=atomsPosition)))
    {
		*firstStart = 0;
		*firstStop  = newStartInterval;

		*secondStart = newStartInterval + 1; // when secondStart = secondStop condition seems to be strange, but thanks to iy
		*secondStop  = newStopInterval;      // some loops will not start
    }

    *startInterval = tmpStartInterval;
    *stopInterval  = tmpStopInterval;

}

static void findGaussNonGaussInterval(unsigned int gaborPosition, 
									  double gaborScale,
									  unsigned int marginalDimension, 
									  unsigned int *intervalCenter, 
									  unsigned int *intervalRange)
{
	*intervalCenter   = marginalDimension + gaborPosition;
	*intervalRange    = (unsigned int)(sqrt(-LOG_EPS_DOT_PRODUCT/M_PI)*gaborScale + 1.5);
	
	*intervalRange  = *intervalRange>marginalDimension ? marginalDimension : *intervalRange;
}

static double findGaussGaussInterval(const Dictionary *dictionary, 
									 const Atom *currentAtom, 
									 const Atom *previousAtom, 
									 unsigned int marginalDimension,
									 unsigned int *intervalCenter,
									 unsigned int *intervalRange)
{
	unsigned int currentAtomPosition = currentAtom->position;
	double currentAtomScale    = *(dictionary->tableOfScalesInOptimalDictionary + currentAtom->scaleIndex);

	unsigned int previousAtomPosition = previousAtom->position;
	double  previousAtomScale    = *(dictionary->tableOfScalesInOptimalDictionary + previousAtom->scaleIndex);

	const double sqrPreviousAtomScale = previousAtomScale*previousAtomScale;
	const double sqrCurrentAtomScale  = currentAtomScale*currentAtomScale;
	const double newScale = sqrPreviousAtomScale + sqrCurrentAtomScale;
	
	double K = -((M_PI*((previousAtomPosition - currentAtomPosition)*(previousAtomPosition - currentAtomPosition)))/newScale);
	
	if(K>=LOG_EPS_DOT_PRODUCT)
	{
		*intervalCenter = marginalDimension + (unsigned int)(0.5 + (currentAtomPosition*sqrPreviousAtomScale + previousAtomPosition*sqrCurrentAtomScale)/newScale);
		*intervalRange  = (unsigned int)(1.5 + sqrt((K - LOG_EPS_DOT_PRODUCT)/(M_PI*newScale/(sqrCurrentAtomScale*sqrPreviousAtomScale))));
		*intervalRange  = *intervalRange>marginalDimension ? marginalDimension : *intervalRange;
	}
	
	return K;
}

static void findRSRCVariables(const MP5Parameters *mp5Parameters,
			      const Dictionary *dictionary,
			      const Atom   *atom,
			      const double *dataTable,
			      unsigned int firstStart,
			      unsigned int firstStop,
			      unsigned int secondStart,
			      unsigned int secondStop,
			      unsigned char mode,
			      double *RS,
			      double *RC)
{

    unsigned             int sample;
    unsigned             int resample;
	
    const unsigned short int scaleIndex = atom->scaleIndex;
    const unsigned       int rifling    = atom->rifling;
    const unsigned       int position   = atom->position;
    const unsigned char      feature    = atom->feature; 
	      unsigned char      located;

    const unsigned       int marginalDimension = mp5Parameters->marginalDimension;
    const unsigned       int startPosition   = marginalDimension + position;

    const double *sinTable = *(mp5Parameters->sinTable + scaleIndex);
    const double *cosTable = *(mp5Parameters->cosTable + scaleIndex);
    const double *expTable = *(mp5Parameters->expTable + scaleIndex);

    const unsigned int commonPeriod = *(dictionary->tableOfPeriodsInOptimalDictionary + scaleIndex);

    double sigLeftVal;
    double sigRightVal;
    double sinVal, cosVal;
    const double *ptrDataTable = dataTable + startPosition;

    if(mode & FULL_RANGE)
		located = (unsigned char)((feature & LEFT_SIDE_POSITION_IN_OFFSET) ? LEFT_SIDE_POSITION : RIGHT_SIDE_POSITION);
    else 
		located = (unsigned char)((feature & LEFT_SIDE_POSITION_IN_RANGE) ? LEFT_SIDE_POSITION : RIGHT_SIDE_POSITION);

    (*RS) = (*RC) = 0.0;

    if(feature & DIRACDELTA)
    {
		*RS = 0.0;
		*RC = *ptrDataTable;           
    }
	else if(feature & GAUSSFUNCTION)
	{
		double expVal;

		for(sample=firstStart;sample<=firstStop;sample++)
		{
			expVal      = *(expTable + sample);
			sigLeftVal  = *(ptrDataTable - sample);
			sigRightVal = *(ptrDataTable + sample);

			(*RC)+=expVal*(sigRightVal + sigLeftVal);
		}

		if(located & LEFT_SIDE_POSITION)
		{
			for(sample=secondStart;sample<=secondStop;sample++)
			{
				expVal      = *(expTable + sample);
				sigRightVal = *(ptrDataTable + sample);

				(*RC)+=expVal*sigRightVal;
			}
		}    
		else if(located & RIGHT_SIDE_POSITION)
		{
			for(sample=secondStart;sample<=secondStop;sample++)
			{
				expVal     = *(expTable + sample);
				sigLeftVal = *(ptrDataTable - sample);

				(*RC)+=expVal*sigLeftVal;
			}
		}
		(*RC)-= (*ptrDataTable);
	}
    else if(feature & SINCOSWAVE)
    {  
		resample = findStartResample(firstStart,rifling,commonPeriod);
		
        for(sample=firstStart;sample<=firstStop;sample++)
		{
			sinVal      = *(sinTable + resample);
            cosVal      = *(cosTable + resample);
			sigLeftVal  = *(ptrDataTable - sample);
    	    sigRightVal = *(ptrDataTable + sample);

            (*RS)+=sinVal*(sigRightVal - sigLeftVal);
			(*RC)+=cosVal*(sigRightVal + sigLeftVal);
   
			resample = findNextResample(resample,rifling,commonPeriod);
		}

		if(located & LEFT_SIDE_POSITION)
        {
			resample = findStartResample(secondStart,rifling,commonPeriod);

			for(sample=secondStart;sample<=secondStop;sample++)
			{
				sinVal      = *(sinTable + resample);
				cosVal      = *(cosTable + resample);
				sigRightVal = *(ptrDataTable + sample);

				(*RS)+=sinVal*sigRightVal;
				(*RC)+=cosVal*sigRightVal;

				resample = findNextResample(resample,rifling,commonPeriod);
			}
		}    
		else if(located & RIGHT_SIDE_POSITION)
		{
			resample = findStartResample(secondStart,rifling,commonPeriod);

			for(sample=secondStart;sample<=secondStop;sample++)
			{

				sinVal     = *(sinTable + resample);
				cosVal     = *(cosTable + resample);

				sigLeftVal = *(ptrDataTable - sample);

				(*RS)+= -sinVal*sigLeftVal;
				(*RC)+=  cosVal*sigLeftVal;
				resample = findNextResample(resample,rifling,commonPeriod);
			}
		}

		(*RC)-= (*ptrDataTable);
    }
    else
    {
		double expVal;

		resample = findStartResample(firstStart,rifling,commonPeriod);

		for(sample=firstStart;sample<=firstStop;sample++)
		{
			expVal      = *(expTable + sample);
			sinVal      = *(sinTable + resample);
			cosVal      = *(cosTable + resample);
			sigLeftVal  = *(ptrDataTable - sample);
			sigRightVal = *(ptrDataTable + sample);

			(*RS)+=(sinVal*expVal)*(sigRightVal - sigLeftVal);
			(*RC)+=(cosVal*expVal)*(sigRightVal + sigLeftVal);

			resample = findNextResample(resample,rifling,commonPeriod);
		}

		if(located & LEFT_SIDE_POSITION)
		{
			resample = findStartResample(secondStart,rifling,commonPeriod);

			for(sample=secondStart;sample<=secondStop;sample++)
			{
				sinVal      = *(sinTable + resample);
				cosVal      = *(cosTable + resample);
				expVal      = *(expTable + sample);
				sigRightVal = *(ptrDataTable + sample);

				(*RS)+=(sinVal*expVal)*sigRightVal;
				(*RC)+=(cosVal*expVal)*sigRightVal;
	
				resample = findNextResample(resample,rifling,commonPeriod);
			}
		}    
		else if(located & RIGHT_SIDE_POSITION)
		{

			resample = findStartResample(secondStart,rifling,commonPeriod);

			for(sample=secondStart;sample<=secondStop;sample++)
			{

				sinVal     = *(sinTable + resample);
				cosVal     = *(cosTable + resample);
				expVal     = *(expTable + sample);
				sigLeftVal = *(ptrDataTable - sample);

				(*RS)+= -(sinVal*expVal)*sigLeftVal;
				(*RC)+=  (cosVal*expVal)*sigLeftVal;

				resample = findNextResample(resample,rifling,commonPeriod);
			}
		}
		(*RC)-= (*ptrDataTable);
    }
}

static void makeOneSinCosTable(double omega, double *sinTable, double *cosTable, unsigned int period)
{

    unsigned int sample;
    double sinOmega, cosOmega;

    for(sample = 0;sample<period;sample++)
    {
	    sincos(omega*sample,&sinOmega,&cosOmega);
		*(sinTable + sample) = sinOmega;
		*(cosTable + sample) = cosOmega;			
    }
}

static void makeOneExpTable(double alpha, double *expTable, unsigned int exponensTableDimension)
{
    unsigned int sample;

    for(sample = 0;sample<exponensTableDimension;sample++)
    {
		*(expTable + sample) = exp(-(alpha*sample)*sample);
    }
}

static void findKSKCKMVariables(const MP5Parameters *mp5Parameters,
								const Dictionary *dictionary,   
								Atom *atom,
								unsigned int firstStart,
								unsigned int firstStop,
								unsigned int secondStart,
								unsigned int secondStop)
{

    register unsigned int sample;
    register unsigned int resample;
	register double sqrSin, sqrExp;
    register double sinVal, cosVal;

	double basicSin, basicCos;
	double newSin, newCos;
    
    double KS = 0.0;
    double KC = 0.0;
    double KM = 0.0;

    if(atom->feature & DIRACDELTA)
    {
		KS = 0.0;
		KC = 1.0;
        KM = 0.0;
    }
	else if(atom->feature & GAUSSFUNCTION)
	{
		const unsigned short int scaleIndex = atom->scaleIndex;
		const double *expTable = *(mp5Parameters->expTable + scaleIndex);
		double expVal;

		for(sample=firstStart;sample<=firstStop;sample++)
		{
			expVal = *(expTable + sample);
			KC+= expVal*expVal; 
		}

		KC = 2.0*KC - 1.0;

		for(sample=secondStart;sample<=secondStop;sample++)
		{
			expVal   = *(expTable + sample);
			KC+= expVal*expVal;
		}

	}
	else if(atom->feature & SINCOSWAVE)
	{
	
		unsigned             int rifling    = atom->rifling;
		const unsigned short int scaleIndex = atom->scaleIndex;
		const unsigned int commonPeriod     = *(dictionary->tableOfPeriodsInOptimalDictionary + scaleIndex);
		const double randomShiftInFrequency = atom->randomShiftInFrequency;
		const double *sinTable = *(mp5Parameters->sinTable + scaleIndex);
		const double *cosTable = *(mp5Parameters->cosTable + scaleIndex);

		resample = findStartResample(firstStart,rifling,commonPeriod);

		sincos(randomShiftInFrequency,&basicSin,&basicCos);	

		newCos = 1.0;
		newSin = 0.0;

		for(sample=firstStart;sample<=firstStop;sample++)
		{
			sinVal = (*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin;
//			cosVal = (*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin;
			sqrSin = sinVal*sinVal;
			KS+= sqrSin;
			KC+= (1 - sqrSin);
//			KC+= cosVal*cosVal;
			/* *KM = 0, because of symtrical sumation of elements of atoms */
			resample = findNextResample(resample,rifling,commonPeriod);
			findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
		}

		KS*= 2.0;
		KC = 2.0*KC - 1.0;

		resample = findStartResample(secondStart,rifling,commonPeriod);
		sincos(secondStart*randomShiftInFrequency,&basicSin,&basicCos);	

		newSin = basicSin;
		newCos = basicCos;

		for(sample=secondStart;sample<=secondStop;sample++)
		{
			sinVal   = (*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin;
			cosVal   = (*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin;

			sqrSin = sinVal*sinVal; 
		
			KS+= sqrSin;
			KC+= (1 - sqrSin);
			KM+= cosVal*sinVal;

			resample = findNextResample(resample,rifling,commonPeriod);
			findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
		}
    }
    else
    {	
		unsigned             int rifling    = atom->rifling;
		const unsigned short int scaleIndex = atom->scaleIndex;
		const unsigned int commonPeriod    = *(dictionary->tableOfPeriodsInOptimalDictionary + scaleIndex);
		const double randomShiftInFrequency = atom->randomShiftInFrequency;
		const double *expTable = *(mp5Parameters->expTable + scaleIndex);
	    const double *sinTable = *(mp5Parameters->sinTable + scaleIndex);
		const double *cosTable = *(mp5Parameters->cosTable + scaleIndex);

		double expVal;

		resample = findStartResample(firstStart,rifling,commonPeriod);

		sincos(randomShiftInFrequency,&basicSin,&basicCos);	

		newCos = 1.0;
		newSin = 0.0;
		
		for(sample=firstStart;sample<=firstStop;sample++)
		{
			sinVal = (*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin;
			cosVal = (*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin;
			expVal = *(expTable + sample);
			sqrSin = sinVal*sinVal;
			sqrExp = expVal*expVal;
			
			KS+= sqrSin*sqrExp;
			KC+= (1-sqrSin)*sqrExp;

			/* *KM = 0, because of symtrical sumation of elements of atoms */
    
			resample = findNextResample(resample,rifling,commonPeriod);
			findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
		}

		KS*= 2.0;
		KC = 2.0*KC - 1.0;

		resample = findStartResample(secondStart,rifling,commonPeriod);
		sincos(secondStart*randomShiftInFrequency,&newSin,&newCos);	

		sincos(secondStart*randomShiftInFrequency,&basicSin,&basicCos);	

		newSin = basicSin;
		newCos = basicCos;

		for(sample=secondStart;sample<=secondStop;sample++)
		{

			sinVal   = (*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin;
			cosVal   = (*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin;
			expVal   = *(expTable + sample);
			
			sqrSin = sinVal*sinVal;
			sqrExp = expVal*expVal;
			
			KS+= sqrSin*sqrExp;
			KC+= (1.0 - sqrSin)*sqrExp;
			KM+= cosVal*sinVal*sqrExp;

			resample = findNextResample(resample,rifling,commonPeriod);
			findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);			
		}
	}

    atom->KS = (float)KS;
    atom->KC = (float)KC;
    atom->KM = (float)KM;
	
}

void setNumberOfAnalysedChannelsAndNumberOfResultsFiles(MP5Parameters *mp5Parameters)
{
    if(mp5Parameters->MPType & SMP)
		mp5Parameters->numberOfAllocatedChannels = mp5Parameters->numberOfAnalysedChannels = 1;
    else if ((mp5Parameters->MPType & MMP1) ||  (mp5Parameters->MPType & MMP2) ||  (mp5Parameters->MPType & MMP3) ||  (mp5Parameters->MPType & MMP4))
    {
		if(mp5Parameters->MPType & MMP2)
			mp5Parameters->numberOfAllocatedChannels = 1;
		else
			mp5Parameters->numberOfAllocatedChannels = mp5Parameters->numberOfChosenChannels;

		mp5Parameters->numberOfAnalysedChannels = mp5Parameters->numberOfChosenChannels;
    }
}

void setMP5Parameters(const Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned short int numberOfChannelsInDataFile = mp5Parameters->numberOfChannelsInDataFile;
    unsigned short int numberOfAnalysedChannels   = mp5Parameters->numberOfAnalysedChannels;
    unsigned short int numberOfStepsInScale       = dictionary->numberOfStepsInScale;
    unsigned       int offsetDimension            = mp5Parameters->offsetDimension;
    unsigned       int offsetExpandedDimension    = mp5Parameters->offsetExpandedDimension;
    unsigned       int exponensTableDimension     = mp5Parameters->exponensTableDimension;
	unsigned       int fftTableDimension          = mp5Parameters->fftTableDimension;

    unsigned int *tableOfPeriodsInOptimalDictionary = dictionary->tableOfPeriodsInOptimalDictionary;
  
    mp5Parameters->rawDataMatrix       = dMatrixAllocate(numberOfChannelsInDataFile,offsetDimension);
    mp5Parameters->processedDataMatrix = dMatrixAllocate(numberOfChannelsInDataFile,offsetExpandedDimension);

    mp5Parameters->zeroSignalTable     = dVectorAllocate(offsetExpandedDimension);
	dSetVectorZero(mp5Parameters->zeroSignalTable,offsetExpandedDimension);
    mp5Parameters->gaussSignalTable    = dVectorAllocate(offsetExpandedDimension);
	
    mp5Parameters->fitted             = createQueue();
    mp5Parameters->bestModulusesTable = dVectorAllocate(numberOfAnalysedChannels);
    mp5Parameters->bestPhasesTable    = fVectorAllocate(numberOfAnalysedChannels);

    mp5Parameters->signalEnergyInEachChannel  = dVectorAllocate(numberOfAnalysedChannels);
    mp5Parameters->residueEnergyInEachChannel = dVectorAllocate(numberOfAnalysedChannels);

    mp5Parameters->meanSignalTable  = dVectorAllocate(offsetExpandedDimension);
    mp5Parameters->meanResidueTable = dVectorAllocate(offsetExpandedDimension);

    mp5Parameters->prevAtomTable    = dMatrixAllocate(numberOfAnalysedChannels,offsetExpandedDimension);

    /* sin/cos function are build for each scale */

    mp5Parameters->sinTable = dVariableMatrixAllocate(numberOfStepsInScale,tableOfPeriodsInOptimalDictionary);
    mp5Parameters->cosTable = dVariableMatrixAllocate(numberOfStepsInScale,tableOfPeriodsInOptimalDictionary);
    mp5Parameters->expTable = dMatrixAllocate(numberOfStepsInScale,exponensTableDimension);

	if(mp5Parameters->FFT)
	{
		mp5Parameters->fftTableInA = dVectorAllocate(fftTableDimension); 
		if(dictionary->typeOfDictionary & OCTAVE_STOCH)
			mp5Parameters->fftTableInB = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*fftTableDimension); // see the fftw documentation		

		if(dictionary->typeOfDictionary & OCTAVE_STOCH)
			mp5Parameters->fftTableOut = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*fftTableDimension);       // see the fftw documentation		
		else
			mp5Parameters->fftTableOut = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(fftTableDimension/2+1)); // see the fftw documentation		

	}

	mp5Parameters->sinAtomTable = dVectorAllocate(offsetExpandedDimension);
	mp5Parameters->cosAtomTable = dVectorAllocate(offsetExpandedDimension);

}

void freeMP5Parameters(const Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
	if(mp5Parameters->prevAtomTable!=NULL)
		dMatrixFree(mp5Parameters->prevAtomTable,mp5Parameters->numberOfAnalysedChannels);

	if(mp5Parameters->sinAtomTable!=NULL)
		dVectorFree(mp5Parameters->sinAtomTable);

	if(mp5Parameters->cosAtomTable!=NULL)
		dVectorFree(mp5Parameters->cosAtomTable);

	if(mp5Parameters->fftwPlan!=NULL)
		fftw_destroy_plan(mp5Parameters->fftwPlan);

	if(mp5Parameters->fftTableInA!=NULL)
		dVectorFree(mp5Parameters->fftTableInA);
		
	if(mp5Parameters->fftTableInB!=NULL)
		fftw_free(mp5Parameters->fftTableInB);
		
	if(mp5Parameters->fftTableOut!=NULL)
		fftw_free(mp5Parameters->fftTableOut);
				
	if(mp5Parameters->signalEnergyInEachChannel!=NULL)
		dVectorFree(mp5Parameters->signalEnergyInEachChannel);

	if(mp5Parameters->residueEnergyInEachChannel!=NULL)
		dVectorFree(mp5Parameters->residueEnergyInEachChannel);

	if(mp5Parameters->meanSignalTable!=NULL)
		dVectorFree(mp5Parameters->meanSignalTable);

	if(mp5Parameters->meanResidueTable!=NULL)
		dVectorFree(mp5Parameters->meanResidueTable);

	if(mp5Parameters->sinTable!=NULL)
		dVariableMatrixFree(mp5Parameters->sinTable,dictionary->numberOfStepsInScale);

	if(mp5Parameters->cosTable!=NULL)
		dVariableMatrixFree(mp5Parameters->cosTable,dictionary->numberOfStepsInScale);

	if(mp5Parameters->expTable!=NULL)
		dMatrixFree(mp5Parameters->expTable,dictionary->numberOfStepsInScale);

	if(mp5Parameters->bestModulusesTable!=NULL)
		dVectorFree(mp5Parameters->bestModulusesTable);
	
	if(mp5Parameters->bestPhasesTable!=NULL)
		fVectorFree(mp5Parameters->bestPhasesTable);
		
    if(mp5Parameters->chosenChannels!=NULL)
		usiVectorFree(mp5Parameters->chosenChannels);

    if(mp5Parameters->chosenOffsets!=NULL)
		usiVectorFree(mp5Parameters->chosenOffsets);

    if(mp5Parameters->rawDataMatrix!=NULL)
		dMatrixFree(mp5Parameters->rawDataMatrix,mp5Parameters->numberOfChannelsInDataFile);

    if(mp5Parameters->processedDataMatrix!=NULL)
		dMatrixFree(mp5Parameters->processedDataMatrix,mp5Parameters->numberOfChannelsInDataFile);		

   if(mp5Parameters->zeroSignalTable!=NULL)
		dVectorFree(mp5Parameters->zeroSignalTable);		

	if(mp5Parameters->gaussSignalTable!=NULL)
		dVectorFree(mp5Parameters->gaussSignalTable);		

	freeQueue(mp5Parameters->fitted,(void (*)(void *))freeAtom);
}

/* prepare Sin, Cos and Exp Tables */

void makeSinCosExpTable(const Dictionary *dictionary, MP5Parameters *mp5Parameters)
{

    unsigned int scaleIndex;
    double  *sinTable, *cosTable, *expTable;
    const unsigned short int numberOfStepsInScale = dictionary->numberOfStepsInScale;
    const unsigned       int exponensTableDimension = mp5Parameters->exponensTableDimension;
    double       omega, alpha;
    double       scale;
    unsigned int period;

    const double  *tableOfScalesInOptimalDictionary       = dictionary->tableOfScalesInOptimalDictionary;
    const double  *tableOfFrequenciesInOptimalDictionary  = dictionary->tableOfFrequenciesInOptimalDictionary;
    const unsigned int *tableOfPeriodsInOptimalDictionary = dictionary->tableOfPeriodsInOptimalDictionary;

	if(applicationMode & PROCESS_USER_MODE)
		printf(" START SIN, COS, EXP, DIRAC TABLES GENERAITING \n\n");

    for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
    {
		omega  = *(tableOfFrequenciesInOptimalDictionary + scaleIndex);
		period = *(tableOfPeriodsInOptimalDictionary     + scaleIndex);
		sinTable = *(mp5Parameters->sinTable + scaleIndex);
		cosTable = *(mp5Parameters->cosTable + scaleIndex);
		makeOneSinCosTable(omega,sinTable,cosTable,period);
    }

    for(scaleIndex=0;scaleIndex<numberOfStepsInScale;scaleIndex++)
    {
		scale = *(tableOfScalesInOptimalDictionary + scaleIndex);

		expTable = *(mp5Parameters->expTable + scaleIndex);

		alpha = M_PI/(scale*scale);

		makeOneExpTable(alpha,expTable,exponensTableDimension);
    }

	if(applicationMode & PROCESS_USER_MODE)
		printf(" END SIN, COS, EXP, DIRAC TABLES GENERAITING \n\n");
}


void makeSinCosExpAtomTable(const Dictionary *dictionary,
						    MP5Parameters *mp5Parameters, 
                            const Atom  *atom)
{

    unsigned       int sample;
    unsigned       int resample;
	double basicSin, basicCos;
	double newSin, newCos;  

    const unsigned short int scaleIndex = atom->scaleIndex;
    const unsigned       int rifling    = atom->rifling;
    const unsigned       int position   = atom->position;
    const unsigned char      feature    = atom->feature;
    unsigned char            located;

    const unsigned       int marginalDimension       = mp5Parameters->marginalDimension;
    const unsigned       int exponensTableDimension  = mp5Parameters->exponensTableDimension;
    const unsigned       int offsetExpandedDimension = mp5Parameters->offsetExpandedDimension;

    const unsigned int startPosition = marginalDimension + position;
    unsigned int firstStart  = 0;
    unsigned int firstStop   = 0;
    unsigned int secondStart = 0;
    unsigned int secondStop  = 0;

    const unsigned int commonPeriod    = *(dictionary->tableOfPeriodsInOptimalDictionary + scaleIndex);
	const double       randomShiftInFrequency = atom->randomShiftInFrequency;

    const double *sinTable = *(mp5Parameters->sinTable + scaleIndex);
    const double *cosTable = *(mp5Parameters->cosTable + scaleIndex);
    const double *expTable = *(mp5Parameters->expTable + scaleIndex);

    double *sinAtomTable = mp5Parameters->sinAtomTable;
    double *cosAtomTable = mp5Parameters->cosAtomTable;

    double *ptrSinAtomTable, *ptrCosAtomTable;

    findStartAndStopConditionsInFullRange(position,marginalDimension,exponensTableDimension,&firstStart,&firstStop,&secondStart,&secondStop,atom->feature);

    ptrSinAtomTable = sinAtomTable + startPosition;
    ptrCosAtomTable = cosAtomTable + startPosition;

    located = (unsigned char)((atom->feature & LEFT_SIDE_POSITION_IN_OFFSET) ? LEFT_SIDE_POSITION : RIGHT_SIDE_POSITION);

    if(feature & DIRACDELTA)
    {	
		dSetVectorZero(sinAtomTable,offsetExpandedDimension);
        dSetVectorZero(cosAtomTable,offsetExpandedDimension);
        *ptrCosAtomTable = 1.0;
    }
	else if(feature & GAUSSFUNCTION)
	{
		double expVal;

		for(sample=firstStart;sample<=firstStop;sample++)
		{
			expVal = *(expTable + sample);
			*(ptrCosAtomTable  - sample) = *(ptrCosAtomTable + sample) = expVal;
		}

		if(located & LEFT_SIDE_POSITION)
		{	          
			for(sample=secondStart;sample<=secondStop;sample++)
			{
				expVal = *(expTable + sample);
				*(ptrCosAtomTable  + sample) = expVal;
			}
		}
		else if(located & RIGHT_SIDE_POSITION)
		{
			for(sample=secondStart;sample<=secondStop;sample++)
			{
				expVal = *(expTable + sample);
				*(ptrCosAtomTable - sample) = expVal;
			}
		}
	}
    else if(feature & SINCOSWAVE)
    {
		resample = findStartResample(firstStart,rifling,commonPeriod);

		sincos(randomShiftInFrequency,&basicSin,&basicCos);	
		
		newCos = 1.0;
		newSin = 0.0;
		
        for(sample=firstStart;sample<=firstStop;sample++)
        {
			*(ptrSinAtomTable  + sample) = (*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin;
			*(ptrSinAtomTable  - sample) = -(*(ptrSinAtomTable + sample));
			*(ptrCosAtomTable  - sample) = *(ptrCosAtomTable + sample) = (*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin;

			resample = findNextResample(resample,rifling,commonPeriod);
			
			findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
		}
        if(located & LEFT_SIDE_POSITION)
		{
			resample = findStartResample(secondStart,rifling,commonPeriod);
              
			sincos(secondStart*randomShiftInFrequency,&basicSin,&basicCos);	

			newSin = basicSin;
			newCos = basicCos;
			  
			for(sample=secondStart;sample<=secondStop;sample++)
			{
				*(ptrSinAtomTable  + sample) = (*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin;
				*(ptrCosAtomTable  + sample) = (*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin;;

				resample = findNextResample(resample,rifling,commonPeriod);
				
				findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
			}
		}
		else if(located & RIGHT_SIDE_POSITION)
		{
			resample = findStartResample(secondStart,rifling,commonPeriod);

			sincos(secondStart*randomShiftInFrequency,&basicSin,&basicCos);	

			newSin = basicSin;
			newCos = basicCos;
		
			for(sample=secondStart;sample<=secondStop;sample++)
			{	
				*(ptrSinAtomTable - sample) = -((*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin);
				*(ptrCosAtomTable - sample) =   (*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin;;

				resample = findNextResample(resample,rifling,commonPeriod);
				
				findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
			}
		}
    }
    else
    {
		double expVal;

		resample = findStartResample(firstStart,rifling,commonPeriod);

		sincos(randomShiftInFrequency,&basicSin,&basicCos);	
		
		newCos = 1.0;
		newSin = 0.0;

		for(sample=firstStart;sample<=firstStop;sample++)
		{
			expVal = *(expTable + sample);

			*(ptrSinAtomTable  + sample) = ((*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin)*expVal;
			*(ptrSinAtomTable  - sample) = -(*(ptrSinAtomTable + sample));
			*(ptrCosAtomTable  - sample) = *(ptrCosAtomTable + sample) = ((*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin)*expVal;
			
			resample = findNextResample(resample,rifling,commonPeriod);
			
			findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
		}
		
		if(located & LEFT_SIDE_POSITION)
		{	          
			resample = findStartResample(secondStart,rifling,commonPeriod);

			sincos(secondStart*randomShiftInFrequency,&basicSin,&basicCos);	

			newSin = basicSin;
			newCos = basicCos;
			
			for(sample=secondStart;sample<=secondStop;sample++)
			{
				expVal = *(expTable + sample);
				*(ptrSinAtomTable  + sample) = ((*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin)*expVal;
				*(ptrCosAtomTable  + sample) = ((*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin)*expVal;

				resample = findNextResample(resample,rifling,commonPeriod);
				
				findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
			}
		}
		else if(located & RIGHT_SIDE_POSITION)
		{
			resample = findStartResample(secondStart,rifling,commonPeriod);

			sincos(secondStart*randomShiftInFrequency,&basicSin,&basicCos);	

			newSin = basicSin;
			newCos = basicCos;
			
			for(sample=secondStart;sample<=secondStop;sample++)
			{
				expVal = *(expTable + sample);
				*(ptrSinAtomTable - sample) = -((*(sinTable + resample))*newCos + (*(cosTable + resample))*newSin)*expVal;
				*(ptrCosAtomTable - sample) = ((*(cosTable + resample))*newCos - (*(sinTable + resample))*newSin)*expVal;

				resample = findNextResample(resample,rifling,commonPeriod);
				
				findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
			}
		}
    }
}

void normAtomTable(const Dictionary *dictionary, const MP5Parameters *mp5Parameters, Atom *atom)
{

    const unsigned int atomPosition      = atom->position;
    const unsigned int marginalDimension = mp5Parameters->marginalDimension;
    const unsigned int exponensTableDimension = mp5Parameters->exponensTableDimension;

    unsigned       int firstStart  = 0;
    unsigned       int secondStart = 0;
    unsigned       int firstStop   = 0;
    unsigned       int secondStop  = 0;

	if(!(atom->feature & DIRACDELTA))
	{
		if(atom->feature & SINCOSWAVE)
		{
			findStartAndStopConditionsInFullRange(atomPosition,marginalDimension,exponensTableDimension,&firstStart,&firstStop,&secondStart,&secondStop,atom->feature);
		}	
		else
		{
			double atomScale  = *(dictionary->tableOfScalesInOptimalDictionary + atom->scaleIndex);

			unsigned int intervalCenter = 0;
			unsigned int intervalRange  = 0;
			
			findGaussNonGaussInterval(atomPosition, 
									  atomScale, 
								      marginalDimension,
									  &intervalCenter, 
							          &intervalRange);
									 
			unsigned int startInterval = intervalCenter - intervalRange;
			unsigned int stopInterval  = intervalCenter + intervalRange;
			
			findStartAndStopConditionsInLimitedRange(atomPosition + marginalDimension,&startInterval,&stopInterval,&firstStart,&firstStop,&secondStart,&secondStop);
		}
	}

    findKSKCKMVariables(mp5Parameters,dictionary,atom,firstStart,firstStop,secondStart,secondStop);                       
} 


void makeAtomTable(MP5Parameters *mp5Parameters, const Atom *atom, unsigned short int channelNumber)
{
    unsigned int sample;
    const unsigned int offsetExpandedDimension = mp5Parameters->offsetExpandedDimension;

    const double *cosAtomTable = mp5Parameters->cosAtomTable;
    const double *sinAtomTable = mp5Parameters->sinAtomTable;
    
    double *prevAtomTable = *(mp5Parameters->prevAtomTable + channelNumber);

    double KS = atom->KS;
    double KC = atom->KC;
    double KM = atom->KM;

    const double phase = *(atom->phase + channelNumber);

    double amplitude;
    double sinPhase, cosPhase;

    sincos(phase,&sinPhase,&cosPhase);

    const double sinPart    = KS*(sinPhase*sinPhase);
    const double cosPart    = KC*(cosPhase*cosPhase);
    const double sinCosPart = 2.0*KM*sinPhase*cosPhase;

    amplitude = sqrt((sinPart + cosPart) - sinCosPart);

    for(sample=0;sample<offsetExpandedDimension;sample++)
		*(prevAtomTable + sample) = ((*(cosAtomTable + sample))*cosPhase - (*(sinAtomTable + sample))*sinPhase)/amplitude;
	
}

double findSignalEnergy(const double *signalTable, unsigned int offsetExpandedDimension)
{
    unsigned int sample;
    double energy = 0.0;

    for(sample=0;sample<offsetExpandedDimension;sample++)
		energy+= (*(signalTable + sample))*(*(signalTable + sample));

    return energy;
}

void findResidue(double *residueTable, const double *atomTable, double modulus, unsigned int offsetExpandedDimension)
{
    unsigned int sample;
    
    for(sample=0;sample<offsetExpandedDimension;sample++)
		*(residueTable + sample) = *(residueTable + sample) - modulus*(*(atomTable + sample));
}

void findAtomDataDotProduct(const Dictionary *dictionary, 
							const MP5Parameters *mp5Parameters,
							Atom *currentAtom,
							const double *dataTable,
							unsigned short int channelNumber,
							unsigned char mode)
{

    const unsigned       int marginalDimension       = mp5Parameters->marginalDimension;
    const unsigned       int exponensTableDimension  = mp5Parameters->exponensTableDimension;

    unsigned       int firstStart  = 0;
    unsigned       int firstStop   = 0;
    unsigned       int secondStart = 0;
    unsigned       int secondStop  = 0;

    double RS, RC;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;

    Atom *previousAtom = mp5Parameters->previousAtom;

    RS = RC = 0.0;
	
	if((previousAtom!=NULL) && (previousAtom->feature & SINCOSWAVE))
	{
		findStartAndStopConditionsInFullRange(currentAtom->position,marginalDimension,exponensTableDimension,&firstStart,&firstStop,&secondStart,&secondStop,currentAtom->feature);
		findRSRCVariables(mp5Parameters,dictionary,currentAtom,dataTable,
						  firstStart,firstStop,secondStart,secondStop,
						  FULL_RANGE,&RS,&RC);	
	}
	else if(currentAtom->feature & DIRACDELTA)
	{
		findRSRCVariables(mp5Parameters,dictionary,currentAtom,dataTable,
						  firstStart,firstStop,secondStart,secondStop,
						  FULL_RANGE,&RS,&RC);
	}
	else if((currentAtom->feature & SINCOSWAVE))
	{
			findStartAndStopConditionsInFullRange(currentAtom->position,marginalDimension,exponensTableDimension,&firstStart,&firstStop,&secondStart,&secondStop,currentAtom->feature);

			findRSRCVariables(mp5Parameters,dictionary,currentAtom,dataTable,
							  firstStart,firstStop,secondStart,secondStop,
							  FULL_RANGE,&RS,&RC);	
	}
	else // GAUSS FUNCTION OR GABOR FUNCTION
	{
		unsigned int currentAtomPosition = currentAtom->position;
		unsigned int  intervalCenter = 0;
		unsigned int  intervalRange  = 0;

		unsigned int startInterval;
		unsigned int stopInterval;
		
		double K = 0.0;

		if(mode & FIRST_ITERATION)
		{
			K = -LOG_EPS_DOT_PRODUCT;

			findGaussNonGaussInterval(currentAtomPosition,
									  *(dictionary->tableOfScalesInOptimalDictionary + currentAtom->scaleIndex), 
									  marginalDimension, 
									  &intervalCenter, 
									  &intervalRange);
									  						  
		}
		else
		{
			K = findGaussGaussInterval(dictionary,
									   currentAtom,
									   previousAtom, 
									   marginalDimension,
									   &intervalCenter,
									   &intervalRange);									   
  		}
		
		if(K>=LOG_EPS_DOT_PRODUCT)
		{
			startInterval = intervalCenter - intervalRange;
			stopInterval  = intervalCenter + intervalRange;
		
			findStartAndStopConditionsInLimitedRange(currentAtomPosition + marginalDimension,&startInterval,&stopInterval,&firstStart,&firstStop,&secondStart,&secondStop);

			if((currentAtomPosition + marginalDimension)<startInterval)
				currentAtom->feature|=LEFT_SIDE_POSITION_IN_RANGE;
			else if((currentAtomPosition + marginalDimension)>stopInterval)
				currentAtom->feature&= (~LEFT_SIDE_POSITION_IN_RANGE);
			else
			{
				const unsigned int sum = stopInterval + startInterval;

				if(2*(currentAtomPosition + marginalDimension)<=sum)
					currentAtom->feature|= LEFT_SIDE_POSITION_IN_RANGE;
				else 
					currentAtom->feature&= ~LEFT_SIDE_POSITION_IN_RANGE;
			}
			
			findRSRCVariables(mp5Parameters,dictionary,currentAtom,dataTable,firstStart,firstStop,secondStart,secondStop,LIMIT_RANGE,&RS,&RC);			

		}
		else
		{
			RS = 0.0;
			RC = 0.0;
		}
	}

	
	if((mode & FIRST_ITERATION) || (mode & MMP1_NEXT_ITERATION))
	{
		*(currentAtom->RS + channelNumber) = (float)RS;
		*(currentAtom->RC + channelNumber) = (float)RC;
	}
	else
	{
		RS = (*(currentAtom->RS + channelNumber)) - (*(bestModulusesTable + channelNumber))*RS;
		RC = (*(currentAtom->RC + channelNumber)) - (*(bestModulusesTable + channelNumber))*RC;
		
		*(currentAtom->RS + channelNumber) = (float)RS;
		*(currentAtom->RC + channelNumber) = (float)RC;
	}
}			

static void makeFFTTable(const double *signalTable,
						 const double *expTable,
						 double *gaussTable, 
						 double *fftTableIn, 
						 const unsigned int intervalCenter,
						 const unsigned int atomsPositionInExpandedSignal,
						 const unsigned int fftDimension,
						 unsigned int *halfOfFFTDimension,
						 const unsigned char feature,
						 const unsigned char mode)
{

	int incx;
	int incy;
	int incz;

	if((fftDimension%2)==0)
		*halfOfFFTDimension = (fftDimension - 1)/2;
	else
		*halfOfFFTDimension = fftDimension/2;
		
	const double *shiftedSignalTable = signalTable + intervalCenter - *halfOfFFTDimension;
			
	if(mode & FIRST_ITERATION)
	{
		double *shiftedGaussTable  = gaussTable + *halfOfFFTDimension;

		if(feature & SINCOSWAVE)
		{
			/* from index table = [0] to [halfOfFFTDimension - 1] */
			incx =  1;
			incy =  1;
			dcopy(fftDimension,shiftedSignalTable,incx,fftTableIn,incy);
		}
		else if(feature & GABORWAVE)
		{
			/* from index table = [0] to [halfOfFFTDimension - 1] */
			incx = 1;
			incy = -1;
			dcopy(*halfOfFFTDimension,expTable+1,incx,shiftedGaussTable - *halfOfFFTDimension,incy);
		
			/* for N = halfOfFFTDimension */
			*shiftedGaussTable = 1.0;

			/* from index table = [halfOfFFTDimension + 1] to [fftDimension - 1] */
			incx = 1;
			incy = 1;
			dcopy(*halfOfFFTDimension,expTable+1,incx,shiftedGaussTable + 1,incy);
		
			if((fftDimension%2)==0)
				*(shiftedGaussTable + *halfOfFFTDimension + 1) = *(expTable + *halfOfFFTDimension + 1);
			//				*(shiftedGaussTable + fftDimension - 1) = *(expTable + *halfOfFFTDimension + 1);

			incz = 1;
			dxypz(fftDimension,shiftedSignalTable,incx,gaussTable,incy,fftTableIn,incz);
		}
	}
	else
	{
		if(feature & SINCOSWAVE)
		{
			/* from index table = [0] to [halfOfFFTDimension - 1] */
			incx =  1;
			incy =  1;
			dcopy(fftDimension,shiftedSignalTable,incx,fftTableIn,incy);
		}
		else
		{
			unsigned int firstRangeLength = 0;
			unsigned int n;
			unsigned char positionSide;
			unsigned int  firstStart  = 0;
			unsigned int  firstStop   = 0;
			unsigned int  secondStart = 0;
			unsigned int  secondStop  = 0;
		
			unsigned int startInterval = 0;
			unsigned int stopInterval  = 0;

			if((fftDimension%2)==0)
			{
				startInterval = intervalCenter - *halfOfFFTDimension;
				stopInterval  = intervalCenter + *halfOfFFTDimension + 1;
			}
			else
			{
				startInterval = intervalCenter - *halfOfFFTDimension;
				stopInterval  = intervalCenter + *halfOfFFTDimension;
			}
			
			findStartAndStopConditionsInLimitedRange(atomsPositionInExpandedSignal,&startInterval,&stopInterval,&firstStart,&firstStop,&secondStart,&secondStop);
			
			if(atomsPositionInExpandedSignal<startInterval)
				positionSide = LEFT_SIDE_POSITION;
			else if(atomsPositionInExpandedSignal>stopInterval)
				positionSide = RIGHT_SIDE_POSITION;
			else
			{
				const unsigned int sum = stopInterval + startInterval;

				if(2*(atomsPositionInExpandedSignal)<=sum)
					positionSide = LEFT_SIDE_POSITION;
				else 
					positionSide = RIGHT_SIDE_POSITION;
			}
			
			if(positionSide & LEFT_SIDE_POSITION)
			{
				if(firstStart<=firstStop)
				{
					n = firstStop - firstStart + 1;
					incx = 1;
					incy = -1;
					dcopy(n,expTable,incx,gaussTable,incy);
					firstRangeLength+=n;

					n = firstStop - firstStart;
					incx = 1;
					incy = 1;
					dcopy(n,expTable+1,incx,gaussTable + n + 1,incy);
					firstRangeLength+=n;
				}
				
				if(secondStart<=secondStop)
				{	
					n = secondStop - secondStart + 1;
					incx = 1;
					incy = 1;
					dcopy(n,expTable + secondStart,incx,gaussTable + firstRangeLength,incy);
				}
			}
			else
			{
				if(firstStart<=firstStop)
				{
				
					n = firstStop - firstStart + 1;
					incx = 1;
					incy = -1;
					dcopy(n,expTable,incx,gaussTable + fftDimension - 2*n + 1,incy);
					firstRangeLength+=n;

					n = firstStop - firstStart;
					incx = 1;
					incy = 1;
					dcopy(n,expTable+1,incx,gaussTable + fftDimension - n,incy);
					firstRangeLength+=n;
				}
			
				if(secondStart<=secondStop)
				{	
					n = secondStop - secondStart + 1;

					incx = 1;
					incy = -1;
					dcopy(n,expTable + secondStart,incx,gaussTable + fftDimension - firstRangeLength - n,incy);
				}
			}
			
			incx = 1;
			incy = 1;
			incz = 1;
				
			dxypz(fftDimension,shiftedSignalTable,incx,gaussTable,incy,fftTableIn,incz);
		}	
	}
}
	
static void findShiftedFFTInFrequency(const double *fftTableInA, 
									  fftw_complex *fftTableInB,
									  const double randomShiftInFrequency,
									  unsigned int fftDimension,
									  unsigned int halfOfFFTDimension,
									  unsigned int atomsFamilyPositionInExpandedSignal,
									  unsigned char mode)
									  
							   
{
	unsigned int i;
	unsigned int forwardIndex, backwardIndex;
	double basicSin, basicCos;
	double newSin, newCos;
	const double phaseCorrection = -M_2PI*randomShiftInFrequency/fftDimension;
	
	
	sincos(phaseCorrection,&basicSin,&basicCos);
	forwardIndex  = halfOfFFTDimension;
	backwardIndex = halfOfFFTDimension;
	
	newCos = 1.0;
	newSin = 0.0;

	fftTableInB[forwardIndex][0] = fftTableInA[forwardIndex];
	fftTableInB[forwardIndex][1] = 0.0;
	
	for(i=1;i<=halfOfFFTDimension;i++)
	{			
		findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
		forwardIndex++;
		backwardIndex--;

		fftTableInB[forwardIndex][0]  = fftTableInA[forwardIndex]*newCos;
		fftTableInB[forwardIndex][1]  = fftTableInA[forwardIndex]*newSin;
		fftTableInB[backwardIndex][0] = fftTableInA[backwardIndex]*newCos;
		fftTableInB[backwardIndex][1] = -fftTableInA[backwardIndex]*newSin;
	}

	if((fftDimension%2)==0)
	{
		findFastlySinCos(&basicSin,&basicCos,&newSin,&newCos);
	
		fftTableInB[fftDimension - 1][0] = fftTableInA[fftDimension - 1]*newCos;
		fftTableInB[fftDimension - 1][1] = fftTableInA[fftDimension - 1]*newSin;
	}	
}

static void findRSRCVariablesWithFFT(fftw_plan fftwPlan,
									 const double *signalTable,
									 const double *expTable,
									 double *gaussSignalTable,
									 double *fftTableInA,
									 fftw_complex *fftTableInB,
									 fftw_complex *fftTableOut,
									 double randomShiftInFrequency,
									 const unsigned int intervalCenter,
									 const unsigned int atomsPositionInExpandedSignal,
									 unsigned int fftDimension,
									 unsigned int *halfOfFFTDimension,
									 unsigned char feature,
									 unsigned char mode,
									 unsigned char typeOfDictionary,
									 unsigned int numberOfThreads)
{


	makeFFTTable(signalTable,
				 expTable,
				 gaussSignalTable,
				 fftTableInA,
				 intervalCenter,
				 atomsPositionInExpandedSignal,
				 fftDimension,
				 halfOfFFTDimension,
				 feature,
				 mode);
	
	if(typeOfDictionary & OCTAVE_STOCH)
	{
		findShiftedFFTInFrequency(fftTableInA,fftTableInB,randomShiftInFrequency,fftDimension,*halfOfFFTDimension,atomsPositionInExpandedSignal,mode);
	
		if(numberOfThreads>1)
			fftw_plan_with_nthreads(numberOfThreads);
			
		fftwPlan = fftw_plan_dft_1d(fftDimension,fftTableInB,fftTableOut,FFTW_FORWARD,FFTW_ESTIMATE);
		fftw_execute(fftwPlan);
	}
	else
	{
		if(numberOfThreads>1)
			fftw_plan_with_nthreads(numberOfThreads);

		fftwPlan = fftw_plan_dft_r2c_1d(fftDimension,fftTableInA,fftTableOut,FFTW_ESTIMATE);
		fftw_execute(fftwPlan);
	}
	
	fftw_destroy_plan(fftwPlan);
}

static void gaborDotProductUpdateFFT(Atom  *atomsFamily,
									 fftw_complex *fftTableOut,
									 double *bestModulusesTable,
									 unsigned int numberOfStepsInFrequencyAtParticularScale,
									 unsigned int multiple,
									 unsigned int atomsFamilyPositionInExpandedSignal,
									 unsigned int atomsFamilyPeriod,
									 double       randomShiftInFrequency,
									 unsigned int intervalCenter,
									 unsigned int halfOfFFTDimension,
									 unsigned int fftDimension,
									 unsigned short int channelNumber,
									 unsigned char mode,
									 unsigned char dictionaryType)										   
{
	Atom *atom = NULL;
	unsigned int frequencyIndex;
	unsigned int frequencyCounter;
	double sinOmega, cosOmega;
	double basicSinOmega, basicCosOmega;
	double frequencyCorrection;
	
	if(mode & FIRST_ITERATION)
		frequencyCorrection  = (M_2PI/atomsFamilyPeriod)*halfOfFFTDimension;
	else
	{
		if((intervalCenter - halfOfFFTDimension)>atomsFamilyPositionInExpandedSignal)
			frequencyCorrection  = -(M_2PI/atomsFamilyPeriod)*(intervalCenter - halfOfFFTDimension - atomsFamilyPositionInExpandedSignal);
		else
			frequencyCorrection  = (M_2PI/atomsFamilyPeriod)*(atomsFamilyPositionInExpandedSignal + halfOfFFTDimension - intervalCenter);

		//		printf(" freq correction: %lf %lf\n",frequencyCorrection,(M_2PI/atomsFamilyPeriod));
	}	

	atom = atomsFamily;
	
	if(mode & FIRST_ITERATION)
	{
		sincos(frequencyCorrection,&basicSinOmega,&basicCosOmega);

		sinOmega = basicSinOmega;
		cosOmega = basicCosOmega;
				
		for(frequencyIndex = multiple, frequencyCounter=0;frequencyCounter<numberOfStepsInFrequencyAtParticularScale;frequencyCounter++,frequencyIndex=(frequencyIndex+multiple))
		{	
			*(atom->RS + channelNumber) = (float)(-(fftTableOut[frequencyIndex][1]*cosOmega + fftTableOut[frequencyIndex][0]*sinOmega));
			*(atom->RC + channelNumber) = (float)(fftTableOut[frequencyIndex][0]*cosOmega - fftTableOut[frequencyIndex][1]*sinOmega);
		
			findFastlySinCos(&basicSinOmega,&basicCosOmega,&sinOmega,&cosOmega);

			atom++;			
		}	
	}
	else if(mode & NEXT_ITERATION)
	{
		double RS, RC;

		sincos(frequencyCorrection,&basicSinOmega,&basicCosOmega);
		
		sinOmega = basicSinOmega;
		cosOmega = basicCosOmega;
		
		for(frequencyIndex = multiple,frequencyCounter=0;frequencyCounter<numberOfStepsInFrequencyAtParticularScale;frequencyCounter++,frequencyIndex=(frequencyIndex+multiple))
		{							
			RS = *(atom->RS + channelNumber) - (*(bestModulusesTable + channelNumber))*(-(fftTableOut[frequencyIndex][1]*cosOmega + fftTableOut[frequencyIndex][0]*sinOmega));
			RC = *(atom->RC + channelNumber) - (*(bestModulusesTable + channelNumber))*(fftTableOut[frequencyIndex][0]*cosOmega - fftTableOut[frequencyIndex][1]*sinOmega);

			*(atom->RS + channelNumber) = (float)RS;
			*(atom->RC + channelNumber) = (float)RC;

			findFastlySinCos(&basicSinOmega,&basicCosOmega,&sinOmega,&cosOmega);

			atom++;
		}
	}
	else // (mode & MMP1_NEXT_ITERATION)
	{
		if(mode & MMP1_IGNORE_RS_RC)
		{
			for(frequencyIndex = multiple, frequencyCounter=0;frequencyCounter<numberOfStepsInFrequencyAtParticularScale;frequencyCounter++,frequencyIndex=(frequencyIndex+multiple))
			{
				*(atom->RS + channelNumber) = 0.0;
				*(atom->RC + channelNumber) = 0.0;
						
				atom++;
			}
		}
		else
		{
			sincos(frequencyCorrection,&basicSinOmega,&basicCosOmega);

			sinOmega = basicSinOmega;
			cosOmega = basicCosOmega;
			
			for(frequencyIndex = multiple, frequencyCounter=0;frequencyCounter<numberOfStepsInFrequencyAtParticularScale;frequencyCounter++,frequencyIndex=(frequencyIndex+multiple))
			{				
				*(atom->RS + channelNumber) = (float)(-(fftTableOut[frequencyIndex][1]*cosOmega + fftTableOut[frequencyIndex][0]*sinOmega));
				*(atom->RC + channelNumber) = (float)(fftTableOut[frequencyIndex][0]*cosOmega - fftTableOut[frequencyIndex][1]*sinOmega);

				findFastlySinCos(&basicSinOmega,&basicCosOmega,&sinOmega,&cosOmega);

				atom++;
			}
		}
	}
}

void findGaborDataDotProductFFT(const Dictionary *dictionary,
							    const MP5Parameters *mp5Parameters,
							    Atom  *atomsFamily,
							    const double *signalTable,
							    unsigned short int channelNumber,
							    unsigned char mode)
{
    const unsigned  int marginalDimension       = mp5Parameters->marginalDimension;

	const unsigned char      atomsFamilyFeature                  = atomsFamily->feature;
	const unsigned short int atomsFamilyScaleIndex               = atomsFamily->scaleIndex;
	const double             atomsFamilyScale                    = *(dictionary->tableOfScalesInOptimalDictionary + atomsFamilyScaleIndex);
	const unsigned int atomsFamilyPositionInExpandedSignal = atomsFamily->position + marginalDimension;
	const unsigned int atomsFamilyPeriod                   = *(dictionary->tableOfPeriodsInOptimalDictionary   + atomsFamily->scaleIndex);
	const unsigned int numberOfStepsInFrequencyAtParticularScale = *(dictionary->numberOfStepsInFrequencyAtParticularScale + atomsFamily->scaleIndex);
		  unsigned int intervalCenter = 0; 
		  unsigned int intervalRange  = 0;
		  
	unsigned int halfOfFFTDimension = 0;
    double       *fftTableInA = mp5Parameters->fftTableInA;
    fftw_complex *fftTableInB = mp5Parameters->fftTableInB;
    fftw_complex *fftTableOut = mp5Parameters->fftTableOut;
	fftw_plan    fftwPlan     = mp5Parameters->fftwPlan;

    double *bestModulusesTable = mp5Parameters->bestModulusesTable;

	const  double randomShiftInFrequency = atomsFamily->randomShiftInFrequency;
	double *expTable         = *(mp5Parameters->expTable + atomsFamilyScaleIndex);
	double *gaussSignalTable = mp5Parameters->gaussSignalTable;
	double K = 0.0;

	unsigned int factor   = 0.0;
	unsigned int multiple = 1;
	unsigned int fftDimension = *(dictionary->tableOfPeriodsInOptimalDictionary + atomsFamily->scaleIndex);

	if((mode & FIRST_ITERATION))
	{
		K = LOG_EPS_DOT_PRODUCT;
		
		if(atomsFamilyFeature & SINCOSWAVE)
		{
			fftDimension   = *(dictionary->tableOfPeriodsInOptimalDictionary + dictionary->numberOfStepsInScale - 1);
			intervalCenter = atomsFamilyPositionInExpandedSignal;
		}
		else
		{	
			findGaussNonGaussInterval(atomsFamily->position,atomsFamilyScale,marginalDimension,&intervalCenter,&intervalRange);
	
			factor   = (unsigned int)((2.0*intervalRange)/fftDimension + 0.5);

			if(factor>1)
			{
				multiple = factor;
				fftDimension = multiple*fftDimension;
			}
		}
	}
	else
	{	
		const unsigned char previousAtomFeature = (mp5Parameters->previousAtom)->feature;

		if(atomsFamilyFeature & SINCOSWAVE)
		{ // dot product of sin/cos with sin/cos  must be estimated in the full range

			K = LOG_EPS_DOT_PRODUCT;
			fftDimension   = *(dictionary->tableOfPeriodsInOptimalDictionary + dictionary->numberOfStepsInScale - 1);
			intervalCenter = atomsFamilyPositionInExpandedSignal;
		}
		else if(previousAtomFeature & SINCOSWAVE)
		{ // dot product of sin/cos with gabor or gauss can be estimated in the shorter range		

			K = LOG_EPS_DOT_PRODUCT;
	
			findGaussNonGaussInterval(atomsFamily->position,atomsFamilyScale,marginalDimension,&intervalCenter,&intervalRange);
			factor   = (unsigned int)((2.0*intervalRange)/fftDimension + 0.5);			   

			if(factor>1)
			{
				multiple = factor;
				fftDimension = multiple*fftDimension;
			}
			
			intervalCenter = atomsFamilyPositionInExpandedSignal;
		}
		else
		{
			// dot product of gabor with gabor or gauss

			K = findGaussGaussInterval(dictionary,
									   atomsFamily,
									   mp5Parameters->previousAtom, 
									   marginalDimension,
									   &intervalCenter, // findGaussGaussInterval return intervalCenter +  marginalDimension
									   &intervalRange);
									   
			factor   = (unsigned int)((2.0*intervalRange)/fftDimension + 0.5);			   

			if(factor>1)
			{
				multiple = factor;
				fftDimension = multiple*fftDimension;
			}
		}
	}	

	if(K>=LOG_EPS_DOT_PRODUCT)
	{
		findRSRCVariablesWithFFT(fftwPlan,
								 signalTable,
		  					     expTable,
								 gaussSignalTable,
								 fftTableInA,
								 fftTableInB,
								 fftTableOut,
								 randomShiftInFrequency,
								 intervalCenter,
								 atomsFamilyPositionInExpandedSignal,
								 fftDimension,
								 &halfOfFFTDimension,
								 atomsFamily->feature,
								 mode,
								 dictionary->typeOfDictionary,
								 mp5Parameters->numberOfThreads);
		
		gaborDotProductUpdateFFT(atomsFamily,
								 fftTableOut,
								 bestModulusesTable,
								 numberOfStepsInFrequencyAtParticularScale,
								 multiple,
								 atomsFamilyPositionInExpandedSignal,
								 atomsFamilyPeriod,
								 randomShiftInFrequency,
								 intervalCenter,
								 halfOfFFTDimension,
								 fftDimension,
								 channelNumber,
								 mode,
								 dictionary->typeOfDictionary);					 
	}
	else
	{	
		if(mode & MMP1_NEXT_ITERATION)
		{
			gaborDotProductUpdateFFT(atomsFamily,
									fftTableOut,
									bestModulusesTable,
									numberOfStepsInFrequencyAtParticularScale,
									multiple,
									atomsFamilyPositionInExpandedSignal,
									atomsFamilyPeriod,
									randomShiftInFrequency,
									intervalCenter,
									halfOfFFTDimension,
									fftDimension,
									channelNumber,
									mode | MMP1_IGNORE_RS_RC,
									dictionary->typeOfDictionary);
		}
	}
}

STATUS findUnknowPhaseDI(Atom *atom, double *modulus, unsigned short int channelNumber)
{

    const double KS = atom->KS;
    const double KC = atom->KC;
    const double KM = atom->KM; 
    const double RS = *(atom->RS + channelNumber);
    const double RC = *(atom->RC + channelNumber);
	double sinSqr;

    double W1, W2;

    double phase     = 0.0;
    double amplitude = 0.0;

    double sinPhase, cosPhase;

    if((!(atom->feature & DIRACDELTA)) && (!(atom->feature & GAUSSFUNCTION)))
    {
		if(solveSystemOfEquestions(&KC,&KM,&KM,&KS,&W1,&W2,&RC,&RS) == ERROR)
		{
			fprintf(stderr,"\n\nAtom with following parameters:");
			fprintf(stderr,"\n POSITION:  %u       \
						    \n RIFLING:   %u       \
		                    \n SCALE:     %hu\n",atom->position,atom->rifling,atom->scaleIndex);							
			fprintf(stderr," Is being marked as a incorrect atom\n");
			fprintf(stderr," Because one can not estimate phase of this atom \n");
    	    fprintf(stderr," according to Dobieslaw Ircha algorithms \n");
			fprintf(stderr," typ atomu: %hu %hu %hu\n",atom->feature & SINCOSWAVE,atom->feature & GAUSSFUNCTION,atom->feature & GABORWAVE);
			fprintf(stderr," KS KC: %1.12lf %1.16lf %1.16lf\n",KS,KC,KM);
			fprintf(stderr," RS RC: %1.12lf %1.16lf\n",RS,RC);
			*(atom->phase + channelNumber) = 0.0;
			fflush(stderr);
			*modulus = 0.0;
	
			atom->feature|= INCORRECTGABOR;
	
			return ERROR;
		}
		else
			phase = atan2(-W2,W1);

	}
    else
		phase = (RC>0 ? -M_2PI : M_PI); // the sufficient condition should be (RC>0 ? 0 : M_PI), but was changed due to mmp1 optimalization
										// from the mathematical point of view, there is no problem, because cos(-2PI) = cos(0)

    sincos(phase,&sinPhase,&cosPhase);

	sinSqr = sinPhase*sinPhase;
	
    const double sinPart    = KS*sinSqr;
    const double cosPart    = KC*(1 - sinSqr);
    const double sinCosPart = 2.0*KM*sinPhase*cosPhase;

    amplitude = sqrt((sinPart + cosPart ) - sinCosPart);
    *modulus  = (RC*cosPhase - RS*sinPhase)/amplitude;

    *(atom->phase + channelNumber) = (float)phase;

    return SUCCESS;
}

STATUS findUnknowPhaseAM(Atom *atom, double *modulusTable, unsigned short int numberOfAnalysedChannels)
{
	if((!(atom->feature & DIRACDELTA)) && (!(atom->feature & GAUSSFUNCTION)))
	{   
		unsigned short int channel;

        double tmpSumX = 0.0;
        double tmpSumY = 0.0;
        double sinPhase, cosPhase;
        double amplitude;

        const double KS = atom->KS;
        const double KC = atom->KC;
        const double KM = atom->KM;
		double modulus;
        double phase;

		for (channel=0;channel<numberOfAnalysedChannels;channel++) 
		{
			modulus =  *(modulusTable + channel);
			phase   =  *(atom->phase + channel);
			tmpSumX += modulus*sin(2.0 * phase);    
			tmpSumY += modulus*cos(2.0 * phase);    
		} 

		const double tmpPhase = atan2(tmpSumX,tmpSumY)/2.0;

		for(channel=0;channel<numberOfAnalysedChannels;channel++) 
		{
			if ((*(atom->RS + channel))>0.0)
			{
				sincos(tmpPhase,&sinPhase,&cosPhase);
			}
			else
			{
				sincos(tmpPhase+M_PI,&sinPhase,&cosPhase);
			}
	    
			amplitude = sqrt(KS*(sinPhase*sinPhase) + KC*(cosPhase*cosPhase) - 2.0*KM*sinPhase*cosPhase);

			*(modulusTable + channel) = ((*(atom->RC + channel))* cosPhase - (*(atom->RS + channel))*sinPhase)/amplitude;
	    
			if ((*(atom->RS + channel))>0.0)
			{
				*(atom->phase + channel) = (float)tmpPhase;
			}
			else
			{
				*(atom->phase + channel) = (float)tmpPhase+M_PI;
			}
	    }
    }
    return SUCCESS;
}

double findLambda(double energyOfResidue)
{
	static double oldEnergyOfResidue = 0.0;
	double lambda = 0.0;

	lambda = sqrt(1.0 - energyOfResidue/oldEnergyOfResidue);
	oldEnergyOfResidue = energyOfResidue;
	
	return lambda;
}
