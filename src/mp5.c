/***********************f****************************************************
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
#include<string.h>
#include<strings.h>
#include<unistd.h>
#include"include/def.h"
#include"include/dic.h"
#include"include/gabor.h"
#include"include/matrix.h"
#include"include/mp5.h"
#include"include/queue.h"
#include"include/tools.h"
#include"include/types.h"
#include"include/vector.h"

#ifdef __MINGW32__
       #define bzero(ptr,size) memset (ptr, 0, size);
       #define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

#define findStartResample(sample,rifling,period)         ((sample*rifling)%period);
#define findNextResample(currentResample,rifling,period) ((currentResample+rifling)%period);


static void findStartAndStopConditionsInFullRange(unsigned short int position,
						  unsigned short int dimOffset,
						  unsigned       int dimExpTable,
						  unsigned int *firstStart,
						  unsigned int *firstStop,
						  unsigned int *secondStart,
						  unsigned int *secondStop,
						  unsigned char feature)
{
/* In this implementation of MP algorithm calculation are performed in the range
 from 0 to N-1, where N is an offset's length.
 However, one can make MP estimation faster, if the fatures of gabors function will be
 used. 'Cosine-Gabor' is a even function and 'Sine-Gabor' is a odd function, therefore
 we don't have to evaluate dot products for the whole offset. Morover, Gabor function deseaper
 very fast, so we can estimate dot products only there, where value of Gabor Function is grater then
 some EPS. This function is looking for the beginning and the end of offset, where calculation
 should be permpormed. */

    const unsigned int startPosition = dimOffset   + position;;
    const unsigned int endOfTable    = dimExpTable - position;
    
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

static void findStartAndStopConditionsInLimitedRange(unsigned short int position,
						     unsigned short int dimOffset,
						     unsigned       int dimExpand,
						     long  int *startRange,
						     long  int *stopRange,
						     unsigned int *firstStart,
						     unsigned int *firstStop,
						     unsigned int *secondStart,
						     unsigned int *secondStop)
{
/* In this implementation of MP algorithm calculation are performed in the range
 from 0 to N-1, where N is an offset's length.
 However, one can make MP estimation faster, if the fatures of gabors function will be
 used. 'Cosine-Gabor' is a even function and 'Sine-Gabor' is a odd function, therefore
 we don't have to evaluate dot products for the whole offset. Morover, Gabor function deseaper
 very fast, so we can estimate dot products only there, where value of Gabor Function is grater then
 some EPS. This function is looking for the beginning and the end of offset, where calculation
 should be permpormed. */

    const long int startPosition = dimOffset + position;
    long int newStartRange, newStopRange, tmpValue;

    *startRange = (*startRange<0                 ? 0             : *startRange);
    *startRange = (*startRange>=(int)(dimExpand) ? dimExpand - 1 : *startRange);

    *stopRange  = (*stopRange<0                 ? 0             : *stopRange);
    *stopRange  = (*stopRange>=(int)(dimExpand) ? dimExpand - 1 : *stopRange);

    const unsigned long int tmpStartRange = *startRange;
    const unsigned long int tmpStopRange  = *stopRange;

    newStartRange = ((startPosition < *startRange) ? (*startRange - startPosition) : (startPosition - *startRange));
    newStopRange  = ((startPosition < *stopRange)  ? (*stopRange  - startPosition) : (startPosition - *stopRange));

    newStartRange = (newStartRange<0 ? -newStartRange : newStartRange);
    newStopRange  = (newStopRange<0  ? -newStopRange  : newStopRange);

    if(newStartRange>newStopRange)
    {
	tmpValue = newStartRange;
	newStartRange = newStopRange;
	newStopRange  = tmpValue;
    }

    if((*startRange<startPosition && *stopRange<startPosition) || (*startRange>startPosition && *stopRange>startPosition))
    {
	
	*firstStart = 1; /* these condition seem to be very strange, but thanks to them */
	*firstStop  = 0; /* some loops will not start, what is needed by us :-) */

	*secondStart = newStartRange;
	*secondStop  = newStopRange;
    }
    if(((*startRange<=startPosition) && (*stopRange>=startPosition)) || ((*startRange>=startPosition) && (*stopRange<=startPosition)))
    {
    
	*firstStart = 0;
	*firstStop  = newStartRange;

	*secondStart = newStartRange + 1; // when secondStart = secondStop condition seems to be strange, but thanks to iy
	*secondStop  = newStopRange;      // some loops will not start
    }

    *startRange = tmpStartRange;
    *stopRange  = tmpStopRange;

}

static void findRSRCVariables(const MP5Parameters *mp5Parameters,
			      const GaborDictionary *gaborDictionary,
			      const Gabor *gabor,
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

    const unsigned short int scaleIndex = gabor->scaleIndex;
    const unsigned       int rifling    = gabor->rifling;
    const unsigned short int position   = gabor->position;
    const unsigned char      feature    = gabor->feature; 
	  unsigned char      located;

    const unsigned short int dimOffset     = mp5Parameters->dimOffset;
    const unsigned       int startPosition = dimOffset + position;

    const double *sinTable = *(mp5Parameters->sinTable + scaleIndex);
    const double *cosTable = *(mp5Parameters->cosTable + scaleIndex);
    const double *expTable = *(mp5Parameters->expTable + scaleIndex);

    const unsigned int commonPeriod = *(gaborDictionary->tableOfPeriodsInOptimalDictionary + scaleIndex);

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
    else if(feature & FFTWAVE)
    {  
	*RC = *ptrDataTable;
        
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

	(*RC)-= 2.0*(*ptrDataTable);
    }
    else
    {
	double expVal;

	*RC = *ptrDataTable;

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
	(*RC)-= 2.0*(*ptrDataTable);
    }
}

static void makeOneSinCosTable(double omega, double *sinTable, double *cosTable, unsigned int period)
{

    unsigned int sample;
    double sinOmega, cosOmega;
    double oldCos, oldSin;
    double newCos;

    sincos(omega,&sinOmega,&cosOmega);

    oldSin = 0.0;
    oldCos = 1.0;

    *sinTable = oldSin;
    *cosTable = oldCos;

    for(sample = 1;sample<period;sample++)
    {
	newCos = oldCos*cosOmega - oldSin*sinOmega;

	*(sinTable + sample) = (oldSin = oldCos*sinOmega + oldSin*cosOmega);
	*(cosTable + sample) = (oldCos = newCos);
    }
}

static void makeOneExpTable(double alpha, double *expTable, unsigned int dimExpTable)
{

    unsigned int sample;
    double oldExpTable;
    double expTwoAlpha, expAlphaConst, expTwoAlphaConst;

    expAlphaConst    = exp(-alpha);
    expTwoAlphaConst = exp(-2.0*alpha);
    expTwoAlpha      = expTwoAlphaConst;

    *expTable = 1.0;

    oldExpTable = expAlphaConst;
    *(expTable + 1) = oldExpTable;

    for(sample = 2;sample<dimExpTable;sample++)
    {
	oldExpTable = oldExpTable*expTwoAlpha*expAlphaConst;
	expTwoAlpha = expTwoAlpha*expTwoAlphaConst;
	*(expTable + sample) = oldExpTable;
    }
}


static void findKSKCKMVariables(const MP5Parameters *mp5Parameters,
				  const GaborDictionary *gaborDictionary,   
				  Gabor *gabor,
				  unsigned int firstStart,
				  unsigned int firstStop,
				  unsigned int secondStart,
				  unsigned int secondStop)
{

    unsigned int sample;
    unsigned int resample;
    double sinVal, cosVal;
    
    double KS = 0.0;
    double KC = 0.0;
    double KM = 0.0;

    unsigned short int scaleIndex = gabor->scaleIndex;
    unsigned       int rifling    = gabor->rifling;

    const double *sinTable = *(mp5Parameters->sinTable + scaleIndex);
    const double *cosTable = *(mp5Parameters->cosTable + scaleIndex);
    const double *expTable = *(mp5Parameters->expTable + scaleIndex);

    const unsigned int commonPeriod = *(gaborDictionary->tableOfPeriodsInOptimalDictionary + scaleIndex);

    if(gabor->feature & DIRACDELTA)
    {
	KS = 0.0;
	KC = 1.0;
        KM = 0.0;
    }
    else if(gabor->feature & FFTWAVE)
    {
	resample = findStartResample(firstStart,rifling,commonPeriod);

	for(sample=firstStart;sample<=firstStop;sample++)
	{
	    sinVal = *(sinTable + resample);
	    cosVal = *(cosTable + resample);

	    KS+= sinVal*sinVal;
	    KC+= cosVal*cosVal;
	    /* *KM = 0, because of symtrical sumation of elements of gabors */

	    resample = findNextResample(resample,rifling,commonPeriod);
	}

	KS*= 2.0;
	KC = 2.0*KC - 1.0;

	resample = findStartResample(secondStart,rifling,commonPeriod);

	for(sample=secondStart;sample<=secondStop;sample++)
	{

	    sinVal   = *(sinTable + resample);
	    cosVal   = *(cosTable + resample);

	    KS+= sinVal*sinVal;
	    KC+= cosVal*cosVal;
	    KM+= cosVal*sinVal;

	    resample = findNextResample(resample,rifling,commonPeriod);
	}
    }
    else
    {
	double expVal;

	resample = findStartResample(firstStart,rifling,commonPeriod);

	for(sample=firstStart;sample<=firstStop;sample++)
	{
	    sinVal = *(sinTable + resample);
	    cosVal = *(cosTable + resample);
	    expVal = *(expTable + sample);

	    KS+= (sinVal*expVal)*(sinVal*expVal);
	    KC+= (cosVal*expVal)*(cosVal*expVal);

	    /* *KM = 0, because of symtrical sumation of elements of gabors */
    
	    resample = findNextResample(resample,rifling,commonPeriod);
	}

	KS*= 2.0;
	KC = 2.0*KC - 1.0;

	resample = findStartResample(secondStart,rifling,commonPeriod);

	for(sample=secondStart;sample<=secondStop;sample++)
	{
	    sinVal   = *(sinTable + resample);
	    cosVal   = *(cosTable + resample);
	    expVal   = *(expTable + sample);

	    KS+= (sinVal*expVal)*(sinVal*expVal);
	    KC+= (cosVal*expVal)*(cosVal*expVal);
	    KM+= (cosVal*expVal)*(sinVal*expVal);

	    resample = findNextResample(resample,rifling,commonPeriod);
	}

    }

    gabor->KS = (float)KS;
    gabor->KC = (float)KC;
    gabor->KM = (float)KM;

}

void setNumberOfAnalysedChannelsAndNumberOfResultsFiles(MP5Parameters *mp5Parameters, DataParameters *dataParameters)
{
    if(mp5Parameters->MPType & SMP)
    {
	mp5Parameters->numberOfAnalysedChannels = 1;
	dataParameters->numberOfResultsFiles = dataParameters->numberOfChosenChannels;
    }
    else if ((mp5Parameters->MPType & MMP1) ||  (mp5Parameters->MPType & MMP2) ||  (mp5Parameters->MPType & MMP3))
    {
	mp5Parameters->numberOfAnalysedChannels = dataParameters->numberOfChosenChannels;
	dataParameters->numberOfResultsFiles = 1;
    }
}

void setMP5Parameters(MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary)
{

    unsigned short int numberOfStepsInScale = gaborDictionary->numberOfStepsInScale;
    unsigned       int dimExpTable          = mp5Parameters->dimExpTable;
    unsigned       int dimExpand            = mp5Parameters->dimExpand;

    unsigned int *tableOfPeriodsInOptimalDictionary = gaborDictionary->tableOfPeriodsInOptimalDictionary;

    unsigned short int numberOfAnalysedChannels = mp5Parameters->numberOfAnalysedChannels;

    mp5Parameters->fitted             = createQueue();
    mp5Parameters->bestModulusesTable = dVectorAllocate(numberOfAnalysedChannels);

    mp5Parameters->signalEnergyInEachChannel  = dVectorAllocate(numberOfAnalysedChannels);
    mp5Parameters->residueEnergyInEachChannel = dVectorAllocate(numberOfAnalysedChannels);

    mp5Parameters->meanSignalTable  = dVectorAllocate(dimExpand);
    mp5Parameters->meanResidueTable = dVectorAllocate(dimExpand);

    mp5Parameters->prevGaborTable   = dMatrixAllocate(numberOfAnalysedChannels,dimExpand);

    /* sin/cos function are build for each scale */

    mp5Parameters->sinTable = dVariableMatrixAllocate(numberOfStepsInScale,tableOfPeriodsInOptimalDictionary);
    mp5Parameters->cosTable = dVariableMatrixAllocate(numberOfStepsInScale,tableOfPeriodsInOptimalDictionary);
    mp5Parameters->expTable = dMatrixAllocate(numberOfStepsInScale,dimExpTable);

    mp5Parameters->sinGaborTable = dVectorAllocate(dimExpand);
    mp5Parameters->cosGaborTable = dVectorAllocate(dimExpand);

    mp5Parameters->memoryAllocated = TRUE;

}

void freeMP5Parameters(MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary)
{
    if(mp5Parameters->memoryAllocated == TRUE)
    {
	dMatrixFree(mp5Parameters->prevGaborTable,mp5Parameters->numberOfAnalysedChannels);
	dVectorFree(mp5Parameters->sinGaborTable);
	dVectorFree(mp5Parameters->cosGaborTable);

	dVectorFree(mp5Parameters->signalEnergyInEachChannel);
	dVectorFree(mp5Parameters->residueEnergyInEachChannel);

	dVectorFree(mp5Parameters->meanSignalTable);
	dVectorFree(mp5Parameters->meanResidueTable);

	dVariableMatrixFree(mp5Parameters->sinTable,gaborDictionary->numberOfStepsInScale);
	dVariableMatrixFree(mp5Parameters->cosTable,gaborDictionary->numberOfStepsInScale);
	dMatrixFree(mp5Parameters->expTable,gaborDictionary->numberOfStepsInScale);

	dVectorFree(mp5Parameters->bestModulusesTable);

	freeQueue(mp5Parameters->fitted,(void (*)(void *))freeGabor);

	mp5Parameters->memoryAllocated = FALSE;
    }
}

/* prepare Sin, Cos and Exp Tables */

void makeSinCosExpTable(MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary)
{

    unsigned int scaleIndex;
    double  *sinTable, *cosTable, *expTable;
    const unsigned short int numberOfStepsInScale = gaborDictionary->numberOfStepsInScale;
    const unsigned       int dimExpTable = mp5Parameters->dimExpTable;
    double  omega, alpha;
    unsigned short int scale;
    unsigned       int period;

    const unsigned short int *tableOfScalesInOptimalDictionary  = gaborDictionary->tableOfScalesInOptimalDictionary;
    const unsigned       int *tableOfPeriodsInOptimalDictionary = gaborDictionary->tableOfPeriodsInOptimalDictionary;
    const double  *tableOfFrequenciesInOptimalDictionary = gaborDictionary->tableOfFrequenciesInOptimalDictionary;

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

	makeOneExpTable(alpha,expTable,dimExpTable);
    }

    printf(" END SIN, COS, EXP, DIRAC TABLES GENERAITING \n\n");
}


void makeSinCosGaborTable(MP5Parameters *mp5Parameters, 
                          const GaborDictionary *gaborDictionary,
                          const Gabor  *gabor)
{

    unsigned       int sample;
    unsigned       int resample;

    const unsigned short int scaleIndex = gabor->scaleIndex;
    const unsigned       int rifling    = gabor->rifling;
    const unsigned short int position   = gabor->position;
    const unsigned char      feature    = gabor->feature;
    unsigned char            located;

    const unsigned short int dimOffset   = mp5Parameters->dimOffset;
    const unsigned       int dimExpTable = mp5Parameters->dimExpTable;
    const unsigned       int dimExpand   = mp5Parameters->dimExpand;

    const unsigned int startPosition = dimOffset + position;
    unsigned int firstStart  = 0;
    unsigned int firstStop   = 0;
    unsigned int secondStart = 0;
    unsigned int secondStop  = 0;

    const unsigned int commonPeriod = *(gaborDictionary->tableOfPeriodsInOptimalDictionary + scaleIndex);

    const double *sinTable = *(mp5Parameters->sinTable + scaleIndex);
    const double *cosTable = *(mp5Parameters->cosTable + scaleIndex);
    const double *expTable = *(mp5Parameters->expTable + scaleIndex);

    double *sinGaborTable = mp5Parameters->sinGaborTable;
    double *cosGaborTable = mp5Parameters->cosGaborTable;

    double *ptrSinGaborTable, *ptrCosGaborTable;

    findStartAndStopConditionsInFullRange(position,dimOffset,dimExpTable,&firstStart,&firstStop,&secondStart,&secondStop,gabor->feature);

    ptrSinGaborTable = sinGaborTable + startPosition;
    ptrCosGaborTable = cosGaborTable + startPosition;

    located = (unsigned char)((gabor->feature & LEFT_SIDE_POSITION_IN_OFFSET) ? LEFT_SIDE_POSITION : RIGHT_SIDE_POSITION);

    if(feature & DIRACDELTA)
    {	
	dSetVectorZero(sinGaborTable,dimExpand);
        dSetVectorZero(cosGaborTable,dimExpand);
        *ptrCosGaborTable = 1.0;
    }
    else if(feature & FFTWAVE)
    {
	resample = findStartResample(firstStart,rifling,commonPeriod);

        for(sample=firstStart;sample<=firstStop;sample++)
        {
	    *(ptrSinGaborTable  + sample) = (*(sinTable + resample));
	    *(ptrSinGaborTable  - sample) = -(*(ptrSinGaborTable + sample));
	    *(ptrCosGaborTable  - sample) = *(ptrCosGaborTable + sample) = (*(cosTable + resample));

	    resample = findNextResample(resample,rifling,commonPeriod);
	}

        if(located & LEFT_SIDE_POSITION)
	{
	    resample = findStartResample(secondStart,rifling,commonPeriod);
                
	    for(sample=secondStart;sample<=secondStop;sample++)
	    {
		*(ptrSinGaborTable  + sample) = (*(sinTable + resample));
		*(ptrCosGaborTable  + sample) = (*(cosTable + resample));

		resample = findNextResample(resample,rifling,commonPeriod);
	    }
	}
	else if(located & RIGHT_SIDE_POSITION)
	{
	    resample = findStartResample(secondStart,rifling,commonPeriod);

	    for(sample=secondStart;sample<=secondStop;sample++)
	    {
		*(ptrSinGaborTable - sample) = -(*(sinTable + resample));
		*(ptrCosGaborTable - sample) = (*(cosTable + resample));

		resample = findNextResample(resample,rifling,commonPeriod);
	    }
	}
    }
    else
    {
	double expVal;

	resample = findStartResample(firstStart,rifling,commonPeriod);

	for(sample=firstStart;sample<=firstStop;sample++)
	{
	    expVal = *(expTable + sample);
	    *(ptrSinGaborTable  + sample) = (*(sinTable + resample))*expVal;
	    *(ptrSinGaborTable  - sample) = -(*(ptrSinGaborTable + sample));
	    *(ptrCosGaborTable  - sample) = *(ptrCosGaborTable + sample) = (*(cosTable + resample))*expVal;

	    resample = findNextResample(resample,rifling,commonPeriod);
	}

	if(located & LEFT_SIDE_POSITION)
	{	          
	    resample = findStartResample(secondStart,rifling,commonPeriod);

	    for(sample=secondStart;sample<=secondStop;sample++)
	    {
		expVal = *(expTable + sample);
		*(ptrSinGaborTable  + sample) = (*(sinTable + resample))*expVal;
		*(ptrCosGaborTable  + sample) = (*(cosTable + resample))*expVal;

		resample = findNextResample(resample,rifling,commonPeriod);
	    }
	}
	else if(located & RIGHT_SIDE_POSITION)
	{
	    resample = findStartResample(secondStart,rifling,commonPeriod);

	    for(sample=secondStart;sample<=secondStop;sample++)
	    {
		expVal = *(expTable + sample);
		*(ptrSinGaborTable - sample) = -(*(sinTable + resample))*expVal;
		*(ptrCosGaborTable - sample) = (*(cosTable + resample))*expVal;

		resample = findNextResample(resample,rifling,commonPeriod);
	    }
	}
    }
}

void normSinCosGaborTable(const MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary, Gabor *gabor)
{

    const unsigned short int position    = gabor->position;
    const unsigned short int dimOffset   = mp5Parameters->dimOffset;
    const unsigned       int dimExpTable = mp5Parameters->dimExpTable;

    unsigned       int firstStart  = 0;
    unsigned       int secondStart = 0;
    unsigned       int firstStop   = 0;
    unsigned       int secondStop  = 0;

    findStartAndStopConditionsInFullRange(position,dimOffset,dimExpTable,&firstStart,&firstStop,&secondStart,&secondStop,gabor->feature);

    findKSKCKMVariables(mp5Parameters,gaborDictionary,gabor,firstStart,firstStop,secondStart,secondStop);
                        
}

void makeGaborTable(MP5Parameters *mp5Parameters, const Gabor *gabor, unsigned short int channelNumber)
{
    unsigned int sample;
    const unsigned int dimExpand = mp5Parameters->dimExpand;

    const double *cosGaborTable = mp5Parameters->cosGaborTable;
    const double *sinGaborTable = mp5Parameters->sinGaborTable;
    
    double *prevGaborTable = *(mp5Parameters->prevGaborTable + channelNumber);

    double KS = gabor->KS;
    double KC = gabor->KC;
    double KM = gabor->KM;

    const double phase = *(gabor->phase + channelNumber);

    double amplitude;
    double sinPhase, cosPhase;

    sincos(phase,&sinPhase,&cosPhase);

    const double sinPart    = KS*(sinPhase*sinPhase);
    const double cosPart    = KC*(cosPhase*cosPhase);
    const double sinCosPart = 2.0*KM*sinPhase*cosPhase;

    amplitude = sqrt((sinPart + cosPart) - sinCosPart);

    for(sample=0;sample<dimExpand;sample++)
	*(prevGaborTable + sample) = ((*(cosGaborTable + sample))*cosPhase - (*(sinGaborTable + sample))*sinPhase)/amplitude;
}

double findSignalEnergy(const double *signalTable, unsigned short int dimOffset)
{
    unsigned short int sample;
    const double *ptrSignalTable = signalTable + dimOffset;
    double energy = 0.0;

    for(sample=0;sample<dimOffset;sample++)
	energy+= (*(ptrSignalTable + sample))*(*(ptrSignalTable + sample));

    return energy;
}

void findResidue(double *residueTable, const double *gaborTable, double modulus, unsigned int dimExpand)
{
    unsigned int sample;
    
    for(sample=0;sample<dimExpand;sample++)
	*(residueTable + sample) = *(residueTable + sample) - modulus*(*(gaborTable + sample));
}

/* Implementation of Dobieslaw Ircha Optimal Phase Algorithm */


void findGaborDataDotProduct(const MP5Parameters *mp5Parameters, const GaborDictionary *gaborDictionary,
			     Gabor *currentGabor,
			     const double *dataTable,
			     unsigned short int channelNumber,
			     unsigned char mode)
{

    const unsigned short int dimOffset   = mp5Parameters->dimOffset;
    const unsigned       int dimExpTable = mp5Parameters->dimExpTable;
    const unsigned       int dimExpand   = mp5Parameters->dimExpand;

    const unsigned short int currentPosition = currentGabor->position;
    const unsigned short int currentScale    = *(gaborDictionary->tableOfScalesInOptimalDictionary + currentGabor->scaleIndex);;

    unsigned       int firstStart  = 0;
    unsigned       int firstStop   = 0;
    unsigned       int secondStart = 0;
    unsigned       int secondStop  = 0;

    double RS, RC;
    double *bestModulusesTable = mp5Parameters->bestModulusesTable;

    Gabor *previousGabor = mp5Parameters->previousGabor;

    RS = RC = 0.0;

    if((mode & FIRST_ITERATION))
    {
	if(!(currentGabor->feature & DIRACDELTA))
	    findStartAndStopConditionsInFullRange(currentPosition,dimOffset,dimExpTable,&firstStart,&firstStop,&secondStart,&secondStop,currentGabor->feature);

	findRSRCVariables(mp5Parameters,gaborDictionary,currentGabor,dataTable,
			  firstStart,firstStop,secondStart,secondStop,
			  FULL_RANGE,&RS,&RC);

	*(currentGabor->RS + channelNumber) = (float)RS;
	*(currentGabor->RC + channelNumber) = (float)RC;
    }
    else if(mode & NEXT_ITERATION)
    {
	if(!(currentGabor->feature & DIRACDELTA))
	{
	    const unsigned short int previousPosition = previousGabor->position;
	    const unsigned short int previousScale    = *(gaborDictionary->tableOfScalesInOptimalDictionary + previousGabor->scaleIndex);

	    long int startRange;
	    long int stopRange;

	    const double sqrPreviousScale = previousScale*previousScale;
	    const double sqrCurrentScale  = currentScale*currentScale;

	    const double newScale = sqrPreviousScale + sqrCurrentScale;

	    double K = -((M_PI*((previousPosition - currentPosition)*(previousPosition - currentPosition)))/newScale);

	    if(K>mp5Parameters->LOG_EPS_DOT_PRODUCT)
	    {
		const unsigned long int  POSITION = dimOffset + (unsigned int)(0.5 + (currentPosition*sqrPreviousScale + previousPosition*sqrCurrentScale)/newScale);
		const unsigned long int  SCALE    = (unsigned int)(1.5 + sqrt((K - mp5Parameters->LOG_EPS_DOT_PRODUCT)/(M_PI*newScale/(sqrCurrentScale*sqrPreviousScale))));

		startRange = POSITION - SCALE;
		stopRange  = POSITION + SCALE;

		findStartAndStopConditionsInLimitedRange(currentPosition,dimOffset,dimExpand,&startRange,&stopRange,&firstStart,&firstStop,&secondStart,&secondStop);

		if((currentPosition + dimOffset)<startRange)
	    	    currentGabor->feature|=LEFT_SIDE_POSITION_IN_RANGE;
		else if((currentPosition + dimOffset)>stopRange)
	    	    currentGabor->feature&= (~LEFT_SIDE_POSITION_IN_RANGE);
		else
		{
		    const unsigned int sum = stopRange + startRange;

		    if(2*(currentPosition + dimOffset)<=sum)
			currentGabor->feature|= LEFT_SIDE_POSITION_IN_RANGE;
		    else 
			currentGabor->feature&= ~LEFT_SIDE_POSITION_IN_RANGE;

		}

		findRSRCVariables(mp5Parameters,gaborDictionary,currentGabor,dataTable,
				  firstStart,firstStop,secondStart,secondStop,LIMIT_RANGE,&RS,&RC);

	    }
    	    else
    	    {
    		RS = 0.0;
		RC = 0.0;
	    }
	}
	else
	{
	    findRSRCVariables(mp5Parameters,gaborDictionary,currentGabor,dataTable,firstStart,firstStop,secondStart,secondStop,LIMIT_RANGE,&RS,&RC);
	}
    
	RS = (*(currentGabor->RS + channelNumber)) - (*(bestModulusesTable + channelNumber))*RS;
	RC = (*(currentGabor->RC + channelNumber)) - (*(bestModulusesTable + channelNumber))*RC;

    }
    *(currentGabor->RS + channelNumber) = (float)RS;
    *(currentGabor->RC + channelNumber) = (float)RC;
}

STATUS findUnknowPhaseDI(Gabor *gabor, double *modulus, unsigned short int channelNumber)
{

    const double KS = gabor->KS;
    const double KC = gabor->KC;
    const double KM = gabor->KM; 
    const double RS = *(gabor->RS + channelNumber);
    const double RC = *(gabor->RC + channelNumber);

    double W1, W2;

    double phase     = 0.0;
    double amplitude = 0.0;

    double sinPhase, cosPhase;

    if(!(gabor->feature & DIRACDELTA))
    {
	if(solveSystemOfEquestions(&KC,&KM,&KM,&KS,&W1,&W2,&RC,&RS) == ERROR)
	{
	    fprintf(stderr,"\n\nGabor with following parameters:");
	    fprintf(stderr,"\n POSITION:  %hu      \
		            \n RIFLING:   %u      \
		            \n SCALE:     %hu\n",gabor->position,gabor->rifling,gabor->scaleIndex);
	    fprintf(stderr," Is being marked as a incorrect gabor\n");
	    fprintf(stderr," Because one can not estimate phase of this gabor \n");
    	    fprintf(stderr," according to Dobieslaw Ircha algorithms \n");

	    *(gabor->phase + channelNumber) = 0.0;
	    *modulus = 0.0;
	
	    gabor->feature|= INCORRECTGABOR;
	
	    return ERROR;
	}
	else
	    phase = atan2(-W2,W1);
    }
    else
	phase = (RC>0 ? 0 : M_PI);

    sincos(phase,&sinPhase,&cosPhase);

    const double sinPart    = KS*(sinPhase*sinPhase);
    const double cosPart    = KC*(cosPhase*cosPhase);
    const double sinCosPart = 2.0*KM*sinPhase*cosPhase;

    amplitude = sqrt((sinPart + cosPart ) - sinCosPart);
    *modulus  = (RC*cosPhase - RS*sinPhase)/amplitude;

    *(gabor->phase + channelNumber) = (float)phase;

    return SUCCESS;

}

STATUS findUnknowPhaseAM(Gabor *gabor, double *modulusTable, unsigned short int numberOfAnalysedChannels)
{
    if(!(gabor->feature & DIRACDELTA))
    {   
	unsigned short int channel;

        double tmpSumX = 0.0;
        double tmpSumY = 0.0;
        double sinPhase, cosPhase;
        double amplitude;

        const double KS = gabor->KS;
        const double KC = gabor->KC;
        const double KM = gabor->KM;
		double modulus;
        double phase;

	for (channel=0;channel<numberOfAnalysedChannels;channel++) 
	{
	    modulus =  *(modulusTable + channel);
	    phase   =  *(gabor->phase + channel);
	    tmpSumX += modulus*sin(2.0 * phase);    
	    tmpSumY += modulus*cos(2.0 * phase);    
	} 

	const double tmpPhase = atan2(tmpSumX,tmpSumY)/2.0;

	for(channel=0;channel<numberOfAnalysedChannels;channel++) 
	{
	    if ((*(gabor->RS + channel))>0.0)
	    {
			sincos(tmpPhase,&sinPhase,&cosPhase);
	    }
	    else
	    {
			sincos(tmpPhase+M_PI,&sinPhase,&cosPhase);
	    }
	    
	    amplitude = sqrt(KS*(sinPhase*sinPhase) + KC*(cosPhase*cosPhase) - 2.0*KM*sinPhase*cosPhase);

	    *(modulusTable + channel) = ((*(gabor->RC + channel))* cosPhase - (*(gabor->RS + channel))*sinPhase)/amplitude;
	    
	     if ((*(gabor->RS + channel))>0.0)
	    {
			*(gabor->phase + channel) = (float)tmpPhase;
	    }
	    else
	    {
			*(gabor->phase + channel) = (float)tmpPhase+M_PI;
	    }
	    
	     
	}

    }
    return SUCCESS;
}
