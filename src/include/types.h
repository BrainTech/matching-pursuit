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

#ifndef __TYPES_H__

	#define __TYPES_H__

	#include<stdio.h>
	#include"def.h"

	struct Node
	{
		void *data;
		struct Node *nextNode;
	};

	typedef struct Node QueueNode;


	typedef struct
	{
		unsigned int size;
		QueueNode *firstNode;
		QueueNode *lastNode;
	}Queue;


	typedef struct
	{
		char text[LENGTH_OF_STRING];
	}String;


	typedef struct
	{
		int  number;
		char text[LENGTH_OF_STRING];
	}Line;

	typedef struct
	{
		char name[LENGTH_OF_NAME_OF_CONFIG_FILE];
		FILE *file;
		Queue *stringQueue;
		Queue *lineQueue;
	}ConfigFile;

	typedef struct
	{
		unsigned short int position;
		unsigned short int scaleIndex;
		unsigned       int rifling;

		float KS;
		float KC;
		float KM;
		float *RS;
		float *RC;
		float *phase;

		/* byte which describes features of the gabor:
		   h g f e d c b a

		   h - a  - bits

		   a       
		          0  - for some reasons incorrect gabor
		          1  - correct gabor

		   b      0 - undefined
		          1 - Dirack Delta

		   c      0 - undefined
			  1 - GaborWave

		   d      0 - undefined
			  1 - FFTWAVE

		   e
		          1 - if position of the gabor is <=dimOffset/2
			  0 - if position of the gabor is >dimOffset/2

                   f      1 - if position of the gabor is <=range/2. We don't have estimate dot gabor in full range of the signal.
			  0 - if position of the gabor is >range/2   Thanks to calculating gabor - gabor dot product, one can find 
		                                                     the range, where dot product is grater the some EPS
		   	    					     bits 'f' and 'g' inlcude information how the center of the 
							             gabor is situated with respect to center of this range
	           g      0  - gabor was not chosen
		          1  - gabor was chosen

		   h 
		          not defined

		   for example:

		   gaborFeature = 00010101 means:

		   uncorrect (will not be taken into further analysis), Gabor Wave, which is
		   situated in the left part of the offset, and 'right side' in the range, where dot product is grater then EPS
		*/

		unsigned char feature;
	}__attribute__((packed)) Gabor;
	
	typedef struct
	{
		unsigned int  sizeOfDictionary;
		unsigned char typeOfDictionary;

		double scaleToPeriodFactor;
		double dilationFactor; /* parameter of the Optimal Dictionary - "density factor" */
		unsigned short int periodDensity;

    		unsigned short int numberOfStepsInScale;

		double       basicStepInPositionInSignal;
		double       basicStepInFrequencyInSignal;
		
		double       basicStepInFrequencyInOptimalDictionary;
		unsigned int basicStepInPeriodInOptimalDictionary;
		double       basicStepInPositionInOptimalDictionary;
				
    		unsigned short int *tableOfScalesInOptimalDictionary;
    		unsigned       int *tableOfPeriodsInOptimalDictionary;
	        double             *tableOfFrequenciesInOptimalDictionary;
		double             *tableOfPositionsInOptimalDictionary;

		unsigned       int *numberOfStepsInFrequencyAtParticularScale;
		unsigned short int *numberOfStepsInPositionAtParticularScale;
		
    		unsigned       int numberSinCosFunctions;

		Gabor *gaborsTable;

		unsigned char memoryAllocated; 
						// 0x01 if tableOfScalesInOptimalDictionary is allocated
						// 0x02 if rest of pointers point to some alocated memory

	}GaborDictionary;

	typedef struct
	{

		char nameOfDataFile[LENGTH_OF_NAME_OF_DATA_FILE];
		char nameOfFileWhereDictionaryWillBeDroped[LENGTH_OF_NAME_OF_DATA_FILE];
		char nameOfFileWhereFitedGaborsWillBeDroped[LENGTH_OF_NAME_OF_DATA_FILE];
	        unsigned short int numberOfResultsFiles;
		char **namesOfResultFiles;
		char extensionOfResultFile[LENGTH_OF_EXTENSION_OF_RESULTS_FILE];
		char nameOfOutputDirectory[LENGTH_OF_OUTPUT_DIRECTORY];

		FILE *dataFile;
		FILE **resultFiles;
		FILE *dictionaryFile;
		FILE *fitedGaborsFile;

		unsigned short int sizeOfHeader;
		unsigned short int sizeOfTail;

		unsigned short int numberOfChannels;
		unsigned int       numberOfPoints;         /* number of samples in data file per channel  */
		unsigned short int numberOfPointsInOffset; /* number of samples in offset per channel     */
		unsigned short int numberOfOffsets;

		double samplingFrequency;

		unsigned short int numberOfChosenChannels;
		unsigned short int *chosenChannels;

		unsigned short int numberOfChosenOffsets;
		unsigned short int *chosenOffsets;

		unsigned int samplesBesideOffsets;

		unsigned char dataFormat;

		double convRate;

		unsigned short int dimOffset; /* number of samples in offset */
		unsigned       int dimExpand; /* 3 * dimOffset =  size of the signal with boundary conditions */

		double **rawDataMatrix;
		double **processedDataMatrix;

		unsigned char writingMode;
		unsigned char verbose;
		unsigned char allocatedElements;

        }DataParameters;

	typedef struct
	{
		unsigned short int numberOfAnalysedChannels;
		unsigned short int dimOffset;      /* number of samples in offset */
		unsigned       int dimExpTable;    /* 2 * dimOffset =  size of the vectors the of the exp values */
		unsigned       int dimExpand;      /* 3 * dimOffset =  size of the signal with boundary conditions */

		double **sinTable;
		double **cosTable;
		double **expTable;

		double *sinGaborTable;
		double *cosGaborTable;

		Gabor  *previousGabor;
		double **prevGaborTable;

		double *singleChannelSignalTable;
		double **multiChannelSignalTable;
		double *singleChannelResidueTable;
		double **multiChannelResidueTable;
		double *meanSignalTable;
		double *meanResidueTable;

		double totalSignalEnergy;           /* inlcudes sum of energies in each channel */
		double totalResidueEnergy;          /* includes sum of residue energy in each channels */
		double meanSignalEnergy;            /* mean signal energy over the channels */
		double meanResidueEnergy;           /* mean residue energy over the channels */

		double *signalEnergyInEachChannel;  /* energy of signal in each channel */
		double *residueEnergyInEachChannel; /* energy of residue in each channel */

		unsigned short int maxNumberOfIterations;
		double energyPercent;

		double *bestModulusesTable;

		unsigned char reinitDictionary;
		unsigned char MPType;

		Queue  *fitted;

		double  LOG_EPS_DOT_PRODUCT;
		BOOLEAN memoryAllocated; /* = FALSE if memory for mp5Parameters (tables, dictionary) was not allocated, otherwise TRUE */

	}MP5Parameters;

#endif
