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


#ifndef __TYPES_H__

	#define __TYPES_H__

	#include<stdio.h>
	#include"def.h"
	#include"fftw3.h"

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
		char text[LENGTH_OF_LINE];
	}Line;

	typedef struct
	{
		char name[LENGTH_OF_NAME_OF_CONFIG_FILE];
		FILE  *file;
		Queue *stringQueue;
		Queue *lineQueue;
	}ConfigFile;

	typedef struct
	{
		unsigned       int position;
		unsigned short int scaleIndex;
		unsigned       int rifling;
		double			   randomShiftInFrequency;

		
		float KS;
		float KC;
		float KM;
		float *RS;
		float *RC;
		float *phase;

		/* byte which describes features of the atom:
			h g f e d c b a

			h - a  - bits

			a      
				0 - undefined
				1 - Dirack Delta

			b
				0 - undefined
				1 - Gauss Functions

			c    
				0 - undefined
				1 - SINCOSWAVE

			d    
				0 - undefined
				1 - GaborWave
			e      
				0  - for some reasons incorrect atom
				1  - correct atom
				
			f
				1 - if position of the atom is <=offsetDimension/2
				0 - if position of the atom is >offsetDimension/2

			g     
				1 - if position of the atom is <=range/2. We don't have estimate dot atom in full range of the signal.
				0 - if position of the atom is >range/2   Thanks to calculating atom - atom dot product, one can find  
				      the range, where dot product is grater the some EPS bits 'f' and 'g' inlcude information how the center of the 	
				      atom is situated with respect to center of this range
	           
			 h      
				0  - atom was not chosen
				1  - atom was chosen

		for example:
			feature = 0000101 means:
				uncorrect (will not be taken into further analysis),  Atom Wave, which i
				situated in the left part of the offset, and 'right side' in the range, where dot product is grater then EPS
		*/

		unsigned char feature;
	}__attribute__((packed)) Atom;
	
	typedef struct
	{
		unsigned int  sizeOfDictionary;
		unsigned int  numberOfCorrectGabors;
		unsigned char typeOfDictionary;

		double scaleToPeriodFactor;
		double dilationFactor; /* parameter of the Optimal Dictionary - "density factor" */
		long int randomSeed; 
		unsigned short int periodDensity;

    	unsigned short int numberOfStepsInScale;

		double       basicStepInPositionInSignal;
		double       basicStepInFrequencyInSignal;
		
		double       basicStepInFrequencyInOptimalDictionary;
		unsigned int basicStepInPeriodInOptimalDictionary;
		double       basicStepInPositionInOptimalDictionary;
				
    	double			   *tableOfScalesInOptimalDictionary;
    	unsigned       int *tableOfPeriodsInOptimalDictionary;
		double             *tableOfFrequenciesInOptimalDictionary;
		double             *tableOfPositionsInOptimalDictionary;

		unsigned       int *numberOfStepsInFrequencyAtParticularScale;
		unsigned       int *numberOfStepsInPositionAtParticularScale;

		unsigned       int numberOfDiracFunctions;	
		unsigned       int numberOfGaussFunctions;	
		unsigned       int numberOfSinCosFunctions;
		unsigned       int numberOfNonFFTAtoms; 

		unsigned char diracInDictionary;
		unsigned char gaussInDictionary;
		unsigned char sinCosInDictionary;

		Atom *atomsTable;

	}Dictionary;

	typedef struct
	{
		char nameOfDataFile[LENGTH_OF_NAME_OF_DATA_FILE];
		char nameOfResultsFile[LENGTH_OF_NAME_OF_RESULTS_FILE];
		char extensionOfResultsFile[LENGTH_OF_EXTENSION_OF_RESULTS_FILE];
		char nameOfOutputDirectory[LENGTH_OF_OUTPUT_DIRECTORY];

		FILE *dataFile;
		FILE *resultsFile;

		unsigned short int sizeOfHeader;
		unsigned short int sizeOfTail;

		unsigned short int numberOfChannelsInDataFile;
		unsigned int       numberOfPoints;         /* number of samples in data file per channel  */
		unsigned short int numberOfOffsets;

		double             samplingFrequency;
		unsigned short int numberOfChosenChannels;
		unsigned short int *chosenChannels;
		unsigned short int numberOfChosenOffsets;
		unsigned short int *chosenOffsets;
		unsigned int       samplesBesideOffsets;
		unsigned char      dataFormat;

		double pointsPerMicrovolt;
		double **rawDataMatrix;
		double **processedDataMatrix;

		unsigned char writingMode;
		unsigned char numberOfThreads;

		unsigned short int numberOfAllocatedChannels;
		unsigned short int numberOfAnalysedChannels;
		unsigned       int offsetDimension;         /* number of samples in offset */
		unsigned       int marginalDimension;       /* marginal condition dimension - how many zeros add to the signal from one side */
		unsigned       int exponensTableDimension;  /* size of the vectors the of the exp values */
		unsigned       int offsetExpandedDimension; /* size of the signal with boundary conditions */
		unsigned       int fftTableDimension;       /* size of fftTable for FFT algorithm */
		
		/* generally:
			marginalDimension offsetDimension marginalDimension
					
					offsetExpandedDimension 
		*/

		double **sinTable;
		double **cosTable;
		double **expTable;
		
		unsigned char FFT;
		fftw_plan     fftwPlan;
		double       *fftTableInA;
		fftw_complex *fftTableInB;
		fftw_complex *fftTableOut;
		
		unsigned char analiticalDotProduct;
		
		double *sinAtomTable;
		double *cosAtomTable;

		Atom   *previousAtom;
		double **prevAtomTable;

		double *zeroSignalTable;
		double *gaussSignalTable;
		double *singleChannelSignalTable;
		double **multiChannelSignalTable;
		double *singleChannelResidueTable;
		double **multiChannelResidueTable;
		double *meanSignalTable;
		double *meanResidueTable;

		double totalSignalEnergy;           /* inlcudes sum of energies over all channels */
		double totalResidueEnergy;          /* includes sum of residue over all channels  */
		double oneChannelSignalEnergy;      /* signal energy in one channels (for SMP or MMP2)    */
		double oneChannelResidueEnergy;     /* residue energy in one channels (for SMP or MMP2) */

		double *signalEnergyInEachChannel;  /* energy of signal in each channel */
		double *residueEnergyInEachChannel; /* energy of residue in each channel */

		unsigned short int maximalNumberOfIterations;
		double             energyPercent;

		double *bestModulusesTable;
		float  *bestPhasesTable;

		unsigned char reinitDictionary;
		unsigned char MPType;

		Queue  *fitted;
		unsigned int  maxGaborScale;
		unsigned char bookWithSignal;
		unsigned char progressBar;
		
		unsigned char accuracy;
		
	}MP5Parameters;

#endif
