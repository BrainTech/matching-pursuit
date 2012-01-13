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


#ifndef __TYPES_H__

	#define __TYPES_H__

	#include<stdio.h>
	#include"def.h"
	#include"fftw3.h"

    #ifdef __MINGW32__
        #define bzero(ptr,size) memset (ptr, 0, size);
        #define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }

        // v3.11 of mingw-runtime deprecates _sleep in favor of Win32 API Sleep
        void __stdcall Sleep(unsigned);
        #define sleep(t) Sleep(1000*t);
    #endif

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

		float KS;
		float KC;
		float KM;
		float *RS;
		float *RC;
		float *phase;

		/* byte which describes features of the atom:
			p o n m l k j i h g f e d c b a

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
				1 - if position of the atom is <=epochSize/2
				0 - if position of the atom is >epochSize/2

			g
				1 - if position of the atom is <=range/2. We don't have estimate dot atom in full range of the signal.
				0 - if position of the atom is >range/2   Thanks to calculating atom - atom dot product, one can find
				      the range, where dot product is grater the some EPS bits 'f' and 'g' inlcude information how the center of the
				      atom is situated with respect to center of this range

			 h
				0  - atom was not seleced
				1  - atom was selected during mp decomposition

			 i
				0  - atom was not selected to stochastic dictionary
				1  - atom was selected to stochastic dictionary

		for example:
			feature = 10000101 means:
				incorrect (will not be taken into further analysis), Atom Wave, which is
				situated in the left part of the epoch, and 'right side' in the range, where dot product is grater then EPS
                atom was taken to the stochastic dictionary
		*/

		unsigned short int feature;
	}__attribute__((packed)) Atom;

	typedef struct
	{
		unsigned char typeOfDictionary;
		double scaleToPeriodFactor; 

		double dilationFactor;                           /* parameter of the Optimal Dictionary - "density factor" */
		double stochasticDictionaryReductionCoefficient; /* parameter of the Optimal Dictionary - "density factor" */
        long int randomSeed;

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

		unsigned       int numberOfInitialDiracFunctions;
		unsigned       int numberOfInitialGaussFunctions;
		unsigned       int numberOfInitialSinCosFunctions;
		unsigned       int numberOfInitialGaborFunctions;
		unsigned       int numberOfInitialCorrectGabors;
		unsigned       int initialNumberOfAtoms;
        
		unsigned       int numberOfFinalDiracFunctions;
		unsigned       int numberOfFinalGaussFunctions;
		unsigned       int numberOfFinalSinCosFunctions;
		unsigned       int numberOfFinalGaborFunctions;
		unsigned       int numberOfFinalCorrectGabors;
        unsigned       int finalNumberOfAtoms;

		unsigned char diracInDictionary;
		unsigned char gaussInDictionary;
		unsigned char sinCosInDictionary;
		unsigned char gaborInDictionary;

		Atom **atomsTable;
		Atom *diracAtomsTable;
		Atom *gaussAtomsTable;
		Atom *sinCosAtomsTable;
		Atom *gaborAtomsTable;

	}Dictionary;
    
	typedef struct
	{
		char nameOfDataFile[LENGTH_OF_NAME_OF_DATA_FILE];
		char nameOfResultsFile[LENGTH_OF_NAME_OF_RESULTS_FILE];
		char nameOfOutputDirectory[LENGTH_OF_OUTPUT_DIRECTORY];

		FILE *dataFile;
		FILE *resultsFile;

		unsigned short int sizeOfHeader;
		unsigned short int sizeOfTail;

		unsigned short int numberOfChannelsInDataFile;
		unsigned int       numberOfPoints;         /* number of samples in data file per channel  */
		unsigned short int numberOfEpochs;

		double              samplingFrequency;
		unsigned short int  numberOfSelectedChannels;
		unsigned short int *selectedChannels;
		unsigned short int  numberOfSelectedEpochs;
		unsigned short int *selectedEpochs;
		unsigned int        samplesBesideEpochs;

		double pointsPerMicrovolt;
		double **rawDataMatrix;
		double **processedDataMatrix;

		unsigned char writingMode;

		unsigned short int numberOfAllocatedChannels;
		unsigned       int numberOfAnalysedChannels;
		unsigned 	   int numberOfReadChannelsAndEpochs;
		unsigned       int epochSize;               /* number of samples in epoch */
		unsigned       int marginalSize;            /* marginal condition size - how many zeros add to the signal from one side */
		unsigned       int exponensTableSize;       /* size of the vectors the of the exp values */
		unsigned       int epochExpandedSize;       /* size of the signal with boundary conditions */
		unsigned       int fftTableSize;            /* size of fftTable for FFT algorithm */

		/* generally:
			marginalSize epochSize marginalSize

					  epochExpandedSize
		*/

		double **sinTable;
		double **cosTable;
		double **expTable;

		unsigned char FFT;
		fftw_plan     fftwPlan;
		double       *fftTableIn;
		fftw_complex *fftTableOut;

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
		double **meanSignalTable;
		double **meanResidueTable;

		double totalSignalEnergy;           /* inlcudes sum of energies over all channels */
		double totalResidueEnergy;          /* includes sum of residue over all channels  */
		double oneChannelSignalEnergy;      /* signal energy in one channels (for SMP or MMP2)  */
		double oneChannelResidueEnergy;     /* residue energy in one channels (for SMP or MMP2) */

		double *signalEnergyInEachChannel;  /* energy of signal in each channel */
		double *residueEnergyInEachChannel; /* energy of residue in each channel */
		double *meanSignalEnergyInEachChannel;  /* energy of signal in each channel for MMP12 MMP21 MMP23 MMP32 */
		double *meanResidueEnergyInEachChannel; /* energy of residue in each channel for MMP12 MMP21 MMP23 MMP32 */

		unsigned short int maximalNumberOfIterations;
		double             energyPercent;

		double *bestModulusesTable;
		float  *bestPhasesTable;

		unsigned char  reinitDictionary;
		unsigned short int MPType;

		Queue  *fitted;

		unsigned char bookWithSignal;
		unsigned char progressBar;

	}MP5Parameters;

	typedef struct
	{
		unsigned int step;
		unsigned int stepInToolbar;
		unsigned char applicationMode;
				
	}Progress;

#endif
