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


#include<stdlib.h>
#include<strings.h>
#include"types.h"

#ifdef __MINGW32__
    #define bzero(ptr,size) memset (ptr, 0, size);
#endif

/* double matrix operation */

double **dMatrixAllocate(unsigned int numberOfRows, unsigned int numberOfColumns)
{
	unsigned int row;
	double **matrix = NULL;
	double  *vector = NULL;

	vector = (double *)malloc(numberOfRows*numberOfColumns*sizeof(double));

	matrix = (double **)malloc(numberOfRows*sizeof(double *));

	for(row=0;row<numberOfRows;row++)
		*(matrix + row) = vector + row*numberOfColumns;

	return matrix;
}

void dMatrixFree(double **matrix)
{
	free(matrix[0]);
	free(matrix);
}

void dSetMatrixZero(double **matrix, unsigned int numberOfRows, unsigned int numberOfColumns)
{
	unsigned int row;

	for(row=0;row<numberOfRows;row++)
		bzero((void *)(*(matrix + row)),numberOfColumns*sizeof(double));
}

double **dVariableMatrixAllocate(unsigned int numberOfRows, unsigned int *tableOfNumbersOfColumns)
{
    unsigned int row;

    double **matrix = (double **)malloc(numberOfRows*sizeof(double *));

    for(row=0;row<numberOfRows;row++)
		*(matrix + row) = (double *)malloc((*(tableOfNumbersOfColumns + row))*sizeof(double));

	return matrix;
}

void dVariableMatrixFree(double **matrix, unsigned int numberOfRows)
{
	unsigned int row;

	for(row=0;row<numberOfRows;row++)
		free(*(matrix + row));

	free(matrix);
}

void dVariableSetMatrixZero(double **matrix, unsigned int numberOfRows, unsigned int *tableOfNumbersOfColumns)
{
	unsigned int row;

	for(row=0;row<numberOfRows;row++)
		bzero((void *)(*(matrix + row)),(*(tableOfNumbersOfColumns + row))*sizeof(double));
}
/* float matrix operation */

float **fMatrixAllocate(unsigned int numberOfRows, unsigned int numberOfColumns)
{
	unsigned int row;
	float **matrix = NULL;
	float  *vector = NULL;

	vector = (float *)malloc(numberOfRows*numberOfColumns*sizeof(float));

	matrix = (float **)malloc(numberOfRows*sizeof(float *));

	for(row=0;row<numberOfRows;row++)
		*(matrix + row) = vector + row*numberOfColumns;

	return matrix;
}

void fMatrixFree(float **matrix)
{
	free(matrix[0]);
	free(matrix);
}

void fSetMatrixZero(float **matrix, unsigned int numberOfRows, unsigned int numberOfColumns)
{
	unsigned int row;

	for(row=0;row<numberOfRows;row++)
		bzero((void *)(*(matrix + row)),numberOfColumns*sizeof(float));
}

unsigned int **iuMatrixAllocate(unsigned int numberOfRows, unsigned int numberOfColumns)
{
	unsigned int row;
	unsigned int **matrix = NULL;
	unsigned int  *vector = NULL;

	vector = (unsigned int *)malloc(numberOfRows*numberOfColumns*sizeof(double));

	matrix = (unsigned int **)malloc(numberOfRows*sizeof(unsigned int *));

	for(row=0;row<numberOfRows;row++)
		*(matrix + row) = vector + row*numberOfColumns;

	return matrix;
}

void uiMatrixFree(unsigned int **matrix)
{
	free(matrix[0]);
	free(matrix);
}

void uiSetMatrixZero(unsigned int **matrix, unsigned int numberOfRows, unsigned int numberOfColumns)
{
	unsigned int row;

	for(row=0;row<numberOfRows;row++)
		bzero((void *)(*(matrix + row)),numberOfColumns*sizeof(unsigned int));
}

void countMeanSignalOverChannelsInSingleEpoch(MP5Parameters *mp5Parameters)
{
	const unsigned int epochExpandedSize = mp5Parameters->epochExpandedSize;
	unsigned int sample;
	double tmpDataValue;

	double **multiChannelSignalTable = mp5Parameters->multiChannelSignalTable;
	double *meanSignalTable          = *mp5Parameters->meanSignalTable;
	unsigned int channel;

	for(sample=0;sample<epochExpandedSize;sample++)
	{
		tmpDataValue = 0.0;

		for(channel=0;channel<mp5Parameters->numberOfAnalysedChannels;channel++)
			tmpDataValue+= (*(*(multiChannelSignalTable + channel) + sample));

		*(meanSignalTable + sample) = tmpDataValue/mp5Parameters->numberOfAnalysedChannels;
	}
}

void countMeanSignalOrResidumOverChannelsAndTrials(MP5Parameters *mp5Parameters)
{
	const unsigned int epochExpandedSize = mp5Parameters->epochExpandedSize;

	unsigned int sample;
	double tmpDataValue;
	double **multiChannelSignalTable = mp5Parameters->multiChannelSignalTable;
	double **meanSignalTable         = mp5Parameters->meanSignalTable;
	unsigned int channel;

	for(sample=0;sample<epochExpandedSize;sample++)
	{
		tmpDataValue = 0.0;

		for(channel=0;channel<mp5Parameters->numberOfAnalysedChannels;channel++)
			tmpDataValue+= (*(*(multiChannelSignalTable + channel) + sample));

		*(*(meanSignalTable) +  sample) = tmpDataValue/mp5Parameters->numberOfAnalysedChannels;
	}
}

void countMeanSignalOrResidumOverChannels(MP5Parameters *mp5Parameters, double **multiChannelSignalTable)
{
	const unsigned int epochExpandedSize = mp5Parameters->epochExpandedSize;

	unsigned int sample;
	double tmpDataValue;
	double **meanSignalTable = mp5Parameters->meanSignalTable;
	unsigned int channelEpochNumber;
	unsigned short int channel;
	unsigned short int epoch;

	for(sample=0;sample<epochExpandedSize;sample++)
	{

		for(epoch=0;epoch<mp5Parameters->numberOfSelectedEpochs;epoch++)
		{
			tmpDataValue = 0.0;

			for(channel=0;channel<mp5Parameters->numberOfSelectedChannels;channel++)
			{
				channelEpochNumber = epoch*mp5Parameters->numberOfSelectedChannels + channel;
				tmpDataValue+= (*(*(multiChannelSignalTable + channelEpochNumber) + sample));

			}

			*(*(meanSignalTable + epoch) + sample)= tmpDataValue/mp5Parameters->numberOfSelectedChannels;

		}
	}

}

void countMeanSignalOrResidumOverEpochs(MP5Parameters *mp5Parameters, double **multiChannelSignalTable)
{
	const unsigned int epochExpandedSize = mp5Parameters->epochExpandedSize;

	unsigned int sample;
	double tmpDataValue;
	double **meanSignalTable = mp5Parameters->meanSignalTable;
	unsigned int channelEpochNumberNumber;
	unsigned short int channel;
	unsigned short int epoch;

	for(sample=0;sample<epochExpandedSize;sample++)
	{
		for(channel=0;channel<mp5Parameters->numberOfSelectedChannels;channel++)
		{
			tmpDataValue = 0.0;

			for(epoch=0;epoch<mp5Parameters->numberOfSelectedEpochs;epoch++)
			{
				channelEpochNumberNumber = epoch*mp5Parameters->numberOfSelectedChannels + channel;
				tmpDataValue+= (*(*(multiChannelSignalTable + channelEpochNumberNumber) + sample));
			}

			*(*(meanSignalTable + channel) + sample)= tmpDataValue/mp5Parameters->numberOfSelectedEpochs;
		}
	}
}
