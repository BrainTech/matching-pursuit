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

#define _GNU_SOURCE

#include<dirent.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<strings.h>
#include<unistd.h>
#include"include/def.h"
#include"include/gabor.h"
#include"include/io.h"
#include"include/matrix.h"
#include"include/mp5.h"
#include"include/new_io.h"
#include"include/queue.h"
#include"include/types.h"
#include"include/vector.h"

#ifdef __MINGW32__
#define bzero(ptr,size) memset (ptr, 0, size);
#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

static void asciiFileSeek(FILE *asciiFile, unsigned long int lineNumber)
{
	unsigned long int line;
	fseek(asciiFile,0L,SEEK_SET);

	for(line=0;line<lineNumber;line++)
	    fscanf(asciiFile,"%*[^\n]\n");
}

static void returnAmplitudeAndModulusForMMP2DI(MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary, Gabor *gabor, float *amplitude, float *modulus, unsigned short int channelNumber)
{
	unsigned int sample;
	double       tmpModulus = 0.0;
	unsigned int dimExpand = mp5Parameters->dimExpand;

        const double phase = (*(gabor->phase + channelNumber));
        const double KS = gabor->KS;
	const double KC = gabor->KC;
	const double KM = gabor->KM;

	double sinPhase, cosPhase;

	sincos(phase,&sinPhase,&cosPhase);

	const double sinPart    = KS*(sinPhase*sinPhase);
	const double cosPart    = KC*(cosPhase*cosPhase);
	const double sinCosPart = 2.0*KM*sinPhase*cosPhase;
	
	double *signalInParticularChannel = *(mp5Parameters->multiChannelSignalTable + channelNumber);
	double *prevGaborTable            = *(mp5Parameters->prevGaborTable);

	makeSinCosGaborTable(mp5Parameters,gaborDictionary,gabor);
	makeGaborTable(mp5Parameters,gabor,0);	
	
	for(sample=0;sample<dimExpand;sample++)
	    tmpModulus+= (*(signalInParticularChannel + sample))*(*(prevGaborTable + sample));

	findResidue(signalInParticularChannel,prevGaborTable,tmpModulus,dimExpand);	

	*amplitude = (float)sqrt((sinPart + cosPart ) - sinCosPart);

	*modulus = (float)(tmpModulus);
}

static void returnAmplitudeAndModulusDI(Gabor *gabor, float *amplitude, float *modulus, unsigned short int channelNumber)
{
    const double RS    = (*(gabor->RS    + channelNumber));
    const double RC    = (*(gabor->RC    + channelNumber));
    const double phase = (*(gabor->phase + channelNumber));

    const double KS = gabor->KS;
    const double KC = gabor->KC;
    const double KM = gabor->KM;

    double sinPhase, cosPhase;

    sincos(phase,&sinPhase,&cosPhase);

    const double sinPart    = KS*(sinPhase*sinPhase);
    const double cosPart    = KC*(cosPhase*cosPhase);
    const double sinCosPart = 2.0*KM*sinPhase*cosPhase;

    *amplitude = (float)sqrt((sinPart + cosPart ) - sinCosPart);
    *modulus   = (float)((RC*cosPhase - RS*sinPhase)/(*amplitude));
}

static void printChosenChannels(DataParameters *dataParameters)
{
    int i;
    const int howMany = 10;
    BOOLEAN breakLine = FALSE;

    printf("                               ");

    for(i=0;i<dataParameters->numberOfChosenChannels;i++)
    {
	breakLine = FALSE;

	if((i%howMany)==0 && i>0)
	{
	    printf("\n");
	    printf("                             ");
	    breakLine = TRUE;
	}

	printf("%-7d ",dataParameters->chosenChannels[i]);
    }

    if(!breakLine)
	printf("\n");
}

static void printChosenOffsets(DataParameters *dataParameters)
{
    int i;
    const int howMany = 10;
    BOOLEAN breakLine = FALSE;

    printf("                               ");

    for(i=0;i<dataParameters->numberOfChosenOffsets;i++)
    {
	breakLine = FALSE;

	if((i%howMany)==0 && i>0)
	{
	    printf("\n");
	    printf("                               ");
	    breakLine = TRUE;
	}

	printf("%-7d ",dataParameters->chosenOffsets[i]);
    }

    if(!breakLine)
	printf("\n");
}

void setDataParameters(DataParameters *dataParameters)
{
    unsigned short int strCounter;
    unsigned short int channelNumber;
    unsigned short int dimOffset            = dataParameters->dimOffset;
    unsigned	   int dimExpand            = dataParameters->dimExpand;
    unsigned short int numberOfChannels     = dataParameters->numberOfChannels;
    unsigned short int numberOfResultsFiles = dataParameters->numberOfResultsFiles;
  
    dataParameters->namesOfResultFiles = (char **)malloc(numberOfResultsFiles*sizeof(char *));
    
    for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
	*(dataParameters->namesOfResultFiles + channelNumber) = charVectorAllocate(LENGTH_OF_NAME_OF_RESULTS_FILE*sizeof(char));

    for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
	for(strCounter=0;strCounter<LENGTH_OF_NAME_OF_RESULTS_FILE;strCounter++)
	    *(*(dataParameters->namesOfResultFiles + channelNumber) + strCounter) = '\0';

    dataParameters->allocatedElements|= NAMES_OF_RESULT_FILES_ALLOCATED;

    dataParameters->resultFiles = (FILE **)malloc(numberOfResultsFiles*sizeof(FILE *));

    for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
	*(dataParameters->resultFiles + channelNumber) = NULL;

    dataParameters->allocatedElements|= RESULT_FILES_ALLOCATED;

    dataParameters->rawDataMatrix = dMatrixAllocate(numberOfChannels,dimOffset);

    dataParameters->allocatedElements|= RAW_DATA_MATRIX_ALLOCATED;

    dataParameters->processedDataMatrix = dMatrixAllocate(numberOfChannels,dimExpand);

    dataParameters->allocatedElements|= PROCESSED_DATA_MATRIX_ALLOCATED;
}

void freeDataParameters(DataParameters *dataParameters)
{
    unsigned short int channelNumber;

    if(dataParameters->allocatedElements & NAMES_OF_RESULT_FILES_ALLOCATED)
    {
	for(channelNumber=0;channelNumber<dataParameters->numberOfResultsFiles;channelNumber++)
	    charVectorFree(*(dataParameters->namesOfResultFiles + channelNumber));

	free(dataParameters->namesOfResultFiles);
    }

    if(dataParameters->allocatedElements & RESULT_FILES_ALLOCATED)
	free(dataParameters->resultFiles);

    if(dataParameters->allocatedElements & CHOSEN_CHANNELS_ALLOCATED)
	usiVectorFree(dataParameters->chosenChannels);

    if(dataParameters->allocatedElements & CHOSEN_OFFSETS_ALLOCATED)
	usiVectorFree(dataParameters->chosenOffsets);

    if(dataParameters->allocatedElements & RAW_DATA_MATRIX_ALLOCATED)
	dMatrixFree(dataParameters->rawDataMatrix,dataParameters->numberOfChannels);

    if(dataParameters->allocatedElements & PROCESSED_DATA_MATRIX_ALLOCATED)
	dMatrixFree(dataParameters->processedDataMatrix,dataParameters->numberOfChannels);
}

STATUS testFilesAndDirectories(DataParameters *dataParameters, const ConfigFile *configFile, char *info)
{
    unsigned short int channelNumber;
    unsigned short int numberOfResultsFiles = dataParameters->numberOfResultsFiles;
    char **namesOfResultFiles = dataParameters->namesOfResultFiles;

    DIR  *directory;
    FILE *file;    
    
    if((directory = opendir(dataParameters->nameOfOutputDirectory))==NULL)
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n CAN NOT OPEN DIRECTORY: %s TO WRITE RESULTS\n CHECK FILE: %s\n",dataParameters->nameOfOutputDirectory,configFile->name);
	return ERROR;
    }
    else
	closedir(directory);
    
    for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
    {
	if(dataParameters->writingMode & CREATE_FILE)
	{
	    file = fopen(*(namesOfResultFiles + channelNumber),"rb");;

	    if(file!=NULL)
	    {
		sprintf(info, "\n DATA PARAMETERS ERROR: \n YOU ARE GOING TO OVERWRITE FILE: %s\n OK, BUT THINK TWICE AND IF NECESSERY REMOVE THIS FILE IT BY HAND\n THE PROGRAM DOES NOT OVERWRITE EXISTING FILES FOR YOU\n",*(dataParameters->namesOfResultFiles + channelNumber));
		fclose(file);
		return ERROR;
	    }
	}
    }

    return SUCCESS;        

}

STATUS testDataParameters(DataParameters *dataParameters, GaborDictionary *gaborDictionary, MP5Parameters *mp5Parameters, char *info)
{
    if(dataParameters->numberOfChosenChannels > dataParameters->numberOfChannels)
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n NUMBER OF CHOSEN CHANELS > NUMBER OF CHANNELS IN FILE: %s\n",dataParameters->nameOfDataFile);
	return ERROR;
    }

    if(gaborDictionary->dilationFactor<=1)
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n DILATIONFACTOR MUST BE > 1 IN FILE: %s\n",dataParameters->nameOfDataFile);
	return ERROR;    
    }

    if((gaborDictionary->typeOfDictionary & OCTAVE_STOCH) && gaborDictionary->periodDensity == 1)
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n IF YOU SET TYPEOFDICTIONARY TO \"OCTAVE_STOCH\", PERIODDENSITY MUST BE > 1 IN FILE: %s\n",dataParameters->nameOfDataFile);
	return ERROR;    
    }

    if(mp5Parameters->energyPercent>=100.0)
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n ENERGYPERCENT SHOULD BE <= 100.0 %% IN FILE: %s\n",dataParameters->nameOfDataFile);
	return ERROR;
    }

    if(!(mp5Parameters->MPType & SMP) && (mp5Parameters->reinitDictionary==REINIT_IN_CHANNEL_DOMAIN || mp5Parameters->reinitDictionary==REINIT_AT_ALL))  
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n THERE IS NO SENSE TO USE MULTICHANNEL MP ALGORITHM AND REINIT DICTIONARY IN CHANNEL DOMAIN: %s\n",dataParameters->nameOfDataFile);
	return ERROR;
    }

    if(!(mp5Parameters->MPType & SMP) && (dataParameters->numberOfChosenChannels==1))  
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n YOU CHOSE ONLY ONE CHANNEL TO MULTICHANNEL ANALYSIS: %s\n",dataParameters->nameOfDataFile);
	return ERROR;
    }

    return SUCCESS;
}

void printInfoAboutData(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary)
{
    unsigned short int channel;
    printf(" \n");
    printf(" THE FOLLOWING PARAMETERS HAS BEEN READ: \n\n");
    printf(" NAME OF DATA FILE:            %s\n",dataParameters->nameOfDataFile);
    printf(" SIZE OF HEADER:               %-5hu\n",dataParameters->sizeOfHeader);
    printf(" TAIL OF HEADER:               %-5hu\n",dataParameters->sizeOfTail);
    printf(" SAMPLE FREQUENCY:             %-4.2f\n",dataParameters->samplingFrequency);
    printf(" FORMAT OF DATA:            ");
    if(dataParameters->dataFormat & FORMAT_ASCII)
	printf("   ASCII\n");
    else if(dataParameters->dataFormat & FORMAT_SHORT)
	printf("   SHORT\n");
    else if(dataParameters->dataFormat & FORMAT_FLOAT)
	printf("   FLOAT\n");
    printf(" NUMBER OF CHANNELS IN FILE:   %-5hu\n",dataParameters->numberOfChannels);
    printf(" NUMBER OF CHOSEN CHANNELS:    %-5hu\n",dataParameters->numberOfChosenChannels);

    printf(" CHOSEN CHANNELS:             \n");
    printChosenChannels(dataParameters);
    printf(" NUMBER OF POINTS IN OFFSET:   %-5hu\n",dataParameters->numberOfPointsInOffset);
    printf(" NUMBER OF CHOSEN OFFSETS:     %-5hu\n",dataParameters->numberOfChosenOffsets);
    printf(" CHOSEN OFFSETS:              \n");
    printChosenOffsets(dataParameters);

    printf(" TYPE OF DICTIONARY:           ");
    if(gaborDictionary->typeOfDictionary & OCTAVE_FIXED)
	printf("OCTAVE_FIXED\n");
    else if(gaborDictionary->typeOfDictionary & OCTAVE_STOCH)
	printf("OCTAVE_STOCH\n");
    printf(" DILATION FACTOR:              %-f\n",gaborDictionary->dilationFactor);
    printf(" PERIOD DENSITY:               %-hu\n",gaborDictionary->periodDensity);
    printf(" REINIT DICTIONATY:            ");
    if(mp5Parameters->reinitDictionary & NO_REINIT_AT_ALL)
	printf("NO REINIT AT ALL \n");
    else if(mp5Parameters->reinitDictionary & REINIT_IN_CHANNEL_DOMAIN)
	printf("REINIT IN CHANNEL DOMAIN \n");
    else if(mp5Parameters->reinitDictionary & REINIT_IN_OFFSET_DOMAIN)
	printf("REINIT IN OFFSET DOMAIN \n");
    else if(mp5Parameters->reinitDictionary & REINIT_AT_ALL)
	printf("REINIT AT ALL\n");

    printf(" SCALE TO PERIOD FACTOR:       %-f\n",gaborDictionary->scaleToPeriodFactor);
    printf(" MAX. GABORS NUMBER:           %-5hu\n",mp5Parameters->maxNumberOfIterations);
    printf(" ENERGY PERCENT:               %-4.2f\n",mp5Parameters->energyPercent);
    printf(" TYPE OF ALGORITHM:            ");
    if(mp5Parameters->MPType & SMP)
	printf("SINLGE CHANNEL MATCHING PURSUIT \n");
    else if(mp5Parameters->MPType & MMP1)
	printf("MULTICHANNEL MATCHING PURSUIT I\n");
    else if(mp5Parameters->MPType & MMP2)
	printf("MULTICHANNEL MATCHING PURSUIT II\n");
    else if(mp5Parameters->MPType & MMP3)
	printf("MULTICHANNEL MATCHING PURSUIT III\n");
    printf(" RESULTS WILL BE WRITTEN TO THE FOLLOWING FILES: \n");
    for(channel=0;channel<dataParameters->numberOfResultsFiles;channel++)
	printf("                                                 %s\n",*(dataParameters->namesOfResultFiles + channel));

    printf(" \n");
    fflush(stdout);
}

STATUS openBinaryDataFile(DataParameters *dataParameters, char *info)
{
    if(!(dataParameters->dataFile = fopen(dataParameters->nameOfDataFile,"rb")))
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n CAN NOT OPEN DATA FILE:  %s\n",dataParameters->nameOfDataFile);
	return ERROR;
    }

    return SUCCESS;
}

STATUS openAsciiDataFile(DataParameters *dataParameters, char *info)
{
    if(!(dataParameters->dataFile = fopen(dataParameters->nameOfDataFile,"rt")))
    {
	sprintf(info,"\n DATA PARAMETERS ERROR: \n CAN NOT OPEN DATA FILE:  %s\n",dataParameters->nameOfDataFile);
	return ERROR;
    }

    return SUCCESS;
}

void createNamesOfResultFiles(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary)
{
	char **namesOfResultFiles = dataParameters->namesOfResultFiles;
	unsigned short int channelNumber;
	const unsigned short int numberOfResultsFiles = dataParameters->numberOfResultsFiles;
	const unsigned short int *chosenChannels = dataParameters->chosenChannels;
	unsigned short int lengthOfDataFileWithOutExpand;

	char *dot = strrchr(dataParameters->nameOfDataFile,'.');
	int len;

	char tmpString[LENGTH_OF_TMP_STRING];

	bzero((void *)tmpString,LENGTH_OF_TMP_STRING);

	strcpy(dataParameters->nameOfFileWhereDictionaryWillBeDroped,dataParameters->nameOfOutputDirectory);
	strcpy(dataParameters->nameOfFileWhereFitedGaborsWillBeDroped,dataParameters->nameOfOutputDirectory);

	if(dot==NULL)
    	{		
		strcpy(dataParameters->nameOfFileWhereDictionaryWillBeDroped,dataParameters->nameOfDataFile);
		strcpy(dataParameters->nameOfFileWhereFitedGaborsWillBeDroped,dataParameters->nameOfDataFile);
	}
	else
	{
		len = strlen(dot);
	    	strncat(dataParameters->nameOfFileWhereDictionaryWillBeDroped,dataParameters->nameOfDataFile,strlen(dataParameters->nameOfDataFile)-len);
	    	strncat(dataParameters->nameOfFileWhereFitedGaborsWillBeDroped,dataParameters->nameOfDataFile,strlen(dataParameters->nameOfDataFile)-len);
	}

	sprintf(tmpString,"_%1.3f_%hu.dic",gaborDictionary->dilationFactor,gaborDictionary->periodDensity);
	strcat(dataParameters->nameOfFileWhereDictionaryWillBeDroped,tmpString);

	sprintf(tmpString,"_%1.3f_%hu.gab",gaborDictionary->dilationFactor,gaborDictionary->periodDensity);
	strcat(dataParameters->nameOfFileWhereFitedGaborsWillBeDroped,tmpString);

	for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
		strcpy(*(namesOfResultFiles + channelNumber),dataParameters->nameOfOutputDirectory);

	if(dot==NULL)
    	{
		for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
		{
			strcpy(*(namesOfResultFiles + channelNumber),dataParameters->nameOfDataFile);
			if(strcmp(dataParameters->extensionOfResultFile,"NONE")==0)
			{
				if(!(mp5Parameters->MPType & SMP))
					sprintf(*(namesOfResultFiles + channelNumber),"_mmp.b");
				else
					sprintf(*(namesOfResultFiles + channelNumber),"_ch_%hu.b",*(chosenChannels + channelNumber));
			}
			else
			{
				if(!(mp5Parameters->MPType & SMP))
					sprintf(*(namesOfResultFiles + channelNumber),"_%s_mmp.b",dataParameters->extensionOfResultFile);
				else
					sprintf(*(namesOfResultFiles + channelNumber),"_%s_ch_%hu.b",dataParameters->extensionOfResultFile,*(chosenChannels + channelNumber));        
			}
		}
    	}
    	else
    	{
		len = strlen(dot);
		for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
		{
	    		strncat(*(namesOfResultFiles + channelNumber),dataParameters->nameOfDataFile,strlen(dataParameters->nameOfDataFile)-len);
	    		lengthOfDataFileWithOutExpand = (unsigned short int)strlen(*(namesOfResultFiles + channelNumber));

	    		if(strcmp(dataParameters->extensionOfResultFile,"NONE")==0)
	    		{
				if(!(mp5Parameters->MPType & SMP))

					sprintf((*(namesOfResultFiles + channelNumber) + lengthOfDataFileWithOutExpand),"_mmp.b");
				else
		    			sprintf((*(namesOfResultFiles + channelNumber) + lengthOfDataFileWithOutExpand),"_ch_%hu.b",*(chosenChannels + channelNumber));
			}	       
	    		else
	    		{
				if(!(mp5Parameters->MPType & SMP))
					sprintf((*(namesOfResultFiles + channelNumber) + lengthOfDataFileWithOutExpand),"_%s_mmp.b",dataParameters->extensionOfResultFile);
				else
		    			sprintf((*(namesOfResultFiles + channelNumber) + lengthOfDataFileWithOutExpand),"_%s_ch_%hu.b",dataParameters->extensionOfResultFile,*(chosenChannels + channelNumber));	    	      
			}
		}
	}

}

STATUS openResultFiles(DataParameters *dataParameters, char *info)
{
	char **namesOfResultFiles = dataParameters->namesOfResultFiles;
	unsigned short int channelNumber;
	const unsigned short int numberOfResultsFiles = dataParameters->numberOfResultsFiles;

	if(dataParameters->verbose & VERBOSE_PRINT_DICTIONARY)
		if(!(dataParameters->dictionaryFile = fopen(dataParameters->nameOfFileWhereDictionaryWillBeDroped,"wt")))
		{
			sprintf(info, "\n DATA PARAMETERS ERROR: \n CAN'T OPEN FILE: %s WHERE DICTIONARY WILL BE DROPED ",dataParameters->nameOfFileWhereDictionaryWillBeDroped);
			return ERROR;
		}

	if(dataParameters->verbose & VERBOSE_PRINT_FITED_GABORS)
		if(!(dataParameters->fitedGaborsFile = fopen(dataParameters->nameOfFileWhereFitedGaborsWillBeDroped,"wt")))
		{
			sprintf(info, "\n DATA PARAMETERS ERROR: \n CAN'T OPEN FILE: %s WHERE FITED GABORS WILL BE DROPED ",dataParameters->nameOfFileWhereFitedGaborsWillBeDroped);
			return ERROR;
		}

	for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
	{
		if(dataParameters->writingMode & CREATE_FILE)
		{
			if(!(*(dataParameters->resultFiles + channelNumber) = fopen(*(namesOfResultFiles + channelNumber),"wb")))
			{
				sprintf(info, "\n DATA PARAMETERS ERROR: \n CAN'T OPEN RESULTS FILE: %s  ",*(dataParameters->namesOfResultFiles + channelNumber));
				return ERROR;
			}
		}
		else if(dataParameters->writingMode & APPEND_FILE)
		{
			if(!(*(dataParameters->resultFiles + channelNumber) = fopen(*(namesOfResultFiles + channelNumber),"a+b")))
			{
				sprintf(info, "\n DATA PARAMETERS ERROR: \n CAN'T OPEN RESULTS FILE: %s  ",*(dataParameters->namesOfResultFiles + channelNumber));
				return ERROR;
			}
		}
	}
	return SUCCESS;
}

void closeFiles(DataParameters *dataParameters)
{
	unsigned short int channelNumber;
	const unsigned short int numberOfResultsFiles = dataParameters->numberOfResultsFiles;

	if(dataParameters->dictionaryFile!=NULL)
		fclose(dataParameters->dictionaryFile);

	if(dataParameters->fitedGaborsFile!=NULL)
		fclose(dataParameters->fitedGaborsFile);

	if(dataParameters->dataFile!=NULL)
		fclose(dataParameters->dataFile);

	for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
	{
		if(dataParameters->allocatedElements & RESULT_FILES_ALLOCATED)
			if((*(dataParameters->resultFiles + channelNumber))!=NULL)
				fclose(*(dataParameters->resultFiles + channelNumber));
	}
}

STATUS analyseBinaryDataFile(DataParameters *dataParameters, char *info)
{
    unsigned long int sizeOfFile;
    unsigned long int sizeOfData;

    fseek(dataParameters->dataFile,0L,SEEK_END);
    sizeOfFile = ftell(dataParameters->dataFile);
    fseek(dataParameters->dataFile,0L,SEEK_SET);

    sizeOfData = sizeOfFile - dataParameters->sizeOfHeader - dataParameters->sizeOfTail;

    if(dataParameters->dataFormat & FORMAT_SHORT)
    {

	if((sizeOfData/sizeof(short))%dataParameters->numberOfChannels!=0)
	{
	    printf("\n DATA WARNING: \n CHANNELS DO NOT CONTAIN THE SAME NUMBER OF SAMPLES IN FILE: %s\n THE ERROR WILL OCCUR ERROR WHILE READING THE LAST OFFSET\n",dataParameters->nameOfDataFile);
	    fflush(stdout);
	    sleep(1);
	}

	if(sizeOfData%sizeof(short)!=0)
	{
	    sprintf(info,"\n DATA ERROR: \n INCORRECT NUMBER OF SAMPLES IN FILE: %s\n",dataParameters->nameOfDataFile);
	    return ERROR;
	}
	else
	{
	    dataParameters->numberOfPoints = sizeOfData/sizeof(short)/dataParameters->numberOfChannels;
	    dataParameters->numberOfOffsets = (unsigned short int)((dataParameters->numberOfPoints/dataParameters->numberOfChannels)/dataParameters->numberOfPointsInOffset);

	    if(dataParameters->numberOfOffsets < dataParameters->numberOfChosenOffsets)
	    {
		printf("\n DATA WARNING: \n TO LESS SAMPLES IN DATA FILE: %s\n NUMBER OF CHOSEN OFFSETS > NUMBER OF OFFSETS IN FILE \n",dataParameters->nameOfDataFile);
		fflush(stdout);
		sleep(1);
	    }

	    dataParameters->samplesBesideOffsets = dataParameters->numberOfPoints - dataParameters->numberOfOffsets*dataParameters->numberOfPoints;
	}
    }
    else if(dataParameters->dataFormat & FORMAT_FLOAT)
    {
	if((sizeOfData/sizeof(float))%dataParameters->numberOfChannels!=0)
	{
	    printf("\n DATA WARNING: \n CHANNELS DO NOT CONTAIN THE SAME NUMBER OF SAMPLES IN FILE: %s\n THE ERROR WILL OCCUR ERROR WHILE READING THE LAST OFFSET\n",dataParameters->nameOfDataFile);
	    fflush(stdout);
	    sleep(1);
	}

	if(sizeOfData%sizeof(float)!=0)
	{
	    sprintf(info,"\n DATA ERROR: \n INCORRECT NUMBER OF SAMPLES IN FILE: %s\n",dataParameters->nameOfDataFile);
	    return ERROR;
	}
	else
        {
	    dataParameters->numberOfPoints  = sizeOfData/sizeof(float)/dataParameters->numberOfChannels;
	    dataParameters->numberOfOffsets = (unsigned short int)((dataParameters->numberOfPoints/dataParameters->numberOfChannels)/dataParameters->numberOfPointsInOffset);

	    if(dataParameters->numberOfOffsets < dataParameters->numberOfChosenOffsets)
	    {
		printf("\n DATA ERROR: \n TO LESS SAMPLES IN DATA FILE: %s\n NUMBER OF CHOSEN OFFSETS > NUMBER OF OFFSETS IN FILE \n",dataParameters->nameOfDataFile);
		fflush(stdout);
		sleep(1);
	    }

	    dataParameters->samplesBesideOffsets = dataParameters->numberOfPoints - dataParameters->numberOfOffsets*dataParameters->numberOfPoints;
	}
    }
    
    return SUCCESS;
}

STATUS analyseAsciiDataFile(DataParameters *dataParameters, char *info)
{
    double value;
    unsigned short int channelNumber;
    unsigned long  int lineNumber;
    unsigned long  int numberOfLinesInFile;
    unsigned long  int numberOfDataLines;
    unsigned short int numberOfChannels = dataParameters->numberOfChannels;

    fseek(dataParameters->dataFile,0L,SEEK_SET);

    numberOfLinesInFile = 0UL;
    while(!feof(dataParameters->dataFile))
    {
	fscanf(dataParameters->dataFile,"%*[^\n]\n");
	numberOfLinesInFile++;
    }

    numberOfDataLines = numberOfLinesInFile - dataParameters->sizeOfHeader - dataParameters->sizeOfTail;

    asciiFileSeek(dataParameters->dataFile,dataParameters->sizeOfHeader);

    for(lineNumber=0;lineNumber<numberOfDataLines;lineNumber++)
    {
	for(channelNumber=0;channelNumber<numberOfChannels-1;channelNumber++)
	{
	    if(fscanf(dataParameters->dataFile,"%lf ",&value)!=1)
	    {
		sprintf(info,"\n DATA ERROR: \n CAN READ NOT SAMPLE IN DATA FILE: %s\n LINE: %lu, CHANNEL: %hu \n",dataParameters->nameOfDataFile,dataParameters->sizeOfHeader+lineNumber,channelNumber);
		return ERROR;
	    }
	}

	if(fscanf(dataParameters->dataFile,"%lf\n",&value)!=1)
	{	
	    sprintf(info,"\n DATA ERROR: \n CAN NOT READ SAMPLE IN DATA FILE: %s\n LINE: %lu, CHANNEL: %hu \n",dataParameters->nameOfDataFile,dataParameters->sizeOfHeader+lineNumber,channelNumber);
	    return ERROR;
	}
    }

    dataParameters->numberOfPoints  = numberOfDataLines;
    dataParameters->numberOfOffsets = (unsigned short int)((dataParameters->numberOfPoints)/dataParameters->numberOfPointsInOffset);

    if(dataParameters->numberOfOffsets < dataParameters->numberOfChosenOffsets)
    {
	printf("\n DATA ERROR: \n TO LESS SAMPLES IN DATA FILE: %s\n NUMBER OF CHOSEN OFFSETS > NUMBER OF OFFSETS IN FILE \n",dataParameters->nameOfDataFile);	
	fflush(stdout);
	sleep(1);
    }

    dataParameters->samplesBesideOffsets = dataParameters->numberOfPoints - dataParameters->numberOfOffsets*dataParameters->numberOfPoints;

    return SUCCESS;
}

void processRawData(DataParameters *dataParameters)
{
    const unsigned short int numberOfChannels  = dataParameters->numberOfChannels;
    const unsigned short int dimOffset         = dataParameters->dimOffset;
    const unsigned       int dimExpand         = dataParameters->dimExpand;
    double               **rawDataMatrix       = dataParameters->rawDataMatrix;
    double               **processedDataMatrix = dataParameters->processedDataMatrix;
    unsigned short int channel;
    unsigned       int sample;

    dSetMatrixZero(processedDataMatrix,numberOfChannels,dimExpand);

    for(channel=0;channel<numberOfChannels;channel++)
	for(sample=0;sample<dimOffset;sample++)
	    *(*(processedDataMatrix + channel) + dimOffset + sample) = *(*(rawDataMatrix + channel) + sample);
}

STATUS writeHeader(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary, char *info)
{
    unsigned short int channelNumber;
    const unsigned short int numberOfResultsFiles = dataParameters->numberOfResultsFiles;

    FILE_HEADER        fileHeader;
    decomposition_info dec_info;
    signal_info        sig_info;

    initField(&fileHeader);

    addTextInfo(&fileHeader,"http://brain.fuw.edu.pl/mp");

    dec_info.energy_percent           =  (float)mp5Parameters->energyPercent;
    dec_info.max_number_of_iterations =  mp5Parameters->maxNumberOfIterations;
    dec_info.dictionary_size          =  gaborDictionary->sizeOfDictionary;
    dec_info.dictionary_type          =  (char)((gaborDictionary->typeOfDictionary & OCTAVE_FIXED) ? 'D' : 'S');
    addDecompInfo(&fileHeader,&dec_info);

    sig_info.sampling_freq             = (float)dataParameters->samplingFrequency;
    sig_info.points_per_microvolt      = (float)dataParameters->convRate;
    sig_info.number_of_chanels_in_file = dataParameters->numberOfChannels;
    addSignalInfo(&fileHeader,&sig_info);

    addDate(&fileHeader);

    for(channelNumber=0;channelNumber<numberOfResultsFiles;channelNumber++)
    {
	if(dataParameters->writingMode & CREATE_FILE)
	{
	    if(WriteFileHeader(&fileHeader,*(dataParameters->resultFiles + channelNumber))==-1)
	    {
		sprintf(info,"\n MP5Parameters ERROR: \n CAN'T WRITE HEADER TO RESULTS FILE: %s \n",*(dataParameters->namesOfResultFiles + channelNumber));
		freeAllFields(&fileHeader);
		return ERROR;
	    }
	    fflush(*(dataParameters->resultFiles + channelNumber));
	}
	else if(dataParameters->writingMode & APPEND_FILE)
	{
	    fseek(*(dataParameters->resultFiles + channelNumber),0,SEEK_END);

	    if(ftell(*(dataParameters->resultFiles + channelNumber))==0)
	    {
		if(WriteFileHeader(&fileHeader,*(dataParameters->resultFiles + channelNumber))==-1)
		{
		    sprintf(info,"\n MP5Parameters ERROR: \n CAN'T WRITE HEADER TO RESULTS FILE: %s \n",*(dataParameters->namesOfResultFiles + channelNumber));
		    freeAllFields(&fileHeader);
		    return ERROR;
		}
		fflush(*(dataParameters->resultFiles + channelNumber));
	    }
	}
    }
    freeAllFields(&fileHeader);

    return SUCCESS;
}


STATUS writeSingleChannelResults(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary, unsigned short int offsetNumber, unsigned short int channelNumber, char *info)
{

    float  amplitude;
    float  phase;
    float  modulus;
    double frequency;

    unsigned short int dimOffset = mp5Parameters->dimOffset;

    SEG_HEADER  head;
    Gabor       *gabor;
    NEW_ATOM	atom;

    head.channel       = *(dataParameters->chosenChannels + channelNumber);
    head.file_offset   = *(dataParameters->chosenOffsets  + offsetNumber);
    head.book_size     = mp5Parameters->fitted->size;
    head.signal_size   = dataParameters->numberOfPointsInOffset;
    head.signal_energy = (float)(*(mp5Parameters->signalEnergyInEachChannel));
    head.book_energy   = (float)(*(mp5Parameters->signalEnergyInEachChannel) - *(mp5Parameters->residueEnergyInEachChannel));

    if(WriteSegmentHeader(&head,*(dataParameters->resultFiles + channelNumber))==-1)
    {
	sprintf(info,"\n mp5Parameters ERROR: \n CAN'T WRITE SEGMENT TO FILE: %s",*(dataParameters->namesOfResultFiles + channelNumber));
	return ERROR;
    }

    do
    {
	gabor = (Gabor *)substractNode(mp5Parameters->fitted);
	phase = *(gabor->phase);

	returnAmplitudeAndModulusDI(gabor,&amplitude,&modulus,0);

	atom.modulus   = modulus;
	atom.position  = (float)(gabor->position);

	if(gabor->feature & DIRACDELTA)
	    atom.scale = 0;
	else
	    atom.scale = (float)(*(gaborDictionary->tableOfScalesInOptimalDictionary + gabor->scaleIndex));

        if(gabor->feature & DIRACDELTA)
	    frequency = dimOffset/2.0;
	else
	    frequency = 0.5*dimOffset*(*(gaborDictionary->tableOfFrequenciesInOptimalDictionary + gabor->scaleIndex))*gabor->rifling/M_PI;

	atom.frequency = (float)(frequency);
	atom.phase     = phase;
	atom.amplitude = (float)(2.0*amplitude);

	if(WriteNewAtom(&atom,*(dataParameters->resultFiles + channelNumber))==-1)
	{
	    sprintf(info,"\n mp5Parameters ERROR: \n CAN'T WRITE SEGMENT TO FILE: %s\n",*(dataParameters->namesOfResultFiles + channelNumber));
	    return ERROR;
	}

	fflush(*(dataParameters->resultFiles + channelNumber));

	freeGabor(gabor);

    }while(mp5Parameters->fitted->size>0);

    return SUCCESS;

}


STATUS writeMultiChannelResults(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary, unsigned short int offsetNumber, char *info)
{
    unsigned short int channelNumber;
    unsigned int       atomNumber;
    float  amplitude;
    float  phase;
    float  modulus;
    double frequency;

    const unsigned short int dimOffset = mp5Parameters->dimOffset;

    SEG_HEADER  head;
    Gabor       *gabor;
    NEW_ATOM	atom;

    for(channelNumber=0;channelNumber<mp5Parameters->numberOfAnalysedChannels;channelNumber++)
    {
	head.channel       = *(dataParameters->chosenChannels + channelNumber);
	head.file_offset   = *(dataParameters->chosenOffsets  + offsetNumber);
	head.book_size     = mp5Parameters->fitted->size;
	head.signal_size   = dataParameters->numberOfPointsInOffset;
	head.signal_energy = (float)(*(mp5Parameters->signalEnergyInEachChannel + channelNumber));
	head.book_energy   = (float)(*(mp5Parameters->signalEnergyInEachChannel + channelNumber) - *(mp5Parameters->residueEnergyInEachChannel + channelNumber));

	if(WriteSegmentHeader(&head,*(dataParameters->resultFiles))==-1)
	{
	    sprintf(info,"\n mp5Parameters ERROR: \n CAN'T WRITE SEGMENT TO FILE: %s",*(dataParameters->namesOfResultFiles + channelNumber));
	    return ERROR;
	}
   
	fflush(*(dataParameters->resultFiles));

	for(atomNumber=0;atomNumber<mp5Parameters->fitted->size;atomNumber++)
	{

	    gabor = (Gabor *)readNNode(mp5Parameters->fitted,atomNumber);

		if(mp5Parameters->MPType & MMP2)
		{
			phase = *(gabor->phase);
			returnAmplitudeAndModulusForMMP2DI(mp5Parameters,gaborDictionary,gabor,&amplitude,&modulus,channelNumber);
		}
		else
		{
			phase = *(gabor->phase+channelNumber);
			returnAmplitudeAndModulusDI(gabor,&amplitude,&modulus,channelNumber);
		}
		
	    atom.modulus   = modulus;
	    atom.position  = (float)(gabor->position);

	    if(gabor->feature & DIRACDELTA)
		atom.scale = 0;
	    else
		atom.scale = (float)(*(gaborDictionary->tableOfScalesInOptimalDictionary + gabor->scaleIndex));
      
    	    if(gabor->feature & DIRACDELTA)
		frequency = dimOffset/2.0;
	    else
		frequency = 0.5*dimOffset*(*(gaborDictionary->tableOfFrequenciesInOptimalDictionary + gabor->scaleIndex))*gabor->rifling/M_PI;
      
	    atom.frequency = (float)(frequency);
	    atom.phase     = phase;
	    atom.amplitude = (float)(2.0*amplitude);

	    if(WriteNewAtom(&atom,*(dataParameters->resultFiles))==-1)
	    {
		sprintf(info,"\n mp5Parameters ERROR: \n CAN'T WRITE SEGMENT TO FILE: %s\n",*(dataParameters->namesOfResultFiles + channelNumber));
		return ERROR;
	    }
	    fflush(*(dataParameters->resultFiles));
	}
    }

    clearQueue(mp5Parameters->fitted,(void (*)(void *))freeGabor);

    return SUCCESS;
}

STATUS readBinaryData(DataParameters *dataParameters, unsigned short int offsetNumber, char *info)
{
    double             **rawDataMatrix = dataParameters->rawDataMatrix;
    const unsigned short int numberOfChannels = dataParameters->numberOfChannels;
    const unsigned short int numberOfPointsInOffset = dataParameters->numberOfPointsInOffset;
    unsigned short int channel;
    unsigned short int sample;
    unsigned short int formatSize = 0;
    unsigned long int  filePosition;

    if(dataParameters->dataFormat & FORMAT_SHORT)
    {
	short int tmpData[numberOfChannels];

	formatSize = sizeof(short int);

	filePosition = dataParameters->sizeOfHeader + (offsetNumber-1)*numberOfPointsInOffset*numberOfChannels*formatSize;

	fseek(dataParameters->dataFile,filePosition,SEEK_SET);

	for(sample=0;sample<numberOfPointsInOffset;sample++)
	{
	    if(fread((void *)tmpData,formatSize,numberOfChannels,dataParameters->dataFile)<numberOfChannels)
	    {
		sprintf(info,"\n ERROR WHILE READING DATA FILE: %s\n CAN NOT READ OFFSET NUMBER: %hu, SAMPLE: %hu \n",dataParameters->nameOfDataFile,offsetNumber,sample+1);
		return ERROR;
	    }

	    for(channel=0;channel<numberOfChannels;channel++)
		*(*(rawDataMatrix + channel) + sample) = (double)(*(tmpData + channel));
	}
    }
    else if(dataParameters->dataFormat & FORMAT_FLOAT)
    {
	float tmpData[numberOfChannels];

	formatSize = sizeof(float);

	filePosition = dataParameters->sizeOfHeader + (offsetNumber-1)*numberOfPointsInOffset*numberOfChannels*formatSize;

	fseek(dataParameters->dataFile,filePosition,SEEK_SET);

	for(sample=0;sample<numberOfPointsInOffset;sample++)
	{
	    if(fread((void *)tmpData,formatSize,numberOfChannels,dataParameters->dataFile)<numberOfChannels)
	    {
		sprintf(info,"\n ERROR WHILE READING DATA FILE: %s\n CAN NOT READ OFFSET NUMBER: %hu, SAMPLE: %hu \n",dataParameters->nameOfDataFile,offsetNumber,sample+1);
		return ERROR;
	    }

	    for(channel=0;channel<numberOfChannels;channel++)
		*(*(rawDataMatrix + channel) + sample) = (double)(*(tmpData + channel));
	}
    }

    return SUCCESS;
}

STATUS readAsciiData(DataParameters *dataParameters, unsigned short int offsetNumber, char *info)
{
    double              **rawDataMatrix = dataParameters->rawDataMatrix;
    const unsigned short int numberOfChannels = dataParameters->numberOfChannels;
    const unsigned short int numberOfPointsInOffset = dataParameters->numberOfPointsInOffset;
    unsigned short int channelNumber;
    unsigned short int lineNumber;

    asciiFileSeek(dataParameters->dataFile,dataParameters->sizeOfHeader + (offsetNumber-1)*numberOfPointsInOffset);

    for(lineNumber=0;lineNumber<numberOfPointsInOffset;lineNumber++)
    {
	for(channelNumber=0;channelNumber<numberOfChannels-1;channelNumber++)
	{
	    if(fscanf(dataParameters->dataFile,"%lf ",(*(rawDataMatrix + channelNumber) + lineNumber))!=1)
	    {
		sprintf(info,"\n ERROR WHILE READING DATA FILE: %s\n CAN NOT SAMPLE IN OFFSET NUMBER: %hu LINE NUMBER: %hu, CHANNEL: %hu \n",dataParameters->nameOfDataFile,offsetNumber,lineNumber+1,channelNumber+1);
		return ERROR;
	    }
	}

	if(fscanf(dataParameters->dataFile,"%lf\n",(*(rawDataMatrix + channelNumber) + lineNumber))!=1)
	{
	    sprintf(info,"\n ERROR WHILE READING DATA FILE: %s\n CAN NOT SAMPLE IN OFFSET NUMBER: %hu LINE NUMBER: %hu, CHANNEL: %hu \n",dataParameters->nameOfDataFile,offsetNumber,lineNumber+1,channelNumber+1);
	    return ERROR;
	}
    }
    return SUCCESS;
}

STATUS readDataFile(DataParameters *dataParameters, unsigned short int offsetNumber, char *info)
{
    if(dataParameters->dataFormat & FORMAT_ASCII)
    {
	if(readAsciiData(dataParameters,offsetNumber,info)==ERROR)
	    return ERROR;
    }
    else if((dataParameters->dataFormat & FORMAT_FLOAT) || (dataParameters->dataFormat & FORMAT_SHORT))
    {
	if(readBinaryData(dataParameters,offsetNumber,info)==ERROR)
	    return ERROR;
    }

    /* copy data matrix (rawDataMatrix) of size numberOfChannels x dimOffset to
       processDataMatrix of size numberOfChannels x dimExpand */

    processRawData(dataParameters);

    return SUCCESS;
}
