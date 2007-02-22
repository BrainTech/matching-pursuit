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
#include<string.h>
#include"include/cmd.h"
#include"include/def.h"
#include"include/queue.h"
#include"include/stringTools.h"
#include"include/vector.h"

static struct CommandsList
{
	char    *command;
	BOOLEAN found;
} commandsList[NUMBER_OF_COMMANDS] =
	{
		{"nameOfDataFile",FALSE},
		{"extensionOfResultFile",FALSE},
		{"nameOfOutputDirectory",FALSE},
		{"writingMode",FALSE},
		{"sizeOfHeader",FALSE},
		{"sizeOfTail",FALSE},
		{"samplingFrequency",FALSE},
		{"formatOfData",FALSE},
		{"numberOfChannels",FALSE},
		{"chosenChannels",FALSE},
		{"numberOfPointsInOffset",FALSE},
		{"chosenOffsets",FALSE},
		{"typeOfDictionary",FALSE},
		{"dilationFactor",FALSE},
		{"periodDensity",FALSE},
		{"reinitDictionary",FALSE},
		{"scaleToPeriodFactor",FALSE},
		{"DOT_EPS",FALSE},
		{"maxNumberOfIterations",FALSE},
		{"energyPercent",FALSE},
		{"convRate",FALSE},
		{"MP",FALSE},
		{"VERBOSE",FALSE}
  };

static STATUS checkIfLineBeginWithCommand(const char *line)
{
	unsigned short int counter;
	char tmpLine[LENGTH_OF_LINE];
	char *command;

	BOOLEAN commandExists = FALSE;

	strcpy(tmpLine,line);

	command = (char *)strtok(tmpLine," ");

	for(counter=0;counter<NUMBER_OF_COMMANDS;counter++)
	{
		if(strcmp((commandsList[counter]).command,command)==0)
			commandExists = TRUE;
	}

	return commandExists==TRUE ? SUCCESS : ERROR;
}

static STATUS isLineBrokenCorrect(const char *line)
{
  char *positionOfBroken = NULL;
  unsigned short int numberOfWhiteBlanks;

  positionOfBroken = strchr(line,'\\');

  numberOfWhiteBlanks = countWhiteBlanks(positionOfBroken);

  if((strlen(positionOfBroken)-1)!=numberOfWhiteBlanks)
    return ERROR;

  return SUCCESS;
}

static void pushString(Queue *stringQueue, char *text)
{
  String *string = (String *)malloc(sizeof(String));

  strcpy(string->text,text);
  addNode(stringQueue,(void *)string);
}

static void popString(Queue *stringQueue, char *text)
{
  String *string;

  string = (String *)substractNode(stringQueue);
  strcpy(text,string->text);

  free(string);
}
/*
  static void updateString(Queue *stringQueue, char *text)
  {
  String *string;

  string = (String *)(stringQueue->lastNode->data);
  strcat(string->text,text);
  }*/

static void pushLine(Queue *lineQueue, char *text, int number)
{
  Line *line = (Line *)malloc(sizeof(Line));
  strcpy(line->text,text);
  line->number = number;
  addNode(lineQueue,(void *)line);
}

static void popLine(Queue *lineQueue, char *text, int *number)
{
  Line *line;

  line = (Line *)substractNode(lineQueue);

  strcpy(text,line->text);
  *number = line->number;
  free(line);
}

static void updateLine(Queue *lineQueue, char *text)
{
  Line *line;

  line = (Line *)(lineQueue->lastNode->data);
  strcat(line->text,text);
}

STATUS openConfigFile(ConfigFile *configFile, char *info)
{

    configFile->file = fopen(configFile->name,"rt");

    if(configFile->file==NULL)
    {
	sprintf(info," FILE OPEN ERROR \n CAN NOT OPEN CONFIG FILE: %s \n",configFile->name);
	return ERROR;
    }
    return SUCCESS;
}

void setConfigFile(ConfigFile *configFile)
{
  configFile->stringQueue = createQueue();
  configFile->lineQueue   = createQueue();
}

void freeConfigFile(ConfigFile *configFile)
{
  freeQueue(configFile->stringQueue,NULL);
  freeQueue(configFile->lineQueue,NULL);
  if(configFile!=NULL)
    fclose(configFile->file);
}

STATUS readConfigFile(ConfigFile *configFile, char *info)
{
    int  lineNumber     = 0;
    char text[LENGTH_OF_LINE];
    char *positionOfBroken = NULL;
    unsigned short int numberOfWhiteBlanks;    
    size_t lengthOfText;

    BOOLEAN addLine = FALSE;
    STATUS  status = SUCCESS;

    do
    {
	lineNumber++;

	if(fgets(text,LENGTH_OF_LINE,configFile->file)==NULL)
	    break;
	
	numberOfWhiteBlanks = countWhiteBlanks(text);
	lengthOfText        = strlen(text); 

	if((strchr(text,'#')!=NULL) || (numberOfWhiteBlanks == lengthOfText))
	    continue;
	else
	{
	    if(addLine)
	    {
		if(strchr(text,'\\')!=NULL)
		{
		    if(!isLineBrokenCorrect(text))
		    {
			sprintf(info,"\n CONFIG FILE COMMANDS ERROR: \
						              \n LINE NUMBER: %d             \
						              \n IS NOT BROKEN CORRECT       \
							      \n CHECK THIS FILE %s          \
							      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",lineNumber,configFile->name);
			status = ERROR;
			break;
		    }

		    positionOfBroken = (char *)strtok(text,"\\");
		    updateLine(configFile->lineQueue,positionOfBroken);
		    addLine = TRUE;
		}
		else
		{
		    updateLine(configFile->lineQueue,positionOfBroken);
		    addLine = FALSE;
		}
	    }
	    else
	    {
		if(strchr(text,'\\')!=NULL)
		{
		    if(!isLineBrokenCorrect(text))
		    {
			sprintf(info,"\n CONFIG FILE COMMANDS ERROR: \
						              \n LINE NUMBER: %d             \
						              \n IS NOT BROKEN CORRECT       \
							      \n CHECK THIS FILE: %s \n      \
							      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",lineNumber,configFile->name);

			status = ERROR;
			break;
		    }

		    positionOfBroken = (char *)strtok(text,"\\");
		    pushLine(configFile->lineQueue,text,lineNumber);
		    addLine = TRUE;
		}
		else
		    pushLine(configFile->lineQueue,text,lineNumber);
	    }
	}
    }while(TRUE);

    return status;

}

STATUS findDataParametersInConfigFile(ConfigFile *configFile, DataParameters *dataParameters, GaborDictionary *gaborDictionary, MP5Parameters *mp5Parameters, char *info)
{

	int i;
	int lineNumber;

	char text[LENGTH_OF_LINE];
	char *parameterPosition = NULL;

	while((configFile->lineQueue->firstNode)!=NULL)
	{
		popLine(configFile->lineQueue,text,&lineNumber);

		if(checkIfLineBeginWithCommand(text)==ERROR)
			goto ERROR_PROCEDURE_1;

		if(strstr(text,"nameOfDataFile")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			commandsList[0].found = TRUE;

			strcpy(dataParameters->nameOfDataFile,parameterPosition);
		}

		if(strstr(text,"extensionOfResultFile")!=NULL)
		{
			
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			commandsList[1].found = TRUE;

			strcpy(dataParameters->extensionOfResultFile,parameterPosition);
		}
		if(strstr(text,"nameOfOutputDirectory")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");
			
			if(*(parameterPosition + (strlen(parameterPosition)-1))!='/')
			    goto ERROR_PROCEDURE_3;

			commandsList[2].found = TRUE;
			strcpy(dataParameters->nameOfOutputDirectory,parameterPosition);
		}
		if(strstr(text,"writingMode")!=NULL)
		{
			dataParameters->writingMode = 0x0;

			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			if(strcmp(parameterPosition,"CREATE")!=0 &&
			strcmp(parameterPosition,"APPEND")!=0) goto ERROR_PROCEDURE_3;

			if(strcmp(parameterPosition,"CREATE")==0)      dataParameters->writingMode|=CREATE_FILE;
			else if(strcmp(parameterPosition,"APPEND")==0) dataParameters->writingMode|=APPEND_FILE;

			commandsList[3].found = TRUE;
		}
		else if(strstr(text,"sizeOfHeader")!=NULL)
		{

			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			if(!isDecimal(parameterPosition,"uint")) goto ERROR_PROCEDURE_4;

			commandsList[4].found = TRUE;

			dataParameters->sizeOfHeader = (unsigned short int)atoi(parameterPosition);
		}
		else if(strstr(text,"sizeOfTail")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uint")) goto ERROR_PROCEDURE_4;

			commandsList[5].found = TRUE;

			dataParameters->sizeOfTail = (unsigned short int)atoi(parameterPosition);
		}
		else if(strstr(text,"samplingFrequency")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz")) goto ERROR_PROCEDURE_7;

			commandsList[6].found = TRUE;

			dataParameters->samplingFrequency = atof(parameterPosition);
		}
		else if(strstr(text,"formatOfData")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			dataParameters->dataFormat = 0x0;

			if(strcmp(parameterPosition,"ASCII")!=0 &&

			strcmp(parameterPosition,"SHORT")!=0 &&
			strcmp(parameterPosition,"FLOAT")!=0) goto ERROR_PROCEDURE_3;

			if(strcmp(parameterPosition,"ASCII")==0)      dataParameters->dataFormat|= FORMAT_ASCII;
			if(strcmp(parameterPosition,"SHORT")==0)      dataParameters->dataFormat|= FORMAT_SHORT;
			else if(strcmp(parameterPosition,"FLOAT")==0) dataParameters->dataFormat|= FORMAT_FLOAT;

			commandsList[7].found = TRUE;
		}
		else if(strstr(text,"numberOfChannels")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz")) goto ERROR_PROCEDURE_5;

			commandsList[8].found = TRUE;

			dataParameters->numberOfChannels = (unsigned short int)atoi(parameterPosition);
		}
		else if(strstr(text,"chosenChannels")!=NULL)
		{
			char backupText[LENGTH_OF_LINE];

			char string[LENGTH_OF_STRING];
			char *firstString;
			char *secondString;

			unsigned short int channel, channelsCounter = 0;
			unsigned short int firstChannel;
			unsigned short int secondChannel;

			strcpy(backupText,text);

			if(countWords(text)==1) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\-\n\t"))!=NULL)
				if(!isDecimal(parameterPosition,"uintgz")) goto ERROR_PROCEDURE_5;
					strcpy(backupText,text);

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\\n\t"))!=NULL)
				pushString(configFile->stringQueue,parameterPosition);

			while((configFile->stringQueue->firstNode)!=NULL)
			{
				popString(configFile->stringQueue,string);

				firstString  = (char *)strtok(string,"-");
				secondString = (char *)strtok(NULL,"-");

				if(!firstString) goto ERROR_PROCEDURE_3;
				else if(strcmp(firstString,"-")==0) goto ERROR_PROCEDURE_3;
				else
				{
					firstChannel = (unsigned short int)atoi(firstString);

					if(!secondString)
						dataParameters->numberOfChosenChannels++;
					else if(strcmp(secondString,"-")==0) goto ERROR_PROCEDURE_2;
					else
					{
						secondChannel = (unsigned short int)atoi(secondString);

						if(secondChannel<firstChannel) goto ERROR_PROCEDURE_3;
						else
							dataParameters->numberOfChosenChannels = (unsigned short int)(dataParameters->numberOfChosenChannels + secondChannel - firstChannel + 1);

					}
				}
			}

			dataParameters->chosenChannels = (unsigned short int *)fVectorAllocate(dataParameters->numberOfChosenChannels);

			strcpy(backupText,text);

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\\n\t"))!=NULL)
				pushString(configFile->stringQueue,parameterPosition);

			while((configFile->stringQueue->firstNode)!=NULL)
			{
				popString(configFile->stringQueue,string);

				firstString  = (char *)strtok(string,"-");
				secondString = (char *)strtok(NULL,"-");

				firstChannel = (unsigned short int)atoi(firstString);

				if(!secondString)
				{
					dataParameters->chosenChannels[channelsCounter] = firstChannel;
					channelsCounter++;
				}
				else
				{
					secondChannel = (unsigned short int)atoi(secondString);

					for(channel=0;channel<secondChannel - firstChannel + 1;channel++,channelsCounter++)
						dataParameters->chosenChannels[channelsCounter] = (unsigned short int)(firstChannel + channel);
				}
			}

			dataParameters->allocatedElements|=CHOSEN_CHANNELS_ALLOCATED;

			commandsList[9].found = TRUE;
		}
		else if(strstr(text,"numberOfPointsInOffset")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz")) goto ERROR_PROCEDURE_5;

			dataParameters->numberOfPointsInOffset = (unsigned short int)atoi(parameterPosition);

			mp5Parameters->dimOffset   = dataParameters->dimOffset = dataParameters->numberOfPointsInOffset;
			mp5Parameters->dimExpTable = 2*mp5Parameters->dimOffset;
			mp5Parameters->dimExpand   = dataParameters->dimExpand = 3*mp5Parameters->dimOffset;

			commandsList[10].found = TRUE;
		}
		else if(strstr(text,"chosenOffsets")!=NULL)
		{
			char backupText[LENGTH_OF_LINE];

			char string[LENGTH_OF_STRING];
			char *firstString;
			char *secondString;

			unsigned short int offset, offsetsCounter = 0;
			unsigned short int firstOffset;
			unsigned short int secondOffset;

			strcpy(backupText,text);

			if(countWords(text)==1) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\-\n\t"))!=NULL)
				if(!isDecimal(parameterPosition,"uintgz")) goto ERROR_PROCEDURE_5;

			strcpy(backupText,text);

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\\n\t"))!=NULL)
				pushString(configFile->stringQueue,parameterPosition);

			while((configFile->stringQueue->firstNode)!=NULL)
			{
				popString(configFile->stringQueue,string);

				firstString  = (char *)strtok(string,"-");
				secondString = (char *)strtok(NULL,"-");

				if(!firstString) goto ERROR_PROCEDURE_3;
				else if(strcmp(firstString,"-")==0) goto ERROR_PROCEDURE_3;
				else
				{
					firstOffset = (unsigned short int)atoi(firstString);

					if(!secondString)
						dataParameters->numberOfChosenOffsets++;
					else if(strcmp(secondString,"-")==0) goto ERROR_PROCEDURE_3;
					else
					{
						secondOffset = (unsigned short int)atoi(secondString);

						if(secondOffset<firstOffset) goto ERROR_PROCEDURE_3;
						else
							dataParameters->numberOfChosenOffsets = (unsigned short int)(dataParameters->numberOfChosenOffsets + secondOffset - firstOffset + 1);
					}
				}
			}

			dataParameters->chosenOffsets = (unsigned short int *)iVectorAllocate(dataParameters->numberOfChosenOffsets);

			strcpy(backupText,text);

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\\n\t"))!=NULL)
				pushString(configFile->stringQueue,parameterPosition);

			while((configFile->stringQueue->firstNode)!=NULL)
			{
				popString(configFile->stringQueue,string);

				firstString  = (char *)strtok(string,"-");
				secondString = (char *)strtok(NULL,"-");

				firstOffset = (unsigned short int)atoi(firstString);

				if(!secondString)
				{
					dataParameters->chosenOffsets[offsetsCounter] = firstOffset;
					offsetsCounter++;
				}
				else
				{
					secondOffset = (unsigned short int)atoi(secondString);

					for(offset=0;offset<secondOffset - firstOffset + 1;offset++,offsetsCounter++)
						dataParameters->chosenOffsets[offsetsCounter] = (unsigned short int)(firstOffset + offset);
				}
			}

			dataParameters->allocatedElements|=CHOSEN_OFFSETS_ALLOCATED;

			commandsList[11].found = TRUE;
		}
		else if(strstr(text,"typeOfDictionary")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"OCTAVE_FIXED")!=0 && strcmp(parameterPosition,"OCTAVE_STOCH")!=0)
				goto ERROR_PROCEDURE_3;

			if(strcmp(parameterPosition,"OCTAVE_FIXED")==0)
				gaborDictionary->typeOfDictionary|= OCTAVE_FIXED;
			else if(strcmp(parameterPosition,"OCTAVE_STOCH")==0)
				gaborDictionary->typeOfDictionary|= OCTAVE_STOCH;

			commandsList[12].found = TRUE;
		}
		else if(strstr(text,"dilationFactor")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz")) goto ERROR_PROCEDURE_7;

			gaborDictionary->dilationFactor = (double)atof(parameterPosition);
			commandsList[13].found = TRUE;
		}
		else if(strstr(text,"periodDensity")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz")) goto ERROR_PROCEDURE_5;

			gaborDictionary->periodDensity = (unsigned short int)atoi(parameterPosition);
			commandsList[14].found = TRUE;
		}
		else if(strstr(text,"reinitDictionary")!=NULL)
		{
			mp5Parameters->reinitDictionary = 0x0;

			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"NO_REINIT_AT_ALL")!=0 &&
			strcmp(parameterPosition,"REINIT_IN_CHANNEL_DOMAIN")!=0 &&
			strcmp(parameterPosition,"REINIT_IN_OFFSET_DOMAIN")!=0  &&
			strcmp(parameterPosition,"REINIT_AT_ALL")!=0) goto ERROR_PROCEDURE_3;

			if(strcmp(parameterPosition,"NO_REINIT_AT_ALL")==0)
				mp5Parameters->reinitDictionary|= NO_REINIT_AT_ALL;
			else if(strcmp(parameterPosition,"REINIT_IN_CHANNEL_DOMAIN")==0)
				mp5Parameters->reinitDictionary|= REINIT_IN_CHANNEL_DOMAIN;
			else if(strcmp(parameterPosition,"REINIT_IN_CHANNEL_DOMAIN")==0)
				mp5Parameters->reinitDictionary|= REINIT_IN_CHANNEL_DOMAIN;
			else if(strcmp(parameterPosition,"REINIT_IN_OFFSET_DOMAIN")==0)
				mp5Parameters->reinitDictionary|= REINIT_IN_OFFSET_DOMAIN;
			else if(strcmp(parameterPosition,"REINIT_AT_ALL")==0)
				mp5Parameters->reinitDictionary|= REINIT_AT_ALL;

			commandsList[15].found = TRUE;
		}
		else if(strstr(text,"scaleToPeriodFactor")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"real")) goto ERROR_PROCEDURE_6;

			gaborDictionary->scaleToPeriodFactor = (double)atof(parameterPosition);
			commandsList[16].found = TRUE;
		}
		else if(strstr(text,"DOT_EPS")!=NULL)
		{
			double tmpDouble;
			char tmpChar[LENGTH_OF_LINE];

			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			sscanf(parameterPosition,"%lE",&tmpDouble);
			sprintf(tmpChar,"%1.17lf",tmpDouble);

			if(!isReal(tmpChar,"realgz")) goto ERROR_PROCEDURE_7;

			if((double)atof(parameterPosition)<EPS_DOUBLE) goto ERROR_PROCEDURE_8;

			mp5Parameters->LOG_EPS_DOT_PRODUCT = log((double)atof(parameterPosition));
			commandsList[17].found = TRUE;
		}
		else if(strstr(text,"maxNumberOfIterations")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz")) goto ERROR_PROCEDURE_5;

			mp5Parameters->maxNumberOfIterations = (unsigned short int)atoi(parameterPosition);

			commandsList[18].found = TRUE;
		}
		else if(strstr(text,"energyPercent")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz")) goto ERROR_PROCEDURE_7;

			mp5Parameters->energyPercent = atof(parameterPosition);

			commandsList[19].found = TRUE;
		}
		else if(strstr(text,"MP")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"SMP")!=0  &&
			strcmp(parameterPosition,"MMP1")!=0 &&
			strcmp(parameterPosition,"MMP2")!=0 &&
			strcmp(parameterPosition,"MMP3")!=0) goto ERROR_PROCEDURE_3;

			if(strcmp(parameterPosition,"SMP")==0)       mp5Parameters->MPType|= SMP;
			else if(strcmp(parameterPosition,"MMP1")==0) mp5Parameters->MPType|= MMP1;
			else if(strcmp(parameterPosition,"MMP2")==0) mp5Parameters->MPType|= MMP2;
			else if(strcmp(parameterPosition,"MMP3")==0) mp5Parameters->MPType|= MMP3;

			commandsList[20].found = TRUE;
		}
		else if(strstr(text,"convRate")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz")) goto ERROR_PROCEDURE_7;

			dataParameters->convRate = atof(parameterPosition);

			commandsList[21].found = TRUE;
		}
		else if(strstr(text,"VERBOSE")!=NULL)
		{
			if(countWords(text)!=2) goto ERROR_PROCEDURE_2;

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz")) goto ERROR_PROCEDURE_5;

			dataParameters->verbose = (unsigned char)atoi(parameterPosition);

			commandsList[22].found = TRUE;
		}
	}


	for(i=0;i<NUMBER_OF_COMMANDS;i++)
	{
		if(commandsList[i].found == FALSE)
		{
			sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s \
		                      \n COMMAND: %s                    \
			              \n WAS NOT FOUND                  \
			              \n CHECK THIS FILE                \
			              \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,commandsList[i].command);

		return ERROR;
		}
	}

	return SUCCESS;

	ERROR_PROCEDURE_1:
		sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s                  \
		              \n LINE: %d                                        \
			      \n DO NOT CONTAIN COMMAND (PROPABLY IT IS COMMENT) \
			      \n CHECK THIS FILE                                 \
			      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,lineNumber);

		return ERROR;

	ERROR_PROCEDURE_2:
		sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s                  \
		              \n INCORRECT NUBER OF ARGUMENTS IN LINE: %d  \
			      \n CHECK THIS FILE                                 \
			      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,lineNumber);

		return ERROR;

	ERROR_PROCEDURE_3:
		sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s        \
		              \n INCORRECT ARGUMENT IN LINE: %d  \
			      \n CHECK THIS FILE                       \
			      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,lineNumber);
		return ERROR;

	ERROR_PROCEDURE_4:
		sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s        \
		              \n INCORRECT ARGUMENT IN LINE: %d  \
			      \n THE ARGUMENT SHOULD BE INTEGER >= 0   \
			      \n CHECK THIS FILE                       \
			      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,lineNumber);

		return ERROR;

	ERROR_PROCEDURE_5:
		sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s        \
		              \n INCORRECT ARGUMENT IN LINE: %d  \
			      \n THE ARGUMENT SHOULD BE INTEGER > 0    \
			      \n CHECK THIS FILE                       \
			      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,lineNumber);

		return ERROR;

	ERROR_PROCEDURE_6:
		sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s        \
		              \n INCORRECT ARGUMENT IN LINE: %d        \
			      \n THE ARGUMENT SHOULD BE FLOAT >= 0     \
			      \n CHECK THIS FILE                       \
			      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,lineNumber);

		return ERROR;

	ERROR_PROCEDURE_7:
		sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s        \
		              \n INCORRECT ARGUMENT IN LINE: %d        \
			      \n THE ARGUMENT SHOULD BE FLOAT > 0      \
			      \n CHECK THIS FILE                       \
			      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,lineNumber);

		return ERROR;

	ERROR_PROCEDURE_8:
		sprintf(info,"\n CONFIG FILE COMMANDS ERROR: %s                                                   \
		              \n DOUBLE ARGUMENT IN LINE %d IS SMALLER THEN REAL ACCURANCY OF DOUBLE TYPE (1E-16) \
			      \n CHECK THIS FILE                                                                  \
			      \n TO GET THE DEFAULT CONFIG FILE TYPE mp5 -g \n",configFile->name,lineNumber);

		return ERROR;

}
