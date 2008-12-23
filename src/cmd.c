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
#include"cmd.h"
#include"def.h"
#include"io_mp5.h"
#include"queue.h"
#include"stringTools.h"
#include"vector.h"

static struct CommandsList
{
	char    *command;
	BOOLEAN found;
} commandsList[NUMBER_OF_ALL_COMMANDS] =
	{
		{"nameOfDataFile",FALSE},
		{"extensionOfResultsFile",FALSE},
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
		{"randomSeed",FALSE},
		{"periodDensity",FALSE},
		{"reinitDictionary",FALSE},
		{"scaleToPeriodFactor",FALSE},
		{"maximalNumberOfIterations",FALSE},
		{"energyPercent",FALSE},
		{"MP",FALSE},
		{"analiticalDotProduct",FALSE},
		{"bookWithSignal",FALSE},
		{"pointsPerMicrovolt",FALSE},
		{"FFT",FALSE},
		{"accuracy",FALSE},	
		{"diracInDictionary",FALSE},
		{"gaussInDictionary",FALSE},
		{"sinCosInDictionary",FALSE},
		{"maxGaborScale",FALSE},
		{"numberOfThreads",FALSE},		
		{"progressBar",FALSE},		
	};

static STATUS checkIfLineBeginWithCommand(const char *line)
{
	unsigned short int counter;
	char tmpLine[LENGTH_OF_LINE];
	char *command;

	BOOLEAN commandExists = FALSE;

	strcpy(tmpLine,line);

	command = (char *)strtok(tmpLine," ");

	for(counter=0;counter<NUMBER_OF_ALL_COMMANDS;counter++)
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

static STATUS pushLine(Queue *lineQueue, char *text, unsigned short int number)
{
	if(strlen(text)>LENGTH_OF_LINE)
		return ERROR;  
    
	Line *line = (Line *)malloc(sizeof(Line));
    strcpy(line->text,text);

    line->number = number;
    addNode(lineQueue,(void *)line);
    
	return SUCCESS;    
}

static void popLine(Queue *lineQueue, char *text, unsigned short int *number)
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

STATUS openConfigFile(ConfigFile *configFile, char *infoMessage)
{

    configFile->file = fopen(configFile->name,"rt");
	if(configFile->file==NULL)
    {
		const char *tmpString[] = {configFile->name};
		printError(infoMessage,CAN_NOT_OPEN_CONFIG_FILE,tmpString,1);
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

STATUS readConfigFile(ConfigFile *configFile, char *infoMessage)
{
    unsigned short int lineNumber     = 0;
	char numberToString[10];
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
						sprintf(numberToString,"%hu",lineNumber);
						const char *tmpString[] = {numberToString};
						printError(infoMessage,LINE_IS_BROKEN_INCORRECTLY,tmpString,1);
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
						sprintf(numberToString,"%hu",lineNumber);
						const char *tmpString[] = {numberToString};
						printError(infoMessage,LINE_IS_BROKEN_INCORRECTLY,tmpString,1);
						status = ERROR;
						break;
					}

					positionOfBroken = (char *)strtok(text,"\\");
					if(!pushLine(configFile->lineQueue,text,lineNumber))
					{
						sprintf(numberToString,"%hu",lineNumber);
						const char *tmpString[] = {numberToString};
						printError(infoMessage,LINE_IS_TOO_LONG,tmpString,1);
						status = ERROR;
						break;
					}
					addLine = TRUE;
				}
				else
				{
					if(!pushLine(configFile->lineQueue,text,lineNumber))
					{
						sprintf(numberToString,"%hu",lineNumber);
						const char *tmpString[] = {numberToString};
						printError(infoMessage,LINE_IS_TOO_LONG,tmpString,1);
						status = ERROR;
						break;
					}
				}
			}
		}
    }while(TRUE);

    return status;

}

STATUS findDataParametersInConfigFile(ConfigFile *configFile, Dictionary *dictionary, MP5Parameters *mp5Parameters, char *infoMessage)
{

	int i;
	unsigned short int lineNumber;

	char text[LENGTH_OF_LINE];
	char numberToString[50];
	
	char *parameterPosition = NULL;

	while((configFile->lineQueue->firstNode)!=NULL)
	{
		popLine(configFile->lineQueue,text,&lineNumber);

		if(!checkIfLineBeginWithCommand(text))
		{
			sprintf(numberToString,"%hu",lineNumber);
			const char *tmpString[] = {numberToString};
			printError(infoMessage,LINE_DOES_NOT_INCLUDE_COMMAND,tmpString,1);
			return ERROR;
		}

		if(strstr(text,"nameOfDataFile")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"nameOfDataFile",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}
			
			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			commandsList[0].found = TRUE;

			strcpy(mp5Parameters->nameOfDataFile,parameterPosition);
		}
		else if(strstr(text,"extensionOfResultsFile")!=NULL)
		{
			
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"extensionOfResultsFile",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			commandsList[1].found = TRUE;

			strcpy(mp5Parameters->extensionOfResultsFile,parameterPosition);
		}
		else if(strstr(text,"nameOfOutputDirectory")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"nameOfOutputDirectory",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");
			
			if(*(parameterPosition + (strlen(parameterPosition)-1))!='/')
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"nameOfOutputDirectory",numberToString};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,2);
				return ERROR;
			}
			commandsList[2].found = TRUE;
			strcpy(mp5Parameters->nameOfOutputDirectory,parameterPosition);
		}
		else if(strstr(text,"writingMode")!=NULL)
		{
			mp5Parameters->writingMode = 0x0;

			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"writingMode",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			if(strcmp(parameterPosition,"CREATE")!=0 &&
			strcmp(parameterPosition,"APPEND")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"nameOfOutputDirectory",numberToString,"CREATE/APPEND"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);
				return ERROR;
			}

			if(strcmp(parameterPosition,"CREATE")==0)      mp5Parameters->writingMode|=CREATE_FILE;
			else if(strcmp(parameterPosition,"APPEND")==0) mp5Parameters->writingMode|=APPEND_FILE;

			commandsList[3].found = TRUE;
		}
		else if(strstr(text,"sizeOfHeader")!=NULL)
		{

			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"sizeOfHeader",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			if(!isDecimal(parameterPosition,"uint"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"sizeOfHeader",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_OR_EQUAL_TO_ZERO,tmpString,2);
				return ERROR;
			}
			commandsList[4].found = TRUE;

			mp5Parameters->sizeOfHeader = (unsigned short int)atoi(parameterPosition);
		}
		else if(strstr(text,"sizeOfTail")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"sizeOfTail",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uint"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"sizeOfTail",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_OR_EQUAL_TO_ZERO,tmpString,2);
				return ERROR;
			}

			commandsList[5].found = TRUE;

			mp5Parameters->sizeOfTail = (unsigned short int)atoi(parameterPosition);
		}
		else if(strstr(text,"samplingFrequency")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"samplingFrequency",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"samplingFrequency",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_FLOAT_GREATER_TO_ZERO,tmpString,2);
				return ERROR;
			}

			commandsList[6].found = TRUE;

			mp5Parameters->samplingFrequency = atof(parameterPosition);
		}
		else if(strstr(text,"formatOfData")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"formatOfData",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);				
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t");

			mp5Parameters->dataFormat = 0x0;

			if(strcmp(parameterPosition,"ASCII")!=0 &&
			strcmp(parameterPosition,"SHORT")!=0 &&
			strcmp(parameterPosition,"FLOAT")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"formatOfData",numberToString,"ASCII/SHORT/FLOAT"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);				
				return ERROR;
			}

			if(strcmp(parameterPosition,"ASCII")==0)      mp5Parameters->dataFormat|= FORMAT_ASCII;
			if(strcmp(parameterPosition,"SHORT")==0)      mp5Parameters->dataFormat|= FORMAT_SHORT;
			else if(strcmp(parameterPosition,"FLOAT")==0) mp5Parameters->dataFormat|= FORMAT_FLOAT;

			commandsList[7].found = TRUE;
		}
		else if(strstr(text,"numberOfChannels")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"numberOfChannels",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);				
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"numberOfChannels",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,tmpString,2);				
				return ERROR;
			}

			commandsList[8].found = TRUE;

			mp5Parameters->numberOfChannelsInDataFile = (unsigned short int)atoi(parameterPosition);
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

			if(countWords(text)==1)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"chosenChannels",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\-\n\t"))!=NULL)
			{
				if(!isDecimal(parameterPosition,"uintgz"))
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"chosenChannels",numberToString};
					printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,tmpString,2);	
					return ERROR;
				}
			}
			
			strcpy(backupText,text);

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\\n\t"))!=NULL)
				pushString(configFile->stringQueue,parameterPosition);

			while((configFile->stringQueue->firstNode)!=NULL)
			{
				popString(configFile->stringQueue,string);

				firstString  = (char *)strtok(string,"-");
				secondString = (char *)strtok(NULL,"-");

				if(!firstString)
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"chosenChannels",numberToString};
					printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);		
					return ERROR;
				}
				else if(strcmp(firstString,"-")==0)
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"chosenChannels",numberToString};
					printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);		
					return ERROR;
				}
				else
				{
					firstChannel = (unsigned short int)atoi(firstString);

					if(!secondString)
						mp5Parameters->numberOfChosenChannels++;
					else if(strcmp(secondString,"-")==0)
					{
						sprintf(numberToString,"%hu",lineNumber);
						const char *tmpString[] = {"chosenChannels",numberToString};
						printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);		
						return ERROR;
					}
					else
					{
						secondChannel = (unsigned short int)atoi(secondString);

						if(secondChannel<firstChannel)
						{
							sprintf(numberToString,"%hu",lineNumber);
							const char *tmpString[] = {"chosenChannels",numberToString};
							printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);		
							return ERROR;
						}
						else
							mp5Parameters->numberOfChosenChannels = (unsigned short int)(mp5Parameters->numberOfChosenChannels + secondChannel - firstChannel + 1);

					}
				}
			}

			mp5Parameters->chosenChannels = (unsigned short int *)fVectorAllocate(mp5Parameters->numberOfChosenChannels);

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
					mp5Parameters->chosenChannels[channelsCounter] = firstChannel;
					channelsCounter++;
				}
				else
				{
					secondChannel = (unsigned short int)atoi(secondString);

					for(channel=0;channel<secondChannel - firstChannel + 1;channel++,channelsCounter++)
						mp5Parameters->chosenChannels[channelsCounter] = (unsigned short int)(firstChannel + channel);
				}
			}

			commandsList[9].found = TRUE;
		}
		else if(strstr(text,"numberOfPointsInOffset")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"numberOfPointsInOffset",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"numberOfPointsInOffset",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,tmpString,2);	
				return ERROR;
			}

			mp5Parameters->offsetDimension = (unsigned int)atoi(parameterPosition);;
			
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
			
			if(countWords(text)==1)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"chosenOffsets",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\-\n\t"))!=NULL)
				if(!isDecimal(parameterPosition,"uintgz"))
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"chosenOffsets",numberToString};
					printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,tmpString,2);	
					return ERROR;
				}

			strcpy(backupText,text);

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\\n\t"))!=NULL)
				pushString(configFile->stringQueue,parameterPosition);

			while((configFile->stringQueue->firstNode)!=NULL)
			{
				popString(configFile->stringQueue,string);

				firstString  = (char *)strtok(string,"-");
				secondString = (char *)strtok(NULL,"-");

				if(!firstString)
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"chosenOffsets",numberToString};
					printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);		
					return ERROR;
				}
				else if(strcmp(firstString,"-")==0)
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"chosenOffsets",numberToString};
					printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);		
					return ERROR;			
				}
				else
				{
					firstOffset = (unsigned short int)atoi(firstString);

					if(!secondString)
						mp5Parameters->numberOfChosenOffsets++;
					else if(strcmp(secondString,"-")==0)
					{
						sprintf(numberToString,"%hu",lineNumber);
						const char *tmpString[] = {"chosenOffsets",numberToString};
						printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);		
						return ERROR;
					}
					else
					{
						secondOffset = (unsigned short int)atoi(secondString);

						if(secondOffset<firstOffset)
						{
							sprintf(numberToString,"%hu",lineNumber);
							const char *tmpString[] = {"chosenOffsets",numberToString};
							printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);		
							return ERROR;
						}
						else
							mp5Parameters->numberOfChosenOffsets = (unsigned short int)(mp5Parameters->numberOfChosenOffsets + secondOffset - firstOffset + 1);
					}
				}
			}

			mp5Parameters->chosenOffsets = (unsigned short int *)iVectorAllocate(mp5Parameters->numberOfChosenOffsets);

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
					mp5Parameters->chosenOffsets[offsetsCounter] = firstOffset;
					offsetsCounter++;
				}
				else
				{
					secondOffset = (unsigned short int)atoi(secondString);

					for(offset=0;offset<secondOffset - firstOffset + 1;offset++,offsetsCounter++)
						mp5Parameters->chosenOffsets[offsetsCounter] = (unsigned short int)(firstOffset + offset);
				}
			}

			commandsList[11].found = TRUE;
		}
		else if(strstr(text,"typeOfDictionary")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"typeOfDictionary",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);		
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"OCTAVE_FIXED")!=0 && strcmp(parameterPosition,"OCTAVE_STOCH")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"typeOfDictionary",numberToString,"OCTAVE_FIXED/OCTAVE_STOCH"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);		
				return ERROR;
			}

			if(strcmp(parameterPosition,"OCTAVE_FIXED")==0)
				dictionary->typeOfDictionary|= OCTAVE_FIXED;
			else if(strcmp(parameterPosition,"OCTAVE_STOCH")==0)
				dictionary->typeOfDictionary|= OCTAVE_STOCH;

			commandsList[12].found = TRUE;
		}
		else if(strstr(text,"dilationFactor")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"typeOfDictionary",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);		
				return ERROR;
			}
			
			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"dilationFactor",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_FLOAT_GREATER_TO_ONE,tmpString,2);		
				return ERROR;
			}

			dictionary->dilationFactor = (double)atof(parameterPosition);
			commandsList[13].found = TRUE;
		}
		else if(strstr(text,"randomSeed")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"randomSeed",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);		
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if((strcmp(parameterPosition,"auto")!=0) && (!isDecimal(parameterPosition,"uint")))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"randomSeed",numberToString,"or string = auto"};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,tmpString,3);		
				return ERROR;
			}

			if(strcmp(parameterPosition,"auto")==0)
				dictionary->randomSeed = AUTO_RANDOM_SEED;
			else
				dictionary->randomSeed = (long int)atoi(parameterPosition);
				
			commandsList[14].found = TRUE;
		}
		else if(strstr(text,"periodDensity")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"periodDensity",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);		
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"periodDensity",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,tmpString,2);		
				return ERROR;
			}

			dictionary->periodDensity = (unsigned short int)atoi(parameterPosition);
			commandsList[15].found = TRUE;
		}
		else if(strstr(text,"reinitDictionary")!=NULL)
		{
			mp5Parameters->reinitDictionary = 0x0;

			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"reinitDictionary",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);		
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"NO_REINIT_AT_ALL")!=0 &&
			strcmp(parameterPosition,"REINIT_IN_CHANNEL_DOMAIN")!=0 &&
			strcmp(parameterPosition,"REINIT_IN_OFFSET_DOMAIN")!=0  &&
			strcmp(parameterPosition,"REINIT_AT_ALL")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"reinitDictionary",numberToString,"NO_REINIT_AT_ALL/REINIT_IN_CHANNEL_DOMAIN/REINIT_IN_OFFSET_DOMAIN/REINIT_AT_ALL"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);		
				return ERROR;
			}

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

			commandsList[16].found = TRUE;
		}
		else if(strstr(text,"scaleToPeriodFactor")!=NULL)
		{
			if(countWords(text)!=2)
			{	
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"scaleToPeriodFactor",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"real"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"scaleToPeriodFactor",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_FLOAT_GREATER_OR_EQUAL_TO_ZERO,tmpString,2);	
				return ERROR;
			}

			dictionary->scaleToPeriodFactor = (double)atof(parameterPosition);
			commandsList[17].found = TRUE;
		}
		else if(strstr(text,"maximalNumberOfIterations")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"maximalNumberOfIterations",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"maximalNumberOfIterations",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_OR_EQUAL_TO_ZERO,tmpString,2);	
				return ERROR;
			}

			mp5Parameters->maximalNumberOfIterations = (unsigned short int)atoi(parameterPosition);

			commandsList[18].found = TRUE;
		}
		else if(strstr(text,"energyPercent")!=NULL)
		{
			if(countWords(text)!=2)			
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"energyPercent",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}
			
			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"energyPercent",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_FLOAT_GREATER_OR_EQUAL_TO_ZERO,tmpString,2);	
				return ERROR;
			}

			mp5Parameters->energyPercent = atof(parameterPosition);

			commandsList[19].found = TRUE;
		}
		else if(strstr(text,"MP")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"MP",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}
			
			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"SMP")!=0  &&
			strcmp(parameterPosition,"MMP1")!=0 &&
			strcmp(parameterPosition,"MMP2")!=0 &&
			strcmp(parameterPosition,"MMP3")!=0 &&
			strcmp(parameterPosition,"MMP4")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"MP",numberToString,"SMP/MMP1/MMP2/MMP3"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);	
				return ERROR;
			}
			
			if(strcmp(parameterPosition,"SMP")==0)       mp5Parameters->MPType|= SMP;
			else if(strcmp(parameterPosition,"MMP1")==0) mp5Parameters->MPType|= MMP1;
			else if(strcmp(parameterPosition,"MMP2")==0) mp5Parameters->MPType|= MMP2;
			else if(strcmp(parameterPosition,"MMP3")==0) mp5Parameters->MPType|= MMP3;
			else if(strcmp(parameterPosition,"MMP4")==0) mp5Parameters->MPType|= MMP4;

			commandsList[20].found = TRUE;
		}
		else if(strstr(text,"analiticalDotProduct")!=NULL)
		{
			if(countWords(text)!=2)	
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"analiticalDotProduct",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"ON")!=0 && strcmp(parameterPosition,"OFF")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"analiticalDotProduct",numberToString,"ON/OFF"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);	
				return ERROR;
			}

			if(strcmp(parameterPosition,"ON")==0)       mp5Parameters->analiticalDotProduct = OFF;
			else if(strcmp(parameterPosition,"OFF")==0) mp5Parameters->analiticalDotProduct = ON;

			commandsList[21].found = TRUE;
		}
		else if(strstr(text,"bookWithSignal")!=NULL)
		{
			if(countWords(text)!=2)	
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"bookWithSignal",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"YES")!=0 && strcmp(parameterPosition,"NO")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"bookWithSignal",numberToString,"YES/NO"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);	
				return ERROR;
			}
			
			if(strcmp(parameterPosition,"NO")==0)       mp5Parameters->bookWithSignal = NO;
			else if(strcmp(parameterPosition,"YES")==0) mp5Parameters->bookWithSignal = YES;

			commandsList[22].found = TRUE;
		}
		else if(strstr(text,"pointsPerMicrovolt")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"pointsPerMicrovolt",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"pointsPerMicrovolt",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_FLOAT_GREATER_TO_ZERO,tmpString,2);	
				return ERROR;
			}

			mp5Parameters->pointsPerMicrovolt = atof(parameterPosition);

			commandsList[23].found = TRUE;
		}
		else if(strstr(text,"FFT")!=NULL)
		{		
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"FFT",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}
			
			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");
			
			if(strcmp(parameterPosition,"OFF")!=0 && strcmp(parameterPosition,"FFT1")!=0 && strcmp(parameterPosition,"FFT2")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"FFT",numberToString,"OFF/FFT1/FFT2"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);
				return ERROR;
			}

			if(strcmp(parameterPosition,"OFF")==0)          mp5Parameters->FFT = OFF;
			else if(strcmp(parameterPosition,"FFT1")==0)   mp5Parameters->FFT = FFT1;
			else if(strcmp(parameterPosition,"FFT2")==0)   mp5Parameters->FFT = FFT2;
		}
		else if(strstr(text,"accuracy")!=NULL)
		{	
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"accuracy",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}
			
			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");
			
			if(strcmp(parameterPosition,"FLOATING")!=0 && strcmp(parameterPosition,"FULL")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"accuracy",numberToString,"FLOATING/FULL"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);
				return ERROR;
			}

			if(strcmp(parameterPosition,"FLOATING")==0)    mp5Parameters->accuracy = FLOATING;
			else if(strcmp(parameterPosition,"FULL")==0)   mp5Parameters->accuracy = FULL;
		}
		else if(strstr(text,"diracInDictionary")!=NULL)
		{
			if(countWords(text)!=2)	
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"diracInDictionary",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"YES")!=0 && strcmp(parameterPosition,"NO")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"diracInDictionary",numberToString,"YES/NO"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);	
				return ERROR;
			}
			
			if(strcmp(parameterPosition,"NO")==0)       dictionary->diracInDictionary = NO;
			else if(strcmp(parameterPosition,"YES")==0) dictionary->diracInDictionary = YES;
		}
		else if(strstr(text,"gaussInDictionary")!=NULL)
		{
			if(countWords(text)!=2)	
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"gaussInDictionary",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"YES")!=0 && strcmp(parameterPosition,"NO")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"gaussInDictionary",numberToString,"YES/NO"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);	
				return ERROR;
			}
			
			if(strcmp(parameterPosition,"NO")==0)       dictionary->gaussInDictionary = NO;
			else if(strcmp(parameterPosition,"YES")==0) dictionary->gaussInDictionary = YES;
		}
		else if(strstr(text,"sinCosInDictionary")!=NULL)
		{
			if(countWords(text)!=2)	
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"sinCosInDictionary",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"YES")!=0 && strcmp(parameterPosition,"NO")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"sinCosInDictionary",numberToString,"YES/NO"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);	
				return ERROR;
			}
			
			if(strcmp(parameterPosition,"NO")==0)       dictionary->sinCosInDictionary = NO;
			else if(strcmp(parameterPosition,"YES")==0) dictionary->sinCosInDictionary = YES;
		}
		else if(strstr(text,"maxGaborScale")!=NULL)
		{		
				if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"maxGaborScale",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uint"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"maxGaborScale",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_OR_EQUAL_TO_ZERO,tmpString,2);	
				return ERROR;
			}

			mp5Parameters->maxGaborScale = atof(parameterPosition);
		}
		else if(strstr(text,"numberOfThreads")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"numberOfThreads",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,3);	
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"numberOfThreads",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,tmpString,2);	
				return ERROR;
			}

			mp5Parameters->numberOfThreads = (unsigned char)atoi(parameterPosition);
		}
		else if(strstr(text,"progressBar")!=NULL)
		{		
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"progressBar",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);	
				return ERROR;
			}
			
			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");
			
			if(strcmp(parameterPosition,"ON")!=0 && strcmp(parameterPosition,"OFF")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"progressBar",numberToString,"ON/OFF"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);
				return ERROR;
			}

			if(strcmp(parameterPosition,"ON")==0)       mp5Parameters->progressBar = ON;
			else if(strcmp(parameterPosition,"OFF")==0) mp5Parameters->progressBar = OFF;
		}
	}

	for(i=0;i<NUMBER_OF_PERMAMENT_COMMANDS;i++)
	{
		if(commandsList[i].found == FALSE)
		{
			const char *tmpString[] = {commandsList[i].command};
			printError(infoMessage,COMMAND_NOT_FOUND,tmpString,1);	
			return ERROR;
		}
	}

	return SUCCESS;

}
