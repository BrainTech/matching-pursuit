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
		{"nameOfOutputDirectory",FALSE},
		{"writingMode",FALSE},
		{"samplingFrequency",FALSE},
		{"numberOfChannels",FALSE},
		{"selectedChannels",FALSE},
		{"numberOfSamplesInEpoch",FALSE},
		{"selectedEpochs",FALSE},
		{"typeOfDictionary",FALSE},
		{"energyError",FALSE},
		{"randomSeed",FALSE},
		{"reinitDictionary",FALSE},
		{"maximalNumberOfIterations",FALSE},
		{"energyPercent",FALSE},
		{"MP",FALSE},
		{"scaleToPeriodFactor",FALSE}, 	
     	{"normType",FALSE},
     	{"pointsPerMicrovolt",FALSE},
		{"diracInDictionary",FALSE},
		{"gaussInDictionary",FALSE},
		{"sinCosInDictionary",FALSE},
		{"gaborInDictionary",TRUE},
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
				const char *tmpString[] = {"nameOfOutputDirectory",numberToString,"/"};
				printError(infoMessage,BAD_NAME_OF_OUTPUT_DIRECTORY,tmpString,3);
				return ERROR;
			}
			commandsList[1].found = TRUE;
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

			commandsList[2].found = TRUE;
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

			commandsList[3].found = TRUE;

			mp5Parameters->samplingFrequency = atof(parameterPosition);
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

			commandsList[4].found = TRUE;

			mp5Parameters->numberOfChannelsInDataFile = (unsigned short int)atoi(parameterPosition);
		}
		else if(strstr(text,"selectedChannels")!=NULL)
		{
			char backupText[LENGTH_OF_LINE];

			char string[LENGTH_OF_STRING];
			char *firstString;
			char *secondString;

			unsigned int channel, channelsCounter = 0;
			unsigned int firstChannel;
			unsigned int secondChannel;

			strcpy(backupText,text);

			if(countWords(text)==1)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"selectedChannels",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\-\n\t"))!=NULL)
			{
				if(!isDecimal(parameterPosition,"uintgz"))
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"selectedChannels",numberToString};
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
					const char *tmpString[] = {"selectedChannels",numberToString};
					printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);
					return ERROR;
				}
				else if(strcmp(firstString,"-")==0)
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"selectedChannels",numberToString};
					printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);
					return ERROR;
				}
				else
				{
					firstChannel = (unsigned short int)atoi(firstString);

					if(!secondString)
						mp5Parameters->numberOfSelectedChannels++;
					else if(strcmp(secondString,"-")==0)
					{
						sprintf(numberToString,"%hu",lineNumber);
						const char *tmpString[] = {"selectedChannels",numberToString};
						printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);
						return ERROR;
					}
					else
					{
						secondChannel = (unsigned short int)atoi(secondString);

						if(secondChannel<firstChannel)
						{
							sprintf(numberToString,"%hu",lineNumber);
							const char *tmpString[] = {"selectedChannels",numberToString};
							printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);
							return ERROR;
						}
						else
							mp5Parameters->numberOfSelectedChannels = (unsigned short int)(mp5Parameters->numberOfSelectedChannels + secondChannel - firstChannel + 1);

					}
				}
			}

			mp5Parameters->selectedChannels = (unsigned short int *)fVectorAllocate(mp5Parameters->numberOfSelectedChannels);

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
					mp5Parameters->selectedChannels[channelsCounter] = firstChannel;
					channelsCounter++;
				}
				else
				{
					secondChannel = (unsigned short int)atoi(secondString);

					for(channel=0;channel<secondChannel - firstChannel + 1;channel++,channelsCounter++)
						mp5Parameters->selectedChannels[channelsCounter] = (unsigned short int)(firstChannel + channel);
				}
			}

			commandsList[5].found = TRUE;
		}
		else if(strstr(text,"numberOfSamplesInEpoch")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"numberOfSamplesInEpoch",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isDecimal(parameterPosition,"uintgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"numberOfSamplesInEpoch",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,tmpString,2);
				return ERROR;
			}

			mp5Parameters->epochSize = (unsigned int)atoi(parameterPosition);

			commandsList[6].found = TRUE;
		}
		else if(strstr(text,"selectedEpochs")!=NULL)
		{
			char backupText[LENGTH_OF_LINE];

			char string[LENGTH_OF_STRING];
			char *firstString;
			char *secondString;

			unsigned short int epoch, epochsCounter = 0;
			unsigned short int firstEpoch;
			unsigned short int secondEpoch;

			strcpy(backupText,text);

			if(countWords(text)==1)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"selectedEpochs",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\-\n\t"))!=NULL)
				if(!isDecimal(parameterPosition,"uintgz"))
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"selectedEpochs",numberToString};
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
					const char *tmpString[] = {"selectedEpochs",numberToString};
					printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);
					return ERROR;
				}
				else if(strcmp(firstString,"-")==0)
				{
					sprintf(numberToString,"%hu",lineNumber);
					const char *tmpString[] = {"selectedEpochs",numberToString};
					printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);
					return ERROR;
				}
				else
				{
					firstEpoch = (unsigned short int)atoi(firstString);

					if(!secondString)
						mp5Parameters->numberOfSelectedEpochs++;
					else if(strcmp(secondString,"-")==0)
					{
						sprintf(numberToString,"%hu",lineNumber);
						const char *tmpString[] = {"selectedEpochs",numberToString};
						printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);
						return ERROR;
					}
					else
					{
						secondEpoch = (unsigned short int)atoi(secondString);

						if(secondEpoch<firstEpoch)
						{
							sprintf(numberToString,"%hu",lineNumber);
							const char *tmpString[] = {"selectedEpochs",numberToString};
							printError(infoMessage,INCORRECT_SYNTAX_OF_ARGUMENT,tmpString,2);
							return ERROR;
						}
						else
							mp5Parameters->numberOfSelectedEpochs = (unsigned short int)(mp5Parameters->numberOfSelectedEpochs + secondEpoch - firstEpoch + 1);
					}
				}
			}

			mp5Parameters->selectedEpochs = (unsigned short int *)iVectorAllocate(mp5Parameters->numberOfSelectedEpochs);

			strcpy(backupText,text);

			parameterPosition = (char *)strtok(backupText," \n\t");

			while((parameterPosition = (char *)strtok(NULL," \\\n\t"))!=NULL)
				pushString(configFile->stringQueue,parameterPosition);

			while((configFile->stringQueue->firstNode)!=NULL)
			{
				popString(configFile->stringQueue,string);

				firstString  = (char *)strtok(string,"-");
				secondString = (char *)strtok(NULL,"-");

				firstEpoch = (unsigned short int)atoi(firstString);

				if(!secondString)
				{
					mp5Parameters->selectedEpochs[epochsCounter] = firstEpoch;
					epochsCounter++;
				}
				else
				{
					secondEpoch = (unsigned short int)atoi(secondString);

					for(epoch=0;epoch<secondEpoch - firstEpoch + 1;epoch++,epochsCounter++)
						mp5Parameters->selectedEpochs[epochsCounter] = (unsigned short int)(firstEpoch + epoch);
				}
			}

			commandsList[7].found = TRUE;
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
			{
				dictionary->typeOfDictionary|= OCTAVE_STOCH;
				mp5Parameters->FFT = OFF;
         	}

			commandsList[8].found = TRUE;
		}
		else if(strstr(text,"energyError")!=NULL)
		{
			if(countWords(text)!=3)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"energyError",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"real"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"energyError",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_FLOAT,tmpString,2);
				return ERROR;
			}

			if(((double)atof(parameterPosition) <= 0.0) || (double)atof(parameterPosition) >= 1.0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"energyError",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_FLOAT_BETWEEN_ZERO_AND_ONE,tmpString,2);
				return ERROR;
			}

			dictionary->energyError = sqrt((double)atof(parameterPosition));
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(!isReal(parameterPosition,"realgz"))
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"stochasticDictionaryReductionCoefficient",numberToString};
				printError(infoMessage,ARGUMENT_SHOUDL_BE_FLOAT_GREATER_TO_ONE,tmpString,2);
				return ERROR;
			}

			dictionary->stochasticDictionaryReductionCoefficient = (double)atof(parameterPosition)/100.0;
			commandsList[9].found = TRUE;
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

			commandsList[10].found = TRUE;
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
			strcmp(parameterPosition,"REINIT_IN_EPOCH_DOMAIN")!=0  &&
			strcmp(parameterPosition,"REINIT_AT_ALL")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"reinitDictionary",numberToString,"NO_REINIT_AT_ALL/REINIT_IN_CHANNEL_DOMAIN/REINIT_IN_EPOCH_DOMAIN/REINIT_AT_ALL"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);
				return ERROR;
			}

			if(strcmp(parameterPosition,"NO_REINIT_AT_ALL")==0)
				mp5Parameters->reinitDictionary|= NO_REINIT_AT_ALL;
			else if(strcmp(parameterPosition,"REINIT_IN_CHANNEL_DOMAIN")==0)
				mp5Parameters->reinitDictionary|= REINIT_IN_CHANNEL_DOMAIN;
			else if(strcmp(parameterPosition,"REINIT_IN_CHANNEL_DOMAIN")==0)
				mp5Parameters->reinitDictionary|= REINIT_IN_CHANNEL_DOMAIN;
			else if(strcmp(parameterPosition,"REINIT_IN_EPOCH_DOMAIN")==0)
				mp5Parameters->reinitDictionary|= REINIT_IN_EPOCH_DOMAIN;
			else if(strcmp(parameterPosition,"REINIT_AT_ALL")==0)
				mp5Parameters->reinitDictionary|= REINIT_AT_ALL;

			commandsList[11].found = TRUE;
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

			commandsList[12].found = TRUE;
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

			commandsList[13].found = TRUE;
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
			strcmp(parameterPosition,"MMP11")!=0 &&
			strcmp(parameterPosition,"MMP12")!=0 &&
			strcmp(parameterPosition,"MMP2")!=0 &&
			strcmp(parameterPosition,"MMP21")!=0 &&
			strcmp(parameterPosition,"MMP22")!=0 &&
			strcmp(parameterPosition,"MMP23")!=0 &&
			strcmp(parameterPosition,"MMP3")!=0 &&
			strcmp(parameterPosition,"MMP32")!=0 &&
			strcmp(parameterPosition,"MMP33")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"MP",numberToString,"SMP/MMP1/MMP11/MMP12/MMP21/MMP2/MMP22/MMP3/MMP23/MMP32/MMP33"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);
				return ERROR;
			}

			if(strcmp(parameterPosition,"SMP")==0)        mp5Parameters->MPType= SMP;
			else if(strcmp(parameterPosition,"MMP1")==0)  mp5Parameters->MPType= MMP1;
			else if(strcmp(parameterPosition,"MMP11")==0) mp5Parameters->MPType= MMP11;
			else if(strcmp(parameterPosition,"MMP12")==0) mp5Parameters->MPType= MMP12;
			else if(strcmp(parameterPosition,"MMP21")==0) mp5Parameters->MPType= MMP21;
			else if(strcmp(parameterPosition,"MMP2")==0)  mp5Parameters->MPType= MMP2;
			else if(strcmp(parameterPosition,"MMP22")==0) mp5Parameters->MPType= MMP22;
			else if(strcmp(parameterPosition,"MMP23")==0) mp5Parameters->MPType= MMP23;
			else if(strcmp(parameterPosition,"MMP3")==0)  mp5Parameters->MPType= MMP3;
			else if(strcmp(parameterPosition,"MMP32")==0) mp5Parameters->MPType= MMP32;
			else if(strcmp(parameterPosition,"MMP33")==0) mp5Parameters->MPType= MMP33;

			commandsList[14].found = TRUE;
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
			commandsList[15].found = TRUE;
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

			mp5Parameters->pointsPerMicrovolt = (double)atof(parameterPosition);
			commandsList[16].found = TRUE;
		}
		else if(strstr(text,"normType")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"normType",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"L1")!=0 && strcmp(parameterPosition,"L2")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"normType",numberToString,"L1/L2"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);
				return ERROR;
			}

			if(strcmp(parameterPosition,"L1")==0)       mp5Parameters->normType = L1;
			else if(strcmp(parameterPosition,"YES")==0) mp5Parameters->normType = L2;
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
		else if(strstr(text,"gaborInDictionary")!=NULL)
		{
			if(countWords(text)!=2)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"gaborInDictionary",numberToString};
				printError(infoMessage,INCORRECT_NUMBER_OF_ARGUMENTS,tmpString,2);
				return ERROR;
			}

			parameterPosition = (char *)strtok(text," \n\t");
			parameterPosition = (char *)strtok(NULL," \n\t\0");

			if(strcmp(parameterPosition,"YES")!=0 && strcmp(parameterPosition,"NO")!=0)
			{
				sprintf(numberToString,"%hu",lineNumber);
				const char *tmpString[] = {"gaborInDictionary",numberToString,"YES/NO"};
				printError(infoMessage,INCORRECT_TYPE_OF_ARGUMENT,tmpString,3);
				return ERROR;
			}

			if(strcmp(parameterPosition,"NO")==0)       dictionary->gaborInDictionary = NO;
			else if(strcmp(parameterPosition,"YES")==0) dictionary->gaborInDictionary = YES;
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
