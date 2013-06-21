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

#define _GNU_SOURCE

#include<dirent.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<strings.h>
#include<unistd.h>
#include"atom.h"
#include"def.h"
#include"io_mp5.h"
#include"matrix.h"
#include"mp5.h"
#include"queue.h"
#include"types.h"
#include"vector.h"
#include"tools.h"

extern unsigned char applicationMode;

static char aTmp[50], bTmp[50];

struct ErrorCodeAndText
{
	unsigned short int errorNumber;
	char               *errorCode;
	char               *errorText;
}errorCodeAndText[] =
{
{CAN_NOT_OPEN_CONFIG_FILE,                           "error.mp5executable.badConfigFile",                             "Failed to open config file {0}"},
{LINE_IS_BROKEN_INCORRECTLY,                         "error.mp5executable.lineIsBrokenIncorrectly",                   "The line #{0} in config file is broken incorrectly"},
{LINE_IS_TOO_LONG,                                   "error.mp5executable.lineIsTooLong",                             "The line #{0} in config file is too long"},
{COMMAND_NOT_FOUND,                                  "error.mp5executable.lackOfCommand",                             "The command {0} has not been found in config file"},
{INCORRECT_NUMBER_OF_ARGUMENTS,                      "error.mp5executable.incorrectNumberOfArgments",                 "The command {0} in line #{1} of config file includes incorrect number of arguments"},
{INCORRECT_TYPE_OF_ARGUMENT,                         "error.mp5executable.incorrectTypeOfArgument",                   "The command {0} in line #{1} of config file has been call with incorrect argments. Allowed arguments are: {3}"},
{INCORRECT_SYNTAX_OF_ARGUMENT,                       "error.mp5executable.incorrectSyntaxOfArgument",                 "The order of arguments of command {0} in line #{1} is incorrect"},
{LINE_DOES_NOT_INCLUDE_COMMAND,                      "error.mp5executable.noCommandInLine",                           "The line #{0} in config file does not contain the command"},
{ARGUMENT_SHOUDL_BE_INTEGER_GREATER_OR_EQUAL_TO_ZERO,"error.mp5executable.argumentShouldBeIntegerGraterOrEqualToZero","The argument {0} in line #{1} of config file should be integer>=0"},
{ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO,         "error.mp5executable.argumentShouldBeIntegerGraterToZero",       "The argument {0} in line #{1} of config file shuold be integer>0"},
{ARGUMENT_SHOUDL_BE_FLOAT_GREATER_OR_EQUAL_TO_ZERO,  "error.mp5executable.argumentShouldBeFloatGraterOrEqualToZero",  "The argument {0} in line #{1} of config file should be float>=0.0"},
{ARGUMENT_SHOUDL_BE_FLOAT,                           "error.mp5executable.argumentShouldBeFloat",                     "The argument {0} in line #{1} of config file should be float"},
{ARGUMENT_SHOUDL_BE_FLOAT_GREATER_TO_ZERO,           "error.mp5executable.argumentShouldBeFloatGraterZero",           "The argument {0} in line #{1} of config file should be float>0.0"},
{ARGUMENT_SHOUDL_BE_FLOAT_GREATER_TO_ONE,            "error.mp5executable.argumentShouldBeFloatGraterToOne",          "The argument {0} in line #{1} of config file should be float>1.0"},
{ARGUMENT_SHOUDL_BE_FLOAT_BETWEEN_ZERO_AND_ONE,      "error.mp5executable.argumentShouldBeBetweenZeroAndOne",         "The argument {0} in line #{1} of config file should be float in range (0.0 1.0)"},
{CAN_NOT_OPEN_DIRECTORY,                             "error.mp5executable.badConfigFile",                             "Failed to open config file {0}"},
{OVERWRITE_RESULTS_ALARM,                            "error.mp5executable.overwriteResultsAlarm",                     "You are trying to overwrite file {0} with results"},
{BAD_NUMBER_OF_SELECTED_CHANNELS,                    "error.mp5executable.badNumberOfSelectedChannels",               "The number of declared selected channel - {0} is greater then number of channels - {1} in file {2}"},
{BAD_NUMBER_OF_SELECTED_EPOCHS,                      "error.mp5executable.badNumberSelectedEpochs",                   "The number of declared selected epochs - {0} is greater then number of epochs - {1} in file {2}"},
{BAD_ENERGY_PERCENT,                                 "error.mp5executable.badEnergyPercent",                          "The value of argument energyPercent can not be greater or equal to 100 percent"},
{BAD_REINIT_ALL,                                     "error.mp5executable.badReinitALL",                              "There is no sense to reinit dictionary with the same seed for eac epoch/channel"},
{BAD_REINIT_MMP,                                     "error.mp5executable.badReinitMMP",                              "There is no sense to use Multichannel MP Algorithms and reinit dictionary in channel domain"},
{BAD_MP_ALGORITHM,                                   "error.mp5executable.badMPAlgorithm",                            "There is no sense to process Multichannel MP Algorithm for one channel only"},
{CAN_NOT_OPEN_DATA_FILE,                             "error.mp5executable.canNotOpenDataFile",                        "Failed to open file {0} with data"},
{CAN_NOT_OPEN_RESULTS_FILE,                          "error.mp5executable.canNotOpenResultsFile",                     "Failed to open file {0} where results will be writting"},
{BAD_NUMBER_OF_SAMPLES_PER_CHANNEL,                  "error.mp5executable.badNumberOfSamplesPerChannel",              "Number of samples per channels is not the same for each channel in data file {0}"},
{BAD_NUMBER_OF_SAMPLES,                              "error.mp5executable.badNumberOfSamples",                        "Incorrect number of samples in data file {0}"},
{CAN_NOT_READ_SAMPLE,                                "error.mp5executable.canNotReadSample",                          "The sample in ascii {0} file in line {1}, channel {2} can not be read"},
{CAN_NOT_WRITE_HEADER,                               "error.mp5executable.canNotWriteHeader",                         "Failed to write file header to results file {0}"},
{CAN_NOT_WRITE_RESULTS,                              "error.mp5executable.canNotWriteResults",                        "Failed to write results to file {0}"},
{CAN_NOT_READ_EPOCH_IN_BINARY_FILE,                  "error.mp5executable.canNotReadEpochInBinaryFile",               "Failed to read sample {0} in epoch {1} of data file {2}"},
{CAN_NOT_READ_EPOCH_IN_ASCII_FILE,                   "error.mp5executable.canNotReadEpochInAsciiFile",                "Failed to read epoch {0}, line {1}, channel {2}, in ascii data file {3}"},
{BAD_NAME_OF_OUTPUT_DIRECTORY,                       "error.mp5executable.badNameOfOutputDirectory",                  "The argument of command {0}, line {1}, should finish with {2}"},
};

STATUS generateErrotTextAndCodeGeneratorFile()
{
	unsigned int counter;
	FILE *errorFile = fopen("errorTranslation.xml","wt");

	if(!errorFile)
	{
		fprintf(stderr,"Can't open the file: \"errorTranslation.xml\" to write translations of errors\n");
		return ERROR;
	}
	else
	{
		fprintf(errorFile,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
		fprintf(errorFile,"<!DOCTYPE properties SYSTEM \"http:\\\\java.sun.com\\dtd\\properties.dtd\">\n");
		fprintf(errorFile,"<properties>\n");

		for(counter=0;counter<NUMBER_OF_ERRORS;counter++)
			fprintf(errorFile,"\t<entry key=\"%s\">%s</entry>\n",errorCodeAndText[counter].errorCode,errorCodeAndText[counter].errorText);

		fprintf(errorFile,"</properties>\n");
		fclose(errorFile);
		return SUCCESS;
	}
}

static void makeErrorText(char *infoMessage, unsigned short int errorNumber, const char *acceptableValues[], unsigned short int sizeOfAcceptableValues)
{
	char tmpInfoMessage[LENGTH_OF_INFO_MESSAGE];
	char *pointer;
	unsigned short int counter = 0;

	bzero(tmpInfoMessage,LENGTH_OF_INFO_MESSAGE);

	for(counter=0;counter<NUMBER_OF_ERRORS;counter++)
	{
		if(errorCodeAndText[counter].errorNumber == errorNumber)
		{
			strcpy(tmpInfoMessage,errorCodeAndText[counter].errorText);
			break;
		}
	}

	pointer = strtok(tmpInfoMessage,strtokMask);
	counter = 0;

	if(pointer!=NULL)
	{
		strcat(infoMessage,pointer);

		if(sizeOfAcceptableValues>0)
			strcat(infoMessage,acceptableValues[counter++]);
	}

	do
	{
		pointer = strtok(NULL,strtokMask);

		if(pointer!=NULL)
		{
			strcat(infoMessage,pointer);
			if(counter<sizeOfAcceptableValues)
				strcat(infoMessage,acceptableValues[counter++]);
		}
	}while(pointer!=NULL);
}

void printError(char *infoMessage, unsigned short int errorNumber, const char *acceptableValues[], unsigned short int sizeOfAcceptableValues)
{
	unsigned short int counter;
	strcat(infoMessage,"ERROR ");

	if((applicationMode & PROCESS_USER_MODE) || (applicationMode & TEST_PARAMETERS_MODE))
		makeErrorText(infoMessage,errorNumber,acceptableValues,sizeOfAcceptableValues);
	else if(applicationMode & PROCESS_SERVER_MODE)
	{
		for(counter=0;counter<NUMBER_OF_ERRORS;counter++)
		{
			if(errorCodeAndText[counter].errorNumber == errorNumber)
			{
				strcat(infoMessage,errorCodeAndText[counter].errorCode);
				break;
			}
		}

		if(acceptableValues!=NULL)
		{
			for(counter=0;counter<sizeOfAcceptableValues;counter++)
			{
				strcat(infoMessage," ");
				strcat(infoMessage,acceptableValues[counter]);
			}
		}
	}

}

#define MAGIC_TEXT         "MPv5.0"
#define WEB_SITE_LINK_TEXT "http:\\\\brain.fuw.edu.pl\\mp"

#define NUMBER_OF_FIELDS_IN_FILE_HEADER_SEGMENT 256   /* number of fields in file header segment, which are other then
													          text fields (for example DateFiled) ot fileds, which require
													          dynamical allocation (are coinsisted of pointers)  */
#define NUMBER_OF_CHARS_IN_SIGNATURE 			10

#define COMMENT_SEGMENT_IDENTITY     ((unsigned char)1)

#define FILE_HEADER_SEGMENT_IDENTITY ((unsigned char)2)
#define WEB_SITE_LINK_FIELD_IDENTITY ((unsigned char)3)
#define DATE_FIELD_IDENTITY 		 ((unsigned char)4)
#define SIGNAL_FIELD_IDENTITY        ((unsigned char)5)
#define DECOMPOSING_FIELD_IDENTITY   ((unsigned char)6)

#define EPOCH_SEGMENT_IDENTITY   	 ((unsigned char)7)
#define SIGNAL_SEGMENT_IDENTITY   	 ((unsigned char)8)
#define ATOMS_SEGMENT_IDENTITY   	 ((unsigned char)9)

#define DIRACDELTA_IDENTITY		     ((unsigned char)10)
#define GAUSSFUNCTION_IDENTITY 	 	 ((unsigned char)11)
#define SINCOSWAVE_IDENTITY    	 	 ((unsigned char)12)
#define GABORWAVE_IDENTITY     	 	 ((unsigned char)13)


#define FIELD_DESCRIPTOR_SIGNATURE "cc"

typedef struct
{
	unsigned char codeOfField;
	unsigned char sizeOfFieldData;
} __attribute__((packed)) FieldDescriptor;

#define SEGMENT_DESCRIPTOR_SIGNATURE "ci"

typedef struct
{
	unsigned char codeOfSegment;
	unsigned int  sizeOfSegmentData;
} __attribute__((packed)) SegmentDescriptor;

typedef struct
{
	FieldDescriptor fieldDescriptor;
	char *webSiteLink;
} __attribute__((packed)) WebSiteLinkField;

typedef struct
{
	FieldDescriptor fieldDescriptor;
	char *date;
} __attribute__((packed)) DateField;

#define SIGNAL_FIELD_SIGNATURE "ffh"

typedef struct
{
	FieldDescriptor    fieldDescriptor;
	float 		  	   samplingFrequency;
	float 			   pointsPerMicrovolt;
	unsigned short int numberOfChannelsInDataFile;
} __attribute__((packed)) SignalField;

#define DECOMPOSING_FIELD_SIGNATURE "fiic"

typedef struct
{
	FieldDescriptor    fieldDescriptor;
	float              energyPercent;
	unsigned int   	   maximalNumberOfIterations;
	unsigned int   	   sizeOfDictionary;
	char           	   typeOfDictionary;
} __attribute__((packed)) DecomposingField;

typedef struct
{
	SegmentDescriptor segmentDescriptor;
	unsigned char numberOfFields;
	void *field[NUMBER_OF_FIELDS_IN_FILE_HEADER_SEGMENT];
	char fieldsSignatures[NUMBER_OF_FIELDS_IN_FILE_HEADER_SEGMENT][NUMBER_OF_CHARS_IN_SIGNATURE];
} __attribute__((packed)) FileHeaderSegmentHeader;

#define EPOCH_SEGMENT_HEADER_SIGNATURE "hi"

typedef struct
{
	SegmentDescriptor 	segmentDescriptor;
	unsigned short int	epochNumber;
	unsigned       int	epochSize;
} __attribute__((packed)) EpochSegmentHeader;

#define SIGNAL_SEGMENT_HEADER_SIGNATURE "h"

typedef struct
{
	SegmentDescriptor  segmentDescriptor;
	unsigned short int channelNumber;
}  __attribute__((packed)) SignalSegmentHeader;

#define ATOMS_SEGMENT_HEADER_SIGNATURE "h"

typedef struct
{
	SegmentDescriptor  segmentDescriptor;
	unsigned short int channelNumber;
}  __attribute__((packed)) AtomsSegmentHeader;

#define DIRACDELTA_SIGNATURE    "fff"
#define GAUSSFUNCTION_SIGNATURE "ffff"
#define SINCOSWAVE_SIGNATURE    "ffff"
#define GABORWAVE_SIGNATURE     "ffffff"

static unsigned char getSizeOf(const char *signature);
static unsigned char countAtomSize(unsigned char type);
static unsigned int  getSizeOfAtomField(const unsigned char type);
static unsigned int  getSizeOfAtomsFields(MP5Parameters *mp5Parameters);
static unsigned int  getSizeOfAtomsSegment(MP5Parameters *mp5Parameters);
static unsigned int  getSizeOfSignalSegment(MP5Parameters *mp5Parameters);
static unsigned int  getSizeOfEpochSegment(MP5Parameters *mp5Parameters);

#ifdef INTELSWP
	#if  !defined(__BORLANDC__) && !defined(WIN32)
		#ifndef __USE_XOPEN
			#define __USE_XOPEN
			#include <unistd.h>
			#undef __USE_XOPEN
		#else
			#include <unistd.h>
		#endif
	#endif

	static void intelShortToJavaShort(unsigned short int *number)
	{
		swab((void *)number,(void *)number,sizeof(unsigned short int));
	}

	static void intelIntToJavaInt(int intelInt, char *javaInt)
	{
		short isav;
		union
		{
			char  fc[4];
			short fsi[2];
			int   ff;
		} gf[2];

		gf[0].ff=intelInt;
		swab((void *)gf[0].fc,(void *)gf[1].fc,4);
		isav=gf[1].fsi[0];
		gf[1].fsi[0]=gf[1].fsi[1]; gf[1].fsi[1]=isav;
		memcpy((void *)javaInt,(void *)&gf[1].ff,4);
	}

	static void intelFloatToJavaFloat(float intelFloat, char *javaFloat)
	{
		short isav;
		union
		{
			char  fc[4];
			short fsi[2];
			float ff;
		} gf[2];

		gf[0].ff=intelFloat;
		swab((void *)gf[0].fc,(void *)gf[1].fc,4);
		isav=gf[1].fsi[0];
		gf[1].fsi[0]=gf[1].fsi[1]; gf[1].fsi[1]=isav;
		memcpy((void *)javaFloat,(void *)&gf[1].ff,4);
	}

	static void intelStructToJavaStruct(const void *intelStruct, void *javaStruct, const char *descriptorSignature, const char *dataSignature)
	{

		unsigned int       counter;
		unsigned short int tmpShort;
		int                tmpIntelInt,   tmpJavaInt;
		float              tmpIntelFloat, tmpJavaFloat;
		unsigned int       position = 0;

		for(counter=0 ; descriptorSignature[counter]!='\0' ; counter++)
		{
			switch(descriptorSignature[counter])
			{
				case 'c':
					memmove((javaStruct + position),(intelStruct + position),sizeof(char));
					position++;
					break;
				case 'h':
					memcpy((void *)&tmpShort,(intelStruct + position),sizeof(short));
					intelShortToJavaShort(&tmpShort);
					memcpy((javaStruct + position),(void *)&tmpShort,sizeof(short));
					position+=sizeof(short);
					break;
				case 'i':
					memcpy((void *)&tmpIntelInt,(intelStruct + position),sizeof(int));
					intelIntToJavaInt(tmpIntelInt,(char *)&tmpJavaInt);
					memcpy((javaStruct + position),(void *)&tmpJavaInt,sizeof(int));
					position+=sizeof(int);
					break;
				case 'f':
					memcpy((void *)&tmpIntelFloat,(intelStruct + position),sizeof(float));
					intelFloatToJavaFloat(tmpIntelFloat,(char *)&tmpJavaFloat);
					memcpy((javaStruct + position),(void *)&tmpJavaFloat,sizeof(float));
					position+=sizeof(float);
					break;
				default: break;
			}
		}

		for(counter=0 ; dataSignature[counter]!='\0' ; counter++)
		{
			switch(dataSignature[counter])
			{
				case 'c':
					memcpy((javaStruct + position),(intelStruct + position),sizeof(char));
					position++;
					break;
				case 'h':
					memcpy((void *)&tmpShort,(intelStruct + position),sizeof(short));
					intelShortToJavaShort(&tmpShort);
					memcpy((javaStruct + position),(void *)&tmpShort,sizeof(short));
					position+=sizeof(short);
					break;
				case 'i':
					memcpy((void *)&tmpIntelInt,(intelStruct + position),sizeof(int));
					intelIntToJavaInt(tmpIntelInt,(char *)&tmpJavaInt);
					memcpy((javaStruct + position),(void *)&tmpJavaInt,sizeof(int));
					position+=sizeof(int);
					break;
				case 'f':
					memcpy((void *)&tmpIntelFloat,(intelStruct + position),sizeof(float));
					intelFloatToJavaFloat(tmpIntelFloat,(char *)&tmpJavaFloat);
					memcpy((javaStruct + position),(void *)&tmpJavaFloat,sizeof(float));
					position+=sizeof(float);
					break;
				default: break;
			}
		}
	}
#endif

static unsigned char countAtomSize(unsigned char type)
{
	if(type & DIRACDELTA)
		return 3*sizeof(float);
	else if(type & GAUSSFUNCTION)
		return 4*sizeof(float);
	else if(type & SINCOSWAVE)
		return 4*sizeof(float);
	else if(type & GABORWAVE)
		return 6*sizeof(float);
	else
		return 0;
}

static unsigned char getSizeOf(const char *signature)
{
	/*  definitions of chars in signature
		'c' - char/unsigned char
		'h' - short int/unsigned short int
		'i' - int/unsigned int
		'f' - float
	*/
	unsigned char counter;
	unsigned char size = 0;

	for(counter=0 ; signature[counter]!='\0' ; counter++)
	{
		switch(signature[counter])
		{
			case 'c': size++;              break;
			case 'h': size+=sizeof(short); break;
			case 'i': size+=sizeof(int);   break;
			case 'f': size+=sizeof(float); break;
			default: break;
		}
	}

	return size;
}

static STATUS writeField(void *field, const char* fieldDataSignature, FILE* resultsFile)
{
	#ifdef INTELSWP
		intelStructToJavaStruct(field,
								field,
								FIELD_DESCRIPTOR_SIGNATURE,
								fieldDataSignature);
	#endif

	if(fwrite(field,getSizeOf(FIELD_DESCRIPTOR_SIGNATURE) + getSizeOf(fieldDataSignature),1,resultsFile)!=1)
		return ERROR;

	return SUCCESS;
}

static STATUS initWebSiteLinkFieldAndDateField(WebSiteLinkField *webSiteLinkField, DateField *dateField)
{
	time_t 		  tmpTime;
	char          tmpBuffer[256];
	unsigned int  lengthOfData;

	time(&tmpTime);
	strcpy(tmpBuffer,(char *)ctime(&tmpTime));
	tmpBuffer[strlen(tmpBuffer)-1]='\0';

	lengthOfData = strlen(tmpBuffer);

	dateField->date = (char *)malloc(lengthOfData*sizeof(char));
	strncpy((char *)(dateField->date),tmpBuffer,lengthOfData);

	(dateField->fieldDescriptor).codeOfField = DATE_FIELD_IDENTITY;
	(dateField->fieldDescriptor).sizeOfFieldData = (unsigned char)(lengthOfData*sizeof(char));

	lengthOfData = strlen(WEB_SITE_LINK_TEXT);

	if(lengthOfData>256)
		return ERROR;

	webSiteLinkField->webSiteLink = (char *)malloc(lengthOfData*sizeof(char));
	strncpy((char *)(webSiteLinkField->webSiteLink),WEB_SITE_LINK_TEXT,lengthOfData);

	(webSiteLinkField->fieldDescriptor).codeOfField = WEB_SITE_LINK_FIELD_IDENTITY;
	(webSiteLinkField->fieldDescriptor).sizeOfFieldData = (unsigned char)(lengthOfData*sizeof(char));

	return SUCCESS;
}

static void initFileHeaderSegmentHeader(FileHeaderSegmentHeader *fileHeaderSegmentHeader)
{
	(fileHeaderSegmentHeader->segmentDescriptor).codeOfSegment = FILE_HEADER_SEGMENT_IDENTITY;
	(fileHeaderSegmentHeader->segmentDescriptor).sizeOfSegmentData    = 0;
	fileHeaderSegmentHeader->numberOfFields = 0;
}

static STATUS writeWebSiteLinkFieldAndDateField(WebSiteLinkField *webSiteLinkField, DateField *dateField, FILE *fileResults)
{
	if(fwrite((void *)&(webSiteLinkField->fieldDescriptor),getSizeOf(FIELD_DESCRIPTOR_SIGNATURE),1,fileResults)!=1)
		return ERROR;
	if(fwrite((void *)(webSiteLinkField->webSiteLink),(webSiteLinkField->fieldDescriptor).sizeOfFieldData,1,fileResults)!=1)
		return ERROR;

	if(fwrite((void *)&(dateField->fieldDescriptor),getSizeOf(FIELD_DESCRIPTOR_SIGNATURE),1,fileResults)!=1)
		return ERROR;
	if(fwrite((void *)(dateField->date),(dateField->fieldDescriptor).sizeOfFieldData,1,fileResults)!=1)
		return ERROR;

	return SUCCESS;
}

static void freeWebSiteLinkFieldAndDateField(WebSiteLinkField *webSiteLinkField, DateField *dateField)
{
	if(webSiteLinkField->webSiteLink!=NULL)
		free(webSiteLinkField->webSiteLink);

	if(dateField->date!=NULL)
		free(dateField->date);
}

static void addFieldToFileHeaderSegment(FileHeaderSegmentHeader *fileHeaderSegmentHeader, void* field, const char *fieldSignature, const unsigned char codeOfField)
{
	FieldDescriptor fieldDescriptor;

	fieldDescriptor.codeOfField      = codeOfField;
	fieldDescriptor.sizeOfFieldData  = getSizeOf(fieldSignature);

	memcpy(field,(void *)&fieldDescriptor,getSizeOf(FIELD_DESCRIPTOR_SIGNATURE));

	fileHeaderSegmentHeader->field[fileHeaderSegmentHeader->numberOfFields] = field;
	strcpy(fileHeaderSegmentHeader->fieldsSignatures[fileHeaderSegmentHeader->numberOfFields],fieldSignature);
	fileHeaderSegmentHeader->numberOfFields = fileHeaderSegmentHeader->numberOfFields + 1;
}

void returnAmplitudeAndModulusForMMP2DI(MP5Parameters *mp5Parameters, Dictionary *dictionary, Atom *atom, float *amplitude, float *modulus, unsigned int channelNumber)
{
	unsigned int sample;
	unsigned int epochExpandedSize = mp5Parameters->epochExpandedSize;

	const double phase = *(atom->phase + channelNumber);
	const double KS = atom->KS;
	const double KC = atom->KC;
	const double KM = atom->KM;

	double sinPhase, cosPhase;

	sincos(phase,&sinPhase,&cosPhase);

	const double sinPart    = KS*(sinPhase*sinPhase);
	const double cosPart    = KC*(cosPhase*cosPhase);
	const double sinCosPart = 2.0*KM*sinPhase*cosPhase;

	double *signalInParticularChannel = *(mp5Parameters->multiChannelSignalTable + channelNumber);
	double *prevAtomTable             = *(mp5Parameters->prevAtomTable + channelNumber);

	makeSinCosExpAtomTable(dictionary,mp5Parameters,atom);
	makeAtomTable(prevAtomTable,mp5Parameters,atom,channelNumber);

	*modulus = 0.0;

	for(sample=0;sample<epochExpandedSize;sample++)
		(*modulus)+= (*(signalInParticularChannel + sample))*(*(prevAtomTable + sample));

	findResidue(signalInParticularChannel,prevAtomTable,*modulus,epochExpandedSize);

	makeAtomTable(prevAtomTable,mp5Parameters,atom,channelNumber);
	*amplitude = (*modulus)/(2.0*sqrt((sinPart + cosPart) - sinCosPart));
	
}

static void returnAmplitudeAndModulusDI(Atom *atom, float *amplitude, float *modulus, unsigned short int channelNumber)
{
    const double RS    = (*(atom->RS    + channelNumber));
    const double RC    = (*(atom->RC    + channelNumber));
    const double phase = (*(atom->phase + channelNumber));

    const double KS = atom->KS;
    const double KC = atom->KC;
    const double KM = atom->KM;

    double sinPhase, cosPhase;

    sincos(phase,&sinPhase,&cosPhase);

    const double sinPart    = KS*(sinPhase*sinPhase);
    const double cosPart    = KC*(cosPhase*cosPhase);
    const double sinCosPart = 2.0*KM*sinPhase*cosPhase;

    *modulus   = (RC*cosPhase - RS*sinPhase)/((float)sqrt((sinPart + cosPart) - sinCosPart));
    *amplitude = (*modulus)/(2.0*sqrt((sinPart + cosPart ) - sinCosPart));

}

static STATUS writeAtom(MP5Parameters *mp5Parameters, Dictionary *dictionary, Atom *atom, unsigned short int channelNumber)
{

	/* Generally, a static allocation of an array inside function's body is forbidden in C/C++ language.
	     The array shold be allocated dynamically, however, then new version of GNU C/C++ let the static allocation
	     of the arrays inside the functions  */

	char  type;

	float modulus;
	float amplitude;
	float position;
	float scale;
	float frequency;
	float phase;

	unsigned char fieldDataSignature[getSizeOf(FIELD_DESCRIPTOR_SIGNATURE) + getSizeOf(GABORWAVE_SIGNATURE)];
	unsigned char field[getSizeOf(FIELD_DESCRIPTOR_SIGNATURE) + getSizeOf(GABORWAVE_SIGNATURE)];
	unsigned char sizeOfFieldDescriptor = getSizeOf(FIELD_DESCRIPTOR_SIGNATURE);
	unsigned char sizeOfFloat = sizeof(float);

	//unsigned short int epochSize = mp5Parameters->epochSize;

	FieldDescriptor fieldDescriptor;

	type = (char)(atom->feature & 0x0F);

	if(mp5Parameters->MPType & SMP)
	{
		returnAmplitudeAndModulusDI(atom,&amplitude,&modulus,0);
		phase      = (float)(*(atom->phase));
	}
	else if((mp5Parameters->MPType & MMP2)  || (mp5Parameters->MPType & MMP22))
	{
		returnAmplitudeAndModulusForMMP2DI(mp5Parameters,dictionary,atom,&amplitude,&modulus,channelNumber);
		phase      = (float)(*(atom->phase));
	}
	else if((mp5Parameters->MPType & MMP12) || (mp5Parameters->MPType & MMP21) ||
 		    (mp5Parameters->MPType & MMP23) || (mp5Parameters->MPType & MMP32))
	{
		returnAmplitudeAndModulusForMMP2DI(mp5Parameters,dictionary,atom,&amplitude,&modulus,channelNumber);
		phase      = (float)(*(atom->phase + channelNumber));
	}
	else
	{
		returnAmplitudeAndModulusDI(atom,&amplitude,&modulus,channelNumber);
		phase      = (float)(*(atom->phase + channelNumber));
	}

	amplitude  = (float)(2.0*amplitude);
	position   = (float)(atom->position);
	scale      = (float)(*(dictionary->tableOfScalesInOptimalDictionary + atom->scaleIndex));
	//frequency  = (float)(0.5*epochSize*(*(dictionary->tableOfFrequenciesInOptimalDictionary + atom->scaleIndex))*atom->rifling/M_PI);
	frequency  = (float)((*(dictionary->tableOfFrequenciesInOptimalDictionary + atom->scaleIndex))*atom->rifling/M_PI);

	if(type & DIRACDELTA)
	{
		amplitude = cos(phase)*amplitude; // we give up phase writing, but in case of Dirac delta and gauss Function the pase codes the sign og Dirac delta and GaussFunction
		fieldDescriptor.codeOfField     = DIRACDELTA_IDENTITY;
		fieldDescriptor.sizeOfFieldData = getSizeOf(DIRACDELTA_SIGNATURE);

		memcpy((void *)&field,(void *)&fieldDescriptor,sizeOfFieldDescriptor);
		memcpy(((void *)&field + sizeOfFieldDescriptor),(void *)&modulus,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 1*sizeOfFloat),(void *)&amplitude,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 2*sizeOfFloat),(void *)&position,sizeof(float));

		strcpy((char *)fieldDataSignature,DIRACDELTA_SIGNATURE);
	}
	else if(type & GAUSSFUNCTION)
	{
		amplitude = cos(phase)*amplitude; // we give up phase writing, but in case of Dirac delta and gauss Function the pase codes the sign of Dirac delta and GaussFunction
		fieldDescriptor.codeOfField     = GAUSSFUNCTION_IDENTITY;
		fieldDescriptor.sizeOfFieldData = getSizeOf(GAUSSFUNCTION_SIGNATURE);

		memcpy((void *)&field,(void *)&fieldDescriptor,sizeOfFieldDescriptor);
		memcpy(((void *)&field + sizeOfFieldDescriptor),(void *)&modulus,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 1*sizeOfFloat),(void *)&amplitude,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 2*sizeOfFloat),(void *)&position,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 3*sizeOfFloat),(void *)&scale,sizeof(float));

		strcpy((char *)fieldDataSignature,GAUSSFUNCTION_SIGNATURE);
	}
	else if(type & SINCOSWAVE)
	{
		fieldDescriptor.codeOfField     = SINCOSWAVE_IDENTITY;
		fieldDescriptor.sizeOfFieldData = getSizeOf(SINCOSWAVE_SIGNATURE);

		memcpy((void *)&field,(void *)&fieldDescriptor,sizeOfFieldDescriptor);
		memcpy(((void *)&field + sizeOfFieldDescriptor),(void *)&modulus,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 1*sizeOfFloat),(void *)&amplitude,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 2*sizeOfFloat),(void *)&frequency,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 3*sizeOfFloat),(void *)&phase,sizeof(float));

		strcpy((char *)fieldDataSignature,SINCOSWAVE_SIGNATURE);
	}
	else if(type & GABORWAVE)
	{
		fieldDescriptor.codeOfField     = GABORWAVE_IDENTITY;
		fieldDescriptor.sizeOfFieldData = getSizeOf(GABORWAVE_SIGNATURE);

		memcpy((void *)&field,(void *)&fieldDescriptor,sizeOfFieldDescriptor);
		memcpy(((void *)&field + sizeOfFieldDescriptor),(void *)&modulus,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 1*sizeOfFloat),(void *)&amplitude,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 2*sizeOfFloat),(void *)&position,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 3*sizeOfFloat),(void *)&scale,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 4*sizeOfFloat),(void *)&frequency,sizeof(float));
		memcpy(((void *)&field + sizeOfFieldDescriptor + 5*sizeOfFloat),(void *)&phase,sizeof(float));

		strcpy((char *)fieldDataSignature,GABORWAVE_SIGNATURE);
	}

	if(writeField((void *)&field,(const char *)fieldDataSignature,mp5Parameters->resultsFile)==ERROR)
		return ERROR;

	fflush(stdout);

	return SUCCESS;
}

static unsigned char getNumberOfDigits(unsigned int number)
{
	if (number == 0)
		return 1;

	unsigned int k = 1, digits = 0;

	while (number > (k - 1))
	{
		k *= 10;
		digits++;
	}

	return digits + 1;
}


static void printSelectedChannels(MP5Parameters *mp5Parameters)
{
    int i;
    const int howMany = 10;
    BOOLEAN breakLine = FALSE;
	unsigned short int numberOfSelectedChannels = mp5Parameters->numberOfSelectedChannels;
	char format[numberOfSelectedChannels + 6];

	unsigned short int numberOfDigits = getNumberOfDigits(numberOfSelectedChannels*numberOfSelectedChannels);
	sprintf(format,"%%-%huu",numberOfDigits);

	printf("                               ");

    for(i=0;i<mp5Parameters->numberOfSelectedChannels;i++)
    {
		breakLine = FALSE;

		if((i%howMany)==0 && i>0)
		{
			printf("\n");
			printf("                             ");
			breakLine = TRUE;
		}

		printf(format,mp5Parameters->selectedChannels[i]);
    }

    if(!breakLine)
		printf("\n");
}

static void printSelectedEpochs(MP5Parameters *mp5Parameters)
{
    int i;
    const int howMany = 10;
    BOOLEAN breakLine = FALSE;
	unsigned short int numberOfSelectedEpochs = mp5Parameters->numberOfSelectedEpochs;
	char format[numberOfSelectedEpochs + 6];

	unsigned short int numberOfDigits = getNumberOfDigits(numberOfSelectedEpochs*numberOfSelectedEpochs);
	sprintf(format,"%%-%huu",numberOfDigits);

    printf("                               ");

    for(i=0;i<mp5Parameters->numberOfSelectedEpochs;i++)
    {
		breakLine = FALSE;

		if((i%howMany)==0 && i>0)
		{
			printf("\n");
			printf("                               ");
			breakLine = TRUE;
		}

		printf("%-7d ",mp5Parameters->selectedEpochs[i]);
    }

    if(!breakLine)
		printf("\n");
}

STATUS testFilesAndDirectories(MP5Parameters *mp5Parameters, const ConfigFile *configFile, char *infoMessage)
{
    char *nameOfResultsFile = mp5Parameters->nameOfResultsFile;

    DIR  *directory;
    FILE *resultsFile;

    if((directory = opendir(mp5Parameters->nameOfOutputDirectory))==NULL)
    {
		const char *tmpString[] = {mp5Parameters->nameOfOutputDirectory};
		printError(infoMessage,CAN_NOT_OPEN_DIRECTORY,tmpString,1);
		return ERROR;
    }
    else
		closedir(directory);

	if(mp5Parameters->writingMode & CREATE_FILE)
	{
		resultsFile = fopen(nameOfResultsFile,"rb");;

		if(resultsFile!=NULL)
		{
			const char *tmpString[] = {mp5Parameters->nameOfResultsFile};
			printError(infoMessage,OVERWRITE_RESULTS_ALARM,tmpString,1);
			fclose(resultsFile);
			return ERROR;
		}
	}

    return SUCCESS;
}

STATUS testMP5Parameters(Dictionary *dictionary, MP5Parameters *mp5Parameters, char *infoMessage)
{

    if(mp5Parameters->numberOfSelectedChannels > mp5Parameters->numberOfChannelsInDataFile)
    {
		sprintf(aTmp,"%hu",mp5Parameters->numberOfSelectedChannels);
		sprintf(bTmp,"%hu",mp5Parameters->numberOfChannelsInDataFile);
		const char *tmpString[] = {aTmp,bTmp,mp5Parameters->nameOfDataFile};
		printError(infoMessage,BAD_NUMBER_OF_SELECTED_CHANNELS,tmpString,3);
		return ERROR;
    }

    if(mp5Parameters->energyPercent>=100.0)
    {
		printError(infoMessage,BAD_ENERGY_PERCENT,NULL,0);
		return ERROR;
    }

    if((mp5Parameters->reinitDictionary!=NO_REINIT_AT_ALL) && (dictionary->randomSeed!=AUTO_RANDOM_SEED))
    {
		printError(infoMessage,BAD_REINIT_ALL,NULL,0);
		return ERROR;
    }

    if(!(mp5Parameters->MPType & SMP) && (mp5Parameters->reinitDictionary==REINIT_IN_CHANNEL_DOMAIN || mp5Parameters->reinitDictionary==REINIT_AT_ALL))
    {
		printError(infoMessage,BAD_REINIT_MMP,NULL,0);
		return ERROR;
    }

    if(!(mp5Parameters->MPType & SMP) && (mp5Parameters->numberOfSelectedChannels==1))
    {
		printError(infoMessage,BAD_MP_ALGORITHM,NULL,0);
		return ERROR;
    }

    return SUCCESS;
}

void printInfoAboutData(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    printf(" \n");
    printf(" THE FOLLOWING PARAMETERS HAS BEEN READ: \n\n");
    printf(" NAME OF DATA FILE:            %s\n",mp5Parameters->nameOfDataFile);
    printf(" SAMPLE FREQUENCY:             %-4.2f\n",mp5Parameters->samplingFrequency);
    printf(" NUMBER OF CHANNELS IN FILE:   %-5hu\n",mp5Parameters->numberOfChannelsInDataFile);
    printf(" NUMBER OF SELECTED CHANNELS:  %-5hu\n",mp5Parameters->numberOfSelectedChannels);

    printf(" SELECTED CHANNELS:            \n");
    printSelectedChannels(mp5Parameters);
    printf(" NUMBER OF SAMPLES IN EPOCH:   %-5u\n",mp5Parameters->epochSize);
    printf(" NUMBER OF SELECTED EPOCHS:    %-5hu\n",mp5Parameters->numberOfSelectedEpochs);
    printf(" SELECTED EPOCHS:               \n");
    printSelectedEpochs(mp5Parameters);

    printf(" TYPE OF DICTIONARY:           ");
    if(dictionary->typeOfDictionary & OCTAVE_FIXED)
		printf("OCTAVE_FIXED\n");
    else if(dictionary->typeOfDictionary & OCTAVE_STOCH)
		printf("OCTAVE_STOCH\n");
    printf(" DILATION FACTOR:              %-f\n",dictionary->dilationFactor);
    printf(" REINIT DICTIONATY:            ");
    if(mp5Parameters->reinitDictionary & NO_REINIT_AT_ALL)
		printf("NO REINIT AT ALL \n");
    else if(mp5Parameters->reinitDictionary & REINIT_IN_CHANNEL_DOMAIN)
		printf("REINIT IN CHANNEL DOMAIN \n");
    else if(mp5Parameters->reinitDictionary & REINIT_IN_EPOCH_DOMAIN)
		printf("REINIT IN EPOCH DOMAIN \n");
    else if(mp5Parameters->reinitDictionary & REINIT_AT_ALL)
		printf("REINIT AT ALL\n");

    printf(" MAXIMAL NUMBER Of ITERATIONS  %-5hu\n",mp5Parameters->maximalNumberOfIterations);
    printf(" ENERGY PERCENT:               %-4.2f\n",mp5Parameters->energyPercent);
    printf(" TYPE OF ALGORITHM:            ");

    if(mp5Parameters->MPType & SMP)
		printf("SINLGE CHANNEL MATCHING PURSUIT \n");
    else if(mp5Parameters->MPType & MMP1)
		printf("MULTICHANNEL MATCHING PURSUIT I\n");
    else if(mp5Parameters->MPType & MMP12)
		printf("MULTICHANNEL MATCHING PURSUIT I-II\n");
    else if(mp5Parameters->MPType & MMP11)
		printf("MULTICHANNEL MULTITRIAL MATCHING PURSUIT I-I\n");
    else if(mp5Parameters->MPType & MMP21)
		printf("MULTICHANNEL MATCHING PURSUIT II-I\n");
	else if(mp5Parameters->MPType & MMP2)
		printf("MULTICHANNEL MATCHING PURSUIT II\n");
	else if(mp5Parameters->MPType & MMP22)
		printf("MULTICHANNEL MULTITRIAL MATCHING PURSUIT II-II\n");
	else if(mp5Parameters->MPType & MMP3)
		printf("MULTICHANNEL MATCHING PURSUIT III\n");
	else if(mp5Parameters->MPType & MMP23)
		printf("MULTICHANNEL MULTITRIAL MATCHING PURSUIT II-III\n");
	else if(mp5Parameters->MPType & MMP32)
		printf("MULTICHANNEL MULTITRIAL MATCHING PURSUIT III-II\n");
    else if(mp5Parameters->MPType & MMP33)
		printf("MULTICHANNEL MULTITRIAL MATCHING PURSUIT III-III\n");

    printf(" RESULTS WILL BE WRITTEN TO THE FOLLOWING FILE: \n");
	printf("                                                %s\n",mp5Parameters->nameOfResultsFile);
	if(mp5Parameters->FFT & ON)
		printf(" FAST FOURIER TRANSFORM:       ON\n");
	else
		printf(" FAST FOURIER TRANSFORM:       OFF\n");

    printf(" \n");
    fflush(stdout);
}

STATUS openBinaryDataFile(MP5Parameters *mp5Parameters, char *infoMessage)
{
    if(!(mp5Parameters->dataFile = fopen(mp5Parameters->nameOfDataFile,"rb")))
    {
		const char *tmpString[] = {mp5Parameters->nameOfDataFile};
		printError(infoMessage,CAN_NOT_OPEN_DATA_FILE,tmpString,1);
		return ERROR;
    }

    return SUCCESS;
}

STATUS openAsciiDataFile(MP5Parameters *mp5Parameters, char *infoMessage)
{
    if(!(mp5Parameters->dataFile = fopen(mp5Parameters->nameOfDataFile,"rt")))
    {
		const char *tmpString[] = {mp5Parameters->nameOfDataFile};
		printError(infoMessage,CAN_NOT_OPEN_DATA_FILE,tmpString,1);
		return ERROR;
    }

    return SUCCESS;
}

void createNamesOfResultFiles(Dictionary *dictionary, MP5Parameters *mp5Parameters)
{
    unsigned short int charCounter;
	char *nameOfResultsFile = mp5Parameters->nameOfResultsFile;
	unsigned short int lengthOfDataFileWithOutExpand;

    for(charCounter=0;charCounter<LENGTH_OF_NAME_OF_RESULTS_FILE;charCounter++)
		*(nameOfResultsFile + charCounter) = '\0';

	char *dot = strrchr(mp5Parameters->nameOfDataFile,'.');
	int len;

	char tmpString[LENGTH_OF_TMP_STRING];

	bzero((void *)tmpString,LENGTH_OF_TMP_STRING);

	strcpy(nameOfResultsFile,mp5Parameters->nameOfOutputDirectory);

	if(dot==NULL)
    {
		strcpy(nameOfResultsFile,mp5Parameters->nameOfDataFile);

        if(!(mp5Parameters->MPType & SMP))
            sprintf(nameOfResultsFile,"_mmp.b");
        else
            sprintf(nameOfResultsFile,"_smp.b");
    }
    else
    {
		len = strlen(dot);
    	strncat(nameOfResultsFile,mp5Parameters->nameOfDataFile,strlen(mp5Parameters->nameOfDataFile)-len);
    	lengthOfDataFileWithOutExpand = (unsigned short int)strlen(nameOfResultsFile);

        if(!(mp5Parameters->MPType & SMP))
            sprintf((nameOfResultsFile + lengthOfDataFileWithOutExpand),"_mmp.b");
        else
            sprintf((nameOfResultsFile+ lengthOfDataFileWithOutExpand),"_smp.b");
	}
}

STATUS openResultFiles(MP5Parameters *mp5Parameters, char *infoMessage)
{
	char *nameOfResultsFile = mp5Parameters->nameOfResultsFile;

	if(mp5Parameters->writingMode & CREATE_FILE)
	{
		if(!(mp5Parameters->resultsFile = fopen(nameOfResultsFile,"wb")))
		{
			const char *tmpString[] = {nameOfResultsFile};
			printError(infoMessage,CAN_NOT_OPEN_RESULTS_FILE,tmpString,1);
			return ERROR;
		}
	}
	else if(mp5Parameters->writingMode & APPEND_FILE)
	{
		if(!(mp5Parameters->resultsFile = fopen(nameOfResultsFile,"a+b")))
		{
			const char *tmpString[] = {nameOfResultsFile};
			printError(infoMessage,CAN_NOT_OPEN_RESULTS_FILE,tmpString,1);
			return ERROR;
		}
	}
	return SUCCESS;
}

void closeFiles(MP5Parameters *mp5Parameters)
{
	if(mp5Parameters->dataFile!=NULL)
		fclose(mp5Parameters->dataFile);

	if(mp5Parameters->resultsFile!=NULL)
		fclose(mp5Parameters->resultsFile);
}

STATUS analyseBinaryDataFile(MP5Parameters *mp5Parameters, char *infoMessage)
{
    unsigned long int sizeOfFile;
    unsigned long int sizeOfData;

    fseek(mp5Parameters->dataFile,0L,SEEK_END);
    sizeOfFile = ftell(mp5Parameters->dataFile);
    fseek(mp5Parameters->dataFile,0L,SEEK_SET);

    sizeOfData = sizeOfFile - mp5Parameters->sizeOfHeader - mp5Parameters->sizeOfTail;

    if((sizeOfData/sizeof(float))%mp5Parameters->numberOfChannelsInDataFile!=0)
    {
        if(applicationMode & PROCESS_USER_MODE)
        {
            printf("\n DATA WARNING: \n CHANNELS DO NOT CONTAIN THE SAME NUMBER OF SAMPLES IN FILE: %s\n THE ERROR WILL OCCUR ERROR WHILE READING THE LAST EPOCH\n",mp5Parameters->nameOfDataFile);
            fflush(stdout);
            sleep(1);
        }
        else if(applicationMode & PROCESS_SERVER_MODE)
		{
		  const char *tmpString[] = {mp5Parameters->nameOfDataFile};
		  printError(infoMessage,BAD_NUMBER_OF_SAMPLES_PER_CHANNEL,tmpString,1);
          return ERROR;
        }
    }

    if(sizeOfData%sizeof(float)!=0)
    {
        const char *tmpString[] = {mp5Parameters->nameOfDataFile};
        printError(infoMessage,BAD_NUMBER_OF_SAMPLES,tmpString,1);
        return ERROR;
    }
	else
    {
        mp5Parameters->numberOfPoints  = sizeOfData/sizeof(float)/mp5Parameters->numberOfChannelsInDataFile;
        mp5Parameters->numberOfEpochs = (unsigned short int)(mp5Parameters->numberOfPoints/mp5Parameters->epochSize);

        if(mp5Parameters->numberOfEpochs < mp5Parameters->numberOfSelectedEpochs)
        {
            sprintf(aTmp,"%hu",mp5Parameters->numberOfSelectedEpochs);
            sprintf(bTmp,"%hu",mp5Parameters->numberOfEpochs);
            const char *tmpString[] = {aTmp,bTmp,mp5Parameters->nameOfDataFile};
            printError(infoMessage,BAD_NUMBER_OF_SELECTED_EPOCHS,tmpString,3);
            return ERROR;
        }

        mp5Parameters->samplesBesideEpochs = mp5Parameters->numberOfPoints - mp5Parameters->numberOfEpochs*mp5Parameters->numberOfPoints;
    }

    return SUCCESS;
}

static void processRawData(MP5Parameters *mp5Parameters)
{
    const unsigned       int numberOfReadChannelsAndEpochs = mp5Parameters->numberOfReadChannelsAndEpochs;
    const unsigned       int epochSize              = mp5Parameters->epochSize;
    const unsigned       int marginalSize           = mp5Parameters->marginalSize;
    const unsigned       int epochExpandedSize = mp5Parameters->epochExpandedSize;
    double               **rawDataMatrix       = mp5Parameters->rawDataMatrix;
    double               **processedDataMatrix = mp5Parameters->processedDataMatrix;
    unsigned       int channel;
    unsigned       int sample;

    dSetMatrixZero(processedDataMatrix,numberOfReadChannelsAndEpochs,epochExpandedSize);

    for(channel=0;channel<numberOfReadChannelsAndEpochs;channel++)
		for(sample=0;sample<epochSize;sample++)
			*(*(processedDataMatrix + channel) + marginalSize + sample) = *(*(rawDataMatrix + channel) + sample);
					
}

static STATUS writeSegmentHeader(void *segmentHeader, const char* segmentDataSignature, FILE* resultsFile)
{
	#ifdef INTELSWP
		intelStructToJavaStruct(segmentHeader,
								segmentHeader,
								SEGMENT_DESCRIPTOR_SIGNATURE,
								segmentDataSignature);
	#endif

	if(fwrite(segmentHeader,getSizeOf(SEGMENT_DESCRIPTOR_SIGNATURE) + getSizeOf(segmentDataSignature),1,resultsFile)!=1)
		return ERROR;

	return SUCCESS;
}

static unsigned int getSizeOfAtomField(const unsigned char type)
{
	return getSizeOf(FIELD_DESCRIPTOR_SIGNATURE) + countAtomSize(type);
}

static unsigned int getSizeOfAtomsFields(MP5Parameters *mp5Parameters)
{
	Atom *atom;
	unsigned int atomNumber;
	unsigned int sizeOfAtomsFields = 0;

	for(atomNumber=0;atomNumber<mp5Parameters->fitted->size;atomNumber++)
	{
		atom = (Atom *)readNNode(mp5Parameters->fitted,atomNumber);

		sizeOfAtomsFields+=getSizeOfAtomField(atom->feature);
	}

	return sizeOfAtomsFields;
}

static unsigned int getSizeOfAtomsSegment(MP5Parameters *mp5Parameters)
{
	unsigned short int numberOfProceesedChannels;

	if(mp5Parameters->MPType == SMP)
		numberOfProceesedChannels = 1;
	else
		numberOfProceesedChannels = mp5Parameters->numberOfSelectedChannels;


	return numberOfProceesedChannels*(getSizeOf(SEGMENT_DESCRIPTOR_SIGNATURE) +
									  getSizeOf(ATOMS_SEGMENT_HEADER_SIGNATURE) +
									  getSizeOfAtomsFields(mp5Parameters));

}

static unsigned int getSizeOfSignalSegment(MP5Parameters *mp5Parameters)
{
	unsigned int sizeOfSignalSegment = 0;

	if(mp5Parameters->bookWithSignal & YES)
	{
		unsigned short int numberOfProceesedChannels;

		if(mp5Parameters->MPType == SMP)
			numberOfProceesedChannels = 1;
		else
			numberOfProceesedChannels = mp5Parameters->numberOfSelectedChannels;


		sizeOfSignalSegment = numberOfProceesedChannels*(getSizeOf(SEGMENT_DESCRIPTOR_SIGNATURE) +
													     getSizeOf(SIGNAL_SEGMENT_HEADER_SIGNATURE) +
														 mp5Parameters->epochSize*sizeof(float));
	}

	return sizeOfSignalSegment;

}

static unsigned int getSizeOfEpochSegment(MP5Parameters *mp5Parameters)
{
	return getSizeOf(EPOCH_SEGMENT_HEADER_SIGNATURE) +
		   getSizeOfSignalSegment(mp5Parameters) +
		   getSizeOfAtomsSegment(mp5Parameters);
}

static STATUS writeEpochSegment(Dictionary *dictionary, MP5Parameters *mp5Parameters, unsigned short int epochNumber, unsigned short int channelNumber, char *infoMessage)
{
	Atom *atom;
    unsigned int atomNumber      = 0;
    const unsigned int epochSize = mp5Parameters->epochSize;
	unsigned int selectedEpochChannel;

	SignalSegmentHeader  signalSegmentHeader;
	AtomsSegmentHeader   atomsSegmentHeader;

	if(mp5Parameters->bookWithSignal & YES)
	{
		(signalSegmentHeader.segmentDescriptor).codeOfSegment     = SIGNAL_SEGMENT_IDENTITY;
		(signalSegmentHeader.segmentDescriptor).sizeOfSegmentData = getSizeOf(SIGNAL_SEGMENT_HEADER_SIGNATURE) + epochSize*sizeof(float);
		signalSegmentHeader.channelNumber                   	  = *(mp5Parameters->selectedChannels + channelNumber);
	}

	(atomsSegmentHeader.segmentDescriptor).codeOfSegment     = ATOMS_SEGMENT_IDENTITY;
	(atomsSegmentHeader.segmentDescriptor).sizeOfSegmentData = getSizeOf(ATOMS_SEGMENT_HEADER_SIGNATURE) + getSizeOfAtomsFields(mp5Parameters);
	atomsSegmentHeader.channelNumber                    	 = *(mp5Parameters->selectedChannels + channelNumber);

	if(mp5Parameters->bookWithSignal & YES)
	{
		unsigned int 	 counter;
		double  *rawDataMatrix;
		float   intelFloat, javaFloat;

		if(writeSegmentHeader((void *)&signalSegmentHeader,SIGNAL_SEGMENT_HEADER_SIGNATURE,mp5Parameters->resultsFile)==ERROR)
			return ERROR;

		if(readDataFileOneTrial(mp5Parameters,mp5Parameters->selectedEpochs[epochNumber],infoMessage)==ERROR)
			return ERROR;

		rawDataMatrix = *(mp5Parameters->rawDataMatrix + channelNumber);

		for(counter = 0;counter<epochSize;counter++)
		{
			intelFloat = (float)(*(rawDataMatrix + counter));

			intelFloatToJavaFloat(intelFloat,(char *)&javaFloat);
			if(fwrite((void *)&javaFloat,sizeof(float),1,mp5Parameters->resultsFile)==0U)
				return ERROR;
		}
	}

	if(writeSegmentHeader((void *)&atomsSegmentHeader,ATOMS_SEGMENT_HEADER_SIGNATURE,mp5Parameters->resultsFile)==ERROR)
		return ERROR;

	fflush(mp5Parameters->resultsFile);

	if((mp5Parameters->MPType & MMP11) || (mp5Parameters->MPType & MMP22) || (mp5Parameters->MPType & MMP33))
		selectedEpochChannel = epochNumber*mp5Parameters->numberOfSelectedChannels + channelNumber;
	else
		selectedEpochChannel = channelNumber;

	for(atomNumber=0;atomNumber<mp5Parameters->fitted->size;atomNumber++)
	{
		atom = (Atom *)readNNode(mp5Parameters->fitted,atomNumber);

		if(writeAtom(mp5Parameters,dictionary,atom,selectedEpochChannel)==ERROR)
			return ERROR;

		fflush(mp5Parameters->resultsFile);
	}

	return SUCCESS;
}

STATUS writeMagic(MP5Parameters *mp5Parameters, char *infoMessage)
{
	char tmpBuffer[6];

	strncpy(tmpBuffer,MAGIC_TEXT,6);

	if(fwrite((void *)tmpBuffer,6*sizeof(char),1,mp5Parameters->resultsFile)!=1)
	{
		const char *tmpString[] = {mp5Parameters->nameOfResultsFile};
		printError(infoMessage,CAN_NOT_WRITE_HEADER,tmpString,1);
		return ERROR;
	}

	return SUCCESS;
}

STATUS writeCommentsSegment(ConfigFile *configFile, MP5Parameters *mp5Parameters, char *infoMessage)
{
	const char *tmpString[] = {mp5Parameters->nameOfResultsFile};

	SegmentDescriptor segmentDescriptor;
	unsigned short int numberOfChars = 0;

    char text[LENGTH_OF_LINE];

	if(fseek(configFile->file,0,0)!=0)
		goto ERROR_PROCEDURE;

    do
    {
		if(fgets(text,LENGTH_OF_LINE,configFile->file)==NULL)
			break;

		if(strstr(text,"##")!=NULL)
		{
			numberOfChars+=strlen(text);
		}
		else
			continue;
	}while(TRUE);

	if(fseek(configFile->file,0,0)!=0)
		goto ERROR_PROCEDURE;

	segmentDescriptor.codeOfSegment     = COMMENT_SEGMENT_IDENTITY;
	segmentDescriptor.sizeOfSegmentData = numberOfChars*sizeof(char);

	#ifdef INTELSWP
		intelStructToJavaStruct((void *)&(segmentDescriptor),
								(void *)&(segmentDescriptor),
								SEGMENT_DESCRIPTOR_SIGNATURE,
								"\0");
	#endif

	if(fwrite((void *)&segmentDescriptor,getSizeOf(SEGMENT_DESCRIPTOR_SIGNATURE),1,mp5Parameters->resultsFile)!=1)
		goto ERROR_PROCEDURE;

	do
    {
		if(fgets(text,LENGTH_OF_LINE,configFile->file)==NULL)
			break;

		if(strstr(text,"##")!=NULL)
		{
			if(fwrite((void *)&text,strlen(text)*sizeof(char),1,mp5Parameters->resultsFile)!=1)
				goto ERROR_PROCEDURE;
		}
		else
			continue;
	}while(TRUE);

	return SUCCESS;

	ERROR_PROCEDURE:
		printError(infoMessage,CAN_NOT_WRITE_RESULTS,tmpString,1);
		return ERROR;

}

STATUS writeFileHeader(Dictionary *dictionary, MP5Parameters *mp5Parameters, char *infoMessage)
{
	const char *tmpString[] = {mp5Parameters->nameOfResultsFile};

	unsigned char numberOfField;
	unsigned char sizeOfFieldDescriptor = getSizeOf(FIELD_DESCRIPTOR_SIGNATURE);
	unsigned char sizeOfData = 0;

	WebSiteLinkField  webSiteLinkField;
	DateField         dateField;
	SignalField       signalField;
	DecomposingField  decomposingField;

	FileHeaderSegmentHeader fileHeaderSegmentHeader;
	initFileHeaderSegmentHeader(&fileHeaderSegmentHeader);

	signalField.samplingFrequency          = (float)mp5Parameters->samplingFrequency;
	signalField.pointsPerMicrovolt         = (float)mp5Parameters->pointsPerMicrovolt;
	signalField.numberOfChannelsInDataFile = mp5Parameters->numberOfChannelsInDataFile;
	addFieldToFileHeaderSegment(&fileHeaderSegmentHeader, (void *)&signalField,SIGNAL_FIELD_SIGNATURE,SIGNAL_FIELD_IDENTITY);

    decomposingField.energyPercent             =  (float)mp5Parameters->energyPercent;
    decomposingField.maximalNumberOfIterations =  mp5Parameters->maximalNumberOfIterations;
    decomposingField.sizeOfDictionary          =  dictionary->finalNumberOfAtoms;
    decomposingField.typeOfDictionary          =  (char)((dictionary->typeOfDictionary & OCTAVE_FIXED) ? 'F' : 'S');
	addFieldToFileHeaderSegment(&fileHeaderSegmentHeader, (void *)&decomposingField,DECOMPOSING_FIELD_SIGNATURE,DECOMPOSING_FIELD_IDENTITY);

	if(!initWebSiteLinkFieldAndDateField(&webSiteLinkField,&dateField))
		goto ERROR_PROCEDURE;

	(fileHeaderSegmentHeader.segmentDescriptor).sizeOfSegmentData = 2*getSizeOf(FIELD_DESCRIPTOR_SIGNATURE) +
																	(webSiteLinkField.fieldDescriptor).sizeOfFieldData +
																	(dateField.fieldDescriptor).sizeOfFieldData;

	for(numberOfField = 0; numberOfField<fileHeaderSegmentHeader.numberOfFields; numberOfField++)
	{
		memcpy((void *)&sizeOfData,fileHeaderSegmentHeader.field[numberOfField] + sizeof(unsigned char),sizeof(unsigned char));
		(fileHeaderSegmentHeader.segmentDescriptor).sizeOfSegmentData+=sizeOfData + sizeOfFieldDescriptor;
	}

	#ifdef INTELSWP
		intelStructToJavaStruct((void *)&(fileHeaderSegmentHeader.segmentDescriptor),
								(void *)&(fileHeaderSegmentHeader.segmentDescriptor),
								SEGMENT_DESCRIPTOR_SIGNATURE,
								"\0");
	#endif

	if(fwrite((void *)&(fileHeaderSegmentHeader.segmentDescriptor),getSizeOf(SEGMENT_DESCRIPTOR_SIGNATURE),1,mp5Parameters->resultsFile)!=1)
		goto ERROR_PROCEDURE;

	if(writeWebSiteLinkFieldAndDateField(&webSiteLinkField,&dateField,mp5Parameters->resultsFile)!=1)
		goto ERROR_PROCEDURE;

	for(numberOfField = 0; numberOfField<fileHeaderSegmentHeader.numberOfFields; numberOfField++)
	{
		if(writeField(fileHeaderSegmentHeader.field[numberOfField],fileHeaderSegmentHeader.fieldsSignatures[numberOfField],mp5Parameters->resultsFile)==ERROR)
			goto ERROR_PROCEDURE;
	}

	freeWebSiteLinkFieldAndDateField(&webSiteLinkField,&dateField);
	fflush(mp5Parameters->resultsFile);

	return SUCCESS;

	ERROR_PROCEDURE:
		printError(infoMessage,CAN_NOT_WRITE_HEADER,tmpString,1);
		return ERROR;
}

STATUS writeSMPResults(Dictionary *dictionary, MP5Parameters *mp5Parameters, unsigned short int epochNumber, unsigned short int channelNumber, char *infoMessage)
{
	const char *tmpString[] = {mp5Parameters->nameOfResultsFile};
    const unsigned int epochSize = mp5Parameters->epochSize;

	EpochSegmentHeader  epochSegmentHeader;

	(epochSegmentHeader.segmentDescriptor).codeOfSegment      = EPOCH_SEGMENT_IDENTITY;
	(epochSegmentHeader.segmentDescriptor).sizeOfSegmentData  = getSizeOfEpochSegment(mp5Parameters);
	epochSegmentHeader.epochNumber                     	      = *(mp5Parameters->selectedEpochs + epochNumber);
	epochSegmentHeader.epochSize                          	  = epochSize;

	if(writeSegmentHeader((void *)&epochSegmentHeader,EPOCH_SEGMENT_HEADER_SIGNATURE,mp5Parameters->resultsFile)==ERROR)
		goto ERROR_PROCEDURE;

	if(writeEpochSegment(dictionary,mp5Parameters,epochNumber,channelNumber,infoMessage)==ERROR)
		goto ERROR_PROCEDURE;

    clearQueue(mp5Parameters->fitted,(void (*)(void *))freeAtom);

    return SUCCESS;

	ERROR_PROCEDURE:
		printError(infoMessage,CAN_NOT_WRITE_RESULTS,tmpString,1);
		return ERROR;
}

STATUS writeMMPResults(Dictionary *dictionary, MP5Parameters *mp5Parameters, unsigned short int epochNumber, char *infoMessage)
{
	const char *tmpString[] = {mp5Parameters->nameOfResultsFile};
	unsigned short int channelNumber;
    const unsigned int epochSize = mp5Parameters->epochSize;

	EpochSegmentHeader  epochSegmentHeader;

	(epochSegmentHeader.segmentDescriptor).codeOfSegment      = EPOCH_SEGMENT_IDENTITY;
	(epochSegmentHeader.segmentDescriptor).sizeOfSegmentData  = getSizeOfEpochSegment(mp5Parameters);
	epochSegmentHeader.epochNumber                     	      = *(mp5Parameters->selectedEpochs + epochNumber);
	epochSegmentHeader.epochSize                     	      = epochSize;

	if(writeSegmentHeader((void *)&epochSegmentHeader,EPOCH_SEGMENT_HEADER_SIGNATURE,mp5Parameters->resultsFile)==ERROR)
		goto ERROR_PROCEDURE;

    for(channelNumber=0;channelNumber<mp5Parameters->numberOfSelectedChannels;channelNumber++)
    {
		if(writeEpochSegment(dictionary,mp5Parameters,epochNumber,channelNumber,infoMessage)==ERROR)
			goto ERROR_PROCEDURE;
	}

	clearQueue(mp5Parameters->fitted,(void (*)(void *))freeAtom);

	return SUCCESS;

	ERROR_PROCEDURE:
		printError(infoMessage,CAN_NOT_WRITE_RESULTS,tmpString,1);
		return ERROR;
}

STATUS writeMMPMultiTrialResults(Dictionary *dictionary, MP5Parameters *mp5Parameters, char *infoMessage)
{
	const char *tmpString[] = {mp5Parameters->nameOfResultsFile};
	unsigned short int channelNumber;
    unsigned short int epochNumber;
	const unsigned int epochSize = mp5Parameters->epochSize;
	EpochSegmentHeader epochSegmentHeader;

	for(epochNumber=0;epochNumber<mp5Parameters->numberOfSelectedEpochs;epochNumber++)
	{
		(epochSegmentHeader.segmentDescriptor).codeOfSegment      = EPOCH_SEGMENT_IDENTITY;
		(epochSegmentHeader.segmentDescriptor).sizeOfSegmentData  = getSizeOfEpochSegment(mp5Parameters);
		epochSegmentHeader.epochNumber                     	      = *(mp5Parameters->selectedEpochs + epochNumber);
		epochSegmentHeader.epochSize                     	      = epochSize;

		if(writeSegmentHeader((void *)&epochSegmentHeader,EPOCH_SEGMENT_HEADER_SIGNATURE,mp5Parameters->resultsFile)==ERROR)
			goto ERROR_PROCEDURE;

		for(channelNumber=0;channelNumber<mp5Parameters->numberOfSelectedChannels;channelNumber++)
		{
			if(writeEpochSegment(dictionary,mp5Parameters,epochNumber,channelNumber,infoMessage)==ERROR)
				goto ERROR_PROCEDURE;
		}
	}

	clearQueue(mp5Parameters->fitted,(void (*)(void *))freeAtom);
	return SUCCESS;

	ERROR_PROCEDURE:
		printError(infoMessage,CAN_NOT_WRITE_RESULTS,tmpString,1);
	return ERROR;
}

static STATUS readBinaryData(MP5Parameters *mp5Parameters, unsigned short int epochNumber, unsigned int epochOnset, char *infoMessage)
{
    double             		 **rawDataMatrix = mp5Parameters->rawDataMatrix + epochOnset;
    const unsigned       int numberOfChannelsInDataFile = mp5Parameters->numberOfChannelsInDataFile;
    const unsigned       int numberOfSelectedChannels   = mp5Parameters->numberOfSelectedChannels;
    const unsigned       int epochSize = mp5Parameters->epochSize;
    unsigned       int channel;
    unsigned       int sample;
    unsigned short int formatSize = 0;
    unsigned long  int filePosition;

    float tmpData[numberOfChannelsInDataFile];

    formatSize = sizeof(float);

    filePosition = mp5Parameters->sizeOfHeader + (epochNumber-1)*epochSize*numberOfChannelsInDataFile*formatSize;
    fseek(mp5Parameters->dataFile,filePosition,SEEK_SET);

	for(sample=0;sample<epochSize;sample++)
	{
        if(fread((void *)tmpData,formatSize,numberOfChannelsInDataFile,mp5Parameters->dataFile)<numberOfChannelsInDataFile)
        {
            sprintf(aTmp,"%u",sample+1);
            sprintf(bTmp,"%hu",epochNumber);
            const char *tmpString[] = {aTmp,bTmp,mp5Parameters->nameOfResultsFile};
            printError(infoMessage,CAN_NOT_READ_EPOCH_IN_BINARY_FILE,tmpString,3);
            return ERROR;
        }

        for(channel=0;channel<numberOfSelectedChannels;channel++)
            *(*(rawDataMatrix + channel) + sample) = (double)(*(tmpData + mp5Parameters->selectedChannels[channel] - 1));
                                    
    }

    return SUCCESS;
}

STATUS readDataFileOneTrial(MP5Parameters *mp5Parameters, unsigned short int epochNumber, char *infoMessage)
{
    if(readBinaryData(mp5Parameters,epochNumber,0,infoMessage)==ERROR)
        return ERROR;

    /* copy data matrix (rawDataMatrix) of size numberOfChannelsInDataFile x epochSize to
	     processDataMatrix of size numberOfChannelsInDataFile x epochExpandedSize */

    processRawData(mp5Parameters);

    return SUCCESS;
}

STATUS readDataFileMultiTrial(MP5Parameters *mp5Parameters, char *infoMessage)
{
	unsigned short int epochNumber;
	unsigned short int selectedEpoch;
	unsigned int epochOnset = 0;
	unsigned short int numberOfSelectedEpochs = mp5Parameters->numberOfSelectedEpochs;
	unsigned short int numberOfSelectedChannels = mp5Parameters->numberOfSelectedChannels;

	for(epochNumber=0;epochNumber<numberOfSelectedEpochs;epochNumber++)
	{
		selectedEpoch = mp5Parameters->selectedEpochs[epochNumber];
		epochOnset    = epochNumber*numberOfSelectedChannels;

        if(readBinaryData(mp5Parameters,selectedEpoch,epochOnset,infoMessage)==ERROR)
            return ERROR;

	}

    processRawData(mp5Parameters);

    return SUCCESS;
}

void printInformationAboutProgress(MP5Parameters *mp5Parameters, Progress *progress, unsigned int atomsCounter)
{
	if(((atomsCounter+1)%(progress->stepInToolbar))==0)
	{
		if(progress->applicationMode & PROCESS_USER_MODE)
		{
			if(mp5Parameters->progressBar)
			{
				toolbar(progress->step);
				progress->step = progress->step + 1;
			}
		}
		else
		{
			printf("TESTED %u\n",atomsCounter);
			fflush(stdout);
		}
	}
		
}
