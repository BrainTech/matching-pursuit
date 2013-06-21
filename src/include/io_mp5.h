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

#ifndef _IO_MP5_H_

	#define _IO_MP5_H_

	#include<ctype.h>
	#include<stdlib.h>
	#include<stdio.h>
	#include<string.h>
	#include<time.h>
	#include"types.h"

	STATUS generateErrotTextAndCodeGeneratorFile();
	void printError(char *infoMessage, unsigned short int errorNumber, const char *acceptableValues[], unsigned short int sizeOfAcceptableValues);
	STATUS testFilesAndDirectories(MP5Parameters *mp5Parameters, const ConfigFile *configFile, char *infoMessage);
	STATUS testMP5Parameters(Dictionary *dictionary, MP5Parameters *mp5Parameters, char *infoMessage);
	STATUS openBinaryDataFile(MP5Parameters *mp5Parameters, char *infoMessage);
	void   createNamesOfResultFiles(Dictionary *dictionary, MP5Parameters *mp5Parameters);
	STATUS openAsciiDataFile(MP5Parameters *mp5Parameters, char *infoMessage);
	STATUS openResultFiles(MP5Parameters *mp5Parameters, char *infoMessage);
	void   closeFiles(MP5Parameters *mp5Parameters);
	STATUS analyseBinaryDataFile(MP5Parameters *mp5Parameters, char *infoMessage);
	STATUS readDataFileOneTrial(MP5Parameters *mp5Parameters, unsigned short int offsetNumber, char *infoMessage);
	STATUS readDataFileMultiTrial(MP5Parameters *mp5Parameters, char *infoMessage);
	void   printInfoAboutData(Dictionary *dictionary, MP5Parameters *mp5Parameters);
	STATUS writeMagic(MP5Parameters *mp5Parameters, char *infoMessage);
	STATUS writeFileHeader(Dictionary *dictionary, MP5Parameters *mp5Parameters, char *infoMessage);
	STATUS writeCommentsSegment(ConfigFile *configFile, MP5Parameters *mp5Parameters, char *infoMessage);
	STATUS writeSMPResults(Dictionary *dictionary, MP5Parameters *mp5Parameters, unsigned short int offsetNumber, unsigned short int channelNumber, char *infoMessage);
	STATUS writeMMPResults(Dictionary *dictionary, MP5Parameters *mp5Parameters, unsigned short int offsetNumber, char *infoMessage);
	STATUS writeMMPMultiTrialResults(Dictionary *dictionary, MP5Parameters *mp5Parameters, char *infoMessage);
	void returnAmplitudeAndModulusForMMP2DI(MP5Parameters *mp5Parameters, Dictionary *dictionary, Atom *atom, float *amplitude, float *modulus, unsigned int channelNumber);
	void printInformationAboutProgress(MP5Parameters *mp5Parameters, Progress *progress, unsigned int atomsCounter);

#endif
