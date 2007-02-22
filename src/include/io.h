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
 

#ifndef _IO_H_

	#define _IO_H_

	#include<stdlib.h>
	#include<stdio.h>
	#include<string.h>
	#include"types.h"

	void   setDataParameters(DataParameters *dataParameters);
	void   freeDataParameters(DataParameters *dataParameters);
	STATUS testFilesAndDirectories(DataParameters *dataParameters, const ConfigFile *configFile, char *info);
	STATUS testDataParameters(DataParameters *dataParameters, GaborDictionary *gaborDictionary, MP5Parameters *mp5Parameters, char *info);
	STATUS openBinaryDataFile(DataParameters *dataParameters, char *info);
	void createNamesOfResultFiles(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary);
	STATUS openAsciiDataFile(DataParameters *dataParameters, char *info);
	STATUS openResultFiles(DataParameters *dataParameters, char *info);
	void   closeFiles(DataParameters *dataParameters);
	STATUS analyseBinaryDataFile(DataParameters *dataParameters, char *info);
	STATUS analyseAsciiDataFile(DataParameters *dataParameters, char *info);
	void   processRawData(DataParameters *dataParameters);
	STATUS readBinaryData(DataParameters *dataParameters, unsigned short int offsetNumber, char *info);
	STATUS readAsciiData(DataParameters *dataParameters, unsigned short int offsetNumber, char *info);
	STATUS readDataFile(DataParameters *dataParameters, unsigned short int offsetNumber, char *info);
	void   printInfoAboutData(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary);
	STATUS writeHeader(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary, char *info);
	STATUS writeSingleChannelResults(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary, unsigned short int offsetNumber, unsigned short int channelNumber, char *info);
	STATUS writeMultiChannelResults(DataParameters *dataParameters, MP5Parameters *mp5Parameters, GaborDictionary *gaborDictionary, unsigned short int offsetNumber, char *info);
	
#endif
