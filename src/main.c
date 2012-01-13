/***********************************************                                                                                                                                                                     **************************************
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
#include<stdio.h>
#include"atom.h"
#include"cmd.h"
#include"def.h"
#include"dic.h"
#include"io_mp5.h"
#include"matrix.h"
#include"mmp.h"
#include"mp5.h"
#include"smp.h"
#include"stringTools.h"
#include"types.h"

#define PROGRAM_VERSION "0.1"

unsigned char applicationMode = 0x00;

#ifdef __MINGW32__
	#define bzero(ptr,size) memset (ptr, 0, size);
	#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

int main(int argc, char *argv[])
{
	char infoMessage[LENGTH_OF_INFO_MESSAGE];
	bzero((void *)infoMessage,LENGTH_OF_INFO_MESSAGE);

    ConfigFile configFile = {"",
							 NULL,
							 NULL,
							 NULL};

    MP5Parameters  mp5Parameters = {"","","",
                                    NULL,NULL,
                                    0,0,0,0,0,
                                    0.0,
                                    0,NULL,0,NULL,
                                    0,
                                    0.0,
                                    NULL,NULL,
                                    0x0,
                                    0,0,0,0,0,0,0,0,
                                    NULL,NULL,NULL,
                                    ON,
                                    NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                    NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                    0.0,0.0,0.0,0.0,
                                    NULL,NULL,NULL,NULL,
                                    0,0.0,
                                    NULL,NULL,
                                    0x0,0,
                                    NULL,
                                    YES,0x0};

    Dictionary dictionary = {0x0,SCALE_TO_PERIOD_FACTOR,
                             0.0,0.0,
                             0L,0,
                             0.0,0.0,0.0,
                             0,
                             0.0,
                             NULL,NULL,NULL,NULL,NULL,NULL,
                             0,0,0,0,0,0,0,0,0,0,0,0,
                             0x0,0x0,0x0,0x0,
                             NULL,NULL,NULL,NULL,NULL};

    int print_help = argc >= 2 && !strcmp(argv[1], "--help");
    int print_version = argc >= 2 && !strcmp(argv[1], "--version");

    if (print_version) {
	    fprintf(stdout, "mp5 " PROGRAM_VERSION "\n");
	    return 0;
    }

    if(argc>6 || argc==1
			  || (argc==2 && ((strcmp(argv[1],"-g")!=0) && (strcmp(argv[1],"-e")!=0)))
			  || ((strcmp(argv[1],"-g")==0) && (argc>2))
			  || ((strcmp(argv[1],"-e")==0) && (argc>2))
			  || (argc==5 && (((strcmp(argv[1],"-m")!=0) && (strcmp(argv[2],"-a")!=0) && (strcmp(argv[4],"-N")!=0)) ||
				              ((strcmp(argv[1],"-m")!=0) && (strcmp(argv[2],"-N")!=0) && (strcmp(argv[4],"-a")!=0)))) ||
	    print_help)
    {
		fprintf(stderr," \n");
		if (!print_help) {
			fprintf(stderr," ERROR: \n");
			fprintf(stderr," INCORRECT CALL OF mp5Parameters PROGRAM \n");
		}
		fprintf(stderr," THE PROPER USE IS AS FOLLOWS: \n");
		fprintf(stderr," mp5 --help | --version         - print help or version\n");
		fprintf(stderr," mp5 -g                         - generate default config file \n");
		fprintf(stderr," mp5 -e                         - generate xml file with errors codes and translations\n");
		fprintf(stderr," mp5 -t [ name of config file ] - test config file \n");
		fprintf(stderr," mp5 -f [ name of config file ] - process of mp decomposition (user mode)\n");
		fprintf(stderr," mp5 -x [ name of config file ] - process of mp decomposition (server mode)\n");
		fprintf(stderr," \n");
		return print_help ? 0 : 1;
    }

    if(strcmp(argv[1],"-g")==0)
		applicationMode = CONFIG_MODE;
	else if(strcmp(argv[1],"-e")==0)
		applicationMode = ERROR_MODE;
	else if(strcmp(argv[1],"-t")==0)
		applicationMode = TEST_PARAMETERS_MODE ;
    else if(strcmp(argv[1],"-f")==0)
		applicationMode = PROCESS_USER_MODE;
    else if(strcmp(argv[1],"-x")==0)
		applicationMode = PROCESS_SERVER_MODE;
    else
    {
		fprintf(stderr," \n");
		fprintf(stderr," ERROR: \n");
		fprintf(stderr," INCORRECT CALL OF mp5Parameters PROGRAM \n");
		fprintf(stderr," UNKNOW PARAMETER: %s \n",argv[1]);
		fprintf(stderr," THE PROPER USE IS AS FOLLOWS:\n");
		fprintf(stderr," mp5 -g                         - generate default config file\n");
		fprintf(stderr," mp5 -e                         - generate xml file with errors codes and translations\n");
		fprintf(stderr," mp5 -t [ name of config file ] - test config file \n");
		fprintf(stderr," mp5 -f [ name of config file ] - process of mp decomposition (USER MODE)\n");
		fprintf(stderr," mp5 -x [ name of config file ] - process of mp decomposition (SERVER MODE)\n");
		fprintf(stderr," \n");
		return 1;
    }

	if((applicationMode & PROCESS_USER_MODE) || (applicationMode & CONFIG_MODE) || (applicationMode & ERROR_MODE))
	{
		printf("\n");
		printf("  Matching Pursuit V (version 2008-01-08)\n");
		printf("  Department of Biomedical Physics at Warsaw University\n");
		printf("  http://brain.fuw.edu.pl, http://eeg.pl\n");
		printf("  Compiled: %s, %s \n", __DATE__ , __TIME__ );

		#ifdef __VERSION__
			printf("  Compiler: %s \n", __VERSION__);
		#endif

		printf("\n");
	}

	if(applicationMode & CONFIG_MODE)
    {
		FILE *exampleFile = fopen("mp5_default.set","wt");
		char exampleTest[] =
		{"# OBLIGATORY PARAMETERS\n\
\n\
nameOfDataFile         test.dat\n\
nameOfOutputDirectory  ./\n\
writingMode            CREATE\n\
samplingFrequency      512.0\n\
numberOfChannels       3\n\
selectedChannels       1 2-3\n\
numberOfSamplesInEpoch 64\n\
selectedEpochs         1\n\
typeOfDictionary       OCTAVE_STOCH\n\
dilationFactor         1.3 90.00\n\
randomSeed             auto\n\
reinitDictionary       NO_REINIT_AT_ALL\n\
maximalNumberOfIterations 3\n\
energyPercent             99.99999999\n\
MP                        SMP\n\
scaleToPeriodFactor       1.0\n\
pointsPerMicrovolt        2.0\n\
\n\
# ADDITIONAL PARAMETERS\n\
\n\
diracInDictionary         YES\n\
gaussInDictionary         NO\n\
sinCosInDictionary        NO\n\
gaborInDictionary         NO\n\
progressBar               ON"
};

		if(exampleFile == NULL)
		{
		    fprintf(stderr," \n");
		    fprintf(stderr," ERROR: \n");
		    fprintf(stderr," Can't open file: mp5_default.set, for default parameters");
		    fprintf(stderr," \n");
		    return 1;
		}

		fprintf(exampleFile,"%s\n",exampleTest);
		fclose(exampleFile);

		printf(" Defaults parameters are being written to mp5_default.set file\n");

		return 0;
    }
	else if(applicationMode & ERROR_MODE)
	{
		if(generateErrotTextAndCodeGeneratorFile()==ERROR)
		{
		    fprintf(stderr," \n");
		    fprintf(stderr," ERROR: Can't write errors codes and text to file: errorTranslation.xml");
		    fprintf(stderr," \n");
			return 1;
		}
		else
		{
			printf(" Errors texts and codes of errors has been written to file: errorsCodes.xml\n");
			return 0;
		}
	}
    else if((applicationMode & TEST_PARAMETERS_MODE ) || (applicationMode & PROCESS_USER_MODE) || (applicationMode & PROCESS_SERVER_MODE))
    {
		/* alocate memory for components of configFile */

		setConfigFile(&configFile);

		if(strlen(argv[2])>LENGTH_OF_NAME_OF_CONFIG_FILE)
		{
			fprintf(stderr," ERROR: Too long name of config file %s. The maximum length is %hu\n",argv[2],LENGTH_OF_NAME_OF_CONFIG_FILE);
			return 1;
		}

		strcpy(configFile.name,argv[2]);

		/* try to open config file */

		if(openConfigFile(&configFile,infoMessage)==ERROR)
		{
			fprintf(stderr,"%s",infoMessage);
			return 1;
		}

		/* read config file */
		if(readConfigFile(&configFile,infoMessage)==ERROR) goto ERROR_PROCEDURE;		

		/* analyse line and command in config file
		"primary" test of config file and parameters included in this file */
		if(findDataParametersInConfigFile(&configFile,&dictionary,&mp5Parameters,infoMessage)==ERROR)
			goto ERROR_PROCEDURE;

		/* "extended" test of parameters read in config file */
		if(testMP5Parameters(&dictionary,&mp5Parameters,infoMessage)==ERROR)
			goto ERROR_PROCEDURE;
			
		if(applicationMode & TEST_PARAMETERS_MODE)		
		{			
			strcpy(infoMessage," NO ERRORS IN CONFIG FILE\n");
			goto ERROR_PROCEDURE;
		}						

		/* try to open data */
        if(openBinaryDataFile(&mp5Parameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;

		/* "primary" test of data file */
        if(analyseBinaryDataFile(&mp5Parameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;

        /* set some constans values for MP5Parameters and DataParameters */
        setNumberOfAnalysedChannelsAndNumberOfResultsFiles(&mp5Parameters);

        /* check dictionary tape and size */
        analyseDictionarySizeAndType(&dictionary,&mp5Parameters);

        /* create names of results files */
        createNamesOfResultFiles(&dictionary,&mp5Parameters);

		/* check wheter directories and results files are existing */
		if(testFilesAndDirectories(&mp5Parameters,&configFile,infoMessage) == ERROR) goto ERROR_PROCEDURE;

		if(applicationMode & PROCESS_USER_MODE)
			printInfoAboutData(&dictionary,&mp5Parameters);

        if((applicationMode & PROCESS_USER_MODE) || (applicationMode & PROCESS_SERVER_MODE))
        {
    	    /* now allocate memory for components of Dictionary */
    	    allocateDictionary(&dictionary,&mp5Parameters);

			/* now allocate memory for components of MP5Parameters */
			setMP5Parameters(&dictionary,&mp5Parameters);

			/* make sin/cos and exp tables */
			makeSinCosExpTable(&dictionary,&mp5Parameters);

			/* results file */
			if(openResultFiles(&mp5Parameters,infoMessage)==ERROR)
				goto ERROR_PROCEDURE;

			unsigned short int epochNumber   = 0;
			unsigned short int channelNumber = 0;

			/* results file is open, header is written into results file */
			{
				unsigned long int currentPositionInResultsFile = ftell(mp5Parameters.resultsFile);
				fseek(mp5Parameters.resultsFile,0L,SEEK_END);
				unsigned long int sizeOfResultsFile = ftell(mp5Parameters.resultsFile);
				fseek(mp5Parameters.resultsFile,currentPositionInResultsFile,SEEK_SET);

				if((mp5Parameters.writingMode & CREATE_FILE) || (sizeOfResultsFile==0))
				{
					if(writeMagic(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					if(writeCommentsSegment(&configFile,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					if(writeFileHeader(&dictionary,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;
				}
			}

			if(mp5Parameters.MPType & SMP)
			{
				if(applicationMode & PROCESS_SERVER_MODE)
				{
					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);
					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);

					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						for(channelNumber=0;channelNumber<mp5Parameters.numberOfSelectedChannels;channelNumber++)
						{
							if(applicationMode & PROCESS_SERVER_MODE)
							{
								printf("CHANNEL %hu\n",channelNumber);
								fflush(stdout);
							}
							else if(applicationMode & PROCESS_USER_MODE)
							{
								printf("\n --EPOCH--: %d, |CHANNEL|: %d\n\n",mp5Parameters.selectedEpochs[epochNumber],mp5Parameters.selectedChannels[channelNumber]);
								fflush(stdout);
							}

							mp5Parameters.singleChannelSignalTable = *(mp5Parameters.processedDataMatrix + channelNumber);

							firstIterationSMPMMP2(&dictionary,&mp5Parameters);
							nextIterationSMPMMP2(&dictionary,&mp5Parameters);

							if(writeSMPResults(&dictionary,&mp5Parameters,epochNumber,channelNumber,infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							resetDictionary(&dictionary);
						}
					}
				}
				else if(mp5Parameters.reinitDictionary & REINIT_IN_CHANNEL_DOMAIN)
				{
					unsigned short int counter;
					time_t       seedTime;
					long int     seed[mp5Parameters.numberOfSelectedChannels];

					for(counter=0;counter<mp5Parameters.numberOfSelectedChannels;counter++)
						if(dictionary.randomSeed==AUTO_RANDOM_SEED)
							seed[counter] = time(&seedTime);
						else
							seed[counter] = dictionary.randomSeed;

					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						for(channelNumber=0;channelNumber<mp5Parameters.numberOfSelectedChannels;channelNumber++)
						{
							dictionary.randomSeed = seed[channelNumber];

							/* create dicionary */
							reinitDictionary(&dictionary,&mp5Parameters);

							/* test atom's feature, for example find INCORRECT atoms */
//							testAtomFeature(&dictionary);

							if(applicationMode & PROCESS_SERVER_MODE)
							{
								printf("CHANNEL  %hu\n",channelNumber);
								fflush(stdout);
							}

							if(applicationMode & PROCESS_USER_MODE)
							{
								printf("\n |CHANNEL|: %d, --EPOCH--: %d\n\n",mp5Parameters.selectedChannels[channelNumber],mp5Parameters.selectedEpochs[epochNumber]);
								fflush(stdout);
							}

							mp5Parameters.singleChannelSignalTable = *(mp5Parameters.processedDataMatrix + channelNumber);

							firstIterationSMPMMP2(&dictionary,&mp5Parameters);
							nextIterationSMPMMP2(&dictionary,&mp5Parameters);

							if(writeSMPResults(&dictionary,&mp5Parameters,epochNumber,channelNumber,infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							resetDictionary(&dictionary);
						}
					}
				}
				else if(mp5Parameters.reinitDictionary & REINIT_IN_EPOCH_DOMAIN)
				{
					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						for(channelNumber=0;channelNumber<mp5Parameters.numberOfSelectedChannels;channelNumber++)
						{
							if(applicationMode & PROCESS_SERVER_MODE)
							{
								printf("CHANNEL  %hu\n",channelNumber);
								fflush(stdout);
							}

							if(applicationMode & PROCESS_USER_MODE)
							{
								printf("\n |CHANNEL|: %d, --EPOCH--: %d\n\n",mp5Parameters.selectedChannels[channelNumber],mp5Parameters.selectedEpochs[epochNumber]);
								fflush(stdout);
							}

							mp5Parameters.singleChannelSignalTable = *(mp5Parameters.processedDataMatrix + channelNumber);

							firstIterationSMPMMP2(&dictionary,&mp5Parameters);
							nextIterationSMPMMP2(&dictionary,&mp5Parameters);

							if(writeSMPResults(&dictionary,&mp5Parameters,epochNumber,channelNumber,infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							resetDictionary(&dictionary);
						}
					}
				}
				else if(mp5Parameters.reinitDictionary & REINIT_AT_ALL)
				{
					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						for(channelNumber=0;channelNumber<mp5Parameters.numberOfSelectedChannels;channelNumber++)
						{
							if(applicationMode & PROCESS_SERVER_MODE)
							{
								printf("CHANNEL  %hu\n",channelNumber);
								fflush(stdout);
							}

							/* create dicionary */
							reinitDictionary(&dictionary,&mp5Parameters);

							/* test atom's feature, for example find INCORRECT atoms */
//							testAtomFeature(&dictionary);

							if(applicationMode & PROCESS_USER_MODE)
							{
								printf("\n --EPOCH--: %d, |CHANNEL|: %d\n\n",mp5Parameters.selectedEpochs[epochNumber],mp5Parameters.selectedChannels[channelNumber]);
								fflush(stdout);
							}

							mp5Parameters.singleChannelSignalTable = *(mp5Parameters.processedDataMatrix + channelNumber);

							firstIterationSMPMMP2(&dictionary,&mp5Parameters);
							nextIterationSMPMMP2(&dictionary,&mp5Parameters);

							if(writeSMPResults(&dictionary,&mp5Parameters,epochNumber,channelNumber,infoMessage)==ERROR)
								goto ERROR_PROCEDURE;
						}
					}
				}
			}
			// writen by Artur Matysiak
			else if(mp5Parameters.MPType & MMP1)
			{
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);

					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --EPOCH--: %d\n\n",mp5Parameters.selectedEpochs[epochNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						firstIterationMMP(&dictionary,&mp5Parameters);
						nextIterationMMP(&dictionary,&mp5Parameters);

						if(writeMMPResults(&dictionary,&mp5Parameters,epochNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
				if(mp5Parameters.reinitDictionary &  REINIT_IN_EPOCH_DOMAIN)
				{

					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --EPOCH--: %d\n\n",mp5Parameters.selectedEpochs[epochNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						firstIterationMMP(&dictionary,&mp5Parameters);
						nextIterationMMP(&dictionary,&mp5Parameters);

						if(writeMMPResults(&dictionary,&mp5Parameters,epochNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
		    }
			//Artur Matysiak end
			else if((mp5Parameters.MPType & MMP12) || (mp5Parameters.MPType & MMP21))
		    {
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
						fflush(stdout);
				}


				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
	                // testAtomFeature(&dictionary);
					// if(applicationMode & PROCESS_SERVER_MODE)
					// {
					// printf("EPOCH  %hu\n",epochNumber);
					// fflush(stdout);
					// }

					if(readDataFileMultiTrial(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					// if(applicationMode & PROCESS_USER_MODE)
					// {
					// printf("\n --EPOCH--: %d\n\n",mp5Parameters.chosenEpochs[epochNumber]);
					// fflush(stdout);
					// }


					mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

					firstIterationMultiChannelMultiTrial(&dictionary,&mp5Parameters);
					nextIterationMultiChannelMultiTrial(&dictionary,&mp5Parameters);

					if(readDataFileMultiTrial(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

					if(writeMMPMultiTrialResults(&dictionary,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

				}
			}
			else if(mp5Parameters.MPType & MMP11)
		    {
				if(applicationMode & PROCESS_SERVER_MODE)
				{
					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);

					// if(applicationMode & PROCESS_SERVER_MODE)
					// {
						// printf("EPOCH  %hu\n",epochNumber);
						// fflush(stdout);
					// }

					if(readDataFileMultiTrial(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					// if(applicationMode & PROCESS_USER_MODE)
					// {
						// printf("\n --EPOCH--: %d\n\n",mp5Parameters.chosenEpochs[epochNumber]);
						// fflush(stdout);
					// }

					mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

					firstIterationMMP(&dictionary,&mp5Parameters);
					nextIterationMMP(&dictionary,&mp5Parameters);

					if(writeMMPMultiTrialResults(&dictionary,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

				}
			}
		    else if(mp5Parameters.MPType & MMP2)
		    {
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);

					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --EPOCH--: %d\n\n",mp5Parameters.selectedEpochs[epochNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;
						countMeanSignalOverChannelsInOneEpoch(&mp5Parameters);
						mp5Parameters.singleChannelSignalTable = *(mp5Parameters.meanSignalTable);

						firstIterationSMPMMP2(&dictionary,&mp5Parameters);
						nextIterationSMPMMP2(&dictionary,&mp5Parameters);

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
						    goto ERROR_PROCEDURE;

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						if(writeMMPResults(&dictionary,&mp5Parameters,epochNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
				if(mp5Parameters.reinitDictionary & REINIT_IN_EPOCH_DOMAIN)
				{

					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{

						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --EPOCH--: %d\n\n",mp5Parameters.selectedEpochs[epochNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						countMeanSignalOverChannelsInOneEpoch(&mp5Parameters);
						mp5Parameters.singleChannelSignalTable = *(mp5Parameters.meanSignalTable);

						firstIterationSMPMMP2(&dictionary,&mp5Parameters);
						nextIterationSMPMMP2(&dictionary,&mp5Parameters);

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
						    goto ERROR_PROCEDURE;

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						if(writeMMPResults(&dictionary,&mp5Parameters,epochNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
		    }
		    else if(mp5Parameters.MPType & MMP22)
		    {				
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);

					// if(applicationMode & PROCESS_SERVER_MODE)
					// {
						// printf("EPOCH  %hu\n",epochNumber);
						// fflush(stdout);
					// }

					if(readDataFileMultiTrial(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					// if(applicationMode & PROCESS_USER_MODE)
					// {
						// printf("\n --EPOCH--: %d\n\n",mp5Parameters.chosenEpochs[epochNumber]);
						// fflush(stdout);
					// }

					mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;
					countMeanSignalOrResidumOverChannelsAndTrials(&mp5Parameters);

					mp5Parameters.singleChannelSignalTable = *(mp5Parameters.meanSignalTable);
					firstIterationSMPMMP2(&dictionary,&mp5Parameters);
					nextIterationSMPMMP2(&dictionary,&mp5Parameters);

					if(readDataFileMultiTrial(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

					if(writeMMPMultiTrialResults(&dictionary,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;
				}
			}
		    else if(mp5Parameters.MPType & MMP3)
		    {
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);
					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);

					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --EPOCH--: %d\n\n",mp5Parameters.selectedEpochs[epochNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						firstIterationMMP(&dictionary,&mp5Parameters);
						nextIterationMMP(&dictionary,&mp5Parameters);

						if(writeMMPResults(&dictionary,&mp5Parameters,epochNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
				if(mp5Parameters.reinitDictionary &  REINIT_IN_EPOCH_DOMAIN)
				{

					for(epochNumber=0;epochNumber<mp5Parameters.numberOfSelectedEpochs;epochNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("EPOCH  %hu\n",epochNumber);
							fflush(stdout);
						}

						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						if(readDataFileOneTrial(&mp5Parameters,mp5Parameters.selectedEpochs[epochNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --EPOCH--: %d\n\n",mp5Parameters.selectedEpochs[epochNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						firstIterationMMP(&dictionary,&mp5Parameters);
						nextIterationMMP(&dictionary,&mp5Parameters);

						if(writeMMPResults(&dictionary,&mp5Parameters,epochNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
			}
			else if((mp5Parameters.MPType & MMP32) || (mp5Parameters.MPType & MMP23))
		    {
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}


				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);
					fflush(stdout);

					/* test atom's feature, for example find INCORRECT atoms */
	                // testAtomFeature(&dictionary);
					// if(applicationMode & PROCESS_SERVER_MODE)
					// {
					// printf("EPOCH  %hu\n",epochNumber);
					// fflush(stdout);
					// }

					if(readDataFileMultiTrial(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					// if(applicationMode & PROCESS_USER_MODE)
					// {
					// printf("\n --EPOCH--: %d\n\n",mp5Parameters.chosenEpochs[epochNumber]);
					// fflush(stdout);
					// }

					mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

					firstIterationMultiChannelMultiTrial(&dictionary,&mp5Parameters);
					nextIterationMultiChannelMultiTrial(&dictionary,&mp5Parameters);

					if(readDataFileMultiTrial(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

					if(writeMMPMultiTrialResults(&dictionary,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

				}
			}
			else if(mp5Parameters.MPType & MMP33)
		    {
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfSelectedEpochs,
															mp5Parameters.numberOfSelectedChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}


				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);

					// if(applicationMode & PROCESS_SERVER_MODE)
					// {
						// printf("EPOCH  %hu\n",epochNumber);
						// fflush(stdout);
					// }

					if(readDataFileMultiTrial(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					// if(applicationMode & PROCESS_USER_MODE)
					// {
						// printf("\n --EPOCH--: %d\n\n",mp5Parameters.chosenEpochs[epochNumber]);
						// fflush(stdout);
					// }

					mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

					firstIterationMMP(&dictionary,&mp5Parameters);
					nextIterationMMP(&dictionary,&mp5Parameters);

					if(writeMMPMultiTrialResults(&dictionary,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

				}
			}
		}
    }

    freeConfigFile(&configFile);
    closeFiles(&mp5Parameters);
    freeMP5Parameters(&dictionary,&mp5Parameters);
    freeDictionary(&dictionary);
	printf(" END\n");
	fflush(stdout);

    return 0;

    ERROR_PROCEDURE:
	fprintf(stderr," \n");
	fprintf(stderr," %s",infoMessage);
	fprintf(stderr,"\n");
	freeConfigFile(&configFile);
	closeFiles(&mp5Parameters);
	freeMP5Parameters(&dictionary,&mp5Parameters);
    freeDictionary(&dictionary);
	fflush(stdout);

    return 1;
}
