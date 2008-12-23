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

#ifndef _DEF_H_

	#define _DEF_H_
    
	#define LENGTH_OF_STRING 64
	#define LENGTH_OF_LINE   4096

	#define LENGTH_OF_NAME_OF_CONFIG_FILE        256
	#define LENGTH_OF_NAME_OF_DATA_FILE          256
	#define LENGTH_OF_NAME_OF_RESULTS_FILE       256
	#define LENGTH_OF_EXTENSION_OF_RESULTS_FILE  256
	#define LENGTH_OF_OUTPUT_DIRECTORY           256
	#define LENGTH_OF_INFO_MESSAGE               4096
	#define LENGTH_OF_TMP_STRING		         256
	#define SIZE_OF_OUTPUT_BUFFOR		         32

	#define NUMBER_OF_ALL_COMMANDS 				 32
	#define NUMBER_OF_PERMAMENT_COMMANDS     	 24
	#define LENGHT_OF_MAGIC_IN_RESULT_FILE		 6

	typedef enum {ERROR,SUCCESS} STATUS;
	typedef enum {FALSE,TRUE}    BOOLEAN;
	typedef enum {NO,YES}        CHOICE;
	typedef enum {OFF,ON}        SWITCH;

	// FFT MODE
	
	#define OFF     0x00
	#define FFT1    0x01 
	#define FFT2	0x02 
	
	#define M_2PI 6.28318530717958623200  /* 2pi */
	#define EXP   2.71828182845905        /* e */

	#define NOT_ALLOCATED 0x01
	#define ALLOCATED     0x02
	
	#define FORMAT_UNKNOW 0x00
	#define FORMAT_ASCII  0x01
	#define FORMAT_SHORT  0x02
	#define FORMAT_FLOAT  0x04

	#define OCTAVE_FIXED 0x01
	#define OCTAVE_STOCH 0x02

	#define SMP  0x01 /* bit 0 */
	#define MMP1 0x02 /* bit 1 */
	#define MMP2 0x04 /* bit 2 */
	#define MMP3 0x08 /* bit 3 */
	#define MMP4 0x10 /* bit 4 */
	
	#define GAUSS_NON_GAUSS 0x01
	#define GAUSS_GAUSS     0x02

//	#define EPS_BAS 1E-16                /* BASIC EPS */
	#define EPS_DET 1E-16	             /* EPS FOR DETRMINANT OF MATRIX */

	/* Dirac's Delta Function  */
	#define DIRACDELTA     0x01 /* bit 1 */

	/* Gauss Function  */
	#define GAUSSFUNCTION  0x02 /* bit 2 */

	/* SIN/COS WAVE  */
	#define SINCOSWAVE     0x04 /* bit 4 */
	
	/* GABOR WAVE */
	#define GABORWAVE      0x08 /* bit 3 */

	/* Wave, which should be removed form dictionary for some reasons */
	#define INCORRECTGABOR 0x10 /* bit 0 */

	/* position of atom>offsetDimension/2  */
	#define LEFT_SIDE_POSITION_IN_OFFSET 0x20 /* bit 5 */

	/* position of atom>range/2  */
	#define LEFT_SIDE_POSITION_IN_RANGE 0x40 /* bit 6 */

	/* Atom was chosen during mp5 processing */
	#define GABOR_WAS_HIT       0x80 /* bit 7 */

	#define LEFT_SIDE_POSITION  0x01 /* bit 1 */
	#define RIGHT_SIDE_POSITION 0x02 /* bit 2 */

	#define FULL_RANGE  0x01 /* bit 1 */
	#define LIMIT_RANGE 0x02 /* bit 2 */

	#define FIRST_ITERATION     0x01 /* bit 1 */
	#define NEXT_ITERATION      0x02 /* bit 2 */
	#define MMP1_NEXT_ITERATION 0x04 /* bit 3 */
	#define MMP1_IGNORE_RS_RC   0x08 /* bit 4 */

	#define FIRST_PART_OF_TABLES_IN_DICTIONARY_ALLOCATED  0x01 /* bit 0 */
	#define SECOND_PART_OF_TABLES_IN_DICTIONARY_ALLOCATED 0x02 /* bit 1 */

	#define NO_REINIT_AT_ALL         0x01 /* bit 0 */
	#define REINIT_IN_CHANNEL_DOMAIN 0x02 /* bit 1 */
	#define REINIT_IN_OFFSET_DOMAIN  0x04 /* bit 2 */
	#define REINIT_AT_ALL            0x08 /* bit 3 */

	#define CREATE_FILE 0x01 /* bit 0 */
	#define APPEND_FILE 0x02 /* bit 1 */
	
	#define SIGNAL 0x01
	#define BOOK   0x02

	#define NUMBER_OF_STEPS_IN_TOOLBAR 10

	#define VERBOSE_PRINT_DICTIONARY   0x01
	#define VERBOSE_PRINT_FITED_ATOMS  0x02
	#define VERBOSE_PRINT_PROCESSBAR   0x04
	
	#define CONFIG_MODE 		   0x01
	#define ERROR_MODE 		       0x02
	#define MEMORY_USAGE_TEST_MODE 0x04
	#define TEST_PARAMETERS_MODE   0x08
	#define PROCESS_USER_MODE      0x10
	#define PROCESS_SERVER_MODE    0x20
	
	#define AUTO_RANDOM_SEED -1

	#define NUMBER_OF_SIGMAS  	  4
	#define LOG_EPS_DOT_PRODUCT   (-((NUMBER_OF_SIGMAS)*(NUMBER_OF_SIGMAS))) // or -log(EPS), where EPS is an accurance of mp5 dot product

	
	#define strtokMask "{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}"

	#define NUMBER_OF_ERRORS			      32
	
	#define CAN_NOT_OPEN_CONFIG_FILE          1
	#define LINE_DOES_NOT_INCLUDE_COMMAND     2
	#define LINE_IS_BROKEN_INCORRECTLY        3
	#define LINE_IS_TOO_LONG                  4
	#define COMMAND_NOT_FOUND                 5
	#define INCORRECT_NUMBER_OF_ARGUMENTS     6
	#define INCORRECT_TYPE_OF_ARGUMENT        7
	#define INCORRECT_SYNTAX_OF_ARGUMENT      8

	#define ARGUMENT_SHOUDL_BE_INTEGER_GREATER_OR_EQUAL_TO_ZERO  9
	#define ARGUMENT_SHOUDL_BE_INTEGER_GREATER_TO_ZERO          10
	#define ARGUMENT_SHOUDL_BE_FLOAT_GREATER_OR_EQUAL_TO_ZERO   11
	#define ARGUMENT_SHOUDL_BE_FLOAT_GREATER_TO_ZERO            12
	#define ARGUMENT_SHOUDL_BE_FLOAT_GREATER_TO_ONE             13

	
	#define CAN_NOT_OPEN_DIRECTORY             14
	#define OVERWRITE_RESULTS_ALARM            15
	#define BAD_NUMBER_OF_CHOSEN_CHANNELS      16
	#define BAD_NUMBER_OF_CHOSEN_OFFSETS       17
	#define BAD_PERIOD_DENSITY_EQUAL           18
	#define BAD_PERIOD_DENSITY_GREATER         19
	#define BAD_PERIOD_DENSITY_GREATER_FFT     20
	#define BAD_ENERGY_PERCENT                 21
	#define BAD_REINIT_ALL                     22
	#define BAD_REINIT_MMP                     23
	#define BAD_MP_ALGORITHM                   24
	
	
	#define CAN_NOT_OPEN_DATA_FILE             25
	#define CAN_NOT_OPEN_RESULTS_FILE          26
	#define BAD_NUMBER_OF_SAMPLES_PER_CHANNEL  27
	#define BAD_NUMBER_OF_SAMPLES              28
	#define CAN_NOT_READ_SAMPLE				   29
	#define CAN_NOT_WRITE_HEADER               30
	#define CAN_NOT_WRITE_RESULTS              31
	#define CAN_NOT_READ_OFFSET_IN_BINARY_FILE 32
	#define CAN_NOT_READ_OFFSET_IN_ASCII_FILE  33
	
	#define WSAVE_DIMENSION 				   1024
	
	#define FLOATING 1
	#define FULL     2
	
#endif


