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
	#define LENGTH_OF_TMP_STRING		     256

	#define NUMBER_OF_COMMANDS 23

	typedef enum {ERROR,SUCCESS} STATUS;
	typedef enum {FALSE,TRUE} BOOLEAN;

	#define M_2PI 6.28318530717958623200  /* 2pi */

	#define EPS_DOUBLE 1E-16

	#define GENERATE 0x01 /* bit 0 */
	#define TEST     0x02 /* bit 1 */
	#define PROCESS  0x04 /* bit 2 */

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

	#define EPS_BAS 1E-16                /* BASIC EPS */
	#define EPS_DET 1E-16	             /* EPS FOR DETRMINANT OF MATRIX */

	/* Wave, which should be removed form dictionary for some reasons */
	#define INCORRECTGABOR 0x01 /* bit 0 */

	/* Dirac's Delta Function */
	#define DIRACDELTA     0x02 /* bit 1*/

	/* GABOR WAVE*/
	#define GABORWAVE      0x04 /* bit 2*/

	/* SIN/COS WAVE */
	#define FFTWAVE        0x08 /* bit 3 */

	/* position of gabor>dimOffset/2 */
	#define LEFT_SIDE_POSITION_IN_OFFSET 0x10 /* bit 4 */

	/* position of gabor>range/2 */
	#define LEFT_SIDE_POSITION_IN_RANGE 0x20 /* bit 5 */

	/* Gabor was chosen during mp5 processing */
	#define GABOR_WAS_HIT       0x40 /* bit 6 */

	#define LEFT_SIDE_POSITION  0x01 /* bit 1 */
	#define RIGHT_SIDE_POSITION 0x02 /* bit 2 */

	#define FULL_RANGE  0x01 /* bit 1 */
	#define LIMIT_RANGE 0x02 /* bit 2 */

	#define FIRST_ITERATION 0x01 /* bit 1 */
	#define NEXT_ITERATION  0x02 /* bit 2 */

	#define NAMES_OF_RESULT_FILES_ALLOCATED 0x01 /* bit 0 */
	#define RESULT_FILES_ALLOCATED          0x02 /* bit 1 */
	#define CHOSEN_CHANNELS_ALLOCATED       0x04 /* bit 3 */
	#define CHOSEN_OFFSETS_ALLOCATED        0x08 /* bit 4 */
	#define RAW_DATA_MATRIX_ALLOCATED       0x10 /* bit 5 */
	#define PROCESSED_DATA_MATRIX_ALLOCATED 0x20 /* bit 6 */

	#define FIRST_PART_OF_TABLES_IN_DICTIONARY_ALLOCATED  0x01 /* bit 0 */
	#define SECOND_PART_OF_TABLES_IN_DICTIONARY_ALLOCATED 0x02 /* bit 1 */

	#define NO_REINIT_AT_ALL         0x01 /* bit 0 */
	#define REINIT_IN_CHANNEL_DOMAIN 0x02 /* bit 1 */
	#define REINIT_IN_OFFSET_DOMAIN  0x04 /* bit 2 */
	#define REINIT_AT_ALL            0x08 /* bit 3 */

	#define CREATE_FILE 0x01 /* bit 0 */
	#define APPEND_FILE 0x02 /* bit 1 */

	#define NUMBER_OF_STEPS_IN_TOOLBAR 10

	#define VERBOSE_PRINT_DICTIONARY   0x01
	#define VERBOSE_PRINT_FITED_GABORS 0x02
	#define VERBOSE_PRINT_PROCESSBAR   0x04

#endif


