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

#ifndef _MATRIX_H_
    
    #define _MATRIX_H_
    
	/* double matrix operation */
    
	double **dMatrixAllocate(unsigned int numberOfRows, unsigned int numberOfColumns);
	void dMatrixFree(double **matrix, unsigned int numberOfRows);
	void dSetMatrixZero(double **matrix, unsigned int numberOfRows, unsigned int numberOfColumns);
        double **dVariableMatrixAllocate(unsigned int numberOfRows, unsigned int *tableOfNumbersOfColumns);
        void dVariableMatrixFree(double **matrix, unsigned int numberOfRows);
        void dVariableSetMatrixZero(double **matrix, unsigned int numberOfRows, unsigned int *tableOfNumbersOfColumns);

	/* float matrix operation */

	float **fMatrixAllocate(unsigned int numberOfRows, unsigned int numberOfColumns);
	void fMatrixFree(float **matrix, unsigned int numberOfRows);
	void fSetMatrixZero(float **matrix, unsigned int numberOfRows, unsigned int numberOfColumns);

#endif

