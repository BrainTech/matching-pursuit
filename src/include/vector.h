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


#ifndef _VEC_MAT_H_

	#define _VEC_MAT_H_

	/* double vector operation */

	double *dVectorAllocate(unsigned int numberOfElements);
	void    dVectorFree(double *vector);
	void    dSetVectorZero(double *vector, unsigned int numberOfElements);

	/* float vector operation */

	float *fVectorAllocate(unsigned int numberOfElements);
	void   fVectorFree(float *vector);
	void   fSetVectorZero(float *vector, unsigned int numberOfElements);

	/* double complex vector operation */

	double *dComplexVectorAllocate(unsigned int numberOfElements);
	void   dComplexVectorFree(double *vector);
	void   dComplexSetVectorZero(double *vector, unsigned int numberOfElements);

	/* float complex vector operation */

	float *fComplexVectorAllocate(unsigned int numberOfElements);
	void   fComplexVectorFree(float *vector);
	void   fComplexSetVectorZero(float *vector, unsigned int numberOfElements);

	/* integer vector operation */

	int *iVectorAllocate(unsigned int numberOfElements);
	void iVectorFree(int *vector);
	void iSetVectorZero(int *vector, unsigned int numberOfElements);


	/* unsigned integer vector operation */

	unsigned int *uiVectorAllocate(unsigned int numberOfElements);
	void uiVectorFree(unsigned int *vector);
	void uiSetVectorZero(unsigned int *vector, unsigned int numberOfElements);

	/* short integer vector operation */

	short int *siVectorAllocate(unsigned int numberOfElements);
	void siVectorFree(short int *vector);
	void siSetVectorZero(short int *vector, unsigned int numberOfElements);

	/* unsigned short integer vector operation */

	unsigned short int *usiVectorAllocate(unsigned short int numberOfElements);
	void usiVectorFree(unsigned short int *vector);
	void usiSetVectorZero(unsigned short int *vector, unsigned int numberOfElements);

	/* unsigned short integer vector operation */

	char *charVectorAllocate(unsigned int numberOfElements);
	void charVectorFree(char *vector);
	void charSetVectorZero(char *vector, unsigned int numberOfElements);

	int dxypz(unsigned int n, const double *x, int incx, const double *y, int incy, double *z, int incz);
	int daxypz(unsigned int n, double alpha, const double *x, int incx, double *y, int incy, double *z, int incz);
	int dxypbz(unsigned int n, const double *x, int incx, double *y, int incy, double beta, double *z, int incz);
	int daxypbz(unsigned int n, double alpha, double *x, int incx, const double *y, int incy, double beta, double *z, int incz);
	double ddot(unsigned int n, const double *x, const int incx, const double *y, int incy);
	int dcopy(unsigned int n, const double *x, int incx, double *y, int incy);
	double dasum(unsigned int n, const double *x, int inc);

#endif
