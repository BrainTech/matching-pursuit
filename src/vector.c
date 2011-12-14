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

#include<stdlib.h>
#include<strings.h>
#include"vector.h"
#include<stdio.h>

#ifdef __MINGW32__
    #define bzero(ptr,size) memset (ptr, 0, size);
#endif
/* double vector operation */

double *dVectorAllocate(unsigned int numberOfElements)
{
	return (double *)malloc(numberOfElements*sizeof(double));
}

void dVectorFree(double *vector)
{
	free(vector);
}

void dSetVectorZero(double *vector, unsigned int numberOfElements)
{
	bzero((void *)vector,numberOfElements*sizeof(double));
}

/* float vector operation */

float *fVectorAllocate(unsigned int numberOfElements)
{
	return (float *)malloc(numberOfElements*sizeof(float));
}

void fVectorFree(float *vector)
{
	free(vector);
}

void fSetVectorZero(float *vector, unsigned int numberOfElements)
{
	bzero((void *)vector,numberOfElements*sizeof(float));
}

/* double complex  vector operation */

double *dComplexVectorAllocate(unsigned int numberOfElements)
{
	return (double *)malloc(2*numberOfElements*sizeof(double));
}

void dComplexVectorFree(double *vector)
{
	free(vector);
}

void dComplexSetVectorZero(double *vector, unsigned int numberOfElements)
{
	bzero((void *)vector,2*numberOfElements*sizeof(double));
}

/* float complex vector operation */

float *fComplexVectorAllocate(unsigned int numberOfElements)
{
	return (float *)malloc(2*numberOfElements*sizeof(float));
}

void fComplexVectorFree(float *vector)
{
	free(vector);
}

void fComplexSetVectorZero(float *vector, unsigned int numberOfElements)
{
	bzero((void *)vector,2*numberOfElements*sizeof(float));
}

/* integer vector operation */

int *iVectorAllocate(unsigned int numberOfElements)
{
	return (int *)malloc(numberOfElements*sizeof(int));
}

void iVectorFree(int *vector)
{
	free(vector);
}

void iSetVectorZero(int *vector, unsigned int numberOfElements)
{
	bzero((int *)vector,numberOfElements*sizeof(int));
}

/* unsigned integer vector operation */

unsigned int *uiVectorAllocate(unsigned int numberOfElements)
{
	return (unsigned int *)malloc(numberOfElements*sizeof(unsigned int));
}

void uiVectorFree(unsigned int *vector)
{
	free(vector);
}

void uiSetVectorZero(unsigned int *vector, unsigned int numberOfElements)
{
	bzero((unsigned int *)vector,numberOfElements*sizeof(unsigned int));
}

/* short integer vector operation */

short int *siVectorAllocate(unsigned int numberOfElements)
{
	return (short int *)malloc(numberOfElements*sizeof(short int));
}

void siVectorFree(short int *vector)
{
	free(vector);
}

void siSetVectorZero(short int *vector, unsigned int numberOfElements)
{
	bzero((short int *)vector,numberOfElements*sizeof(short int));
}

/* unsigned short integer vector operation */

unsigned short int *usiVectorAllocate(unsigned short int numberOfElements)
{
	return (unsigned short int *)malloc(numberOfElements*sizeof(unsigned short int));
}

void usiVectorFree(unsigned short int *vector)
{
	free(vector);
}

void usiSetVectorZero(unsigned short int *vector, unsigned int numberOfElements)
{
	bzero((unsigned short int *)vector,numberOfElements*sizeof(unsigned short int));
}

/* char vector operation */

char *charVectorAllocate(unsigned int numberOfElements)
{
	return (char *)malloc(numberOfElements*sizeof(char));
}

void charVectorFree(char *vector)
{
	free(vector);
}

void charSetVectorZero(char *vector, unsigned int numberOfElements)
{
  unsigned short int element;

	for(element=0;element<numberOfElements;element++)
		*(vector + element) = '\0';
}

int dxypz(unsigned int n, const double *x, int incx, const double *y, int incy, double *z, int incz)
{

/* the procedure performs fast operation:
	z(i) = x(i)*y(i)
*/
	int i, ix, iy, iz, m;

	if((incx==1) && (incy==1) && (incz==1))
	{
		m = n%4;

		if(m!=0)
			for(i=0;i<m;i++)
			{
					z[i] = x[i]*y[i];
			}

		if(n<4)
			return 0;

		for(i=m;i<n;i=i+4)
		{
			z[i] = x[i]*y[i];
			z[i+1] = x[i+1]*y[i+1];
			z[i+2] = x[i+2]*y[i+2];
			z[i+3] = x[i+3]*y[i+3];
		}
		return 0;
	}
	else
	{
		ix = 0;
		iy = 0;
		iz = 0;

		if(incx<0)
			ix = (-n+1)*incx;

		if(incy<0)
			iy = (-n+1)*incy;

		if(incz<0)
			iz = (-n+1)*incz;

		for(i=0;i<n;i++)
		{
			z[iz] = x[ix]*y[iy];

			ix = ix + incx;
			iy = iy + incy;
			iz = iz + incz;
		}

		return 0;
	}
}

int daxypz(unsigned int n, double alpha, const double *x, int incx, double *y, int incy, double *z, int incz)
{

/* the procedure performs fast operation:
	z(i) = alpha*x(i)*y(i)
*/
	double temp;
	int i, ix, iy, iz, m;

	if((incx==1) && (incy==1) && (incz==1))
	{
		m = n%4;

		if(m!=0)
		{
			for(i=0;i<m;i++)
			{
				temp = x[i]*y[i];
				z[i] = alpha*temp;
			}
		}

		if(n<4)
			return 0;

		for(i=m;i<n;i=i+4)
		{
			temp = x[i]*y[i];
			z[i] = alpha*temp ;

			temp   = x[i+1]*y[i+1];
			z[i+1] = alpha*temp;

			temp   = x[i+2]*y[i+2];
			z[i+2] = alpha*temp;

			temp   = x[i+3]*y[i+3];
			z[i+3] = alpha*temp;
		}
		return 0;
	}
	else
	{
		ix = 0;
		iy = 0;
		iz = 0;

		if(incx<0)
			ix = (-n+1)*incx;

		if(incy<0)
			iy = (-n+1)*incy;

		if(incz<0)
			iz = (-n+1)*incz;

		for(i=0;i<n;i++)
		{
			temp = x[ix]*y[iy];
			z[iz] = alpha*temp;

			ix = ix + incx;
			iy = iy + incy;
			iz = iz + incz;
		}

		return 0;
	}
}

int dxypbz(unsigned int n, const double *x, int incx, double *y, int incy, double beta, double *z, int incz)
{

/* the procedure performs fast operation:
	z(i) = x(i)*y(i) +  beta*z(i)
*/
	int i, ix, iy, iz, m;

	if((incx==1) && (incy==1) && (incz==1))
	{
		m = n%4;

		if(m!=0)
		{
			for(i=0;i<m;i++)
			{
				z[i] = beta*z[i];
				z[i] = x[i]*y[i] + z[i];
			}
		}

		if(n<4)
			return 0;

		for(i=m;i<n;i=i+4)
		{
			z[i] = beta*z[i];
			z[i] = x[i]*y[i] + z[i];

			z[i+1] = beta*z[i+1];
			z[i+1] = x[i+1]*y[i+1] + z[i+1];

			z[i+2] = beta*z[i+2];
			z[i+2] = x[i+2]*y[i+2] + z[i+2];

			z[i+3] = beta*z[i+3];
			z[i+3] = x[i+3]*y[i+3] + z[i+3];
		}
		return 0;
	}
	else
	{
		ix = 0;
		iy = 0;
		iz = 0;

		if(incx<0)
			ix = (-n+1)*incx;

		if(incy<0)
			iy = (-n+1)*incy;

		if(incz<0)
			iz = (-n+1)*incz;

		for(i=0;i<n;i++)
		{
			z[iz] = beta*z[iy];
			z[iz] = x[ix]*y[iy] + z[iz];

			ix = ix + incx;
			iy = iy + incy;
			iz = iz + incz;
		}

		return 0;
	}
}

int daxypbz(unsigned int n, double alpha, double *x, int incx, const double *y, int incy, double beta, double *z, int incz)
{

/* the procedure performs fast operation:
	z(i) = alpha*x(i)*y(i) +  beta*z(i)
*/
	double temp;
	int i, ix, iy, iz, m;

	if((incx==1) && (incy==1) && (incz==1))
	{
		m = n%4;

		if(m!=0)
		{
			for(i=0;i<m;i++)
			{
				z[i] = beta*z[i];
				temp = x[i]*y[i];
				z[i] = alpha*temp + z[i];
			}
		}

		if(n<4)
			return 0;

		for(i=m;i<n;i=i+4)
		{
			z[i] = beta*z[i];
			temp = x[i]*y[i];
			z[i] = alpha*temp + z[i];

			z[i+1] = beta*z[i+1];
			temp   = x[i+1]*y[i+1];
			z[i+1] = alpha*temp + z[i+1];

			z[i+2] = beta*z[i+2];
			temp   = x[i+2]*y[i+2];
			z[i+2] = alpha*temp + z[i+2];

			z[i+3] = beta*z[i+3];
			temp   = x[i+3]*y[i+3];
			z[i+3] = alpha*temp + z[i+3];
		}

		return 0;
	}
	else
	{
		ix = 0;
		iy = 0;
		iz = 0;

		if(incx<0)
			ix = (-n+1)*incx;

		if(incy<0)
			iy = (-n+1)*incy;

		if(incz<0)
			iz = (-n+1)*incz;

		for(i=0;i<n;i++)
		{
			z[iz] = beta*z[iz];
			temp = x[ix]*y[iy];
			z[iz] = alpha*temp + z[iz];

			ix = ix + incx;
			iy = iy + incy;
			iz = iz + incz;
		}

		return 0;
	}
}

double ddot(unsigned int n, const double *x, const int incx, const double *y, int incy)
{
/* the procedure performs fast operation:
	dotProduct = x*y ^{T}
*/

	double temp;
	int i, ix, iy, m;

	temp = 0.0;

	if((incx==1) && (incy==1))
	{
		m = n%5;

		if(m!=0)
		{
			for(i=0;i<m;i++)
				temp = temp + x[i]*y[i];
		}

		if(n<5)
			return temp;

		for(i=m;i<n;i=i+5)
			temp = temp + x[i]*y[i] + x[i+1]*y[i+1] +
						  x[i+2]*y[i+2] + x[i+3]*y[i+3] + x[i+4]*y[i+4];

		return temp;
	}
	else
	{
		ix = 0;
		iy = 0;

		if(incx<0)
			ix = (-n+1)*incx;

		if(incy<0)
			iy = (-n+1)*incy;

		for(i=0;i<n;i++)
		{
			temp = temp + x[ix]*y[iy];

			ix = ix + incx;
			iy = iy + incy;
		}

		return temp;
	}
}

int dcopy(unsigned int n, const double *x, int incx, double *y, int incy)
{
/* the procedure performs fast operation:
	y(i) = x(i)
*/
	int i, ix, iy, m;

	if((incx==1) && (incy==1))
	{
		m = n%7;

		if(m!=0)
		{
			for(i=0;i<m;i++)
				y[i] = x[i];
		}

		if(n<7)
			return 0;

		for(i=m;i<n;i=i+7)
		{
			y[i] = x[i];
			y[i + 1] = x[i + 1];
			y[i + 2] = x[i + 2];
			y[i + 3] = x[i + 3];
			y[i + 4] = x[i + 4];
			y[i + 5] = x[i + 5];
			y[i + 6] = x[i + 6];
		}
		return 0;
	}
	else
	{
		ix = 0;
		iy = 0;

		if(incx<0)
			ix = (-n+1)*incx;

		if(incy<0)
			iy = (-n+1)*incy;

		for(i=0;i<n;i++)
		{
			y[iy] = x[ix];
			ix = ix + incx;
			iy = iy + incy;
		}

		return 0;
	}

	return 0;
}

double dasum(unsigned int n, const double *x, int incx)
{
/* the procedure performs fast operation:
	a = sum(x)
*/
	double temp = 0.0;
	int i, m, ix;

	if(incx!=1)
	{
		ix = 0;

		if(incx<0)
			ix = (-n+1)*incx;

		for(i=0;i<n;i++)
		{
			temp = temp + x[ix];
			ix = ix + incx;
		}
		return temp;
	}
	else
	{
		m = n%6;

		if(m!=0)
		{
			for(i=0;i<m;i++)
				temp = temp + x[i];
		}

		if(n<6)
			return temp;

		for(i=m;i<n;i=i+6)
		{
			temp = temp + x[i] +  x[i + 1] + x[i+2]+ + x[i+3] + x[i+4] + x[i+5];
		}

		return temp;
	}
}

