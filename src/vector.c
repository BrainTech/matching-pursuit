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

#include<stdlib.h>
#include<strings.h>
#include"include/vector.h"

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
