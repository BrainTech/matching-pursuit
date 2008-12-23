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


#include<stdio.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include"def.h"
#include"tools.h"

#ifdef INLINE

	inline void sincos(double __x, double *__sinx, double *__cosx)
	{
		register double __cosr,__sinr;
		__asm __volatile__
		(
			"fsincos\n\t"
			"fnstsw    %%ax\n\t"
			"testl     $0x400, %%eax\n\t"
			"jz        1f\n\t"
			"fldpi\n\t"
			"fadd      %%st(0)\n\t"
			"fxch      %%st(1)\n\t"
			"2: fprem1\n\t"
			"fnstsw    %%ax\n\t"
			"testl     $0x400, %%eax\n\t"
			"jnz       2b\n\t"
			"fstp      %%st(1)\n\t"
			"fsincos\n\t"
			"1:"
			: "=t" (__cosr), "=u" (__sinr) : "0" (__x));
			*__sinx = __sinr;
			*__cosx = __cosr;
	}

	inline double in_atan2(double __y, double __x)
	{
		register long double __value;
		__asm __volatile__("fpatan" : "=t" (__value) : "0" (__x), "u" (__y) : "st(1)");
		return __value;
	}

	#define atan2(X,Y) in_atan2( (X), (Y))

#else
	#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif


double Clock(void)
{
    #ifdef LINUX
	struct tms buff;
	return 0.01*times(&buff);
    #else
	return clock()/((double)CLOCKS_PER_SEC);
    #endif
}

STATUS solveSystemOfEquestions(const double *A11, const double *A12, const double *A21, const double *A22, double *X1, double *X2, const double *B1, const double *B2)
{

	/* This function solves system of linear equestions:

		| A11    A12 || X1 |    | B1 |
		|            ||    | =  |    |
		| A21    A22 || X2 |    | B2 | 
	*/

	const double detA = (*A11)*(*A22) - (*A12)*(*A21);

	if(detA<=EPS_DET)
		return ERROR;

	const double Y1 = (*A22)*(*B1) - (*A12)*(*B2);
	const double Y2 = (*A11)*(*B2) - (*A21)*(*B1);

	*X1 = Y1/detA;
	*X2 = Y2/detA;

	return SUCCESS;
}

void toolbar(unsigned long int step)
{
//    printf(" step: %lu \n",step);
	char bar[NUMBER_OF_STEPS_IN_TOOLBAR+2];
    int i;
	
    for(i=0;i<NUMBER_OF_STEPS_IN_TOOLBAR+1;i++)
	bar[i]=' ';
		    
    bar[NUMBER_OF_STEPS_IN_TOOLBAR+1]='\0';
			
    for(i=0;i<=step;i++)
    {
		bar[i]='=';
		fprintf(stdout,"\r|%s| %d%%",bar,NUMBER_OF_STEPS_IN_TOOLBAR*i);
		fflush(stdout);
    }
							    
    if(step==NUMBER_OF_STEPS_IN_TOOLBAR)
		fprintf(stdout, "\n");
									
}

double findMax(double firstNumber, double secondNumber)
{
	return (firstNumber>=secondNumber) ? firstNumber : secondNumber;
}
