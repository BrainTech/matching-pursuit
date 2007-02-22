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

#include<stdio.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include"include/def.h"
#include"include/tools.h"

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

double ddot(const int *n, const double *dx, const int *incx, const double *dy, const int *incy)
{
	/* System generated locals */
	int i__1;
	double ret_val;

	/* Local variables */
	static int i, m;
	static double dtemp;
	static int ix, iy, mp1;

/*      forms the dot product of two vectors.   
        uses unrolled loops for increments equal to one.   
        jack dongarra, linpack, 3/11/78.   
        modified 12/3/93, array(1) declarations changed to array(*)   
   
        Parameter adjustments   
        Function Body */

	#define DY(I) dy[(I)-1]
	#define DX(I) dx[(I)-1]

	ret_val = 0.;
	dtemp = 0.;

	if (*n <= 0)
		return ret_val;

	if (*incx == 1 && *incy == 1)
		goto L20;

	/* code for unequal increments or equal increments not equal to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0)
		ix = (-(*n) + 1) * *incx + 1;

	if (*incy < 0)
		iy = (-(*n) + 1) * *incy + 1;

	i__1 = *n;

	for (i = 1; i <= *n; ++i)
	{
		dtemp += DX(ix) * DY(iy);
		ix += *incx;
		iy += *incy;
		/* L10: */
	}

	ret_val = dtemp;
	return ret_val;

	/* code for both increments equal to 1   
           clean-up loop */

	L20:
		m = *n % 5;
		if (m == 0)
			goto L40;

		i__1 = m;
		for (i = 1; i <= m; ++i)
		{
			dtemp += DX(i) * DY(i);
			/* L30: */
		}
		if (*n < 5)
			goto L60;

	L40:
		mp1 = m + 1;
		i__1 = *n;
		for (i = mp1; i <= *n; i += 5)
		{
			dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) * 
			DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
			/* L50: */
		}

	L60:
		ret_val = dtemp;
		return ret_val;

} /* ddot_ */

