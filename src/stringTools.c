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

#include<ctype.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"include/def.h"
#include"include/stringTools.h"

unsigned short int countWhiteBlanks(const char *line)
{
  unsigned short int numberOfWhiteBlanks = 0;

  while(*line!='\0')
    {
      if(*line==' ' || *line=='\n' || *line=='\t')
	numberOfWhiteBlanks++;
      line++;
    }
	
  return numberOfWhiteBlanks;

}

unsigned short int  countWords(const char *line)
{
  BOOLEAN out = TRUE;
  unsigned short int numberOfWords = 0;
	
  while(*line!='\0')
    {
      if(*line==' ' || *line=='\n' || *line=='\t')
	out = TRUE;
      else if(out)
	{
	  out = FALSE;
	  numberOfWords++;
	}
      line++;
    }
	
  return numberOfWords;
}

STATUS isDecimal(const char *string, const char *mode)
{
  const char *c = string;

  while(*c!='\0')
    {
      if(!isdigit(*c))
	return ERROR;
		
      c++;
    }

  if(strcmp(mode,"uintgz")==0 && atoi(string)>0) /* uintgtz - unsigned integer greater then zero */
    return SUCCESS;
  else if(strcmp(mode,"uint")==0 && atoi(string)>=0) /* uint - unsigned integer */
    return SUCCESS;
  else if(strcmp(mode,"int")==0) /* uint - integer */
    return SUCCESS;
  else
    fprintf(stderr," INCORRECT USE OF isDecimal PROCEDURE \n");


  return ERROR;
}

STATUS isReal(const char *string, const char *mode)
{
  const char *c = string;
  BOOLEAN foundDot = FALSE;

  while(*c!='\0')
    {
      if(!isdigit(*c) && *c!='.')
	return ERROR;
      if(*c=='.')
	foundDot = TRUE;
      c++;
    }

  if(!foundDot)
    return ERROR;

  if(strcmp(mode,"realgz")==0 && (double)atof(string)>0) /* realgz - real, greater then zero */
    return SUCCESS;
  else if(strcmp(mode,"realgez")==0 && atof(string)>=0) /* realgez - real greater or equal zero */
    return SUCCESS;
  else if(strcmp(mode,"real")==0) /* real - ... :-) */
    return SUCCESS;
  else
    fprintf(stderr," INCORRECT USE OF isReal PROCEDURE \n");
	
  return ERROR;
}

void reverseString(char *string)
{
  unsigned char c;
  int i,j;
	
  for(i=0,j=strlen(string)-1;i<j;i++,j--)
    {
      c = string[i];
      string[i] = string[j];
      string[j] = c;
    }
}

char *decimalToString(char *string,int number)
{
  int i, sign;
	
  if((sign=number)<0)
    number = -number;
  i = 0;
	
  do
    {
      string[i++] = (char)(number % 10 + '0');
    }while((number/=10)>0);
	
  if(sign<0)
    string[i++] = '-';
	
  string[i] = '\0';
	
  reverseString(string);

  return string;
}
