/* 1999 06 11/17/20 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "include/new_io.h"

#ifdef INTELSWP

#if  !defined(__BORLANDC__) && !defined(WIN32)
#ifndef __USE_XOPEN
#define __USE_XOPEN
#include <unistd.h>
#undef __USE_XOPEN 
#else
#include <unistd.h>
#endif
#endif

static void INTELfloatToIBM(float f,char *ret) {
  short isav;
  union {
    char  fc[4];
    short fsi[2];
    float ff;
  } gf[2];

  gf[0].ff=f;
  swab((void *)gf[0].fc,(void *)gf[1].fc,4);
  isav=gf[1].fsi[0];
  gf[1].fsi[0]=gf[1].fsi[1]; gf[1].fsi[1]=isav;
  memcpy((void *)ret,(void *)&gf[1].ff,4);
}

static float IBMfloatToINTEL(const char *f) {
  short isav;
  union	{
    char  fc[4];
    short fsi[2];
    float ff;
  } gf[2];

  memcpy((void *)gf[0].fc,(void *)f,4);
  swab((void *)gf[0].fc,(void *)gf[1].fc,4);
  isav=gf[1].fsi[0];
  gf[1].fsi[0]=gf[1].fsi[1]; gf[1].fsi[1]=isav;
  return gf[1].ff;
}

static int IBMintToINTEL(const char *f) {
  short isav;
  union	{
    char fc[4];
    short fsi[2];
    int ff;
  } gf[2];
  
  memcpy((void *)gf[0].fc,(void *)f,4);
  swab((void *)gf[0].fc,(void *)gf[1].fc,4);
  isav=gf[1].fsi[0];
  gf[1].fsi[0]=gf[1].fsi[1]; gf[1].fsi[1]=isav;
  return gf[1].ff;
}

static void INTELintToIBM(int f,char *ret) {
  short isav;
  union {
    char  fc[4];
    short fsi[2];
    int   ff;
  } gf[2];

  gf[0].ff=f;
  swab((void *)gf[0].fc,(void *)gf[1].fc,4);
  isav=gf[1].fsi[0];
  gf[1].fsi[0]=gf[1].fsi[1]; gf[1].fsi[1]=isav;
  memcpy((void *)ret,(void *)&gf[1].ff,4);
}

static char *IBMStructToINTEL(const char *format,const char *input,
			      char *output) {
  register int i,index=0;
  float ftmp;
  int   itmp;

  for(i=0 ; format[i]!='\0' ; i++)
    switch(format[i]) {
    case 'f': 
      ftmp=IBMfloatToINTEL(input+index);
      memcpy((void *)(output+index),(void *)&ftmp,sizeof(float));
      index+=sizeof(float);
      break;
    case 'd': 
      itmp=IBMintToINTEL(input+index);
      memcpy((void *)(output+index),(void *)&itmp,sizeof(int));
      index+=sizeof(int);
      break;
    case 'c':
      output[index]=input[index];
      index++;
      break;
    case 'h':
      swab((void *)(input+index),(void *)(output+index),sizeof(short));
      index+=sizeof(short);
      break;
    }
  return output;
}

static char *INTELStructToIBM(const char *format,const char *input,
			      char *output) {
  register int i,index=0;

  for(i=0 ; format[i]!='\0'; i++)
    switch(format[i]) {
    case 'f':
      INTELfloatToIBM(*((float *)(input+index)),output+index);
      index+=sizeof(float);
      break;
    case 'd':
      INTELintToIBM(*((int *)(input+index)),output+index);
      index+=sizeof(int);
      break;
    case 'c':
      output[index]=input[index];
      index++;
      break;
    case 'h':
      swab((void *)(input+index),(void *)(output+index),sizeof(short));
      index+=sizeof(short);
      break;
    }
  return output;
}

#endif                  

int WriteSegmentHeader(SEG_HEADER *head,FILE *plik) {
#ifdef INTELSWP
  INTELintToIBM(head->channel,(char *)&head->channel);
  INTELintToIBM(head->file_offset,(char *)&head->file_offset);
  INTELintToIBM(head->book_size,(char *)&head->book_size);
  INTELintToIBM(head->signal_size,(char *)&head->signal_size);
  INTELfloatToIBM(head->signal_energy,(char *)&head->signal_energy);
  INTELfloatToIBM(head->book_energy,(char *)&head->book_energy);
#endif
  if(fwrite((void *)head,sizeof(SEG_HEADER),1,plik)==0U) 
    return -1;
  return 0;
}

int WriteNewAtom(NEW_ATOM *atom,FILE *plik) {
#ifdef INTELSWP
  INTELfloatToIBM(atom->scale,(char *)&atom->scale);
  INTELfloatToIBM(atom->frequency,(char *)&atom->frequency);
  INTELfloatToIBM(atom->position,(char *)&atom->position);
  INTELfloatToIBM(atom->modulus,(char *)&atom->modulus);
  INTELfloatToIBM(atom->amplitude,(char *)&atom->amplitude);
  INTELfloatToIBM(atom->phase,(char *)&atom->phase);
#endif
  if(fwrite((void *)atom,sizeof(NEW_ATOM),1,plik)==0U) 
    return -1;
  return 0;
}

int ReadSegmentHeader(SEG_HEADER *head,FILE *plik) {
  if(fread((void *)head,sizeof(SEG_HEADER),1,plik)==0U)
    return -1;
#ifdef INTELSWP
  head->channel=      IBMintToINTEL((char *)&head->channel);
  head->file_offset=  IBMintToINTEL((char *)&head->file_offset);
  head->signal_size=  IBMintToINTEL((char *)&head->signal_size);
  head->book_size=    IBMintToINTEL((char *)&head->book_size);
  head->signal_energy=IBMfloatToINTEL((char *)&head->signal_energy);
  head->book_energy  =IBMfloatToINTEL((char *)&head->book_energy);
#endif
  return 0;
}

int ReadNewAtom(NEW_ATOM *atom,FILE *plik) {
  if(fread((void *)atom,sizeof(NEW_ATOM),1,plik)==0U)
    return -1;
#ifdef INTELSWP
  atom->scale=    IBMfloatToINTEL(  (char *)&atom->scale);
  atom->frequency=IBMfloatToINTEL((char *)&atom->frequency);
  atom->position= IBMfloatToINTEL((char *)&atom->position);
  atom->modulus=  IBMfloatToINTEL((char *)&atom->modulus);
  atom->amplitude=IBMfloatToINTEL((char *)&atom->amplitude);
  atom->phase  =  IBMfloatToINTEL((char *)&atom->phase);
#endif
  return 0;
}

int checkBookVersion(FILE *file) {
  const long fpos=ftell(file);
  int ok=0;
  char magic[5];

  fseek(file,0L,SEEK_SET);
  if(fread((void *)magic,4,1,file)!=1) 
    return -1;
  magic[4]='\0';
  if(strcmp(magic,"MPv4")!=0)
    ok=-1;
  fseek(file,fpos,SEEK_SET);
  return ok;
}

int setBookPosition(int where,FILE *file) {
  SEG_HEADER  seg_head;	
  long pos;
  int i;
 
  if(skipHeader(file)==-1)
    return -1;

  pos=ftell(file);
  for(i=0 ; i<where ; i++) {
    if(ReadSegmentHeader(&seg_head,file)==-1)
      return -1;
    pos+=(long)seg_head.book_size*sizeof(NEW_ATOM)+sizeof(SEG_HEADER);
    if(fseek(file,pos,SEEK_SET)!=0)
      return -1;
  }
  return 0;
}

unsigned long NumberOfAllWaveForms;

static char *trim(char *string) {
  register int i;
  
  for(i=strlen(string)-1 ; i>=0 ; i--)
    if(!isspace(string[i])) {
      string[i+1]='\0';
      break;
    }
  return string;
}

int countBook(char *filename) {                           
  SEG_HEADER  seg_head;
  FILE *file;
  int count;
  long pos=0L;
  
  NumberOfAllWaveForms=0UL;
  if((file=fopen(trim(filename),"rb"))==0)
    return -1;

  if(skipHeader(file)==-1) {
    fclose(file);
    return -1;
  }

  pos=ftell(file);
  for(count=0 ; ; count++) {
    if(ReadSegmentHeader(&seg_head,file)==-1)
      break;
    pos+=(long)seg_head.book_size*sizeof(NEW_ATOM)+sizeof(SEG_HEADER);
    NumberOfAllWaveForms+=(unsigned long)seg_head.book_size;
    if(fseek(file,pos,SEEK_SET)!=0)
      break;
  }
  
  fprintf(stdout,"LICZBA ZESTAWOW ATOMOW W PLIKU    : %d\n"
	  "LICZBA WSZYSTKICH ATOMOW W PLIKU  : %lu\n",
	  count,NumberOfAllWaveForms);

  fclose(file);
  return count;
}

void initField(FILE_HEADER *head) {
  register int i;

  head->numOfFields=0;
  for(i=0 ; i<MAX_FIELDS ; i++) {
    head->field[i].code=ZERO;
    head->field[i].size=ZERO;
    head->field[i].field=NULL;
  }
}

void freeAllFields(FILE_HEADER *head) {
  int i;
  for(i=0 ; i<head->numOfFields ; i++)
    if(head->field[i].field!=NULL)
      free((void *)head->field[i].field);
  head->numOfFields=0;
}

void deleteField(FILE_HEADER *head,Byte code) {
  register int i,j;

  for(i=0 ; i<head->numOfFields ; i++)
    if(head->field[i].code==code) {
      free(head->field[i].field);
      head->numOfFields--;
      for(j=i ; j<head->numOfFields ; j++)
	head->field[j]=head->field[j+1];
      break;
    }
}

static int getSizeOf(const char *sign) {
  register int i,size=0;
  
  for(i=0 ; sign[i]!='\0' ; i++)
    switch(sign[i]) {
    case 'f': size+=sizeof(float); break;
    case 'd': size+=sizeof(int);   break;
    case 'c': size++;              break;
    case 'h': size+=sizeof(short); break;
    }
  return size;
}

int addField(FILE_HEADER *head,Byte code,void *ptr) {
  int index,size=0;
#ifdef INTELSWP
  char *sign=NULL;
#endif

  if(head->numOfFields>=MAX_FIELDS)
    return -1;

  deleteField(head,code);
  head->field[index=head->numOfFields].code=code;
  switch(code) {
  case DATE_INFO:
  case TEXT_INFO: 
    size=strlen((const char *)ptr)+1;
    if(size>=256)
      return -1;
    head->field[index].size=(Byte)size;
    if((head->field[index].field=malloc(size*sizeof(char)))==NULL)
      return -1;
    strcpy((char *)head->field[index].field,(const char *)ptr);
    break;
  case SIGNAL_INFO:
  case DECOMP_INFO:
    if(code==SIGNAL_INFO) 
      size=getSizeOf(SIGNAL_INFO_SIGNATURE);
    else if(code==DECOMP_INFO)
      size=getSizeOf(DECOMP_INFO_SIGNATURE);

    head->field[index].size=(Byte)size;
    if((head->field[index].field=malloc(size*sizeof(char)))==NULL)
      return -1;
#ifdef INTELSWP    
    if(code==SIGNAL_INFO) 
      sign=SIGNAL_INFO_SIGNATURE;
    else 
      sign=DECOMP_INFO_SIGNATURE;

    INTELStructToIBM(sign,(char *)ptr,head->field[index].field);
#else
    memcpy(head->field[index].field,ptr,size);
#endif
    break;
  default: 
    return -1;
  }

  head->numOfFields++;
  return 0;
}

static int readField(FILE_HEADER *head,Byte *block,FILE *file) {
  int index;
  Byte code;

  if(head->numOfFields>=MAX_FIELDS)
    return -1;
  
  if(fread((void *)&code,sizeof(Byte),1,file)!=1)
    return -1;
  if(fread((void *)block,sizeof(Byte),1,file)!=1)
    return -1;

  head->field[index=head->numOfFields].code=code;
  head->field[index].size=*block;

  if((head->field[index].field=malloc(*block))==NULL)
    return -1;

  if(fread(head->field[index].field,*block,1,file)!=1) {
    free((void *)head->field[index].field);
    return -1;
  }

  head->numOfFields++;
  return 0;
}

void *getField(FILE_HEADER *head,Byte code) {
  static char buff[256];
  int i;
#ifdef INTELSWP
  char *sign=NULL;
#endif

  for(i=0 ; i<head->numOfFields ; i++)
    if(head->field[i].code==code) {
      switch(code) {
      case DATE_INFO:
      case TEXT_INFO:
	memcpy(buff,head->field[i].field,head->field[i].size);
	return (void *)buff;
      case SIGNAL_INFO:
      case DECOMP_INFO:
#ifdef INTELSWP
	if(code==SIGNAL_INFO) 
	  sign=SIGNAL_INFO_SIGNATURE;
	else 
	  sign=DECOMP_INFO_SIGNATURE;
	
	IBMStructToINTEL(sign,(char *)head->field[i].field,buff);
	return (void *)buff;
#else 
	memcpy(buff,head->field[i].field,head->field[i].size);
	return (void *)buff;
#endif
      }
    }
	    
  return NULL;
}

static Word sizeOfHead(FILE_HEADER *head) {
  Word size=4+sizeof(Word);
  register int i;

  for(i=0 ; i<head->numOfFields ; i++)
    size+=(Word)(2+head->field[i].size);
  return size;
}

int ReadFileHeader(FILE_HEADER *head,FILE *file) {
  char magic[5];
  Byte block;
  Word size=4;

  fseek(file,0L,SEEK_SET);
  if(fread((void *)magic,4,1,file)!=1)
    return -1;

  magic[4]='\0';
  if(strcmp(magic,"MPv4")!=0)
    return -1;

  freeAllFields(head);
  if(fread((void *)&head->HeaderSize,sizeof(Word),1,file)!=1)
    return -1;
  size+=sizeof(Word);
#ifdef INTELSWP
  swab((void *)&head->HeaderSize,(void *)&head->HeaderSize,2);
#endif

  while(size<head->HeaderSize) {
    if(readField(head,&block,file)==-1)
      return -1;
    size+=block+2*sizeof(Byte);
  }
  fseek(file,head->HeaderSize,SEEK_SET);
  return 0;
}

int WriteFileHeader(FILE_HEADER *head,FILE *file) {
  char magic[]={ 'M', 'P', 'v', '4' };
  Word size;
  int i;

  head->HeaderSize=(size=sizeOfHead(head));
  if(fwrite((void *)magic,4,1,file)!=1)
    return -1;

#ifdef INTELSWP
  swab((void *)&size,(void *)&size,2);
#endif
  if(fwrite((void *)&size,2,1,file)!=1)
    return -1;
  for(i=0 ; i<head->numOfFields ; i++) {
    if(fwrite((void *)&head->field[i].code,1,1,file)!=1)
      return -1;
    if(fwrite((void *)&head->field[i].size,1,1,file)!=1)
      return -1;
    if(fwrite(head->field[i].field,head->field[i].size,1,file)!=1)
      return -1;
  }
  return 0;
}

int skipHeader(FILE *file) {
  Word size;
  if(checkBookVersion(file)==-1) 
    return -1;
  fseek(file,4L,SEEK_SET);
  if(fread((void *)&size,2,1,file)!=1)
    return -1;
#ifdef INTELSWP
  swab((void *)&size,(void *)&size,2);
#endif
  fseek(file,size,SEEK_SET);
  return 0;
}

int addDate(FILE_HEADER *head) {
  time_t tim;
  char buff[128];

  time(&tim);
  strcpy(buff,ctime(&tim));
  buff[strlen(buff)-1]='\0';
  return addField(head,DATE_INFO,(void *)buff);
}


