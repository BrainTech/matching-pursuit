/* 1999 06 11/20 */

#ifndef NEW_IOBOOK_H
#define NEW_IOBOOK_H

#define MAX_FIELDS 256

#define ZERO          ((Byte)0)
#define TEXT_INFO     ((Byte)1)
#define DATE_INFO     ((Byte)2) 
#define SIGNAL_INFO   ((Byte)3)
#define DECOMP_INFO   ((Byte)4)

typedef unsigned char  Byte;
typedef unsigned short Word;

typedef struct {
  float energy_percent;
  int   max_number_of_iterations;
  int   dictionary_size;
  char  dictionary_type;
  char  dummy[3];
} decomposition_info;

#define DECOMP_INFO_SIGNATURE "fddcccc"

typedef struct {
  float sampling_freq;
  float points_per_microvolt;
  int   number_of_chanels_in_file;
} signal_info;

#define SIGNAL_INFO_SIGNATURE "ffd"

typedef struct {
  Byte code;
  Byte size;
  void *field;
} TFIELD;

typedef struct {
  int    numOfFields;
  Word   HeaderSize;
  TFIELD field[MAX_FIELDS];
} FILE_HEADER;

typedef struct {
  int   channel;
  int   file_offset;
  int   book_size;
  int   signal_size;
  float signal_energy;
  float book_energy;
} SEG_HEADER;

typedef struct {		
  float scale;
  float frequency;
  float position;
  float modulus;
  float amplitude;
  float phase;
} NEW_ATOM;

#ifdef __cplusplus
extern "C" {
#endif

int   WriteFileHeader(FILE_HEADER *,FILE *);	
int   WriteSegmentHeader(SEG_HEADER *,FILE *);
int   ReadFileHeader(FILE_HEADER *,FILE *);
int   ReadSegmentHeader(SEG_HEADER *,FILE *);	
int   WriteNewAtom(NEW_ATOM *,FILE *);	
int   ReadNewAtom(NEW_ATOM *,FILE *);  
int   checkBookVersion(FILE *);
int   setBookPosition(int,FILE *);
int   countBook(char *);
int   skipHeader(FILE *);
void  initField(FILE_HEADER *);
void *getField(FILE_HEADER *,Byte);
int   addField(FILE_HEADER *,Byte,void *);
void  freeAllFields(FILE_HEADER *);
void  deleteField(FILE_HEADER *,Byte);
int   addDate(FILE_HEADER *);

#ifdef __cplusplus
}
#endif

#define getDecompPtr(head) ((decomposition_info *)getField(head,DECOMP_INFO))
#define getSignalPtr(head) ((signal_info *)getField(head,SIGNAL_INFO))
#define getTextPtr(head)   ((char *)getField(head,TEXT_INFO))
#define getDatePtr(head)   ((char *)getField(head,DATE_INFO))

#define addDecompInfo(head,data) addField(head,DECOMP_INFO,(void *)(data))
#define addSignalInfo(head,data) addField(head,SIGNAL_INFO,(void *)(data))
#define addTextInfo(head,data)   addField(head,TEXT_INFO,(void *)(data))

#endif




