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

#include<math.h>
#include<stdio.h>
#include"atom.h"
#include"cmd.h"
#include"def.h"
#include"dic.h"
#include"io_mp5.h"
#include"mmp.h"
#include"mp5.h"
#include"smp.h"
#include"stringTools.h"
#include"types.h"

unsigned char applicationMode = 0x00;

#ifdef __MINGW32__
	#define bzero(ptr,size) memset (ptr, 0, size);
	#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

static void countMeanSignalOverChannels(MP5Parameters *mp5Parameters)
{
	const unsigned int offsetDimension   = mp5Parameters->offsetDimension;
	const unsigned int marginalDimension = mp5Parameters->marginalDimension;

	unsigned int sample;
	double tmpDataValue;
	
	double **multiChannelSignalTable = mp5Parameters->multiChannelSignalTable;
	double *meanSignalTable          = mp5Parameters->meanSignalTable;
	int channel;

	for(sample=0;sample<offsetDimension;sample++)
	{
		tmpDataValue = 0.0;
		
		for(channel=0;channel<mp5Parameters->numberOfAnalysedChannels;channel++)
		{
			tmpDataValue+= (*(*(multiChannelSignalTable + mp5Parameters->chosenChannels[channel] - 1) + marginalDimension + sample));
		}
		*(meanSignalTable + marginalDimension + sample) = tmpDataValue/mp5Parameters->numberOfAnalysedChannels;
	}
}

int main(int argc, char *argv[])
{
	char infoMessage[LENGTH_OF_INFO_MESSAGE];
	bzero((void *)infoMessage,LENGTH_OF_INFO_MESSAGE);

    ConfigFile configFile = {"",
							 NULL,
							 NULL,
							 NULL};
   
    MP5Parameters  mp5Parameters = {"","","","",
									NULL,NULL,
									0,0,
									0,0,0,
									0.0,
									0,NULL,
									0,NULL,
									0,
									0x0,
									0.0,
									NULL,NULL,
									0x0,0x0,
									0,0,0,
									0,0,0,0,
									NULL,NULL,NULL,
									0x00,
									NULL,
									NULL,
									NULL,
									NULL,
									0x00,
									NULL,
									NULL,
									NULL,
									NULL,
									NULL,
									NULL,
									NULL,NULL,NULL,NULL,NULL,NULL,
									0.0,0.0,0.0,0.0,
									NULL,NULL,
									0,
									0.0,
									NULL,NULL,
									0x0,0x0,
									NULL,
									0x00,
									0,
									0x00,
									FULL
								};

    Dictionary dictionary = {0,0,0,
				   		     0.0,0.0,0,0,
							 0,
							 0.0,0.0,
							 0.0,
							 0,
							 0.0,
							 NULL,NULL,NULL,
							 NULL,NULL,NULL,
							 0,0,0,0,
							 YES,YES,YES,
							 NULL};
		     	
    if(argc>6 || argc==1
			  || (argc==2 && ((strcmp(argv[1],"-g")!=0) && (strcmp(argv[1],"-e")!=0)))
			  || ((strcmp(argv[1],"-g")==0) && (argc>2))
			  || ((strcmp(argv[1],"-e")==0) && (argc>2))
			  || (argc==5 && (((strcmp(argv[1],"-m")!=0) && (strcmp(argv[2],"-a")!=0) && (strcmp(argv[4],"-N")!=0)) ||
				              ((strcmp(argv[1],"-m")!=0) && (strcmp(argv[2],"-N")!=0) && (strcmp(argv[4],"-a")!=0)))))
    {
		fprintf(stderr," \n");
		fprintf(stderr," ERROR: \n");
		fprintf(stderr," INCORRECT CALL OF mp5Parameters PROGRAM \n");
		fprintf(stderr," THE PROPER USE IS AS FOLLOWS: \n");
		fprintf(stderr," mp5 -g                         - generate default config file \n");
		fprintf(stderr," mp5 -e                         - generate xml file with errors codes and translations\n");
		fprintf(stderr," mp5 -t [ name of config file ] - test config file \n");
		fprintf(stderr," mp5 -f [ name of config file ] - process of mp decomposition (user mode)\n");
		fprintf(stderr," mp5 -x [ name of config file ] - process of mp decomposition (server mode)\n");
		fprintf(stderr," \n");
		return ERROR;
    }

    if(strcmp(argv[1],"-g")==0)
		applicationMode = CONFIG_MODE;
	else if(strcmp(argv[1],"-e")==0)
		applicationMode = ERROR_MODE;						
	else if(strcmp(argv[1],"-t")==0)
		applicationMode = TEST_PARAMETERS_MODE ;
    else if(strcmp(argv[1],"-f")==0)
		applicationMode = PROCESS_USER_MODE;
    else if(strcmp(argv[1],"-x")==0)
		applicationMode = PROCESS_SERVER_MODE;		
    else
    {
		fprintf(stderr," \n");
		fprintf(stderr," ERROR: \n");
		fprintf(stderr," INCORRECT CALL OF mp5Parameters PROGRAM \n");
		fprintf(stderr," UNKNOW PARAMETER: %s \n",argv[1]);
		fprintf(stderr," THE PROPER USE IS AS FOLLOWS:\n");
		fprintf(stderr," mp5 -g                         - generate default config file\n");
		fprintf(stderr," mp5 -e                         - generate xml file with errors codes and translations\n");
		fprintf(stderr," mp5 -t [ name of config file ] - test config file \n");
		fprintf(stderr," mp5 -f [ name of config file ] - process of mp decomposition (USER MODE)\n");
		fprintf(stderr," mp5 -x [ name of config file ] - process of mp decomposition (SERVER MODE)\n");
		fprintf(stderr," \n");
		return ERROR;
    }

	if((applicationMode & PROCESS_USER_MODE) || (applicationMode & CONFIG_MODE) || (applicationMode & ERROR_MODE))
	{
		printf("\n");
		printf("  Matching Pursuit V (version 2008-01-08)\n");
		printf("  Department of Biomedical Physics at Warsaw University\n");
		printf("  http://brain.fuw.edu.pl, http://eeg.pl\n");
		printf("  Compiled: %s, %s \n", __DATE__ , __TIME__ );

		#ifdef __VERSION__
			printf("  Compiler: %s \n", __VERSION__);
		#endif

		printf("\n");
	}
    
	if(applicationMode & CONFIG_MODE)
    {
		FILE *exampleFile = fopen("example.set","wt");
		char exampleTest[] =
		{"#\n\
# PLIK KONFIGURACYJNY PROGRAMU mp5 -- SPECYFIKACJA\n\
#\n\
# CONFIGURATION FILE FOR THE mp5 PROGRAMME -- SPECIFICATION\n\
#\n\
#\n\
# Plik konfiguracyjny sklada sie z argumentow i ich wartosci oraz komentarzy.\n\
# Komentarzem jest linia rozpoczynajaca sie od znaku '#'. Komentarzy moze byc dowolna ilosc.\n\
# Lista argumentow jest scisle okreslona. Brak lub dopisanie jakiegos argumentu\n\
# powoduje przerwanie dzialania programu mp5Parameters.\n\
# Obecnie program mp5Parameters rozpoznaje 14 nastepujacych argumentow:\n\
#\n\
# The configuration file consists of arguments, their values, and comments.\n\
# Each line starting with the # character is treated as a comment. There is no space restriction on comments.\n\
# The list of arguments is strictly defined. The lack of an argument or an extra argument\n\
# results in termination of the mp5Parameters programme.\n\
# Currently, the mp5Parameters programme recognizes the 19 following arguments:\n\
#\n\
# nameOfDataFile\n\
# extensionOfResultsFile\n\
# writingMode\n\
# sizeOfHeader\n\
# sizeOfTail\n\
# samplingFrequency\n\
# formatOfData\n\
# numberOfChannels\n\
# chosenChannels\n\
# numberOfPointsInOffset\n\
# chosenOffsets\n\
# typeOfDictionary\n\
# dilationFactor\n\
# reinitDictionary\n\
# maximalNumberOfIterations\n\
# energyPercent\n\
# MP\n\
# analiticalDotProduct\n\
# pointsPerMicrovolt\n\
#\n\
#\n\
# Argumenty oraz ich wartosci moga byc podane w dowolnej kolejnosci.\n\
#\n\
# Arguments, followed by their values, can be given in any order.\n\
#\n\
#\n\
# W skladni pliku konfiguracyjnego wyrozniono nastepujace cztery znaki specjalne:\n\
# The syntax of the configuration file treats the following four characters as special:\n\
#\n\
# '#' -- Komentarz            Comment\n\
# '\' -- Lamanie linii        Linebreak\n\
# '-' -- Zakres wartosci      Range of values\n\
# ' ' -- Znak przestankowy    Pause\n\
#\n\
# Ad. '#':\n\
# Wszystkie linie zaczynajace sie od '#' sa traktowane jak komentarze\n\
# i nie sa poddawane dalszej analizie.\n\
# Komentarzy moze byc w pliku konfiguracyjmym dowolna ilosc\n\
# Mozna je takze na wlasny uzytek dopisywac. Aby jednak uchronic\n\
# uzytkownikow przed pozniejszym krazeniem po sieci plikow konfiguracujnych\n\
# z roznymi wersjami komentarzy, zawsze mozna wygenerowac domyslny\n\
# plik konfiguracyjny wydajac polecenie mp5 -g.\n\
#\n\
# All lines starting with # are treated as comments\n\
# and are thus not analysed any further.\n\
# There is no space restriction on comments.\n\
# One can add new comments for his/her own needs. However, to protect\n\
# future users from multiple versions, it is possible to generate the default\n\
# configuration file by typping mp5 -g command.\n\
#\n\
# Ad. '\':\n\
# Polecenia wydawane dla programu mp5Parameters moga byc dosyc dlugie,\n\
# o czym przekonamy sie dalej. Aby nie pisac nieskonczenie dlugich linii,\n\
# dopusczono ich lamanie -- tak samo jak w jezyku C.\n\
# Commands for mp5Parameters programme can be relatively long.\n\
#\n\
# In order not to write infinitely long lines,\n\
# one can break them, as it is in the C language.\n\
#\n\
# Ad. '-':\n\
# Analizie programem mp5Parameters moga byc poddane dane zawierajace wiele kanalow i skladek.\n\
# Kanaly lub skladki, ktore maja byc analizowane przez program, nalezy wymienic w pliku konfiguracyjnym.\n\
# W przypadku gdy lista kanalow lub skladek bylaby dluga, aby uniknac zmudnego ich wypisywania,\n\
# mozna podac zakres kanalow 'od-do', np. 1-23.\n\
# UWAGA I:  Miedzy liczbami nie moze wystapic spacja.\n\
# UWAGA II:	Lamanie linii ('\') oraz uzywanie znaku zakres ('-') jest dopuszczalne tylko dla argumentow\n\
#           chosenChannels i chosenOffsets. Wartosci wszystkich pozostalych parametrow powinny byc jednoargumentowe.\n\
#\n\
# mp5Parameters programme can analyse data containing many channels and trials/epochs.\n\
# Channels and epochs to be analysed should be listed in the configuration file.\n\
# When the list of channels or epochs is about to be long, one can specify the range 'from-to', e.g. 1-23.\n\
# NOTE I:  There can be no space in between the numbers.\n\
# NOTE II: The special characters '\' (linebreak)  and '-' range of values can be applied\n\
#          only to chosenChannels and chosenOffset arguments. The values of all other parameters\n\
#          should be single arguments.\n\
#\n\
#\n\
# Nazwa pliku z danymi.\n\
#\n\
# Name of data set.\n\
#\n\
\n\
nameOfDataFile transien.dat\n\
\n\
#\n\
# Nazwa rozszerzenia pliku z wynikami.\n\
# Nazwy plikow wyjsciowych sa generowane automatycznie.\n\
# Do nazwy pliku z danymi dodawane jest rozszerzenie postaci:\n\
# _smp.b dla algorytmu jednokanalowego (SMP) lub _mmp.b dla algorytmow wielokanalowych (MMP1, MMP2, MMP3).\n\
# W takim przypadku pelna nazwa pliku wyjsciowego jest postaci:\n\
# nameOfDataFile_extensionOfResultsFile_smp.b lub nameOfDataFile_extensionOfResultsFile_mmp.b.\n\
# W przypadku, gdy taki sposob tworzenia nazwy pliku wyjsciowego jest niewystarczajacy,\n\
# mozna dodac do nazwy wlasne rozszerzenie.\n\
# Jezeli wartosc argumentu jest rowna NONE, nazwa pliku wyjsciowego jest generowana\n\
# bez zadnego dodatkowego rozszerzenia.\n\
#\n\
# Name of the extension of the file containg results.\n\
# Names of the output files are generated automatically.\n\
# The name of the file with the data is appended with\n\
# _smp.b in the case of the single channel mp algorithm (SMP), or _mmp.b in the case of multichannel mp algorithms (MMP1, MMP2, MMP3).\n\
# Then, the full name looks as follows:\n\
# nameOfDataFile_extensionOfResultsFile_smp.b or nameOfDataFile_extensionOfResultsFile_mmp.b\n\
# When the resulting file name is not sufficient, one can add his/her own extension.\n\
# If the argument's value is set to NONE,\n\
# the file name is generated without any additional extension.\n\
#\n\
\n\
extensionOfResultsFile NONE\n\
\n\
#\n\
# Mozna podac nazwe katalogu, do ktorego maja byc zrzucane wyniki.\n\
#\n\
# One can specify a directory to write the results to.\n\
#\n\
\n\
nameOfOutputDirectory ./\n\
\n\
#\n\
# Tryb zapisu wynikow.\n\
# Dowzolone sa nastepujace wartosci:\n\
# CREATE -- Utworz nowy plik do zapisu wynikow. Jesli juz plik o takiej nazwie istnieje, to go zamaz.\n\
# APPEND -- Dolacz wyniki do istniejacego juz pliku, lub utworz nowy plik, jesli plik o danej nazwie nie istnieje.\n\
#\n\
# Results writing mode.\n\
# The following values are allowed:\n\
# CREATE -- Create a new file for results writing. Replace an already existing file with same name.\n\
# APPEND -- Append the results to an already existing file or create a new file if there is no file to append to.\n\
#\n\
\n\
writingMode CREATE\n\
\n\
\n\
#\n\
# Wielkosc naglowka w pliku z danymi.\n\
# W przypadku, gdy plik z danymi jest plikiem tesktowym,\n\
# sizeOfHeader oznacza liczbe linii w naglowku.\n\
# W przypadku gdy dane sa zapisane w pliku binarnym, rozmiar naglowka\n\
# powinien byc wyrazony w bajtach.\n\
# Wartosc argumentu: liczba calkowita dodatnia lub zero.\n\
#\n\
# Size of the header in the data file.\n\
# If the data file is a text file, sizeOfHeader denotes the number of lines in the header.\n\
# If the data file is a binary file, the size should be given in bytes.\n\
# Argument's value: positive integer number or zero.\n\
#\n\
\n\
sizeOfHeader 0\n\
\n\
#\n\
# Czasami do danych EEG dolaczany jest jeszcze \"ogon\".\n\
# Jego rozmiar powinien byc podany w taki sam sposob, jak zostalo okreslone\n\
# dla rozmiaru naglowka.\n\
# Wartosc argumentu: liczba calkowita dodatnia lub zero.\n\
#\n\
# Sometimes, EEG data contain the so-called \"tail\". In such a case,\n\
# its size should be given in the same way as for the header.\n\
# Argument's value: positive integer or zero.\n\
#\n\
\n\
sizeOfTail 0\n\
\n\
#\n\
# Czestosc probkowania.\n\
# Wartosc argumentu: liczba zmiennoprzecinkowa, dodatnia.\n\
# Obowiazkowo w wartosci musi wystapic kropka dziesietna.\n\
#\n\
# Sampling frequency.\n\
# Argument's value: floating-point positive number.\n\
# Decimal dot is obligatory.\n\
#\n\
\n\
samplingFrequency 1017.25\n\
\n\
#\n\
# Format danych.\n\
# W przypadku formatu ASCII zaklada sie, ze w pliku nie ma pustych linii.\n\
# Wartosc argumentu: ASCII/SHORT/FLOAT.\n\
#\n\
# Data format.\n\
# In the case of ASCII format, it is assumed that the file does not contain empty lines.\n\
# Argument's value: ASCII/SHORT/FLOAT.\n\
#\n\
\n\
formatOfData FLOAT\n\
\n\
#\n\
# Liczba kanalow w pliku.\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
#\n\
# Number of channels in a file.\n\
# Argument's value: positive integer.\n\
#\n\
\n\
numberOfChannels 1\n\
\n\
#\n\
# Kanaly, ktore maja zostac poddane analizie.\n\
# Liczby calkowite dodatanie, oddzielone przecinkami lub laczone znakiem '-'.\n\
# Kanaly sa numerowane od 1.\n\
# Przyklad skladni:\n\
# Plik zawiera 128 kanalow. Analizie maja zostac poddane nastepujace kanaly:\n\
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 25, 27, 29, 31, 33, 35, 37, 45, 65, 66, 67, 121, 122, 123, 124, 125, 126, 127, 128.\n\
# Mozna je wymienic w pliku w nastepujacy sposob:\n\
#\n\
# Channels to be analysed.\n\
# Positive integers, separated by commas or joined with '-'.\n\
# Channel numeration starts with 1.\n\
# An example of the syntax:\n\
# A file contains 128 channels. The following channels are to be analysed:\n\
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 25, 27, 29, 31, 33, 35, 37, 45, 65, 66, 67, 121, 122, 123, 124, 125, 126, 127, 128.\n\
# They can be listed in the following way:\n\
#\n\
# chosenChannels 1-9 25 27 29 31 \\n\
#                33 35 37 45 65-67 \\n\
#                121-128\n\
#\n\
\n\
chosenChannels 1\n\
\n\
#\n\
# Dlugosc skladki sygnalu, mierzona jako liczba probek,\n\
# ktora jednorazowo ma zostac poddana analizie.\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
#\n\
# Length of an epoch, measured in the number of samples to be analysed.\n\
# Argument's value: positive integer.\n\
#\n\
\n\
numberOfPointsInOffset 512\n\
\n\
#\n\
# Numery skladek, ktore chcemy analizowac.\n\
# Nomenklatura taka sama jak przy chosenChannels.\n\
#\n\
# Numbers of epochs to be analysed.\n\
# Same nomenclature as for chosenChannels.\n\
#\n\
\n\
chosenOffsets 1\n\
\n\
#\n\
# Rodzaj slownika.\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
# OCTAVE_FIXED: Slownik optymalizowany. Wartosci parametery atomow {u, f, s} sa obliczane w taki sposob,\n\
#               by zbudowany slownik jak najefektywniej pokrywal przestrzen parametrow Atoma.\n\
# OCTAVE_STOCH: Slownik optymalny obliczany taj jak wyzej, z ta roznica, ze parametry {u, f, s} sa losowane\n\
#               wokol obliczonych optymalnych wartosci.\n\
#\n\
# Type of the dictionary.\n\
# Argument's value: positive integer number.\n\
# OCTAVE_FIXED: Optimised dictionary. The values of atom parameters {u, f, s} are computed such that the\n\
#               built dictionary covers the space of atom parameters most efficiently.\n\
# OCTAVE_STOCH: Optimal dictionary computed as above, yet, the parameters {u, f, s} are randomized in the\n\
#               neighbourhood of the computed optimal values.\n\
#\n\
\n\
typeOfDictionary OCTAVE_FIXED\n\
\n\
#\n\
# \"Gestosc slownika\".\n\
# Wartosc argumentu: liczba zmiennoprzecinkowa dodatnia.\n\
#\n\
# \"Dictionary's density\".\n\
# Argument's value: positive floating-point number.\n\
#\n\
\n\
dilationFactor 1.7\n\
\n\
randomSeed auto\n\
\n\
#\n\
# Obecna wersja programu MP -- mp5 -- ma zaprogramowany Slownik Optymalny. Stwarza to mozliwosc dokonania pewnych optymalizacji\n\
# numerycznych, m.in. zapamietywania talbic z wartosciami funkcji sin/cos/exp. Niestety, takie tablice zabieraja bardzo duza\n\
# ilosc pamieci, co czesto uniemozliwia dokonywanie obliczen. Aby zmniejszyc zapotrzebowanie na RAM, postanowiono zapamietywac\n\
# wartosci sin/cos w specjalny sposob. Po zapoznaniu sie ze slownikiem optymalnym mozna zauwazyc, ze dla danej sakli DS\n\
# zostaje obliczony krok w czestosci, DF. Atomy dla zadanej skali powinny miec czestosci rozlozone co DF w zakrecie (0 PI).\n\
# W zwiazku z tym dla skali DS czestosc atoma jest wielokrotnoscia pewnej podstawowej czestosci, dalej nazywanej DF0.\n\
# Pakowanie tablic sin/cos w pamieci komputera polega na tym, ze majac stablicowana funkcje sin o okresie np. 10 probek\n\
# mozna poprzez wybieranie z opdowiednim skokiem tych probek szybko wygenerowac funkcje sin o okresach 9x, 8x, ... , 2x mniejszych.\n\
# Podobna sytucacje mamy w slowniku optymalnym, gdzie czestosci sa wielokrotnoscia pewnej czestosci bazowej,\n\
# a zatem ich okresy sa dzielnikami pewnego okresu podstawowego. W programie MP5 zapamietywane sa zatem dla zadanej skali DS\n\
# tylko sinusy o podstawowym okresie DF0. Niestety, powstaja pewne problemy, jezeli chcemy miec slownik stochastyczny.\n\
# Nie mamy jak losowac czestosci, poniewaz uniemozliwia to siatka, na ktorej sa zapisane atomy. Parametr periodDenstiy umozliwia\n\
# zageszczenie tej siatki. W przypadku slownika \"FIXED\", dla zadnej czestosci DS, czestosci atomow sa rozlozone w nastepujacy\n\
# sposob:\n\
# DF0   2DF0   3DF0   ...   PI\n\
# W przypadku slownika stochastycznego podstawowa czestosc to DF0' = DF0 / periodDensity i wtedy czestosci sa rozlozone gesciej:\n\
# DF0/periodDenstiy   2*DF0/periodDensity   ...   PI\n\
# Mozna zauwazyc, ze czestosci dokladnie wyliczone ze wzorow na Slownik Optymalny znajduja sie w pozycjach\n\
# k * periodDensity * DF, gdzie k = 0, ... , zas pomiedzy tymi czestosciami mamy periodDensity posrednich czestosci,\n\
# wsrod ktorych mozemy dokonac losowania. Im wieksze jest periodDensity, tym slownik bedzie bardziej stochastyczny,\n\
# ale zuzycie pamieci bedzie wieksze, a dokladnie periodDensity razy wieksze niz w przypadku slownika zafiksowanego.\n\
#\n\
# The current version of MP -- mp5 -- contains implementaion of the so-called Optimal Dictionary. This creates an\n\
# opportunity to perform certain numerical optimizations, e.g. with respect to storing of the arrays containing the values of sin/cos/exp functions.\n\
# Unfortunately, such arrays consume lots of memory, which often makes computation impossible. In order to decrease RAM consumption,\n\
# sin/cos values are packed in a special way. After getting acquainted with the optimal dictionary, one can note that for a given\n\
# scale DS, an optimal step/increment/delta in frequency, DF, is computed. For a given, required scale, atom frequencies should be\n\
# distributed every DF within the range (0 PI). Therefore, for a scale DS, the frequency of a atom is a multiple of a certain basic\n\
# frequency, later called DF0. Packing of sin/cos arrays in PC's memory is performed such that for an exemplary sin function\n\
# with 10-sample period, one can easily generate a sine with a period which is 9x, 8x, ... , 2x smaller,\n\
# by proper selection of the 10-sample sine. Similarly, the optimal dictionary contains frequencies which are multiples of a certain\n\
# basic frequency, and thus, their periods are inverse multiples of a certain basic period. For these reasons,\n\
# for a required scale DS, only sines of the basic FREQUENCY DF0 are stored in mp5 implementation. Unfortunately, certain problems\n\
# arise if one wants to have a stochastic dictionary. There is no way to randomize frequency, due to the grid on which\n\
# atoms are placed. periodDenstiy parameter enables making this grid denser. In the case of the FIXED dictionary, for a required frequency DS,\n\
# atom frequencies are distributed as follows:\n\
# DF0   2DF0   3DF0   ...   PI\n\
# In the case of the stochastic dictionary, basic frequency is equal to DF0' = DF0 / periodDensity, an then the frequencies are distributed denser:\n\
# DF0/periodDenstiy   2*DF0/periodDensity   ...   PI\n\
# One can note that the frequencies computed with the formulas for the Optimal Dictionary are located at the positions\n\
# k * periodDensity * DF, where k = 0, ... , while in between these frequencies we have periodDensity of intermediate frequencies\n\
# to randomize from.\n\
#\n\
\n\
periodDensity 1\n\
\n\
#\n\
# Reinicjalizacja slownika.\n\
# Dozwolone wartosci argumentu:\n\
# NO_REINIT_AT_ALL         -- Do dekompozcji sygnalow w poszczegolnych kanalach i skladkach\n\
#                             bedzie uzywany jeden, ten sam slownik.\n\
# REINIT_IN_CHANNEL_DOMAIN -- Dane beda analizowane kanalami. Wygnerowany slownik zostnie uzyty do dekompozycji\n\
#                             wszystkich skladek (offsetow) wchodzacych w sklad danego kanalu.\n\
#                             Przed dekompozycja sygnalu w nastepnym kanale zostanie wygnerowany nowy slownik.\n\
#                             Ta opcja jest dozwolona tylko w przypadku metody jednokanalowej (single matching pursuit -- SMP).\n\
# REINIT_IN_OFFSET_DOMAIN  -- Dane beda analizowane skladkami. Wygnerowany slownik zostnie uzyty do dekompozycji\n\
#                             wszystkich kanalow wchodzacych w sklad danej skladki (offsetu).\n\
#                             Przed dekompozycja sygnalu w nastepnym offsecie, zostanie wygnerowany nowy slownik.\n\
# REINIT_AT_ALL            -- Przed dekompozycja jakiegokolwiek fragmentu sygnalu (kanalu/skladki) bedzie generowany nowy slownik.\n\
#\n\
# Re-initialisation of the dictionary.\n\
# Argument's values:\n\
# NO_REINIT_AT_ALL         -- One and same dictionary for decompositions in all channels and all epochs.\n\
# REINIT_IN_CHANNEL_DOMAIN -- Data will be analysed channel by channel. A generated dictionary will be used to decompose all epochs\n\
#                             within a channel. Before the analysis of another channel, a new dictionary will be generated.\n\
#                             This option in only allowed in the case of the single-channel matching pursuit -- SMP.\n\
# REINIT_IN_OFFSET_DOMAIN  -- Data will be analysed epoch by epoch. A generated dictionary will be used to decompose all channels\n\
#                             within an epoch. Before the analysis of another epoch, new dictionary will be generated.\n\
# REINIT_AT_ALL            -- Before the analysis of any piece of data (channel/epoch), a new dictionary will be generated.\n\
#\n\
\n\
reinitDictionary NO_REINIT_AT_ALL\n\
\n\
#\n\
# Parametr umozliwiajacy budowe slownika oscylacyjnego.\n\
# Wszystkie atomy, dla ktorych w ramach danej skali krotnosc okresu oscylacji okreslona przez ten parametr sie nie miesci,\n\
# zostana odrzucone ze slownika.\n\
#\n\
# The parameter which enables building of an oscillatory dictionary.\n\
# All atoms for which the number of oscillatory periods, given by this parameter, exceeds the length of a scale\n\
# will be removed from the dictionary.\n\
#\n\
\n\
scaleToPeriodFactor 0.0\n\
\n\
#\n\
# Maksymalna liczba atomow, ktora moze byc uzyta do rozlozenia sygnalu.\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
#\n\
# Maximal number of atoms which can be used to decompose a signal.\n\
# Argument's value: positive integer.\n\
#\n\
\n\
maximalNumberOfIterations 100\n\
\n\
#\n\
# Zatrzymaj, jesli dokompozycja opisala conajmniej x%  energii sygnalu.\n\
# Wartosc argumentu: liczba zmiennoprzecinkowa, dodatnia.\n\
# Obowiazkowo w wartosci musi wystapic kropka dziesietna.\n\
#\n\
# Stop if the decomposition explained at least x% of the signal's energy.\n\
# Argument's value: positive floating-point number.\n\
# Decimal dot is obligatory.\n\
#\n\
\n\
energyPercent 99.0\n\
\n\
#\n\
# Rodzaj algorytmu MP.\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
# Do wyboru sa nastepujace algorytmy (zobacz Piotr Durka \"Matching pursuit and unification in EEG analysis\"):\n\
#\n\
# MP strategy.\n\
# Argument's value: positive integer.\n\
# The following algorithms are available (see Piotr Durka \"Matching pursuit and unification in EEG analysis\"):\n\
# SMP  -- SINGLE CHANNEL MATCHING PURSUIT.\n\
# MMP1 -- MULTICHANNEL MATCHING PURSUIT (constant phase, maximum sum of energies).\n\
# MMP2 -- MULTICHANNEL MATCHING PURSUIT (constant phase, maximum sum of products).\n\
# MMP3 -- MULTICHANNEL MATCHING PURSUIT (variable phase, maximum sum of energies).\n\
#\n\
\n\
MP MMP2\n\
\n\
#\n\
# Czy do plikow z ksiazkami zawierajacymi wynik obliczen (dopasowane atomy) \n\
# dolaczac dane, na ktorych wykonywano obliczenia?.\n\
# Wartosc argumentu: YES/NO.\n\
#\n\
# Eanglish version.\n\
#\n\
\n\
analiticalDotProduct OFF\n\
\n\
#\n\
\n\
bookWithSignal NO\n\
\n\
#\n\
# Stala konwersji probka -> napiecie w uV.\n\
# Obowiazkowo w wartosci musi wystapic kropka dziesietna.\n\
#\n\
# Sample -> voltage [uV] conversion constant.\n\
# Decimal dot is obligatory.\n\
#\n\
\n\
pointsPerMicrovolt 100.0\n\
\n"};

		if(exampleFile == NULL)
		{
		    fprintf(stderr," \n");
		    fprintf(stderr," ERROR: \n");
		    fprintf(stderr," Can't open file: example.set, for default parameters");
		    fprintf(stderr," \n");
		    return 1;
		}

		fprintf(exampleFile,"%s\n",exampleTest);
		fclose(exampleFile);
    
		printf(" Defaults parameters are being written to example.set file\n");

		return 0;
    }
	else if(applicationMode & ERROR_MODE)
	{
		if(generateErrotTextAndCodeGeneratorFile()==ERROR)
		{
		    fprintf(stderr," \n");
		    fprintf(stderr," ERROR: Can't write errors codes and text to file: errorTranslation.xml");
		    fprintf(stderr," \n");
			return 1;
		}
		else
		{
			printf(" Errors texts and codes of errors has been written to file: errorsCodes.xml\n");
			return 0;
		}
	}
    else if((applicationMode & TEST_PARAMETERS_MODE ) || (applicationMode & PROCESS_USER_MODE) || (applicationMode & PROCESS_SERVER_MODE))
    {
		/* alocate memory for components of configFile */

		setConfigFile(&configFile);
		if(strlen(argv[2])>LENGTH_OF_NAME_OF_CONFIG_FILE)
		{
			fprintf(stderr," ERROR: Too long name of config file %s. The maximum length is %hu\n",argv[2],LENGTH_OF_NAME_OF_CONFIG_FILE);
			return 1;
		}

		strcpy(configFile.name,argv[2]);

		/* try to open config file */

		if(openConfigFile(&configFile,infoMessage)==ERROR)
		{
			fprintf(stderr,"%s",infoMessage);
			return 1;
		}

		/* read config file */
		if(readConfigFile(&configFile,infoMessage)==ERROR) goto ERROR_PROCEDURE;

		/* analyse line and command in config file
		"primary" test of config file and parameters included in this file */
		if(findDataParametersInConfigFile(&configFile,&dictionary,&mp5Parameters,infoMessage)==ERROR)
			goto ERROR_PROCEDURE;
		
		/* "extended" test of parameters read in config file */
		if(testMP5Parameters(&dictionary,&mp5Parameters,infoMessage)==ERROR)
			goto ERROR_PROCEDURE;

		/* try to open data */
		if(mp5Parameters.dataFormat & FORMAT_ASCII)
		{
			if(openAsciiDataFile(&mp5Parameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;
		}
		else if((mp5Parameters.dataFormat & FORMAT_FLOAT) || (mp5Parameters.dataFormat & FORMAT_SHORT))
		{
			if(openBinaryDataFile(&mp5Parameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;
		}

		/* "primary" test of data file */
		if(mp5Parameters.dataFormat & FORMAT_ASCII)
		{
			if(analyseAsciiDataFile(&mp5Parameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;
		}
		else if((mp5Parameters.dataFormat & FORMAT_FLOAT) || (mp5Parameters.dataFormat & FORMAT_SHORT))
		{
			if(analyseBinaryDataFile(&mp5Parameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;
		}

        /* set some constans values for MP5Parameters and DataParameters */
        setNumberOfAnalysedChannelsAndNumberOfResultsFiles(&mp5Parameters);

        /* check dictionary tape and size */
        analyseDictionarySizeAndType(&dictionary,&mp5Parameters);

        /* create names of results files */
        createNamesOfResultFiles(&dictionary,&mp5Parameters);

		/* check wheter directories and results files are existing */
		if(testFilesAndDirectories(&mp5Parameters,&configFile,infoMessage) == ERROR) goto ERROR_PROCEDURE;

		if(applicationMode & PROCESS_USER_MODE)
			printInfoAboutData(&dictionary,&mp5Parameters);

		if(mp5Parameters.numberOfThreads>1)
			fftw_init_threads();				
			
        if((applicationMode & PROCESS_USER_MODE) || (applicationMode & PROCESS_SERVER_MODE))
        {
    	    /* now allocate memory for components of Dictionary */
    	    allocateDictionary(&dictionary,&mp5Parameters);
  	    
			/* now allocate memory for components of MP5Parameters */
			setMP5Parameters(&dictionary,&mp5Parameters);

			/* make sin/cos and exp tables */
			makeSinCosExpTable(&dictionary,&mp5Parameters);

			/* results file */
			if(openResultFiles(&mp5Parameters,infoMessage)==ERROR)
				goto ERROR_PROCEDURE;

			unsigned short int offsetNumber  = 0;
			unsigned short int channelNumber = 0;

			/* results file is open, header is written into results file */
			{
				unsigned long int currentPositionInResultsFile = ftell(mp5Parameters.resultsFile);
				fseek(mp5Parameters.resultsFile,0L,SEEK_END);
				unsigned long int sizeOfResultsFile = ftell(mp5Parameters.resultsFile);
				fseek(mp5Parameters.resultsFile,currentPositionInResultsFile,SEEK_SET); 

				if((mp5Parameters.writingMode & CREATE_FILE) || (sizeOfResultsFile==0))
				{
					if(writeMagic(&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					if(writeCommentsSegment(&configFile,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;

					if(writeFileHeader(&dictionary,&mp5Parameters,infoMessage)==ERROR)
						goto ERROR_PROCEDURE;
				}
			}
		
			if(mp5Parameters.MPType & SMP)
			{
				if(applicationMode & PROCESS_SERVER_MODE)
				{
					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfChosenOffsets,
															mp5Parameters.numberOfChosenChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);
					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);

					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						for(channelNumber=0;channelNumber<mp5Parameters.numberOfChosenChannels;channelNumber++)
						{	
							if(applicationMode & PROCESS_SERVER_MODE)
							{
								printf("CHANNEL %hu\n",channelNumber);
								fflush(stdout);
							}
							else if(applicationMode & PROCESS_USER_MODE)
							{
								printf("\n --OFFSET--: %d, |CHANNEL|: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber],mp5Parameters.chosenChannels[channelNumber]);
								fflush(stdout);
							}
						
							mp5Parameters.singleChannelSignalTable = *(mp5Parameters.processedDataMatrix + mp5Parameters.chosenChannels[channelNumber]-1);
							firstIterationSMPMMP2(&dictionary,&mp5Parameters);
							nextIterationSMPMMP2(&dictionary,&mp5Parameters);

							if(writeSMPResults(&dictionary,&mp5Parameters,offsetNumber,channelNumber,infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							resetDictionary(&dictionary);
						}
					}      
				}
				else if(mp5Parameters.reinitDictionary & REINIT_IN_CHANNEL_DOMAIN)
				{
					unsigned short int counter;
					time_t       seedTime;
					long int     seed[mp5Parameters.numberOfChosenChannels];
					
					for(counter=0;counter<mp5Parameters.numberOfChosenChannels;counter++)
						if(dictionary.randomSeed==AUTO_RANDOM_SEED)
							seed[counter] = time(&seedTime);
						else
							seed[counter] = dictionary.randomSeed;

					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)					
					{	
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						for(channelNumber=0;channelNumber<mp5Parameters.numberOfChosenChannels;channelNumber++)
						{
							dictionary.randomSeed = seed[channelNumber];

							/* create dicionary */
							reinitDictionary(&dictionary,&mp5Parameters);

							/* test atom's feature, for example find INCORRECT atoms */
//							testAtomFeature(&dictionary);

							if(applicationMode & PROCESS_SERVER_MODE)
							{
								printf("CHANNEL  %hu\n",channelNumber);
								fflush(stdout);
							}

							if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							if(applicationMode & PROCESS_USER_MODE)
							{
								printf("\n |CHANNEL|: %d, --OFFSET--: %d\n\n",mp5Parameters.chosenChannels[channelNumber],mp5Parameters.chosenOffsets[offsetNumber]);
								fflush(stdout);
							}

							mp5Parameters.singleChannelSignalTable = *(mp5Parameters.processedDataMatrix + mp5Parameters.chosenChannels[channelNumber]-1);

							firstIterationSMPMMP2(&dictionary,&mp5Parameters);
							nextIterationSMPMMP2(&dictionary,&mp5Parameters);

							if(writeSMPResults(&dictionary,&mp5Parameters,offsetNumber,channelNumber,infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							resetDictionary(&dictionary);
						}
					}
				}
				else if(mp5Parameters.reinitDictionary & REINIT_IN_OFFSET_DOMAIN)
				{
					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						for(channelNumber=0;channelNumber<mp5Parameters.numberOfChosenChannels;channelNumber++)
						{
							if(applicationMode & PROCESS_SERVER_MODE)
							{
								printf("CHANNEL  %hu\n",channelNumber);
								fflush(stdout);
							}

							if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							if(applicationMode & PROCESS_USER_MODE)
							{
								printf("\n |CHANNEL|: %d, --OFFSET--: %d\n\n",mp5Parameters.chosenChannels[channelNumber],mp5Parameters.chosenOffsets[offsetNumber]);
								fflush(stdout);
							}

							mp5Parameters.singleChannelSignalTable = *(mp5Parameters.processedDataMatrix + mp5Parameters.chosenChannels[channelNumber]-1);

							firstIterationSMPMMP2(&dictionary,&mp5Parameters);
							nextIterationSMPMMP2(&dictionary,&mp5Parameters);
 
							if(writeSMPResults(&dictionary,&mp5Parameters,offsetNumber,channelNumber,infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							resetDictionary(&dictionary);
						}
					}
				}
				else if(mp5Parameters.reinitDictionary & REINIT_AT_ALL)
				{
					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}
	
						for(channelNumber=0;channelNumber<mp5Parameters.numberOfChosenChannels;channelNumber++)
						{
							if(applicationMode & PROCESS_SERVER_MODE)
							{
								printf("CHANNEL  %hu\n",channelNumber);
								fflush(stdout);
							}
	
							/* create dicionary */
							reinitDictionary(&dictionary,&mp5Parameters);

							/* test atom's feature, for example find INCORRECT atoms */
//							testAtomFeature(&dictionary);

							if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
								goto ERROR_PROCEDURE;

							if(applicationMode & PROCESS_USER_MODE)
							{
								printf("\n --OFFSET--: %d, |CHANNEL|: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber],mp5Parameters.chosenChannels[channelNumber]);
								fflush(stdout);
							}

							mp5Parameters.singleChannelSignalTable = *(mp5Parameters.processedDataMatrix + mp5Parameters.chosenChannels[channelNumber]-1);
		
							firstIterationSMPMMP2(&dictionary,&mp5Parameters);
							nextIterationSMPMMP2(&dictionary,&mp5Parameters);

							if(writeSMPResults(&dictionary,&mp5Parameters,offsetNumber,channelNumber,infoMessage)==ERROR)
								goto ERROR_PROCEDURE;
						}
					}
				}
			}
			// writen by Artur Matysiak
			else if(mp5Parameters.MPType & MMP1)
			{
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfChosenOffsets,
															mp5Parameters.numberOfChosenChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}
														
				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);
                
					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --OFFSET--: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						firstIterationMMP(&dictionary,&mp5Parameters);
						nextIterationMMP(&dictionary,&mp5Parameters);

						if(writeMMPResults(&dictionary,&mp5Parameters,offsetNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
				if(mp5Parameters.reinitDictionary &  REINIT_IN_OFFSET_DOMAIN)
				{

					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --OFFSET--: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;
		
						firstIterationMMP(&dictionary,&mp5Parameters);
						nextIterationMMP(&dictionary,&mp5Parameters);

						if(writeMMPResults(&dictionary,&mp5Parameters,offsetNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
		    }    
			//Artur Matysiak end
		    else if(mp5Parameters.MPType & MMP2)
		    {
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfChosenOffsets,
															mp5Parameters.numberOfChosenChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);
	                
					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --OFFSET--: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						countMeanSignalOverChannels(&mp5Parameters);
						mp5Parameters.singleChannelSignalTable = mp5Parameters.meanSignalTable;
				
						firstIterationSMPMMP2(&dictionary,&mp5Parameters);
						nextIterationSMPMMP2(&dictionary,&mp5Parameters);

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
						    goto ERROR_PROCEDURE;

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;
					
						if(writeMMPResults(&dictionary,&mp5Parameters,offsetNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
				if(mp5Parameters.reinitDictionary &  REINIT_IN_OFFSET_DOMAIN)
				{

					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{

						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}
	
						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --OFFSET--: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						countMeanSignalOverChannels(&mp5Parameters);
						mp5Parameters.singleChannelSignalTable = mp5Parameters.meanSignalTable;

						firstIterationSMPMMP2(&dictionary,&mp5Parameters);
						nextIterationSMPMMP2(&dictionary,&mp5Parameters);

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
						    goto ERROR_PROCEDURE;

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;
		
						if(writeMMPResults(&dictionary,&mp5Parameters,offsetNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
		    }
		    else if(mp5Parameters.MPType & MMP3)
		    {
					if(applicationMode & PROCESS_SERVER_MODE)
					{

						printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfChosenOffsets,
																mp5Parameters.numberOfChosenChannels,
																mp5Parameters.maximalNumberOfIterations,
																mp5Parameters.energyPercent);
						fflush(stdout);
					}

														
				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);
	                
					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --OFFSET--: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						firstIterationMMP(&dictionary,&mp5Parameters);
						nextIterationMMP(&dictionary,&mp5Parameters);

						if(writeMMPResults(&dictionary,&mp5Parameters,offsetNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
				if(mp5Parameters.reinitDictionary &  REINIT_IN_OFFSET_DOMAIN)
				{

					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --OFFSET--: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;

						firstIterationMMP(&dictionary,&mp5Parameters);
						nextIterationMMP(&dictionary,&mp5Parameters);

						if(writeMMPResults(&dictionary,&mp5Parameters,offsetNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
			}
			else if(mp5Parameters.MPType & MMP4)
		    {
				if(applicationMode & PROCESS_SERVER_MODE)
				{

					printf("\nSTART\t%hu\t%hu\t%u\t%6.2f\n",mp5Parameters.numberOfChosenOffsets,
															mp5Parameters.numberOfChosenChannels,
															mp5Parameters.maximalNumberOfIterations,
															mp5Parameters.energyPercent);
					fflush(stdout);
				}

				if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
				{
					/* create dicionary */
					makeDictionary(&dictionary,&mp5Parameters);

					/* test atom's feature, for example find INCORRECT atoms */
//					testAtomFeature(&dictionary);
	                
					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{
						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --OFFSET--: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;			
						firstAndNextIterationMMP4(&dictionary,&mp5Parameters);
					
						if(writeMMPResults(&dictionary,&mp5Parameters,offsetNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
				if(mp5Parameters.reinitDictionary &  REINIT_IN_OFFSET_DOMAIN)
				{
					for(offsetNumber=0;offsetNumber<mp5Parameters.numberOfChosenOffsets;offsetNumber++)
					{

						if(applicationMode & PROCESS_SERVER_MODE)
						{
							printf("OFFSET  %hu\n",offsetNumber);
							fflush(stdout);
						}
	
						/* create dicionary */
						reinitDictionary(&dictionary,&mp5Parameters);

						/* test atom's feature, for example find INCORRECT atoms */
//						testAtomFeature(&dictionary);

						if(readDataFile(&mp5Parameters,mp5Parameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						if(applicationMode & PROCESS_USER_MODE)
						{
							printf("\n --OFFSET--: %d\n\n",mp5Parameters.chosenOffsets[offsetNumber]);
							fflush(stdout);
						}

						mp5Parameters.multiChannelSignalTable = mp5Parameters.processedDataMatrix;
						firstAndNextIterationMMP4(&dictionary,&mp5Parameters);
		
						if(writeMMPResults(&dictionary,&mp5Parameters,offsetNumber,infoMessage)==ERROR)
							goto ERROR_PROCEDURE;

						resetDictionary(&dictionary);
					}
				}
		    }
		} 
    }
	
    freeConfigFile(&configFile);
    closeFiles(&mp5Parameters);
    freeMP5Parameters(&dictionary,&mp5Parameters);
    freeDictionary(&dictionary);
	if(mp5Parameters.numberOfThreads>1)
		fftw_cleanup_threads();
	printf(" END\n");
	fflush(stdout);

    return 0;

    ERROR_PROCEDURE:
	fprintf(stderr," \n");
	fprintf(stderr," %s",infoMessage);
	fprintf(stderr,"\n");
	freeConfigFile(&configFile);
	closeFiles(&mp5Parameters);
	freeMP5Parameters(&dictionary,&mp5Parameters);
    freeDictionary(&dictionary);
	if(mp5Parameters.numberOfThreads>1)
		fftw_cleanup_threads();
	fflush(stdout);

    return 1;
}
