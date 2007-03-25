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

#include<math.h>
#include<stdio.h>
#include"include/cmd.h"
#include"include/def.h"
#include"include/dic.h"
#include"include/io.h"
#include"include/mmp1.h"
#include"include/mmp3.h"
#include"include/mp5.h"
#include"include/smp.h"
#include"include/stringTools.h"
#include"include/types.h"

#ifdef __MINGW32__
#define bzero(ptr,size) memset (ptr, 0, size);
#define sincos(th,x,y) { (*(x))=sin(th); (*(y))=cos(th); }
#endif

static void countMeanSignalOverChannels(MP5Parameters *mp5Parameters,DataParameters *dataParameters )
{
	unsigned short int sample;
	double tmpDataValue;
	
	double **multiChannelSignalTable = mp5Parameters->multiChannelSignalTable;
	double *meanSignalTable         = mp5Parameters->meanSignalTable;
	int channel;
	for(sample=0;sample<mp5Parameters->dimOffset;sample++)
	{
		tmpDataValue = 0.0;
		
		for(channel=0;channel<dataParameters->numberOfChosenChannels;channel++)
			tmpDataValue+= *(*(multiChannelSignalTable + channel) + sample);

		*(meanSignalTable + sample) = tmpDataValue/mp5Parameters->numberOfAnalysedChannels;
		mp5Parameters->meanSignalEnergy+= (*(meanSignalTable + sample))*(*(meanSignalTable + sample));
	}
}

int main(int argc, char *argv[])
{
    unsigned char mode;

    ConfigFile configFile = {"",
			     NULL,
			     NULL,
			     NULL};

  
    DataParameters dataParameters = {"","","",
				     0,
				     NULL,
				     "","",
				     NULL,NULL,NULL,NULL,
				     0,0,
				     0,0,0,0,
				     0.0,
				     0,NULL,
				     0,NULL,
				     0,
				     0x0,
				     0.0,
				     0,0,
				     NULL,NULL,
				     0x0,0x0,0x0};
 
    MP5Parameters  mp5Parameters = {0,0,0,0,
				    NULL,NULL,NULL,
				    NULL,NULL,
				    NULL,NULL,
				    NULL,NULL,NULL,NULL,NULL,NULL,
				    0.0,0.0,0.0,0.0,
				    NULL,NULL,
				    0,0.0,
				    NULL,
				    0x0,0x0,
				    NULL,
				    0.0,
				    FALSE};

    GaborDictionary gaborDictionary = {0,0,
				       0.0,0.0,0,
				       0,
				       0.0,0.0,
				       0.0,
				       0,
				       0.0,
				       NULL,NULL,NULL,NULL,
				       NULL,NULL,
				       0,NULL,
				       0x0};
		     
    char infoMessage[LENGTH_OF_INFO_MESSAGE];

    printf("\n");
    printf("  Matching Pursuit V (version 2006-08-03)\n");
    printf("  Department of Biomedical Physics at Warsaw University\n");
    printf("  http://brain.fuw.edu.pl, http://eeg.pl\n");
    printf("  Compiled: %s, %s \n", __DATE__ , __TIME__ );

    #ifdef __VERSION__
	printf("  Compiler: %s \n", __VERSION__);
    #endif

    printf("\n");

    if(argc>3 || argc==1 || (argc==2 && (strcmp(argv[1],"-g")!=0)))
    {
	fprintf(stderr," \n");
	fprintf(stderr," ERROR: \n");
	fprintf(stderr," INCORRECT CALL OF mp5Parameters PROGRAM \n");
	fprintf(stderr," THE PROPER USE IS AS FOLLOWS: \n");
	fprintf(stderr," mp5 -g                         - generate default config file \n");
	fprintf(stderr," mp5 -t [ name of config file ] - test config file \n");
	fprintf(stderr," mp5 -f [ name of config file ] - process of mp decomposition \n");
	fprintf(stderr," \n");
	return ERROR;
    }

    if(strcmp(argv[1],"-g")==0)
	mode = GENERATE;
    else if(strcmp(argv[1],"-t")==0)
	mode = TEST;
    else if(strcmp(argv[1],"-f")==0)
	mode = PROCESS;
    else
    {
	fprintf(stderr," \n");
	fprintf(stderr," ERROR: \n");
	fprintf(stderr," INCORRECT CALL OF mp5Parameters PROGRAM \n");
	fprintf(stderr," UNKNOW PARAMETER: %s \n",argv[1]);
	fprintf(stderr," THE PROPER USE IS AS FOLLOWS: \n");
	fprintf(stderr," mp5 -g                         - generate default config file \n");
	fprintf(stderr," mp5 -t [ name of config file ] - test config file \n");
	fprintf(stderr," mp5 -f [ name of config file ] - process of mp decomposition \n");
	fprintf(stderr," \n");
	return ERROR;
    }

    if(mode & GENERATE)
    {
	FILE *exampleFile = fopen("example.set","wt");
	char exampleTest[] =
{"#\n\
# PLIK KONFIGURACYJNY PROGRAMU mp5 - SPECYFIKACJA\n\
#\n\
# Plik konfiguracyjny sklada sie z argumentow i ich wartosci oraz komentarzy\n\
# Komentarzem jest linia rozpoczynajaca sie od znaku '#'. Komentarzy moze byc dowolna ilosc.\n\
# Lista argumentow jest scisle okreslona. Brak lub dopisanie jakiegos argumentu\n\
# powoduje przerwanie dzialalani programu mp5Parameters. \n\
# Obecnie program mp5Parameters rozpoznaje 14 nastepujacych argumentow: \n\
#\n\
#	nameOfDataFile\n\
#	extensionOfResultFile\n\
#       nameOfOutputDirectory\n\
#	writingMode\n\
#	sizeOfHeader\n\
#	sizeOfTail\n\
#	samplingFrequency\n\
#	formatOfData\n\
#	numberOfChannels\n\
#	chosenChannels\n\
#	numberOfPointsInOffset\n\
#	chosenOffsets\n\
#	typeOfDictionary\n\
#	dilationFactor\n\
#	periodDensity\n\
#	reinitDictionary\n\
#       scaleToPeriodFactor\n\
#       DOT_EPS\n\
#	maxNumberOfIterations\n\
#	energyPercent\n\
#	MP\n\
#	convRate\n\
#	VERBOSE\n\
#\n\
# Argumenty oraz ich wartosci moga byc podane w dowolnej kolejnosci.\n\
#\n\
# W skladni pliku konfiguracyjnego wyrozniono nastepujace znaki specjalne:\n\
#\n\
# 1. '#' - komentarz\n\
# 2. '\\' - lamanie linii\n\
# 3. '-' - zakres wartosci\n\
# 4. ' ' - znak przestankowy\n\
#\n\
# Adn. 1 '#'\n\
# Wszytkie linie zaczynajace sie od '#' sa traktowane jak komentarze\n\
# i nie sa poddawane dalszej analizie.\n\
# Komentarzy moze byc w pliku konfiguracyjmym dowolna ilosc\n\
# Mozna je takze na wlasny uzytek dopisywac. Aby jednak uchronic\n\
# uzytkownikow przed pozniejszym krazeniem po sieci plikow konfiguracujnych\n\
# z roznymi wersjami komentarzy, zawsze mozna wygenerowac domyslny\n\
# plik konfiguracyjny wydajac polecenie mp5 -g.\n\
#\n\
# Adn. 2 '\\'\n\
# Polecenia wydawane dla programu mp5Parameters moga byc dosyc dlugie,\n\
# o czym przekonamy sie dalej. Aby nie pisac nieskonczenie dlugich linii\n\
# dopusczono ich lamanie tak samo jak w jezyku C.\n\
#\n\
# And. 3 '-'\n\
# Analizie programem mp5Parameters moga byc poddane dane zawierajace wiele kanalow i skladek.\n\
# Kanaly lub skladki, ktore maja byc analizowane przez program nalezy wymienic w pliku konfiguracyjnym.\n\
# W przypadku gdy lista kanalow lub skladek bylaby dluga, aby uniknac zmudnego ich wypisywania\n\
# mozna podac zakres kanalow 'od-do', np. 1-23.\n\
# UWAGA I:\n\
#	Miedzy liczbami nie moze wystapic spacja !!!\n\
# UWAGA II.:\n\
#	Lamanie linii ('-') oraz uzywanie znaku zakres ('\') jest dopuszczone tylko dla argumentow\n\
# 	chosenChannels i chosenOffsets\n\
# 	Wartosci wszystkich pozostalych parametrow powinny byc jednoargumentowe !!\n\
#\n\
\n\
\n\
#\n\
# Nazwa pliku z danymi.\n\
#\n\
\n\
nameOfDataFile transien.dat\n\
\n\
#\n\
# Nazwa rozszerzenia pliku z wynikami.\n\
# Nazwy plikow wyjsciowych sa generowane automatycznie.\n\
# Do nazwy pliku z danymi dodawane jest rozszerzenie\n\
# postaci _ch_x.b, gdzie x - to numer kanalu, wymieniony dalej w pliku konfiguracyjnym\n\
# Pelna nazwa pliku wyjsciowego jest zatem nastepujaca:\n\
# nameOfDataFile_ch_x.b\n\
# W przypadku, gdy taki sposob tworzenia nazwy pliku wyjsciowego jest niewystarczajacy\n\
# mozna dodac do nazwy wlasne rozszerzenie. W takim przypadku pelna nazwa pliku wyjsciowego\n\
# jest postaci: \n\
# nameOfDataFile_extensionOfResultFile_ch_x.b.\n\
# Jezeli wartosc argumentu jest rowna NONE, nazwa pliku wyjsciowego jest generowana\n\
# bez rzadnego dodatkowego rozszerzenia.\n\
#\n\
\n\
extensionOfResultFile NONE\n\
\n\
#\n\
# Mozna podac nazwe katalogu, do ktorego maja byc zrzucane wyniki\n\
# Nazea katalogu musie sie konczyc znakiem '/'\n\
#\n\
\n\
nameOfOutputDirectory ./\n\
\n\\
#\n\
# Tryb zapisu wynikow.\n\
# Dowzolone sa nastepujace wartosci:\n\
# CREATE - utworz nowy plik do zapisu wynikow. Jesli juz plik o takiej nazwie istnieje,\n\
#          to go zamaz.\n\
# APPEND - dolacz wyniki do istniejacego juz pliku/utworz nowy plik, jesli plik o danej nazwie\n\
#\n\
\n\
writingMode CREATE\n\
\n\
#\n\
# Wielkosc naglowka w pliku z danymi.\n\
# W przypadku, gdy plik z danymi jest plikiem tesktowym\n\
# sizeOfHeader oznacza liczbe linii w naglowku.\n\
# W przypadku gdy dane sa zapisane w pliku binarnym, rozmiar naglowka\n\
# powinna byc wyrazona w bajtach.\n\
# Wartosc argumentu: liczba calkowita dodatnia lub zero.\n\
#\n\
\n\
sizeOfHeader 0\n\
\n\
#\n\
# Czasami do danych EEG dolaczany jest jeszcze \"ogon\".\n\
# Jego rozmiar powinien byc podany w taki sam sposob jak zostalo okreslone \n\
# dla rozmiaru naglowka\n\
# Wartosc argumentu: liczba calkowita dodatnia lub zero.\n\
#\n\
\n\
sizeOfTail 0\n\
\n\
#\n\
# Czestosc probkowania.\n\
# Wartosc argumentu: liczba zmiennoprzecinkowa, dodatnia.\n\
# Obowiazkowo w wartosci musi wystapic kropka dziesietna.\n\
#\n\
\n\
samplingFrequency 128.0\n\
\n\
#\n\
# Format danych.\n\
# W przypadku formatu ASCII zaklada sie, ze w pliku nie ma pustych linii\n\
# Wartosc argumentu: ASCII/SHORT/FLOAT\n\
#\n\
\n\
formatOfData FLOAT\n\
\n\
#\n\
# Liczba kanalow w pliku.\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
#\n\
\n\
numberOfChannels 2\n\
\n\
#\n\
# Kanaly, ktore maja zostac poddane analizie.\n\
# Liczby calkowite dodatanie, oddzielone przecinkami lub laczone znakiem '-'.\n\
# Kanaly sa numerowane od 1.\n\
# Przyklad skladni:\n\
# Plik zawiera 128 kanalow. Analizie maja zostac poddane nastepujace kanaly\n\
# 1 2 3 4 5 6 7 8 9 25 27 29 31 33 35 37 45 65 66 67 121 122 123 124 125 126 127 128\n\
# Mozna je wymienic w pliku w nastepujacy sposob:\n\
#\n\
# chosenChannels 1-9 25 27 29 31 \\\n\
#		  33 35 37 45 65-67 \\\n\
#		  121-128\n\
#\n\
\n\
chosenChannels 1 2\n\
\n\
#\n\
# Dlugosc skladki sygnalu, mierzona w ilosci probek, ktora\n\
# jednorazowo ma zostac poddana analizie.\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
#\n\
\n\
numberOfPointsInOffset 512\n\
\n\
#\n\
# Numery skladek, ktore chcemy analizowac.\n\
# Nomenklatura taka sama jak przy chosenChannels\n\
#\n\
\n\
chosenOffsets 1\n\
\n\
#\n\
# Rodzaj slownika\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
#\n\
# OCTAVE_FIXED: Slownik optymalizowany. Wartosci parametery gaborow {u,f,s} sa obliczane w taki sposob, by  \n\
#                   zbudowany slownik jak najefektywniej pokrywal przestrzen parametrow Gabora \n\
# OCTAVE_STOCH  Slownik optymalny obliczany taj jak wyzej, z ta roznica, ze parametry {u,f,s} sa losowane wokol \n\
#		    obliczonych optymalnych wartosci.\n\
#\n\
\n\
typeOfDictionary OCTAVE_FIXED\n\
\n\
#\n\
#\"Gestosc slownika\" \n\
# Wartosc argumentu: liczba zmiennoprzecinkowa dodatnia.\n\
#\n\
\n\
dilationFactor 2.0\n\
\n\
#\n\
# Obecna wersja programu MP - MP5 ma zaprogramowany Slownik Optymalny. Stwarza to mozliwosc dokonania pewnych optymalizacji\n\ 
# numerycznych m.in zapamietywania talbic z wartosciami funkcji sin/cos/exp. Niestety takie tablice zabieraja bardzo duza\n\
# ilosc pamieci, co czesto uniemozliwia dokonywanie obliczen. Aby zmniejszyc zapotrzebowanie na RAM, postanowiono zapaktywac\n\
# wartosci sin/cos w specjalny sposob. Po zapozmniu sie ze slownikiem optymalnym mozna zauwazyc, ze dla danej sakli DS\n\
# zostaje obliczony krok w czestosci DF. Gabory dla zadanej skali powinny miec czestosci rozlozone co DF w zakrecie (0 PI).\n\
# w zwiazku z tym dla skali DS czestosc gabora jest wielokrotnoscia pewnej podstawoehj czestosci, dalej nazywanej DF0.\n\
# Pakowanie tablic sin/cos w pamieci komputera polega na tym, ze majac stablicowana funkcje sin o okresie np. 10 probek\n\
# mozna poprzez wybieranie z opdowiednim skokiem tych problek szybko wygenerowac funkcje sin o okresach 9x, 8x, .... 2x\n\
# mniejszych !. Podobna sytucacje mamy w slowniku optymalnym, gdzie czestosci sa wielokrotnoscia pewnej czestosci bazowej\n\
# a zatem ich okresy sa dzielnikami pewnego okresu podstawowego. W programie MP5 zapoamietywane sa zatem dla zadanej skali DS\n\
# tylko sinusy o podstawowym okresie DF0. Niestety powstaja pewne problemy jezeli chcemy miec slownik stochastyczny.\n\
# Nie mamy jak losowac czestosci, poniewaz unimowliwia to siatka na ktorej sa zapisane gabory. Parametr periodDenstiy umozliwia\n\
# zageszczenie tej siatki. W przypadku slownika \"FIXED\" dla zadnej czestosci DS czestosci gaborow sa rozlozone w nastepujacy\n\
# sposob:\n\
# DF0 2DF0 3DF0 ... PI\n\
# W przypadku slownika stochastycznego podstawowy czestosc to DF0' = DF0/periodDensity i wtedy czestosci sa rozlozone gesciej:\n\
# DF0/periodDenstiy 2*DF0/periodDensity ... PI\n\
# Mozna zauwazyc, ze czestosc dokladnie wyliczone ze wzorow na Optymalny Slownik znajduja sie w pozycjach\n\
# k*periodDensity*DF, gdzie k = 0 ..., zas pomiedzy tymi czestosciami mamy periodDensity posrednich czestosci\n\
# wsrod ktorych mozemy dokonac losowania. Im wieksze jest periodDensity, tym slosnik bedzie bardziej stochastyczny,\n\
# ale zuzycie pamieci bedzie wieksze, a dokladnie periodDensity razy wieksze niz w przypadku slownika \"FIXED\"\n\
#\n\
\n\
periodDensity 1\n\
\n\
#\n\
# Reinicjalizacja slownika.\n\
# Dozwolone wartosci argumentu:\n\
#\n\
# NO_REINIT_AT_ALL         - do dekompozcji sygnalow w poszczegolnych kanalach i skladkach\n\
#                            bedzie uzywany jeden, ten sam slownik\n\
# REINIT_IN_CHANNEL_DOMAIN - dane beda analizowane kanalami. Wygnerowany slownik zostnie uzyty do dekompozycji\n\
#			     wszystkich skladek (offsetow) wchodzacych w sklad danego kanalu.\n\
#			     Przed dekompozycja sygnalu w nastepnym kanale, zostanie wygnerowany nowy slownik.\n\
#			     Ta opcja jest dozwolona tylko w przypadku metody jednokanalowej (single matching pursuit - SMP)\n\
# REINIT_IN_OFFSET_DOMAIN  - dane beda analizowane skladkami. Wygnerowany slownik zostnie uzyty do dekompozycji\n\
#			     wszystkich kanalow wchodzacych w sklad danej skladki (offsetu).\n\
#			     Przed dekompozycja sygnalu w nastepnym offsecie, zostanie wygnerowany nowy slownik.\n\
# REINIT_AT_ALL            - Przed dekompozycja jakiegokolwiek fragmentu sygnalu (kanalu/skladki) bedzie generowany nowy\n\
#			     slownik.\n\
#\n\
\n\
reinitDictionary NO_REINIT_AT_ALL\n\
\n\
#\n\
# Parametr umozliwiajacy budowe slownika oscylacyjnego \n\
# Wszystkie gabory, dla ktorych w ramach jednej skali nie miesci sie okres oscylacji \n\
# zostana odrzucone ze slownika. \n\
#\n\
\n\
scaleToPeriodFactor 0.5 \n\
\n\
#\n\
# Dokladnosc obliczania iloczynu skalarnego dwoch gaborow.\n\
# Obliczanie iloczynu skalarnego dwoch gaborow, mozna rozbic na dwie operacje:\n\
# 1. mnozenie punkt po punkcie obydwu funkcji.\n\
# 2. sumowanie otrzymanych punktow.\n\
# Po wykonaniu kroku pierwszego dostajemy gabora (iloczyn 2 gaborow daje gabora).\n\
# Polozenie i skale tego gabora, mozna wyliczyc analitycznie. Dzieki temu mozna oszacowac\n\
# przedzial, na ktorym  krok drugi - sumowanie elementow do iloczynu skalarnego,\n\
# jest najbardziej efektywne numerycznie. Przedzial ten jest oczywiscie krotszy niz dlugosc skladki\n\
# wraz wraz z warunkami brzegowymi i dobrany tak, by dokladnosc wyliczenia iloczynu skalarnego\n\
# byla wieksza niz DOT_EPS. Zastosowana metoda oszacowania iloczynu skalarnego nie tylko pomaga\n\
# w skroceniu obliczen lecz takze umozliwia znajdowanie ortogonalnych gaborow, dla ktorych obliczanie\n\
# iloczynu skalarnego nie ma sensu. Jak zostalo to opisane wyzej, iloczyn dwoch gaborow daje gabora.\n\
# We wzorze tego gabora znajduje sie pewna stala, ktora moze byc zastosowana do szacowania ortogonalnosci\n\
# gaborow. Gabory, dla ktorych iloczyn daje te stala mniejsza niz DOT_EPS sa traktowane jako funkcje ortogonalne.\n\
# UWAGI:\n\
# a) w pierwszym kroku obliczen nie stosuje sie zadnej z powyzszcyh optymalizacji\n\
# b) w MP4 DOT_EPS domyslnie ustawione jest na wartosc rowna 1E-8.\n\
#\n\
\n\
DOT_EPS 1E-16\n\
\n\
#\n\
# Maks. ilosc gaborow, ktora moze byc uzyta do rozlozenia sygnalu\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
#\n\
\n\
maxNumberOfIterations 100\n\
\n\
#\n\
# Zatrzymaj, jesli dokompozycja opisala conajmniej x%  energii sygnalu\n\
# Wartosc argumentu: liczba zmiennoprzecinkowa, dodatnia.\n\
# Obowiazkowo w wartosci musi wystapic kropka dziesietna.\n\
#\n\
\n\
energyPercent 99.0\n\
\n\
#\n\
# Rodzaj algorytmu MP\n\
# Wartosc argumentu: liczba calkowita dodatnia.\n\
# Do wyboru sa nastepujace algorytmy: \n\
#  SMP  - SINGLE CHANNEL MATCHING PURSUIT\n\
#  MMP1 - MULTICHANNEL MATCHING PURSUIT, first algorithm\n\
#  MMP2 - MULTICHANNEL MATCHING PURSUIT, second algorithm\n\
#  MMP3 - MULTICHANNEL MATCHING PURSUIT, third algorithm (not existing now)\n\
#\n\
\n\
MP SMP \n\
\n\
#\n\
# Stala konwersji probka ->napiecie w uV\n\
# Obowiazkowo w wartosci musi wystapic kropka dziesietna.\n\
#\n\
\n\
convRate 100.0\n\
\n\
#\n\
# VERBOSE:\n\
# 1 - zapisuje do pliku tekstowego wygenerowany slownik\n\
# 2 - zapisuje do pliku tekstowego dopasowane gabory\n\
# 4 - wyswietlaj progressbar\n\
# Uwaga argumenty opcji verbose mozna ze soba laczyc, np. wpisujac 5 wypiszemy na ekran infromacje o kolejnych iteracjach\n\
# oraz zapiszemy do pliku tekstowego dopasowane gabory\n\
#\n\
\n\
VERBOSE 4\n\
\n\
"};

	if(exampleFile == NULL)
	{
	    fprintf(stderr," \n");
	    fprintf(stderr," ERROR: \n");
	    fprintf(stderr," CAN'T OPEN EXAMPLE FILE: example.set");
	    fprintf(stderr," \n");
	    return ERROR;
	}

        fprintf(exampleFile,"%s\n",exampleTest);
	fclose(exampleFile);
    
        printf(" Example file has been generated and written to file example.set \n");

	return SUCCESS;
    }
    else if((mode & TEST) || (mode & PROCESS))
    {
	/* alocate memory for components of configFile */

	setConfigFile(&configFile);

	if(strlen(argv[2])>LENGTH_OF_NAME_OF_CONFIG_FILE)
	{
	    char numberToString[5];

	    fprintf(stderr," FILE OPEN ERROR \n");
	    fprintf(stderr," TO LONG NAME OF CONFIG FILE: ");
	    fprintf(stderr," %s \n",argv[2]);
	    fprintf(stderr,"\nTHE MAXIMUM LENGTH IS: ");
	    fprintf(stderr," %s ",decimalToString(numberToString,LENGTH_OF_NAME_OF_CONFIG_FILE));
	    fprintf(stderr,"\n");
	    return ERROR;
	}

	strcpy(configFile.name,argv[2]);

	/* try to open config file */

	if(openConfigFile(&configFile,infoMessage)==ERROR)
	{
	    fprintf(stderr," %s ",infoMessage);
	    return ERROR;
	}

	/* read config file */
	if(readConfigFile(&configFile,infoMessage)==ERROR) goto ERROR_PROCEDURE;

	/* analyse line and command in config file
	 "primary" test of config file and parameters included in this file */
	if(findDataParametersInConfigFile(&configFile,&dataParameters,&gaborDictionary,&mp5Parameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;

	/* "extended" test of parameters read in config file */
	if(testDataParameters(&dataParameters,&gaborDictionary,&mp5Parameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;

	/* try to open data */
	if(dataParameters.dataFormat & FORMAT_ASCII)
	{
	    if(openAsciiDataFile(&dataParameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;
	}
	else if((dataParameters.dataFormat & FORMAT_FLOAT) || (dataParameters.dataFormat & FORMAT_SHORT))
	{
	    if(openBinaryDataFile(&dataParameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;
	}

	/* "primary" test of data file */
	if(dataParameters.dataFormat & FORMAT_ASCII)
	{
	    if(analyseAsciiDataFile(&dataParameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;
	}
	else if((dataParameters.dataFormat & FORMAT_FLOAT) || (dataParameters.dataFormat & FORMAT_SHORT))
	{
	    if(analyseBinaryDataFile(&dataParameters,infoMessage)==ERROR) goto ERROR_PROCEDURE;
	}

        /* set some constans values for MP5Parameters and DataParameters */
        setNumberOfAnalysedChannelsAndNumberOfResultsFiles(&mp5Parameters,&dataParameters);

        /* allocate memory for components of dataParameters */
        setDataParameters(&dataParameters);

        /* check dictionary tape and size */
        analyseDictionarySizeAndType(&mp5Parameters,&gaborDictionary);

        printSizeOfDictionaryAndSizeOfSinCosExpTables(&mp5Parameters,&gaborDictionary);

        /* create names of results files */
        createNamesOfResultFiles(&dataParameters,&mp5Parameters,&gaborDictionary);

	/* check wheter directories and results files are existing */
	if(testFilesAndDirectories(&dataParameters,&configFile,infoMessage) == ERROR) goto ERROR_PROCEDURE;

        printInfoAboutData(&dataParameters,&mp5Parameters,&gaborDictionary);

        if(mode & PROCESS)
        {
    	    /* now allocate memory for components of GaborDictionary */
    	    allocateDictionary(&mp5Parameters,&gaborDictionary);

    	    /* now allocate memory for components of MP5Parameters */
	    setMP5Parameters(&mp5Parameters, &gaborDictionary);

	    /* make sin/cos and exp tables */
	    makeSinCosExpTable(&mp5Parameters,&gaborDictionary);

	    /* results file */
	    if(openResultFiles(&dataParameters,infoMessage)==ERROR)
		goto ERROR_PROCEDURE;

	    unsigned short int offsetNumber  = 0;
	    unsigned short int channelNumber = 0;

	    /* results file is open, header is written into results file */
	    if(writeHeader(&dataParameters,&mp5Parameters,&gaborDictionary,infoMessage)==ERROR)
		goto ERROR_PROCEDURE;

	    printf("\n START MP DECOMPOSITION \n");

	    if(mp5Parameters.MPType & SMP)
	    {
		if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
		{
		    /* create dicionary */
		    makeDictionary(&mp5Parameters,&gaborDictionary);

		    /* test gabor's feature, for example find INCORRECT gabors */
		    testGaborFeature(&gaborDictionary);

		    if(dataParameters.verbose & VERBOSE_PRINT_DICTIONARY)
			printDictionaryToAsciFile(&dataParameters,&gaborDictionary);

		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			for(channelNumber=0;channelNumber<dataParameters.numberOfChosenChannels;channelNumber++)
			{
			    printf("\n --OFFSET--: %d, |CHANNEL|: %d\n\n",dataParameters.chosenOffsets[offsetNumber],dataParameters.chosenChannels[channelNumber]);

			    mp5Parameters.singleChannelSignalTable = *(dataParameters.processedDataMatrix + dataParameters.chosenChannels[channelNumber]-1);

			    firstIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);
			    nextIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);

			    if(writeSingleChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,channelNumber,infoMessage)==ERROR)
				goto ERROR_PROCEDURE;

			    resetDictionary(&gaborDictionary);
			}
		    }      
		}
		else if(mp5Parameters.reinitDictionary & REINIT_IN_CHANNEL_DOMAIN)
		{
		    for(channelNumber=0;channelNumber<dataParameters.numberOfChosenChannels;channelNumber++)
		    {	
		      /* create dicionary */
			reinitDictionary(&mp5Parameters,&gaborDictionary);

		      /* test gabor's feature, for example find INCORRECT gabors */
			testGaborFeature(&gaborDictionary);

			for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
			{
			    if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
				goto ERROR_PROCEDURE;

			    printf("\n |CHANNEL|: %d, --OFFSET--: %d\n\n",dataParameters.chosenChannels[channelNumber],dataParameters.chosenOffsets[offsetNumber]);

			    mp5Parameters.singleChannelSignalTable = *(dataParameters.processedDataMatrix + dataParameters.chosenChannels[channelNumber]-1);

			    firstIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);
			    nextIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);

			    if(writeSingleChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,channelNumber,infoMessage)==ERROR)
				goto ERROR_PROCEDURE;

			    resetDictionary(&gaborDictionary);
			}
		    }
		}
		else if(mp5Parameters.reinitDictionary & REINIT_IN_OFFSET_DOMAIN)
		{
		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			/* create dicionary */
			reinitDictionary(&mp5Parameters,&gaborDictionary);

			/* test gabor's feature, for example find INCORRECT gabors */
			testGaborFeature(&gaborDictionary);

			for(channelNumber=0;channelNumber<dataParameters.numberOfChosenChannels;channelNumber++)
			{
			    if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
				goto ERROR_PROCEDURE;

			    printf("\n |CHANNEL|: %d, --OFFSET--: %d\n\n",dataParameters.chosenChannels[channelNumber],dataParameters.chosenOffsets[offsetNumber]);

			    mp5Parameters.singleChannelSignalTable = *(dataParameters.processedDataMatrix + dataParameters.chosenChannels[channelNumber]-1);

			    firstIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);
			    nextIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);
 
			    if(writeSingleChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,channelNumber,infoMessage)==ERROR)
				goto ERROR_PROCEDURE;

			    resetDictionary(&gaborDictionary);
			}
		    }
		}
		else if(mp5Parameters.reinitDictionary & REINIT_AT_ALL)
		{
		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			for(channelNumber=0;channelNumber<dataParameters.numberOfChosenChannels;channelNumber++)
			{
			    /* create dicionary */
			    reinitDictionary(&mp5Parameters,&gaborDictionary);

			    /* test gabor's feature, for example find INCORRECT gabors */
			    testGaborFeature(&gaborDictionary);

			    if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
				goto ERROR_PROCEDURE;

			    printf("\n --OFFSET--: %d, |CHANNEL|: %d\n\n",dataParameters.chosenOffsets[offsetNumber],dataParameters.chosenChannels[channelNumber]);

			    mp5Parameters.singleChannelSignalTable = *(dataParameters.processedDataMatrix + dataParameters.chosenChannels[channelNumber]-1);

			    firstIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);
			    nextIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);

			    if(writeSingleChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,channelNumber,infoMessage)==ERROR)
				goto ERROR_PROCEDURE;
			}
		    }
		}
	    }
	    // writen by Artur Matysiak
	    else if(mp5Parameters.MPType & MMP1)
	    {
		if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
		{
		    /* create dicionary */
		    makeDictionary(&mp5Parameters,&gaborDictionary);

		    /* test gabor's feature, for example find INCORRECT gabors */
		    testGaborFeature(&gaborDictionary);
                
		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			printf("\n --OFFSET--: %d\n\n",dataParameters.chosenOffsets[offsetNumber]);

			mp5Parameters.multiChannelSignalTable = dataParameters.processedDataMatrix;

			firstIterationMMP1(&mp5Parameters,&dataParameters,&gaborDictionary);
			nextIterationMMP1(&mp5Parameters,&dataParameters,&gaborDictionary);

			if(writeMultiChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			resetDictionary(&gaborDictionary);
		    }
		}
		if(mp5Parameters.reinitDictionary &  REINIT_IN_OFFSET_DOMAIN)
		{

		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			/* create dicionary */
			reinitDictionary(&mp5Parameters,&gaborDictionary);

			/* test gabor's feature, for example find INCORRECT gabors */
			testGaborFeature(&gaborDictionary);

			if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			printf("\n --OFFSET--: %d\n\n",dataParameters.chosenOffsets[offsetNumber]);

			mp5Parameters.multiChannelSignalTable = dataParameters.processedDataMatrix;

			firstIterationMMP1(&mp5Parameters,&dataParameters,&gaborDictionary);
			nextIterationMMP1(&mp5Parameters,&dataParameters,&gaborDictionary);

			if(writeMultiChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			resetDictionary(&gaborDictionary);
		    }
		}
	    }    
	    //Artur Matysiak end
	    else if(mp5Parameters.MPType & MMP2)
	    {
		if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
		{
		    /* create dicionary */
		    makeDictionary(&mp5Parameters,&gaborDictionary);

		    /* test gabor's feature, for example find INCORRECT gabors */
		    testGaborFeature(&gaborDictionary);
                
		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			printf("\n --OFFSET--: %d\n\n",dataParameters.chosenOffsets[offsetNumber]);

			mp5Parameters.multiChannelSignalTable = dataParameters.processedDataMatrix;

			countMeanSignalOverChannels(&mp5Parameters,&dataParameters);
			mp5Parameters.singleChannelSignalTable = mp5Parameters.meanSignalTable;
				
			firstIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);
			nextIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);

			if(writeMultiChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			resetDictionary(&gaborDictionary);
		    }
		}
		if(mp5Parameters.reinitDictionary &  REINIT_IN_OFFSET_DOMAIN)
		{

		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			/* create dicionary */
			reinitDictionary(&mp5Parameters,&gaborDictionary);

			/* test gabor's feature, for example find INCORRECT gabors */
			testGaborFeature(&gaborDictionary);

			if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			printf("\n --OFFSET--: %d\n\n",dataParameters.chosenOffsets[offsetNumber]);

			mp5Parameters.multiChannelSignalTable = dataParameters.processedDataMatrix;

			countMeanSignalOverChannels(&mp5Parameters,&dataParameters);
			mp5Parameters.singleChannelSignalTable = mp5Parameters.meanSignalTable;

			firstIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);
			nextIterationSMP(&mp5Parameters,&dataParameters,&gaborDictionary);

			if(writeMultiChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			resetDictionary(&gaborDictionary);
		    }
		}
	    }
	    else if(mp5Parameters.MPType & MMP3)
	    {
		if(mp5Parameters.reinitDictionary & NO_REINIT_AT_ALL)
		{
		    /* create dicionary */
		    makeDictionary(&mp5Parameters,&gaborDictionary);

		    /* test gabor's feature, for example find INCORRECT gabors */
		    testGaborFeature(&gaborDictionary);
                
		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			printf("\n --OFFSET--: %d\n\n",dataParameters.chosenOffsets[offsetNumber]);

			mp5Parameters.multiChannelSignalTable = dataParameters.processedDataMatrix;

			firstIterationMMP3(&mp5Parameters,&dataParameters,&gaborDictionary);
			nextIterationMMP3(&mp5Parameters,&dataParameters,&gaborDictionary);

			if(writeMultiChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			resetDictionary(&gaborDictionary);
		    }
		}
		if(mp5Parameters.reinitDictionary &  REINIT_IN_OFFSET_DOMAIN)
		{

		    for(offsetNumber=0;offsetNumber<dataParameters.numberOfChosenOffsets;offsetNumber++)
		    {
			/* create dicionary */
			reinitDictionary(&mp5Parameters,&gaborDictionary);

			/* test gabor's feature, for example find INCORRECT gabors */
			testGaborFeature(&gaborDictionary);

			if(readDataFile(&dataParameters,dataParameters.chosenOffsets[offsetNumber],infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			printf("\n --OFFSET--: %d\n\n",dataParameters.chosenOffsets[offsetNumber]);

			mp5Parameters.multiChannelSignalTable = dataParameters.processedDataMatrix;

			firstIterationMMP3(&mp5Parameters,&dataParameters,&gaborDictionary);
			nextIterationMMP3(&mp5Parameters,&dataParameters,&gaborDictionary);

			if(writeMultiChannelResults(&dataParameters,&mp5Parameters,&gaborDictionary,offsetNumber,infoMessage)==ERROR)
			    goto ERROR_PROCEDURE;

			resetDictionary(&gaborDictionary);
		    }
		}
	    }
		} 
    }

    freeConfigFile(&configFile);
    closeFiles(&dataParameters);
    freeDataParameters(&dataParameters);
    freeMP5Parameters(&mp5Parameters,&gaborDictionary);
    freeDictionary(&gaborDictionary);
    return SUCCESS;

    ERROR_PROCEDURE:
	fprintf(stderr," \n");
        fprintf(stderr," ERROR: \n");
	fprintf(stderr,"%s",infoMessage);
        fprintf(stderr," \n");
	freeConfigFile(&configFile);
        closeFiles(&dataParameters);
        freeDataParameters(&dataParameters);
	freeMP5Parameters(&mp5Parameters,&gaborDictionary);
        freeDictionary(&gaborDictionary);

    return ERROR;

}
