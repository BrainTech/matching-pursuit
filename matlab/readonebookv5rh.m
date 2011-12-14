function [book header epochSize] = readonebookv5rh(nameOfResultsFile,epochNumber,channelNumber)

sizeOfFloat = 4;

% atom = struct('type',0,...
%               'parameters',0);
atom = zeros(0,7); % atom(n,1) - atoms #n type; atom(n,2:7) - atoms parameters

book = struct('signal',[],...
    'atoms',atom,...
    'atomsType',char([]));

header = struct('samplingFrequency',[],...
    'pointsPerMicrovolt',[],...
    'numberOfChannelsInDataFile',[],...
    'energyPercent',[],...
    'maximalNumberOfIterations',[],...
    'sizeOfDictionary',[],...
    'typeOfDictionary',char([]),...
    'dateField',char([]));


FIELD_DESCRIPTOR_SIGNATURE   = 'cc';
SEGMENT_DESCRIPTOR_SIGNATURE = 'ci';

COMMENT_SEGMENT_IDENTITY      = 1;

FILE_HEADER_SEGMENT_IDENTITY  = 2;
WEB_SITE_LINK_FIELD_IDENTITY  = 3;
DATE_FIELD_IDENTITY 		  = 4;
SIGNAL_FIELD_IDENTITY         = 5;
DECOMPOSING_FIELD_IDENTITY    = 6;

OFFSET_SEGMENT_IDENTITY   	  = 7;
SIGNAL_SEGMENT_IDENTITY   	  = 8;
ATOMS_SEGMENT_IDENTITY   	  = 9;

OFFSET_SEGMENT_HEADER_SIGNATURE = 'hi';
SIGNAL_SEGMENT_HEADER_SIGNATURE = 'h';
ATOMS_SEGMENT_HEADER_SIGNATURE  = 'h';

resultsFile = fopen(nameOfResultsFile,'rb','ieee-be');

fseek(resultsFile,0,'eof');
sizeOfResultsFile = ftell(resultsFile);
fseek(resultsFile,0,'bof');

readMagicAndComments(resultsFile);
[webSiteLinkField,dateField,signalField,decomposingField] = readFileHeader(resultsFile);

% fprintf(' \n');
% fprintf(' web site link:                %s\n',webSiteLinkField.webSiteLink);
% fprintf(' \n');
% fprintf(' SIGNAL INFORMATION:\n');
% fprintf(' sampling frequency:           %f\n',signalField.samplingFrequency);
% fprintf(' pointsPerMicrovolt:           %f\n',signalField.pointsPerMicrovolt);
% fprintf(' numberOfChannelsInDataFile:   %u\n',signalField.numberOfChannelsInDataFile);
% fprintf(' \n');
% fprintf(' DECOMPOSING INFORMATION:\n');
% fprintf(' energy percent of the signal: %f \n',decomposingField.energyPercent);
% fprintf(' maximal number of iterations: %u \n',decomposingField.maximalNumberOfIterations);
% fprintf(' size of the dictionary:       %u \n',decomposingField.sizeOfDictionary);
% if(decomposingField.typeOfDictionary=='F')
%     fprintf(' type of the dictionary:       fixed\n');
% elseif(decomposingField.typeOfDictionary=='S')
%     fprintf(' type of the dictionary:       stochastic\n');
% end
% fprintf(' \n');
% fprintf(' date:                         %s\n',dateField.date);
% fprintf(' \n');

header.samplingFrequency          = signalField.samplingFrequency;
header.pointsPerMicrovolt         = signalField.pointsPerMicrovolt;
header.numberOfChannelsInDataFile = signalField.numberOfChannelsInDataFile;
header.energyPercent              = decomposingField.energyPercent;
header.maximalNumberOfIterations  = decomposingField.maximalNumberOfIterations;
header.sizeOfDictionary           = decomposingField.sizeOfDictionary;
header.typeOfDictionary           = decomposingField.typeOfDictionary;
header.dateField                  = dateField.date;

tmpEpochNumber        = 0;
tmpChannelNumber      = 0;
positionInResultsFile = ftell(resultsFile);
positionInSegment     = 0;
atomCounter           = 0;

while(positionInResultsFile<sizeOfResultsFile)

    codeOfPrimarySegment     = fread(resultsFile,1,'uchar');
    sizeOfPrimaryDataSegment = fread(resultsFile,1,'uint32');

    if(codeOfPrimarySegment==COMMENT_SEGMENT_IDENTITY)
        fseek(resultsFile,sizeOfPrimaryDataSegment,'cof');
    end

    if(codeOfPrimarySegment==OFFSET_SEGMENT_IDENTITY)

        tmpEpochNumber = fread(resultsFile,1,'uint16');
        epochSize      = fread(resultsFile,1,'uint32');

        if(tmpEpochNumber == epochNumber)

            positionInPrimarySegment = getSizeOf(OFFSET_SEGMENT_HEADER_SIGNATURE);

            while(positionInPrimarySegment<sizeOfPrimaryDataSegment)

                codeOfSecondarySegment     = fread(resultsFile,1,'uchar');
                sizeOfSecondaryDataSegment = fread(resultsFile,1,'uint32');

                if(codeOfSecondarySegment==SIGNAL_SEGMENT_IDENTITY)
                    tmpChannelNumber = fread(resultsFile,1,'uint16');

                    if(tmpChannelNumber == channelNumber)
                        book.signal = fread(resultsFile,epochSize,'float32');
                    else
                        fseek(resultsFile,sizeOfSecondaryDataSegment - getSizeOf(SIGNAL_SEGMENT_HEADER_SIGNATURE),'cof');
                    end
                elseif(codeOfSecondarySegment==ATOMS_SEGMENT_IDENTITY)
                    tmpChannelNumber = fread(resultsFile,1,'uint16');
                    if(tmpChannelNumber == channelNumber)
                        atomsCounter = 1;
                        positionInSecondarySegment = getSizeOf(ATOMS_SEGMENT_HEADER_SIGNATURE);
                        while(positionInSecondarySegment<sizeOfSecondaryDataSegment)
                            [book.atoms(atomsCounter,:), sizeOfField] = readOneAtom(resultsFile);
                            positionInSecondarySegment = positionInSecondarySegment + sizeOfField;
                            atomsCounter = atomsCounter + 1;
                        end
                    else
                        fseek(resultsFile,sizeOfSecondaryDataSegment - getSizeOf(ATOMS_SEGMENT_HEADER_SIGNATURE),'cof');
                    end
                else
                    fprintf(' unknow SEGMENT [size: %hu, code: %hu] !\n',sizeOfSecondaryDataSegment, codeOfSecondarySegment);
                    fseek(resultsFile,sizeOfSecondaryDataSegment,'cof');
                end

                positionInPrimarySegment = positionInPrimarySegment + sizeOfSecondaryDataSegment + getSizeOf(SEGMENT_DESCRIPTOR_SIGNATURE);

            end
        else
            fseek(resultsFile,sizeOfPrimaryDataSegment - getSizeOf(OFFSET_SEGMENT_HEADER_SIGNATURE),'cof');
        end
    else
        fprintf(' unknow SEGMENT [size: %hu, code: %hu] !\n',sizeOfPrimaryDataSegment, codeOfPrimarySegment);
        fseek(resultsFile,sizeOfPrimaryDataSegment,'cof')
    end

    positionInResultsFile = ftell(resultsFile);

end

%ConversionCoeff=[3.5 -123.5 1445 -5532]; %conversion coefficients HAK
%(without shift)
ConversionCoeff = [3.5 13 8.5 71]; %conversion coefficients HAK
%(with shift 13)
book.atomsType = char(polyval(ConversionCoeff,book.atoms(:,1)-13));
fclose(resultsFile);

%[sortedModuluses, indexes] = sort(book.atoms(:,2),1,'descend');
%book.atoms = book.atoms(indexes,:);

function readMagicAndComments(resultsFile)

COMMENT_SEGMENT_IDENTITY      = 1;

magic = fread(resultsFile,6,'char');

if((strcmp(char(magic'),'MPv5.0')~=1) & (strcmp(char(magic'),'MPv6.0')~=1))
    fclose(resultsFile);
    error(' bad file format');
end

codeOfPrimarySegment     = fread(resultsFile,1,'uchar');
sizeOfPrimaryDataSegment = fread(resultsFile,1,'uint32');

if(codeOfPrimarySegment==COMMENT_SEGMENT_IDENTITY)
    fseek(resultsFile,sizeOfPrimaryDataSegment,'cof');
end

function [webSiteLinkField,dateField,signalField,decomposingField] = readFileHeader(resultsFile)

COMMENT_SEGMENT_IDENTITY 	  = 1;
FILE_HEADER_SEGMENT_IDENTITY  = 2;
WEB_SITE_LINK_FIELD_IDENTITY  = 3;
DATE_FIELD_IDENTITY 		  = 4;
SIGNAL_FIELD_IDENTITY         = 5;
DECOMPOSING_FIELD_IDENTITY    = 6;

OFFSET_SEGMENT_IDENTITY       = 7;
CHANNEL_SEGMENT_IDENTITY   	  = 8;

webSiteLinkField = struct('webSiteLink',char([]));
dateField        = struct('date',char([]));

signalField = struct('samplingFrequency',0,...
    'pointsPerMicrovolt',0.0,...
    'numberOfChannelsInDataFile',0);

decomposingField = struct('energyPercent',0.0,...
    'maximalNumberOfIterations',0,...
    'sizeOfDictionary',0,...
    'typeOfDictionary',char([]));

codeOfSegment     = fread(resultsFile,1,'uchar');
sizeOfDataSegment = fread(resultsFile,1,'uint32');

if(codeOfSegment==FILE_HEADER_SEGMENT_IDENTITY)

    positionInFileHeader = 0;

    while(positionInFileHeader<sizeOfDataSegment)

        codeOfField     = fread(resultsFile,1,'uchar');
        sizeOfFieldData = fread(resultsFile,1,'uchar');

        if(codeOfField == WEB_SITE_LINK_FIELD_IDENTITY)

            webSiteLinkField.webSiteLink = char(fread(resultsFile,sizeOfFieldData,'char'));

        elseif(codeOfField == DATE_FIELD_IDENTITY)

            dateField.date = char(fread(resultsFile,sizeOfFieldData,'char'));

        elseif(codeOfField == SIGNAL_FIELD_IDENTITY)
            
            signalField.samplingFrequency          = fread(resultsFile,1,'float32');
            signalField.pointsPerMicrovolt         = fread(resultsFile,1,'float32');
            signalField.numberOfChannelsInDataFile = fread(resultsFile,1,'uint16');

        elseif(codeOfField == DECOMPOSING_FIELD_IDENTITY)

            decomposingField.energyPercent             = fread(resultsFile,1,'float32');
            decomposingField.maximalNumberOfIterations = fread(resultsFile,1,'uint32');
            decomposingField.sizeOfDictionary          = fread(resultsFile,1,'uint32');
            decomposingField.typeOfDictionary          = char(fread(resultsFile,1,'char'));

        else

            fread(resultsFile,sizeOfFieldData,'uchar');
            fprintf(' unknow block [size: %hu code: %hu] !\n',sizeOfFieldData, codeOfField);

        end

        positionInFileHeader = positionInFileHeader + sizeOfFieldData + 2;

    end
end

function [atom, sizeOfAtomsField] = readOneAtom(resultsFile)

DIRACDELTA_IDENTITY		      = 10;
GAUSSFUNCTION_IDENTITY 	 	  = 11;
SINCOSWAVE_IDENTITY    	 	  = 12;
GABORWAVE_IDENTITY     	 	  = 13;

% atom(1,1) -> modulus
% atom(1,2) -> amplitude
% atom(1,3) -> position
% atom(1,4) -> scale
% atom(1,5) -> frequency
% atom(1,6) -> phase

% atom = struct('type',0,...
%               'parameters',0);
atom = zeros(1,7); % atom(n,1) - atoms #n type; atom(n,2:7) - atoms parameters

type       = fread(resultsFile,1,'char');
sizeOfAtom = fread(resultsFile,1,'char');

sizeOfAtomsField = 2 + sizeOfAtom;

if(type==DIRACDELTA_IDENTITY) % DIRAC DELTA
%     atom.type = 'D';
%     atom.parameters(1,1) = fread(resultsFile,1,'float32');
%     atom.parameters(1,2) = fread(resultsFile,1,'float32');
%     atom.parameters(1,3) = fread(resultsFile,1,'float32');
    [x num] = fread(resultsFile,[1,3],'float32');
    %if num ... %test integrity
    atom(2:4) = x;
elseif(type==GAUSSFUNCTION_IDENTITY) % GAUSS FUNCTION
%     atom.type = 'N';
%     atom.parameters(1,1) = fread(resultsFile,1,'float32');
%     atom.parameters(1,2) = fread(resultsFile,1,'float32');
%     atom.parameters(1,3) = fread(resultsFile,1,'float32');
%     atom.parameters(1,4) = fread(resultsFile,1,'float32');
    [x num] = fread(resultsFile,[1,4],'float32');
    %if num ... %test integrity
    atom(2:5) = x;
elseif(type==SINCOSWAVE_IDENTITY) % HARMONIC (SIN/COS) WAVE
%     atom.type = 'H';
%     atom.parameters(1,1) = fread(resultsFile,1,'float32');
%     atom.parameters(1,2) = fread(resultsFile,1,'float32');
%     atom.parameters(1,5) = fread(resultsFile,1,'float32');
%     atom.parameters(1,6) = fread(resultsFile,1,'float32');
    [x num] = fread(resultsFile,[1,4],'float32');
    %if num ... %test integrity
    atom(2:3) = x(1:2);
    atom(6:7) = x(3:4);
elseif(type==GABORWAVE_IDENTITY) % GABOR WAVE
%     atom.type = 'G';
%     atom.parameters(1,1) = fread(resultsFile,1,'float32');
%     atom.parameters(1,2) = fread(resultsFile,1,'float32');
%     atom.parameters(1,3) = fread(resultsFile,1,'float32');
%     atom.parameters(1,4) = fread(resultsFile,1,'float32');
%     atom.parameters(1,5) = fread(resultsFile,1,'float32');
%     atom.parameters(1,6) = fread(resultsFile,1,'float32');
    [x num] = fread(resultsFile,[1,6],'float32');
    %if num ... %test integrity
    atom(2:7) = x;
else
    %atom.type = 'U';
    error(' unknow ATOM [size: %hu, code: %hu] !\n',sizeOfAtom,type);
    %fseek(resultsFile,sizeOfAtom,'cof');
end

atom(1,1) = type;

function size = getSizeOf(signature)

% definitions of chars in signature
%	'c' - char/unsigned char
%	'h' - short int/unsigned short int
%	'i' - int/unsigned int
%	'f' - float

sizeOfChar  = 1;
sizeOfShort = 2;
sizeOfInt   = 4;
sizeOfFloat = 5;

size = 0;

for(counter=1:length(signature))
    switch(signature(counter))
        case 'c'
            size = size + sizeOfChar;
        case 'h'
            size = size + sizeOfShort;
        case 'i'
            size = size + sizeOfInt;
        case 'f'
            size = size + sizeOfFloat;
        otherwise
            break;
    end
end
