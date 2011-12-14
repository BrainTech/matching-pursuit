function [wigXY,xx,yy]=mp2tfv5(book, header, epochSize, Dt, Df, minF, maxF, minT, maxT)
% minF, maxF, minT, maxT, dt, df -- in samples

dimBase = maxT;
Fsamp   = header.samplingFrequency;

if nargin<5
    minF = 0;
    maxF = dimBase/2;
    minT = 0;
    maxT = dimBase;
end

if nargin<3
    DT = 1;
    DF = 1;
end

t = 1:maxT - minT; % skala czasu w punktach
f = 1:maxF - minF; % skala czestosci w punktach

Tsize = length(t);
Fsize = length(f);
Xsize = (maxT-minT)/Dt;
Ysize = (maxF-minF)/Df;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CAUTION: the values below are correct only if the MP decomposition
% was calculated with border conditions as zeros outside 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
book_min_t_points = 0;
book_max_t_points = dimBase;
book_min_f_points = 0;
book_max_f_points = dimBase/2;

[wigXY, xx, yy] = deal([],[],[]);

% sprawdzamy czy wyjsciowe piksle zawieraja calkowita liczbe pikseli dokladnej mapy

if mod(Tsize,Xsize)~=0

    disp(sprintf('wrong Xsize, possible values are: \n'));
    
    for i=1:Tsize
        if mod(Tsize,i)==0
            disp(i)
        end
    end
    error('... so the end');
elseif mod(Fsize,Ysize)~=0
    disp(sprintf('wrong Ysize, possible values are: \n'));
    for i=1:Fsize
        if mod(Fsize,i)==0
            disp(i)
        end
    end
    error('... so the end');
end

DX = Tsize/Xsize;
DY = Fsize/Ysize;
tx = minT + t(1:DX:Tsize); %brzegi duzych pikseli na osi czasu
fy = minF + f(1:DY:Fsize); %brzegi duzych pikseli na osi czestosci

tlowx = ((tx-1)); %granice calkowania - duze piksele
tuppx = ( tx-1 + DX); %granice calkowania - duze piksele

flowy = ((fy-1));  %granice calkowania - duze piksele
fuppy = ( fy-1 + DY);  %granice calkowania - duze piksele

xx = tx+DX/2; % srodki duzych pikseli na osi czasu (T==>X)
yy = fy+DY/2; % srodki duzych pikseli na osi czestosci (F==>Y)

df=dimBase/Fsamp;% konwersja Hz na punkty

PI4=4.0*pi;
SQRT_PI = sqrt(pi);
BI_SQRT_PI = 2*SQRT_PI;
num_atom=size(book.atoms, 1);
DX_DY = DX*DY;

wigXY=zeros(Xsize,Ysize);
wig_gabXY=sparse([]);

hm = floor(num_atom/10);

for k=1:num_atom
   
    modulus   = book.atoms(k,2)./header.pointsPerMicrovolt;
    amplitude = book.atoms(k,3)./header.pointsPerMicrovolt;
    scale     = book.atoms(k,5);
    position  = book.atoms(k,4);
    frequency = ((book.atoms(k,6))/pi)*(epochSize/2.0);
 
    if book.atomsType(k)=='H'% sinus
        if(frequency < Fsize)
            floor((frequency - minF)/DY) + 1;          
            freqXY = fy(floor((frequency - minF)/DY) + 1);
        
            if (freqXY>=1) & (freqXY<=Ysize)
                wigXY(:,freqXY)=wigXY(:,freqXY)+(modulus^2)/Tsize; % normalizacja HAK
            end
        end
    elseif book.atomsType(k)=='G' | book.atomsType(k)=='N' %gabor lub gauss
        %liczenie w osiach XxY
        g1=(BI_SQRT_PI/scale)*(tlowx - position);
        g2=(BI_SQRT_PI/scale)*(tuppx - position);
        hgab_tx=(erf(g2)-erf(g1))'; %normowanie ponizej
        gf1=(SQRT_PI*scale/dimBase)*(flowy - frequency);
        gf2=(SQRT_PI*scale/dimBase)*(fuppy - frequency);
        hgab_fy=(erf(gf2)-erf(gf1));%normowanie ponizej     
        wig_gabXY=kron(sparse(hgab_tx), sparse(hgab_fy));     
        g1_book=(BI_SQRT_PI/scale)*(book_min_t_points - position);
        g2_book=(BI_SQRT_PI/scale)*(book_max_t_points - position);
        gf1_book=(SQRT_PI*scale/dimBase)*(book_min_f_points - frequency);
        gf2_book=(SQRT_PI*scale/dimBase)*(book_max_f_points - frequency);

        %calka z czesci atomu ktora jest w granicach liczenia ksiazki MP
        NORM_XY = (erf(g2_book)-erf(g1_book))*(erf(gf2_book)-erf(gf1_book));
         
        if NORM_XY~=0
            wigXY=wigXY+modulus.^2*full(wig_gabXY)./NORM_XY;
        end
    elseif book.atomsType(k)=='D' % dirac
          transXY=tx(1+floor((position-minT)/DX));% przesuwamy czestosc - floor -zeby pasowala do siatki wig, 
          if transXY>=1 & transXY<=Xsize
             wigXY(transXY,:)=wigXY(transXY,:)+(modulus^2)/Fsize; %normalizacja HAK
          end
    end

    if mod(k,hm) == 0
        fprintf(1,'*');
    end
end;

%PJD
wigXY=sqrt(wigXY); %U
wigXY=wigXY'; 

fprintf(1,'\n');
xx=xx/Fsamp;
yy=yy*Fsamp/dimBase;