% samplingFrequency = 256;
% 
% sizeOfSignal      = 1024;
% samplingFrequency = 256;
%                                        % amp   pos[sec.]  width[sec.]   freq [Hz]    phase
% g = gabor(sizeOfSignal,samplingFrequency,4,      1.5,        0.5,          5.0,      0.1,'G'); % GABOR FUNCTION
% g2= gabor(sizeOfSignal,samplingFrequency,10,     3.5,        1.0,          12.0,     0.7,'G'); % GABOR FUNCTION
% g3= gabor(sizeOfSignal,samplingFrequency,20,     3.5,        0.5,          2.0,     -0.7,'G'); % GABOR FUNCTION
% 
% n = gabor(sizeOfSignal,samplingFrequency,4,      2.5,        0.1,          0.0,      0.0,'N'); % GAUSS FUNCTION
% d = gabor(sizeOfSignal,samplingFrequency,15.0,   3.0,        0.0,          0.0,      0.0,'D'); % DIRAC'S DELTA
% h = gabor(sizeOfSignal,samplingFrequency,1.5,    0.0,        0.0,          20,       0.0,'H'); % HARMONIC WAVE
% 
% signal = g+g2+g3 + n + d + h;
% plot(signal)
% % signal = g;
% % 
% file = fopen('mp5_demo.dat','wb');
% fwrite(file,signal,'float32');
% fclose(file);

[book header offsetDimension]=readonebookv5rh('S1_testing_fs512_30_0_45Hz_r1_smp.b', 1,1);

% dimBase=offsetDimension;
% maxF=dimBase/4; %max freq on the t-f map
% c_f=1; %stala kalibracji

% figure
% [map,xx,yy]=mp2tfv5(book, header, offsetDimension, 1, 1, 0,maxF,0,dimBase); 
%mp2tfv5(book, header, offsetDimension, Dt, Df, minF, maxF, minT, maxT)

% subplot('position',[.07 .3 .9 .65])
% imagesc(xx,yy,map); set(gca,'ydir', 'normal')

TYPE      = 1;
MODULUS   = 2;
AMPLITUDE = 3;
POSITION  = 4;
SCALE     = 5;
FREQUENCY = 6;
PHASE     = 7;

samplingFrequency = header.samplingFrequency;

reconstruction = zeros(1,offsetDimension);
original       = book.signal;

for atom=1:size(book.atoms,1)
    amplitude = book.atoms(atom,AMPLITUDE);
    position  = book.atoms(atom,POSITION);
    width     = book.atoms(atom,SCALE)/samplingFrequency;
    frequency = book.atoms(atom,FREQUENCY)*(0.5*samplingFrequency);
    phase     = book.atoms(atom,PHASE);
    reconstruction = reconstruction + gabor(offsetDimension,samplingFrequency,amplitude,position,width,frequency,phase,book.atomsType(atom));
end
time = (0:1:offsetDimension-1)/samplingFrequency;
% subplot('position',[.07 .06 .9 .17])
plot(time,original./header.pointsPerMicrovolt);
xlim([0 offsetDimension/samplingFrequency])
hold
plot(time,reconstruction./header.pointsPerMicrovolt,'red');
hold off

