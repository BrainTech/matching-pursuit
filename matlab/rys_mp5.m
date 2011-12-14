function book = rys_mp5(nameOfBookFile,epochNumber,channelNumber,frequencyFactor)

%plot the t-f map, original signal and its mp reconstruction
figure
[book header epochSize] = readonebookv5rh(nameOfBookFile,epochNumber,channelNumber);

dimBase = epochSize;
maxF    = floor((dimBase/2)/frequencyFactor); %max freq on the t-f map
c_f     = header.pointsPerMicrovolt;          %stala kalibracji

[map,xx,yy] = mp2tfv5(book, header, epochSize, 1, 1, 0,maxF,0,dimBase); 
%mp2tfv5(book, header, epochSize, Dt, Df, minF, maxF, minT, maxT)

subplot('position',[.07 .3 .9 .65])
imagesc(xx,yy,map+1.0); set(gca,'ydir', 'normal')
ylabel('[Hz]');
xlabel('[s]');

TYPE      = 1;
MODULUS   = 2;
AMPLITUDE = 3;
POSITION  = 4;
SCALE     = 5;
FREQUENCY = 6;
PHASE     = 7;

samplingFrequency = header.samplingFrequency;

reconstruction = zeros(1,epochSize);
original       = book.signal;

for atom=1:size(book.atoms,1)
    amplitude = book.atoms(atom,AMPLITUDE);
    position  = book.atoms(atom,POSITION);
    width     = book.atoms(atom,SCALE)/samplingFrequency;
    frequency = book.atoms(atom,FREQUENCY)*(0.5*samplingFrequency);
    phase     = book.atoms(atom,PHASE);
    reconstruction = reconstruction + gabor(epochSize,samplingFrequency,amplitude,position,width,frequency,phase,book.atomsType(atom));
end

time = (0:1:epochSize-1)/samplingFrequency;
subplot('position',[.07 .06 .9 .17])

plot(time,original./c_f);
xlim([0 dimBase/samplingFrequency])
hold
plot(time,reconstruction./c_f,'red');
xlabel('s');
hold off

% figure
% %[px,py] = gradient(map,0.1,0.1);
% %quiver(px,py);
% 
% J = entropyfilt(map,ones(3,3));
% imshow(J,[]);