function [signal, time] = gabor(sizeOfSignal,sampleFrequency,atomAmplitude,atomPosition,atomWidth,atomFrequency,atomPhase,atomType)

% the input parameters of the function should be specified in SI units:
% sampleFrequency in Hz
% width in sec.
% atomPosition in sec.
% atomAmplitude in uV, V etc (depends on the signal)
% atomPhase in radians (?)

time = 0:1:sizeOfSignal-1;
%parameters in samples (points):
position  = atomPosition;
width     = atomWidth*sampleFrequency;
frequency = (atomFrequency/(0.5*sampleFrequency))*pi;

if atomType=='H'
    signal = atomAmplitude*cos(frequency*time + atomPhase);
elseif atomType=='D'
    signal = zeros(size(time));
    signal(position + 1) = atomAmplitude;
elseif atomType=='N'
    signal = atomAmplitude*exp(-pi.*((time-position)/width).^2);
elseif atomType=='G'
    signal = atomAmplitude*exp(-pi.*((time-position)/width).^2).*cos(frequency.*(time-position) + atomPhase);
end

