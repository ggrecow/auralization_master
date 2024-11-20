function [resPhase] = correctPhaseDelay(atmoDelay,dirPathLength,freqVector,c)
%CORRECTPHASEDELAY calculates the phase delay correction factors that can
%be applied to an urban TF

if nargin < 4
    c = 340;
end

deltaDelay = atmoDelay - dirPathLength/c;

resPhase = [0;exp(1j*2*pi*deltaDelay.*freqVector(2:end))];

end

