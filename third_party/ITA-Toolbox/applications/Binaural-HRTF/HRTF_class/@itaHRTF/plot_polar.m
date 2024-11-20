function plot_polar(this, varargin)
% itaHRTF.plot_polar ... by default plots the HRTF magnitude of a given ear 
% in the horizontal plane 
%
% if you want to plot a different plane use hrtf.sphericalSlice() to get
% the coorect plane before plotting
%
%   options (default):
%       'freq'  ([125, 250, 500, 1000, 2000, 4000, 8000, 16000]) in Hz
%       'ear'   ('L'): which ear to plot


% Hark Braren -- hbr@akustik.rwth-aachen.de

defaultFreqs = [125, 250, 500, 1000, 2000, 4000, 8000, 16000];

sArgs.freq  = defaultFreqs;
sArgs.ear   = 'L';

sArgs = ita_parse_arguments(sArgs,varargin);

if this.nPointsAzimuth == 1 || this.nPointsAzimuth == 2
    %special case: saggital plane -> multiple theta angles at one phi angle
    %workaroun: rotate into xz plane and use existing plot function
    this.channelCoordinates.phi = this.channelCoordinates.phi-min(this.channelCoordinates.phi);
    plotFunction = @(freq) ita_plot_polar(this.getEar(sArgs.ear),freq,'newFigure',false,'plotPlane','xz');
else
    if this.nPointsElevation == 1
        %move data into horizontal plane
        this.channelCoordinates.z = 0;
    end
   %default case horziontal plane
   plotFunction = @(freq) ita_plot_polar(this.getEar(sArgs.ear),freq,'newFigure',false,'plotPlane','xy');
end


% figure()
for iFreq = 1:numel(sArgs.freq)
    currFreq = sArgs.freq (iFreq);
    plotFunction(currFreq);
    hold on
end

title('HRTF Magnitude in [dB]')
set(gca,'ThetaZeroLocation','top')

legend(num2str(sArgs.freq'))