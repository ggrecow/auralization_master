function hrtfOut = ita_individualizeITDByMeasurements(hrtfIn,headWidth,headHeight,headFrontCircumference)
%ITA_INDIVIDUALIZEITDBYMEASUREMENTS - +++ Short Description here +++
%  This function uses the ITA-toolbox to measure and generates an
%  HRTF with an individualisation of the head, which can be saved.
%  The following parameters are ued to specify the head size: 
%  width, height and frontal circumference.
%
%  Syntax:
%   hrtfObjOut = ita_individualizeITDByMeasurements(hrtfObjIn, options)
%
%   Options (default):
%           'hrtfObjIn'         : original HRTF of which the ITD is to be adapted
%           'headWitdh'  (148)   : head width of subject, measured fromt
%                                  tragus to tragus with caliper [mm] 
%           'headHeight' (180)   : head height, measured vertically from
%                                  top of head to ear canal entrance [mm] 
%           'headFrontCircumference'  (270)  : head circumference, measured from tragus to tragus over bend of nose  shortest distance) [mm] 
%
%   How to measure: see HRTF_ITD_individualization_measurements.png 
%       image(imread('HRTF_ITD_individualization_measurements.png'));
%   
%  Example:
%   hrtfObjOut = ita_individualizeITDByMeasurements(hrtfObjIn)
%
%  See also:
%   itaHRTF, itaAnthroHRTF
%
%   Reference page in Help browser 
%        <a href="matlab:doc ita_individualizeITDByMeasurements">doc ita_individualizeITDByMeasurements</a>
%
% Related Publication(s)
%    Bomhardt, R., & Fels, J. (2014, October). Analytical Interaural Time Difference Model for the Individualization of Arbitrary
%    Head-Related Impulse Responses. Audio Engineering Society Convention 137. http://www.aes.org/e-lib/browse.cfm?elib=17454
%

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>

% Authors: 2021 Hark Braren -- Email: hark.braren@akustik.rwth-aachen.de


%% Input Check
if nargin < 1
    ita_verbose_info('Please supply a HRTF to individualize!!',0)
    return
end

%% When parameters are missing: use default parameters from ITA Artificial Head
if nargin < 4
    ita_verbose_info('Incomplete set of measurements, using default values for the rest',1)
end

% default parameters
if nargin < 2
    ita_verbose_info('Using headWidth = 148 [mm]',2)
    headWidth = 148; %mm
end
if nargin < 3
    ita_verbose_info('Using headHeight = 180 [mm]',2)
    headHeight = 180; %mm
end
if nargin < 4
    ita_verbose_info('Using headFrontCircumference = 270 [mm]',2)
    headFrontCircumference = 270; %mm
end

%% calculate measurements in accordance with Bomhardt, Fels (2014) AES
ind_w = headWidth/2;                                            % in mm
ind_h = headHeight;                                             % in mm
ind_d = ita_HRTFindividualization_circ2radius(ind_w, headFrontCircumference);         % in mm


%% calculate phase/ITD on elipsoid
indHRTF_subj = itaAnthroHRTF(hrtfIn, 'w', ind_w/1000 ,'d', ind_d/1000, 'h', ind_h/1000); %function wants unis in [m]
indHRTF_subj.calcEllipsoid = true;


%% use elipsoid ITD phase component on original magnitude data - HBR

%remove ITD information
hrtfOut = ita_minimumphase(hrtfIn); 
meanDelay = mean(hrtfIn.meanTimeDelay);

% add ITD from elipsoid model
phaseITD = angle(indHRTF_subj.freqData)-angle(ita_minimumphase(indHRTF_subj).freqData);
hrtfOut.freqData = hrtfOut.freqData.*exp(-1j*phaseITD);

% add mean time delay back in
hrtfOut = ita_time_shift(hrtfOut,meanDelay,'time');

%% Add history line
hrtfOut = ita_metainfo_add_historyline(hrtfOut,mfilename,varargin);

%end function
end