%% Visualize spatial audio data with itaCoordinates
%
% <<../../pics/ita_toolbox_logo_wbg.png>>
% 
% This tutorial demonstrates how to work with an HRTF. At first, an HRTF is
% created, then certain directions are chosen and plotted. Finally, the HRTF is modified and stored. 
% 
% This file is part of the application HRTF_class for the ITA-Toolbox. All rights reserved.
% You can find the license for this m-file in the application folder.
% It is recommended to open this file as a Live Script.

%% Init object

% create sphere
coord   =  ita_sph_sampling_equiangular(11,36);
scatter(coord)

% calculate sphere's HRTF 
coord.r = 0.1;
pSphere = test_rbo_pressureSphere('sph',coord,'fftDeg',8);
HRTF_sphere    = ita_time_shift(itaHRTF(pSphere),pSphere.trackLength/2,'time');

%% Find Functions 
% Note regarding the angles: 
% As spherical coordinates are used, the angle theta = 0 degree is not in
% the horizontal  plane. Therefore it is defined for the interval [0, pi].

% define coordinates to look for two specific HRTFs
coordF          = itaCoordinates([1 pi/2 pi/2; 1 pi/2 pi/4],'sph');
% find the HRTF for the defined coordinates
HRTF_find       = HRTF_sphere.findnearestHRTF(coordF); 
% faster workaround for the first HRTF's calculation
HRTF_find1 = HRTF_sphere.findnearestHRTF(90,90);

% Slice of the HRTF function for the horizontal and the vertical plane
horizontalPlaneHRTF        = HRTF_sphere.sphericalSlice('theta_deg',90);
verticalPlaneHRTF          = HRTF_sphere.sphericalSlice('phi_deg',0);

%% Plot sliced HRTFs 
% plot frequency domain in dependence of the angle (elevation or azimuth)
figure
verticalPlaneHRTF.plot_freqSlice
pause(5)
horizontalPlaneHRTF.plot_freqSlice('earSide','R')

%% Plot ITD and HRTF in 
% ITD Plot for horinzontal plane
horizontalPlaneHRTF.plot_ITD('method','xcorr','plot_type','line')
title('ITDs for the horinzontal plane')

% plot the two HRTFs for the coordinates chosen above 
% plot in time domain
HRTF_find.pt 
title('HRTFs in time domain')

% plot in frequency domain
HRTF_find.pf
ylim([-20 10])
title('HRTFs in frequency domain')

% plot only right ear
HRTF_find.getEar('R').pf
ylim([-20 10])
title('Right ear`s HRTFs in frequency domain')
%% Play gui
pinkNoise = ita_generate('pinknoise',1,44100,12)*10;
HRTF_find.play_gui(pinkNoise); 

%% Binaural parameters
ITD = horizontalPlaneHRTF.ITD;  % different methods are available: see method in itaHRTF
%ILD = slicePhi.ILD;
 
%% Modifications
% calculate DTF
[DTF_sphere, common] = HRTF_sphere.DTF;

%plot HRTF and DTF for left and right ear for one position
HRTFvsDTF = ita_merge(DTF_sphere.findnearestHRTF(90,90),HRTF_sphere.findnearestHRTF(90,90));
HRTFvsDTF.pf
legend('DTF left','DTF right','HRTF left','HRTF right')
title('HRTF and DTF for position (90,90)')

% interpolate HRTF
phiI     = deg2rad(0:5:355);
thetaI   = deg2rad(15:15:90);
[THETA_I, PHI_I] = meshgrid(thetaI,phiI);
rI       = ones(numel(PHI_I),1);
coordI   = itaCoordinates([rI THETA_I(:) PHI_I(:)],'sph'); % itaCoordinates object

HRTF_interp = HRTF_sphere.interp(coordI);


%% Write and init
nameDaff_file = 'HRTF_sphere.daff';
HRTF_sphere.writeDAFFFile(nameDaff_file);

%HRTF_daff = itaHRTF('daff',nameDaff_file);

nameDaff_file2 = 'yourHRTF.daff';
if ~strcmp(nameDaff_file2,'yourHRTF.daff')
    HRTF_daff2 = itaHRTF('daff',nameDaff_file2);
    HRTF_daff2.plot_freqSlice
else
   ita_disp('use an existing daff-file') 
end