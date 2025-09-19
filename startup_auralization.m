function startup_auralization(bp)
% function startup_auralization(bp)
%
% This scripts initialises the auralization framework. It includes all
% necessary codes to the local MATLAB path and add any required third party
% software required for its use. All required dependencies are distributed 
% (in original form) and are located within the /third_party folder.
% 
% As recommended by the author to avoid any conflicts, the added paths will 
% be removed from the MATLAB directories when MATLAB is closed. This means 
% that startup needs to be run once after MATLAB has started.
%
% -------------------------------
% Author: Gil Felix Greco (ggrecow@gmail.com)
% Institution: Technische Universität Braunschweig 
%
% Date created: 20.11.2024
% Date last modified: 19.09.2025
% MATLAB version: 2024b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bp = [fileparts(mfilename('fullpath')) filesep]; % obtains the folder where this script is
end

addpath(bp);

main_dirs = { 'auralization', 'utilities' };

%% add main_dirs

for i_dir = 1:length(main_dirs)
    main_dir_i = [bp main_dirs{i_dir} filesep];
    if exist(main_dir_i,'dir')
        bAdd = 1;
        if bAdd
            addpath(main_dirs{i_dir});
            fprintf('%s.m: Main dir added to path:\n \t%s\n', mfilename, main_dir_i);
        end
    end
end

%% add ITA toolbox

ITA_folder        = [bp 'third_party' filesep 'ITA-Toolbox' filesep];

bAdd = ~exist('ita_toolbox_setup.m','file');
if bAdd
addpath(ITA_folder);
run ita_toolbox_setup;
end


%% add ART

ART        = [bp 'third_party' filesep 'ARTMatlab_v2023a' filesep];

bAdd = ~exist('AtmosphericPropagation.m','file');
if bAdd
    addpath(ART);
end

%% add Faddeeva codes for numerical error calculation
Faddeeva = [bp 'third_party' filesep 'Faddeeva_MATLAB' filesep];

bAdd = ~exist('Faddeeva_w.m','file');
if bAdd
    addpath(Faddeeva);
end

%% add AKtools
AKtools = [bp 'third_party' filesep 'AKtools' filesep];

bAdd = ~exist('AKphaseManipulation.m','file');
if bAdd
    addpath(AKtools);
end

%% add AKtools - plot tools 
AKtools_plotting = [bp 'third_party' filesep 'AKtools' filesep 'Plotting'];

bAdd = ~exist('AKp.m','file');
if bAdd
    addpath(AKtools_plotting);
end

%% add AKtools - plot tools brewer
AKtools_brewer = [bp 'third_party' filesep 'AKtools' filesep 'Plotting' filesep 'cbrewer'];

bAdd = ~exist('cbrewer.m','file');
if bAdd
    addpath(AKtools_brewer);
end

%% add data from FABIAN HRTF database 
HRTF= [bp 'third_party' filesep 'FABIAN_HRTF_DATABASE_v4' filesep];

bAdd = ~exist('FABIAN_HRIR_measured_HATO_0.m','file');
if bAdd
    addpath(HRTF);
end
