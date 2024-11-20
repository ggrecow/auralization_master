function startup_auralization(bp)
% function startup_auralization(bp)
%
% This scripts initialises the toolbox. As recommended by the authors, the
%   added paths will be removed from the MATLAB directories when MATLAB is
%   closed. This means that startup needs to be run once after MATLAB
%   has started.
%
% Author: Alejandro Osses
% modified by Gil F. Greco, 20.11.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bp = [fileparts(mfilename('fullpath')) filesep]; % obtains the folder where this script is
end

addpath(bp);

main_dirs = { 'auralization', ...
                      'utilities'};

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
