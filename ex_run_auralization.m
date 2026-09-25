% script ex_run_auralization.m
%
% Gives an example on how to call the auralization framework
%
% INPUT DATA:
%
%   - exemplary case available under the 'input_data' folder
%
%   - Generic flyover of an A319 during approach procedure
%
%   - has only one receiver, positioned at (x,y,z) = (0, 0, 1.2) meters
%
% FOLDER STRUCTURE:
%
%   input_data/
%     auralization_input.dat           <- shared by all examples
%     geschw_hoehe_verlauf.dat         <- shared by all examples
%     example_homogeneous_atmosphere/
%       input_file_auralization.ini
%     example_sounding_atmosphere/
%       input_file_auralization_import_atmosphere.ini
%       atmosphere_soundings/
%
%   auralization_input.dat - sound emissions per time step, per sound
%   source. File needs to have this name.
%
%   geschw_hoehe_verlauf.dat - flight trajectory and associated operational
%   conditions. File needs to have this name.
%
%   example_homogeneous_atmosphere/input_file_auralization.ini -
%   initializing file containing inputs related to signal processing and
%   setup of the atmospheric conditions used for the sound propagation
%   simulation. Here, a homogeneous atmosphere is assumed.
%
%   example_sounding_atmosphere/input_file_auralization_import_atmosphere.ini -
%   same as above, but the atmosphere is imported from atmospheric soundings
%   (.txt files stored in the 'atmosphere_soundings' subfolder), as provided by
%   http://weather.uwyo.edu/upperair/sounding_legacy.html. The sounding
%   (i.e. real measurement of atmosphere parameters over height, or a
%   measured inhomogeneous atmosphere) to be used for sound propagation
%   simulation is defined inside the .ini file.
%
% POSSIBLE CALLS:
%
%   auralization_master(core_path)
%   auralization_master(core_path, tag)
%   auralization_master(core_path, tag, input_file)
%   auralization_master(core_path, tag, input_file, results_path)
%
%   core_path    - folder containing the .dat input files (required)
%   tag          - name of the results folder; also used in all output
%                  file names (default: 'auralization_results')
%   input_file   - full path to the .ini file (default: the .ini file
%                  inside <core_path>; an error is thrown if <core_path>
%                  does not contain exactly ONE .ini file)
%   results_path - folder where results are stored; created if non
%                  existent (default: <core_path>)
%
%   NOTE: with the example folder structure above, <core_path> contains
%   no .ini file, so <input_file> must be provided. The shorter calls are
%   meant for case folders holding both the .dat files and a single .ini.
%
%   NOTE: all paths are built from the location of this script, so it can
%   be run regardless of MATLAB's current folder. Run the whole script
%   (F5 / Run), not single sections, otherwise the script location is unknown.
%
% Assumption: all required inputs are available
%
% -------------------------------
% Author: Gil Felix Greco (ggrecow@gmail.com)
% Institution: Technische Universität Braunschweig
%
% Date created: 13.03.2025
% Date last modified: 25.09.2026
% MATLAB version: 2024b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% select example case

% 'homogeneous' - homogeneous atmosphere
% atmosphere = 'homogeneous'; % <-- uncomment here for homogeneous atmosphere

% 'sounding'    - atmosphere imported from atmospheric soundings
atmosphere = 'sounding'; % <-- uncomment here for inhomogeneous atmosphere

atmosphere = validatestring(atmosphere, {'homogeneous', 'sounding'});

%% setup

% folder of this script (paths are built from here, independent of MATLAB's current folder)
script_folder = fileparts(mfilename('fullpath'));

% input data folder (contains the .dat files)
core_path = fullfile(script_folder, 'input_data');

% initializing file of the selected example case
switch atmosphere
    case 'homogeneous'
        ini_name = 'input_file_auralization.ini';
    case 'sounding'
        ini_name = 'input_file_auralization_import_atmosphere.ini';
end

input_file = fullfile(core_path, ['example_' atmosphere '_atmosphere'], ini_name);

% case tag. The 'tag' is used to name the folder where results and plots will be saved.
% All output files (plots, data) are also saved using this tag in their filenames.
tag = ['VR_approach_' atmosphere '_atmosphere'];

% results folder
results_path = [fullfile(script_folder, ['output_data_test_' atmosphere '_atmosphere']) filesep];

%% run auralization (see POSSIBLE CALLS in the header for other options)

auralization_master(core_path, tag, input_file, results_path);
