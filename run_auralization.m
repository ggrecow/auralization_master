% script run_auralization.m
%
% Gives an example on how to call the auralization framework
%
% INPUT DATA:
%
%   - exemplary case available under the '\input_data' folder 
%
%   - Generic flyover of an A319 during approach procedure
%
%   - has only one receiver, positioned at (x,y,z) = (0, 0, 1.2) meters
%
%   auralization_input.dat - sound emissions per time step, per sound
%   source. File needs to have this name.
%
%   geschw_hoehe_verlauf.dat - flight trajectory and associated operational
%   conditions. File needs to have this name.
%   
%   input_file_auralization.ini - initializing file containing some inputs
%   related to signal processing and setup of atmospherical conditions used
%   for sound propagation simulation 
%
%   input_file_auralization_import_atmosphere.ini - initializing file, same
%   as above, but is provided as an extra example to show how to import 
%   atmospheric soundings as .txt file from 
%   http://weather.uwyo.edu/upperair/sounding_legacy.html
%
% Assumption: all required inputs are available  
%
% -------------------------------
% Author: Gil Felix Greco (ggrecow@gmail.com)
% Institution: Technische Universität Braunschweig 
%
% Date created: 13.03.2025
% Date last modified: 29.04.2026
% MATLAB version: 2024b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% setup 

% input data folder
core_path = [pwd '\input_data'];

% input_file folder
input_file = [ core_path '\input_file_auralization.ini' ] ; 

% case tag. The 'tag' is used to name the folder where results and plots will be saved.
% All output files (plots, data) are also saved using this tag in their filenames.
tag = 'VR_approach';  

%% Ex. 1 - call auralization_master(core_path)

% Only <core_path> is provided. Then, tag = 'auralization_results', 
% and the code looks for an .ini inside <core_path>. Results are stored at <core_path> 

% WARNING: this setup only works if you have only one .ini file inside <core_path>

% auralization_master(core_path);

%% Ex. 2 - call auralization_master(core_path, tag)

% Here, <core_path> and <tag> are provided. The code will look for an 
% input_file inside <core_path>. Results are stored at <core_path> 

% WARNING: this setup only works if you have only one .ini file inside <core_path>

% auralization_master(core_path, tag); 

%% Ex. 3 - call auralization_master(core_path, tag, input_file)

% Here, <core_path>, <tag>, and <input_file> are provided. 
% Results are stored at <core_path> 

% auralization_master(core_path, tag, input_file); 

%% Ex. 4 -call auralization_master(core_path, tag, input_file, results_path)

% In this case, all inputs are provided. Results are stored at <results_path>. 
% HINT: <results_path> will be automatically created if non existent  

results_path = [pwd '\output_data_test\'];
auralization_master(core_path, tag, input_file, results_path); 