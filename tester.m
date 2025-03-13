clear all; close all; clc;
% script tester.m
%
%   Gives and example on how to call the auralization framework
%
%   - Input data: generic flyover of an A319 during approach procedure
%   (therefore, tehre are no (buzzsaw and fan) tonal components from the engine noise) 
%
%   - Only one receiver, positioned at (x,y,z) = (0, 0, 1.2) meters
%
% Gil Felix Greco, Braunschweig 13.03.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% input data
core_path = [pwd '\TEST_input_data'];

% input_file
input_file = [ core_path '\input_file_auralization.ini' ] ;

% case tag
tag = 'VR_app'; % tag for plots and saving results/figs

%% call <auralization_framework>

auralization_master(core_path); % no tag is provided, results folder will be called <auralization_results>
% % auralization_master(core_path, tag); % no 'input_file' provided, the code will look for an input_file inside <core_path>
% % auralization_master(core_path, tag, input_file); % no <results_path> provided, results will be saved inside <core_path>
