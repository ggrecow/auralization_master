clear all; close all; clc;

% input data
core_path = [pwd '\TEST_input_data'];

% input_file
input_file = [ core_path '\input_file_auralization.ini' ] ;

% case tag
tag = 'VR_app'; % tag for plots and saving results/figs

auralization_master(core_path);
% % auralization_master(core_path, tag);
% auralization_master(core_path, tag, input_file);
