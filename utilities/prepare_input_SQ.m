function [output_1, output_2, output_3, output_4] = prepare_input_SQ(input_1, input_2, input_3, input_4, TrimTime, tag_auralization)
%
% function [output_1, output_2] = prepare_input_SQ(input_1, input_2, TrimTime)
%
% this function prepares the input data for both objective SQ calculation 
% and auralizations. This is necessary when we want to compare metrics calculated
% from 1) directly from PANAM and 2) from the auralized sound files. Also may be necessary 
% when using ray-tracing based auralization, because shadow zones may occur because of 
% sound refraction. Thus, trimming the aircraft trajectory around the max. SPL may avoid this issue  
%
% Procedure:
%
%   1) check for doubled values on the time vector and remove them
%
%   2) uses function <trim_data_TimeBased> to trim the databased on a desired <TrimTime>.
%   type <help trim_data_TimeBased> for more details
%
% INPUTS:
%
% input_1 - input data to be trimmed. Only use a <data> output from the 
%           <PANAM_SQAT_data_conversion> because this function was
%           designed to work with this data
%
% input_2 - input data to be trimmed. Only use a <SPECTROGRAM> output from the 
%           <PANAM_SQAT_data_conversion> because this function was
%           designed to work with this data
%
% input_3 - input data to be trimmed. Only use a <SPECTROGRAM> output from the 
%           <PANAM_SQAT_data_conversion> because this function was
%           designed to work with this data
%
% input_4 - Flight profile (see <get_flight_profile> function)
%
% TrimTime - time span to trim the data before and after the max spl
%
% OUTPUTS:
%
% output_1 - <data> output already prepared to use on for SQ analysis and auralization
%
% output_2 - <SPECTROGRAM> output already prepared to use on for SQ analysis and auralization
% 
% output_3 - <SPECTROGRAM> output already prepared to use on for SQ analysis and auralization
%
% output_4 - trimmed flight profile
%
% Gil Felix Greco - Braunschweig 22.09.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check and correct input data for doubled input values on time vector

% procedure on input_1
time_PANAM = zeros(1,length(input_1)); % declare variable for memory allocation

for i = 1:size(input_1,1)
    time_PANAM(i) = input_1{i}.source_time; % get time vector; obs: based on the source time, change if for some reason you need receiver time
end

% it may happens that PANAM inputs have repeated time (and thus noise) bins, this needs to be corrected (i.e. repeated idx removed)
[~, ind] = unique(time_PANAM); % xind = index of first occurrence of a repeated value 

input_unique = cell(length(ind),1); % declare variable for memory allocation
for i = 1:length(ind)
    input_unique{i,1} = input_1{ind(i)}; % keep only not repeated time instances on <input> variable
end

for i = 1:length(ind)

    % procedure on input_2
    input_2{1,1}.source_time(:,i) = input_2{1,1}.source_time(:,ind(i));
    input_2{1,1}.retarded_time(:,i) = input_2{1,1}.retarded_time(:,ind(i));
    
    input_2{1,1}.overall(:,i) = input_2{1,1}.overall(:,ind(i));
    input_2{1,1}.overall_broadband(:,i) = input_2{1,1}.overall_broadband(:,ind(i));
    input_2{1,1}.airframe_toc(:,i) = input_2{1,1}.airframe_toc(:,ind(i));
    input_2{1,1}.engine(:,i) = input_2{1,1}.engine(:,ind(i));
    input_2{1,1}.engine_without_fan_harmonics(:,i) = input_2{1,1}.engine_without_fan_harmonics(:,ind(i));
    input_2{1,1}.engine_broadband_toc(:,i) = input_2{1,1}.engine_broadband_toc(:,ind(i));
    input_2{1,1}.engine_buzzsaw_toc(:,i) = input_2{1,1}.engine_buzzsaw_toc(:,ind(i));
    input_2{1,1}.fan_harmonics_toc(:,i) = input_2{1,1}.fan_harmonics_toc(:,ind(i));
    
    % procedure on input_3
    input_3{1,1}.source_time(:,i) = input_3{1,1}.source_time(:,ind(i));
    input_3{1,1}.retarded_time(:,i) = input_3{1,1}.retarded_time(:,ind(i));
    
    input_3{1,1}.overall(:,i) = input_3{1,1}.overall(:,ind(i));
    input_3{1,1}.overall_broadband(:,i) = input_3{1,1}.overall_broadband(:,ind(i));
    input_3{1,1}.airframe(:,i) = input_3{1,1}.airframe(:,ind(i));
    input_3{1,1}.engine(:,i) = input_3{1,1}.engine(:,ind(i));
    input_3{1,1}.engine_without_fan_harmonics(:,i) = input_3{1,1}.engine_without_fan_harmonics(:,ind(i));
    input_3{1,1}.engine_broadband(:,i) = input_3{1,1}.engine_broadband(:,ind(i));
    input_3{1,1}.engine_buzzsaw(:,i) = input_3{1,1}.engine_buzzsaw(:,ind(i));
    input_3{1,1}.fan_harmonics(:,i) = input_3{1,1}.fan_harmonics(:,ind(i));
       
end

% procedure on input_4
input_4.x = input_4.x(ind);
input_4.y = input_4.y(ind);
input_4.z = input_4.z(ind);
input_4.TAS = input_4.TAS(ind);
input_4.thrust = input_4.thrust(ind);
                       
%% trim data to a specific time duration

% TrimTime = 20; % seconds, time after and before the max SPL on the input data

% [idx_lower, idx_upper] = trim_data_TimeBased(input_unique, TrimTime, input_4, tag_auralization);
[idx_lower, idx_upper] = trim_data_TimeDistanceBased(input_unique, TrimTime, input_4, tag_auralization);
 
% trim input_1
output_1 = input_unique(idx_lower:idx_upper,:);

% trim input_2
output_2{1,1}.source_time = input_2{1,1}.source_time(:,idx_lower:idx_upper);
output_2{1,1}.retarded_time = input_2{1,1}.retarded_time(:,idx_lower:idx_upper);

output_2{1,1}.overall = input_2{1,1}.overall(:,idx_lower:idx_upper);
output_2{1,1}.overall_broadband = input_2{1,1}.overall_broadband(:,idx_lower:idx_upper);
output_2{1,1}.airframe_toc = input_2{1,1}.airframe_toc(:,idx_lower:idx_upper);
output_2{1,1}.engine = input_2{1,1}.engine(:,idx_lower:idx_upper);
output_2{1,1}.engine_without_fan_harmonics = input_2{1,1}.engine_without_fan_harmonics(:,idx_lower:idx_upper);
output_2{1,1}.engine_broadband_toc = input_2{1,1}.engine_broadband_toc(:,idx_lower:idx_upper);
output_2{1,1}.engine_buzzsaw_toc = input_2{1,1}.engine_buzzsaw_toc(:,idx_lower:idx_upper);
output_2{1,1}.fan_harmonics_toc = input_2{1,1}.fan_harmonics_toc(:,idx_lower:idx_upper);

% trim input_3
output_3{1,1}.source_time = input_3{1,1}.source_time(:,idx_lower:idx_upper);
output_3{1,1}.retarded_time = input_3{1,1}.retarded_time(:,idx_lower:idx_upper);

output_3{1,1}.overall = input_3{1,1}.overall(:,idx_lower:idx_upper);
output_3{1,1}.overall_broadband = input_3{1,1}.overall_broadband(:,idx_lower:idx_upper);
output_3{1,1}.airframe = input_3{1,1}.airframe(:,idx_lower:idx_upper);
output_3{1,1}.engine = input_3{1,1}.engine(:,idx_lower:idx_upper);
output_3{1,1}.engine_without_fan_harmonics = input_3{1,1}.engine_without_fan_harmonics(:,idx_lower:idx_upper);
output_3{1,1}.engine_broadband = input_3{1,1}.engine_broadband(:,idx_lower:idx_upper);
output_3{1,1}.engine_buzzsaw = input_3{1,1}.engine_buzzsaw(:,idx_lower:idx_upper);
output_3{1,1}.fan_harmonics = input_3{1,1}.fan_harmonics(:,idx_lower:idx_upper);

% trim input_4
output_4.x = input_4.x(idx_lower:idx_upper);
output_4.y = input_4.y(idx_lower:idx_upper);
output_4.z = input_4.z(idx_lower:idx_upper);
output_4.TAS = input_4.TAS(idx_lower:idx_upper);
output_4.thrust = input_4.thrust(idx_lower:idx_upper);

end
