function [time_PANAM_retarded, time_PANAM_auralization, time]  = getAuralizationTime(input, input_type)
%
% function [time_PANAM_retarded, time_PANAM_auralization, time]  = getAuralizationTime(input, input_type)
%
% auralization based on receiver data needs to be carefully processed
% because retarded time does not have an fixed dt. Moreover, needs to start
% from zero. This function deals with:
%
% 1) getting retarded time and normalizing to start in t = zero
% output: <time_PANAM_retarded>
%
% 2) creating a time vector  with the same size as <time_PANAM_retarded>
% and a desired fixed dt
% output: <time_PANAM_auralization>
%
% 3) creating a time vector for auralization based on <time_PANAM_auralization>
% and an input sampling frequency
%
% OBS: for auralization purposes, we need to interpolate the data to be
% auralized later on based on the original <time_PANAM_retarded> and <time_PANAM_auralization>
%
% Author: Gil Felix Greco, Braunschweig 26.09.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global fs

global dt_panam

%% normalize PANAM time vector (based on retarded time - normalized to start in zero)

switch input_type

    % get time vector from PANAM results
    case 'immission'
        for i = 1:size(input,1)
            time_PANAM(i) = input{i}.retarded_time; % get time vector; obs: based on the source time, change if for some reason you need receiver time
        end

    case 'emission'
        for i = 1:size(input,1)
            time_PANAM(i) = input{i}.source_time; % get time vector; obs: based on the source time, change if for some reason you need receiver time
        end

end

% get retarded time vector from PANAM results
for i = 1:length(time_PANAM)
    if i==1
        dt_panam_retarded(i) = time_PANAM(i) - time_PANAM(i);
    else
        dt_panam_retarded(i) = time_PANAM(i) - time_PANAM(i-1);
    end
end

time_PANAM_retarded = cumsum(dt_panam_retarded);

%% create vector with uniform dt

% for i = 1:size(input,1)
%     time_PANAM(i) = input{i}.source_time; % get time vector; obs: based on the source time, change if for some reason you need receiver time
% end
%
% dt_panam =  time_PANAM(2) - time_PANAM(1); % time resolution from panam. Assumes it is constant (dt from PANAM is usually .5 sec but not always constant). this has to be truncated, because if dt is not constant, we would have problems to reconstruct the broadband signal using the overlap and add method
%

% dt_panam = 0.5; % desired time resolution

time_PANAM_auralization = 0:dt_panam:(length(time_PANAM_retarded)*dt_panam)-dt_panam; % make time vector initiating in zero, with same length as original input and <dt_panam> timesteps

% time_PANAM_retarded = time_PANAM_auralization;

%% Signal temporal characteristics (for auralization)

dt = 1/fs;    % time discretization of the (to be) syntesized signal
% fmax = fs/2;  % max freq of the output auralized signal

time = time_PANAM_auralization(1):dt:time_PANAM_auralization(end);      % Time vector

if mod(length(time),2)~=0   % we need an even number of samples
    time = time_PANAM_auralization(1):dt:time_PANAM_auralization(end)+dt;       % add another sample to time vector
end

end
