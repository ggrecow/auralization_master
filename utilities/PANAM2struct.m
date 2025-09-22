function [data,OASPL,OASPL_dBA,SPECTROGRAM,SPECTROGRAM_dBA] = PANAM2struct(PATH)
% function [data,OASPL,OASPL_dBA,SPECTROGRAM,SPECTROGRAM_dBA] = PANAM2struct(PATH)
%
% READ and manage freq dependent emission/immission data per time step 
% from PANAM 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  INPUT: PATH = string containing the directory path to .dat input data 
%                from PANAM - all data is provided in SPL (dB(Z)) values.
%                Values in dB(A) are computed within this code.
%
%  OUTPUTS: data -> struct containing all processed data per time step in   
%                   the format - data{nTimeSteps,nObserver} - contains
%                   db(Z) and dB(A) data
%
%           OASPL -> struct containing OASPL (time-dependent) values for 
%                    each individual noise sources - dB(Z) and dB(A) data 
%                    are separate outputs
%
%           SPECTROGRAM -> struct containing time, freq and SPL(nFreq,nTime) 
%                          vectors for each individual noise sources -
%                          dB(Z) and dB(A) data are separate outputs
%
%  Main steps:
%
%  1) PANAM immission data is loaded and divided in HEADER and DATA using 
%     function mhdrload
%     
%     HEADER: char containing info about data 
%     DATA: (31,2)=(nFreq,SPL(z)) - contribution for 4 different sound sources
%                                 1- Airframe broadband (1/3-OB)   
%                                 2- Engine broadband (no buzzsaw) (1/3-OB)   
%                                 3- Engine Buzzsaw noise (tones - actual freq.)
%                                   3.1 - Engine Buzzsaw noise (1/3-OB) 
%                                 4- Fan harmonics (tones - actual freq.)
%                                   4.1- Fan harmonics (1/3-OB) 
%
%  Important: 1 + 2 + 3.1 + 4.1 = overall
%                 1 + 2                  = overall broadband
%                       2 + 3.1 + 4.1 = Engine 
%              
%  Important: - the dimension z of HEADER and DATA is a sequence of 
%               time-steps 
%
%  2) Get data from HEADER(char), transform into double and store info of
%     each imission data as struct variable in the rows of a cell(data)
%
%  3) Get DATA from each noise source and their corresponding freq vectors and
%     stores together with their correspondent emission/immission (i.e. SPL) info in each cell(data)
%     rows. For the first three sources (airframe, engine bbn and engine buuzsaw), 
%     the following steps are carried out:
%
%   - Step 1) struct in the format data{nTimeStep,obs}.SourceName(nFreq,nLevel)
%   - Step 2) freq vector may not correspond to the correct central freq of
%             1/3 octave bands, so this is corrected
%   - Step 3) initial freq and final freq may not be the same, so all data 
%             is corrected to have the same 28 nFreqBands 
%             i.e., 1st band (25Hz) to 28th band (12.5kHz)
%
%   The fan harmonics is a special guy because they are given in their 
%   actual corresponding narrow band freq. The following steps are adopted:
%
%       - Step 1) struct in the format data{nTimeStep,obs}.SourceName(nFreq,nLevel) 
%                 (freq bands correspond to narrow bands)
%       - Step 2) truncate negative/NaN/inf level values to be zero
%       - Step 3) an extra variable is created to include the tons on the 
%                 nearest corresponding 1/3 octave band
%
%  4) Organize the cell(data), which contains n rows with sequential info
%  about emission, per observers, so that the output is the cell variable
%  data_observer(nTime,nobserver) 
%
%  Gil Greco, Braunschweig, 16/04/2020
%  updated: 21.04.2021
%  udpated: 09.06.2023 (input data is not restricted anymore, just needs to have a valid format .dat)
%  udpated: 11.10.2023 (buzzsaw is now treated as tonal and 1/3-OB separetly. overall_broadband noise does not include buzzsaw anymore)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%  udpated: 19.06.2024 (fundamental sources (i.e. airframe BBN, engine BBN,
%  engine buzssaw and engine fan tonal) are now given as raw (i.e. as from
%  PANAM) and in TOC (from 25 Hz to 12.5 kHz). For those fundamental
%  sources, the TOC or nothing is on the variable. All other metrics
%  composed by these fundamental sources are already in TOC, including all
%  dBA value, so the TOC nomenclature is not included but inferred).
%  SPECTROGRAMS ONLY TOC !!!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) load data
%
% Header --> contains char data info about each time step (see function load_header)
% DATA   --> contains emission data (1st column = freq bands; 2nd column = emission SPL)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read input data;
[HEADER,DATA] = mhdrload(PATH);

%% 2) Get data from header and store in a cell 
% (each cell contains a structure with all data emission/immission for a time-step)

i = 1:4:size(HEADER,3);     % 3 in 3 because there is 3 different noise source contributions
data = cell(length(i),1);   % pre-alocatting cell size for memory
clear i

aux = 1;   % initiating auxiliar counter

% each cell row --> one time step with corresponding emission/immission info and data 

for i = 1:4:size(HEADER,3)
    
    if i==1      % use other function because 1st header is special (has 2 extra rows in the beginning)
        
        a = HEADER(1:size(HEADER,1),1:size(HEADER,2),i);
        [data{aux}] = load_header_first(a);
        aux = aux+1;
        clear a
        
    else
        
        a = HEADER(1:size(HEADER,1),1:size(HEADER,2),i);
        [data{aux}] = load_header(a);
        aux = aux+1;
        clear a
        
    end
end

%% 1/3 octave central bands 
% (freq vector data from panam may not be already in the correct central freq of 1/3 octrave bands)   

% USE this (entire toc bands from 16 Hz - 25 khz) - PANAM input will be checked in the following and cutted only to 25 Hz - 12.5 kHz  
toc = [16 20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000 25000];

% toc = [25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 ];

%% A weighting factor in 1/3 octave central bands    

A_weighting_factor = [20,-50.5000000000000;25,-44.7000000000000;32,-39.4000000000000;40,-34.6000000000000;50,-30.2000000000000;63,-26.2000000000000;80,-22.5000000000000;100,-19.1000000000000;125,-16.1000000000000;160,-13.4000000000000;200,-10.9000000000000;250,-8.60000000000000;315,-6.60000000000000;400,-4.80000000000000;500,-3.20000000000000;630,-1.90000000000000;800,-0.800000000000000;1000,0;1250,0.600000000000000;1600,1;2000,1.20000000000000;2500,1.30000000000000;3150,1.20000000000000;4000,1;5000,0.500000000000000;6300,-0.100000000000000;8000,-1.10000000000000;10000,-2.50000000000000;12500,-4.30000000000000;16000,-6.60000000000000;20000,-9.30000000000000];

% load A_weighting_factor         % (1st column -> tob from 20-20kHz)
                                  % (2nd column -> A-weight factor for each tob)

A_weighting_factor = A_weighting_factor(2:29,:); % only from 25 Hz tpo 12.5 kHz                                         
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) get data from each noise source and its corresponding freq vector
%  3.1) airframe data
%       Step 1 - read raw input data and store in a struct of the format:
%       data(nTimeSteps,nObserver).NoiseSource(nFreq,nLevel)  

aux = 1;
for i = 1:4:size(DATA,3)

    data{aux}.airframe = DATA(:,:,i);
    data{aux}.airframe(isnan(data{aux}.airframe)) = 0; % truncate NaN values
    data{aux}.airframe(isinf(data{aux}.airframe)) = 0; % truncate inf values
    data{aux}.airframe(data{aux}.airframe<0) = 0; % truncate negative values
    aux = aux+1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify if values are the same (error: Orignal data (Lothar) - processed data)
% error_1=data{1}.airframe(:,2)-DATA(:,2,1);                             % error: all values should be zero
% error_2=data{size(data,1)}.airframe(:,2)-DATA(:,2,size(DATA,3)-3);     % error: all values should be zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.1) airframe data
%       Step 2 - checking/correcting toc freq bands

% brute force approach 

for j = 1:size(data,1)  % sweeps dat struct rows (each row correspond to a time steps)
    
    % 1) find which idx corresponds to the 25 Hz freq band
    % fmin=25 Hz
    b = abs( toc(3)-data{j}.airframe(:,1) );  % error vector
    [idx_min] = find(b==min(b));              % find idx of min error
    clear b
    
    % 2) cut to have first idx always related to 25 Hz
    a = data{j}.airframe(idx_min:idx_min+27,:);
    data{j}.airframe_toc = a;
    
    % 3) impose 1/3 octave bands till 1.25 kHz
    for k=1:28
        data{j}.airframe_toc(k,1) = transpose(toc(1,k+2));
    end
    
    clear idx_max idx_min a
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.1) airframe data
%       Step 3 - A-weighting

for j = 1:size(data,1)  % sweeps dat struct rows (each row correspond to a time steps)
    
    data{j}.airframe_dBA = zeros(size(data{j}.airframe_toc,1),2);  % pre allocating memory
    
    for k = 1:28
        
        data{j}.airframe_dBA(k,1) = A_weighting_factor(k,1);
        
        aux(j) = data{j}.airframe_toc(k,2) + A_weighting_factor(k,2);
        
        if aux(j)<0
            data{j}.airframe_dBA(k,2) = 0;
        else
            data{j}.airframe_dBA(k,2) = aux(j);
        end
        
    end
    
    clear aux
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Get data from each noise source and its corresponding freq vector
%  3.2) engine broadband data
%       Step 1 - read raw input data and store in a struct in the format:
%       data(nTimeSteps,nObserver).NoiseSource(nFreq,nLevel)  

aux = 1;

for i = 2:4:size(DATA,3)

    data{aux}.engine_broadband = DATA(:,:,i);
    data{aux}.engine_broadband(isnan(data{aux}.engine_broadband)) = 0; % truncate NaN values
    data{aux}.engine_broadband(isinf(data{aux}.engine_broadband)) = 0; % truncate inf values
    data{aux}.engine_broadband(data{aux}.engine_broadband<0) = 0; % truncate negative values
    aux = aux+1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify if values are the same (error: Orignal data (Lothar) - processed data)
% error_1=data{1}.engine_broadband(:,2)-DATA(:,2,2);                             % error: all values should be zero
% error_2=data{size(data,1)}.engine_broadband(:,2)-DATA(:,2,size(DATA,3)-2);     % error: all values should be zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.2) engine broadband data
%       Step 2 - checking/correcting toc freq bands

for j = 1:size(data,1)  % sweeps dat struct rows (each row correspond to a time steps)
    
    % 1) find which idx corresponds to the 25 Hz freq band
    % fmin=25 Hz
    b = abs( toc(3)-data{j}.engine_broadband(:,1) );  % error vector
    [idx_min] = find(b==min(b));            % find idx of min error
    clear b
    
    % 2) cut to have first idx always related to 25 Hz
    a = data{j}.engine_broadband(idx_min:idx_min+27,:);
    data{j}.engine_broadband_toc = a;
    
    % 3) impose 1/3 octave bands till 1.25 kHz
    for k = 1:28
        data{j}.engine_broadband_toc(k,1) = transpose(toc(1,k+2));
    end
    
    clear idx_max idx_min a
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.2) engine broadband data
%       Step 3 -  A weighting

for j=1:size(data,1)  % sweeps dat struct rows (each row correspond to a time steps)
    
    data{j}.engine_broadband_dBA = zeros(size(data{j}.engine_broadband_toc,1),2);  % pre allocating memory
    
    for k = 1:28
        
        data{j}.engine_broadband_dBA(k,1) = A_weighting_factor(k,1);
        
        aux(j) = data{j}.engine_broadband_toc(k,2) + A_weighting_factor(k,2);
        
        if aux(j)<0
            data{j}.engine_broadband_dBA(k,2) = 0;
        else
            data{j}.engine_broadband_dBA(k,2) = aux(j);
        end
        
    end
    
    clear aux
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) get data from each noise source and its corresponding freq vector
%  3.3) engine Buzzsaw data
%       Step 1 - read raw input data and store in a struct data(nTimeSteps,nObserver).<noise_source>(nFreq,nLevel)  

aux = 1;
for i = 3:4:size(DATA,3)

    data{aux}.engine_buzzsaw = DATA(:,:,i);
    data{aux}.engine_buzzsaw(isnan(data{aux}.engine_buzzsaw)) = 0; % truncate NaN values
    data{aux}.engine_buzzsaw(isinf(data{aux}.engine_buzzsaw)) = 0; % truncate inf values
    data{aux}.engine_buzzsaw(data{aux}.engine_buzzsaw<0) = 0; % truncate negative values
    aux = aux + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify if values are the same (error: Orignal data (Lothar) - processed data)
% error_1=data{1}.engine_buzzsaw(:,2)-DATA(:,2,3);                             % error: all values should be zero
% error_2=data{size(data,1)}.engine_buzzsaw(:,2)-DATA(:,2,size(DATA,3)-1);     % error: all values should be zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.3) engine Buzzsaw data
%       Step 2 - checking/correcting toc freq bands

for j = 1:size(data,1)  % sweeps dat struct rows (each row correspond to a time steps)
    
    % 1) find which idx corresponds to the 25 Hz freq band
    % fmin=25 Hz
    b = abs( toc(3)-data{j}.engine_buzzsaw(:,1) );  % error vector
    [idx_min] = find(b==min(b));            % find idx of min error
    clear b

    % 2) cut to have first idx always related to 25 Hz
    a = data{j}.engine_buzzsaw(idx_min:idx_min+27,:);
    data{j}.engine_buzzsaw_toc = a;
    
    % 3) impose 1/3 octave bands till 1.25 kHz
    for k = 1:28
        data{j}.engine_buzzsaw_toc(k,1) = transpose(toc(1,k+2));
    end
    
    clear idx_max idx_min a
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.3) engine Buzzsaw data
%       Step 3 - A weighting

for j = 1:size(data,1)  % sweeps dat struct rows (each row correspond to a time steps)
    
    data{j}.engine_buzzsaw_dBA = zeros(size(data{j}.engine_buzzsaw_toc,1),2);  % pre allocating memory
    
    for k = 1:28
        
        data{j}.engine_buzzsaw_dBA(k,1) = A_weighting_factor(k,1);
        
        aux(j) = data{j}.engine_buzzsaw_toc(k,2) + A_weighting_factor(k,2);
        
        if aux(j)<0
            data{j}.engine_buzzsaw_dBA(k,2) = 0;
        else
            data{j}.engine_buzzsaw_dBA(k,2) = aux(j);
        end
        
    end
    
    clear aux
end
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3) get data from each noise source and its corresponding freq vector
%  3.4) fan harmonics
%       Step 1 - read raw input data and store in a struct of the format:
%                data(nTimeSteps,nObserver).NoiseSource(nFreq,nLevel)  

aux = 1;

for i = 4:4:size(DATA,3)

    data{aux}.fan_harmonics = DATA(1:28,:,i);
    aux = aux + 1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify if values are the same (error: Orignal data (Lothar) - processed data)
% clear error_1 error_2
% error_1=data{1}.fan_harmonics(:,2)-DATA(:,2,4);                          % error: all values should be zero
% error_2=data{size(data,1)}.fan_harmonics(:,2)-DATA(:,2,size(DATA,3));    % error: all values should be zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.4) fan harmonics
%       Step 2 - force negative values to be zero

for j = 1:size(data,2)   % observer loop
    for i = 1:size(data,1)   % time loop
        
        a = data{i,j}.fan_harmonics(:,2);
        
        a(a<0) = 0; % truncate negative values
        a(isinf(a)) = 0; % truncate inf values
        a(isnan(a)) = 0; % truncate NaN values
        
        data{i,j}.fan_harmonics(:,2) = a;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clear a;
        % a =  data{i,j}.fan_harmonics(:,1)~=0 & ... % find idx where freq is not zero
        %      data{i,j}.fan_harmonics(:,2)==0;      % find idx where SPL is zero
        %
        % data{i,j}.fan_harmonics(a,1) = 0;  % truncate freq of tones with 0dB from vector (14.06.2023)
        %
        % %not good approach for immission data, since tones may be zero
        % %at a certain time step and not zero after. This would lead to a
        % %vector of tones with different sizes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end

clear a

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.4) fan harmonics
%       Step 3.1 - create fan harmonics vector in 1/3 octave bands

freq = data{1,1}.airframe_toc(:,1);  % toc bands
nFreq_toc = length(data{1,1}.airframe_toc(:,1)); % n toc bands

% create variable fan_harmonics_toc
fan_harmonics_toc = zeros(length(freq),2);
fan_harmonics_toc(:,1) = freq;

for j = 1:size(data,2)   % observer loop
    for i = 1:size(data,1)   % time loop
        
        data{i,j}.fan_harmonics_toc = fan_harmonics_toc;
        
    end
end

clear fan_harmonics_toc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.4) fan harmonics
%       Step 3.2 - adding fan harmonics to the corresponding (nearest) 1/3 octave band

for j = 1:size(data,2)   % observer loop
    for i = 1:size(data,1)   % time loop
        
        n_harmonics = size(nonzeros(data{i,j}.fan_harmonics(:,1)),1); % find number of harmonics from data
        
        for w = 1:n_harmonics    % loop harmonics
            
            f_target = data{i,j}.fan_harmonics(w,1);
            A(:) = abs(freq(:,1)-f_target);
            B = min(A);
            
            decide = find(A==B);  % avoid singular values that are right in the middle of two TOC bands
            
            if length(decide)>1
                decide = decide(1);
            else
            end
            
            [idx(:,w)] = decide;
            
        end
        
        s = 1; % aux counter
        
        for k = 1:nFreq_toc   % loop freq tob
            
            for y = 1:length(idx)
                
                % if k==idx(1) || k==idx(2) || k==idx(3) || k==idx(4) || k==idx(5)
                if k==idx(y)
                    
                    data{i,j}.fan_harmonics_toc(k,2) = data{i,j}.fan_harmonics(s,2);
                    
                    s = s + 1;
                    
                else
                    
                end
                
            end
        end
        
    end
end

clear idx A B s w i j y f_target n_harmonics k a 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  3.4) fan harmonics
%       Step 4 - A weighting

for j = 1:size(data,1)  % sweeps dat struct rows (each row correspond to a time steps)
    
    data{j}.fan_harmonics_dBA = zeros(size(data{j}.engine_broadband_toc,1),2);  % pre allocating memory
    
    for k = 1:28
        
        data{j}.fan_harmonics_dBA(k,1) = A_weighting_factor(k,1);
        
        aux(j) = data{j}.fan_harmonics_toc(k,2);
        
        if aux(j)>0
            data{j}.fan_harmonics_dBA(k,2) = aux(j) + A_weighting_factor(k,2);
        else
            data{j}.fan_harmonics_dBA(k,2) = 0;
        end
        
    end
    
    clear aux
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  checkpoint (all original input data is contained inside the output <data> struct.
%               Till here, the output <data> struct is merely a huge row vector,
%               not discriminatating if there is data for more than one receiver position. 
%               From now onwards, additional SPL metrics and vectors are 
%               computed using the original input data, and the output <data> struct
%               is re-organized if the number of receiver positions ir greater than 1) 

clear HEADER DATA aux toc A_weighting_factor j k   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% summing noise sources to obtain emission/immission levels(freq) per time step 

for j = 1:size(data,2)   % observer loop
    for i = 1:size(data,1)   % time loop
        
        % declaring variables for memory allocation
        data{i,j}.engine = zeros(length(nFreq_toc),2);
        data{i,j}.engine_without_fan_harmonics = zeros(length(nFreq_toc),2);
        data{i,j}.overall_broadband = zeros(length(nFreq_toc),2);
        data{i,j}.overall = zeros(length(nFreq_toc),2);
        
        data{i,j}.engine_dBA = zeros(length(nFreq_toc),2);
        data{i,j}.engine_without_fan_harmonics_dBA = zeros(length(nFreq_toc),2);
        data{i,j}.overall_broadband_dBA = zeros(length(nFreq_toc),2);
        data{i,j}.overall_dBA = zeros(length(nFreq_toc),2);
        
        for k = 1:nFreq_toc
            
            % filling freq column with toc bands
            data{i,j}.engine(k,1) = freq(k,1);                         
            data{i,j}.engine_without_fan_harmonics(k,1) = freq(k,1);  
            data{i,j}.overall_broadband(k,1) = freq(k,1);              
            data{i,j}.overall(k,1) = freq(k,1);                        
                        
            data{i,j}.engine_dBA(k,1)= freq(k,1);                                  
            data{i,j}.engine_without_fan_harmonics_dBA(k,1)= freq(k,1);            
            data{i,j}.overall_broadband_dBA(k,1)= freq(k,1);                       
            data{i,j}.overall_dBA(k,1)= freq(k,1);                       
            
            % ENGINE (total)
            data{i,j}.engine(k,2) = ...
                10.*log10( 10.^(data{i,j}.fan_harmonics_toc(k,2)./10)...
                         + 10.^(data{i,j}.engine_broadband_toc(k,2)./10)...
                         + 10.^(data{i,j}.engine_buzzsaw_toc(k,2)./10) );
            
            data{i,j}.engine_dBA(k,2) = ...
                10.*log10( 10.^(data{i,j}.fan_harmonics_dBA(k,2)./10)...
                         + 10.^(data{i,j}.engine_broadband_dBA(k,2)./10)...
                         + 10.^(data{i,j}.engine_buzzsaw_dBA(k,2)./10) );
            
            % ENGINE (BBN only, i.e. without fan harmonics) - deprecated
            % but still here because aures tonality on the panam module
            % uses it (24.11.2023)
            % because buzzsaw is also tonal
            data{i,j}.engine_without_fan_harmonics(k,2) = ...
                10.*log10( 10.^(data{i,j}.engine_broadband_toc(k,2)./10)...
                         + 10.^(data{i,j}.engine_buzzsaw_toc(k,2)./10) );
            
            data{i,j}.engine_without_fan_harmonics_dBA(k,2) = ...
                10.*log10( 10.^(data{i,j}.engine_broadband_dBA(k,2)./10)...
                         + 10.^(data{i,j}.engine_buzzsaw_dBA(k,2)./10) );
           
            % OVERALL (BROADBAND NOISE ONLY)
            data{i,j}.overall_broadband(k,2) = ...
                10.*log10( 10.^(data{i,j}.airframe_toc(k,2)./10)...
                         + 10.^(data{i,j}.engine_broadband_toc(k,2)./10) );
            
            data{i,j}.overall_broadband_dBA(k,2) = ...
                10.*log10( 10.^(data{i,j}.airframe_dBA(k,2)./10)...
                         + 10.^(data{i,j}.engine_broadband_dBA(k,2)./10) );

%             data{i,j}.overall_broadband(k,2) = ...
%                          10.*log10( 10.^(data{i,j}.airframe(k,2)./10)...
%                          + 10.^(data{i,j}.engine_broadband(k,2)./10)...
%                          + 10.^(data{i,j}.engine_buzzsaw_toc(k,2)./10) );
%                      
%             data{i,j}.overall_broadband_dBA(k,2) = ...
%                 10.*log10( 10.^(data{i,j}.airframe_dBA(k,2)./10)...
%                          + 10.^(data{i,j}.engine_broadband_dBA(k,2)./10)...
%                          + 10.^(data{i,j}.engine_buzzsaw_toc_dBA(k,2)./10) );
                         
            % OVERALL 
            data{i,j}.overall(k,2) = ...
                10.*log10( 10.^(data{i,j}.engine(k,2)./10)...
                         + 10.^(data{i,j}.airframe_toc(k,2)./10) );
            
            data{i,j}.overall_dBA(k,2) = ...
                10.*log10( 10.^(data{i,j}.engine_dBA(k,2)./10)...
                         + 10.^(data{i,j}.airframe_dBA(k,2)./10) );
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OASPL of different noise sources per each time step

for j = 1:size(data,2)   % observer loop
    for i = 1:size(data,1)   % time loop
        
        data{i,j}.OASPL_airframe = 10.*log10(sum(10.^(data{i,j}.airframe(:,2)./10)));
        data{i,j}.OASPL_airframe_toc = 10.*log10(sum(10.^(data{i,j}.airframe_toc(:,2)./10)));
        data{i,j}.OASPL_fan_harmonics= 10.*log10(sum(10.^(data{i,j}.fan_harmonics(:,2)./10)));
        data{i,j}.OASPL_fan_harmonics_toc = 10.*log10(sum(10.^(data{i,j}.fan_harmonics_toc(:,2)./10)));
        data{i,j}.OASPL_engine_broadband = 10.*log10(sum(10.^(data{i,j}.engine_broadband(:,2)./10)));
        data{i,j}.OASPL_engine_broadband_toc = 10.*log10(sum(10.^(data{i,j}.engine_broadband_toc(:,2)./10)));
        data{i,j}.OASPL_engine_buzzsaw = 10.*log10(sum(10.^(data{i,j}.engine_buzzsaw(:,2)./10)));
        data{i,j}.OASPL_engine_buzzsaw_toc = 10.*log10(sum(10.^(data{i,j}.engine_buzzsaw_toc(:,2)./10)));
        data{i,j}.OASPL_engine = 10.*log10(sum(10.^(data{i,j}.engine(:,2)./10)));
        data{i,j}.OASPL_engine_without_fan_harmonics = 10.*log10(sum(10.^(data{i,j}.engine_without_fan_harmonics(:,2)./10)));
        data{i,j}.OASPL_overall_broadband = 10.*log10(sum(10.^(data{i,j}.overall_broadband(:,2)./10)));
        data{i,j}.OASPL_overall = 10.*log10(sum(10.^(data{i,j}.overall(:,2)./10)));
        
        data{i,j}.OASPL_airframe_dBA = 10.*log10(sum(10.^(data{i,j}.airframe_dBA(:,2)./10)));
        data{i,j}.OASPL_fan_harmonics_dBA = 10.*log10(sum(10.^(data{i,j}.fan_harmonics_dBA(:,2)./10)));
        data{i,j}.OASPL_engine_broadband_dBA = 10.*log10(sum(10.^(data{i,j}.engine_broadband_dBA(:,2)./10)));
        data{i,j}.OASPL_engine_buzzsaw_dBA = 10.*log10(sum(10.^(data{i,j}.engine_buzzsaw_dBA(:,2)./10)));
        data{i,j}.OASPL_engine_dBA = 10.*log10(sum(10.^(data{i,j}.engine_dBA(:,2)./10)));
        data{i,j}.OASPL_engine_without_fan_harmonics_dBA = 10.*log10(sum(10.^(data{i,j}.engine_without_fan_harmonics_dBA(:,2)./10)));
        data{i,j}.OASPL_overall_broadband_dBA = 10.*log10(sum(10.^(data{i,j}.overall_broadband_dBA(:,2)./10)));
        data{i,j}.OASPL_overall_dBA = 10.*log10(sum(10.^(data{i,j}.overall_dBA(:,2)./10)));
        
    end
end

clear i j k nFreq_toc freq
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% organizing data per observer position 
% (checks how many receiver positions are available in the input file and 
% and organize it. Data for each receiver position is allocated in each 
% column of the output <data> struct, while each time step is allocated in the rows)

% source time approach
source_time = zeros(size(data,1),size(data,2)); % pre allocating for memory

for i = 1:size(data,1)   % time loop
    source_time(i) = data{i}.source_time;
end

% [idx]=find(source_time==0);       % find idx when source time is equal to zero
[idx] = find(source_time==source_time(1));       % find idx when source time is equal to zero
n_observer = size(idx,1);                        % estimating the total number of observers 
clear source_time;

if n_observer>1 % check if all observers have the same number of time-steps

    for i = 1:(n_observer-2)       % got the number of time-steps per observer-1

        nTime_a(:) = idx(i+1)-idx(i);
        nTime_b(:) = idx(i+2)-idx(i+1);

    end

    if nTime_a==nTime_b % this only works if all observers have the same number of time steps

        nTime = nTime_a;
        clear nTime_a nTime_b

        data_observer=cell(nTime,n_observer);      % declaring cell (nTimeSteps,nObserver)

        % starts a nasty logic to re-organized the output <data> struct according to the number of receiver positions!!!
        for j = 1:size(data_observer,2)              % column loop (nObservers)

            if j==1
                a = 1;
                b = nTime;
            else
                a = (j-1)*(nTime+1);
                b = j*nTime;
            end

            k = round(linspace(a,b,nTime));

            for i = 1:nTime      % time steps loop

                aux_data = data{k(i)};

                data_observer{i,j} = aux_data;

                clear aux_data
            end

        end

        data = data_observer; % the final output is the struct data(nTime,nObserver)

    else
        disp('WARNING: Observer positions dont have the same number of time-steps !!!');
        return
    end

    clear a b i j k nObs nTime data_observer

else

    %   return   % if only one observer is available, outputs already computed struct data(nTime,1)

end

%% check for input data for doubled input values on time vector (only for checking and displaying, correction cant be made now in order to keep the same time vector for all receiver positions)

for j = 1:size(data,2)     % observer loop
    
    for i = 1:size(data,1) % time loop
        time_PANAM(i,j) = data{i,j}.source_time; % get time vector; obs: based on the source time, change if for some reason you need receiver time
    end
    
    % it may happens that PANAM inputs have repeated time (and thus noise) bins, this needs to be corrected (i.e. repeated idx removed)
    [~, ind] = unique(time_PANAM(:,j)); % xind = index of first occurrence of a repeated value
    unique_length(j) = length(ind); % store number of unique time values per j observer
    clear ind
    
    if unique_length(j) == size(data(:,j),1)
        doubled_time(j) = 0; % <double_time> receives 0 if no repeated time values were found for the j-th receiver 
    else
        doubled_time(j) = 1; % <double_time> receives 1 if repeated time values were found for the j-th receiver 
    end
 
end

 [~,idx_nonzero] = find(doubled_time==1); % find which receivers have doubled value 

clear time_PANAM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% screen info about observers founded in the file

fprintf( '\n*--------------------------------------------------------------------------*' );
fprintf( '\nLog from <PANAM_SQAT_data_conversion> function' );
fprintf( '\n\n- Input file path: %s' ,char(PATH));

if isempty(idx_nonzero)
    fprintf( '\n\n- Repeated values on the time vector were not found for any receiver.' );
else
    fprintf( '\n\n- WARNING: Repeated values on the time vector were found in receiver no. %i ', idx_nonzero);
end

fprintf( '\n\n- Number of receiver(s) found: %i \n\n',n_observer );

for i=1:n_observer
    
    fprintf('Receiver position %i: x = %.5g\t| y = %.5g \t| z = %.5g\n',i, data{1,i}.xobs,data{1,i}.yobs,data{1,i}.zobs );
    
end
fprintf('*--------------------------------------------------------------------------*\n' );

% clear i myFolders file show

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT: OASPL struct with observers in the columns 
%          OASPL - time-dependent metrics (freq is integrated for each time-step)

OASPL = cell(1,n_observer);      % declaring cell (nTimeSteps,nObserver)

OASPL_dBA = cell(1,n_observer);      % declaring cell (nTimeSteps,nObserver)

for j = 1:size(data,2)   % observer loop
    for i = 1:size(data,1)   % time loop
        
        OASPL{1,j}.airframe(i)=data{i,j}.OASPL_airframe;
        OASPL{1,j}.airframe_toc(i)=data{i,j}.OASPL_airframe_toc;
        OASPL{1,j}.fan_harmonics(i)=data{i,j}.OASPL_fan_harmonics;
        OASPL{1,j}.fan_harmonics_toc(i)=data{i,j}.OASPL_fan_harmonics_toc;
        OASPL{1,j}.engine_broadband(i)=data{i,j}.OASPL_engine_broadband;
        OASPL{1,j}.engine_broadband_toc(i)=data{i,j}.OASPL_engine_broadband_toc;
        OASPL{1,j}.engine_buzzsaw(i)=data{i,j}.OASPL_engine_buzzsaw;
        OASPL{1,j}.engine_buzzsaw_toc(i)=data{i,j}.OASPL_engine_buzzsaw_toc;
        OASPL{1,j}.engine_without_fan_harmonics(i)=data{i,j}.OASPL_engine_without_fan_harmonics;
        OASPL{1,j}.engine(i)=data{i,j}.OASPL_engine;
        OASPL{1,j}.overall_broadband(i)=data{i,j}.OASPL_overall_broadband;
        OASPL{1,j}.overall(i)=data{i,j}.OASPL_overall;
        
        OASPL{1,j}.source_time(i)=data{i,j}.source_time;
        OASPL{1,j}.retarded_time(i)=data{i,j}.retarded_time;
        
        OASPL_dBA{1,j}.airframe(i)=data{i,j}.OASPL_airframe_dBA;
        OASPL_dBA{1,j}.fan_harmonics(i)=data{i,j}.OASPL_fan_harmonics_dBA;
        OASPL_dBA{1,j}.engine_broadband(i)=data{i,j}.OASPL_engine_broadband_dBA;
        OASPL_dBA{1,j}.engine_buzzsaw(i)=data{i,j}.OASPL_engine_buzzsaw_dBA;
        OASPL_dBA{1,j}.engine_without_fan_harmonics(i)=data{i,j}.OASPL_engine_without_fan_harmonics_dBA;
        OASPL_dBA{1,j}.engine(i)=data{i,j}.OASPL_engine_dBA;
        OASPL_dBA{1,j}.overall_broadband(i)=data{i,j}.OASPL_overall_broadband_dBA;
        OASPL_dBA{1,j}.overall(i)=data{i,j}.OASPL_overall_dBA;
        
        OASPL_dBA{1,j}.source_time(i)=data{i,j}.source_time;
        OASPL_dBA{1,j}.retarded_time(i)=data{i,j}.retarded_time;
    end
end

aux = 1;

for i = 1:length(OASPL{1,j}.retarded_time) - 1
    dt(aux) = OASPL{1,j}.retarded_time(i) - OASPL{1,j}.retarded_time(i+1);
    aux = aux + 1;
end

dt = abs(dt);
dt = [dt dt(end)];
    
% SEL - LEq - Lmax
for j = 1:size(data,2)   % observer loop

% SEL - deprecated
%     OASPL{1,j}.SEL.airframe=10.*log10( sum(10.^(OASPL{1,j}.airframe./10)) );
%     OASPL{1,j}.SEL.engine_broadband=10.*log10( sum(10.^(OASPL{1,j}.engine_broadband./10)) );
%     OASPL{1,j}.SEL.engine_buzzsaw=10.*log10( sum(10.^(OASPL{1,j}.engine_buzzsaw./10)) );
%     OASPL{1,j}.SEL.engine_without_fan_harmonics=10.*log10( sum(10.^(OASPL{1,j}.engine_without_fan_harmonics./10)) );
%     OASPL{1,j}.SEL.engine=10.*log10( sum(10.^(OASPL{1,j}.engine./10)) );   
%     OASPL{1,j}.SEL.overall_broadband=10.*log10( sum(10.^(OASPL{1,j}.overall_broadband./10)) );        
%     OASPL{1,j}.SEL.overall=10.*log10( sum(10.^(OASPL{1,j}.overall./10)) ); 
%     
%     OASPL_dBA{1,j}.SEL_dBA.airframe=10.*log10( sum(10.^(OASPL_dBA{1,j}.airframe./10)) );
%     OASPL_dBA{1,j}.SEL_dBA.engine_broadband=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_broadband./10)) );
%     OASPL_dBA{1,j}.SEL_dBA.engine_buzzsaw=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_buzzsaw./10)) );
%     OASPL_dBA{1,j}.SEL_dBA.engine_without_fan_harmonics=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_without_fan_harmonics./10)) );
%     OASPL_dBA{1,j}.SEL_dBA.engine=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine./10)) );   
%     OASPL_dBA{1,j}.SEL_dBA.overall_broadband=10.*log10( sum(10.^(OASPL_dBA{1,j}.overall_broadband./10)) );        
%     OASPL_dBA{1,j}.SEL_dBA.overall=10.*log10( sum(10.^(OASPL_dBA{1,j}.overall./10)) ); 
    
   % SEL - match PANAM
    OASPL{1,j}.SEL.airframe=10.*log10( sum(10.^(OASPL{1,j}.airframe./10).*dt) );
    OASPL{1,j}.SEL.airframe_toc=10.*log10( sum(10.^(OASPL{1,j}.airframe_toc./10).*dt) );
    OASPL{1,j}.SEL.fan_harmonics=10.*log10( sum(10.^(OASPL{1,j}.fan_harmonics./10).*dt) );
    OASPL{1,j}.SEL.fan_harmonics_toc=10.*log10( sum(10.^(OASPL{1,j}.fan_harmonics_toc./10).*dt) );
    OASPL{1,j}.SEL.engine_broadband=10.*log10( sum(10.^(OASPL{1,j}.engine_broadband./10).*dt) );
    OASPL{1,j}.SEL.engine_broadband_toc=10.*log10( sum(10.^(OASPL{1,j}.engine_broadband_toc./10).*dt) );
    OASPL{1,j}.SEL.engine_buzzsaw=10.*log10( sum(10.^(OASPL{1,j}.engine_buzzsaw./10).*dt) );
    OASPL{1,j}.SEL.engine_buzzsaw_toc=10.*log10( sum(10.^(OASPL{1,j}.engine_buzzsaw_toc./10).*dt) );
    OASPL{1,j}.SEL.engine_without_fan_harmonics=10.*log10( sum(10.^(OASPL{1,j}.engine_without_fan_harmonics./10).*dt) );
    OASPL{1,j}.SEL.engine=10.*log10( sum(10.^(OASPL{1,j}.engine./10).*dt) );   
    OASPL{1,j}.SEL.overall_broadband=10.*log10( sum(10.^(OASPL{1,j}.overall_broadband./10).*dt) );        
    OASPL{1,j}.SEL.overall=10.*log10( sum(10.^(OASPL{1,j}.overall./10).*dt) ); 
    
    OASPL_dBA{1,j}.SEL_dBA.airframe=10.*log10( sum(10.^(OASPL_dBA{1,j}.airframe./10).*dt) );
    OASPL_dBA{1,j}.SEL_dBA.fan_harmonics_toc=10.*log10( sum(10.^(OASPL_dBA{1,j}.fan_harmonics./10).*dt) );
    OASPL_dBA{1,j}.SEL_dBA.engine_broadband=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_broadband./10).*dt) );
    OASPL_dBA{1,j}.SEL_dBA.engine_buzzsaw=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_buzzsaw./10).*dt) );
    OASPL_dBA{1,j}.SEL_dBA.engine_without_fan_harmonics=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_without_fan_harmonics./10).*dt) );
    OASPL_dBA{1,j}.SEL_dBA.engine=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine./10).*dt) );   
    OASPL_dBA{1,j}.SEL_dBA.overall_broadband=10.*log10( sum(10.^(OASPL_dBA{1,j}.overall_broadband./10).*dt) );        
    OASPL_dBA{1,j}.SEL_dBA.overall=10.*log10( sum(10.^(OASPL_dBA{1,j}.overall./10).*dt) ); 
    
% Leq - deprecated
%     OASPL{1,j}.Leq.airframe=10.*log10(1/length(OASPL{1,j}.retarded_time)*sum(10.^(OASPL{1,j}.airframe./10)) );
%     OASPL{1,j}.Leq.engine_broadband=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL{1,j}.engine_broadband./10)) );
%     OASPL{1,j}.Leq.engine_buzzsaw=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL{1,j}.engine_buzzsaw./10)) );
%     OASPL{1,j}.Leq.engine_without_fan_harmonics=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL{1,j}.engine_without_fan_harmonics./10)) );
%     OASPL{1,j}.Leq.engine=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL{1,j}.engine./10)) );   
%     OASPL{1,j}.Leq.overall_broadband=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL{1,j}.overall_broadband./10)) );        
%     OASPL{1,j}.Leq.overall=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL{1,j}.overall./10)) ); 
%     
%     OASPL_dBA{1,j}.Leq_dBA.airframe=10.*log10(1/length(OASPL{1,j}.retarded_time)*sum(10.^(OASPL_dBA{1,j}.airframe./10)) );
%     OASPL_dBA{1,j}.Leq_dBA.engine_broadband=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL_dBA{1,j}.engine_broadband./10)) );
%     OASPL_dBA{1,j}.Leq_dBA.engine_buzzsaw=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL_dBA{1,j}.engine_buzzsaw./10)) );
%     OASPL_dBA{1,j}.Leq_dBA.engine_without_fan_harmonics=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL_dBA{1,j}.engine_without_fan_harmonics./10)) );
%     OASPL_dBA{1,j}.Leq_dBA.engine=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL_dBA{1,j}.engine./10)) );   
%     OASPL_dBA{1,j}.Leq_dBA.overall_broadband=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL_dBA{1,j}.overall_broadband./10)) );        
%     OASPL_dBA{1,j}.Leq_dBA.overall=10.*log10(1/length(OASPL{1,j}.retarded_time)* sum(10.^(OASPL_dBA{1,j}.overall./10)) ); 
    
   % Leq - match PANAM
    OASPL{1,j}.Leq.airframe=10.*log10( sum(10.^(OASPL{1,j}.airframe./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.airframe_toc=10.*log10( sum(10.^(OASPL{1,j}.airframe_toc./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.fan_harmonics=10.*log10( sum(10.^(OASPL{1,j}.fan_harmonics./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.fan_harmonics_toc=10.*log10( sum(10.^(OASPL{1,j}.fan_harmonics_toc./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.engine_broadband=10.*log10( sum(10.^(OASPL{1,j}.engine_broadband./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.engine_broadband_toc=10.*log10( sum(10.^(OASPL{1,j}.engine_broadband_toc./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.engine_buzzsaw=10.*log10( sum(10.^(OASPL{1,j}.engine_buzzsaw./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.engine_buzzsaw_toc=10.*log10( sum(10.^(OASPL{1,j}.engine_buzzsaw_toc./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.engine_without_fan_harmonics=10.*log10( sum(10.^(OASPL{1,j}.engine_without_fan_harmonics./10).*dt)/sum(dt) );
    OASPL{1,j}.Leq.engine=10.*log10( sum(10.^(OASPL{1,j}.engine./10).*dt)/sum(dt) );   
    OASPL{1,j}.Leq.overall_broadband=10.*log10( sum(10.^(OASPL{1,j}.overall_broadband./10).*dt)/sum(dt) );        
    OASPL{1,j}.Leq.overall=10.*log10( sum(10.^(OASPL{1,j}.overall./10).*dt)/sum(dt) ); 
    
    OASPL_dBA{1,j}.Leq_dBA.airframe=10.*log10( sum(10.^(OASPL_dBA{1,j}.airframe./10).*dt)/sum(dt) );
    OASPL_dBA{1,j}.Leq_dBA.fan_harmonics=10.*log10( sum(10.^(OASPL_dBA{1,j}.fan_harmonics./10).*dt)/sum(dt) );
    OASPL_dBA{1,j}.Leq_dBA.engine_broadband=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_broadband./10).*dt)/sum(dt) );
    OASPL_dBA{1,j}.Leq_dBA.engine_buzzsaw=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_buzzsaw./10).*dt)/sum(dt) );
    OASPL_dBA{1,j}.Leq_dBA.engine_without_fan_harmonics=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine_without_fan_harmonics./10).*dt)/sum(dt) );
    OASPL_dBA{1,j}.Leq_dBA.engine=10.*log10( sum(10.^(OASPL_dBA{1,j}.engine./10).*dt) );   
    OASPL_dBA{1,j}.Leq_dBA.overall_broadband=10.*log10( sum(10.^(OASPL_dBA{1,j}.overall_broadband./10).*dt)/sum(dt) );        
    OASPL_dBA{1,j}.Leq_dBA.overall=10.*log10( sum(10.^(OASPL_dBA{1,j}.overall./10).*dt)/sum(dt) ); 
    
    % Lmax
    OASPL{1,j}.Lmax.airframe=max(OASPL{1,j}.airframe);
    OASPL{1,j}.Lmax.airframe_toc=max(OASPL{1,j}.airframe_toc);
    OASPL{1,j}.Lmax.fan_harmonics=max(OASPL{1,j}.fan_harmonics);
    OASPL{1,j}.Lmax.fan_harmonics_toc=max(OASPL{1,j}.fan_harmonics_toc);
    OASPL{1,j}.Lmax.engine_broadband=max(OASPL{1,j}.engine_broadband);
    OASPL{1,j}.Lmax.engine_broadband_toc=max(OASPL{1,j}.engine_broadband_toc);
    OASPL{1,j}.Lmax.engine_buzzsaw=max(OASPL{1,j}.engine_buzzsaw);
    OASPL{1,j}.Lmax.engine_buzzsaw_toc=max(OASPL{1,j}.engine_buzzsaw_toc);
    OASPL{1,j}.Lmax.engine_without_fan_harmonics=max(OASPL{1,j}.engine_without_fan_harmonics);
    OASPL{1,j}.Lmax.engine=max(OASPL{1,j}.engine);
    OASPL{1,j}.Lmax.overall_broadband=max(OASPL{1,j}.overall_broadband);
    OASPL{1,j}.Lmax.overall=max(OASPL{1,j}.overall);
    
    OASPL_dBA{1,j}.LAmax.airframe=max(OASPL_dBA{1,j}.airframe);
    OASPL_dBA{1,j}.LAmax.fan_harmonics=max(OASPL_dBA{1,j}.fan_harmonics);
    OASPL_dBA{1,j}.LAmax.engine_broadband=max(OASPL_dBA{1,j}.engine_broadband);
    OASPL_dBA{1,j}.LAmax.engine_buzzsaw=max(OASPL_dBA{1,j}.engine_buzzsaw);
    OASPL_dBA{1,j}.LAmax.engine_without_fan_harmonics=max(OASPL_dBA{1,j}.engine_without_fan_harmonics);
    OASPL_dBA{1,j}.LAmax.engine=max(OASPL_dBA{1,j}.engine);
    OASPL_dBA{1,j}.LAmax.overall_broadband=max(OASPL_dBA{1,j}.overall_broadband);
    OASPL_dBA{1,j}.LAmax.overall=max(OASPL_dBA{1,j}.overall);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT EPNL (loaded from PANAM output folder)

% EPNL=cell(1,n_observer);      % declaring cell (nTimeSteps,nObserver)
% 
% [HEADER_PNLT,DATA_PNLT] = mhdrload([PATH 'level_time_history.dat']);
% a = importdata([PATH 'EPNL.plt']);
% 
% for j=1:size(data,2)   % observer loop
% 
%     EPNL{1,j}.PNLT=DATA_PNLT(:,13,j);
%     EPNL{1,j}.time=DATA_PNLT(:,1,j);
%     
%     EPNL{1,j}.EPNL=a.data(j,4);
%     EPNL{1,j}.TCPNL=a.data(j,5);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT: SPECTROGRAM struct with observers in the columns 
%          struct contains all processed data per time step in   
%          the format - data{nTimeSteps,nObserver} 

SPECTROGRAM = cell(1,n_observer);      % declaring cell (nTimeSteps,nObserver)
SPECTROGRAM_dBA = cell(1,n_observer);      % declaring cell (nTimeSteps,nObserver)

for j = 1:size(data,2)   % observer loop
    
    SPECTROGRAM{1,j}.freq=data{1,1}.airframe_toc(:,1);  % toc bands
    SPECTROGRAM_dBA{1,j}.freq=data{1,1}.airframe_dBA(:,1);  % toc bands
    
    for i=1:size(data,1)   % time loop
        
        SPECTROGRAM{1,j}.source_time(i)=data{i,j}.source_time;
        SPECTROGRAM{1,j}.retarded_time(i)=data{i,j}.retarded_time;
        
        SPECTROGRAM_dBA{1,j}.source_time(i)=data{i,j}.source_time;
        SPECTROGRAM_dBA{1,j}.retarded_time(i)=data{i,j}.retarded_time;
        
        %     SPECTROGRAM{1,j}.fan_harmonics(k,i)=data{i,j}.fan_harmonics(k,2);
        
        for k=1:length(SPECTROGRAM{1,1}.freq)   % freq loop
            
            %dB
            SPECTROGRAM{1,j}.overall(k,i)=data{i,j}.overall(k,2);
            
            SPECTROGRAM{1,j}.overall_broadband(k,i)=data{i,j}.overall_broadband(k,2);
            
            SPECTROGRAM{1,j}.airframe_toc(k,i)=data{i,j}.airframe_toc(k,2);

            SPECTROGRAM{1,j}.engine(k,i)=data{i,j}.engine(k,2);

            SPECTROGRAM{1,j}.engine_without_fan_harmonics(k,i)=data{i,j}.engine_without_fan_harmonics(k,2);

            SPECTROGRAM{1,j}.engine_broadband_toc(k,i)=data{i,j}.engine_broadband_toc(k,2);
            SPECTROGRAM{1,j}.engine_buzzsaw_toc(k,i)=data{i,j}.engine_buzzsaw_toc(k,2);
            SPECTROGRAM{1,j}.fan_harmonics_toc(k,i)=data{i,j}.fan_harmonics_toc(k,2);

            %dBA
            SPECTROGRAM_dBA{1,j}.overall(k,i)=data{i,j}.overall_dBA(k,2);
            
            SPECTROGRAM_dBA{1,j}.overall_broadband(k,i)=data{i,j}.overall_broadband_dBA(k,2);
            
            SPECTROGRAM_dBA{1,j}.airframe(k,i)=data{i,j}.airframe_dBA(k,2);
            SPECTROGRAM_dBA{1,j}.engine(k,i)=data{i,j}.engine_dBA(k,2);
            
            SPECTROGRAM_dBA{1,j}.engine_without_fan_harmonics(k,i)=data{i,j}.engine_without_fan_harmonics_dBA(k,2);
            SPECTROGRAM_dBA{1,j}.engine_broadband(k,i)=data{i,j}.engine_broadband_dBA(k,2);
            SPECTROGRAM_dBA{1,j}.engine_buzzsaw(k,i)=data{i,j}.engine_buzzsaw_dBA(k,2);
            SPECTROGRAM_dBA{1,j}.fan_harmonics(k,i)=data{i,j}.fan_harmonics_dBA(k,2);
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save SQ_data data data_observer

% save emission_data data_observer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% embedded functions for header data extraction

    function [data] = load_header_first(a)   % first header is special (has 2 extra rowa for TITLE and Variables)
        
        data.retarded_time = str2double(a(4,18:30));       % retarded time
        data.source_time = str2double(a(5,16:30));         % source time
        data.sound_speed = str2double(a(6,18:30));         % sound speed
        data.xpos = str2double(a(7,9:30));                 % x pos
        data.ypos = str2double(a(8,9:30));                 % y pos
        data.alt = str2double(a(9,8:30));                  % altitude
        data.xobs = str2double(a(10,10:30));               % x observador
        data.yobs = str2double(a(11,10:30));               % y observador
        data.zobs = str2double(a(12,10:30));               % z observador
        data.phix = str2double(a(13,10:30));               % phi x
        data.phiy = str2double(a(14,10:30));               % phi y
        data.dist = str2double(a(15,10:30));               % distance
        data.distxy = str2double(a(16,10:30));             % distance xy
        data.vel = str2double(a(17,10:30));                % velocity
        data.n1 = str2double(a(18,10:30));                 % N1
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% embedded functions for header data extraction

    function [data] = load_header(a)
        
        data.retarded_time = str2double(a(2,18:30));       % retarded time
        data.source_time = str2double(a(3,16:30));         % source time
        data.sound_speed = str2double(a(4,18:30));         % sound speed
        data.xpos = str2double(a(5,9:30));                 % x pos
        data.ypos = str2double(a(6,9:30));                 % y pos
        data.alt = str2double(a(7,8:30));                  % altitude
        data.xobs = str2double(a(8,10:30));                % x observador
        data.yobs = str2double(a(9,10:30));                % y observador
        data.zobs = str2double(a(10,10:30));               % z observador
        data.phix = str2double(a(11,10:30));               % phi x
        data.phiy = str2double(a(12,10:30));               % phi y
        data.dist = str2double(a(13,10:30));               % distance
        data.distxy = str2double(a(14,10:30));             % distance xy
        data.vel = str2double(a(15,10:30));                % velocity
        data.n1 = str2double(a(16,10:30));                 % N1
        
    end

end