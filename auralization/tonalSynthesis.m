function TonalSignal = tonalSynthesis(input, time_PANAM, time, show, tag_auralization, tag_source, input_type)
% function TonalSignal = tonalSynthesis(input, time_PANAM, time, show, tag_auralization, tag_source, input_type)
%
% INPUTS:
%
%   input : struct 
%   contains two vectors describing the tones' freq and SPL evolution over
%   time/source position
%
%       - <tonesFreqTime> [nFreq x nTime] : which gives the evolution of the tones' freq over
%           time/source position 
%
%       - <tonesSPLTime> [nFreq x nTime] : which gives the evolution of the tones' SPL over
%           time/source position 
%
%   time_PANAM : array [1xN]
%   time vector from PANAM (retarded time, starting in t=0) with a fixed dt_panam
%
%   time : array [1xN]
%   time vector of the signal to be synthesised/auralized
%
%   show : logical (boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable)
%
%   tag_auralization : string
%   containg the main path (where files should be saved) and their corresponding basic info 
%   (i.e. emission/immission & aircraft & procedure & Receiver #) to be
%   included in the name of the saved files related to auralization processes
%
%   tag_source : string
%   tells which tonal noise is being auralized to the function. It only
%   affects plots and info display though.
%   possible strings are:  'fan_harmonics' or 'buzzsaw'
%
%   input_type : string
%   tells the function whether the input being auralized is from emission or immission. It only
%   affects plots and info display though.
%   possible strings are:  'emission' or 'immission'
%
% OUTPUT
%
%   TonalSignal : row vector [Nx1]
%   contains the synthesised signal composed by the sum of the individual
%   tones predicted by panam, in Pascal values over time
%
% This code is largely based on the work of:
% M. P. Allen, Analysis and synthesis of aircraft engine fan noise for use 
% in psychoacoustic studies. Master thesis, Virginia Polytechnic Institute 
% and State University, 2012
% 
% Author: Gil Felix Greco, Braunschweig 20.06.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global pref   % specified in SQAT_auralization_master()
global fs
global save_mat_fig

%% get tonal content over time as from PANAM

% number of tones (assumes they dont change over time)
tones = input.tones;  

% freq of the tones over time
tonesFreqTime = input.tonesFreqTime;  

% SPL of the tones overtime
tonesSPLTime = input.tonesSPLTime;  

%% check plot (input data) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show == true

    figure('name',['INPUT - Tonal content from PANAM - ' input_type ' - ' tag_source]);
    h  = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    for nTones = 1:tones % loop over the number of tones i-th source_time steps

        subplot(2,1,1); plot(time_PANAM,tonesSPLTime(nTones,:)); hold on;
        ylabel('SPL, $L_{\mathrm{p,Z}}$ (dB)','Interpreter','Latex');

    end

    switch tag_source
        case  'fan_harmonics'
            l=legend('Tone 1','Tone 2','Tone 3','Tone 4','Tone 5'); set(l,'Orientation','horizontal','Location','northoutside');
        case  'buzzsaw'
            l=legend( sprintf('Buzzsaw tones - K = %g',tones)); set(l,'Orientation','horizontal','Location','northoutside');
    end


    legend box off;

    for nTones = 1:tones % loop over the number of tones i-th source_time steps

        subplot(2,1,2); plot(time_PANAM,tonesFreqTime(nTones,:)./1000); hold on;
        xlabel('PANAM receiver time (trimmed), $t^*_{\mathrm{P,i}}$ (s)','Interpreter','Latex');
        ylabel('Frequency, $f$ (kHz)','Interpreter','Latex');

    end

    set(gcf,'color','w');

    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization,  '_tonalContent_PANAM_', tag_source);
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

else
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% translate amp and freq of the tones over time from panam to auralization time sampling
% approach 1 (build vectors per blocks)

% dt_panam = 0.5; % time discretization from PANAM
% LengthBlock = ceil(dt_panam/dt); % number of samples I need to get one dt_panam
% nBlocks = floor(nSamples/LengthBlock); % number of blocks i need to build the whole time vector
%
% EnvelopeSPL = zeros(length(tones),nSamples);
% EnvelopeFreq = zeros(length(tones),nSamples);
%
% for nTones = 1:length(tones)
%
%     j = 1;
%     EnvelopeSPL_aux = zeros(LengthBlock,nBlocks); % initialize/clean vector
%     EnvelopeFreq_aux = zeros(LengthBlock,nBlocks); % initialize/clean vector
%
%     for i = 1:nBlocks
%
%         EnvelopeSPL_aux(:,i) = tonesAmpTime(nTones,j)*ones(LengthBlock,1); % get SPL over time in auralization time
%         EnvelopeFreq_aux(:,i) = tonesFreqTime(nTones,j)*ones(LengthBlock,1); % get freq over time in auralization time
%         j=j+1;
%
%     end
%
%     EnvelopeSPL(nTones,:) = reshape(EnvelopeSPL_aux,1,(LengthBlock*nBlocks));
%     EnvelopeFreq(nTones,:) = reshape(EnvelopeFreq_aux,1,(LengthBlock*nBlocks));
% end
%
% EnvelopePressure = sqrt(2).*pref.*10.^(EnvelopeSPL./20); % get sound pressure over time in auralization time
%
% clear EnvelopeSPL_aux EnvelopeFreq_aux
%%%% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on;
% subplot(2,1,1); plot(time,EnvelopeSPL,'o');
% subplot(2,1,2); plot(time,EnvelopeFreq,'o');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% translate amp and freq of the tones over time from panam to auralization time sampling
% approach 2 interpolation (use this)

for nTones = 1:tones % loop over the number of tones i-th source_time steps
    
    % 1) convert SPL to pressure --> 2) interpolate
%     EnvelopeAmplitude(nTones,:) = interp1(time_PANAM,...
%                                           sqrt(2).*pref.*10.^(tonesSPLTime(nTones,:)./20),...
%                                           time,'linear','extrap'); % get sound pressure (amplitude) over time in auralization time from SPL (amplitude) provided by PANAM

    % 1) interpolate --> 2) convert SPL to pressure
    EnvelopeSPL(nTones,:) = interp1( time_PANAM, tonesSPLTime(nTones,:), time, 'linear', 'extrap' ); % get SPL over time in auralization time
%     EnvelopeAmplitude(nTones,:) = sqrt(2).*pref.*10.^(EnvelopeSPL(nTones,:)./20); % get sound pressure (amplitude) over time in auralization time from SPL (amplitude) provided by PANAM
    EnvelopeAmplitude(nTones,:) = sqrt(2)*sqrt(pref^2.*10.^(EnvelopeSPL(nTones,:)/10));
    
    EnvelopeFreq(nTones,:) = interp1(time_PANAM,tonesFreqTime(nTones,:),time,'linear','extrap'); % get freq over time in auralization time
    
end

%%%% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on;
% subplot(2,1,1); plot(time,EnvelopeSPL,'o');
% subplot(2,1,2); plot(time,EnvelopeFreq,'o');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% sinthesize tonal components using additive synthesis
tic;

a=-pi; b=pi; % random number btw [-pi pi], use for a cos function
% a=0; b=2*pi; % random number btw [0 2*pi], use for a sin function

phase = (b-a).*rand(1,tones) + a; % create random phase vector

kf = (2*pi)/fs; % modulation sensitivity for 1 Hz/Hz variation

nSamples = length(time);           % number of samples to synthesise
TonalSignal = zeros(1,nSamples);   % initialize output

for nTones = 1:tones % loop over the number of tones i-th source_time steps
    
    currentTone = EnvelopeAmplitude(nTones,:).*... % amplitude of each tone, p(t) in Pascal
                  cos(kf*cumsum(EnvelopeFreq(nTones,:)) + phase(nTones)); % freq of each tone
    
    TonalSignal = TonalSignal + currentTone; % sum all tones
    
    %%%% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure; plot(time,currentTone);
%         a = yline(rms(currentTone),'k--');
%         legend(a,sprintf('$p_{\\mathrm{rms}}=$%g (Pa)',rms(currentTone)),'Location','NorthEast','Interpreter','Latex'); %legend boxoff
%         xlabel('Time, $t$ (s)','Interpreter','Latex');
%         ylabel('Sound pressure, $p$ (Pa)','Interpreter','Latex'); %grid on;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     sound(currentTone);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear currentTone
end

% convention: sound file always row vector
TonalSignal = transpose(TonalSignal);  

switch tag_source
    case  'fan_harmonics'
        fprintf('\n*--------------------------------------------------------------------------*');
        fprintf( '\nLog from <SQAT_auralization_master> function (%s-based) \n' , input_type);
end

fprintf('\n- TonalSynthesis (%s): Synthesis time:\t%f sec\n', tag_source,toc)

%%%% check plot (spectrogram) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show == true
    
    tag_title = ['OUTPUT - Spectrogram of auralized tonal noise - ' input_type ' - ' tag_source];
    tag_save = [ '_tonalSignal_Spectrogram_' tag_source];
    PLOT_spectrogram(TonalSignal, fs, tag_title, tag_auralization, tag_save);
    
else 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
