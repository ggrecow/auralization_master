function BroadbandSignal = broadbandSynthesis_smooth(input, source, time_PANAM, time, smoothing, show, tag_auralization, input_type)
%function BroadbandSignal = broadbandSynthesis_smooth(input, source, time_PANAM, time, smoothing, show, tag_auralization, input_type)
%
% INPUTS:
%
%   input : struct
%   struct containing all processed data per time step in the format - data{nTimeSteps,nObserver}. 
%   The broadband nosie description is taken automatically according to the chosen sound 
%    source to be auralized, which is indicated by the <source> string (see below).  
%
%   source : string
%   string containing the BROADBAND noise source to be auralized
%   options are: 'engine', 'airframe', 'overall'
%
%   time_PANAM : array [1xN]
%   time vector from PANAM
%
%   time : array [1xN]
%   time vector of the signal to be synthesised
%
%   smoothing : scalar
%   number of spectral smoothing iterations during the conversion from toc to
%   narrowband (can be ZERO)
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable)
%
%   tag_auralization : string
%   containg the main path (where files should be saved) and their corresponding basic info 
%   (i.e. emission/immission & aircraft & procedure & Receiver #) to be
%   included in the name of the saved files related to auralization processes
%
%   input_type : string
%   tells the function whether the input being auralized is from emission or immission. It only
%   affects plots and info display though.
%   possible strings are:  'emission' or 'immission'
%
% OUTPUT
%
%   BroadbandSignal : row vector [Nx1]
%   synthesised signal, in Pascal values over time, which is only composed by broadband noise
%
% This code is largely based on the work of:
% M. P. Allen, Analysis and synthesis of aircraft engine fan noise for use 
% in psychoacoustic studies. Master thesis, Virginia Polytechnic Institute 
% and State University, 2012
%
% Author: Gil Felix Greco, Braunschweig 30.06.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global pref;     % specified in SQAT_auralization_master()
global fs
global dt_panam
global save_mat_fig

WinFlag = 1;
WinCorr = 1;
OvlpCorr = 1;

%% Signal temporal characteristics (from PANAM)

% dt_panam =  time_PANAM(2) - time_PANAM(1); % time resolution from panam. Assumes it is constant (source noise prediction dt from PANAM is usally .5 seconds and cte) 

%% Signal temporal characteristics (for auralization)

BlockLen = round(fs*dt_panam);               % define synthesis block length [samples]
HopSizeFactor = 5;
HopSize = BlockLen/HopSizeFactor;         % synthesis hop size [samples]
nBlocks = floor(fs*time(end)/HopSize);      % number of synthesis blocks (round down)
% nBlocks = size(input,1); % number of time blocks from PANAM (should be the same as above) 

BlockTime = (0:1:nBlocks-1)'*HopSize/fs;    % time at each block [seconds]
L = (HopSize*nBlocks + BlockLen - HopSize);   % number of output time vector [samples]
OutTime = (0:1:L-1)'/fs;                    % time at each sample [seconds]

df = fs/BlockLen; % freq resolution based on PANAM time-steps (blocks), the BBN auralization is fixed to the same df
% df = 1/dt_panam; % freq resolution based on PANAM time-steps (blocks), the BBN auralization is fixed to the same df

%% divide third octave SPL into narrowband bins

tic;
plot_figs_conversion = 1;
% smoothings = 1; % number of smoothing iterationr during the conversion from toc to narrowband

for i = 1:size(input, 1) % get freq vector from panam
    switch source  
        case 'engine'
            freq_PANAM(:,i) = input{i,1}.engine_broadband(:,1);  
        case 'airframe'
            freq_PANAM(:,i) = input{i,1}.airframe(:,1); 
        case 'overall'
            freq_PANAM(:,i) = input{i,1}.overall_broadband(:,1); 
    end
end

nBands = size(freq_PANAM,1);

BandsVsTime = zeros(size(input,1),nBands);

for i = 1:size(input, 1)
    switch source
        case 'engine'
            BandsVsTime(i,:) = input{i,1}.engine_broadband(:,2); % (nBlocks,nBands), i.e. get spectrogram from PANAM
        case 'airframe'
            BandsVsTime(i,:) = input{i,1}.airframe(:,2); % (nBlocks,nBands), i.e. get spectrogram from PANAM
        case 'overall'
            BandsVsTime(i,:) = input{i,1}.overall_broadband(:,2); % (nBlocks,nBands), i.e. get spectrogram from PANAM  
    end
end

freq_PANAM = freq_PANAM';

if HopSizeFactor==1
    
    [BinsVsTime, ~] = ConvertFreqBands(BandsVsTime,freq_PANAM,df,time_PANAM,BlockLen,smoothing,plot_figs_conversion);
    
else
    
    RepBandsVsTime = zeros(nBlocks, nBands);
    freq_PANAM_rep = zeros(nBlocks, nBands);

    j = 1;
    for i = 1:HopSizeFactor:nBlocks
        
        RepBandsVsTime(i,:) = BandsVsTime(j,:);
        freq_PANAM_rep(i,:) = freq_PANAM(j,:);

        if j ~= nBlocks/HopSizeFactor
            RepBandsVsTime(i+1:i+(HopSizeFactor-1),:) = 0.5.*(BandsVsTime(j,:)+BandsVsTime(j+1,:)).*ones((HopSizeFactor-1),nBands);
            freq_PANAM_rep(i+1:i+(HopSizeFactor-1),:) = freq_PANAM(j,:).*ones((HopSizeFactor-1),nBands);
        else
            RepBandsVsTime(i+1:i+(HopSizeFactor-1),:) = BandsVsTime(j,:).*ones((HopSizeFactor-1),nBands);
            freq_PANAM_rep(i+1:i+(HopSizeFactor-1),:) = freq_PANAM(j,:).*ones((HopSizeFactor-1),nBands);
        end
        
        j = j + 1;
        
        if i == nBlocks
            return
        end
    end
    
    input_time = linspace(0, time_PANAM(end), size(RepBandsVsTime,1));

    [BinsVsTime, ~] = ConvertFreqBands(RepBandsVsTime, freq_PANAM_rep, df, input_time, BlockLen, smoothing, plot_figs_conversion);
end

BinsVsTime = BinsVsTime';

fprintf('\n- BroadbandSynthesis (%s): 1/3-OB to narrowband conversion time:\t%f sec\n', source, toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthesize BBN at each block
tic;

BroadbandSignal = zeros(L,1);        % initialize output [Pa]

% Hann window of appropriate length.  This is a column, so we'll need to
% transpose the rows that come out of the iFFT stage.
win = hann(BlockLen);

% correction factors for overlap and window energy
nOvlp = BlockLen/HopSize;           % number of overlapping pieces
WinEnergy = sum(win.^2/BlockLen);   % window energy [p_rms^2]

for iBlock = 1:nBlocks
    
    % indices of current block
    BlockStartIndex = 1+(iBlock-1)*HopSize;
    BlockEndIndex = BlockStartIndex + BlockLen - 1;
    
    %*****************************************************************
    % no matter how the spectrum at each block was generated, it
    % needs to end up here in the form of pressure amplitude in each bin
    % at each block. Phases get assigned as random no matter what.
    %*****************************************************************
    
    % pressure amplitudes in each bin for the current synthesis block.
    % all work is done as a row vector for now, but later in the loop they
    % get transposed to columns to combine with output.  phases at DC and
    % Nyquist are zero (same thing as zero imaginary component).
    filtAmplitude = BinsVsTime(iBlock,:);
    filtPhase = 2*pi*rand(size(filtAmplitude));
    filtPhase(1) = 0;
    filtPhase(end) = 0;
    
    % generate real and imaginary parts of complex frequency-domain filter
    % **up to the nyquist only**! Remember, we still have to
    % combine/conjugate/flip these to get a usable frequency-domain filter.
    Re_OneSide = filtAmplitude.*cos(filtPhase);
    Im_OneSide = filtAmplitude.*sin(filtPhase);
    
    % form one-sided complex filter from real and imaginary parts
    Filter_OneSide = Re_OneSide + 1i*Im_OneSide;
    
    % scale by block length since output of ifft is amplitude/BlockLen
    Filter_OneSide(2:end) = Filter_OneSide(2:end)*BlockLen;
    
    % scale by two for symmetry around DC (DC doesn't get scaled).  we
    % divide by two since we're about to add the other side of the
    % filter, doubling the amplitudes.  we cut the amplitudes in half now
    % in anticipation.
    Filter_OneSide(2:end) = Filter_OneSide(2:end)/2;
    
    % form the complete complex trace (in frequency domain).
    Filter_Total = [Filter_OneSide ...
        fliplr(conj(Filter_OneSide(2:end-1)))];

    % inverse transform filter to get time-domain block (unwindowed!).
    % the real() command removes tiny imaginary parts.  making sure we
    % constructed a symmetric complex trace earlier ensures that they're
    % on the order of 10^(-14).
    BlockOut = real(ifft(Filter_Total));

    % transpose output block to a column, and window
    if WinFlag == 1
        Grain = BlockOut'.*win;
    else
        Grain = BlockOut';
    end
    
    if WinCorr == 1     % correct for window energy loss
        Grain = Grain./sqrt(WinEnergy);
    end
    
    if OvlpCorr == 1    % correct for overlap
        % this correction corresponds to removing 3dB for each doubling of
        % the number of overlapping windows.
        % 2 windows = 3dB  -->  10^(3/20) = 1.41 times the amplitude
        % 4 windows = 6dB  -->  10^(6/20) = 2    times the amplitude
        % 8 windows = 9dB  -->  10^(9/20) = 2.82 times the amplitude
        % etc...
        Grain = Grain./sqrt(nOvlp);
    end
    
    % add grain to total synthesized signal
    BroadbandSignal(BlockStartIndex:BlockEndIndex) = ...
        BroadbandSignal(BlockStartIndex:BlockEndIndex) + Grain;
    
end

fprintf('- BroadbandSynthesis (%s): Synthesis time:\t%f sec\n', source, toc);
% fprintf('*--------------------------------------------------------------------------*\n');

%%
% due to the block length size, it may occur that the output time L of the 
% auralized BBN does not match the time vector used to auralize the tonal part.
% This is only acceptable if they are off by few samples. here we correct that offset
% by adding (duplicating) or subtracting the last sample(s) of the auralized signal 
if L < length(time) % auralized time is smaller than required time vector
    
    add = BroadbandSignal(end).*ones(length(time)-L,1); % get last value of the auralized BBN signal and generate a vector with the number of samples that we need to add
    BroadbandSignal = [BroadbandSignal;add];
    
elseif L > length(time) % auralized time is larger than required time vector
    
    BroadbandSignal(L-(L-length(time)-1):L)=[]; % remove (L-length(time)) samples from the BroadbandSignal vector      
    
end

%%%% check plot (spectrogram) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show == true    
    tag_title =  ['OUTPUT - Spectrogram of auralized BBN - ' input_type ' - ' source];
    tag_save = [ '_broadbandSignal_spectrogram_' source];
    PLOT_spectrogram(BroadbandSignal, fs, tag_title, tag_auralization, tag_save);    
else
end
% sound(BroadbandSignal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function p = ConvertFreqBands(pMspl,fM,df,t)

    function [pp, f] = ConvertFreqBands(pMspl, fM, df, t, BlockLen, smoothing, plot)
        %
        % INPUTS:
        %   pMspl - spectrogram
        %   fM - original freq vector from panam
        %   df - freq step
        %   t - original time vector from panam
        %
        %   OUTPUT
        %   fMs - narrowband spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         smoothings = 5; % number of smoothing iterations
                
        pM = pref^2*10.^(pMspl'/10); % convert from spl to pascal - % p_rms^2 vs time in current band
        % clear pMspl;
        
        % % f = 0:df:round(fL(end));  % frequency array bounded to max freq predicted by panam and to the giving df
        f = fs/2*linspace(0,1,BlockLen/2+1);
        % f = f1(1):df:f2(end);  % frequency array bounded to max freq predicted by panam and to the giving df

        p = zeros(length(f), length(t));
        pEqual = zeros(length(f), length(t));

        b = 3;  % bands (1: octave, 3: third octave, ...)
        G   = 10^(3/10);
             
        fM = fM';

        for tt = 1:length(t)

            fL(tt,:) = [ (fM(:,tt)./G.^(1./(2.*b)))', fM(end,tt)*G.^(1./(2.*b))]; % Define frequency bands
            n = diff(fL(tt,:))/df;  % number of bins in a band

            % Equal distribution
            for i = 1:size(fM,1)
                idx = f >= fL(tt,i) & f <= fL(tt,i+1);
                p(idx, tt) = sqrt( pM(i, tt) / n(i) )*sqrt(2);
                % p(f > fL(i), tt) = sqrt( pM(i, tt) / n(i) )*sqrt(2);
            end           

            pEqual(:, tt) = p(:, tt);
            
            % Iterative band smoothing with equal distributions
            fMs = 1.0 * fM(:,tt);
            fMtmp = zeros(1, 4*length(fMs));
            
            for j = 1:smoothing
                
                % Iterate over frequency bands
                for i = 1:length(fMs)-1
                    indices = (f > fMs(i)*G.^(1./(4.^j.*b))) & (f < fMs(i+1)/G.^(1./(4.^j.*b)));
                    p(indices, tt) = mean(p(indices, tt));
                end
                
                % Update frequency bands with memory allocation
                k = 1;
                for i = 1:length(fMs)
                    if fMs(i) > fL(tt,5*j)
                        fMtmp(k:k+3) = [fMs(i)./G.^(3./(2.*4.^j.*b)), fMs(i)./G.^(1./(2.*4.^j.*b)), fMs(i).*G.^(1./(2.*4.^j.*b)), fMs(i).*G.^(3./(2.*4.^j.*b))];
                        k = k + 4;
                    end
                end
                fMs = fMtmp(1:k-1);
            end
        end
        
        ff = fs/2*linspace(0,1,BlockLen/2+1);    % frequency bins for auralization need to go till fs/2
        pp = [p; zeros( length(ff)-length(f),length(t)) ]; % correct p to proportional size of frquency bins up to fs/2 by completing the matrix with zeros 

                % sum_pEqual_squared = sum(sum(pEqual.^2));
                % sum_p_squared = sum(sum(p.^2));
                % 
                % disp(sum_pEqual_squared);
                % disp(sum_p_squared);
        
        %% plot spl vs freq (1/3-OB and NB for a time-step)
        
        if plot==1
            
            figure('name','PROCESSING - BBN auralization: TOC to Narrowband conversion (1st block)');
            h  = gcf;
            set(h,'Units','Inches');
            pos = get(h,'Position');
            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

            % semilogx(fM,pMspl(1, :), 'k'); hold all;  % original input 1/3-OB

            % idx_t = 20;
             idx_t = round(size(pMspl,1)/2);
            % idx_t = 141; % overhead for Case 6, centerline receiver

            for i=1:size(fL,2)-1
                a = semilogx( [ fL(idx_t,i) fL(idx_t,i+1) ], [pMspl(idx_t,i) pMspl(idx_t,i)], 'Color', [ 0, 0, 0]);hold on; % create horizontal lines
                semilogx( linspace( fL(idx_t, i), fL(idx_t, i), 24), linspace(0,pMspl(idx_t,i) ,24), 'Color', [ 0, 0, 0]); % create vertical lines (flow)
                % semilogx(linspace(TOC(i,2),TOC(i,2),24),linspace(0,BandsVsTime(1,i),24),'b:'); % create vertical lines (fcenter)
                semilogx( linspace( fL(idx_t, i+1), fL(idx_t, i+1), 24 ), linspace(0,pMspl(idx_t,i),24), 'Color', [ 0, 0, 0]); % create vertical lines (fupper)
            end

            b = semilogx(f, p2spl(pEqual(:, 1)), 'b'); hold all;  % Plotting equal distribution
            % c = semilogx(f, p2spl(p(:, 1)), 'r'); % iterative band smoothing
            
            % legend([a,b,c], {'input: 1/3-OB',...
            %              sprintf('Equal Distribution: $\\Delta f=%g$ Hz',df),...
            %              sprintf('Iterative Band Smoothing: $\\Delta f=%g$ Hz',df)},...
            %              'Location','SW','Interpreter','Latex');

            legend([a,b], {'input: 1/3-OB',...
                sprintf('Equal Distribution: $\\Delta f=%g$ Hz',df)},...
                'Location','SW','Interpreter','Latex');

            xlabel('Frequency, $f$ (kHz)','Interpreter','Latex');
            ylabel('SPL, $L_{\mathrm{p,Z}}$ (dB re 20$~\mu$Pa)','Interpreter','Latex');
            xlim([20, 15000]);
            grid on;
            set(gcf,'color','w');

            if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
            else
                filename = strcat(tag_auralization, '_BBN_TOC_to_ narrowband_', source);
                save_pdf = 1; save_png = 0;
                export_figures( filename, save_mat_fig, save_png, save_pdf );
            end

            %% plot original and smoothed spectrograms
            
            % figure('name','PROCESSING - BBN auralization: TOC to Narrowband conversion - spectrogram ');
            % h  = gcf;
            % set(h,'Units','Inches');
            % pos = get(h,'Position');
            % set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            % 
            % subplot(1, 3, 1);         % Plotting original spectrogram
            % pcolor(t, fM./1000, pMspl');
            % shading 'flat';
            % colormap jet;
            % title('Original spectrogram');
            % ylabel('Frequency, $f$ (kHz)','Interpreter','Latex');
            % xlabel('Time, $t$ (s)','Interpreter','Latex');
            % ylabel(colorbar, 'SPL, $L_{\mathrm{p,Z}}$ (dB re 20$~\mu$Pa)','Interpreter','Latex');
            % 
            % ylim([0 15]);
            % caxis([0 90]);
            % %         caxis([0 max(max(p2spl(pM)))]);
            % 
            % % Plotting spectrogram with iterative band smoothing
            % subplot(1, 3, 2);
            % pcolor(t, f./1000, p2spl(pEqual));
            % shading 'flat';
            % colormap jet;
            % title('Equal distr. spectrogram');
            % ylabel('Frequency, $f$ (kHz)','Interpreter','Latex');
            % xlabel('Time, $t$ (s)','Interpreter','Latex');
            % ylabel(colorbar, 'SPL, $L_{\mathrm{p,Z}}$ (dB re 20$~\mu$Pa)','Interpreter','Latex');
            % 
            % ylim([0 15]);
            % caxis([0 90]);
            % %         caxis([0 max(max(p2spl(pEqual)))]);
            % 
            % % Assuming tNew and pNew are defined earlier
            % subplot(1, 3, 3);
            % pcolor(t, f./1000, p2spl(p));
            % shading 'flat';
            % colormap jet;
            % title('Iterative Band Smoothing');
            % ylabel('Frequency, $f$ (kHz)','Interpreter','Latex');
            % xlabel('Time, $t$ (s)','Interpreter','Latex');
            % set(gcf,'color','w');
            % 
            % ylim([0 15]);
            % caxis([0 90]);
            % %         caxis([0 max(max(p2spl(p)))]);
            % 
            % ylabel(colorbar, 'SPL, $L_{\mathrm{p,Z}}$ (dB re 20$~\mu$Pa)','Interpreter','Latex');
            % 
            % set(gcf,'color','w');
            % 
            % if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
            % else
            %     filename = strcat(tag_auralization, '_BBN_TOC_to_ narrowband_spectrogram_', source);
            %     save_pdf = 0; save_png = 1;
            %     export_figures( filename, save_mat_fig, save_png, save_pdf );
            % end

        else
        end
        
    end % end function p = ConvertFreqBands(pMspl,fM,df,t)

end