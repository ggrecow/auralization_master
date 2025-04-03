function PLOT_spectrogram(input, fs, tag_title, tag_auralization, tag_save)

global save_mat_fig

%% main custom <PLOT_spectrogram> function 
figure( 'name', tag_title);

windowsize = 1024;
pref = 20e-6;
% pref = 1;

tiledlayout('vertical')

for nChannels = 1:size( input, 2 )

    ax = nexttile(nChannels);

    [P,F,T] = myspecgram(input(:,nChannels), fs, windowsize , 0.75); % overlap is 75% of nfft here

    SPL = 20.*log10(abs(P)./pref);
    imagesc(T,F./1000,SPL) ; 
    set(gca,'YDir','Normal')

    % jet_white = load('utilities\COLORMAP_jet_white');
    % colormap('jet'); % colormap(jet_white.jet_white); % colormap(artemis_colormap);
    colormap(ax, 'jet');

    if max(max(SPL))>110 % probably emission at the source
        cMin  = 100;
    else
        cMin = 0;
    end

    ylim([0 15]); % freq axis
    clim([cMin max(max(SPL))]);
    ylabel('Frequency, $f$ (kHz)','Interpreter','Latex');
    xlabel('Time, $t$ (s)', 'Interpreter', 'Latex');

    if size( input, 2 ) >1
        if nChannels==1
            ax.Title.String = 'Left ear';
            ax.Title.FontWeight = 'normal';
            ax.Title.Interpreter = 'latex';
        else
            ax.Title.String = 'Right ear';
            ax.Title.FontWeight = 'normal';
            ax.Title.Interpreter = 'latex';
        end
    end

end

% Common features
cb = colorbar;
cb.Layout.Tile = 'east';

zString = 'SPL, $L_{\mathrm{Z}}$ (dB re 20$~\mu$Pa)';
set( get(cb,'label'), 'string', zString, 'fontsize', 16, 'Interpreter','Latex');

set(gcf,'color','w');

if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
else
    filename = strcat(tag_auralization, tag_save);
    save_pdf = 0; save_png = 1;
    export_figures( filename, save_mat_fig, save_png, save_pdf );
end

%% myspecgram function

    function  [fft_specgram,freq_vector,time] = myspecgram(signal, Fs, nfft, Overlap)
        %
        % FFT peak spectrogram of signal  (example sinus amplitude 1   = 0 dB after fft).
        %   signal - input signal,
        %   Fs - Sampling frequency (Hz).
        %   nfft - FFT window size
        %   Overlap - buffer overlap % (between 0 and 0.95)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % to check amplitude, generate a sine with know rms amplitude of 40 dB SPL
        % fc=1000;    % Center frequency (Hz)
        % LevelIn=40; % Level (dBSPL)
        % L=10;       % Length (seconds)
        % fs=48000;   % Sampling frequency
        % dt=1/fs;    % Time step
        % t = 0:dt:L; % Time vector
        %
        % %% Make signal
        %
        % A=20e-6*10^(LevelIn/20)*sqrt(2); % Amplitude
        % RefSignal=A*sin(2*pi*fc*t)';     % Generate signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % source: https://de.mathworks.com/matlabcentral/answers/775437-why-does-spectrogram-of-stft-shows-different-magnitude-while-changing-the-window-size
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        signal = signal(:);
        samples = length(signal);
        % fill signal with zeros if its length is lower than nfft
        if samples<nfft
            s_tmp = zeros(nfft,1);
            s_tmp((1:samples)) = signal;
            signal = s_tmp;
            samples = nfft;
        end
        % window : hanning
        window = hanning(nfft);
        window = window(:);
        %    compute fft with overlap
        offset = fix((1-Overlap)*nfft);
        spectnum = 1+ fix((samples-nfft)/offset); % Number of windows
        %     % for info is equivalent to :
        %     noverlap = Overlap*nfft;
        %     spectnum = fix((samples-noverlap)/(nfft-noverlap));	% Number of windows
        % main loop
        fft_specgram = [];
        for ci=1:spectnum
            start = (ci-1)*offset;
            sw = signal((1+start):(start+nfft)).*window;
            fft_specgram = [fft_specgram abs(fft(sw))*4/nfft];     % X=fft(x.*hanning(N))*4/N; % hanning only
        end
        % one sidded fft spectrum  % Select first half
        if rem(nfft,2)    % nfft odd
            select = (1:(nfft+1)/2)';
        else
            select = (1:nfft/2+1)';
        end
        fft_specgram = fft_specgram(select,:);
        freq_vector = (select - 1)*Fs/nfft;
        % time vector
        % time stamps are defined in the middle of the buffer
        time = ((0:spectnum-1)*offset + round(offset/2))/Fs;
    end


end