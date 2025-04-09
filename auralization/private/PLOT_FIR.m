function PLOT_FIR( impulseResponse_FIR, FRF_ART, tBlock, fs, tag_auralization)
% function PLOT_FIR( impulseResponse_FIR, FRF_ART, tBlock, fs, tag_auralization)
%
% This function plots the freq reponse of the FIR filter and compares to
% the original one. Important to double check how good the FIR filters are.
%
% INPUTS
% <impulseResponse> [nFreq x nTime] to be ploted
% FRF_ART : single sided spectrum of the FRF from ART
% tBlock : nTime, single value only, meaning we plot the FRF of <impulseResponse> only for one fixed nTime value
% fs: sampling freq. Necessary  to create a freq vector
%
% Gil Felix Greco, Braunschweig 12.03.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global save_mat_fig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nTaps

%  compute freq vector for <impulseResponse_FIR>
nTaps = size( impulseResponse_FIR, 1 );
nBins = ceil( (nTaps+1)/2) ;
frequencyVector_FIR = linspace(0, fs/2, nBins);

% get single-sided spectrum of <impulseResponse_FIR>
FR_cut = fft(impulseResponse_FIR);
FR_cut_single_side = FR_cut(1:ceil((size(FR_cut,1)+1)/2), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZERO-PADDED (this is the actuall size of the vector due to circular
% convolution

%  compute freq vector for <impulseResponseCut>
% nfreq = size(transferFunction_doubleSideSpectrum,1);
% nBins = ceil( (nfreq +1)/2) ;
% frequencyVector_cut = linspace(0, fs/2, nBins);
%
% % get single-sided spectrum of <impulseResponseCut>
% FR_cut = fft(impulseResponseCut, nfreq);
% FR_cut_single_side = FR_cut(1:ceil((size(FR_cut,1)+1)/2), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tBlock = round(length(transferFunction)/2);
% tBlock = 21;

%  compute freq vector for <FRF_ART>
frequencyVector_FRF_ART = linspace(0, fs/2, size( FRF_ART, 1 ));

h  = figure;
set(h, 'name', 'PROCESSING - Atmospheric transfer function (FIR filter)');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% original FRF from FIR filter
b = semilogx(frequencyVector_FIR, 10.*log10(abs(FR_cut_single_side(:, tBlock)).^2) );  hold on;

% original FRF from ART
a = semilogx(frequencyVector_FRF_ART, 10.*log10(abs(FRF_ART(:, tBlock)).^2), '--', 'LineWidth', 1.5 );  

xlabel( 'Frequency, $f$ (Hz)', 'Interpreter', 'Latex' );

ylabel('$10\cdot\log_{10}\left(|A_{\mathrm{atm}}|^2 \right)$','Interpreter','Latex');
% ylabel( 'Transmission loss, TL (dB re 20$~\mu$Pa) ', 'Interpreter', 'Latex' );

xlim( [ 20 25e3 ] );
ylim( [ -200 0 ] );

xticks([ 20 1e2 1e3 10e3 20e3 ]);
xticklabels({'20','100','1k', '10k', '20k'});

h = gca;
h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
h.YAxis.MinorTickValues = -200:10:0; % Minor ticks which don't line up with majors

if length(frequencyVector_FRF_ART) == length(frequencyVector_FIR)

yyaxis right
semilogx(frequencyVector_FRF_ART, 10.*log10(abs(FR_cut_single_side(:, tBlock)).^2) - 10.*log10(abs(FRF_ART(:, tBlock)).^2), '-', 'LineWidth', 0.5 );  
ylim( [ -20 20 ] );
ylabel('FIR-ART','Interpreter','Latex');

end

legend([a, b], 'ART', sprintf('FIR - %.5g Taps', nTaps), 'Location', 'southwest');
grid on;

set( gcf, 'color', 'w' );

if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
else
    filename = strcat(tag_auralization, '_atmospheric_transfer_function_FIR');
    save_pdf = 1; save_png = 0;
    export_figures( filename, save_mat_fig, save_png, save_pdf );
end

end