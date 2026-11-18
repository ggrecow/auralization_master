function PLOT_amplitude_fluctuations( AF, t_ctrl, tag_auralization )
% function PLOT_amplitude_fluctuations( AF, t_ctrl, tag_auralization )
%
% Plots the standard deviation of the turbulence-induced level fluctuations
% used in the auralization (semi-empirical model of Lincke & Pieren 2024, see
% <get_log_amplitude_fluctuation>), over time and frequency.
%
% Figure (saved if save_figs = 1):
%   <tag>_amplitude_fluctuations_std.pdf
%     top   : sigma(f, t) in dB, 1/3-octave bands over time
%     bottom: sigma(t) for selected 1/3-octave bands, and the limit frequency
%             f_lim between weak and strong scattering (right axis)
%
% INPUTS
%   AF : struct (from <get_propagation>, OUT.amplitudeFluctuations)
%       AF.f     : [1 x nF] 1/3-octave band centre frequencies [Hz]
%       AF.sigma : [nT x nF] standard deviation of the level fluctuations [dB]
%       AF.flim  : [nT x 1] limit frequency [Hz] (NaN if above 20 kHz)
%   t_ctrl : [nT x 1] time of each row of AF.sigma on the receiver time axis [s]
%   tag_auralization : string, path + name tag used to save the figure
%       (empty: figure is not saved)
%
% Author: Gil Felix Greco, Braunschweig (29.10.2026)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global save_mat_fig

fontSize = 12;

h = figure( 'name', 'PROCESSING - Turbulence-induced amplitude fluctuations' );
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'Position', [pos(1), pos(2), pos(3), pos(4)*1.5]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)*1.5])

% top: sigma over time and frequency
subplot(2,1,1)
imagesc( t_ctrl, 1:numel(AF.f), AF.sigma.' ); set(gca,'YDir','Normal');
cb = colorbar('vert');
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Std. of level fluctuations, $\hat{\sigma}_\chi$ (dB)';
cb.Label.FontSize = fontSize;
fTicks = [63 250 1000 4000 16000];
[~, iTicks] = min( abs( log(AF.f(:)) - log(fTicks) ), [], 1 );
yticks( iTicks ); yticklabels( {'63', '250', '1k', '4k', '16k'} );
xlabel( 'Time, $t$ (s)', 'Interpreter', 'Latex' );
ylabel( 'Frequency, $f$ (Hz)', 'Interpreter', 'Latex' );

% bottom: selected bands and limit frequency
subplot(2,1,2)
fSel = [250 1000 4000];
hold on; box on; grid on
lab = cell(1, numel(fSel));
for k = 1:numel(fSel)
    [~, i] = min( abs( log(AF.f) - log(fSel(k)) ) );
    plot( t_ctrl, AF.sigma(:,i), 'LineWidth', 1.2 );
    lab{k} = sprintf( '%g Hz', fSel(k) );
end
ylabel( '$\hat{\sigma}_\chi$ (dB)', 'Interpreter', 'Latex' );
yyaxis right
plot( t_ctrl, AF.flim, 'k--', 'LineWidth', 1 );
set( gca, 'YScale', 'log' );
ylabel( '$f_{\mathrm{lim}}$ (Hz)', 'Interpreter', 'Latex' );
legend( [lab, {'$f_{\mathrm{lim}}$'}], 'Interpreter', 'Latex', 'Location', 'best' );
xlabel( 'Time, $t$ (s)', 'Interpreter', 'Latex' );

set( gcf, 'color', 'w' );

if ~isempty( tag_auralization )
    filename = strcat( tag_auralization, '_amplitude_fluctuations_std' );
    save_pdf = 1; save_png = 0;
    export_figures( filename, save_mat_fig, save_png, save_pdf );
end

fprintf( '\n- Amplitude fluctuations: standard deviation over time and frequency plotted\n' );

end
