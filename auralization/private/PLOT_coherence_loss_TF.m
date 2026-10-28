function PLOT_coherence_loss_TF( TF, w, fs, dt_block, idx_closest, tag_auralization )
% function PLOT_coherence_loss_TF( TF, w, fs, dt_block, idx_closest, tag_auralization )
%
% Compares the atmospheric transfer functions (direct + ground-reflected path)
% WITHOUT and WITH turbulence-induced coherence loss. The transfer function
% with coherence loss is the realization actually used in the auralization,
% i.e. the reflected path filtered with the random decorrelation filters of
% <get_FIR> (same <coherence_seed>, redrawn every 0.1 s). Both transfer
% functions are built as H_direct + H_reflected, so they differ ONLY by the
% decorrelation filters.
%
% Two figures are created (and saved, if save_figs = 1):
%
%   1) <tag>_coherence_loss_TF_difference_spectrogram.png
%      spectrogram of the level difference (time x frequency)
%        dL = 10*log10(|H_d + H_r*B|^2) - 10*log10(|H_d + H_r|^2)   (dB)
%      positive: interference dips filled in; negative: interference peaks lowered.
%      The dashed line (see legend) marks the closest point to the receiver.
%      Colorbar limits are adapted to each case (symmetric, 99th percentile
%      of |dL|, see <il_symmetric_clim>)
%
%   2) <tag>_coherence_loss_TF_comparison.pdf
%      transfer functions with and without coherence loss at two time steps:
%      the closest point to the receiver (shortest direct-path propagation
%      time) and the time step with the strongest coherence loss (largest w)
%
% INPUTS
%   TF : struct (from <get_FIR>, OUT.TF_coherenceLoss)
%       TF.coherent     : [nBins x nTimeSteps] single-sided, without coherence loss
%       TF.decorrelated : [nBins x nTimeSteps] single-sided, with coherence loss
%   w : [nTimeSteps x 1] coherence parameter, Ccoh(f) = exp(-w*f^2) [s^2]
%   fs : sampling frequency [Hz]
%   dt_block : duration of each time step (FIR block) [s]
%   idx_closest : index of the time step at the closest point to the receiver
%       (e.g. shortest direct-path propagation time). For a receiver below the
%       flight path, this is the overhead position
%   tag_auralization : string, path + name tag used to save the figures
%       (empty: figures are not saved)
%
% Gil Felix Greco, Braunschweig 28.10.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global save_mat_fig

fontSize = 12;

[nBins, nTimeSteps] = size( TF.coherent );
freq = linspace( 0, fs/2, nBins );        % [Hz]
time = ( 0:nTimeSteps-1 ) .* dt_block;    % [s], start of each time step

% levels (eps avoids -Inf in perfect interference dips)
L_coherent     = 10.*log10( abs( TF.coherent ).^2     + eps );
L_decorrelated = 10.*log10( abs( TF.decorrelated ).^2 + eps );
dL = L_decorrelated - L_coherent;

% coherence factor at 1 kHz, for the figure titles
Ccoh_1kHz = exp( -w(:) .* 1e6 );

%% figure 1: spectrogram of the level difference

h = figure( 'name', 'PROCESSING - Transfer function: with minus without coherence loss' );
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

imagesc( time, freq, dL ); set(gca,'YDir','Normal');
colormap( gca, il_diverging_colormap(256) );
% colorbar limits adapted to each case: symmetric around 0 dB (white = no change),
% set by the 99th percentile of |dL| within the plotted frequency range, rounded up
% to the next full dB (up to 5 dB) or to the next multiple of 5 dB (above 5 dB).
% A percentile instead of the maximum is used, because nearly perfect interference
% dips of the coherent case give isolated, very large values (limited only by <eps>)
% that would wash out the rest of the plot.
f_plot_max = min( 8e3, fs/2 );
idx_f_plot = freq >= 20 & freq <= f_plot_max;
cLim = il_symmetric_clim( dL(idx_f_plot, :), 99 );
clim( [-cLim cLim] );

cb = colorbar('vert');
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\Delta L$ (dB), with minus without coherence loss';
cb.Label.FontSize = fontSize;

hold on;
% dashed line: closest point to the receiver, i.e. time step with the shortest
% direct-path propagation time (same time step as the upper panel of figure 2)
hx = xline_local( time(idx_closest), 'k--', 'LineWidth', 1 );
if ~isempty( hx ) % xline not available in MATLAB <= R2018a
    legend( hx, 'closest point to the receiver', 'Location', 'northeast' );
end

xlabel( 'Time, $t$ (s)', 'Interpreter', 'Latex' );
ylabel( 'Frequency, $f$ (Hz)', 'Interpreter', 'Latex' );
ylim( [ 20 f_plot_max ] );
yticks( [ 20 2e3 4e3 6e3 8e3 ] );
yticklabels( {'20', '2k', '4k', '6k', '8k'} );

set( gcf, 'color', 'w' );

if ~isempty( tag_auralization )
    filename = strcat( tag_auralization, '_coherence_loss_TF_difference_spectrogram' );
    save_pdf = 0; save_png = 1; % color plot
    export_figures( filename, save_mat_fig, save_png, save_pdf );
end

%% figure 2: transfer functions at the closest point and at the strongest coherence loss

[~, idx_max_w] = max( w(:) );

if idx_max_w == idx_closest
    idx_plot = idx_closest;
    tag_plot = {'closest point and strongest coherence loss'};
else
    idx_plot = [idx_closest, idx_max_w];
    tag_plot = {'closest point', 'strongest coherence loss'};
end

h = figure( 'name', 'PROCESSING - Transfer function with and without coherence loss' );
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'Position', [pos(1), pos(2), pos(3), pos(4)*numel(idx_plot)*0.75]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)*numel(idx_plot)*0.75])
movegui(h, 'center');

t = tiledlayout( numel(idx_plot), 1 );
t.TileSpacing = 'compact';

for k = 1:numel(idx_plot)

    ax = nexttile;
    n = idx_plot(k);

    % with coherence loss first (solid), without coherence loss on top (dotted),
    % so that both curves stay visible where they overlap
    b = semilogx( freq, L_decorrelated(:,n), '-', 'Color', [0.80 0.15 0.15], 'LineWidth', 1 ); hold on;
    a = semilogx( freq, L_coherent(:,n), 'k:', 'LineWidth', 1.5 );

    xlim( [ 20 min(20e3, fs/2) ] );
    xticks( [ 20 1e2 1e3 10e3 20e3 ] );
    xticklabels( {'20', '100', '1k', '10k', '20k'} );
    ylabel( '$10\cdot\log_{10}\left(|H_{\mathrm{atm}}|^2 \right)$', 'Interpreter', 'Latex' );
    grid on;

    ax.Title.String = sprintf( '%s: $t$ = %.1f s, $C_{\\mathrm{coh}}$(1 kHz) = %.2f', ...
        tag_plot{k}, time(n), Ccoh_1kHz(n) );
    ax.Title.FontWeight = 'normal';
    ax.Title.Interpreter = 'latex';

    if k == numel(idx_plot)
        xlabel( 'Frequency, $f$ (Hz)', 'Interpreter', 'Latex' );
    else
        xticklabels( '' );
    end

    if k == 1
        legend( [b, a], 'with coherence loss', 'without coherence loss', 'Location', 'southwest' );
    end

end

set( gcf, 'color', 'w' );

if ~isempty( tag_auralization )
    filename = strcat( tag_auralization, '_coherence_loss_TF_comparison' );
    save_pdf = 1; save_png = 0; % 2D line plot
    export_figures( filename, save_mat_fig, save_png, save_pdf );
end

fprintf( ['\n- Coherence loss: transfer functions with/without coherence loss plotted ' ...
          '(closest point at t = %.1f s, strongest coherence loss at t = %.1f s, colorbar +/- %g dB)\n'], ...
    time(idx_closest), time(idx_max_w), cLim );

end

%% local function : il_symmetric_clim

function cLim = il_symmetric_clim(data, pct)
% symmetric colorbar limit [dB] from the <pct>-th percentile of |data|, rounded up
% to the next full dB (<= 5 dB) or the next multiple of 5 dB (> 5 dB), at least 1 dB.
% No toolbox needed (no <prctile>).
values = sort( abs( data(:) ) );
values = values( isfinite(values) );
if isempty( values )
    cLim = 1;
    return;
end
idx = max( 1, ceil( pct/100 * numel(values) ) );
x = max( values(idx), 1 ); % at least +/- 1 dB, even if the effect is negligible
if x <= 5
    cLim = ceil( x );
else
    cLim = 5 * ceil( x/5 );
end
end

%% local function : il_diverging_colormap

function cmap = il_diverging_colormap(n)
% blue - white - red colormap, symmetric around zero (no toolbox needed)
n2 = floor( n/2 );
blue  = [0.13 0.40 0.67];
red   = [0.70 0.09 0.17];
white = [1 1 1];
lower = [ linspace(blue(1), white(1), n2)', linspace(blue(2), white(2), n2)', linspace(blue(3), white(3), n2)' ];
upper = [ linspace(white(1), red(1), n2)',  linspace(white(2), red(2), n2)',  linspace(white(3), red(3), n2)' ];
cmap = [ lower; upper ];
end