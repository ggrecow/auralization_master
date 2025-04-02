function PLOT_HRIR_spectrogram(input_direct, input_reflected, dt_panam, tag_title, tag_auralization, tag_save)

global save_mat_fig

fontSize = 16;

stretchY = 1; % stretch plot in the vertical direction

h = figure( 'name', tag_title);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'Position', [pos(1), pos(2), pos(3), pos(4)*stretchY]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)*stretchY])
movegui(h, 'center');

% t = tiledlayout('horizontal');
t = tiledlayout(2,2);
t.TileSpacing = 'compact';

% direct - left ear
ax1 = nexttile;

data = input_direct.l_a;

timeVec = linspace(0, size(data,2).*dt_panam,size(data,2)) ;

% Your function for plotting
AKp(data, 'm3d', 'y', timeVec, 'dr', [-20 20]);

XLabelStr =  'Time, $t$ (s)';

% Set axis labels and titles with LaTeX interpreter
ax1.YLabel.String = 'Frequency, $f$ (Hz)';
ax1.YLabel.Interpreter = 'latex';
ax1.XTick = [0, floor(timeVec(round(length(timeVec)/2))), floor(timeVec(end))];
% ax1.XTickLabel =  timeVec;
ax1.XTickLabelRotation = 0;
ax1.XLabel.Interpreter = 'latex';
ax1.XLabel.String =  XLabelStr;
ax1.XLabel.Rotation = 0;
ax1.Title.String = 'Direct - left ear';
ax1.Title.FontWeight = 'normal';
ax1.Title.Interpreter = 'latex';

% xline(90, 'k--');

xticklabels(""); % Remove y axis numbering
xlabel(""); % Remove y axis label
colorbar('off'); % Remove the colorbar from the figure

% direct - right ear
ax2 = nexttile;

data = input_direct.r_a;

% Your function for plotting
% Your function for plotting
AKp(data, 'm3d', 'y', timeVec, 'dr', [-20 20]);

XLabelStr =  'Time, $t$ (s)';

% Set axis labels and titles with LaTeX interpreter
ax2.YLabel.String = 'Frequency, $f$ (Hz)';
ax2.YLabel.Interpreter = 'latex';
ax2.XTick = [0, floor(timeVec(round(length(timeVec)/2))), floor(timeVec(end))];
% ax2.XTickLabel = ["-90", "0", "90", "0", "-90"];
ax2.XTickLabelRotation = 0;
ax2.XLabel.Interpreter = 'latex';
ax2.XLabel.String =  XLabelStr;
ax2.XLabel.Rotation = 0;
ax2.Title.String = 'Direct - right ear';
ax2.Title.FontWeight = 'normal';
ax2.Title.Interpreter = 'latex';

yticklabels(""); % Remove y axis numbering
ylabel(""); % Remove y axis label
xticklabels(""); % Remove y axis numbering
xlabel(""); % Remove y axis label
colorbar('off'); % Remove the colorbar from the figure

% xline(90, 'k--');

% reflected - left ear
ax3 = nexttile;

data = input_reflected.l_a;

timeVec = linspace(0, size(data,2).*dt_panam,size(data,2)) ;

% Your function for plotting
AKp(data, 'm3d', 'y', timeVec, 'dr', [-20 20]);

XLabelStr =  'Time, $t$ (s)';

% Set axis labels and titles with LaTeX interpreter
ax3.YLabel.String = 'Frequency, $f$ (Hz)';
ax3.YLabel.Interpreter = 'latex';
ax3.XTick = [0, floor(timeVec(round(length(timeVec)/2))), floor(timeVec(end))];
% ax3.XTickLabel =  timeVec;
ax3.XTickLabelRotation = 0;
ax3.XLabel.Interpreter = 'latex';
ax3.XLabel.String =  XLabelStr;
ax3.XLabel.Rotation = 0;
ax3.Title.String = ' Reflected - left ear';
ax3.Title.FontWeight = 'normal';
ax3.Title.Interpreter = 'latex';

colorbar('off'); % Remove the colorbar from the figure

% xline(90, 'k--');

colorbar('off'); % Remove the colorbar from the figure

% reflected - right ear
ax4 = nexttile;

data = input_reflected.r_a;

timeVec = linspace(0, size(data,2).*dt_panam,size(data,2)) ;

% Your function for plotting
AKp(data, 'm3d', 'y', timeVec, 'dr', [-20 20]);

XLabelStr =  'Time, $t$ (s)';

% Set axis labels and titles with LaTeX interpreter
ax4.YLabel.String = 'Frequency, $f$ (Hz)';
ax4.YLabel.Interpreter = 'latex';
ax4.XTick = [0, floor(timeVec(round(length(timeVec)/2))), floor(timeVec(end))];
% ax1.XTickLabel =  timeVec;
ax4.XTickLabelRotation = 0;
ax4.XLabel.Interpreter = 'latex';
ax4.XLabel.String =  XLabelStr;
ax4.XLabel.Rotation = 0;
ax4.Title.String = 'Reflected - right ear';
ax4.Title.FontWeight = 'normal';
ax4.Title.Interpreter = 'latex';

% xline(90, 'k--');

yticklabels(""); % Remove y axis numbering
ylabel(""); % Remove y axis label

% common settings

% Create the colorbar and set the label
cb = colorbar; % Colorbar for this tile
zString = ('SPL, $L_{\mathrm{z}}$ (dB re 20 $\mu$Pa)');
set(get(cb,'label'),'string', zString, 'fontsize', fontSize, 'Interpreter', 'latex');

if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
else
    filename = strcat(tag_auralization, tag_save);
    save_pdf = 0; save_png = 1;
    export_figures( filename, save_mat_fig, save_png, save_pdf );
end

end