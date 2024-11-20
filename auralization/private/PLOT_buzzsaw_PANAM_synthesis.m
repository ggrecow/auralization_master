
function PLOT_buzzsaw_PANAM_synthesis(freq_plot_panam, SPL_plot_panam, freq_synthesis, SPL_synthesis, f1, f2, source, tag_auralization)

global save_mat_fig

idx_plot_panam = freq_plot_panam >= f1 & freq_plot_panam <= f2;

nTonesPerToc_plot_panam = sum( idx_plot_panam~=0 );

idxx =  find( nTonesPerToc_plot_panam ~= 0);

figure('name','PROCESSING - Buzzsaw noise synthesis- overhead position' );
h  = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% colorplt = [.4, .4, .4];
colorplt = [ 0, 0, 0];

if length(idxx)<length( SPL_plot_panam )
    return;
else

% plot 1/3-OB as bars - PANAM predicted results for overhead position
for i = 1:length( SPL_plot_panam ) % this is not going to work if PANAM has more than 1 tone within the same band
    a = semilogx( [ f1(idxx(i)) f2(idxx(i)) ], [SPL_plot_panam(i) SPL_plot_panam(i)], 'Color', colorplt);hold on; % create horizontal lines
    semilogx( linspace( f1(idxx(i)), f1(idxx(i)), 24), linspace(0,SPL_plot_panam(i) ,24), 'Color', colorplt); % create vertical lines (flow)
    % semilogx(linspace(TOC(i,2),TOC(i,2),24),linspace(0,BandsVsTime(1,i),24),'b:'); % create vertical lines (fcenter)
    semilogx( linspace( f2(idxx(i)), f2(idxx(i)), 24 ), linspace(0,SPL_plot_panam(i),24), 'Color', colorplt); % create vertical lines (fupper)
end

end

colorstem = [ 0.49,0.18,0.56 ];
% semilogx( input{t_bin}.engine_buzzsaw(:,1), input{t_bin}.engine_buzzsaw(:,2) ); hold on;
b = stem( freq_synthesis, SPL_synthesis, 'Color', colorstem, 'LineWidth', 0.5, 'Marker', '*', MarkerSize=4);

ylabel('SPL, $L_{\mathrm{Z}}$ (dB re 20$~\mu$Pa)','Interpreter','Latex');
xlabel('Frequency, $f$ (Hz)', 'Interpreter' , 'latex');
% xlabel('1/3 octave nominal center frequency, $f$ (Hz)', 'Interpreter' , 'latex');
set(gcf,'color','w');
legend( [a,b], 'PANAM', 'Synthesis');

fnom = [ 25 31.5 40, 50 63 80, 100 125 160, 200 250 315, 400 500 630, ...
    800 1000 1250, 1600 2000 2500, 3150 4000 5000, 6300 8000 10000, ...
    12500 16000 20000 ] ;    % nominal center freq - preferred for freq. labeling

% xticks( fnom ) ;
% xtickangle(45);

h = gca; % Get axis to modify
h.XAxis.MinorTick = 'off'; % Must turn on minor ticks if they are off

% xlim([161 3811]);
% ylim([101 122]);

if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
else
    filename = strcat(tag_auralization, '_synthesised_overhead_', source);
    save_pdf = 1; save_png = 0;
    export_figures( filename, save_mat_fig, save_png, save_pdf );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% another plot, stored here to sabe space in the main code. to run this plot,
% please insert it just below this function's caller in the main code and run it
%
% plot predicted Buzzsaw SPL (1/3-OB) and distributed SPL over calculated BPF+nHarmonics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% for i = 1:nTimes
%     SplPlot(i) = input{i}.engine_buzzsaw(30,1);
% end
% plot( SplPlot ); hold on;
%
% % plot( tonesSPLTime(200, : ) ); hold on;
% ylabel('SPL (dB)', 'Interpreter' , 'latex'); xlabel('Frequency, $f$ (Hz)', 'Interpreter' , 'latex'); set(gcf,'color','w');
% legend('PANAM', 'Synthesis');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end