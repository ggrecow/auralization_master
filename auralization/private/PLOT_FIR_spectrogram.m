function PLOT_FIR_spectrogram( impulseResponse_FIR, FRF_ART, fs, tag_auralization )

global save_mat_fig

% ART
% P_ART= 10.*log10( abs( FRF_ART ).^2./ (2e-5)^2 );
FRF_ART= 10.*log10( abs( FRF_ART ).^2);
T = 0:size(FRF_ART,2);

%  compute freq vector for <FRF_ART>
frequencyVector_FRF_ART = linspace(0, fs/2, size( FRF_ART, 1 ));

% % original - ART
% h  = figure;
% set(h, 'name', 'PROCESSING - Atmospheric transfer function spectrogram (ART)');
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% imagesc(T, frequencyVector_FRF_ART, FRF_ART); colorbar('vert'); set(gca,'YDir','Normal')
% 
% ylabel('Frequency, $f$ (Hz)','Interpreter','Latex');
% xlabel('Time bins','Interpreter','Latex');
% ylabel(colorbar, '$10\cdot\log_{10}\left(|A_{\mathrm{atm,ART}}|^2 \right)$','Interpreter','Latex');
% 
% colormap('jet'); clim([-150 0]);
% % title('Original from ART');
% 
% set(gcf,'color','w');
% 
% if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
% else
%     if save_mat_fig==1
%         saveas(gcf,strcat(tag_auralization, '_atmospheric_transfer_function_spectrogram_ART'), 'fig');
%     else
%     end
%     % saveas(gcf,strcat(tag_auralization, '_atmospheric_transfer_function_spectrogram_ART'), 'pdf');
%     saveas(gcf,strcat(tag_auralization, '_atmospheric_transfer_function_spectrogram_ART'), 'png');
% end

% FIR

% get single-sided spectrum of <impulseResponse_FIR>
FR_cut = fft(impulseResponse_FIR);
FR_cut_single_side = FR_cut(1:ceil((size(FR_cut,1)+1)/2), :);

P_FIR= 10.*log10( abs( FR_cut_single_side ).^2.);

%  compute freq vector for <impulseResponse_FIR>
frequencyVector_FIR = linspace(0, fs/2, size( P_FIR, 1 ));

% h  = figure;
% set(h, 'name', 'PROCESSING - Atmospheric transfer function spectrogram (FIR filter)');
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% imagesc(T, frequencyVector_FIR, P_FIR); colorbar('vert'); set(gca,'YDir','Normal')
% colormap('jet'); clim([-150 0]);
% 
% ylabel('Frequency, $f$ (Hz)','Interpreter','Latex');
% xlabel('Time bins','Interpreter','Latex');
% ylabel(colorbar, '$10\cdot\log_{10}\left(|A_{\mathrm{atm,FIR}}|^2 \right)$','Interpreter','Latex');
% 
% % title('FIR filter');
% set(gcf,'color','w');
% 
% if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
% else
%     if save_mat_fig==1
%         saveas(gcf,strcat(tag_auralization, '_atmospheric_transfer_function_spectrogram_FIR'), 'fig');
%     else
%     end
%     % saveas(gcf,strcat(tag_auralization, '_atmospheric_transfer_function_spectrogram_FIR'), 'pdf');
%     saveas(gcf,strcat(tag_auralization, '_atmospheric_transfer_function_spectrogram_FIR'), 'png');
% end

% DELTA
if size( P_FIR,1 ) == size( FRF_ART,1 ) % only plots delta if sizes (nTaps) are equal

    D = P_FIR - FRF_ART;

    h  = figure;
    set(h, 'name', 'PROCESSING - Atmospheric transfer function spectrogram (FIR-ART)');
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    imagesc(T, frequencyVector_FIR, D); colorbar('vert'); set(gca,'YDir','Normal')
    colormap('jet'); %clim([0 50]);

    ylabel('Frequency, $f$ (Hz)','Interpreter','Latex');
    xlabel('Time bins','Interpreter','Latex');
    ylabel(colorbar, '$10\cdot\log_{10}\left(\frac{|A_{\mathrm{atm,FIR}}}{A_{\mathrm{atm,ART}}}|^2 \right)$','Interpreter','Latex');

    % title('ART minus FIR filter');

    if isempty(tag_auralization) % if tag_auralization is empty, dont save anything
    else
        filename = strcat(tag_auralization, '_atmospheric_transfer_function_spectrogram_FIR_minus_ART');
        save_pdf = 0; save_png = 1;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end

else
end

set(gcf,'color','w');

end