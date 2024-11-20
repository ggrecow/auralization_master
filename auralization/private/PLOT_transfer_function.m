function PLOT_transfer_function(xx, yy, data, tag_title)

figure( 'name', tag_title);
imagesc(xx, yy, data) ; colorbar('vert'); set(gca,'YDir','Normal')
colormap('jet'); 

% ylim([0 max(yy)]);
% caxis([0 max(max(data))]);
caxis([-100 -40]);

xlabel('Time, $t$ (s)', 'Interpreter', 'Latex');
ylabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel(colorbar, '$10\cdot\log_{10}\left(|H_{\mathrm{atm}}|^2 \right) $ (dB SPL)','Interpreter','Latex');
    
% load C:\Users\greco\Desktop\auralization_verification_thesis\private\COLORMAP_jet_white.mat
% colormap(jet_white);

ylim( [ 20 8e3 ] );

% yticks([ 20 1e2 1e3 10e3 20e3]);
% yticklabels({'20','100','1k', '10k', '20k'});

yticks([ 20 2e3 4e3 6e3 8e3 ]);
yticklabels({'20', '2k', '4k', '6k', '8k'});

% ax = gca;
% ax.XAxis.MinorTick = 'on';
% ax.XAxis.MinorTickValues =  0:0.5:50;

set(gcf,'color','w');

end