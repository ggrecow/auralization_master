
function export_figures( filename, save_mat_fig, save_png, save_pdf )
% function export_figures( filename, save_mat_fig, save_png, save_pdf )
%
% Function to export figures with a fixed resolution
%
% INPUTS: 
%       filename : string
%       contains the name of the file (without <.format>).
%
%       save_mat_fig, save_mat_png, save_mat_pdf : boolean
%       inform the function whether .fig, .png, and/or .pdf shoud be saved - 0=NO; 1; yes 
%
% Gil Felix Greco, Braunschweig 02.08.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    % save_mat_fig = 0;
    save_png = 1; % always save png
    save_pdf = 0;
end

resolution = '-r600'; % Default resolution of 600 DPI

if save_mat_fig == 1 % always good to save, but it may use a LOT OF SPACE is used everytime

    saveas( gcf, filename, 'fig' );

elseif save_png == 1 % preferred format for color plots like spectrograms

    format = '-dpng'; % png format
    print( gcf, filename, format, resolution );

elseif save_pdf == 1 % preferred format for normal 2D plots

    format = '-dpdf'; % pdf format
    print( gcf, filename, format, resolution );

end

end